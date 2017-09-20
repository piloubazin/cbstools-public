package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainIntensityBasedSkullStripping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume mainImage;
	private ParamVolume extraImage;
	
	private ParamOption	 bgParam;
	private String		bgType = "exponential";
	private 	static final String[]	bgTypes = {"exponential","half-normal","exp+log-normal","half+log-normal"};
	private ParamBoolean	iterateParam;
	
	private		static final byte	HNORM = 101;
	private		static final byte	EXP = 102;
	private		static final byte	UNI = 21;
	private		static final byte	LOGN = 22;
	
	private ParamBoolean	skip0Param;
	private ParamBoolean	slab2dParam;
	private ParamInteger	dilateParam;
	
	private ParamVolume maskImage;
	private ParamVolume probaImage;
	private ParamVolume mainmaskedImage;
	private ParamVolume extramaskedImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(mainImage = new ParamVolume("Main Intensity Image"));
		
		inputParams.add(extraImage = new ParamVolume("extra image with high intensity at brain boundary (opt)"));
		extraImage.setMandatory(false);
		
		inputParams.add(bgParam = new ParamOption("Background noise model", bgTypes));
		bgParam.setValue(bgType);
		
		inputParams.add(iterateParam = new ParamBoolean("Iterative estimation", false));
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", false));
		inputParams.add(slab2dParam = new ParamBoolean("2D slab data", false));
		inputParams.add(dilateParam = new ParamInteger("Additional mask dilation", -5,5, 0));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("Intensity-based Skull Stripping");
		inputParams.setName("IntensityBasedSkullStripping");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimate a brain mask for a dataset with good brain/background intensity separation (e.g. PD-weighted). "
							+"An extra image can be used to ensure high intensities are preserved (e.g. T1 map, T2-weighted data or a probability map for a ROI");
		
		info.setVersion("3.1.1");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(maskImage = new ParamVolume("Brain Mask Image",VoxelType.UBYTE));
		outputParams.add(probaImage = new ParamVolume("Foreground Probability Image",VoxelType.FLOAT));
		outputParams.add(mainmaskedImage = new ParamVolume("Masked Main Image",VoxelType.FLOAT));
		outputParams.add(extramaskedImage = new ParamVolume("Masked Extra Image",VoxelType.FLOAT));
		
		extramaskedImage.setMandatory(false);
		
		outputParams.setName("skull stripped images");
		outputParams.setLabel("skull stripped images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		int main = 0;
		int nimg = 1;
		
		int extra = -1;
		if (extraImage.getImageData() != null) { extra = nimg; nimg++; }
		
		ImageDataFloat mainImg = new ImageDataFloat(mainImage.getImageData());
		ImageDataFloat extraImg = null;
		if (extra>-1) extraImg = new ImageDataFloat(extraImage.getImageData());
		
		float[][][] buffer = mainImg.toArray3d();
		int nx = mainImg.getRows();
		int ny = mainImg.getCols();
		int nz = mainImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = mainImg.getHeader().getDimResolutions()[0];
		float ry = mainImg.getHeader().getDimResolutions()[1];
		float rz = mainImg.getHeader().getDimResolutions()[2];
		
		float[][] image = new float[nimg][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[main][xyz] = buffer[x][y][z];
		}
		mainImg = null;
		buffer = null;
		
		if (extra>-1) {
			buffer = extraImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[extra][xyz] = buffer[x][y][z];
			}
			extraImg = null;
			buffer = null;
		}			

		// main algorithm

		boolean is2d = slab2dParam.getValue().booleanValue();
		/*
		byte off = 4;
		if (is2d) {
			// pad all images to avoid boundary problems	
			int nxyzp = (nx+2*off)*(ny+2*off)*(nz+2*off);
			float[][] tmp = new float[nimg][nxyzp];
			for (int i=0;i<nimg;i++) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x + nx*y + nx*ny*z;
					int xyzp = x+off + (nx+2*off)*(y+off) + (nx+2*off)*(ny+2*off)*(z+off);
					tmp[i][xyzp] = image[i][xyz];
				}
			}
			image = tmp;
			nx = nx+2*off;
			ny = ny+2*off;
			nz = nz+2*off;
			nxyz = nxyzp;
		}
		*/
		boolean skip0 = skip0Param.getValue().booleanValue();
		
		int ndata=0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[main][xyz]!=0) ndata++;
		
		double[] data = new double[ndata];
		int n=0;
		double max=0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[main][xyz]!=0) {
			data[n] = image[main][xyz];
			if (data[n]>max) max = data[n];
			n++;
		}
		Percentile measure = new Percentile();
		measure.setData(data, 0, ndata);
		max = measure.evaluate(99.8);
		
		double mean = 0.0;		
		double pi = Math.PI;
		
		byte model = EXP;
		if (bgParam.getValue().startsWith("half")) model = HNORM;
		byte fgmodel = UNI;
		if (bgParam.getValue().endsWith("log-normal")) fgmodel = LOGN;
		boolean iterate = iterateParam.getValue().booleanValue();
		
		if (model==EXP) mean = ImageStatistics.robustExponentialFit(data, iterate, ndata);
		else if (model==HNORM) mean = ImageStatistics.robustHalfGaussianFit(data, iterate, ndata);
		
		BasicInfo.displayMessage("background mean: "+mean+"\n");
		
		double fmean = 0.0;
		double fdev = 0.0;
		if (fgmodel==LOGN) {
			// model the filter response as something more interesting, e.g. log-normal (removing the bg samples)
			int nb=0;
			for (int xyz=0;xyz<nxyz;xyz++) if (image[main][xyz]>0) {
				// only count positive values
				nb++;
			}
			double[] response = new double[nb];
			n=0;
			for (int xyz=0;xyz<nxyz;xyz++) if (image[main][xyz]>0) {
				response[n] = image[main][xyz];
				n++;
			}
			double[] weights = new double[nb];
			for (int b=0;b<nb;b++) { 
				if (model==EXP) weights[b] = (1.0-FastMath.exp( -response[b]/mean));
				else if (model==HNORM)  weights[b] = (1.0-Math.exp(-response[b]*response[b]/(pi*mean*mean)));
				else weights[b] = 1.0;
				response[b] = FastMath.log(response[b]);
			}
		
			fmean = ImageStatistics.weightedPercentile(response,weights,50.0,nb);
		
			// stdev: 50% +/- 1/2*erf(1/sqrt(2)) (~0.341344746..)
			double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
			fdev = 0.5*(ImageStatistics.weightedPercentile(response,weights,50.0+dev,nb) - ImageStatistics.weightedPercentile(response,weights,50.0-dev,nb));
		
			BasicInfo.displayMessage("Log-normal parameter estimates: mean = "+FastMath.exp(fmean)+", stdev = "+FastMath.exp(fdev)+",\n");
		}
		
		// re-normalized probability map
		float[] proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-image[main][xyz]*image[main][xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(Math.exp(-image[main][xyz]/mean)/mean);
				
				
			// bg vs. outlier test : proba is p(xyz is background)
			if (fgmodel==UNI) proba[xyz] = (float)max*proba[xyz]/(1.0f+(float)max*proba[xyz]);
			else if (fgmodel==LOGN) {
				double plg = FastMath.exp(-Numerics.square(FastMath.log(image[main][xyz])-fmean)/(2.0*fdev*fdev))/FastMath.sqrt(2.0*FastMath.PI*fdev*fdev);
				proba[xyz] = (float)(proba[xyz]/(plg+proba[xyz]));
			}
		}

		BasicInfo.displayMessage("background-based skull stripping");
		
		// start from the bg mask
		
		MinMaxFiltering minmax = new MinMaxFiltering(proba, nx,ny,nz, rx,ry,rz);
		
		float[] brain = minmax.growRegion(new float[]{0.0f}, new float[]{0.9f}, new float[]{0.9f}, 16, 10, is2d);
		
		BinaryTopology topo = null;
		Gdm3d gdm = null;
		if (extra>-1) {
			// final refinement using extra contrast
			float extramax = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (extra>-1 && image[extra][xyz]>extramax) extramax = image[extra][xyz];
			}
			// keep CSF inside + regions away from boundary (=> limit number of iterations)
			int[] mask = new int[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (brain[xyz]>0) {
					mask[xyz] = 1;
				} else {
					mask[xyz] = 0;
				}
			}
			
			float[] balloon = new float[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				float force = -0.1f;
				if (extra>-1) force = Numerics.max(force, Numerics.bounded( (image[extra][xyz]/extramax-0.5f)/0.5f, -1.0f, 1.0f));
				balloon[xyz] = force;
			}
			
			
			if (is2d) {
				// cheap option: just add mirrored values and an extra boundary in the slab direction...
				int[] maskb = new int[nx*ny*(nz+12)];
				float[] balloonb = new float[nx*ny*(nz+12)];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					maskb[x+nx*y+nx*ny*(z+6)] = mask[x+nx*y+nx*ny*z];
					balloonb[x+nx*y+nx*ny*(z+6)] = balloon[x+nx*y+nx*ny*z];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<6;z++) {
					maskb[x+nx*y+nx*ny*z] = mask[x+nx*y+nx*ny*0];
					balloonb[x+nx*y+nx*ny*z] = balloon[x+nx*y+nx*ny*0];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=nz+6;z<nz+12;z++) {
					maskb[x+nx*y+nx*ny*z] = mask[x+nx*y+nx*ny*(nz-1)];
					balloonb[x+nx*y+nx*ny*z] = balloon[x+nx*y+nx*ny*(nz-1)];
				}

				// topology correction for the mask?
				topo = new BinaryTopology(maskb, nx, ny, nz+12, rx, ry, rz, "wcs");
			
				topo.outsideSphericalTopology();
			
				int[] toposeg = topo.exportIntSegmentation();

				gdm = new Gdm3d(toposeg, nx, ny, nz+12, rx, ry, rz, null, balloonb, 0.0f, 0.1f, 0.9f, "wcs");
						
				gdm.evolveNarrowBand(200, 0.0005f);
				
				float[] brainb = gdm.exportSegmentation();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					brain[x+nx*y+nx*ny*z] = brainb[x+nx*y+nx*ny*(z+6)];
				}
			} else {
				// topology correction for the mask?
				topo = new BinaryTopology(mask, nx, ny, nz, rx, ry, rz, "wcs");
				
				topo.outsideSphericalTopology();
				
				int[] toposeg = topo.exportIntSegmentation();
				
				gdm = new Gdm3d(toposeg, nx, ny, nz, rx, ry, rz, null, balloon, 0.0f, 0.1f, 0.9f, "wcs");
						
				gdm.evolveNarrowBand(200, 0.0005f);

				brain = gdm.exportSegmentation();
			}
		}
		
		byte[][][] brainmask = new byte[nx][ny][nz];
		float[][][] probamask = new float[nx][ny][nz];	
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			probamask[x][y][z] = 1.0f-proba[xyz];
			if (brain[xyz]>0) {
				brainmask[x][y][z] = 1;
			} else {
				brainmask[x][y][z] = 0;
			}
		}
		
		int dilate = dilateParam.getValue().intValue();
		if (dilate>0) {
			brainmask = Morphology.dilateObject(brainmask, nx,ny,nz, dilate,dilate,dilate);
		} else if (dilate<0) {
			brainmask = Morphology.erodeObject(brainmask, nx,ny,nz, -dilate,-dilate,-dilate);
		}

		/*		
		if (is2d) {
			// un-pad all images
			int nxyzu = (nx-2*off)*(ny-2*off)*(nz-2*off);
			float[][] tmp = new float[nimg][nxyzu];
			float[][][] tmp3d = new float[nx-2*off][ny-2*off][nz-2*off];
			byte[][][] tmp3db = new byte[nx-2*off][ny-2*off][nz-2*off];
			for (int i=0;i<nimg;i++) {
				for (int x=off;x<nx-2*off;x++) for (int y=off;y<ny-2*off;y++) for (int z=off;z<nz-2*off;z++) {
					int xyz = x + nx*y + nx*ny*z;
					int xyzu = x-off + (nx-2*off)*(y-off) + (nx-2*off)*(ny-2*off)*(z-off);
					tmp[i][xyzu] = image[i][xyz];
				}
			}
			for (int x=off;x<nx-2*off;x++) for (int y=off;y<ny-2*off;y++) for (int z=off;z<nz-2*off;z++) {
				tmp3d[x-off][y-off][z-off] = probamask[x][y][z];
				tmp3db[x-off][y-off][z-off] = brainmask[x][y][z];
			}
			image = tmp;
			probamask = tmp3d;
			brainmask = tmp3db;
			nx = nx-2*off;
			ny = ny-2*off;
			nz = nz-2*off;
			nxyz = nxyzu;
		}
		*/
		
		ImageDataUByte resultData = new ImageDataUByte(brainmask);		
		resultData.setHeader(mainImage.getImageData().getHeader());
		resultData.setName(mainImage.getImageData().getName()+"_stripmask");
		maskImage.setValue(resultData);
		resultData = null;
			
		ImageDataFloat bufferData = new ImageDataFloat(probamask);
		bufferData.setHeader(mainImage.getImageData().getHeader());
		bufferData.setName(mainImage.getImageData().getName()+"_probamask");
		probaImage.setValue(bufferData);
		bufferData = null;
		
		buffer = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (brainmask[x][y][z]>0) buffer[x][y][z] = image[main][xyz];
			else buffer[x][y][z] = 0.0f;
		}
		bufferData = new ImageDataFloat(buffer);
		bufferData.setHeader(mainImage.getImageData().getHeader());
		bufferData.setName(mainImage.getImageData().getName()+"_strip");
		mainmaskedImage.setValue(bufferData);
		bufferData = null;
		buffer = null;
		
		if (extra>-1) {
			buffer = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (brainmask[x][y][z]>0) buffer[x][y][z] = image[extra][xyz];
				else buffer[x][y][z] = 0.0f;
			}
			bufferData = new ImageDataFloat(buffer);
			bufferData.setHeader(extraImage.getImageData().getHeader());
			bufferData.setName(extraImage.getImageData().getName()+"_strip");
			extramaskedImage.setValue(bufferData);
			bufferData = null;
			buffer = null;
		}
	}


}
