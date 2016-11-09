package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
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


/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainIntensityBasedSkullStripping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume mainImage;
	private ParamVolume extraImage;
	
	private ParamOption	 bgParam;
	private String		bgType = "exponential";
	private 	static final String[]	bgTypes = {"exponential","half-normal"};
	private ParamBoolean	iterateParam;
	
	private		static final byte	HNORM = 101;
	private		static final byte	EXP = 102;
	
	private ParamBoolean	skip0Param;
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
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("Intensity-based Skull Stripping");
		inputParams.setName("IntensityBasedSkullStripping");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimate a brain mask for a dataset with good brain/background intensity separation (e.g. PD-weighted). "
							+"An extra image can be used to ensure high intensities are preserved (e.g. T1 map, T2-weighted data or a probability map for a ROI");
		
		info.setVersion("3.0.6");
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
		double mean = 0.0;		
		double pi = Math.PI;
		
		byte model = EXP;
		if (bgParam.getValue().equals("half-normal")) model = HNORM;
		boolean iterate = iterateParam.getValue().booleanValue();
		
		if (model==EXP) mean = ImageStatistics.robustExponentialFit(data, iterate, ndata);
		else if (model==HNORM) mean = ImageStatistics.robustHalfGaussianFit(data, iterate, ndata);
		
		Interface.displayMessage("background mean: "+mean);
		
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
			proba[xyz] = proba[xyz]/(1.0f+(float)max*proba[xyz]);
		}

		Interface.displayMessage("background-based skull stripping");
		
		// start from the bg mask
		MinMaxFiltering minmax = new MinMaxFiltering(proba, nx,ny,nz, rx,ry,rz);
		
		float[] brain = minmax.growRegion(new float[]{0.0f}, new float[]{0.9f}, new float[]{0.9f}, 16, 10);
		
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
			
			// topology correction for the mask?
			topo = new BinaryTopology(mask, nx, ny, nz, rx, ry, rz, "wcs");
			
			topo.outsideSphericalTopology();
			
			
			gdm = new Gdm3d(topo.exportIntSegmentation(), nx, ny, nz, rx, ry, rz, null, balloon, 0.0f, 0.1f, 0.9f, "wcs");
						
			gdm.evolveNarrowBand(100, 0.001f);

			brain = gdm.exportSegmentation();
		}
		
		byte[][][] brainmask = new byte[nx][ny][nz];
		float[][][] probamask = new float[nx][ny][nz];	
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (brain[xyz]>0) {
				brainmask[x][y][z] = 1;
				probamask[x][y][z] = 1.0f-proba[xyz];
			} else {
				brainmask[x][y][z] = 0;
				probamask[x][y][z] = 0.0f;
			}
		}

		ImageDataUByte resultData = new ImageDataUByte(brainmask);		
		resultData.setHeader(mainImage.getImageData().getHeader());
		resultData.setName(mainImage.getImageData().getName()+"_stripmask");
		maskImage.setValue(resultData);
		resultData = null;
		brainmask = null;
		
		ImageDataFloat bufferData = new ImageDataFloat(probamask);
		bufferData.setHeader(mainImage.getImageData().getHeader());
		bufferData.setName(mainImage.getImageData().getName()+"_probamask");
		probaImage.setValue(bufferData);
		bufferData = null;
		
		buffer = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (brain[xyz]>0) buffer[x][y][z] = image[main][xyz];
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
				if (brain[xyz]>0) buffer[x][y][z] = image[extra][xyz];
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
