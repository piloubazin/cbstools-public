package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
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

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMp2rageDuraEstimation extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	pvtypes = {"dura_region","boundary","dura_prior","bg_prior", "intens_prior"};
	private		String		pvtype = "dura_region";
	
	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private ParamOption pvParam;
	private ParamDouble distParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) Image"));
		inputParams.add(maskImage = new ParamVolume("Skull Stripping Mask"));
		inputParams.add(distParam = new ParamDouble("Distance to background (mm)", -1E10, 1E10, 5.0));
		inputParams.add(pvParam = new ParamOption("output type", pvtypes));
		pvParam.setDescription("Outputs an estimate of the dura / CSF boundary or an estimate of the entire dura region.");
		pvParam.setValue(pvtype);

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("MP2RAGE Dura Estimation");
		inputParams.setName("Mp2rageDuraEstimation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Filters a MP2RAGE brain image to obtain a probability map of dura matter.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Dura Image",VoxelType.FLOAT));
		outputParams.setName("dura_image");
		outputParams.setLabel("Dura Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat inv2 = new ImageDataFloat(inv2Image.getImageData());
		float[][][] inv2img = inv2.toArray3d();
		int nx = inv2.getRows();
		int ny= inv2.getCols();
		int nz = inv2.getSlices();
		float rx = inv2.getHeader().getDimResolutions()[0];
		float ry = inv2.getHeader().getDimResolutions()[1];
		float rz = inv2.getHeader().getDimResolutions()[2];
		
		ImageDataUByte mask = new ImageDataUByte(maskImage.getImageData());
		byte[][][] maskimg = mask.toArray3d();
		
		// main algorithm
		
		// 1. Estimate regions of potential partial volume
		float[][][] pvmap = new float[nx][ny][nz];
		float pvavg = 0.0f;
		float pvsig = 0.0f;
		float pvden = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (!zeroNeighbor(inv2img,x,y,z,2)) {
				// check for zero-valued neighbors as well
				pvmap[x][y][z] = 0.0f;
				float best = 0.0f;
				float dmax = 13;
				for (int d=0;d<dmax;d++) {
					float val = planeScore(inv2img, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				pvmap[x][y][z] = best;
				pvavg += best;
				pvden++;
			}
		}
		pvavg /= pvden;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			//if (!zeroNeighbor(inv2img,x,y,z,2)) {
			pvsig += Numerics.square(pvmap[x][y][z]-pvavg);
		}
		pvsig /= (pvden-1.0f);
		// probability
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (pvmap[x][y][z]>pvavg) {
				pvmap[x][y][z] = 1.0f - (float)Math.exp( -0.5f*(pvmap[x][y][z]-pvavg)*(pvmap[x][y][z]-pvavg)/pvsig );
			} else {
				pvmap[x][y][z] = 0.0f;
			}
		}

		// 2. Expand background region
		float scale = distParam.getValue().floatValue()/rx;
		float[][][] bgmap = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			bgmap[x][y][z] = maskimg[x][y][z];
		}
		float[][] kernel = ImageFilters.separableGaussianKernel(scale, scale, scale);
		int kx = (kernel[0].length-1)/2;
		bgmap = ImageFilters.separableConvolution(bgmap, nx, ny, nz, kernel, kx, kx, kx);
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			bgmap[x][y][z] = 2.0f*(1.0f-bgmap[x][y][z])*maskimg[x][y][z];
		}

		// 3. Find regions of low intensity
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 20;
		
		// min max for t1 map fixing
		float max = -1e10f;
		float min = 1e10f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (inv2img[x][y][z]>max) max = inv2img[x][y][z];
			if (inv2img[x][y][z]<min) min = inv2img[x][y][z];
		}
		// re-normalize
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			inv2img[x][y][z] = (inv2img[x][y][z]-min)/(max-min);
		}
		System.out.println("image range: "+min+", "+max);

		double mean = 0.0;
		double den = 0.0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (inv2img[x][y][z]!=0) {
				mean += inv2img[x][y][z];
				den++;
			}
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		// re-normalized probability map
		float[][][] intensmap = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			intensmap[x][y][z] = (float)(Math.exp(-inv2img[x][y][z]/mean)/mean);
			intensmap[x][y][z] = intensmap[x][y][z]/(1.0f+intensmap[x][y][z]);
		}
		// loop
		double diff = 1.0;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			diff = mean;
			mean = 0.0;
			den = 0.0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (inv2img[x][y][z]!=0) {
					mean += intensmap[x][y][z]*inv2img[x][y][z];
					den += intensmap[x][y][z];
				}
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(diff-mean);
			System.out.println("diff parameters: "+diff);
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				intensmap[x][y][z] = (float)(Math.exp(-(inv2img[x][y][z])/mean)/mean);
				intensmap[x][y][z] = intensmap[x][y][z]/(1.0f+intensmap[x][y][z]);
			}
		}

		// 4. Combine information sources, expand PV in regions of higher BG
		float[][][] result = new float[nx][ny][nz];
		ImageDataFloat resultData;
		if (pvParam.getValue().equals("dura_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x][y][z] = pvmap[x][y][z];
			}
		} else 
		if (pvParam.getValue().equals("bg_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x][y][z] = Numerics.min(1.0f, bgmap[x][y][z]);
			}
		} else 
		if (pvParam.getValue().equals("intens_prior")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x][y][z] = intensmap[x][y][z];
			}
		} else 
		if (pvParam.getValue().equals("boundary")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				//result[x][y][z] = (float)Math.sqrt(bgmap[x][y][z]*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]));
				result[x][y][z] = Numerics.min(1.0f, bgmap[x][y][z])*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]);
			}
		} else 
		if (pvParam.getValue().equals("dura_region")) {
			// take the product
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				//result[x][y][z] = (float)Math.sqrt(bgmap[x][y][z]*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]));
				result[x][y][z] = Numerics.min(1.0f, bgmap[x][y][z])*Numerics.max(pvmap[x][y][z],intensmap[x][y][z]);
			}
			// propagate maximum values toward background
			float[][][] tmp = new float[nx][ny][nz];
			for (int t=0;t<2*scale;t++) {
				for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
					tmp[x][y][z] = result[x][y][z];
					if (bgmap[x+1][y][z]<bgmap[x][y][z] && result[x+1][y][z]>tmp[x][y][z]) tmp[x][y][z] = result[x+1][y][z];
					if (bgmap[x-1][y][z]<bgmap[x][y][z] && result[x-1][y][z]>tmp[x][y][z]) tmp[x][y][z] = result[x-1][y][z];
					if (bgmap[x][y+1][z]<bgmap[x][y][z] && result[x][y+1][z]>tmp[x][y][z]) tmp[x][y][z] = result[x][y+1][z];
					if (bgmap[x][y-1][z]<bgmap[x][y][z] && result[x][y-1][z]>tmp[x][y][z]) tmp[x][y][z] = result[x][y-1][z];
					if (bgmap[x][y][z+1]<bgmap[x][y][z] && result[x][y][z+1]>tmp[x][y][z]) tmp[x][y][z] = result[x][y][z+1];
					if (bgmap[x][y][z-1]<bgmap[x][y][z] && result[x][y][z-1]>tmp[x][y][z]) tmp[x][y][z] = result[x][y][z-1];
				}
				for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
					result[x][y][z] = tmp[x][y][z];
				}
			}
			/*
			// square the result to dampen lower values?
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				result[x][y][z] = result[x][y][z]*result[x][y][z];
			}
			*/
		}
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(inv2.getHeader());
		resultData.setName(inv2.getName()+"_dura");
		resultImage.setValue(resultData);
		resultData = null;
		result = null;
		
	}

	boolean zeroNeighbor(float[][][] image, int x, int y, int z, int d) {
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			if (image[x+i][y+j][z+l]!=image[x][y][z] && i*i+j*j+l*l<=2*d*d) return false;
		}
		return true;
	}
	
	float planeScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {		
			val =(2.0f*image[x][y][z]-image[x-1][y][z]-image[x+1][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z]-image[x+1][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y][z-1]-image[x+1][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y][z+1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z-1]-image[x-1][y-1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x][y-1][z+1]-image[x-1][y-1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x-1][y+1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x][y+1][z+1]-image[x-1][y+1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==1) {
			val =(2.0f*image[x][y][z]-image[x][y-1][z]-image[x][y+1][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z]-image[x-1][y+1][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x][y-1][z-1]-image[x][y+1][z-1]
				 +2.0f*image[x][y][z+1]-image[x][y-1][z+1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y][z-1]-image[x-1][y-1][z-1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y][z+1]-image[x-1][y-1][z+1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x+1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y][z+1]-image[x+1][y-1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==2) { 			
			val =(2.0f*image[x][y][z]-image[x][y][z-1]-image[x][y][z+1]
				 +2.0f*image[x-1][y][z]-image[x-1][y][z-1]-image[x-1][y][z+1]
				 +2.0f*image[x+1][y][z]-image[x+1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z]-image[x][y-1][z-1]-image[x][y-1][z+1]
				 +2.0f*image[x][y+1][z]-image[x][y+1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y-1][z]-image[x-1][y-1][z-1]-image[x-1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z]-image[x-1][y+1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z]-image[x+1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y+1][z]-image[x+1][y+1][z-1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==3) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x-1][y+1][z]-image[x-2][y][z]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]-image[x][y+2][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y-1][z+1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]-image[x][y+2][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==4) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x][y-1][z+1]-image[x][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]-image[x-1][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]-image[x-1][y][z+2]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]-image[x+1][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]-image[x+1][y][z+2])/18.0f;
		} else if (d==5) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x][y][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y][z+1]-image[x-2][y][z]-image[x][y][z+2]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]-image[x+2][y-1][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]-image[x][y-1][z+2]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]-image[x+2][y+1][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]-image[x][y+1][z+2])/18.0f;
		} else if (d==6) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z]-image[x+1][y-1][z]
				 +2.0f*image[x-1][y-1][z]-image[x-2][y][z]-image[x][y-2][z]
				 +2.0f*image[x+1][y+1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y+1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]-image[x][y-2][z-1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y+1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]-image[x][y-2][z+1]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==7) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z+1]-image[x][y+1][z-1]
				 +2.0f*image[x][y-1][z-1]-image[x][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x][y+1][z+1]-image[x][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z+1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]-image[x-1][y][z-2]
				 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]-image[x-1][y+2][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]-image[x+1][y][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]-image[x+1][y+2][z])/18.0f;
		} else if (d==8) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z+1]-image[x+1][y][z-1]
				 +2.0f*image[x-1][y][z-1]-image[x-2][y][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y][z+1]-image[x][y][z+2]-image[x+2][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z+1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]-image[x][y-1][z-2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]-image[x+2][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]-image[x][y+1][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]-image[x+2][y+1][z])/18.0f;
		} else if (d==9) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]-image[x+2][y][z]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]-image[x+2][y+2][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]-image[x][y+2][z+2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==10) { 			
			val =(2.0f*image[x][y][z]-image[x+1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]-image[x-2][y][z]
				 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]-image[x-2][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]-image[x-2][y][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]-image[x][y+2][z+2])/14.0f;
		} else if (d==11) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]-image[x][y-2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]-image[x+2][y-2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]-image[x][y-2][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==12) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]-image[x+2][y][z-2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]-image[x][y+2][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]-image[x+2][y+2][z])/14.0f;
		}
		return 0.5f*val;
	}

}
