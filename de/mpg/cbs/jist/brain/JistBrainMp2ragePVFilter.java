package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;


/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMp2ragePVFilter extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputT1mapImage;
	private ParamVolume inputT1weiImage;
	private ParamVolume resultImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputT1mapImage = new ParamVolume("T1 map (T1_images) Image (opt)"));
		inputParams.add(inputT1weiImage = new ParamVolume("T1-weighted (UNI) Image (opt)"));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("MP2RAGE Partial Volume Filter");
		inputParams.setName("Mp2ragePVFilter");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Filters a MP2RAGE image for regions of partial voluming assuming a ridge-like model of intensity.");
		
		info.setVersion("3.0.1");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Partial Volume Image",VoxelType.FLOAT));
		outputParams.setName("pv_image");
		outputParams.setLabel("PV Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat in = null;
		String pvType = null;
		String outType = "probability";
		if (inputT1mapImage.getImageData()!=null) {
			in = new ImageDataFloat(inputT1mapImage.getImageData());
			pvType = "bright";
		} else if (inputT1weiImage.getImageData()!=null) {
			in = new ImageDataFloat(inputT1weiImage.getImageData());
			pvType = "dark";
		}			
		float[][][] image = in.toArray3d();
		int nx = in.getRows();
		int ny= in.getCols();
		int nz = in.getSlices();
		
		// main algorithm
		float[][][] result = new float[nx][ny][nz];
		ImageDataFloat resultData;
		float avg = 0.0f;
		float sig = 0.0f;
		float den = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
			// check for zero-valued neighbors as well
			result[x][y][z] = 0.0f;
			float best = 0.0f;
			float dmax = 13;
			for (int d=0;d<dmax;d++) {
				float val = planeScore(image, x,y,z,d);
				if (val*val>best*best) best = val;
			}
			result[x][y][z] = best;
			avg += best;
			den++;
		}
		avg /= den;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
			sig += (result[x][y][z]-avg)*(result[x][y][z]-avg);
			
			if ( (pvType.equals("bright") && result[x][y][z]<0) || (pvType.equals("dark") && result[x][y][z]>0) ) {
				result[x][y][z] = 0.0f;
			} else if  (pvType.equals("dark")) {
				result[x][y][z] = -result[x][y][z];
			}
		}
		sig /= (den-1.0f);
		
		if (outType.equals("probability")) {
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				result[x][y][z] = 1.0f - (float)Math.exp( -0.5f*(result[x][y][z]-avg)*(result[x][y][z]-avg)/sig );
			}
		}
		
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(in.getHeader());
		resultData.setName(in.getName()+"_pv");
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
			val =(2.0f*image[x][y][z]		-image[x-1][y][z]		-image[x+1][y][z]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x][y-1][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x][y+1][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x][y+1][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==1) {
			val =(2.0f*image[x][y][z]		-image[x][y-1][z]		-image[x][y+1][z]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1]
				 +2.0f*image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1]
				 +2.0f*image[x-1][y][z-1]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y][z+1]	-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y][z-1]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y][z+1]	-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==2) { 			
			val =(2.0f*image[x][y][z]		-image[x][y][z-1]		-image[x][y][z+1]
				 +2.0f*image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1]
				 +2.0f*image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1]
				 +2.0f*image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1]
				 +2.0f*image[x-1][y-1][z]	-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z]	-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z]	-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y+1][z]	-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==3) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x-1][y+1][z]	-image[x-2][y][z]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z]	-image[x][y-2][z]		-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]		-image[x][y+2][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]		-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]		-image[x][y+2][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]		-image[x+2][y][z+1])/18.0f;
		} else if (d==4) { 			
			val =(2.0f*image[x][y][z]		-image[x][y-1][z-1]		-image[x][y+1][z+1]
				 +2.0f*image[x][y+1][z-1]	-image[x][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x][y-1][z+1]	-image[x][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]		-image[x-1][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]		-image[x-1][y][z+2]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]		-image[x+1][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]		-image[x+1][y][z+2])/18.0f;
		} else if (d==5) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y][z-1]		-image[x+1][y][z+1]
				 +2.0f*image[x+1][y][z-1]	-image[x][y][z-2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y][z+1]	-image[x-2][y][z]		-image[x][y][z+2]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]		-image[x+2][y-1][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]		-image[x][y-1][z+2]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]		-image[x+2][y+1][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]		-image[x][y+1][z+2])/18.0f;
		} else if (d==6) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y+1][z]		-image[x+1][y-1][z]
				 +2.0f*image[x-1][y-1][z]	-image[x-2][y][z]		-image[x][y-2][z]
				 +2.0f*image[x+1][y+1][z]	-image[x][y-2][z]		-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]		-image[x][y-2][z-1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]		-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]		-image[x][y-2][z+1]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]		-image[x+2][y][z+1])/18.0f;
		} else if (d==7) { 			
			val =(2.0f*image[x][y][z]		-image[x][y-1][z+1]		-image[x][y+1][z-1]
				 +2.0f*image[x][y-1][z-1]	-image[x][y-2][z]		-image[x][y][z-2]
				 +2.0f*image[x][y+1][z+1]	-image[x][y][z+2]		-image[x][y+2][z]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]		-image[x-1][y][z-2]
				 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]		-image[x-1][y+2][z]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]		-image[x+1][y][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]		-image[x+1][y+2][z])/18.0f;
		} else if (d==8) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y][z+1]		-image[x+1][y][z-1]
				 +2.0f*image[x-1][y][z-1]	-image[x-2][y][z]		-image[x][y][z-2]
				 +2.0f*image[x+1][y][z+1]	-image[x][y][z+2]		-image[x+2][y][z]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]		-image[x][y-1][z-2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]		-image[x+2][y-1][z]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]		-image[x][y+1][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]		-image[x+2][y+1][z])/18.0f;
		} else if (d==9) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]		-image[x+2][y][z]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]		-image[x+2][y+2][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]		-image[x][y+2][z+2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]		-image[x+2][y][z+2])/14.0f;
		} else if (d==10) { 			
			val =(2.0f*image[x][y][z]		-image[x+1][y-1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]		-image[x-2][y][z]
				 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]		-image[x-2][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]		-image[x-2][y][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]		-image[x][y+2][z+2])/14.0f;
		} else if (d==11) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]		-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]		-image[x][y-2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]		-image[x+2][y-2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]		-image[x][y-2][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]		-image[x+2][y][z+2])/14.0f;
		} else if (d==12) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]		-image[x][y][z-2]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]		-image[x+2][y][z-2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]		-image[x][y+2][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]		-image[x+2][y+2][z])/14.0f;
		}
		return 0.5f*val;
	}

}
