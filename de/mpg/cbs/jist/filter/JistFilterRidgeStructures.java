package de.mpg.cbs.jist.filter;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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

import de.mpg.cbs.utilities.Numerics;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistFilterRidgeStructures extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	pvtypes = {"bright","dark","both"};
	private		String		pvtype = "bright";
	private		static final String[]	outtypes = {"probability","intensity"};
	private		String		outtype = "probability";
	
	// jist containers
	private ParamVolume inputImage;
	private ParamVolume resultImage;
	private ParamOption pvParam;
	private ParamOption outParam;
	private ParamBoolean strictParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(pvParam = new ParamOption("Structure intensity", pvtypes));
		pvParam.setDescription("Filter specifically for bright or dark intensity structures.");
		pvParam.setValue(pvtype);
		inputParams.add(outParam = new ParamOption("output", outtypes));
		pvParam.setDescription("Outputs the raw intensity values or a probability score for the detected structures.");
		pvParam.setValue(outtype);
		inputParams.add(strictParam = new ParamBoolean("use strict min/max filter", true));
		strictParam.setDescription("Using either a broader response filter or a very strict min/max filter");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Filter");
		inputParams.setLabel("Ridge Structures");
		inputParams.setName("RidgeStructures");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Filters an image for 3D ridge-like structures.");
		
		info.setVersion("3.0.4");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Ridge Structure Image",VoxelType.FLOAT));
		outputParams.setName("ridge_image");
		outputParams.setLabel("Ridge Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat in = new ImageDataFloat(inputImage.getImageData());
		float[][][] image = in.toArray3d();
		int nx = in.getRows();
		int ny= in.getCols();
		int nz = in.getSlices();
		
		String pvType = pvParam.getValue();
		String outType = outParam.getValue();

		// main algorithm
		System.out.println("compute directional filter");
		float[][][] filter = new float[nx][ny][nz];
		ImageDataFloat resultData;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
			// check for zero-valued neighbors as well
			int dmax = 13;
			if (strictParam.getValue()) {
				filter[x][y][z] = minmaxplaneScore(image, x,y,z, dmax);
			} else {
				float best = 0.0f;
				for (int d=0;d<dmax;d++) {
					float val = planeScore(image, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				filter[x][y][z] = best;
			}
		}
		
		float[][][] result = new float[nx][ny][nz];
		if (outType.equals("intensity")) {
			System.out.println("output raw filter response");
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				if ( (pvType.equals("bright") && filter[x][y][z]<0) || (pvType.equals("dark") && filter[x][y][z]>0) ) {
					result[x][y][z] = 0.0f;
				} else if  (pvType.equals("dark")) {
					result[x][y][z] = -filter[x][y][z];
				} else {
					result[x][y][z] = filter[x][y][z];
				}
			}
		}
		else if (outType.equals("probability")) {
			System.out.println("output filter response probability");
			boolean[][][] mask = new boolean[nx][ny][nz];
			double mean = 0.0;
			double den = 0.0;
			double max = 0.0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				// keep only the proper sign
				if ( (pvType.equals("bright") && filter[x][y][z]<0) || (pvType.equals("dark") && filter[x][y][z]>0) ) {
					filter[x][y][z] = 0.0f;
				} else if  (pvType.equals("dark")) {
					filter[x][y][z] = -filter[x][y][z];
				}
				mask[x][y][z] = (filter[x][y][z]!=0);
				
				if (mask[x][y][z]) {
					mean += filter[x][y][z];
					max = Numerics.max(Numerics.abs(filter[x][y][z]),max);
					den++;
				}
			}
			mean /= den;
			double normg = 2.0/(mean*FastMath.PI);
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
				double pg = FastMath.exp(-filter[x][y][z]*filter[x][y][z]/(FastMath.PI*mean*mean));
				result[x][y][z] = (float)( (1.0-pg)/max/( (1.0-pg)/max+pg*normg));
				//double pe = FastMath.exp( -0.5*filter[x][y][z]/mean);
				//result[x][y][z] = (float)( (1.0-pe)/max/( (1.0-pe)/max+pe/mean));
			}
			// loop to stable results? no: may converge to zero too easily...
		}
		/*
		else if (outType.equals("probability") && (!strictParam.getValue())) {
			System.out.println("output filter response probability");
			float avg = 0.0f;
			float sig = 0.0f;
			float den = 0.0f;
			float max = 0.0f;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				avg += filter[x][y][z];
				max = Numerics.max(Numerics.abs(filter[x][y][z]), max);
				den++;
			}
			avg /= den;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				sig += (filter[x][y][z]-avg)*(filter[x][y][z]-avg);
			}
			sig /= (den-1.0);
			
			float normg = (float)FastMath.sqrt(2.0*FastMath.PI*sig);
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				// keep only the proper sign
				if ( (pvType.equals("bright") && filter[x][y][z]<0) || (pvType.equals("dark") && filter[x][y][z]>0) ) {
					filter[x][y][z] = 0.0f;
				} else if  (pvType.equals("dark")) {
					filter[x][y][z] = -filter[x][y][z];
				}
				float pg = (float)FastMath.exp( -0.5f*(filter[x][y][z]-avg)*(filter[x][y][z]-avg)/sig );
				result[x][y][z] = (1.0f-pg)/max/( (1.0f-pg)/max+pg/normg);
			}
			
		} 
		else if (outType.equals("probability") && (strictParam.getValue())) {
			System.out.println("output filter response probability");
			boolean[][][] mask = new boolean[nx][ny][nz];
			double mean = 0.0;
			double den = 0.0;
			double max = 0.0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				// keep only the proper sign
				if ( (pvType.equals("bright") && filter[x][y][z]<0) || (pvType.equals("dark") && filter[x][y][z]>0) ) {
					filter[x][y][z] = 0.0f;
				} else if  (pvType.equals("dark")) {
					filter[x][y][z] = -filter[x][y][z];
				}
				mask[x][y][z] = (filter[x][y][z]!=0);
				
				if (mask[x][y][z]) {
					mean += filter[x][y][z];
					max = Numerics.max(Numerics.abs(filter[x][y][z]),max);
					den++;
				}
			}
			mean /= den;
			double normg = 2.0/(mean*FastMath.PI);
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
				double pg = FastMath.exp(-filter[x][y][z]*filter[x][y][z]/(FastMath.PI*mean*mean));
				result[x][y][z] = (float)( (1.0-pg)/max/( (1.0-pg)/max+pg*normg));
				//double pe = FastMath.exp( -0.5*filter[x][y][z]/mean);
				//result[x][y][z] = (float)( (1.0-pe)/max/( (1.0-pe)/max+pe/mean));
			}
			// loop to stable results? no: may converge to zero too easily...
		}
		*/
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(in.getHeader());
		resultData.setName(in.getName()+"_rdg");
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

	float minmaxplaneScore(float[][][] image, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		for (int d=0;d<dmax;d++) {
			float val1 = 0.0f, val2 = 0.0f;
			if (d==0) {		
				val1=(image[x][y][z]		-image[x-1][y][z]
					 +image[x][y-1][z]		-image[x-1][y-1][z]
					 +image[x][y+1][z]		-image[x-1][y+1][z]
					 +image[x][y][z-1]		-image[x-1][y][z-1]
					 +image[x][y][z+1]		-image[x-1][y][z+1]
					 +image[x][y-1][z-1]	-image[x-1][y-1][z-1]
					 +image[x][y-1][z+1]	-image[x-1][y-1][z+1]
					 +image[x][y+1][z-1]	-image[x-1][y+1][z-1]
					 +image[x][y+1][z+1]	-image[x-1][y+1][z+1])/9.0f;
				val2=(image[x][y][z]		-image[x+1][y][z]
					 +image[x][y-1][z]		-image[x+1][y-1][z]
					 +image[x][y+1][z]		-image[x+1][y+1][z]
					 +image[x][y][z-1]		-image[x+1][y][z-1]
					 +image[x][y][z+1]		-image[x+1][y][z+1]
					 +image[x][y-1][z-1]	-image[x+1][y-1][z-1]
					 +image[x][y-1][z+1]	-image[x+1][y-1][z+1]
					 +image[x][y+1][z-1]	-image[x+1][y+1][z-1]
					 +image[x][y+1][z+1]	-image[x+1][y+1][z+1])/9.0f;
			} else if (d==1) {
				val1=(image[x][y][z]		-image[x][y-1][z]
					 +image[x-1][y][z]		-image[x-1][y-1][z]
					 +image[x+1][y][z]		-image[x+1][y-1][z]	
					 +image[x][y][z-1]		-image[x][y-1][z-1]	
					 +image[x][y][z+1]		-image[x][y-1][z+1]
					 +image[x-1][y][z-1]	-image[x-1][y-1][z-1]
					 +image[x-1][y][z+1]	-image[x-1][y-1][z+1]
					 +image[x+1][y][z-1]	-image[x+1][y-1][z-1]
					 +image[x+1][y][z+1]	-image[x+1][y-1][z+1])/9.0f;
				val2=(image[x][y][z]		-image[x][y+1][z]
					 +image[x-1][y][z]		-image[x-1][y+1][z]
					 +image[x+1][y][z]		-image[x+1][y+1][z]
					 +image[x][y][z-1]		-image[x][y+1][z-1]
					 +image[x][y][z+1]		-image[x][y+1][z+1]
					 +image[x-1][y][z-1]	-image[x-1][y+1][z-1]
					 +image[x-1][y][z+1]	-image[x-1][y+1][z+1]
					 +image[x+1][y][z-1]	-image[x+1][y+1][z-1]
					 +image[x+1][y][z+1]	-image[x+1][y+1][z+1])/9.0f;
			} else if (d==2) { 			
				val1=(image[x][y][z]		-image[x][y][z-1]
					 +image[x-1][y][z]		-image[x-1][y][z-1]
					 +image[x+1][y][z]		-image[x+1][y][z-1]
					 +image[x][y-1][z]		-image[x][y-1][z-1]
					 +image[x][y+1][z]		-image[x][y+1][z-1]	
					 +image[x-1][y-1][z]	-image[x-1][y-1][z-1]
					 +image[x-1][y+1][z]	-image[x-1][y+1][z-1]
					 +image[x+1][y-1][z]	-image[x+1][y-1][z-1]
					 +image[x+1][y+1][z]	-image[x+1][y+1][z-1])/9.0f;
				val2=(image[x][y][z]		-image[x][y][z+1]
					 +image[x-1][y][z]		-image[x-1][y][z+1]
					 +image[x+1][y][z]		-image[x+1][y][z+1]
					 +image[x][y-1][z]		-image[x][y-1][z+1]
					 +image[x][y+1][z]		-image[x][y+1][z+1]
					 +image[x-1][y-1][z]	-image[x-1][y-1][z+1]
					 +image[x-1][y+1][z]	-image[x-1][y+1][z+1]
					 +image[x+1][y-1][z]	-image[x+1][y-1][z+1]
					 +image[x+1][y+1][z]	-image[x+1][y+1][z+1])/9.0f;
			} else if (d==3) { 			
				val1=(image[x][y][z]		-image[x-1][y-1][z]
					 +image[x-1][y+1][z]	-image[x-2][y][z]
					 +image[x+1][y-1][z]	-image[x][y-2][z]
					 +image[x][y][z-1]		-image[x-1][y-1][z-1]
					 +image[x-1][y+1][z-1]	-image[x-2][y][z-1]
					 +image[x+1][y-1][z-1]	-image[x][y-2][z-1]
					 +image[x][y][z+1]		-image[x-1][y-1][z+1]
					 +image[x-1][y+1][z+1]	-image[x-2][y][z+1]
					 +image[x+1][y-1][z+1]	-image[x][y-2][z+1])/9.0f;
				val2=(image[x][y][z]		-image[x+1][y+1][z]
					 +image[x-1][y+1][z]	-image[x][y+2][z]
					 +image[x+1][y-1][z]	-image[x+2][y][z]
					 +image[x][y][z-1]		-image[x+1][y+1][z-1]
					 +image[x-1][y+1][z-1]	-image[x][y+2][z-1]
					 +image[x+1][y-1][z-1]	-image[x+2][y][z-1]
					 +image[x][y][z+1]		-image[x+1][y+1][z+1]
					 +image[x-1][y+1][z+1]	-image[x][y+2][z+1]
					 +image[x+1][y-1][z+1]	-image[x+2][y][z+1])/9.0f;
			} else if (d==4) { 			
				val1=(image[x][y][z]		-image[x][y-1][z-1]
					 +image[x][y+1][z-1]	-image[x][y][z-2]
					 +image[x][y-1][z+1]	-image[x][y-2][z]
					 +image[x-1][y][z]		-image[x-1][y-1][z-1]
					 +image[x-1][y+1][z-1]-image[x-1][y][z-2]
					 +image[x-1][y-1][z+1]-image[x-1][y-2][z]
					 +image[x+1][y][z]		-image[x+1][y-1][z-1]
					 +image[x+1][y+1][z-1]-image[x+1][y][z-2]
					 +image[x+1][y-1][z+1]-image[x+1][y-2][z])/9.0f;
				val2=(image[x][y][z]		-image[x][y+1][z+1]
					 +image[x][y+1][z-1]	-image[x][y+2][z]
					 +image[x][y-1][z+1]	-image[x][y][z+2]
					 +image[x-1][y][z]		-image[x-1][y+1][z+1]
					 +image[x-1][y+1][z-1]	-image[x-1][y+2][z]
					 +image[x-1][y-1][z+1]	-image[x-1][y][z+2]
					 +image[x+1][y][z]		-image[x+1][y+1][z+1]
					 +image[x+1][y+1][z-1]	-image[x+1][y+2][z]
					 +image[x+1][y-1][z+1]	-image[x+1][y][z+2])/9.0f;
			} else if (d==5) { 			
				val1=(image[x][y][z]		-image[x-1][y][z-1]
					 +image[x+1][y][z-1]	-image[x][y][z-2]
					 +image[x-1][y][z+1]	-image[x-2][y][z]
					 +image[x][y-1][z]		-image[x-1][y-1][z-1]
					 +image[x+1][y-1][z-1]-image[x][y-1][z-2]
					 +image[x-1][y-1][z+1]-image[x-2][y-1][z]
					 +image[x][y+1][z]		-image[x-1][y+1][z-1]
					 +image[x+1][y+1][z-1]-image[x][y+1][z-2]
					 +image[x-1][y+1][z+1]-image[x-2][y+1][z])/9.0f;
				val2=(image[x][y][z]		-image[x+1][y][z+1]
					 +image[x+1][y][z-1]	-image[x+2][y][z]
					 +image[x-1][y][z+1]	-image[x][y][z+2]
					 +image[x][y-1][z]		-image[x+1][y-1][z+1]
					 +image[x+1][y-1][z-1]	-image[x+2][y-1][z]
					 +image[x-1][y-1][z+1]	-image[x][y-1][z+2]
					 +image[x][y+1][z]		-image[x+1][y+1][z+1]
					 +image[x+1][y+1][z-1]	-image[x+2][y+1][z]
					 +image[x-1][y+1][z+1]	-image[x][y+1][z+2])/9.0f;
			} else if (d==6) { 			
				val1=(image[x][y][z]		-image[x-1][y+1][z]
					 +image[x-1][y-1][z]	-image[x-2][y][z]
					 +image[x+1][y+1][z]	-image[x][y-2][z]
					 +image[x][y][z-1]		-image[x-1][y+1][z-1]
					 +image[x-1][y-1][z-1]-image[x-2][y][z-1]
					 +image[x+1][y+1][z-1]-image[x][y-2][z-1]
					 +image[x][y][z+1]		-image[x-1][y+1][z+1]
					 +image[x-1][y-1][z+1]-image[x-2][y][z+1]
					 +image[x+1][y+1][z+1]-image[x][y-2][z+1])/9.0f;
				val2=(image[x][y][z]		-image[x+1][y-1][z]
					 +image[x-1][y-1][z]	-image[x][y-2][z]
					 +image[x+1][y+1][z]	-image[x+2][y][z]
					 +image[x][y][z-1]		-image[x+1][y-1][z-1]
					 +image[x-1][y-1][z-1]	-image[x][y-2][z-1]
					 +image[x+1][y+1][z-1]	-image[x+2][y][z-1]
					 +image[x][y][z+1]		-image[x+1][y-1][z+1]
					 +image[x-1][y-1][z+1]	-image[x][y-2][z+1]
					 +image[x+1][y+1][z+1]	-image[x+2][y][z+1])/9.0f;
			} else if (d==7) { 			
				val1=(image[x][y][z]		-image[x][y-1][z+1]
					 +image[x][y-1][z-1]	-image[x][y-2][z]
					 +image[x][y+1][z+1]	-image[x][y][z+2]
					 +image[x-1][y][z]		-image[x-1][y-1][z+1]
					 +image[x-1][y-1][z-1]	-image[x-1][y-2][z]
					 +image[x-1][y+1][z+1]	-image[x-1][y][z+2]
					 +image[x+1][y][z]		-image[x+1][y-1][z+1]
					 +image[x+1][y-1][z-1]	-image[x+1][y-2][z]
					 +image[x+1][y+1][z+1]	-image[x+1][y][z+2])/9.0f;
				val2=(image[x][y][z]		-image[x][y+1][z-1]
					 +image[x][y-1][z-1]	-image[x][y][z-2]
					 +image[x][y+1][z+1]	-image[x][y+2][z]
					 +image[x-1][y][z]		-image[x-1][y+1][z-1]
					 +image[x-1][y-1][z-1]	-image[x-1][y][z-2]
					 +image[x-1][y+1][z+1]	-image[x-1][y+2][z]
					 +image[x+1][y][z]		-image[x+1][y+1][z-1]
					 +image[x+1][y-1][z-1]	-image[x+1][y][z-2]
					 +image[x+1][y+1][z+1]	-image[x+1][y+2][z])/9.0f;
			} else if (d==8) { 			
				val1=(image[x][y][z]		-image[x-1][y][z+1]
					 +image[x-1][y][z-1]	-image[x-2][y][z]
					 +image[x+1][y][z+1]	-image[x][y][z+2]
					 +image[x][y-1][z]		-image[x-1][y-1][z+1]
					 +image[x-1][y-1][z-1]-image[x-2][y-1][z]
					 +image[x+1][y-1][z+1]-image[x][y-1][z+2]
					 +image[x][y+1][z]		-image[x-1][y+1][z+1]
					 +image[x-1][y+1][z-1]-image[x-2][y+1][z]
					 +image[x+1][y+1][z+1]-image[x][y+1][z+2])/9.0f;
				val2=(image[x][y][z]		-image[x+1][y][z-1]
					 +image[x-1][y][z-1]	-image[x][y][z-2]
					 +image[x+1][y][z+1]	-image[x+2][y][z]
					 +image[x][y-1][z]		-image[x+1][y-1][z-1]
					 +image[x-1][y-1][z-1]	-image[x][y-1][z-2]
					 +image[x+1][y-1][z+1]	-image[x+2][y-1][z]
					 +image[x][y+1][z]		-image[x+1][y+1][z-1]
					 +image[x-1][y+1][z-1]	-image[x][y+1][z-2]
					 +image[x+1][y+1][z+1]	-image[x+2][y+1][z])/9.0f;
			} else if (d==9) { 			
				val1=(image[x][y][z]		-image[x-1][y-1][z-1]
					 +image[x-1][y-1][z+1]	-image[x-2][y-2][z]
					 +image[x-1][y+1][z-1]	-image[x-2][y][z-2]
					 +image[x+1][y-1][z-1]	-image[x][y-2][z-2]
					 +image[x+1][y+1][z-1]	-image[x][y][z-2]
					 +image[x-1][y+1][z+1]	-image[x-2][y][z]
					 +image[x+1][y-1][z+1]	-image[x][y-2][z])/7.0f;
				val2=(image[x][y][z]		-image[x+1][y+1][z+1]
					 +image[x-1][y-1][z+1]	-image[x][y][z+2]
					 +image[x-1][y+1][z-1]	-image[x][y+2][z]
					 +image[x+1][y-1][z-1]	-image[x+2][y][z]
					 +image[x+1][y+1][z-1]	-image[x+2][y+2][z]
					 +image[x-1][y+1][z+1]	-image[x][y+2][z+2]
					 +image[x+1][y-1][z+1]	-image[x+2][y][z+2])/7.0f;
			} else if (d==10) { 			
				val1=(image[x][y][z]		-image[x+1][y-1][z-1]
					 +image[x+1][y-1][z+1]	-image[x+2][y-2][z]
					 +image[x+1][y+1][z-1]	-image[x+2][y][z-2]
					 +image[x-1][y-1][z-1]	-image[x][y-2][z-2]
					 +image[x-1][y+1][z-1]	-image[x][y][z-2]
					 +image[x-1][y-1][z+1]	-image[x][y-2][z]
					 +image[x+1][y+1][z+1]	-image[x+2][y][z])/7.0f;
				val2=(image[x][y][z]		-image[x-1][y+1][z+1]
					 +image[x+1][y-1][z+1]	-image[x][y][z+2]
					 +image[x+1][y+1][z-1]	-image[x][y+2][z]
					 +image[x-1][y-1][z-1]	-image[x-2][y][z]
					 +image[x-1][y+1][z-1]	-image[x-2][y+2][z]
					 +image[x-1][y-1][z+1]	-image[x-2][y][z+2]
					 +image[x+1][y+1][z+1]	-image[x][y+2][z+2])/7.0f;
			} else if (d==11) { 			
				val1=(image[x][y][z]		-image[x-1][y+1][z-1]
					 +image[x-1][y+1][z+1]	-image[x-2][y+2][z]
					 +image[x+1][y+1][z-1]	-image[x][y+2][z-2]
					 +image[x-1][y-1][z-1]	-image[x-2][y][z-2]
					 +image[x+1][y-1][z-1]	-image[x][y][z-2]
					 +image[x-1][y-1][z+1]	-image[x-2][y][z]
					 +image[x+1][y+1][z+1]	-image[x][y+2][z])/7.0f;
				val2=(image[x][y][z]		-image[x+1][y-1][z+1]
					 +image[x-1][y+1][z+1]	-image[x][y][z+2]
					 +image[x+1][y+1][z-1]	-image[x+2][y][z]
					 +image[x-1][y-1][z-1]	-image[x][y-2][z]
					 +image[x+1][y-1][z-1]	-image[x+2][y-2][z]
					 +image[x-1][y-1][z+1]	-image[x][y-2][z+2]
					 +image[x+1][y+1][z+1]	-image[x+2][y][z+2])/7.0f;
			} else if (d==12) { 			
				val1=(image[x][y][z]		-image[x-1][y-1][z+1]
					 +image[x-1][y+1][z+1]	-image[x-2][y][z+2]
					 +image[x+1][y-1][z+1]	-image[x][y-2][z+2]
					 +image[x-1][y-1][z-1]	-image[x-2][y-2][z]
					 +image[x+1][y-1][z-1]	-image[x][y-2][z]
					 +image[x-1][y+1][z-1]	-image[x-2][y][z]
					 +image[x+1][y+1][z+1]	-image[x][y][z+2])/7.0f;
				val2=(image[x][y][z]		-image[x+1][y+1][z-1]
					 +image[x-1][y+1][z+1]	-image[x][y+2][z]
					 +image[x+1][y-1][z+1]	-image[x+2][y][z]
					 +image[x-1][y-1][z-1]	-image[x][y][z-2]
					 +image[x+1][y-1][z-1]	-image[x+2][y][z-2]
					 +image[x-1][y+1][z-1]	-image[x][y+2][z-2]
					 +image[x+1][y+1][z+1]	-image[x+2][y+2][z])/7.0f;
			}
			// find the strongest gradient direction, then estimate the corresponding filter response
			if (val1*val1+val2*val2>maxgrad) {
				maxgrad = val1*val1+val2*val2;
				if (val1*val1<val2*val2) minval = val1;
				else minval = val2;
				sign = val1*val2;
			}
		}
		if (sign>0) return minval;
		else return 0.0f;
	}
}
