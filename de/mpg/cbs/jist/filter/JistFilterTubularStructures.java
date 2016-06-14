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
import de.mpg.cbs.structures.Histogram;

/*
 * @author Pierre-Louis Bazin
 */
public class JistFilterTubularStructures extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	pvtypes = {"bright","dark","both"};
	private		String		pvtype = "bright";
	private		static final String[]	outtypes = {"intensity","probability","squared_proba"};
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
		pvParam.setDescription("Outputs the raw intensity values or a probability score for the structures.");
		pvParam.setValue(outtype);
		inputParams.add(strictParam = new ParamBoolean("use strict min/max filter", true));
		strictParam.setDescription("Using either a broader response filter or a very strict min/max filter");
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Filter");
		inputParams.setLabel("Tubular Structures");
		inputParams.setName("TubularStructures");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Filters an image for 3D tubular structures.");
		
		info.setVersion("3.0.4");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Tubular Structure Image",VoxelType.FLOAT));
		outputParams.setName("tubular_image");
		outputParams.setLabel("Tubular Image");
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
		float[][][] filter = new float[nx][ny][nz];
		float[][][] planescore = new float[nx][ny][nz];
		byte[][][] planedir = new byte[nx][ny][nz];
		ImageDataFloat resultData;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			//System.out.println("("+x+", "+y+", "+z+")");
			if (!zeroNeighbor(image,x,y,z,2)) {
				// check for zero-valued neighbors as well
				int dmax = 13;
				if (strictParam.getValue()) {
					minmaxplaneScore(image, planescore, planedir, x,y,z, 13);
					// sign issue: remove all that have different sign, keep global sign
					float linescore = minmaxlineScore(image, planescore, planedir, x,y,z, 4);
					if (planescore[x][y][z]*linescore>0) {
						filter[x][y][z] = Numerics.sign(linescore)*Numerics.sqrt(planescore[x][y][z]*linescore);
					} else {
						filter[x][y][z] = 0.0f;
					}
				} else {
					float best = 0.0f;
					for (int d=0;d<dmax;d++) {
						float val = tubularScore(image, x,y,z,d);
						if (val*val>best*best) best = val;
					}
					filter[x][y][z] = best;
				}
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
			double mean = 0.0f;
			double den = 0.0f;
			double max = 0.0f;
			boolean[][][] mask = new boolean[nx][ny][nz];
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (!zeroNeighbor(image,x,y,z,2)) {
				// keep only the proper sign
				if ( (pvType.equals("bright") && filter[x][y][z]<0) || (pvType.equals("dark") && filter[x][y][z]>0) ) {
					filter[x][y][z] = 0.0f;
				} else if  (pvType.equals("dark")) {
					filter[x][y][z] = -filter[x][y][z];
				}
				mask[x][y][z] = (filter[x][y][z]!=0);
				if (mask[x][y][z]) {
					// fit exp only to non-zero data: much nicer!
					mean += filter[x][y][z];
					den++;
					max = Numerics.max(Numerics.abs(filter[x][y][z]), max);
				}
			}
			mean /= den;
			
			// exponential fit better than gaussian fit
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
				double pe = FastMath.exp( -0.5*filter[x][y][z]/mean)/mean;
				result[x][y][z] = (float)(1.0/max/( 1.0/max+pe));
			}
		}
		
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(in.getHeader());
		resultData.setName(in.getName()+"_tub");
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
	public static float tubularScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z])
						-image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z]
						-image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z]
						-image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1]
						-image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1]
						-image[x][y-1][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1]
						-image[x][y-1][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1]
						-image[x][y+1][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1]
						-image[x][y+1][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1])/24.0f;		
		} else if (d==1) {
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z])
						-image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z]
						-image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z]
						-image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1]
						-image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1]
						-image[x-1][y][z-1]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1]
						-image[x-1][y][z+1]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1]
						-image[x+1][y][z-1]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1]
						-image[x+1][y][z+1]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1])/24.0f;
		} else if (d==2) {
			val =(8.0f*(image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1])
						-image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1]
						-image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1]
						-image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1]
						-image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1]
						-image[x-1][y-1][z]		-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1]
						-image[x-1][y+1][z]		-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1]
						-image[x+1][y-1][z]		-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1]
						-image[x+1][y+1][z]		-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1])/24.0f;
		} else if (d==3) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z])
						-image[x-1][y+1][z]		-image[x-2][y][z]		-image[x][y+2][z]
						-image[x+1][y-1][z]		-image[x][y-2][z]		-image[x+2][y][z]
						-image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-1]		-image[x][y+2][z-1]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
						-image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1]
						-image[x-1][y+1][z+1]	-image[x-2][y][z+1]		-image[x][y+2][z+1]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
		} else if (d==4) {		
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1])
						-image[x][y+1][z-1]		-image[x][y][z-2]		-image[x][y+2][z]
						-image[x][y-1][z+1]		-image[x][y-2][z]		-image[x][y][z+2]
						-image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1]
						-image[x-1][y+1][z-1]	-image[x-1][y][z-2]		-image[x-1][y-2][z]
						-image[x-1][y-1][z+1]	-image[x-1][y-2][z]		-image[x-1][y][z+2]
						-image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x+1][y][z-2]		-image[x+1][y][z+2]
						-image[x+1][y-1][z+1]	-image[x+1][y-2][z]		-image[x+1][y][z+2])/24.0f;
		} else if (d==5) {	
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1])
						-image[x+1][y][z-1]		-image[x][y][z-2]		-image[x+2][y][z]
						-image[x-1][y][z+1]		-image[x-2][y][z]		-image[x][y][z+2]
						-image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1]
						-image[x+1][y-1][z-1]	-image[x][y-1][z-2]		-image[x+2][y-1][z]
						-image[x-1][y-1][z+1]	-image[x-2][y-1][z]		-image[x][y-1][z+2]
						-image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x][y+1][z-2]		-image[x+2][y+1][z]
						-image[x-1][y+1][z+1]	-image[x-2][y+1][z]		-image[x][y+1][z+2])/24.0f;
		} else if (d==6) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z])
						-image[x-1][y-1][z]		-image[x-2][y][z]		-image[x][y-2][z]
						-image[x+1][y+1][z]		-image[x][y-2][z]		-image[x+2][y][z]
						-image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y][z-1]		-image[x][y-2][z-1]
						-image[x+1][y+1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
						-image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1]
						-image[x-1][y-1][z+1]	-image[x-2][y][z+1]		-image[x][y-2][z+1]
						-image[x+1][y+1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
		} else if (d==7) {
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1])
						-image[x][y-1][z-1]		-image[x][y-2][z]		-image[x][y][z-2]
						-image[x][y+1][z+1]		-image[x][y][z+2]		-image[x][y+2][z]
						-image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1]
						-image[x-1][y-1][z-1]	-image[x-1][y-2][z]		-image[x-1][y][z-2]
						-image[x-1][y+1][z+1]	-image[x-1][y][z+2]		-image[x-1][y+2][z]
						-image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1]
						-image[x+1][y-1][z-1]	-image[x+1][y-2][z]		-image[x+1][y][z-2]
						-image[x+1][y+1][z+1]	-image[x+1][y][z+2]		-image[x+1][y+2][z])/24.0f;
		} else if (d==8) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1])
						-image[x-1][y][z-1]		-image[x-2][y][z]		-image[x][y][z-2]
						-image[x+1][y][z+1]		-image[x][y][z+2]		-image[x+2][y][z]
						-image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y-1][z]		-image[x][y-1][z-2]
						-image[x+1][y-1][z+1]	-image[x][y-1][z+2]		-image[x+2][y-1][z]
						-image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y+1][z]		-image[x][y+1][z-2]
						-image[x+1][y+1][z+1]	-image[x][y+1][z+2]		-image[x+2][y+1][z])/24.0f;
		} else if (d==9) {
			val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1])
						-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z]
						-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z]
						-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2]
						-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2])/18.0f;
		} else if (d==10) {
			val =(6.0f*(image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1])
						-image[x+1][y-1][z+1]	-image[x+2][y-2][z]		-image[x][y][z+2]
						-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z]
						-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z]
						-image[x-1][y+1][z-1]	-image[x][y][z-2]		-image[x-2][y+2][z]
						-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2]
						-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2])/18.0f;
		} else if (d==11) {				
			val =(6.0f*(image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1])
						 -image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2]
						 -image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z]
						 -image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z]
						 -image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z]
						 -image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2]
						 -image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2])/18.0f;
		} else if (d==12) {				
			val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1])
						-image[x-1][y+1][z+1]	-image[x-2][y][z-2]		-image[x][y+2][z]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+2]		-image[x+2][y][z]
						-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2]
						-image[x+1][y-1][z-1]	-image[x][y-2][z]		-image[x+2][y][z-2]
						-image[x-1][y+1][z-1]	-image[x-2][y][z]		-image[x][y+2][z-2]
						-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z])/18.0f;
		}
		return 0.5f*val;
	}

	/*
	public static float minmaxtubularScore(float[][][] image, int x, int y, int z, int dmax) {
		// same as ridge filter:
		// find region with strongest tubular difference (almost as before, with abs values)
		// use the minimum difference around the tube
		float maxgrad = 0.0f;
		float minval = 0.0f;
		float sign = 0.0f;
		for (int d=0;d<dmax;d++) {
			float[] val = new float[8];
			for (int n=0;n<8;n++) val[n] = 0.0f;
			if (d==0) {
				val[0]=  image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z];
				val[1]=  image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z];
				val[2]=  image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1];
				val[3]=	 image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1];
				val[4]=	 image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y-1][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1];
				val[5]=	 image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y+1][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1];
				val[6]=	 image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y-1][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1];
				val[7]=	 image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z]
						-image[x][y+1][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1];
			} else if (d==1) {
				val[0]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z];
				val[1]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z];
				val[2]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1];
				val[3]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1];
				val[4]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x-1][y][z-1]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1];
				val[5]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x+1][y][z+1]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1];
				val[6]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x-1][y][z+1]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1];
				val[7]=	 image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z]
						-image[x+1][y][z-1]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1];
			} else if (d==2) {
				val[0]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1];
				val[1]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1];
				val[2]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1];
				val[3]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1];
				val[4]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x-1][y-1][z]		-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1];
				val[5]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x+1][y+1][z]		-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1];
				val[6]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x-1][y+1][z]		-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1];
				val[7]=	 image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1]
						-image[x+1][y-1][z]		-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1];
			} else if (d==3) {
				val[0]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x-1][y+1][z]		-image[x-2][y][z]		-image[x][y+2][z];
				val[1]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x+1][y-1][z]		-image[x][y-2][z]		-image[x+2][y][z];
				val[2]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1];
				val[3]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1];
				val[4]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-1]		-image[x][y+2][z-1];
				val[5]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1];
				val[6]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1];
				val[7]=	 image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z]
						-image[x-1][y+1][z+1]	-image[x-2][y][z+1]		-image[x][y+2][z+1];
			} else if (d==4) {		
				val[0]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x][y+1][z-1]		-image[x][y][z-2]		-image[x][y+2][z];
				val[1]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x][y-1][z+1]		-image[x][y-2][z]		-image[x][y][z+2];
				val[2]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1];
				val[3]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1];
				val[4]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x-1][y+1][z-1]	-image[x-1][y][z-2]		-image[x-1][y-2][z];
				val[5]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x+1][y-1][z+1]	-image[x+1][y-2][z]		-image[x+1][y][z+2];
				val[6]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x-1][y-1][z+1]	-image[x-1][y-2][z]		-image[x-1][y][z+2];
				val[7]=	 image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x+1][y][z-2]		-image[x+1][y][z+2];
			} else if (d==5) {	
				val[0]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x+1][y][z-1]		-image[x][y][z-2]		-image[x+2][y][z];
				val[1]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x-1][y][z+1]		-image[x-2][y][z]		-image[x][y][z+2];
				val[2]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1];
				val[3]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1];
				val[4]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x+1][y-1][z-1]	-image[x][y-1][z-2]		-image[x+2][y-1][z];
				val[5]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x-1][y+1][z+1]	-image[x-2][y+1][z]		-image[x][y+1][z+2];
				val[6]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x-1][y-1][z+1]	-image[x-2][y-1][z]		-image[x][y-1][z+2];
				val[7]=	 image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1]
						-image[x+1][y+1][z-1]	-image[x][y+1][z-2]		-image[x+2][y+1][z];
			} else if (d==6) {
				val[0]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x-1][y-1][z]		-image[x-2][y][z]		-image[x][y-2][z];
				val[1]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x+1][y+1][z]		-image[x][y-2][z]		-image[x+2][y][z];
				val[2]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1];
				val[3]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1];
				val[4]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x-1][y-1][z-1]	-image[x-2][y][z-1]		-image[x][y-2][z-1];
				val[5]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x+1][y+1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1];
				val[6]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x+1][y+1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1];
				val[7]=	 image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z]
						-image[x-1][y-1][z+1]	-image[x-2][y][z+1]		-image[x][y-2][z+1];
			} else if (d==7) {
				val[0]=  image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x][y-1][z-1]		-image[x][y-2][z]		-image[x][y][z-2];
				val[1]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x][y+1][z+1]		-image[x][y][z+2]		-image[x][y+2][z];
				val[2]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1];
				val[3]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1];
				val[4]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x-1][y-1][z-1]	-image[x-1][y-2][z]		-image[x-1][y][z-2];
				val[5]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x+1][y+1][z+1]	-image[x+1][y][z+2]		-image[x+1][y+2][z];
				val[6]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x-1][y+1][z+1]	-image[x-1][y][z+2]		-image[x-1][y+2][z];
				val[7]=	 image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1]
						-image[x+1][y-1][z-1]	-image[x+1][y-2][z]		-image[x+1][y][z-2];
			} else if (d==8) {
				val[0]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x-1][y][z-1]		-image[x-2][y][z]		-image[x][y][z-2];
				val[1]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x+1][y][z+1]		-image[x][y][z+2]		-image[x+2][y][z];
				val[2]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1];
				val[3]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1];
				val[4]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y-1][z]		-image[x][y-1][z-2];
				val[5]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x+1][y+1][z+1]	-image[x][y+1][z+2]		-image[x+2][y+1][z];
				val[6]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x+1][y-1][z+1]	-image[x][y-1][z+2]		-image[x+2][y-1][z];
				val[7]=	 image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y+1][z]		-image[x][y+1][z-2];
			} else if (d==9) {
				val[0]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2];
				val[1]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z];
				val[2]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z];
				val[3]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2];
				val[4]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z];
				val[5]=	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
						-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2];
			} else if (d==10) {
				val[0]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x+1][y-1][z+1]	-image[x+2][y-2][z]		-image[x][y][z+2];
				val[1]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x-1][y+1][z-1]	-image[x][y][z-2]		-image[x-2][y+2][z];
				val[2]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z];
				val[3]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2];
				val[4]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z];
				val[5]=	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
						-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2];
			} else if (d==11) {				
				val[0]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2];
				val[1]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z];
				val[2]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z];
				val[3]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2];
				val[4]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z];
				val[5]=	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
						-image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2];
			} else if (d==12) {				
				val[0]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x-1][y+1][z+1]	-image[x-2][y][z-2]		-image[x][y+2][z];
				val[1]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x+1][y-1][z-1]	-image[x][y-2][z]		-image[x+2][y][z-2];
				val[2]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+2]		-image[x+2][y][z];
				val[3]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y][z]		-image[x][y+2][z-2];
				val[4]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2];
				val[5]=	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
						-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z];
			}
			
			// select the direction of strongest variation from center
			float grad = (Numerics.abs(val[0])+Numerics.abs(val[1])+Numerics.abs(val[2])+Numerics.abs(val[3])
						 +Numerics.abs(val[4])+Numerics.abs(val[5])+Numerics.abs(val[6])+Numerics.abs(val[7]) );
			/*				
			float grad = (Numerics.abs(val1+val2)+Numerics.abs(val3+val4)
						 +Numerics.abs(val5+val6)+Numerics.abs(val7+val8) );
			//float grad = (Numerics.abs(val[0]+val[1]+val[2]+val[3]+val[4]+val[5]+val[6]+val[7]) );
			float grad = (Numerics.abs(val[0]+val[1]+val[2]+val[3]) );
			*/	
			/*
			if (d>8) grad /= 6.0f;
			else grad /= 8.0f;
			
			if (grad>maxgrad) {
				maxgrad = grad;
				
				// check opposite signs for gradient-like changes (strictly negative)
				sign = Numerics.min(val[0]*val[1], val[2]*val[3], val[4]*val[5], val[6]*val[7]);
				
				// choose the smallest
				if (d>8) minval = Numerics.minmag(val, 6);
				else minval = Numerics.minmag(val, 8);
				/*
				// combine opposites?
				val[0] = val[0]+val[1];
				val[1] = val[2]+val[3];
				val[2] = val[4]+val[5];
				val[3] = val[6]+val[7];
				
				// choose the smallest
				if (d>8) minval = Numerics.minmag(val, 3);
				else minval = Numerics.minmag(val, 4);
				
				/*
				//if (d>8) minval = Numerics.minmag(val1+val2,val3+val4,val5+val6);
				//else minval = Numerics.minmag(val1+val2,val3+val4,val5+val6,val7+val8);
				//if (d>8) minval = Numerics.minmag(val, 6);
				//else minval = Numerics.minmag(val, 8);
				
				// use the  smallest values?
				if (d>8) {
					Numerics.sortmag(val, 3);
					//minval = val[0]+val[1]+val[2];
					//sign = minval*(val[3]+val[4]+val[5]);
					//minval = val[0]+val[1];
					//sign = minval*(val[4]+val[5]);
					minval = val[0];
					sign = val[0]*val[2];
				} else {
					Numerics.sortmag(val, 4);
					//minval = val[0]+val[1]+val[2]+val[3];
					//sign = minval*(val[4]+val[5]+val[6]+val[7]);
					//minval = val[0]+val[1];
					//sign = minval*(val[6]+val[7]);
					minval = val[0];
					sign = val[0]*val[3];		
				}
				*/
				//all the same sign? if so, any pair is >=0
				/*if (val1+val2>=0 && val3+val4>=0 && val5+val6>=0 && val7+val8>=0) sign = +1.0f;
				else if (val1+val2<=0 && val3+val4<=0 && val5+val6<=0 && val7+val8<=0) sign = +1.0f;
				
				if (val1>=0 && val2>=0 && val3>=0 && val4>=0 && val5>=0 && val6>=0 && val7>=0 && val8>=0) sign = +1.0f;
				else if (val1<=0 && val2<=0 && val3<=0 && val4<=0 && val5<=0 && val6<=0 && val7<=0 && val8<=0) sign = +1.0f;
				else sign = -1.0f;
				*/
				/*
				//sign = 1.0f;
			}
		}
		if (sign<0) return 0.0f;
		else return minval;
	}
	*/
	void minmaxplaneScore(float[][][] image, float[][][] plane, byte[][][] dir, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		byte direction = -1;
		for (byte d=0;d<dmax;d++) {
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
				direction = d;
			}
		}
		if (sign>0) {
			plane[x][y][z] = minval;
			dir[x][y][z] = direction;
		} else {
			plane[x][y][z] = 0.0f;
			dir[x][y][z] = -1;
		}
		return;
	}
	float minmaxlineScore(float[][][] image, float[][][] plane, byte[][][] dir, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		byte direction = -1;
		
		float val1 = 0.0f, val2 = 0.0f;
		if (dir[x][y][z]==0) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x][y][z-1]	-image[x][y-1][z-1]	-image[x][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x][y][z+1]	-image[x][y-1][z+1]	-image[x][y+1][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x][y-1][z]	-image[x][y-1][z-1]	-image[x][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x][y+1][z+1]	-image[x][y+1][z-1]	-image[x][y+1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z-1]	+image[x][y+1][z+1]
							-image[x][y-1][z+1]	-image[x][y-2][z]	-image[x][y][z+2];
					val2 = 	 image[x][y][z]		+image[x][y-1][z-1]	+image[x][y+1][z+1]
							-image[x][y+1][z-1]	-image[x][y][z-2]	-image[x][y+2][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z+1]	+image[x][y+1][z-1]
							-image[x][y-1][z-1]	-image[x][y-2][z]	-image[x][y][z-2];
					val2 = 	 image[x][y][z]		+image[x][y+1][z-1]	+image[x][y-1][z+1]
							-image[x][y+1][z+1]	-image[x][y+2][z]	-image[x][y][z+2];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==1) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x-1][y][z]	-image[x-1][y][z-1]	-image[x-1][y][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x+1][y][z]	-image[x+1][y][z-1]	-image[x+1][y][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y][z-1]	-image[x-1][y][z-1]	-image[x+1][y][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y][z+1]	-image[x-1][y][z+1]	-image[x+1][y][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z-1]	+image[x+1][y][z+1]
							-image[x-1][y][z+1]	-image[x-2][y][z]	-image[x][y][z+2];
					val2 = 	 image[x][y][z]		+image[x-1][y][z-1]	+image[x+1][y][z+1]
							-image[x+1][y][z-1]	-image[x][y][z-2]	-image[x+2][y][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z+1]	+image[x+1][y][z-1]
							-image[x-1][y][z-1]	-image[x-2][y][z]	-image[x][y][z-2];
					val2 = 	 image[x][y][z]		+image[x-1][y][z+1]	+image[x+1][y][z-1]
							-image[x+1][y][z+1]	-image[x][y][z+2]	-image[x+2][y][z];
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==2) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y-1][z]	-image[x-1][y-1][z]	-image[x+1][y-1][z];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y+1][z]	-image[x-1][y+1][z]	-image[x+1][y+1][z];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x-1][y][z]	-image[x-1][y-1][z]	-image[x-1][y+1][z];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x+1][y][z]	-image[x+1][y-1][z]	-image[x+1][y+1][z];
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x-1][y-1][z]	+image[x+1][y+1][z]
							-image[x-1][y+1][z]	-image[x-2][y][z]	-image[x][y+2][z];
					val2 = 	 image[x][y][z]		+image[x-1][y-1][z]	+image[x+1][y+1][z]
							-image[x+1][y-1][z]	-image[x][y-2][z]	-image[x+2][y][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x-1][y+1][z]	+image[x+1][y-1][z]
							-image[x-1][y-1][z]	-image[x-2][y][z]	-image[x][y-2][z];
					val2 = 	 image[x][y][z]		+image[x-1][y+1][z]	+image[x+1][y-1][z]
							-image[x+1][y+1][z]	-image[x][y+2][z]	-image[x+2][y][z];
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==3) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y+1][z]		+image[x+1][y-1][z]
							-image[x][y][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y+1][z]		+image[x+1][y-1][z]
							-image[x][y][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x-1][y+1][z]	-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x+1][y-1][z]	-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
							-image[x-1][y+1][z-1]	-image[x-2][y+2][z]		-image[x][y][z-2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
							-image[x+1][y-1][z+1]	-image[x][y][z+2]		-image[x+2][y-2][z];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==4) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y+1][z-1]		+image[x][y-1][z+1]
							-image[x-1][y][z]	-image[x-1][y+1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y+1][z-1]		+image[x][y-1][z+1]
							-image[x+1][y][z]	-image[x+1][y+1][z-1]	-image[x+1][y-1][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y+1][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y-1][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x+1][y-1][z+1]	-image[x+2][y][z]		-image[x][y-2][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x-1][y+1][z-1]	-image[x][y+2][z-2]		-image[x-2][y][z];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==5) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x+1][y][z-1]		+image[x-1][y][z+1]
							-image[x][y-1][z]	-image[x+1][y-1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x+1][y][z-1]		+image[x-1][y][z+1]
							-image[x][y+1][z]	-image[x+1][y+1][z-1]	-image[x-1][y+1][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x+1][y][z-1]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x-1][y][z+1]	-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x-1][y+1][z+1]	-image[x][y+2][z]		-image[x-2][y][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x+1][y-1][z-1]	-image[x+2][y][z-2]		-image[x][y-2][z];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==6) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y-1][z]		+image[x+1][y+1][z]
							-image[x][y][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y-1][z]		+image[x+1][y+1][z]
							-image[x][y][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1];
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x-1][y-1][z]	-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x+1][y+1][z]	-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z];
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
							-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
							-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==7) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]		+image[x][y-1][z-1]		+image[x][y+1][z+1]
							-image[x-1][y][z]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1];
					val2 =	 image[x][y][z]		+image[x][y-1][z-1]		+image[x][y+1][z+1]
							-image[x+1][y][z]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1];
				} else if (d==1) {
					val1 =	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y-1][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1];
					val2 =	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y+1][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1];
				} else if (d==2) {
					val1 =   image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z];
					val2 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2];
				} else if (d==3) {
					val1 =   image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z];
					val2 =	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2];
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==8) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]		+image[x-1][y][z-1]		+image[x+1][y][z+1]
							-image[x][y-1][z]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1];
					val2 =	 image[x][y][z]		+image[x-1][y][z-1]		+image[x+1][y][z+1]
							-image[x][y+1][z]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1];
				} else if (d==1) {
					val1 =	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x-1][y][z-1]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1];
					val2 =	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x+1][y][z+1]	-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1];
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2];
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2];
				}					
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==9) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x-1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x+1][y-1][z-1]);
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
											  -1.5f*(image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]);
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x+1][y-1][z+1]);
					val2 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x-1][y+1][z-1]);
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==10) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x-1][y-1][z-1]);
				} else 	if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x-1][y+1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x+1][y-1][z+1]);
				} else 	if (d==2) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x+1][y+1][z-1]);
				}		
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}			
		} else if (dir[x][y][z]==11) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x-1][y-1][z-1]);
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x-1][y+1][z+1]);
				} else 	if (d==2) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]);
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		} else if (dir[x][y][z]==12) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x-1][y-1][z-1]);
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x-1][y+1][z+1]);
				} else if (d==2) {
					val1 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y+1][z-1]);
					val2 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y+1][z-1]);
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					direction = d;
				}
			}
		}
		if (sign>0) return minval;
		else return 0.0f;
	}
}
