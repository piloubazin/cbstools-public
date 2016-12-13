package de.mpg.cbs.jist.filter;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataDouble;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.analysis.*;

import java.io.*;

public class JistFilterRecursiveRidgeDiffusion extends ProcessingAlgorithm {
	
	// jist containers
	private ParamVolume inputImage;
	private ParamOption brightParam;
	private static final String[] brightTypes = {"bright","dark"};
	
	private ParamVolume surfaceImage;
	private ParamOption	orientationParam;
	private static final String[] orientationTypes = {"parallel","orthogonal","undefined"};
	private ParamFloat	angleParam;
	
	private ParamFloat scalingParam;
	private ParamInteger nscalesParam;
	
	//private ParamFloat thresholdParam;
	
	private ParamOption filterParam;
	private static final String[] filterTypes = {"2D","1D"};
	
	private ParamOption propagationParam;
	private static final String[] propagationTypes = {"diffusion","belief"};
	private ParamFloat factorParam;
	private ParamFloat diffParam;
	private ParamInteger iterParam;
	
	private ParamVolume filterImage;
	private ParamVolume probaImage;
	private ParamVolume scaleImage;
	private ParamVolume directionImage;
	
	// global variables
	int nx, ny, nz, nc, nxyz;
	int nscales;
	
	// numerical quantities
	private static final	float	INVSQRT2 = (float)(1.0/FastMath.sqrt(2.0));
	private static final	float	INVSQRT3 = (float)(1.0/FastMath.sqrt(3.0));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0);
	private static final	float	SQRT2PI = (float)FastMath.sqrt(2.0*(float)Math.PI);
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private static final	float   	L2N2=2.0f*(float)(FastMath.sqrt(2.0*(float)(FastMath.log(2.0))));
	
	// direction labeling		
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;
	public	static 	final 	byte 	XpY = 3;
	public	static 	final 	byte 	YpZ = 4;
	public	static 	final 	byte 	ZpX = 5;
	public	static 	final 	byte 	XmY = 6;
	public	static 	final 	byte 	YmZ = 7;
	public	static 	final 	byte 	ZmX = 8;
	public	static 	final 	byte 	XpYpZ = 9;
	public	static 	final 	byte 	XmYmZ = 10;
	public	static 	final 	byte 	XmYpZ = 11;
	public	static 	final 	byte 	XpYmZ = 12;
	public	static	final	byte	mX = 13;
	public	static	final	byte	mY = 14;
	public	static	final	byte	mZ = 15;
	public	static 	final 	byte 	mXpY = 16;
	public	static 	final 	byte 	mYpZ = 17;
	public	static 	final 	byte 	mZpX = 18;
	public	static 	final 	byte 	mXmY = 19;
	public	static 	final 	byte 	mYmZ = 20;
	public	static 	final 	byte 	mZmX = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	
	public	static 	final 	byte 	NC = 13;
	public	static 	final 	byte 	NC2 = 26;
	
	private static final boolean debug=true;
	
	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(brightParam = new ParamOption("Ridge intensities", brightTypes));
		inputParams.add(filterParam = new ParamOption("Ridge filter", filterTypes));
		inputParams.add(surfaceImage = new ParamVolume("Surface(s) level sets"));
		surfaceImage.setMandatory(false);
		inputParams.add(orientationParam = new ParamOption("Orientation relative to surface", orientationTypes));
		inputParams.add(scalingParam = new ParamFloat("Angular factor", 0.0f, 10.0f, 2.0f));
		
		
		//inputParams.add(thresholdParam = new ParamFloat("Probability threshold", 0.0, 1.0, 0.5));
		//inputParams.add(scalingParam = new ParamFloat("Scale factor ", 0.25f, 2.0f, 1.0f));
		inputParams.add(nscalesParam = new ParamInteger("Number of scales ", 0, 10, 3));

		inputParams.add(propagationParam = new ParamOption("Propagation model", propagationTypes));
		inputParams.add(factorParam = new ParamFloat("Diffusion factor", 0.0f, 100.0f, 0.5f));
		inputParams.add(diffParam = new ParamFloat("Max difference", 0.0f, 1.0f, 0.001f));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000, 100));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Filter.devel");
		inputParams.setLabel("Recursive Ridge Diffusion");
		inputParams.setName("RecursiveRidgeDiffusion");
		
		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Extracts planar of tubular structures across multiple scales, with an optional directional bias.");
		
		info.setVersion("3.0.9");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}
	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(filterImage = new ParamVolume("Filter response",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(probaImage = new ParamVolume("Probability response",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(scaleImage = new ParamVolume("Detection scale",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(directionImage = new ParamVolume("Ridge direction",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("Recursive Ridge Diffusion");
		outputParams.setLabel("Recursive Ridge Diffusion");
	}
	
	@Override
	protected void execute(CalculationMonitor monitor){
		// import the image data into 1D arrays
		float[] image = Interface.getFloatImage3D(inputImage);
		int[] dim = Interface.getDimensions(inputImage);
		nx = dim[X]; ny = dim[Y]; nz = dim[Z];
		nxyz = nx*ny*nz;
		
		// load surface, if used
		float[] surface = null;
		nc = 0;
		String orientation = orientationParam.getValue();
		if (orientation.equals("parallel") || orientation.equals("orthogonal")) {
			nc = Interface.getComponents(surfaceImage);
			surface = Interface.getFloatImage4D(surfaceImage);
		}
		
		// load scale values
		nscales = nscalesParam.getValue().intValue();
		float scaling = scalingParam.getValue().floatValue();
		
		boolean[] mask = new boolean[nxyz];
		float minI = 1e9f;
		float maxI = -1e9f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			// mask
			if (image[id]==0) mask[id] = false;
			else mask[id] = true;
			// remove border from computations
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[id] = false;
			// normalize
			if (image[id]>maxI) maxI = image[id];
			if (image[id]<minI) minI = image[id];
		}
		
		// normalize, invert image if looking for dark features
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (brightParam.getValue().equals("bright"))
				image[xyz] = (image[xyz]-minI)/(maxI-minI);
			else if (brightParam.getValue().equals("dark"))
				image[xyz] = (maxI-image[xyz])/(maxI-minI);
		}
		
		// Compute filter at different scales
		// new filter response from raw image		
		float[] maxresponse = new float[nxyz];
		byte[] maxdirection = new byte[nxyz];
		byte[] maxscale = new byte[nxyz];
		
		directionFromRecursiveRidgeFilter(image, mask, maxresponse, maxdirection);
		
		for(int i=0;i<nscales;i++){
			
			float scale = 1.0f+i;
			float[] smoothed = new float[nxyz];

			// Gaussian Kernel
			float[][] G = ImageFilters.separableGaussianKernel(scale/L2N2,scale/L2N2,scale/L2N2);
			int gx = (G[0].length-1)/2;
			int gy = (G[1].length-1)/2;
			int gz = (G[2].length-1)/2;
				
			// smoothed image
			smoothed = ImageFilters.separableConvolution(image,nx,ny,nz,G,gx,gy,gz); 

			byte[] direction = new byte[nxyz];
			float[] response = new float[nxyz];
			directionFromRecursiveRidgeFilter(smoothed, mask, response, direction);
			
			//Combine scales: keep maximum response
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				if (response[xyz]>maxresponse[xyz]) {
					maxresponse[xyz] = response[xyz];
					maxdirection[xyz] = direction[xyz];
					maxscale[xyz] = (byte)(1+i);
				}
			}
		}
		
		// Attenuate response according to orientation
		if (orientation.equals("parallel")|| orientation.equals("orthogonal")) {
			float factor = scalingParam.getValue().floatValue();
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float[] vec = directionVector(maxdirection[xyz]);
				// find the closest surfaces; weighted mean
				int closest = 0;
				float dist = 1e9f;
				for (int c=1;c<nc;c++) {
					if (Numerics.abs(surface[xyz+c*nxyz])<dist) {
						closest = c;
						dist = Numerics.abs(surface[xyz+c*nxyz]);
					}
				}
				int second = closest;
				if (closest<nc-1 && surface[xyz+closest*nxyz]*surface[xyz+(closest+1)*nxyz]<0) second = closest+1; 
				if (closest>0 && surface[xyz+closest*nxyz]*surface[xyz+(closest-1)*nxyz]<0) second = closest-1; 
				
				float w1 = 1.0f, w2 = 0.0f;
				if (second!=closest) {
					w1 = surface[xyz+second*nxyz]/(surface[xyz+second*nxyz]-surface[xyz+closest*nxyz]);
					w2 = surface[xyz+closest*nxyz]/(surface[xyz+closest*nxyz]-surface[xyz+second*nxyz]);
				}
				double[] grad = new double[3];
				// central differences
				grad[X] = w1*0.5f*(surface[xyz+closest*nxyz+1]-surface[xyz+closest*nxyz-1]);
				grad[Y] = w1*0.5f*(surface[xyz+closest*nxyz+nx]-surface[xyz+closest*nxyz-nx]);
				grad[Z] = w1*0.5f*(surface[xyz+closest*nxyz+nx*ny]-surface[xyz+closest*nxyz-nx*ny]);
				if (second!=closest) {
					grad[X] += w2*0.5f*(surface[xyz+second*nxyz+1]-surface[xyz+second*nxyz-1]);
					grad[Y] += w2*0.5f*(surface[xyz+second*nxyz+nx]-surface[xyz+second*nxyz-nx]);
					grad[Z] += w2*0.5f*(surface[xyz+second*nxyz+nx*ny]-surface[xyz+second*nxyz-nx*ny]);
				}					
				double norm = grad[X]*grad[X]+grad[Y]*grad[Y]+grad[Z]*grad[Z];
				if (norm>0) {
					norm = FastMath.sqrt(norm);
					grad[X] /= norm;
					grad[Y] /= norm;
					grad[Z] /= norm;
				}
				double correlation = FastMath.pow(grad[X]*vec[X]+grad[Y]*vec[Y]+grad[Z]*vec[Z],factor);
				if (orientation.equals("parallel") && filterParam.getValue().equals("1D")) {
					correlation = 1.0-correlation;
				} else if (orientation.equals("orthogonal") && filterParam.getValue().equals("2D")) {
					correlation = 1.0-correlation;
				}
				maxresponse[xyz] *= (float)correlation;
			}
		}

		// Equalize histogram (Exp Median)
		float[] proba = new float[nxyz];
		probabilityFromRecursiveRidgeFilter(maxresponse, proba);	
		
		// generate a direction vector
		float[][] direction = new float[3][nxyz];
		for (int xyz=0;xyz<nxyz;xyz++){
			float[] vec = directionVector(maxdirection[xyz]);
			direction[X][xyz] = vec[X];
			direction[Y][xyz] = vec[Y];
			direction[Z][xyz] = vec[Z];
		}
		
		// 3. diffuse the data to neighboring structures
		// compute neighborhood structure
		Interface.displayMessage("compute neighborhood structure\n");
		float[] similarity1 = new float[nxyz];
		float[] similarity2 = new float[nxyz];
		byte[] neighbor1 = new byte[nxyz];
		byte[] neighbor2 = new byte[nxyz];
		//estimateDiffusionSimilarity(direction, mask, neighbor1, neighbor2, similarity1, similarity2);
		estimateSimpleDiffusionSimilarity(direction, mask, neighbor1, neighbor2, similarity1, similarity2);
		
		// original score
		float[] diffused = new float[nxyz];
		float shapeScale = 1.0f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (!mask[xyz]) proba[xyz] = 0.0f;

			//diffused[xyz] = proba[xyz];
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}
		
		// pre-compute different similarity models; default is angle only
		for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
			// value more neighbors that are well aligned and have strong tubular response
			similarity1[xyz] = (float)FastMath.sqrt(similarity1[xyz]*proba[neighborIndex(neighbor1[xyz],xyz)]);
			similarity2[xyz] = (float)FastMath.sqrt(similarity2[xyz]*proba[neighborIndex(neighbor2[xyz],xyz)]);
		}
		
		// diffuse along neighborhood
		int maxiter = iterParam.getValue().intValue();
		float maxdiff = diffParam.getValue().floatValue();
		float factor = factorParam.getValue().floatValue();
		
		if (propagationParam.getValue().equals("diffusion")) {
			// diffuse along the tensor
			Interface.displayMessage("diffusion smoothing ("+factorParam.getValue().floatValue()+")\n");
		
			proba = probaDiffusion(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff, factor);
		} else if (propagationParam.getValue().equals("belief")) {
			// belief propagation
			Interface.displayMessage("belief propagation\n");
		
			proba = beliefPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
		}
				
		// Output
		Interface.displayMessage("Prepare output images\n");
		String outname = Interface.getName(inputImage);
		ImageHeader outheader = Interface.getHeader(inputImage);
		
		if (filterParam.getValue().equals("1D")) outname+="_rrf1d";
		else outname+="_rrf2d";
		
		if (propagationParam.getValue().equals("diffusion") && iterParam.getValue().intValue()>0) outname += "d";
		else if (propagationParam.getValue().equals("belief") && iterParam.getValue().intValue()>0) outname += "bp";
		
		Interface.displayMessage("filter out\n");
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = maxresponse[id];
		}					
		ImageDataFloat resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_filter");
		filterImage.setValue(resData);
		
		Interface.displayMessage("proba out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = proba[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_proba");
		probaImage.setValue(resData);
		
		Interface.displayMessage("scale out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = maxscale[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_scale");
		scaleImage.setValue(resData);
		
		Interface.displayMessage("direction out\n");
		float[][][][] result4 = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result4[x][y][z][X] = direction[X][id];
			result4[x][y][z][Y] = direction[Y][id];
			result4[x][y][z][Z] = direction[Z][id];
		}					
		resData = new ImageDataFloat(result4);	
		resData.setHeader(outheader);
		resData.setName(outname+"_dir");
		directionImage.setValue(resData);
		
		return;
	}
	
	private final void directionFromRecursiveRidgeFilter(float[] img, boolean[] mask, float[] filter,byte[] direction) {
			
			// get the tubular filter response
			float[][][] planescore = new float[nx][ny][nz];
			byte[][][] planedir = new byte[nx][ny][nz];
			byte[][][] linedir = new byte[nx][ny][nz];
			float[][][] image = new float[nx][ny][nz];
			//float[] filter = new float[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				image[x][y][z] = img[xyz];
			}
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				filter[xyz] = 0.0f;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					// check for zero-valued neighbors as well
					minmaxplaneScore(image, planescore, planedir, x,y,z, 13);
					// sign issue: remove all that have different sign, keep global sign
					float linescore = minmaxlineScore(image, planedir, linedir, x,y,z, 4);
					if (planescore[x][y][z]*linescore>0) {
						filter[xyz] = Numerics.sign(linescore)*Numerics.sqrt(planescore[x][y][z]*linescore);
						direction[xyz] = linedir[x][y][z];
						if(filter[xyz]<0) { filter[xyz]=0; direction[xyz] = -1; }
					} else {
						filter[xyz] = 0.0f;
						direction[xyz] = -1;
					}
				}
			}
			planescore = null;
			planedir = null;
			linedir = null;
			image = null;
			return;
	}
	private final void probabilityFromRecursiveRidgeFilter( float[] filter, float[] shape) {
		// normalization: best is the iterative robust exponential (others are removed)
		int nb = 0;
		double max = 0.0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
				// keep only the proper sign
				if (filter[xyz]<=0) filter[xyz] = 0.0f;
				else {
					// fit exp only to non-zero data
					nb++;
					max = Numerics.max(filter[xyz], max);
				}
		}
		
		// robust measures
		double[] response = new double[nb];
		int n=0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
				if (filter[xyz]>0) {
					response[n] = filter[xyz];
					n++;
				}
		}
		Percentile measure = new Percentile();
		double median = measure.evaluate(response, 50.0);
		double beta = median/FastMath.log(2.0);
		
		Interface.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			int xyz = x + nx*y + nx*ny*z;
			if (filter[xyz]>0) {
					double pe = FastMath.exp( -filter[xyz]/beta)/beta;
					shape[xyz] = (float)(1.0/max/( 1.0/max+pe));
			}
		}
		return;
	}
	boolean zeroNeighbor(float[] image, boolean[] mask, int x, int y, int z, int d) {
			for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
				if (image[x+i+nx*(y+j)+nx*ny*(z+l)]!=image[x+nx*y+nx*ny*z] && i*i+j*j+l*l<=2*d*d) return false;
			}
			return true;
		}
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
	
	float minmaxlineScore(float[][][] image, byte[][][] planedir, byte[][][] linedir, int x, int y, int z, int dmax) {
		float maxgrad = 0.0f; 
		float minval = 0.0f; 
		float sign = 0.0f;
		byte direction = -1;
		
		float val1 = 0.0f, val2 = 0.0f;
		if (planedir[x][y][z]==0) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x][y][z-1]	-image[x][y-1][z-1]	-image[x][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x][y][z+1]	-image[x][y-1][z+1]	-image[x][y+1][z+1];
					direction = Y;
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x][y-1][z]	-image[x][y-1][z-1]	-image[x][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x][y+1][z+1]	-image[x][y+1][z-1]	-image[x][y+1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z-1]	+image[x][y+1][z+1]
							-image[x][y-1][z+1]	-image[x][y-2][z]	-image[x][y][z+2];
					val2 = 	 image[x][y][z]		+image[x][y-1][z-1]	+image[x][y+1][z+1]
							-image[x][y+1][z-1]	-image[x][y][z-2]	-image[x][y+2][z];
					direction = YpZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z+1]	+image[x][y+1][z-1]
							-image[x][y-1][z-1]	-image[x][y-2][z]	-image[x][y][z-2];
					val2 = 	 image[x][y][z]		+image[x][y+1][z-1]	+image[x][y-1][z+1]
							-image[x][y+1][z+1]	-image[x][y+2][z]	-image[x][y][z+2];
					direction = YmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==1) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x-1][y][z]	-image[x-1][y][z-1]	-image[x-1][y][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]	+image[x][y][z+1]
							-image[x+1][y][z]	-image[x+1][y][z-1]	-image[x+1][y][z+1];
					direction = Z;		
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y][z-1]	-image[x-1][y][z-1]	-image[x+1][y][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y][z+1]	-image[x-1][y][z+1]	-image[x+1][y][z+1];
					direction = X;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z-1]	+image[x+1][y][z+1]
							-image[x-1][y][z+1]	-image[x-2][y][z]	-image[x][y][z+2];
					val2 = 	 image[x][y][z]		+image[x-1][y][z-1]	+image[x+1][y][z+1]
							-image[x+1][y][z-1]	-image[x][y][z-2]	-image[x+2][y][z];
					direction = ZpX;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z+1]	+image[x+1][y][z-1]
							-image[x-1][y][z-1]	-image[x-2][y][z]	-image[x][y][z-2];
					val2 = 	 image[x][y][z]		+image[x-1][y][z+1]	+image[x+1][y][z-1]
							-image[x+1][y][z+1]	-image[x][y][z+2]	-image[x+2][y][z];
					direction = ZmX;		
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==2) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y-1][z]	-image[x-1][y-1][z]	-image[x+1][y-1][z];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]	+image[x+1][y][z]
							-image[x][y+1][z]	-image[x-1][y+1][z]	-image[x+1][y+1][z];
					direction = X;		
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x-1][y][z]	-image[x-1][y-1][z]	-image[x-1][y+1][z];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]	+image[x][y+1][z]
							-image[x+1][y][z]	-image[x+1][y-1][z]	-image[x+1][y+1][z];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]		+image[x-1][y-1][z]	+image[x+1][y+1][z]
							-image[x-1][y+1][z]	-image[x-2][y][z]	-image[x][y+2][z];
					val2 = 	 image[x][y][z]		+image[x-1][y-1][z]	+image[x+1][y+1][z]
							-image[x+1][y-1][z]	-image[x][y-2][z]	-image[x+2][y][z];
					direction = XpY;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]		+image[x-1][y+1][z]	+image[x+1][y-1][z]
							-image[x-1][y-1][z]	-image[x-2][y][z]	-image[x][y-2][z];
					val2 = 	 image[x][y][z]		+image[x-1][y+1][z]	+image[x+1][y-1][z]
							-image[x+1][y+1][z]	-image[x][y+2][z]	-image[x+2][y][z];
					direction = XmY;		
				}	
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==3) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y+1][z]		+image[x+1][y-1][z]
							-image[x][y][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y+1][z]		+image[x+1][y-1][z]
							-image[x][y][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1];
					direction = XmY;		
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x-1][y+1][z]	-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x+1][y-1][z]	-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z];
					direction = XmYpZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
							-image[x-1][y+1][z-1]	-image[x-2][y+2][z]		-image[x][y][z-2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
							-image[x+1][y-1][z+1]	-image[x][y][z+2]		-image[x+2][y-2][z];
					direction = XmYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==4) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x][y+1][z-1]		+image[x][y-1][z+1]
							-image[x-1][y][z]	-image[x-1][y+1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y+1][z-1]		+image[x][y-1][z+1]
							-image[x+1][y][z]	-image[x+1][y+1][z-1]	-image[x+1][y-1][z+1];
					direction = YmZ;		
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y+1][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y-1][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1];
					direction = X;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z];
					direction = XmYpZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x+1][y-1][z+1]	-image[x+2][y][z]		-image[x][y-2][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x-1][y+1][z-1]	-image[x][y+2][z-2]		-image[x-2][y][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==5) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x+1][y][z-1]		+image[x-1][y][z+1]
							-image[x][y-1][z]	-image[x+1][y-1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x+1][y][z-1]		+image[x-1][y][z+1]
							-image[x][y+1][z]	-image[x+1][y+1][z-1]	-image[x-1][y+1][z+1];
					direction = ZmX;			
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x+1][y][z-1]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x-1][y][z+1]	-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z];
					direction = XmYmZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x-1][y+1][z+1]	-image[x][y+2][z]		-image[x-2][y][z+2];
					val2 = 	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
							-image[x+1][y-1][z-1]	-image[x+2][y][z-2]		-image[x][y-2][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==6) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 = 	 image[x][y][z]		+image[x-1][y-1][z]		+image[x+1][y+1][z]
							-image[x][y][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1];
					val2 = 	 image[x][y][z]		+image[x-1][y-1][z]		+image[x+1][y+1][z]
							-image[x][y][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1];
					direction = XpY;		
				} else if (d==1) {
					val1 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x-1][y-1][z]	-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1];
					val2 = 	 image[x][y][z]		+image[x][y][z-1]		+image[x][y][z+1]
							-image[x+1][y+1][z]	-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1];
					direction = Z;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
							-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
							-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z];
					direction = XpYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==7) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]		+image[x][y-1][z-1]		+image[x][y+1][z+1]
							-image[x-1][y][z]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1];
					val2 =	 image[x][y][z]		+image[x][y-1][z-1]		+image[x][y+1][z+1]
							-image[x+1][y][z]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1];
					direction = YpZ;		
				} else if (d==1) {
					val1 =	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y-1][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1];
					val2 =	 image[x][y][z]		+image[x-1][y][z]		+image[x+1][y][z]
							-image[x][y+1][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1];
					direction = X;		
				} else if (d==2) {
					val1 =   image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z];
					val2 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 =   image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z];
					val2 =	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
							-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2];
					direction = XmYmZ;		
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==8) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]		+image[x-1][y][z-1]		+image[x+1][y][z+1]
							-image[x][y-1][z]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1];
					val2 =	 image[x][y][z]		+image[x-1][y][z-1]		+image[x+1][y][z+1]
							-image[x][y+1][z]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1];
					direction = ZpX;		
				} else if (d==1) {
					val1 =	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x-1][y][z-1]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1];
					val2 =	 image[x][y][z]		+image[x][y-1][z]		+image[x][y+1][z]
							-image[x+1][y][z+1]	-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1];
					direction = Y;		
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z];
					val2 = 	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
							-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2];
					direction = XpYpZ;		
				} else if (d==3) {
					val1 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z];
					val2 = 	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
							-image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2];
					direction = XmYpZ;		
				}					
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==9) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x-1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x+1][y-1][z-1]);
					direction = XpYmZ;						  
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1]
											  -1.5f*(image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]);
					direction = XmYpZ;						  
				} else if (d==2) {
					val1 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x+1][y-1][z+1]);
					val2 = 	 image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x-1][y+1][z-1]);
					direction = XmYmZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==10) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x-1][y-1][z-1]);
					direction = XmYpZ;						  
				} else 	if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x-1][y+1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x+1][y-1][z+1]);
					direction = XpYmZ;						  
				} else 	if (d==2) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x+1][y+1][z-1]);
					direction = XpYpZ;						  
				}		
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}			
		} else if (planedir[x][y][z]==11) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x+1][y+1][z-1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x-1][y-1][z+1]	+image[x-1][y-1][z-1]);
					direction = XmYmZ;						  
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z-1]	+image[x-1][y-1][z+1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x-1][y+1][z+1]);
					direction = XpYmZ;						  
				} else 	if (d==2) {
					val1 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y-1][z+1]);
					val2 =	 image[x][y][z]			+image[x+1][y+1][z+1]	+image[x-1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]);
					direction = XpYpZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		} else if (planedir[x][y][z]==12) {
			for (byte d=0;d<dmax;d++) {
				if (d==0) {
					val1 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x+1][y-1][z+1]	+image[x+1][y+1][z+1]);
					val2 =	 image[x][y][z]			+image[x-1][y+1][z+1]	+image[x+1][y-1][z-1]
											  -1.5f*(image[x-1][y+1][z-1]	+image[x-1][y-1][z-1]);
					direction = XmYmZ;						  
				} else if (d==1) {
					val1 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]);
					val2 =	 image[x][y][z]			+image[x+1][y-1][z+1]	+image[x-1][y+1][z-1]
											  -1.5f*(image[x+1][y+1][z+1]	+image[x-1][y+1][z+1]);
					direction = XmYpZ;						  
				} else if (d==2) {
					val1 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y+1][z-1]);
					val2 =	 image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1]
											  -1.5f*(image[x-1][y+1][z+1]	+image[x-1][y+1][z-1]);
					direction = XpYpZ;						  
				}
				// find the strongest gradient direction, then estimate the corresponding filter response
				if (val1*val1+val2*val2>maxgrad) {
					maxgrad = val1*val1+val2*val2;
					if (val1*val1<val2*val2) minval = val1;
					else minval = val2;
					sign = val1*val2;
					linedir[x][y][z] = direction;
				}
			}
		}
		if (sign>0) return minval;
		else return 0.0f;
	}
	private final float[] directionVector(int d) {
		if (d==X) return new float[]{1.0f, 0.0f, 0.0f};
		else if (d==Y) return new float[]{0.0f, 1.0f, 0.0f};
		else if (d==Z) return new float[]{0.0f, 0.0f, 1.0f};
		else if (d==XpY) return new float[]{INVSQRT2, INVSQRT2, 0.0f};
		else if (d==YpZ) return new float[]{0.0f, INVSQRT2, INVSQRT2};
		else if (d==ZpX) return new float[]{INVSQRT2, 0.0f, INVSQRT2};
		else if (d==XmY) return new float[]{INVSQRT2, -INVSQRT2, 0.0f};
		else if (d==YmZ) return new float[]{0.0f, INVSQRT2, -INVSQRT2};
		else if (d==ZmX) return new float[]{-INVSQRT2, 0.0f, INVSQRT2};
		else if (d==XpYpZ) return new float[]{INVSQRT3, INVSQRT3, INVSQRT3};
		else if (d==XmYmZ) return new float[]{INVSQRT3, -INVSQRT3, -INVSQRT3};
		else if (d==XmYpZ) return new float[]{INVSQRT3, -INVSQRT3, INVSQRT3};
		else if (d==XpYmZ) return new float[]{INVSQRT3, INVSQRT3, -INVSQRT3};
		else return new float[]{0.0f, 0.0f, 0.0f};
	}
	private final byte[] directionNeighbor(int d) {
		if (d==X) return new byte[]{1, 0, 0};
		else if (d==Y) return new byte[]{0, 1, 0};
		else if (d==Z) return new byte[]{0, 0, 1};
		else if (d==XpY) return new byte[]{1, 1, 0};
		else if (d==YpZ) return new byte[]{0, 1, 1};
		else if (d==ZpX) return new byte[]{1, 0, 1};
		else if (d==XmY) return new byte[]{1, -1, 0};
		else if (d==YmZ) return new byte[]{0, 1, -1};
		else if (d==ZmX) return new byte[]{-1, 0, 1};
		else if (d==XpYpZ) return new byte[]{1, 1, 1};
		else if (d==XmYmZ) return new byte[]{1, -1, -1};
		else if (d==XmYpZ) return new byte[]{1, -1, 1};
		else if (d==XpYmZ) return new byte[]{1, 1, -1};
		else return new byte[]{0, 0, 0};
	}
	
	private final float[] beliefPropagation(float[] proba, float[] similarity1, float[] similarity2, byte[] ngb1, byte[] ngb2, boolean[] mask, int nxyz, int iter, float maxdiff) {
		
		float[] message1fg = new float[nxyz];
		float[] message2fg = new float[nxyz];
		float[] message1bg = new float[nxyz];
		float[] message2bg = new float[nxyz];
		
		float[] newmsg1fg = new float[nxyz];
		float[] newmsg2fg = new float[nxyz];
		float[] newmsg1bg = new float[nxyz];
		float[] newmsg2bg = new float[nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) {
			message1fg[xyz] = 1.0f;
			message2fg[xyz] = 1.0f;
			message1bg[xyz] = 1.0f;
			message2bg[xyz] = 1.0f;
			newmsg1fg[xyz] = 1.0f;
			newmsg2fg[xyz] = 1.0f;
			newmsg1bg[xyz] = 1.0f;
			newmsg2bg[xyz] = 1.0f;
		}
		
		// min value for probas, similarities
		float minp = 0.00f;
		float mins = 0.00f;
		
		// message passing
		for (int t=0;t<iter;t++) {
			Interface.displayMessage("iteration "+(t+1)+": ");
			float diff = 0.0f;
			//float prev;
			int ngb;
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) if (proba[xyz]>0 && proba[xyz]<1) {
				float simfgfg1 = similarity1[xyz];
				float simbgfg1 = 0.0f; //(1.0f-similarity1[xyz]); //similarity1[xyz];
				float simfgbg1 = (1.0f-similarity1[xyz]); // good results
				float simbgbg1 = 1.0f; //(1.0f-similarity1[xyz]); //1.0f; //similarity1[xyz];
				
				float simfgfg2 = similarity2[xyz];
				float simbgfg2 = 0.0f;//(1.0f-similarity2[xyz]); //similarity2[xyz];
				float simfgbg2 = (1.0f-similarity2[xyz]); // good results
				float simbgbg2 = 1.0f; //(1.0f-similarity2[xyz]); //1.0f; //similarity2[xyz];
				
				// first neighbor message, fg label
				ngb = neighborIndex(ngb1[xyz],xyz);
				float mfgfg1 = Numerics.max(mins,simfgfg1)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1fg[ngb]>0) mfgfg1 *= message1fg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2fg[ngb]>0) mfgfg1 *= message2fg[ngb];
				
				float mfgbg1 = Numerics.max(mins,simfgbg1)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1bg[ngb]>0) mfgbg1 *= message1bg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2bg[ngb]>0) mfgbg1 *= message2bg[ngb];
				
				newmsg1fg[xyz] = Numerics.max(mfgfg1, mfgbg1);
				diff = Numerics.max(diff, Numerics.abs(newmsg1fg[xyz]-message1fg[xyz]));
				
				// first neighbor message, bg label
				ngb = neighborIndex(ngb1[xyz],xyz);
				float mbgbg1 = Numerics.max(mins,simbgbg1)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1bg[ngb]>0) mbgbg1 *= message1bg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2bg[ngb]>0) mbgbg1 *= message2bg[ngb];
				
				float mbgfg1 = Numerics.max(mins,simbgfg1)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1fg[ngb]>0) mbgfg1 *= message1fg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2fg[ngb]>0) mbgfg1 *= message2fg[ngb];
				
				newmsg1bg[xyz] = Numerics.max(mbgbg1, mbgfg1);
				diff = Numerics.max(diff, Numerics.abs(newmsg1bg[xyz]-message1bg[xyz]));
				
				// second neighbor message, fg label
				ngb = neighborIndex(ngb2[xyz],xyz);
				float mfgfg2 = Numerics.max(mins,simfgfg2)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1fg[ngb]>0) mfgfg2 *= message1fg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2fg[ngb]>0) mfgfg2 *= message2fg[ngb];
				
				float mfgbg2 = Numerics.max(mins,simfgbg2)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1bg[ngb]>0) mfgbg2 *= message1bg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2bg[ngb]>0) mfgbg2 *= message2bg[ngb];
				
				newmsg2fg[xyz] = Numerics.max(mfgfg2, mfgbg2);
				diff = Numerics.max(diff, Numerics.abs(newmsg2fg[xyz]-message2fg[xyz]));
				
				// second neighbor message, bg label
				ngb = neighborIndex(ngb2[xyz],xyz);
				float mbgbg2 = Numerics.max(mins,simbgbg2)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1bg[ngb]>0) mbgbg2 *= message1bg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2bg[ngb]>0) mbgbg2 *= message2bg[ngb];
				
				float mbgfg2 = Numerics.max(mins,simbgfg2)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz && message1fg[ngb]>0) mbgfg2 *= message1fg[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz && message2fg[ngb]>0) mbgfg2 *= message2fg[ngb];
				
				newmsg2bg[xyz] = Numerics.max(mbgbg2, mbgfg2);
				diff = Numerics.max(diff, Numerics.abs(newmsg2bg[xyz]-message2bg[xyz]));
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) if (proba[xyz]>0 && proba[xyz]<1){
				message1fg[xyz] = newmsg1fg[xyz];
				message1bg[xyz] = newmsg1bg[xyz];
				message2fg[xyz] = newmsg2fg[xyz];
				message2bg[xyz] = newmsg2bg[xyz];
			}
			Interface.displayMessage("diff "+diff+"\n");
			if (diff < maxdiff) t = iter;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz];
			
			if (mask[xyz]) if (proba[xyz]>0 && proba[xyz]<1) {
				if (message1fg[xyz]>0) belief[xyz] *= message1fg[xyz];
				if (message2fg[xyz]>0) belief[xyz] *= message2fg[xyz];
			}
			
			bgbelief = (1.0f-proba[xyz]);
			if (mask[xyz]) if (proba[xyz]>0 && proba[xyz]<1) {
				if (message1bg[xyz]>0) bgbelief *= message1bg[xyz];
				if (message2bg[xyz]>0) bgbelief *= message2bg[xyz];
			}
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
			belief[xyz] = belief[xyz]/(belief[xyz]+bgbelief);
		}
		return belief;
	}
	
	private final float[] probaDiffusion(float[] proba, float[] similarity1, float[] similarity2, byte[] neighbor1, byte[] neighbor2, boolean[] mask, int nxyz, int maxiter, float maxdiff, float factor) {
		
		float[] diffused = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}
		
		for (int t=0;t<maxiter;t++) {
			Interface.displayMessage("iteration "+(t+1)+": ");
			float diff = 0.0f;
			/*
			float globalw = 0.0f;
			float den = 0.0f;
			float ratio = 1.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float ngb1 = diffused[neighborIndex(neighbor1[xyz],xyz)];
				float ngb2 = diffused[neighborIndex(neighbor2[xyz],xyz)];
				float w1 = similarity1[xyz]*probaWeight(ngb1);
				float w2 = similarity2[xyz]*probaWeight(ngb2);
				globalw += w1+w2;
				den += 1.0f;
			}
			ratio = globalw/den;
			*/
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float prev = diffused[xyz];
				// weight with distance to 0 or 1
				float ngb1 = diffused[neighborIndex(neighbor1[xyz],xyz)];
				float ngb2 = diffused[neighborIndex(neighbor2[xyz],xyz)];
				//float w1 = similarity1[xyz]*probaWeight(ngb1);
				//float w2 = similarity2[xyz]*probaWeight(ngb2);
				
				// remap neighbors?
				ngb1 = (float)FastMath.exp(ngb1)-1.0f;
				ngb2 = (float)FastMath.exp(ngb2)-1.0f;
				
				// now this depends on similarity choices
				// value more neighbors that are well aligned and have strong tubular response
				//float w1 = similarity1[xyz]*proba[neighborIndex(neighbor1[xyz],xyz)];
				//float w2 = similarity2[xyz]*proba[neighborIndex(neighbor2[xyz],xyz)];
				float w1 = similarity1[xyz];
				float w2 = similarity2[xyz];
				
				// stable normalization	 (works OK, but not much different from the raw filter?)			 
				//diffused[xyz] = Numerics.bounded((proba[xyz] + factor*w1*ngb1 + factor*w2*ngb2)/(1.0f + factor*globalw), 0, 1);
				// integration over the whole vessel (log version is more stable (??) )
				diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz] + factor*w1*ngb1 + factor*w2*ngb2);
				//diffused[xyz] = proba[xyz] + factor*w1*ngb1 + factor*w2*ngb2;
				
				if (Numerics.abs(prev-diffused[xyz])>diff) diff = Numerics.abs(prev-diffused[xyz]);
			}
			
			Interface.displayMessage("diff "+diff+"\n");
			if (diff<maxdiff) t=maxiter;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (maxiter==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		return diffused;
	}
    private final void estimateSimpleDiffusionSimilarity(float[][] imdir, boolean[] mask, byte[] neighbor1, byte[] neighbor2, 
    													float[] similarity1, float[] similarity2) {
    	
		float[]	ds1 = new float[NC2];
		float[]	ds2 = new float[NC2];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (mask[id]) {
				// find the aligned discrete directions
				byte d1 = -1; float dv1 = 0.0f;
				byte d2 = -1; float dv2 = 0.0f;
				for (byte d=0;d<NC2;d++) {
					int idn = neighborIndex(d,id);
					
					if (mask[idn]) {
						float dv = directionProduct(d, id, imdir);
						if (dv>dv1) {
							dv1 = dv;
							d1 = d;
						}
						if (dv<dv2) {
							dv2 = dv;
							d2 = d;
						}
					}
				}
				if (d1!=-1) {
					neighbor1[id] = d1;
					similarity1[id] = simpleAngleSimilarityGain(imdir, d1,id,neighborIndex(d1,id));
				} else {
					neighbor1[id] = 0;
					similarity1[id] = 0.0f;
				}
				if (d2!=-1) {
					neighbor2[id] = d2;
					similarity2[id] = simpleAngleSimilarityGain(imdir, d2,id,neighborIndex(d2,id));
				} else {
					neighbor2[id] = 0;
					similarity2[id] = 0.0f;
				}
			}
		}
		
		return;
    

    }
    
	private final float simpleAngleSimilarityGain(float[][] imdir, int dir, int id1, int id2) {
		float prod = (imdir[0][id1]*imdir[0][id2]+imdir[1][id1]*imdir[1][id2]+imdir[2][id1]*imdir[2][id2]);
		float dataAngle = (float)FastMath.acos(Numerics.min(1.0f,Numerics.abs(prod)))/PI2;
				
		return (1.0f-dataAngle);
	}	

	private final int neighborIndex(byte d, int id) {
		int idn=id;
		
			 if (d==X) 	idn+=1; 		
		else if (d==mX)	idn-=1;
		else if (d==Y) 	idn+=nx;
		else if (d==mY)	idn-=nx;
		else if (d==Z) 	idn+=nx*ny;
		else if (d==mZ)	idn-=nx*ny;
		else if (d==XpY) 	idn+=1+nx;
		else if (d==mXpY) 	idn-=1+nx;
		else if (d==YpZ) 	idn+=nx+nx*nx;
		else if (d==mYpZ)	idn-=nx+nx*ny;
		else if (d==ZpX) 	idn+=1+nx*ny;	
		else if (d==mZpX)	idn-=1+nx*ny;
		else if (d==XmY) 	idn+=1-nx;	
		else if (d==mXmY)	idn-=1-nx;
		else if (d==YmZ) 	idn+=nx-nx*ny;
		else if (d==mYmZ)	idn-=nx-nx*ny;
		else if (d==ZmX) 	idn+=1-nx*ny;
		else if (d==mZmX)	idn-=1-nx*ny;
		else if (d==XpYpZ) 		idn+=1+nx+nx*ny;
		else if (d==mXpYpZ)		idn-=1+nx+nx*ny;
		else if (d==XmYmZ) 		idn+=1-nx-nx*ny; 
		else if (d==mXmYmZ)		idn-=1-nx-nx*ny;
		else if (d==XmYpZ) 		idn+=1-nx+nx*ny;
		else if (d==mXmYpZ)		idn-=1-nx+nx*ny;
		else if (d==XpYmZ) 		idn+=1+nx-nx*ny; 
		else if (d==mXpYmZ)		idn-=1+nx-nx*ny;

		return idn;
	}
	private final float directionProduct(int dir, int id, float[][] imdir) {
		float dv=0.0f;
			
			 if (dir==X) dv = imdir[0][id];
		else if (dir==Y) dv = imdir[1][id];
		else if (dir==Z) dv = imdir[2][id];
		else if (dir==XpY) dv = imdir[0][id]/SQRT2+imdir[1][id]/SQRT2;
		else if (dir==YpZ) dv = imdir[1][id]/SQRT2+imdir[2][id]/SQRT2;
		else if (dir==ZpX) dv = imdir[2][id]/SQRT2+imdir[0][id]/SQRT2;
		else if (dir==XmY) dv = imdir[0][id]/SQRT2-imdir[1][id]/SQRT2;
		else if (dir==YmZ) dv = imdir[1][id]/SQRT2-imdir[2][id]/SQRT2;
		else if (dir==ZmX) dv = imdir[2][id]/SQRT2-imdir[0][id]/SQRT2;
		else if (dir==XpYpZ) dv = imdir[0][id]/SQRT3+imdir[1][id]/SQRT3+imdir[2][id]/SQRT3;
		else if (dir==XmYpZ) dv = imdir[0][id]/SQRT3-imdir[1][id]/SQRT3+imdir[2][id]/SQRT3;
		else if (dir==XpYmZ) dv = imdir[0][id]/SQRT3+imdir[1][id]/SQRT3-imdir[2][id]/SQRT3;
		else if (dir==XmYmZ) dv = imdir[0][id]/SQRT3-imdir[1][id]/SQRT3-imdir[2][id]/SQRT3;
		else if (dir==mX) dv = -imdir[0][id];
		else if (dir==mY) dv = -imdir[1][id];
		else if (dir==mZ) dv = -imdir[2][id];
		else if (dir==mXpY) dv = -(imdir[0][id]/SQRT2+imdir[1][id]/SQRT2);
		else if (dir==mYpZ) dv = -(imdir[1][id]/SQRT2+imdir[2][id]/SQRT2);
		else if (dir==mZpX) dv = -(imdir[2][id]/SQRT2+imdir[0][id]/SQRT2);
		else if (dir==mXmY) dv = -(imdir[0][id]/SQRT2-imdir[1][id]/SQRT2);
		else if (dir==mYmZ) dv = -(imdir[1][id]/SQRT2-imdir[2][id]/SQRT2);
		else if (dir==mZmX) dv = -(imdir[2][id]/SQRT2-imdir[0][id]/SQRT2);
		else if (dir==mXpYpZ) dv = -(imdir[0][id]/SQRT3+imdir[1][id]/SQRT3+imdir[2][id]/SQRT3);
		else if (dir==mXmYpZ) dv = -(imdir[0][id]/SQRT3-imdir[1][id]/SQRT3+imdir[2][id]/SQRT3);
		else if (dir==mXpYmZ) dv = -(imdir[0][id]/SQRT3+imdir[1][id]/SQRT3-imdir[2][id]/SQRT3);
		else if (dir==mXmYmZ) dv = -(imdir[0][id]/SQRT3-imdir[1][id]/SQRT3-imdir[2][id]/SQRT3);

		return dv;
	}	
	
	private final void directionFromHessian(float[] img, boolean[] mask, float[] shape, byte[] direction) {
		double[][] hessian = new double[3][3];
			   
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int id = x + nx*y + nx*ny*z;
			if (mask[id] && !zeroNeighbor(img, mask, x,y,z,2)) {
				ImageGeometry.computeHessianOrder2At(img, id, hessian, nx, ny, nz);
				
				// compute eigenvalues
				double[] vals = Matrix3D.eigenvalues(hessian);
				double[][] dirs = Matrix3D.eigenvectors(hessian, vals);
				
				// shape: basic score is not good enough (too much background noise)
				if (vals[0]!=0) shape[id] = (float)(-vals[0]*(Numerics.abs(vals[0])-Numerics.abs(vals[2]))/Numerics.abs(vals[0]));
				
				// main direction: lowest eigenvalue : approx with closest neighbor
				double maxcorr = 0.0;
				byte maxd = -1;
				for (byte d=0;d<NC;d++) {
					float[] dird = directionVector(d);
					double corr = Numerics.abs(dirs[X][2]*dird[X]+dirs[Y][2]*dird[Y]+dirs[Z][2]*dird[Z]);
					if (corr>maxcorr) {
						maxcorr = corr;
						maxd = d;
					}
				}
				direction[id] = maxd;
			}
		}

		// remove differences of incorrect sign
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (shape[xyz]<0) shape[xyz] = 0.0f;
			}
		}

		return;
    }
    
	private final float tubularScore(float[] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y)+nx*ny*(z)]		+image[x+1+nx*(y)+nx*ny*(z)])
						-image[x+nx*(y-1)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z)]		-image[x+1+nx*(y-1)+nx*ny*(z)]
						-image[x+nx*(y+1)+nx*ny*(z)]		-image[x-1+nx*(y+1)+nx*ny*(z)]		-image[x+1+nx*(y+1)+nx*ny*(z)]
						-image[x+nx*(y)+nx*ny*(z-1)]		-image[x-1+nx*(y)+nx*ny*(z-1)]		-image[x+1+nx*(y)+nx*ny*(z-1)]
						-image[x+nx*(y)+nx*ny*(z+1)]		-image[x-1+nx*(y)+nx*ny*(z+1)]		-image[x+1+nx*(y)+nx*ny*(z+1)]
						-image[x+nx*(y-1)+nx*ny*(z-1)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y-1)+nx*ny*(z-1)]
						-image[x+nx*(y-1)+nx*ny*(z+1)]		-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y-1)+nx*ny*(z+1)]
						-image[x+nx*(y+1)+nx*ny*(z-1)]		-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z-1)]
						-image[x+nx*(y+1)+nx*ny*(z+1)]		-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)])/24.0f;		
		} else if (d==1) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x+nx*(y-1)+nx*ny*(z)]		+image[x+nx*(y+1)+nx*ny*(z)])
						-image[x-1+nx*(y)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z)]		-image[x-1+nx*(y+1)+nx*ny*(z)]
						-image[x+1+nx*(y)+nx*ny*(z)]		-image[x+1+nx*(y-1)+nx*ny*(z)]		-image[x+1+nx*(y+1)+nx*ny*(z)]
						-image[x+nx*(y)+nx*ny*(z-1)]		-image[x+nx*(y-1)+nx*ny*(z-1)]		-image[x+nx*(y+1)+nx*ny*(z-1)]
						-image[x+nx*(y)+nx*ny*(z+1)]		-image[x+nx*(y-1)+nx*ny*(z+1)]		-image[x+nx*(y+1)+nx*ny*(z+1)]
						-image[x-1+nx*(y)+nx*ny*(z-1)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-1+nx*(y+1)+nx*ny*(z-1)]
						-image[x-1+nx*(y)+nx*ny*(z+1)]		-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-1+nx*(y+1)+nx*ny*(z+1)]
						-image[x+1+nx*(y)+nx*ny*(z-1)]		-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z-1)]
						-image[x+1+nx*(y)+nx*ny*(z+1)]		-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)])/24.0f;
		} else if (d==2) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x+nx*(y)+nx*ny*(z-1)]		+image[x+nx*(y)+nx*ny*(z+1)])
						-image[x-1+nx*(y)+nx*ny*(z)]		-image[x-1+nx*(y)+nx*ny*(z-1)]		-image[x-1+nx*(y)+nx*ny*(z+1)]
						-image[x+1+nx*(y)+nx*ny*(z)]		-image[x+1+nx*(y)+nx*ny*(z-1)]		-image[x+1+nx*(y)+nx*ny*(z+1)]
						-image[x+nx*(y-1)+nx*ny*(z)]		-image[x+nx*(y-1)+nx*ny*(z-1)]		-image[x+nx*(y-1)+nx*ny*(z+1)]
						-image[x+nx*(y+1)+nx*ny*(z)]		-image[x+nx*(y+1)+nx*ny*(z-1)]		-image[x+nx*(y+1)+nx*ny*(z+1)]
						-image[x-1+nx*(y-1)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-1+nx*(y-1)+nx*ny*(z+1)]
						-image[x-1+nx*(y+1)+nx*ny*(z)]		-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-1+nx*(y+1)+nx*ny*(z+1)]
						-image[x+1+nx*(y-1)+nx*ny*(z)]		-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y-1)+nx*ny*(z+1)]
						-image[x+1+nx*(y+1)+nx*ny*(z)]		-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)])/24.0f;
		} else if (d==3) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y-1)+nx*ny*(z)]		+image[x+1+nx*(y+1)+nx*ny*(z)])
						-image[x-1+nx*(y+1)+nx*ny*(z)]		-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x+1+nx*(y-1)+nx*ny*(z)]		-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x+nx*(y)+nx*ny*(z-1)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z-1)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-2+nx*(y)+nx*ny*(z-1)]		-image[x+nx*(y+2)+nx*ny*(z-1)]
						-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y-2)+nx*ny*(z-1)]		-image[x+2+nx*(y)+nx*ny*(z-1)]
						-image[x+nx*(y)+nx*ny*(z+1)]		-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)]
						-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-2+nx*(y)+nx*ny*(z+1)]		-image[x+nx*(y+2)+nx*ny*(z+1)]
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+nx*(y-2)+nx*ny*(z+1)]		-image[x+2+nx*(y)+nx*ny*(z+1)])/24.0f;
		} else if (d==4) {		
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x+nx*(y-1)+nx*ny*(z-1)]		+image[x+nx*(y+1)+nx*ny*(z+1)])
						-image[x+nx*(y+1)+nx*ny*(z-1)]		-image[x+nx*(y)+nx*ny*(z-2)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x+nx*(y-1)+nx*ny*(z+1)]		-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z+2)]
						-image[x-1+nx*(y)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-1+nx*(y+1)+nx*ny*(z+1)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-1+nx*(y)+nx*ny*(z-2)]		-image[x-1+nx*(y-2)+nx*ny*(z)]
						-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-1+nx*(y-2)+nx*ny*(z)]		-image[x-1+nx*(y)+nx*ny*(z+2)]
						-image[x+1+nx*(y)+nx*ny*(z)]		-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)]
						-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+1+nx*(y)+nx*ny*(z-2)]		-image[x+1+nx*(y)+nx*ny*(z+2)]
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y-2)+nx*ny*(z)]		-image[x+1+nx*(y)+nx*ny*(z+2)])/24.0f;
		} else if (d==5) {	
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y)+nx*ny*(z-1)]		+image[x+1+nx*(y)+nx*ny*(z+1)])
						-image[x+1+nx*(y)+nx*ny*(z-1)]		-image[x+nx*(y)+nx*ny*(z-2)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x-1+nx*(y)+nx*ny*(z+1)]		-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z+2)]
						-image[x+nx*(y-1)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y-1)+nx*ny*(z+1)]
						-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y-1)+nx*ny*(z-2)]		-image[x+2+nx*(y-1)+nx*ny*(z)]
						-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-2+nx*(y-1)+nx*ny*(z)]		-image[x+nx*(y-1)+nx*ny*(z+2)]
						-image[x+nx*(y+1)+nx*ny*(z)]		-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x+1+nx*(y+1)+nx*ny*(z+1)]
						-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+nx*(y+1)+nx*ny*(z-2)]		-image[x+2+nx*(y+1)+nx*ny*(z)]
						-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-2+nx*(y+1)+nx*ny*(z)]		-image[x+nx*(y+1)+nx*ny*(z+2)])/24.0f;
		} else if (d==6) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y+1)+nx*ny*(z)]		+image[x+1+nx*(y-1)+nx*ny*(z)])
						-image[x-1+nx*(y-1)+nx*ny*(z)]		-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y-2)+nx*ny*(z)]
						-image[x+1+nx*(y+1)+nx*ny*(z)]		-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x+nx*(y)+nx*ny*(z-1)]		-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x+1+nx*(y-1)+nx*ny*(z-1)]
						-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-2+nx*(y)+nx*ny*(z-1)]		-image[x+nx*(y-2)+nx*ny*(z-1)]
						-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+nx*(y-2)+nx*ny*(z-1)]		-image[x+2+nx*(y)+nx*ny*(z-1)]
						-image[x+nx*(y)+nx*ny*(z+1)]		-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x+1+nx*(y-1)+nx*ny*(z+1)]
						-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-2+nx*(y)+nx*ny*(z+1)]		-image[x+nx*(y-2)+nx*ny*(z+1)]
						-image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+nx*(y-2)+nx*ny*(z+1)]		-image[x+2+nx*(y)+nx*ny*(z+1)])/24.0f;
		} else if (d==7) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x+nx*(y-1)+nx*ny*(z+1)]		+image[x+nx*(y+1)+nx*ny*(z-1)])
						-image[x+nx*(y-1)+nx*ny*(z-1)]		-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z-2)]
						-image[x+nx*(y+1)+nx*ny*(z+1)]		-image[x+nx*(y)+nx*ny*(z+2)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x-1+nx*(y)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-1+nx*(y+1)+nx*ny*(z-1)]
						-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-1+nx*(y-2)+nx*ny*(z)]		-image[x-1+nx*(y)+nx*ny*(z-2)]
						-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-1+nx*(y)+nx*ny*(z+2)]		-image[x-1+nx*(y+2)+nx*ny*(z)]
						-image[x+1+nx*(y)+nx*ny*(z)]		-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y+1)+nx*ny*(z-1)]
						-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+1+nx*(y-2)+nx*ny*(z)]		-image[x+1+nx*(y)+nx*ny*(z-2)]
						-image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+1+nx*(y)+nx*ny*(z+2)]		-image[x+1+nx*(y+2)+nx*ny*(z)])/24.0f;
		} else if (d==8) {
			val =(8.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y)+nx*ny*(z+1)]		+image[x+1+nx*(y)+nx*ny*(z-1)])
						-image[x-1+nx*(y)+nx*ny*(z-1)]		-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z-2)]
						-image[x+1+nx*(y)+nx*ny*(z+1)]		-image[x+nx*(y)+nx*ny*(z+2)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x+nx*(y-1)+nx*ny*(z)]		-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x+1+nx*(y-1)+nx*ny*(z-1)]
						-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-2+nx*(y-1)+nx*ny*(z)]		-image[x+nx*(y-1)+nx*ny*(z-2)]
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+nx*(y-1)+nx*ny*(z+2)]		-image[x+2+nx*(y-1)+nx*ny*(z)]
						-image[x+nx*(y+1)+nx*ny*(z)]		-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x+1+nx*(y+1)+nx*ny*(z-1)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-2+nx*(y+1)+nx*ny*(z)]		-image[x+nx*(y+1)+nx*ny*(z-2)]
						-image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+nx*(y+1)+nx*ny*(z+2)]		-image[x+2+nx*(y+1)+nx*ny*(z)])/24.0f;
		} else if (d==9) {
			val =(6.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y-1)+nx*ny*(z-1)]	+image[x+1+nx*(y+1)+nx*ny*(z+1)])
						-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-2+nx*(y-2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z+2)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-2+nx*(y)+nx*ny*(z-2)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y-2)+nx*ny*(z-2)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+nx*(y)+nx*ny*(z-2)]		-image[x+2+nx*(y+2)+nx*ny*(z)]
						-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y+2)+nx*ny*(z+2)]
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+2+nx*(y)+nx*ny*(z+2)])/18.0f;
		} else if (d==10) {
			val =(6.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x+1+nx*(y-1)+nx*ny*(z-1)]	+image[x-1+nx*(y+1)+nx*ny*(z+1)])
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+2+nx*(y-2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z+2)]
						-image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+2+nx*(y)+nx*ny*(z-2)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y-2)+nx*ny*(z-2)]		-image[x-2+nx*(y)+nx*ny*(z)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x+nx*(y)+nx*ny*(z-2)]		-image[x-2+nx*(y+2)+nx*ny*(z)]
						-image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x+nx*(y-2)+nx*ny*(z)]		-image[x-2+nx*(y)+nx*ny*(z+2)]
						-image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y+2)+nx*ny*(z+2)])/18.0f;
		} else if (d==11) {				
			val =(6.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y+1)+nx*ny*(z-1)]	+image[x+1+nx*(y-1)+nx*ny*(z+1)])
						 -image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-2+nx*(y+2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z+2)]
						 -image[x+1+nx*(y+1)+nx*ny*(z-1)]	-image[x+nx*(y+2)+nx*ny*(z-2)]		-image[x+2+nx*(y)+nx*ny*(z)]
						 -image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-2+nx*(y)+nx*ny*(z-2)]		-image[x+nx*(y-2)+nx*ny*(z)]
						 -image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y)+nx*ny*(z-2)]		-image[x+2+nx*(y-2)+nx*ny*(z)]
						 -image[x-1+nx*(y-1)+nx*ny*(z+1)]	-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y-2)+nx*ny*(z+2)]
						 -image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+nx*(y+2)+nx*ny*(z)]		-image[x+2+nx*(y)+nx*ny*(z+2)])/18.0f;
		} else if (d==12) {				
			val =(6.0f*(image[x+nx*(y)+nx*ny*(z)]			+image[x-1+nx*(y-1)+nx*ny*(z+1)]	+image[x+1+nx*(y+1)+nx*ny*(z-1)])
						-image[x-1+nx*(y+1)+nx*ny*(z+1)]	-image[x-2+nx*(y)+nx*ny*(z-2)]		-image[x+nx*(y+2)+nx*ny*(z)]
						-image[x+1+nx*(y-1)+nx*ny*(z+1)]	-image[x+nx*(y-2)+nx*ny*(z+2)]		-image[x+2+nx*(y)+nx*ny*(z)]
						-image[x-1+nx*(y-1)+nx*ny*(z-1)]	-image[x-2+nx*(y-2)+nx*ny*(z)]		-image[x+nx*(y)+nx*ny*(z-2)]
						-image[x+1+nx*(y-1)+nx*ny*(z-1)]	-image[x+nx*(y-2)+nx*ny*(z)]		-image[x+2+nx*(y)+nx*ny*(z-2)]
						-image[x-1+nx*(y+1)+nx*ny*(z-1)]	-image[x-2+nx*(y)+nx*ny*(z)]		-image[x+nx*(y+2)+nx*ny*(z-2)]
						-image[x+1+nx*(y+1)+nx*ny*(z+1)]	-image[x+nx*(y)+nx*ny*(z+2)]		-image[x+2+nx*(y+2)+nx*ny*(z)])/18.0f;
		}
		return 0.5f*val;
	}
	
	private final void directionFromTubularFilter(float[] img, boolean[] mask, float[] filter,byte[] direction) {
			
			// get the tubular filter response
			float[][][] image = new float[nx][ny][nz];
			//float[] filter = new float[nx*ny*nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				image[x][y][z] = img[xyz];
			}
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				filter[xyz] = 0.0f;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					// check for zero-valued neighbors as well
					float best = 0.0f;
					int dmax = 13;
					byte dbest = -1;
					for (byte d=0;d<dmax;d++) {
						float val = tubularScore(img, x,y,z,d);
						if (val*val>best*best) {
							best = val;
							dbest = d;
						}
					}
					filter[xyz] = best;
					if (filter[xyz]<0) {
						filter[xyz]=0;
					}
					direction[xyz] = dbest;
				}
			}
			image = null;
			return;
	}

}
