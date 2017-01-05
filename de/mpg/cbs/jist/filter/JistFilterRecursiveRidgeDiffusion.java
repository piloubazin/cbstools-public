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
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
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
	
	private ParamVolume locationImage;
	
	private ParamFloat scalingParam;
	private ParamInteger nscalesParam;
	
	//private ParamFloat thresholdParam;
	
	private ParamOption filterParam;
	private static final String[] filterTypes = {"2D","1D"};
	
	private ParamOption propagationParam;
	private static final String[] propagationTypes = {"diffusion","belief","SUR"};
	private ParamFloat factorParam;
	private ParamFloat diffParam;
	private ParamInteger ngbParam;
	private ParamInteger iterParam;
	private ParamFloat maxdiffParam;
	
	private ParamVolume filterImage;
	private ParamVolume probaImage;
	private ParamVolume scaleImage;
	private ParamVolume directionImage;
	private ParamVolume correctImage;
	
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
		inputParams.add(scalingParam = new ParamFloat("Angular factor", 0.0f, 10.0f, 1.0f));
		
		inputParams.add(locationImage = new ParamVolume("Location prior (opt)"));
		locationImage.setMandatory(false);
		
		//inputParams.add(thresholdParam = new ParamFloat("Probability threshold", 0.0, 1.0, 0.5));
		//inputParams.add(scalingParam = new ParamFloat("Scale factor ", 0.25f, 2.0f, 1.0f));
		inputParams.add(nscalesParam = new ParamInteger("Number of scales ", 0, 10, 3));

		inputParams.add(propagationParam = new ParamOption("Propagation model", propagationTypes));
		inputParams.add(factorParam = new ParamFloat("Diffusion factor", 0.0f, 10.0f, 1.0f));
		inputParams.add(diffParam = new ParamFloat("Similarity scale", 0.0f, 1.0f, 0.1f));
		inputParams.add(ngbParam = new ParamInteger("Neighborhood size", 0, 26, 4));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000, 100));
		inputParams.add(maxdiffParam = new ParamFloat("Max difference", 0.0f, 1.0f, 0.001f));
		
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
		outputParams.add(scaleImage = new ParamVolume("Detection scale",VoxelType.BYTE,-1,-1,-1,-1));
		outputParams.add(directionImage = new ParamVolume("Ridge direction",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(correctImage = new ParamVolume("Directional correction",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("Recursive Ridge Diffusion");
		outputParams.setLabel("Recursive Ridge Diffusion");
	}
	
	@Override
	protected void execute(CalculationMonitor monitor){
		BasicInfo.displayMessage("recursive ridge diffusion:\n");
		
		// import the image data into 1D arrays
		float[] image = Interface.getFloatImage3D(inputImage);
		int[] dim = Interface.getDimensions(inputImage);
		nx = dim[X]; ny = dim[Y]; nz = dim[Z];
		nxyz = nx*ny*nz;
		BasicInfo.displayMessage("...load data\n");
		
		// load surface, if used
		float[] surface = null;
		nc = 0;
		String orientation = orientationParam.getValue();
		if (orientation.equals("parallel") || orientation.equals("orthogonal")) {
			BasicInfo.displayMessage("...load surface orientation(s)\n");
			nc = Interface.getComponents(surfaceImage);
			surface = Interface.getFloatImage4D(surfaceImage);
		}
		
		// load masking, if used
		float[] location = null;
		if (Interface.isValid(locationImage)) {
			location = Interface.getFloatImage3D(locationImage);
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
		
		BasicInfo.displayMessage("...first filter response\n");
		
		if (filterParam.getValue().equals("1D")) directionFromRecursiveRidgeFilter1D(image, mask, maxresponse, maxdirection);
		else if (filterParam.getValue().equals("2D")) directionFromRecursiveRidgeFilter2D(image, mask, maxresponse, maxdirection);
			
		for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
			if (maxresponse[xyz]>0) {
				maxscale[xyz] = (byte)0;
			} else {
				maxscale[xyz] = (byte)-1;
			}
		}
		for (int i=0;i<nscales;i++) {
			float scale = 1.0f+i;
			float[] smoothed = new float[nxyz];

			BasicInfo.displayMessage("...filter response at scale "+scale+"\n");
		
			// Gaussian Kernel
			float[][] G = ImageFilters.separableGaussianKernel(scale/L2N2,scale/L2N2,scale/L2N2);
			int gx = (G[0].length-1)/2;
			int gy = (G[1].length-1)/2;
			int gz = (G[2].length-1)/2;
				
			// smoothed image
			smoothed = ImageFilters.separableConvolution(image,nx,ny,nz,G,gx,gy,gz); 

			byte[] direction = new byte[nxyz];
			float[] response = new float[nxyz];
			if (filterParam.getValue().equals("1D")) directionFromRecursiveRidgeFilter1D(smoothed, mask, response, direction);
			else if (filterParam.getValue().equals("2D")) directionFromRecursiveRidgeFilter2D(smoothed, mask, response, direction);
			
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
		float[] correction = new float[nxyz];
		if (orientation.equals("parallel") || orientation.equals("orthogonal")) {
			BasicInfo.displayMessage("...orientation response filtering\n");
			float factor = scalingParam.getValue().floatValue();
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float[] vec = directionVector(maxdirection[xyz]);
				// find the closest surfaces; weighted mean
				int closest = 0;
				float dist = Numerics.abs(surface[xyz]);
				for (int c=1;c<nc;c++) {
					if (Numerics.abs(surface[xyz+c*nxyz])<dist) {
						closest = c;
						dist = Numerics.abs(surface[xyz+c*nxyz]);
					}
				}
				int second = closest;
				if (closest<nc-1 && surface[xyz+closest*nxyz]*surface[xyz+(closest+1)*nxyz]<0) second = closest+1; 
				if (closest>0 && surface[xyz+closest*nxyz]*surface[xyz+(closest-1)*nxyz]>0) second = closest-1; 
				
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
				// angle projected back into [0,1], 0: parallel, 1:orthogonal
				double angle = 2.0*FastMath.acos(grad[X]*vec[X]+grad[Y]*vec[Y]+grad[Z]*vec[Z])/FastMath.PI;
				angle = Numerics.min(angle,2.0-angle);
				if (orientation.equals("parallel") && filterParam.getValue().equals("1D")) {
					correction[xyz] = (float)FastMath.pow(angle, factor);
				} else if (orientation.equals("parallel") && filterParam.getValue().equals("2D")) {
					correction[xyz] = (float)FastMath.pow(1.0-angle, factor);
				} else if (orientation.equals("orthogonal") && filterParam.getValue().equals("1D")) {
					correction[xyz] = (float)FastMath.pow(1.0-angle, factor);
				} else if (orientation.equals("orthogonal") && filterParam.getValue().equals("2D")) {
					correction[xyz] = (float)FastMath.pow(angle, factor);
				} else {
					correction[xyz] = 1.0f;
				}
				maxresponse[xyz] *= correction[xyz];
			}
		}
		
		// mask response with location prior before equalization
		if (location!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) maxresponse[xyz] *= location[xyz];
		}

		// Equalize histogram (Exp Median)
		BasicInfo.displayMessage("...normalization into probabilities\n");
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
		
		// 3. diffuse the data to neighboring structures with SUR
		BasicInfo.displayMessage("...diffusion\n");
		int ngbsize = ngbParam.getValue().intValue();
		float scale = diffParam.getValue().floatValue();
		float factor = factorParam.getValue().floatValue();
		int iter = iterParam.getValue().intValue();
		float maxdiff = maxdiffParam.getValue().floatValue();
		
		if (propagationParam.getValue().equals("SUR")) {
			if (filterParam.getValue().equals("1D"))  spatialUncertaintyRelaxation1D(proba, maxdirection, ngbsize, scale, factor, iter);
			else if (filterParam.getValue().equals("2D"))  spatialUncertaintyRelaxation2D(proba, maxdirection, ngbsize, scale, factor, iter);
		} else if  (propagationParam.getValue().equals("diffusion")) {
			if (filterParam.getValue().equals("1D"))  probabilisticDiffusion1D(proba, maxdirection, ngbsize, maxdiff, scale, factor, iter);
			else if (filterParam.getValue().equals("2D"))  probabilisticDiffusion2D(proba, maxdirection, ngbsize, maxdiff, scale, factor, iter);
		} else if  (propagationParam.getValue().equals("belief")) {
			if (filterParam.getValue().equals("1D"))  beliefPropagation1D(proba, maxdirection, ngbsize, scale, iter, maxdiff);
			else if (filterParam.getValue().equals("2D"))  beliefPropagation2D(proba, maxdirection, ngbsize, scale, iter, maxdiff);
		} 				
		// Output
		BasicInfo.displayMessage("...output images\n");
		String outname = Interface.getName(inputImage);
		ImageHeader outheader = Interface.getHeader(inputImage);
		
		if (filterParam.getValue().equals("1D")) outname+="_rrf1d";
		else outname+="_rrf2d";
		
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = maxresponse[id];
		}					
		ImageDataFloat resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_filter");
		filterImage.setValue(resData);
		
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = proba[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_proba");
		probaImage.setValue(resData);
		
		byte[][][] resultb = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			resultb[x][y][z] = maxscale[id];
		}					
		ImageDataByte resDataB = new ImageDataByte(resultb);	
		resDataB.setHeader(outheader);
		resDataB.setName(outname+"_scale");
		scaleImage.setValue(resDataB);
		
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
		
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = correction[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(outheader);
		resData.setName(outname+"_correct");
		correctImage.setValue(resData);
		
		return;
	}
	
	private final void directionFromRecursiveRidgeFilter1D(float[] img, boolean[] mask, float[] filter,byte[] direction) {
			
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
	private final void directionFromRecursiveRidgeFilter2D(float[] img, boolean[] mask, float[] filter,byte[] direction) {
			
			// get the tubular filter response
			float[][][] planescore = new float[nx][ny][nz];
			byte[][][] planedir = new byte[nx][ny][nz];
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
					filter[xyz] = planescore[x][y][z];
					direction[xyz] = planedir[x][y][z];
					if (filter[xyz]<0) { filter[xyz]=0; direction[xyz] = -1; }
				}
			}
			planescore = null;
			planedir = null;
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

	private final void spatialUncertaintyRelaxation1D(float[] proba, byte[] dir, int ngbsize, float scale, float factor, int iter) {
		// run SUR with weighting against normal directions
		//mix with the neighbors?
		float[] certainty = new float[nx*ny*nz];   	
		float[] newproba = new float[nx*ny*nz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nx*ny*nz][26];  
		float[][] angleweight1D = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			angleweight1D[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		
		// SOR-scheme?
		//float sorfactor = 1.95f;
		float sorfactor = 1.0f;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyzi = x+nx*y+nx*ny*z;
			for (byte dj=0;dj<26;dj++) {
				int xyzj = Ngb.neighborIndex(dj, xyzi, nx, ny, nz);
				imgweight[xyzi][dj] = sorfactor*diffusionWeightFunction(proba, dir, xyzi, xyzj, dj, angleweight1D, scale)/ngbsize;
			}
		}
		for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
			newproba[xyzi] = proba[xyzi];
		}
		
		float maxdiff = 1.0f;
		float meandiff = 1.0f;
		float t0diff = 1.0f;
		for (int t=0;t<iter && meandiff>0.001f*t0diff;t++) {
		//for (int t=0;t<iter && maxdiff>0.1f;t++) {
			BasicInfo.displayMessage("iter "+(t+1));
			maxdiff = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			float nflip = 0.0f;
			float nproc = 0.0f;
							
			// mask: remove regions 0,1 away from ]0,1[ regions
			boolean[] diffmask = new boolean[nxyz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (proba[xyz]>0 && proba[xyz]<1) diffmask[xyz] = true;
				else diffmask[xyz] = false;
				for (int i=-1;i<=1 && !diffmask[xyz];i++) for (int j=-1;j<=1 && !diffmask[xyz];j++) for (int k=-1;k<=1 && !diffmask[xyz];k++) {
					if (proba[xyz+i+j*nx+k*nx*ny]>0 && proba[xyz+i+j*nx+k*nx*ny]<1) diffmask[xyz] = true;
				}
			}
			
			// main loop: label-per-label
			//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			//BasicInfo.displayMessage(".");

			// get the gain ; normalize
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				certainty[xyzi] = Numerics.abs(2.0f*proba[xyzi]-1.0f);
			}
			// propagate the values : diffusion
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				nproc++;
						
				float den = certainty[xyzi];
				float num = den*proba[xyzi];
				float prev = proba[xyzi];
					
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
					ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
				}
				byte[] rank = Numerics.argmax(ngbweight, ngbsize);
				for (int l=0;l<ngbsize;l++) {
					int xyzl = Ngb.neighborIndex(rank[l], xyzi, nx, ny, nz);
					num += ngbweight[rank[l]]*proba[xyzl];
					den += ngbweight[rank[l]];
				}
				if (den>1e-9f) num /= den;
					
				newproba[xyzi] = num;
						
				meandiff += Numerics.abs(num-prev);
				ndiff++;
				maxdiff = Numerics.max(maxdiff, Numerics.abs(num-prev));
				if (prev<0.5f && num>0.5f) nflip++;
				if (prev>0.5f && num<0.5f) nflip++;
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				proba[xyzi] = newproba[xyzi];
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiff+"\n");
			//BasicInfo.displayMessage("n processed "+nproc+", n flipped "+nflip+"\n");
			
			//BasicInfo.displayMessage("n resorted"+nresort+"\n");
			
		}
		newproba = null;
	}
	
	private final void spatialUncertaintyRelaxation2D(float[] proba, byte[] dir, int ngbsize, float scale, float factor, int iter) {
		// run SUR with weighting against normal directions
		//mix with the neighbors?
		float[] certainty = new float[nx*ny*nz];   	
		float[] newproba = new float[nx*ny*nz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nx*ny*nz][26];  
		float[][] angleweight2D = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			angleweight2D[d1][d2] = (float)FastMath.pow(2.0f*FastMath.acos(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		
		// SOR-scheme?
		//float sorfactor = 1.95f;
		float sorfactor = 1.0f;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyzi = x+nx*y+nx*ny*z;
			for (byte dj=0;dj<26;dj++) {
				int xyzj = Ngb.neighborIndex(dj, xyzi, nx, ny, nz);
				imgweight[xyzi][dj] = sorfactor*diffusionWeightFunction(proba, dir, xyzi, xyzj, dj, angleweight2D, scale)/ngbsize;
			}
		}
		for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
			newproba[xyzi] = proba[xyzi];
		}
		
		float maxdiff = 1.0f;
		float meandiff = 1.0f;
		float t0diff = 1.0f;
		for (int t=0;t<iter && meandiff>0.001f*t0diff;t++) {
		//for (int t=0;t<iter && maxdiff>0.1f;t++) {
			BasicInfo.displayMessage("iter "+(t+1));
			maxdiff = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			float nflip = 0.0f;
			float nproc = 0.0f;
							
			// mask: remove regions 0,1 away from ]0,1[ regions
			boolean[] diffmask = new boolean[nxyz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (proba[xyz]>0 && proba[xyz]<1) diffmask[xyz] = true;
				else diffmask[xyz] = false;
				for (int i=-1;i<=1 && !diffmask[xyz];i++) for (int j=-1;j<=1 && !diffmask[xyz];j++) for (int k=-1;k<=1 && !diffmask[xyz];k++) {
					if (proba[xyz+i+j*nx+k*nx*ny]>0 && proba[xyz+i+j*nx+k*nx*ny]<1) diffmask[xyz] = true;
				}
			}
			
			// main loop: label-per-label
			//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			//BasicInfo.displayMessage(".");

			// get the gain ; normalize
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				certainty[xyzi] = Numerics.abs(2.0f*proba[xyzi]-1.0f);
			}
			// propagate the values : diffusion
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				nproc++;
						
				float den = certainty[xyzi];
				float num = den*proba[xyzi];
				float prev = proba[xyzi];
					
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
					ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
				}
				byte[] rank = Numerics.argmax(ngbweight, ngbsize);
				for (int l=0;l<ngbsize;l++) {
					int xyzl = Ngb.neighborIndex(rank[l], xyzi, nx, ny, nz);
					num += ngbweight[rank[l]]*proba[xyzl];
					den += ngbweight[rank[l]];
				}
				if (den>1e-9f) num /= den;
					
				newproba[xyzi] = num;
						
				meandiff += Numerics.abs(num-prev);
				ndiff++;
				maxdiff = Numerics.max(maxdiff, Numerics.abs(num-prev));
				if (prev<0.5f && num>0.5f) nflip++;
				if (prev>0.5f && num<0.5f) nflip++;
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (diffmask[xyzi]) {
				proba[xyzi] = newproba[xyzi];
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiff+"\n");
			//BasicInfo.displayMessage("n processed "+nproc+", n flipped "+nflip+"\n");
			
			//BasicInfo.displayMessage("n resorted"+nresort+"\n");
			
		}
		newproba = null;
	}
	
	private final float diffusionWeightFunction(float[] proba, byte[] dir, int xyz, int ngb, byte dngb, float[][] corr, float scale) {
		float diff = Numerics.abs((proba[xyz] - proba[ngb])/scale);
		float weight = 1.0f;
		if (dir[xyz]>-1 && dngb>-1) weight = corr[dir[xyz]][dngb];
		return weight/(1.0f+Numerics.square(diff/scale));
	}
	
	
	private final float[] probabilisticDiffusion1D(float[] proba, byte[] dir, int ngbsize, float maxdiff, float angle, float factor, int iter) {
		// build a similarity function from aligned neighbors
		float[][] similarity = null;
		byte[][] neighbor = null;
		
		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbsize, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}

		for (int t=0;t<iter;t++) {
			Interface.displayMessage("iteration "+(t+1)+": ");
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				float prev = diffused[xyz];
				diffused[xyz] = proba[xyz];
				// weight with distance to 0 or 1
				for (int n=0;n<ngbsize;n++) {
					float ngb = diffused[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				
				if (Numerics.abs(prev-diffused[xyz])>diff) diff = Numerics.abs(prev-diffused[xyz]);
			}
			
			Interface.displayMessage("diff "+diff+"\n");
			if (diff<maxdiff) t=iter;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iter==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		return diffused;
	}

    private final void estimateSimpleDiffusionSimilarity1D(byte[] dir, float[] proba, int ngbsize, byte[][] neighbor, float[][] similarity, float factor) {
    	
    	// initialize
    	neighbor = new byte[ngbsize][nxyz];
    	similarity = new float[ngbsize][nxyz];
    	
    	float[][] angleweight1D = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			angleweight1D[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		float[] weight = new float[26];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (proba[id]>0) {
				// find the N best aligned discrete directions
				for (byte d=0;d<NC2;d++) {
					int idn = neighborIndex(d,id);
					
					if (proba[idn]>0) {
						if (ngbsize==2) weight[d] = angleweight1D[d][dir[id]];
						else weight[d] = angleweight1D[d][dir[id]]*proba[idn];
					} else {
						weight[d] = 0.0f;
					}
				}
				byte[] ngb = Numerics.argmax(weight, ngbsize);
				for (int n=0;n<ngbsize;n++) {
					neighbor[n][id] = ngb[n];
					int idn = neighborIndex(ngb[n],id);
					if (proba[idn]>0) {
						similarity[n][id] = angleweight1D[dir[id]][dir[idn]];
					} else {
						similarity[n][id] = 0.0f;
					}
				}
			}
		}
		
		return;
    }
    
	private final float[] probabilisticDiffusion2D(float[] proba, byte[] dir, int ngbsize, float maxdiff, float angle, float factor, int iter) {
		// build a similarity function from aligned neighbors
		float[][] similarity = null;
		byte[][] neighbor = null;
		
		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbsize, neighbor, similarity, angle);
		
		// run the diffusion process
		float[] diffused = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}

		for (int t=0;t<iter;t++) {
			Interface.displayMessage("iteration "+(t+1)+": ");
			float diff = 0.0f;
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0) {
				float prev = diffused[xyz];
				diffused[xyz] = proba[xyz];
				// weight with distance to 0 or 1
				for (int n=0;n<ngbsize;n++) {
					float ngb = diffused[neighborIndex(neighbor[n][xyz],xyz)];
					// remap neighbors?
					ngb = (float)FastMath.exp(ngb)-1.0f;
				
					// integration over the whole vessel (log version is more stable (??) )
					diffused[xyz] += factor*similarity[n][xyz]*ngb;
				}
				diffused[xyz] = (float)FastMath.log(1.0f + diffused[xyz]);
				
				if (Numerics.abs(prev-diffused[xyz])>diff) diff = Numerics.abs(prev-diffused[xyz]);
			}
			
			Interface.displayMessage("diff "+diff+"\n");
			if (diff<maxdiff) t=iter;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (iter==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0)/(1.0f+2.0f*factor), 0.0f, 1.0f);
		}
		return diffused;
	}

    private final void estimateSimpleDiffusionSimilarity2D(byte[] dir, float[] proba, int ngbsize, byte[][] neighbor, float[][] similarity, float factor) {
    	
    	// initialize
    	neighbor = new byte[ngbsize][nxyz];
    	similarity = new float[ngbsize][nxyz];
    	
    	float[][] angleweight1D = new float[26][26];
		float[][] angleweight2D = new float[26][26];
		for (int d1=0;d1<26;d1++) for (int d2=0;d2<26;d2++) {
			float[] dir1 = directionVector(d1);
			float[] dir2 = directionVector(d2);
			angleweight1D[d1][d2] = (float)FastMath.pow(2.0f*FastMath.asin(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
			angleweight2D[d1][d2] = (float)FastMath.pow(2.0f*FastMath.acos(Numerics.abs(dir1[X]*dir2[X] + dir1[Y]*dir2[Y] + dir1[Z]*dir2[Z]))/FastMath.PI,factor);
		}
		float[] weight = new float[26];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (proba[id]>0) {
				// find the N best planar discrete directions
				for (byte d=0;d<NC2;d++) {
					int idn = neighborIndex(d,id);
					
					if (proba[idn]>0) {
						if (ngbsize==2) weight[d] = angleweight2D[d][dir[id]];
						else weight[d] = angleweight2D[d][dir[id]]*proba[idn];
					} else {
						weight[d] = 0.0f;
					}
				}
				byte[] ngb = Numerics.argmax(weight, ngbsize);
				for (int n=0;n<ngbsize;n++) {
					neighbor[n][id] = ngb[n];
					int idn = neighborIndex(ngb[n],id);
					if (proba[idn]>0) {
						// similarity comes from the normal direction
						similarity[n][id] = angleweight1D[dir[id]][dir[idn]];
					} else {
						similarity[n][id] = 0.0f;
					}
				}
			}
		}
		
		return;
    }
    
    private final float[] beliefPropagation1D(float[] proba, byte[] dir, int ngbsize, float angle, int iter, float maxdiff) {
		// build a similarity function from aligned neighbors
		float[][] similarity = null;
		byte[][] neighbor = null;
		
		estimateSimpleDiffusionSimilarity1D(dir, proba, ngbsize, neighbor, similarity, angle);
		
		float[][] messagefg = new float[ngbsize][nxyz];
		float[][] messagebg = new float[ngbsize][nxyz];
		
		float[][] newmsgfg = new float[ngbsize][nxyz];
		float[][] newmsgbg = new float[ngbsize][nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) for (int n=0;n<ngbsize;n++) {
			messagefg[n][xyz] = 1.0f;
			messagebg[n][xyz] = 1.0f;
			newmsgfg[n][xyz] = 1.0f;
			newmsgbg[n][xyz] = 1.0f;
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
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					float simfgfg = similarity[n][xyz];
					float simbgfg = 0.0f; //(1.0f-similarity1[xyz]); //similarity1[xyz];
					float simfgbg = (1.0f-similarity[n][xyz]); // good results
					float simbgbg = 1.0f; //(1.0f-similarity1[xyz]); //1.0f; //similarity1[xyz];
				
					ngb = neighborIndex(neighbor[n][xyz],xyz);
					// neighbor message, fg label
					float mfgfg = Numerics.max(mins,simfgfg)*Numerics.max(minp,proba[ngb]);
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mfgfg *= messagefg[m][ngb];
					}
					
					float mfgbg = Numerics.max(mins,simfgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mfgbg *= messagebg[m][ngb];
					}
					newmsgfg[n][xyz] = Numerics.max(mfgfg, mfgbg);
					diff = Numerics.max(diff, Numerics.abs(newmsgfg[n][xyz]-messagefg[n][xyz]));
				
					// neighbor message, bg label
					float mbgbg = Numerics.max(mins,simbgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mbgbg *= messagebg[m][ngb];
					}
				
					float mbgfg = Numerics.max(mins,simbgfg)*Numerics.max(minp,proba[ngb]);
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mbgfg *= messagefg[m][ngb];
					}
					newmsgbg[n][xyz] = Numerics.max(mbgbg, mbgfg);
					diff = Numerics.max(diff, Numerics.abs(newmsgbg[n][xyz]-messagebg[n][xyz]));
				}
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					messagefg[n][xyz] = newmsgfg[n][xyz];
					messagebg[n][xyz] = newmsgbg[n][xyz];
				}
			}
			Interface.displayMessage("diff "+diff+"\n");
			if (diff < maxdiff) t = iter;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz];
			
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					if (messagefg[n][xyz]>0) belief[xyz] *= messagefg[n][xyz];
				}
			}
			
			bgbelief = (1.0f-proba[xyz]);
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					if (messagebg[n][xyz]>0) bgbelief *= messagebg[n][xyz];
				}
			}
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
			belief[xyz] = belief[xyz]/(belief[xyz]+bgbelief);
		}
		return belief;
	}
	
    private final float[] beliefPropagation2D(float[] proba, byte[] dir, int ngbsize, float angle, int iter, float maxdiff) {
		// build a similarity function from aligned neighbors
		float[][] similarity = null;
		byte[][] neighbor = null;
		
		estimateSimpleDiffusionSimilarity2D(dir, proba, ngbsize, neighbor, similarity, angle);
		
		float[][] messagefg = new float[ngbsize][nxyz];
		float[][] messagebg = new float[ngbsize][nxyz];
		
		float[][] newmsgfg = new float[ngbsize][nxyz];
		float[][] newmsgbg = new float[ngbsize][nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) for (int n=0;n<ngbsize;n++) {
			messagefg[n][xyz] = 1.0f;
			messagebg[n][xyz] = 1.0f;
			newmsgfg[n][xyz] = 1.0f;
			newmsgbg[n][xyz] = 1.0f;
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
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					float simfgfg = similarity[n][xyz];
					float simbgfg = 0.0f; //(1.0f-similarity1[xyz]); //similarity1[xyz];
					float simfgbg = (1.0f-similarity[n][xyz]); // good results
					float simbgbg = 1.0f; //(1.0f-similarity1[xyz]); //1.0f; //similarity1[xyz];
				
					ngb = neighborIndex(neighbor[n][xyz],xyz);
					// neighbor message, fg label
					float mfgfg = Numerics.max(mins,simfgfg)*Numerics.max(minp,proba[ngb]);
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mfgfg *= messagefg[m][ngb];
					}
					
					float mfgbg = Numerics.max(mins,simfgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mfgbg *= messagebg[m][ngb];
					}
					newmsgfg[n][xyz] = Numerics.max(mfgfg, mfgbg);
					diff = Numerics.max(diff, Numerics.abs(newmsgfg[n][xyz]-messagefg[n][xyz]));
				
					// neighbor message, bg label
					float mbgbg = Numerics.max(mins,simbgbg)*Numerics.max(minp,(1.0f-proba[ngb]));
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagebg[m][ngb]>0) mbgbg *= messagebg[m][ngb];
					}
				
					float mbgfg = Numerics.max(mins,simbgfg)*Numerics.max(minp,proba[ngb]);
					for (int m=0;m<ngbsize;m++) {
						if (neighborIndex(neighbor[m][ngb],ngb)!=xyz && messagefg[m][ngb]>0) mbgfg *= messagefg[m][ngb];
					}
					newmsgbg[n][xyz] = Numerics.max(mbgbg, mbgfg);
					diff = Numerics.max(diff, Numerics.abs(newmsgbg[n][xyz]-messagebg[n][xyz]));
				}
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					messagefg[n][xyz] = newmsgfg[n][xyz];
					messagebg[n][xyz] = newmsgbg[n][xyz];
				}
			}
			Interface.displayMessage("diff "+diff+"\n");
			if (diff < maxdiff) t = iter;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz];
			
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					if (messagefg[n][xyz]>0) belief[xyz] *= messagefg[n][xyz];
				}
			}
			
			bgbelief = (1.0f-proba[xyz]);
			if (proba[xyz]>0 && proba[xyz]<1) {
				for (int n=0;n<ngbsize;n++) {
					if (messagebg[n][xyz]>0) bgbelief *= messagebg[n][xyz];
				}
			}
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
			belief[xyz] = belief[xyz]/(belief[xyz]+bgbelief);
		}
		return belief;
	}
	

}
