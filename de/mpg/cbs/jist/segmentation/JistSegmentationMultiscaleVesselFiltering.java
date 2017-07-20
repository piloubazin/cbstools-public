package de.mpg.cbs.jist.segmentation;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
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

import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import    org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.*;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

import java.io.*;

public class JistSegmentationMultiscaleVesselFiltering extends ProcessingAlgorithm {
	
	// jist containers
	private ParamVolume inputImage;
	private ParamOption brightParam;
	private static final String[] brightTypes = {"bright","dark"};
	
	private ParamFloat thresholdParam;
	
	private ParamOption filterParam;
	private static final String[] filterTypes = {"RRF","Hessian"};
	
	private ParamOption propagationParam;
	private static final String[] propagationTypes = {"diffusion","belief"};
	private ParamFloat factorParam;
	private ParamFloat diffParam;
	private ParamInteger iterParam;
	
	private ParamVolume vesselImage;
	private ParamVolume filterImage;
	private ParamVolume probaImage;
	private ParamVolume scaleImage;
	private ParamVolume diameterImage;
	private ParamVolume lengthImage;
	private ParamVolume pvImage;
	
	private ParamFloat scaleStepParam;
	private ParamInteger nbrScaleParam;
	
	private ParamVolume scaleMap;
	
	/*
	private ParamVolume valIn;
	private ParamVolume valOut;
	private ParamVolume probaNew;
	
	private ParamFloat thresholdParamProba;
	private ParamVolume probaSeg;
	
	private ParamVolume diamRange0;
	private ParamVolume diamRange;
	
	private ParamVolume diam;	
		
	private ParamFloat thresholdParamPv;
	private ParamVolume pvEst;
	private ParamVolume largeSeg;
	
	private ParamInteger thresholdReg;
	private ParamVolume cleanPvEst;
	private ParamVolume cleanLargeSeg;
	*/
	
	
	//private ParamFile 	outputParam;
	// global variables
	int nx, ny, nz, nc, nxyz;
	int scaleNbr;
	
	// numerical quantities
	private static final	float	INVSQRT2 = (float)(1.0f/FastMath.sqrt(2.0f));
	private static final	float	INVSQRT3 = (float)(1.0f/FastMath.sqrt(3.0f));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0f);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0f);
	private static final	float	SQRT2PI = (float)FastMath.sqrt(2.0f*(float)Math.PI);
	private static final	float	PI2 = (float)(Math.PI/2.0f);
	private static final	float   	L2N2=2.0f*(float)(FastMath.sqrt(2.0f*(float)(FastMath.log(2.0f))));
	
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
		inputParams.add(brightParam = new ParamOption("Vessel intensities", brightTypes));
		inputParams.add(filterParam = new ParamOption("Shape filter", filterTypes));
		inputParams.add(thresholdParam = new ParamFloat("Probability threshold", 0.0f, 1.0f, 0.5f));
		//inputParams.add(thresholdParamProba = new ParamFloat("In Probability threshold", 0.0f, 1.0f, 0.5f));
		//inputParams.add(thresholdReg = new ParamInteger("Clean Estimation threshold", 0, 50, 15));
		//inputParams.add(thresholdParamPv = new ParamFloat("Partial Voluming threshold", 0.0f, 1.0f, 0.7));
		inputParams.add(scaleStepParam = new ParamFloat("Scale Step ", 0.25f, 2.0f, 1.0f));
		inputParams.add(nbrScaleParam = new ParamInteger("Scale number ", 1, 10, 4));

		inputParams.add(propagationParam = new ParamOption("Propagation model", propagationTypes));
		inputParams.add(factorParam = new ParamFloat("Diffusion factor", 0.0f, 100.0f, 0.5f));
		inputParams.add(diffParam = new ParamFloat("Max difference", 0.0f, 1.0f, 0.001f));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000, 100));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation.devel");
		inputParams.setLabel("Multiscale Vessel Filtering");
		inputParams.setName("MultiscaleVesselFiltering");
		
		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Victoire Plessis", "plessis@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Segments vasculature at multiple scales.");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}
	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(vesselImage = new ParamVolume("Segmented Vessels",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(filterImage = new ParamVolume("Filtered Vessels",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(probaImage = new ParamVolume("Vessel Probability",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(scaleImage = new ParamVolume("Vessel Scale",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(diameterImage = new ParamVolume("Vessel Diameter",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(pvImage = new ParamVolume("Vessel Partial Voluming",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(lengthImage = new ParamVolume("Vessel Length",VoxelType.FLOAT,-1,-1,-1,-1));
		/*
		outputParams.add(largeSeg = new ParamVolume("large Segmentation",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(cleanLargeSeg = new ParamVolume("Clean large Segmentation",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(probaSeg = new ParamVolume("proba Segmentation",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(labelImage = new ParamVolume("Labeled Vessels",VoxelType.INT,-1,-1,-1,-1));
	
		outputParams.add(probaNew = new ParamVolume("Vessel In Probability",VoxelType.FLOAT,-1,-1,-1,-1));
		if (debug) {
			outputParams.add(scaleMap = new ParamVolume("Scale Map",VoxelType.FLOAT,-1,-1,-1,-1));
			outputParams.add(valIn = new ParamVolume("In values",VoxelType.FLOAT,-1,-1,-1,-1));
			outputParams.add(valOut = new ParamVolume("Out values",VoxelType.FLOAT,-1,-1,-1,-1));
		}	
		outputParams.add(diamRange0 = new ParamVolume("Diameter estime",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(diamRange = new ParamVolume("Diameter estime 2",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(diam = new ParamVolume("Diameter",VoxelType.DOUBLE,-1,-1,-1,-1));
		outputParams.add(pvEst = new ParamVolume("Partial Voluming Estimation",VoxelType.DOUBLE,-1,-1,-1,-1));
		outputParams.add(cleanPvEst = new ParamVolume("CleanPartial Voluming Estimation",VoxelType.DOUBLE,-1,-1,-1,-1));
		*/
		
		outputParams.setName("Multiscale vessel filtering");
		outputParams.setLabel("Multiscale vessel filtering");
	}
	
	@Override
	protected void execute(CalculationMonitor monitor){
		// load scale values
		scaleNbr=nbrScaleParam.getValue().intValue();
		float scaleStep=scaleStepParam.getValue().floatValue();
		float scaleFactors[]=new float[scaleNbr];
		for(int i=0;i<scaleNbr;i++){
		scaleFactors[i]=1.0f+i*scaleStep;
		}
		// import the image data into 1D arrays
		ImageDataFloat inputImg = new ImageDataFloat(inputImage.getImageData());
		nx = inputImg.getRows();
		ny = inputImg.getCols();
		nz = inputImg.getSlices();
		nc = inputImg.getComponents();
		nxyz = nx*ny*nz;
		
		float[] image= new float[nxyz];	
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			image[id] = inputImg.getFloat(x,y,z);		
		}
		
		boolean[] maskCort = new boolean[nxyz];
		float minI = 1e9f;
		float maxI = -1e9f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			// mask
			if (image[id]==0) maskCort[id] = false;
			else maskCort[id] = true;
			// remove border from computations
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) maskCort[id] = false;
			// normalize
			if (image[id]>maxI) maxI = image[id];
			if (image[id]<minI) minI = image[id];
		}
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			image[xyz] = (image[xyz]-minI)/(maxI-minI);
		}
		
		// Compute filter at different scales
		// new filter response: empty
		float[][] response = new float[scaleNbr][nxyz];
		//new direction empty
		byte[][] scaleDirection = new byte[scaleNbr][nxyz];
		for(int i=0;i<scaleNbr;i++){
			
			float scale=scaleFactors[i];
			float[] smoothed = new float[nxyz];
			if(scale==1.0f){
				smoothed=image;	
			}
			else {
				// Gaussian Kernel
				float[][] G = ImageFilters.separableGaussianKernel(scale/L2N2,scale/L2N2,scale/L2N2);
				int gx = (G[0].length-1)/2;
				int gy = (G[1].length-1)/2;
				int gz = (G[2].length-1)/2;
				
				// smoothed image
				smoothed = ImageFilters.separableConvolution(image,nx,ny,nz,G,gx,gy,gz); 
			}
			
			//Begin Filter
			boolean[] mask = new boolean[nxyz];
			float min = 1e9f;
			float max = -1e9f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				// mask
				if (smoothed[id]==0) mask[id] = false;
				else mask[id] = true;
				// remove border from computations
				if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[id] = false;
				// normalize
				if (smoothed[id]>max) max = smoothed[id];
				if (smoothed[id]<min) min = smoothed[id];
			}
			
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (brightParam.getValue().equals("bright")) {
					smoothed[xyz] = (smoothed[xyz]-min)/(max-min);
				} else if (brightParam.getValue().equals("dark")) {
					smoothed[xyz] = (max-smoothed[xyz])/(max-min);
				}
			}		
			
			byte[] direction = new byte[nxyz];
			float[] filter= new float[nxyz];
			
			if (filterParam.getValue().equals("Hessian")) {
				BasicInfo.displayMessage("estimate diffusion tensor from Hessian\n");
				directionFromHessian(smoothed, mask, filter, direction);
			} else {		
				BasicInfo.displayMessage("estimate diffusion tensor from RRF\n");
				directionFromRecursiveRidgeFilter(smoothed, mask, filter,direction);
			}
			response[i]=filter;	
			scaleDirection[i]=direction;
		}

		
		//Combine scales
		float[] finalResponse = response[0];	
		float[] scaleMax = new float[nxyz];
		for(int id=0;id<nxyz;id++){
			scaleMax[id]=scaleFactors[0];
			for(int i=1;i<scaleNbr;i++){
				finalResponse[id]=Numerics.max(finalResponse[id],response[i][id]);
				if(finalResponse[id]==response[i][id]){
					scaleMax[id]=scaleFactors[i];
				}
			}
		}
		// Adjust Direction
		byte[] finalDirection= scaleDirection[0];
		for(int id=0;id<nxyz;id++){
			for(int i=1;i<scaleNbr;i++){
				if(scaleMax[id]==scaleFactors[i]){
					finalDirection[id]=scaleDirection[i][id];
				}
			}
		}
		// Egalise histogram (Exp Median)
		float[] proba = new float[nxyz];
		if (filterParam.getValue().equals("Hessian")) {
			shapeFromHessianFilter( finalResponse, proba);	
		} else {
			shapeFromRecursiveRidgeFilter( finalResponse, proba);	
		}
		
		// generate a vector
		float[][] direction = new float[3][nxyz];
		for (int xyz=0;xyz<nxyz;xyz++){
			float[] vec = directionVector(finalDirection[xyz]);
			direction[X][xyz] = vec[X];
			direction[Y][xyz] = vec[Y];
			direction[Z][xyz] = vec[Z];
		}
		
		// 2. threshold each original element
		float threshold = Numerics.bounded(thresholdParam.getValue().floatValue(), 0.001f, 0.999f);
		boolean[] objinit = ObjectExtraction.objectFromImage(proba, nx,ny,nz, threshold, ObjectExtraction.SUPERIOR);
		
		// 3. diffuse the data to neighboring structures
		// compute neighborhood structure
		BasicInfo.displayMessage("compute neighborhood structure\n");
		float[] similarity1 = new float[nxyz];
		float[] similarity2 = new float[nxyz];
		byte[] neighbor1 = new byte[nxyz];
		byte[] neighbor2 = new byte[nxyz];
		//estimateDiffusionSimilarity(direction, mask, neighbor1, neighbor2, similarity1, similarity2);
		estimateSimpleDiffusionSimilarity(direction, maskCort, neighbor1, neighbor2, similarity1, similarity2);
		
		// original score
		float[] diffused = new float[nxyz];
		float shapeScale = 1.0f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (!maskCort[xyz]) proba[xyz] = 0.0f;

			//diffused[xyz] = proba[xyz];
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}
		
		// pre-compute different similarity models; default is angle only
		for (int xyz=0;xyz<nxyz;xyz++) if (maskCort[xyz]) {
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
			BasicInfo.displayMessage("diffusion smoothing ("+factorParam.getValue().floatValue()+")\n");
		
			proba = probaDiffusion(proba, similarity1, similarity2, neighbor1, neighbor2, maskCort, nxyz, maxiter, maxdiff, factor);
		} else if (propagationParam.getValue().equals("belief")) {
			// belief propagation
			BasicInfo.displayMessage("belief propagation\n");
		
			proba = beliefPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, maskCort, nxyz, maxiter, maxdiff);
			//proba = basicBeliefPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
			//proba = integralPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
		}
		
		// threshold the new diffused components		
		boolean[] objdiff = ObjectExtraction.objectFromImage(proba, nx,ny,nz, threshold, ObjectExtraction.SUPERIOR);
		
		// 4. estimate diameters
		BasicInfo.displayMessage("Vessel Intensity \n");
		// in and out vessels intensity initialisation
		float Iin=0.0f;
		int nIn=0;
		float[][][] inVal=new float[nx][ny][nz];	
		float minIn = 1e9f;
		float maxIn = -1e9f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			inVal[x][y][z]=0.0f;
			if(objinit[id] && scaleMax[id]>1){
				Iin+=image[id];
				inVal[x][y][z]=image[id];
				nIn++;
				if (image[id] > maxIn) maxIn = image[id];
				if (image[id] < minIn) minIn = image[id];
			}
		}
		Iin=Iin/(float)nIn;
		
		BasicInfo.displayMessage("Ring \n");
		boolean[] maskOut=new boolean[nxyz];
		boolean[] maskVote=new boolean[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			maskOut[id]=false;
			if(objinit[id]) maskVote[id]=true;
			else maskVote[id]=false;
		}		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if(objinit[xyz]){
				float[] finalDir=directionVector(finalDirection[xyz]);
				for(int d=0;d<13;d++){
					int rIn=3;
					int rOut =5;
					byte[] orthVector = directionNeighbor(d);					
					float orthTest = finalDir[X]*orthVector[X]+finalDir[Y]*orthVector[Y]+finalDir[Z]*orthVector[Z];
					if(orthTest*orthTest<1.0e-9f){
						float dist= (float)FastMath.sqrt( orthVector[X]*orthVector[X]+orthVector[Y]*orthVector[Y]+orthVector[Z]*orthVector[Z]);
						int s=1;
						if(scaleMax[xyz]>1){
							rIn*=2;
							rOut*=2;
						}
						while(dist<=rOut){	
							if(x+s*orthVector[X]>=0 && x+s*orthVector[X]<nx && (y+s*orthVector[Y])>=0 && (y+s*orthVector[Y])<ny && (z+s*orthVector[X])>=0 && (z+s*orthVector[X])<nz){
								if(dist>rIn){
									//if(!objinit[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])])	maskOut[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])]=true;
									if(!objdiff[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])])	maskOut[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])]=true;
								}
								else {
									//if(!objinit[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])])	maskVote[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])]=true;
									if(!objdiff[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])])	maskVote[x+s*orthVector[X]+nx*(y+s*orthVector[Y])+nx*ny*(z+s*orthVector[X])]=true;
								}
							}
							if(x-s*orthVector[X]>=0 && x+s*orthVector[X]<nx && (y-s*orthVector[Y])>=0 && (y-s*orthVector[Y])<ny && (z-s*orthVector[X])>=0 && (z-s*orthVector[X])<nz){
								if(dist>rIn){
									//if(!objinit[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])])	maskOut[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])]=true;
									if(!objdiff[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])])	maskOut[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])]=true;
								}
								else {
									//if(!objinit[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])])	maskVote[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])]=true;
									if(!objdiff[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])])	maskVote[x-s*orthVector[X]+nx*(y-s*orthVector[Y])+nx*ny*(z-s*orthVector[X])]=true;
								}
							}
							s++;
							dist= (float)FastMath.sqrt( orthVector[X]*orthVector[X]+orthVector[Y]*orthVector[Y]+orthVector[Z]*orthVector[Z])*(float)s;
						}
					}
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if(maskVote[xyz]){
				maskOut[xyz]=false;
			}
		}
				
		BasicInfo.displayMessage("Background Intensity \n");
		float Iout=0.0f;
		int nOut=0;
		float minOut = 1e9f;
		float maxOut = -1e9f;
		float[][][] outVal=new float[nx][ny][nz];		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			outVal[x][y][z]=0.0f;
			if(maskOut[id]){
				Iout+=image[id];
				outVal[x][y][z]=image[id];
				nOut++;
				if (image[id] > maxOut) maxOut = image[id];
				if (image[id] < minOut) minOut = image[id];
			}
		}
		Iout=Iout/(float)nOut;
		
		if(debug){
			// mean and dev calculus
			float[] vectIn=new float[nIn];
			BasicInfo.displayMessage("nIn ="+nIn+" \n");
			float[] vectOut=new float[nOut];
			BasicInfo.displayMessage("nOut ="+nOut+" \n");
	
			int no=0;       
			int ni=0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
				if(objinit[id] && scaleMax[id]>1){
					vectIn[ni]=image[id];
					ni++;
				}
				if(maskOut[id]){
					vectOut[no]=image[id];
					no++;
				}
			}
			
			//float medOut=vectOut[nOut/2];
			float medOut=Iout;
			float sigmaOut=0.0f;
			for(int n=0; n<nOut;n++){
				sigmaOut+=(vectOut[n]-medOut)*(vectOut[n]-medOut);
			}
			sigmaOut/=(float)(nOut-1);
			sigmaOut=(float)FastMath.sqrt(sigmaOut);
			BasicInfo.displayMessage("med = "+medOut+"\n");
			BasicInfo.displayMessage("sigma = "+sigmaOut+"\n");
			
			float medIn=Iin;
			float sigmaIn=0.0f;
			for(int n=0; n<nIn;n++){
				sigmaIn+=(vectIn[n]-medIn)*(vectIn[n]-medIn);
			}
			sigmaIn/=(float)(nIn-1);
			sigmaIn=(float)FastMath.sqrt(sigmaIn);
			BasicInfo.displayMessage("med = "+medIn+"\n");
			BasicInfo.displayMessage("sigma = "+sigmaIn+"\n");
		}
				
		// in probability map
		float[][][] probaInImg=new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if(maskCort[xyz]){
				//}
				if(image[xyz]>Iin){
					probaInImg[x][y][z]=1.0f;
				}
				else if (image[xyz]<Iout){
					probaInImg[x][y][z]=0.0f;
				}
				else {
					probaInImg[x][y][z]=(image[xyz]-Iout)/(Iin-Iout);
				}
			}
			else {
				probaInImg[x][y][z]=0.0f;
			}
			
		}
		
		
		// Estime diameter 
			// Estime window size
		BasicInfo.displayMessage("Diameter profile\n");
		float[][][] diamEst= new float[nx][ny][nz];
		int[][][] voisNbr= new int[nx][ny][nz];
		boolean[][] obj2=new boolean[nxyz][5]; 
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			obj2[id][0]=false;
			if(objdiff[id]){
				float[] finalDir=directionVector(finalDirection[id]);
				diamEst[x][y][z]=0.0f;
				int rIn=1;
				int vois=0;
				for (int i=-rIn;i<=rIn;i++) for (int j=-rIn;j<=rIn;j++) for (int k=-rIn;k<=rIn;k++){
					float dist= (float)FastMath.sqrt((float)i*i +(float)j*j +(float)k*k);
					if(dist<=rIn){
						if(x+i<nx &&  x+i>=0 && y+j<ny &&  y+j>=0 && z+k<nz &&  z+k>=0 ){
							float orthTest=(float)i*finalDir[X]+(float)j*finalDir[Y]+(float)k*finalDir[Z];
							orthTest=orthTest*orthTest;
							if(orthTest<1.0e-9f){
								vois++;
								diamEst[x][y][z]+=probaInImg[x+i][y+j][z+k]; 	
							}
						}
					}
				}
				diamEst[x][y][z]/=(float)vois; 
				voisNbr[x][y][z]=vois;
				if(diamEst[x][y][z]>0.5f) obj2[id][0]=true;
			}		          
		}
		
		float[][][][] diamEst2= new float[nx][ny][nz][4];
		for(int s=1;s<5;s++){
		BasicInfo.displayMessage("Initial search "+s+"\n");
		int rIn=1+s;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			obj2[id][s]=false;
			if(obj2[id][s-1]){
				float[] finalDir=directionVector(finalDirection[id]);
				diamEst2[x][y][z][s-1]=0.0f;

				int vois=0;
				for (int i=-rIn;i<=rIn;i++) for (int j=-rIn;j<=rIn;j++) for (int k=-rIn;k<=rIn;k++){
					float dist= (float)FastMath.sqrt((float)i*i +(float)j*j +(float)k*k);
					if(dist<=rIn){
						if(x+i<nx &&  x+i>=0 && y+j<ny &&  y+j>=0 && z+k<nz &&  z+k>=0 ){
							float orthTest=(float)i*finalDir[X]+(float)j*finalDir[Y]+(float)k*finalDir[Z];
							orthTest=orthTest*orthTest;
							if(orthTest<1.0e-9f){
								vois++;
								diamEst2[x][y][z][s-1]+=probaInImg[x+i][y+j][z+k]; 	
							}
						}
					}
				}
				diamEst2[x][y][z][s-1]/=(float)vois; 
				voisNbr[x][y][z]=vois;
				if(diamEst2[x][y][z][s-1]>0.5f) obj2[id][s]=true;
			}		          
		}
		}
			// Fit to model	
		double[][][] firstEst=new double[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			firstEst[x][y][z]=0.0f;
		}
		double[][][] pv=new double[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			pv[x][y][z]=0.0f;	
		}			
		BasicInfo.displayMessage("Optimization\n");
		int wtOptim=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id=x+nx*y+nx*ny*z;
			if(objdiff[id]){
				int rIn=1;
				for(int s=1;s<5;s++){
					if(obj2[id][s-1] ){			
						rIn++;	
					}			
				}
				
				float[] finalDir=directionVector(finalDirection[id]);
				float[] tan1=new float[3];
			
				boolean btan=true;
				for(int d=0;d<13;d++){
					if(btan){
						tan1=directionVector(d);
						float orthTest=(float)tan1[X]*finalDir[X]+(float)tan1[Y]*finalDir[Y]+(float)tan1[Z]*finalDir[Z];
						orthTest=orthTest*orthTest;
						if(orthTest<1.0e-9f){
							btan=false;
						}
					}
				}
				float[] tan2=new float[3];
				tan2[X]=tan1[Y]*finalDir[Z]-tan1[Z]*finalDir[Y];
				tan2[Y]=tan1[Z]*finalDir[X]-tan1[X]*finalDir[Z];
				tan2[Z]=tan1[X]*finalDir[Y]-tan1[Y]*finalDir[X];
				
				BasicInfo.displayMessage(".");
		
				VesselDiameterCostFunction jfct= new VesselDiameterCostFunction();
				jfct.setImagesData(probaInImg,nx,ny,nz);
				jfct.setPointData(x,y,z, rIn,finalDir,  tan1, tan2);
				
				
				double[] init={0.0f,0.0f,(double)rIn} ;
				SimplexOptimizer optimizer= new SimplexOptimizer(1e-5,1e-10);
				MaxEval max=new MaxEval(10000);
				ObjectiveFunction g=new ObjectiveFunction(jfct);
				
				InitialGuess start=new InitialGuess(init);
				MaxIter mIt = new MaxIter(200); 
				double[] step={0.01f,0.01f,0.01f};
				NelderMeadSimplex simplex= new NelderMeadSimplex(step);
				PointValuePair resultPair=optimizer.optimize(max,g, GoalType.MINIMIZE,start,simplex);
				double[] result= resultPair.getPoint();
				
				firstEst[x][y][z]=FastMath.abs(result[2]);	
				float maxDist= (float)rIn+0.5f;
				if(firstEst[x][y][z]>maxDist){
					firstEst[x][y][z]=0.0f;	
					wtOptim++;
				}else{				
					float xc=FastMath.round(result[0]*tan1[X]+result[1]*tan2[X]);
					float yc=FastMath.round(result[0]*tan1[Y]+result[1]*tan2[Y]);
					float zc=FastMath.round(result[0]*tan1[Z]+result[1]*tan2[Z]);
					for(int i=-rIn;i<=rIn;i++)for(int j=-rIn;j<=rIn;j++)for(int k=-rIn;k<=rIn;k++){	
						if(x+i>=0 && x+i<nx && (y+j)>=0 && (y+j)<ny && (z+k)>=0 && (z+k)<nz){
								float orthTest = finalDir[X]*(float)i+finalDir[Y]*(float)j+finalDir[Z]*(float)k;
								if(orthTest*orthTest<=1e-6){
								float dist= (float)FastMath.sqrt( ((float)i-xc)*((float)i-xc) +((float)j-yc)*((float)j-yc) +((float)k-zc)*((float)k-zc)  );	
								if(dist<=firstEst[x][y][z]-0.5f){
									pv[x+i][y+j][z+k]=1;
								}							
								else if(dist<=firstEst[x][y][z]+0.5f){
									pv[x+i][y+j][z+k]=FastMath.max(0.5f+firstEst[x][y][z]-dist,pv[x+i][y+j][z+k]);
								}
							}
						}					
					}					
				}											
			}
		}
		/*
		BasicInfo.displayMessage("Post-processing\n");
		
		boolean[] objLarge=new boolean[nxyz];
					
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x+nx*y+nx*ny*z;
			if (pv[x][y][z]>0) objLarge[id]=true;
			else objLarge[id]=false;
		}
			
		int[] vesselsPv = ObjectLabeling.connected26Object3D(objLarge, nx,ny,nz);
		int lbNbr=ObjectLabeling.countLabels(vesselsPv,nx,ny,nz);
		
		BasicInfo.displayMessage("Number of labeled region : "+lbNbr+"\n");
		
		int[] labelRegionSize=new int[lbNbr];
		float[] regionVote = new float[lbNbr];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id=x+nx*y+nx*ny*z;
			labelRegionSize[vesselsPv[id]]++;
			regionVote[vesselsPv[id]]+=probaInImg[x][y][z];
		}
		for(int i=0;i<lbNbr;i++){
			//regionVote[i]/=labelRegionSize[i];
			BasicInfo.displayMessage("Region : "+i+" score "+regionVote[i]+"\n");
		}
		
		int thresholdVote = thresholdReg.getValue().intValue();
		float thresholdProba = Numerics.bounded(thresholdParamProba.getValue().floatValue(), 0.001f, 0.999f);
		double[][][] cleanPv=new double[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x+nx*y+nx*ny*z;
			if(objLarge[id] && regionVote[vesselsPv[id]]>thresholdProba*thresholdVote)	cleanPv[x][y][z]=pv[x][y][z];
			else cleanPv[x][y][z]=0;
		}
		*/
		
		// 4. grow the region around the vessel centers to fit the scale (pv modeling)
		float[] pvol = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (objdiff[xyz]) {
				// when the diameter is smaller than the voxel, linear approx
				pvol[xyz] = Numerics.bounded((float)firstEst[x][y][z], 0.0f, 1.0f);
				int nsc = Numerics.ceil(scaleMax[xyz]/2.0f); // scale is diameter
				for (int i=-nsc;i<=nsc;i++) for (int j=-nsc;j<=nsc;j++) for (int k=-nsc;k<=nsc;k++) {
					if (x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+k>=0 && z+k<nz) {
						float dist = (float)FastMath.sqrt(i*i+j*j+k*k);
						// update the neighbor value for diameters above the voxel size
						//pvol[xyz+i+j*nx+k*nx*ny] = Numerics.max(pvol[xyz+i+j*nx+k*nx*ny], Numerics.bounded(scaleMax[xyz]/2.0f+0.5f-dist, 0.0f, 1.0f));
						pvol[xyz+i+j*nx+k*nx*ny] = Numerics.max(pvol[xyz+i+j*nx+k*nx*ny], Numerics.bounded((float)firstEst[x][y][z]/2.0f+0.5f-dist, 0.0f, 1.0f));
					}
				}
			}
		}
		
		// 5. measure the length of each segmented vessel element?
		boolean[] pvseg = ObjectExtraction.objectFromImage(pvol, nx,ny,nz, threshold, ObjectExtraction.SUPEQUAL);
		int[] vessels = ObjectLabeling.connected6Object3D(pvseg, nx,ny,nz);
		int[] vesselsizes = ObjectStatistics.volumes(vessels, ObjectLabeling.countLabels(vessels, nx,ny,nz), nx,ny,nz);
		
		// Output
		BasicInfo.displayMessage("Prepare output images\n");
		String outname = inputImg.getName();
		
		if (filterParam.getValue().equals("Hessian")) outname+="_hsn";
		else outname+="_rrf";
		
		if (propagationParam.getValue().equals("diffusion") && iterParam.getValue().intValue()>0) outname += "d";
		else if (propagationParam.getValue().equals("belief") && iterParam.getValue().intValue()>0) outname += "bp";
		
		/*
		byte[][][] pSeg = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (probaInImg[x][y][z]>thresholdProba) pSeg[x][y][z] = 1;
			else pSeg[x][y][z] = 0;
		}
		
		ImageDataUByte pSegData = new ImageDataUByte(pSeg);	
		pSegData.setHeader(inputImg.getHeader());
		pSegData.setName(outname+"_pSeg");
		probaSeg.setValue(pSegData);
		
		float thresholdPv = Numerics.bounded(thresholdParamPv.getValue().floatValue(), 0.001f, 0.999f);
		byte[][][] large = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x+nx*y+nx*ny*z;
			if (pv[x][y][z]>thresholdPv) large[x][y][z]=1;
			else	large[x][y][z]=0;
		}
				
		
		ImageDataUByte largeVesselData = new ImageDataUByte(large);	
		largeVesselData.setHeader(inputImg.getHeader());
		largeVesselData.setName(outname+"_"+scaleNbr+"_largeSeg");
		largeSeg.setValue(largeVesselData);
		
		
		ImageDataDouble  diamData = new ImageDataDouble(firstEst);	
		diamData.setHeader(inputImg.getHeader());
		diamData.setName(outname+"_diam");
		diam.setValue(diamData);
		
		ImageDataDouble  dataPv = new ImageDataDouble(pv);	
		dataPv.setHeader(inputImg.getHeader());
		dataPv.setName(outname+"_pv");
		pvEst.setValue(dataPv);
		
		ImageDataDouble  cleanDataPv = new ImageDataDouble(cleanPv);	
		cleanDataPv.setHeader(inputImg.getHeader());
		cleanDataPv.setName(outname+"_cleanPv");
		cleanPvEst.setValue(cleanDataPv);
		
		byte[][][] largeClean = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x+nx*y+nx*ny*z;
			if (cleanPv[x][y][z]>thresholdPv) largeClean[x][y][z]=1;
			else	largeClean[x][y][z]=0;
		}
		
		ImageDataUByte cleanLargeVesselData = new ImageDataUByte(largeClean);	
		cleanLargeVesselData.setHeader(inputImg.getHeader());
		cleanLargeVesselData.setName(outname+"_"+scaleNbr+"_CleanLargeSeg");
		cleanLargeSeg.setValue(cleanLargeVesselData);
		
		BasicInfo.displayMessage("diam est \n");		
				
		ImageDataFloat windowData = new ImageDataFloat(diamEst);	
		windowData.setHeader(inputImg.getHeader());
		windowData.setName(outname+"_diamEst");
		diamRange0.setValue(windowData);
		
		ImageDataFloat windowData2 = new ImageDataFloat(diamEst2);	
		windowData2.setHeader(inputImg.getHeader());
		windowData2.setName(outname+"_diamEst2");
		diamRange.setValue(windowData2);
		
		BasicInfo.displayMessage("Proba in \n");
	
			ImageDataFloat imProbaIn = new ImageDataFloat(probaInImg);
			imProbaIn.setHeader(inputImg.getHeader());
			imProbaIn.setName(outname+"_ProbaIn");
			probaNew.setValue(imProbaIn);	
		
		if(debug){
			BasicInfo.displayMessage("In values \n");
	
			ImageDataFloat imIn = new ImageDataFloat(inVal);
			imIn.setHeader(inputImg.getHeader());
			imIn.setName(outname+"_In");
			valIn.setValue(imIn);	
			
			BasicInfo.displayMessage("Out values \n");
	
			ImageDataFloat imOut = new ImageDataFloat(outVal);
			imOut.setHeader(inputImg.getHeader());
			imOut.setName(outname+"_Out");
			valOut.setValue(imOut);
			
			BasicInfo.displayMessage("Scale Map \n");
			
			float[][][] imSc= new float[nx][ny][nz]; 
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				if(obj[id]){
					imSc[x][y][z] = scaleMax[id];
				}
				else imSc[x][y][z]=0;
			}
			ImageDataFloat imSData = new ImageDataFloat(imSc);
			imSData.setHeader(inputImg.getHeader());
			imSData.setName(outname+"_scaleMap");
			scaleMap.setValue(imSData);
		}
		*/
		BasicInfo.displayMessage("seg \n");		
		byte[][][] seg = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (objdiff[id]) seg[x][y][z] = 1;
			else seg[x][y][z] = 0;
		}
		
		ImageDataUByte vesselData = new ImageDataUByte(seg);	
		vesselData.setHeader(inputImg.getHeader());
		vesselData.setName(outname+"_"+scaleNbr+"_seg");
		vesselImage.setValue(vesselData);
		
		BasicInfo.displayMessage("filter out\n");
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = finalResponse[id];
		}					
		ImageDataFloat resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_filter");
		filterImage.setValue(resData);
		
		BasicInfo.displayMessage("proba out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = proba[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_proba");
		probaImage.setValue(resData);
		
		BasicInfo.displayMessage("partial volume out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = pvol[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_pv");
		pvImage.setValue(resData);
		
		BasicInfo.displayMessage("scale out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (objdiff[id]) result[x][y][z] = scaleMax[id];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_scale");
		scaleImage.setValue(resData);
		
		BasicInfo.displayMessage("diameter out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			//if (obj[id]) result[x][y][z] = scaleMax[id];
			if (objdiff[id]) result[x][y][z] = (float)firstEst[x][y][z];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_diam");
		diameterImage.setValue(resData);
		
		BasicInfo.displayMessage("length out\n");
		result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (vessels[id]>0) result[x][y][z] = vesselsizes[vessels[id]-1];
		}					
		resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_length");
		lengthImage.setValue(resData);
		
		/*
		BasicInfo.displayMessage("label out\n");
		int[][][] vasc = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			vasc[x][y][z] = vessels[id];
		}				
		ImageDataInt labelData = new ImageDataInt(vasc);	
		labelData.setHeader(inputImg.getHeader());
		labelData.setName(outname+"_label");
		labelImage.setValue(labelData);
		*/
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
	private final void shapeFromRecursiveRidgeFilter( float[] filter, float[] shape) {
		// normalization: best is the iterative robust exponential (others are removed)
		int nb = 0;
		double max = 0.0f;
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
		double median = measure.evaluate(response, 50.0f);
		double beta = median/FastMath.log(2.0f);
		
		BasicInfo.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
			int xyz = x + nx*y + nx*ny*z;
			if (filter[xyz]>0) {
					double pe = FastMath.exp( -filter[xyz]/beta)/beta;
					shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
			}
		}
		return;
	}
	private final void shapeFromHessianFilter( float[] filter, float[] shape) {
		// normalization: best is the iterative robust exponential (others are removed)
		int nb = 0;
		double max = 0.0f;
		double sqmean = 0.0f;
		double mean = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			// keep only the proper sign
			if (filter[xyz]<=0) filter[xyz] = 0.0f;
			else {
				// fit exp only to non-zero data
				mean += filter[xyz];
				sqmean += filter[xyz]*filter[xyz];
				nb++;
				max = Numerics.max(filter[xyz], max);
			}
		}
		mean /= nb;
		sqmean /= nb;
		
		// half-normal estimates
		double sigma2 = sqmean;
		double normg = FastMath.sqrt(2.0f/(sigma2*FastMath.PI));
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (filter[xyz]>0) {
				double pg = normg*FastMath.exp(-filter[xyz]*filter[xyz]/(2.0f*sigma2));
				shape[xyz] = (float)(1.0f/max/( 1.0f/max+pg));
			}
		}
		
		// iterate
		double sigma2prev;
		for (int t=0;t<20;t++) {
			BasicInfo.displayMessage("iteration "+(t+1)+"\n");

			sigma2prev = sigma2;
				
			mean = 0.0f;
			sqmean = 0.0f;
			double den = 0.0f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (filter[xyz]>0) {
					mean += shape[xyz]*filter[xyz];
					sqmean += shape[xyz]*filter[xyz]*filter[xyz];
					den += shape[xyz];
				}
			}
			mean /= den;
			sqmean /= den;
				
			BasicInfo.displayMessage("parameter estimates: mean "+mean+", sq mean "+FastMath.sqrt(sqmean)+"\n");
				
			// half-normal estimates
			sigma2 = sqmean;
			normg = FastMath.sqrt(2.0f/(sigma2*FastMath.PI));
				
			BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
				
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)  {
				int xyz = x + nx*y + nx*ny*z;
				if (filter[xyz]>0) {
					double pg = normg*FastMath.exp(-filter[xyz]*filter[xyz]/(2.0f*sigma2));
					shape[xyz] = (float)(1.0f/max/( 1.0f/max+pg));
				}
			}
			if (2.0f*Numerics.abs(FastMath.sqrt(sigma2)-FastMath.sqrt(sigma2prev))
				/(FastMath.sqrt(sigma2)+FastMath.sqrt(sigma2prev))<0.01f) t=1000;
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
			BasicInfo.displayMessage("iteration "+(t+1)+": ");
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
			BasicInfo.displayMessage("diff "+diff+"\n");
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
			BasicInfo.displayMessage("iteration "+(t+1)+": ");
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
			
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff<maxdiff) t=maxiter;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (maxiter==0) diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0f), 0.0f, 1.0f);
			else diffused[id] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0f)/(1.0f+2.0f*factor), 0.0f, 1.0f);
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
				double maxcorr = 0.0f;
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
