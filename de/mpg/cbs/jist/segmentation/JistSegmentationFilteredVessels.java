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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
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

/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationFilteredVessels extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	private ParamVolume inputDirections;
	private ParamOption brightParam;
	private static final String[] brightTypes = {"bright","dark"};
	private ParamOption filterParam;
	private static final String[] filterTypes = {"tubular","strict_tubular","hessian","Frangi","pre-filtered"};
	
	//private ParamFloat origParam;
	//private ParamOption probaParam;
	//private static final String[] probaTypes = {"GaussianM","GaussianS","ExponentialM","ExponentialS"};
	private ParamOption similarityParam;
	private static final String[] similarityTypes = {"full","recomputed","angle","neighbor","constant"};
	private ParamOption propagationParam;
	private static final String[] propagationTypes = {"diffusion","belief"};
	private ParamFloat shapeParam;
	private ParamFloat thresholdParam;
	private ParamFloat factorParam;
	private ParamFloat diffParam;
	private ParamInteger iterParam;
	
	private ParamVolume shapeImage;
	private ParamVolume labelImage;
	private ParamVolume vesselImage;
	private ParamVolume probaImage;
	private ParamVolume similarityImage;
	private ParamVolume neighborImage;
	private ParamVolume dirImage;
	
	// global variables
	int nx, ny, nz, nc, nxyz;
	
	// numerical quantities
	private static final	float   INF=1e15f;
	private static final	float   ZERO=1e-15f;
	private static final	float	PI2 = (float)(Math.PI/2.0f);
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0f);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0f);
	private static final	float	INVSQRT2 = (float)(1.0f/FastMath.sqrt(2.0f));
	private static final	float	INVSQRT3 = (float)(1.0f/FastMath.sqrt(3.0f));
	private static final	byte	UNKNOWN = -1;
	
	// direction labeling
	private	static	final	byte	EMPTY = -1;
			
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
	
	private static final boolean debug=false;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(brightParam = new ParamOption("Vessel intensities", brightTypes));
		
		inputParams.add(filterParam = new ParamOption("Shape filter", filterTypes));
		inputParams.add(inputDirections = new ParamVolume("Input Directions (opt)"));
		inputDirections.setMandatory(false);
		
		//inputParams.add(origParam = new ParamFloat("Original image scaling", 0.0f, 10.0f, 0.5f));
		//inputParams.add(probaParam = new ParamOption("Filter probability model", probaTypes));
		inputParams.add(similarityParam = new ParamOption("Filter similarity", similarityTypes));
		inputParams.add(shapeParam = new ParamFloat("Shape filter scaling", 0.0f, 10.0f, 1.0f));
		inputParams.add(thresholdParam = new ParamFloat("Probability threshold", 0.0f, 1.0f, 0.5f));
		inputParams.add(propagationParam = new ParamOption("Propagation model", propagationTypes));
		inputParams.add(factorParam = new ParamFloat("Diffusion factor", 0.0f, 100.0f, 0.5f));
		inputParams.add(diffParam = new ParamFloat("Max difference", 0.0f, 1.0f, 0.0001f));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000, 100));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation.devel");
		inputParams.setLabel("Filtered Vessels");
		inputParams.setName("FilteredVessels");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Segments vasculature with a MRF diffusion technique.");
		
		info.setVersion("3.0.6");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(vesselImage = new ParamVolume("Segmented Vessels",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(probaImage = new ParamVolume("Vessel Probability",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(labelImage = new ParamVolume("Labeled Vessels",VoxelType.INT,-1,-1,-1,-1));
		if (debug) {
			outputParams.add(shapeImage = new ParamVolume("Vessel Shape Filter",VoxelType.FLOAT,-1,-1,-1,-1));
			outputParams.add(similarityImage = new ParamVolume("Vessel Similarity",VoxelType.FLOAT,-1,-1,-1,-1));
			outputParams.add(neighborImage = new ParamVolume("Vessel Neighbors",VoxelType.UBYTE,-1,-1,-1,-1));
			outputParams.add(dirImage = new ParamVolume("Vessel Direction",VoxelType.FLOAT,-1,-1,-1,-1));
		}		
		outputParams.setName("filtered vessel images");
		outputParams.setLabel("filtered vessel images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat inputImg = new ImageDataFloat(inputImage.getImageData());
		nx = inputImg.getRows();
		ny = inputImg.getCols();
		nz = inputImg.getSlices();
		nc = inputImg.getComponents();
		nxyz = nx*ny*nz;
		BasicInfo.displayMessage("Image dims: "+nx+", "+ny+", "+nz+"; "+nc+"\n");
		float rx = inputImg.getHeader().getDimResolutions()[0];
		float ry = inputImg.getHeader().getDimResolutions()[1];
		float rz = inputImg.getHeader().getDimResolutions()[2];
		
		int orient = inputImg.getHeader().getImageOrientation().ordinal();
		int orx = inputImg.getHeader().getAxisOrientation()[0].ordinal();
		int ory = inputImg.getHeader().getAxisOrientation()[1].ordinal();
		int orz = inputImg.getHeader().getAxisOrientation()[2].ordinal();
		
		// new proba: empty
		float[] proba = new float[nxyz];
		float min = 1e9f;
		float max = -1e9f;
		boolean[] mask = new boolean[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			proba[id] = inputImg.getFloat(x,y,z);
			// mask
			if (proba[id]==0) mask[id] = false;
			else mask[id] = true;
			// remove border from computations
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[id] = false;
			// normalize
			if (proba[id]>max) max = proba[id];
			if (proba[id]<min) min = proba[id];
		}
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (brightParam.getValue().equals("bright")) {
				proba[xyz] = (proba[xyz]-min)/(max-min);
			} else if (brightParam.getValue().equals("dark")) {
				proba[xyz] = (max-proba[xyz])/(max-min);
			}
		}
		
		// compute hessian tensor direction and shape
		float[] shape = new float[nxyz];
		float[][] direction = new float[3][nxyz];
		
		if (filterParam.getValue().equals("hessian")) {
			BasicInfo.displayMessage("estimate diffusion tensor from Hessian\n");
			shapeDirectionFromHessian(proba, mask, shape, direction);
		} else if (filterParam.getValue().equals("Frangi")) {
			BasicInfo.displayMessage("estimate diffusion tensor from Frangi's filter\n");
			shapeDirectionFrangi(proba, mask, shape, direction, shapeParam.getValue().floatValue());
		} else if (filterParam.getValue().equals("tubular")) {
			BasicInfo.displayMessage("estimate diffusion tensor from tubular filter\n");
			shapeDirectionFromTubularFilter(proba, mask, shape, direction);
		} else if (filterParam.getValue().equals("strict_tubular")) {
			BasicInfo.displayMessage("estimate diffusion tensor from tubular filter\n");
			shapeDirectionFromStrictTubularFilter(proba, mask, shape, direction);
		} else if (filterParam.getValue().equals("pre-filtered")) {
			BasicInfo.displayMessage("estimate diffusion tensor from pre-computed filter\n");
			boolean recomputeDir = false;
			if (inputDirections.getImageData()==null) recomputeDir = true;
			else {
				ImageDataFloat inputDir = new ImageDataFloat(inputDirections.getImageData());
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int id = x + nx*y + nx*ny*z;
					direction[X][id] = inputDir.getFloat(x,y,z, X);
					direction[Y][id] = inputDir.getFloat(x,y,z, Y);
					direction[Z][id] = inputDir.getFloat(x,y,z, Z);
					// renormalize
					float norm = (float)FastMath.sqrt(direction[X][id]*direction[X][id]+direction[Y][id]*direction[Y][id]+direction[Z][id]*direction[Z][id]);
					direction[X][id] /= norm;
					direction[Y][id] /= norm;
					direction[Z][id] /= norm;
				}
			}		
			shapeDirectionFromPrecomputedFilter(proba, mask, shape, direction, recomputeDir);
		}
		
		// compute neighborhood structure
		BasicInfo.displayMessage("compute neighborhood structure\n");
		float[] similarity1 = new float[nxyz];
		float[] similarity2 = new float[nxyz];
		byte[] neighbor1 = new byte[nxyz];
		byte[] neighbor2 = new byte[nxyz];
		//estimateDiffusionSimilarity(direction, mask, neighbor1, neighbor2, similarity1, similarity2);
		estimateSimpleDiffusionSimilarity(direction, mask, neighbor1, neighbor2, similarity1, similarity2);
		
		// 1. original score
		float[] diffused = new float[nxyz];
		float shapeScale = 1.0f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (mask[xyz]) {
				proba[xyz] = shape[xyz];
			} else {
				proba[xyz] = 0.0f;
			}
			//diffused[xyz] = proba[xyz];
			diffused[xyz] = (float)FastMath.log(1.0f + proba[xyz]);
		}
		
		// pre-compute different similarity models; default is angle only
		if (similarityParam.getValue().equals("full")) {
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				// value more neighbors that are well aligned and have strong tubular response
				similarity1[xyz] = (float)FastMath.sqrt(similarity1[xyz]*proba[neighborIndex(neighbor1[xyz],xyz)]);
				similarity2[xyz] = (float)FastMath.sqrt(similarity2[xyz]*proba[neighborIndex(neighbor2[xyz],xyz)]);
			}
		} else if (similarityParam.getValue().equals("neighbor")) {
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				// value more neighbors that are well aligned and have strong tubular response
				similarity1[xyz] = proba[neighborIndex(neighbor1[xyz],xyz)];
				similarity2[xyz] = proba[neighborIndex(neighbor2[xyz],xyz)];
			}
		} else if (similarityParam.getValue().equals("constant")) {
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				// value more neighbors that are well aligned and have strong tubular response
				similarity1[xyz] = 1.0f;
				similarity2[xyz] = 1.0f;
			}
		}
		
		// 2. diffuse along neighborhood
		int maxiter = iterParam.getValue().intValue();
		float maxdiff = diffParam.getValue().floatValue();
		float factor = factorParam.getValue().floatValue();
		
		if (propagationParam.getValue().equals("diffusion")) {
			// diffuse along the tensor
			BasicInfo.displayMessage("diffusion smoothing ("+factorParam.getValue().floatValue()+")\n");
		
			proba = probaDiffusion(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff, factor);
		} else if (propagationParam.getValue().equals("belief")) {
			// belief propagation
			BasicInfo.displayMessage("belief propagation\n");
		
			proba = beliefPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
			//proba = basicBeliefPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
			//proba = integralPropagation(proba, similarity1, similarity2, neighbor1, neighbor2, mask, nxyz, maxiter, maxdiff);
		}
		
		// 3. threshold, then label each connected component		
		float threshold = Numerics.bounded(thresholdParam.getValue().floatValue(), 0.001f, 0.999f);
		boolean[] obj = ObjectExtraction.objectFromImage(proba, nx,ny,nz, threshold, ObjectExtraction.SUPERIOR);
		int[] vessels = ObjectLabeling.connected26Object3D(obj, nx,ny,nz);
		
		// outputs
		byte[][][] seg = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (obj[id]) seg[x][y][z] = 1;
			else seg[x][y][z] = 0;
		}
		obj = null;

		String outname = inputImg.getName()+"_"+filterParam.getValue().substring(0,3);
		if (maxiter>0) outname+="d";
		
		ImageDataUByte vesselData = new ImageDataUByte(seg);	
		vesselData.setHeader(inputImg.getHeader());
		vesselData.setName(outname+"_seg");
		vesselImage.setValue(vesselData);
		
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			result[x][y][z] = proba[id];
		}
		proba = null;
		
		ImageDataFloat resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_proba");
		probaImage.setValue(resData);
		
		int[][][] vasc = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			vasc[x][y][z] = vessels[id];
		}
		vessels = null;
		
		ImageDataInt labelData = new ImageDataInt(vasc);	
		labelData.setHeader(inputImg.getHeader());
		labelData.setName(outname+"_label");
		labelImage.setValue(labelData);
		
		/*
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
			*//*
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
				if (similarityParam.getValue().equals("recomputed")) {
					w1 = (float)FastMath.sqrt(w1*ngb1);
					w2 = (float)FastMath.sqrt(w2*ngb2);
					// need to avoid saturation?
					w1 = Numerics.min(w1,1.0f);
					w2 = Numerics.min(w2,1.0f);
				}
				
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
		proba = null;
		
		// 3. threshold, then label each connected component		
		//float threshold = (1.0f+2.0f*factor)*Numerics.bounded(thresholdParam.getValue().floatValue(), 0.001f, 0.999f);
		float threshold = 0;
		if (maxiter==0) threshold = (float)FastMath.log(1.0f+Numerics.bounded(thresholdParam.getValue().floatValue(), 0.0f00001f, 0.999999f));
		else threshold = (float)FastMath.log(1.0f+(1.0f+2.0f*factor)*Numerics.bounded(thresholdParam.getValue().floatValue(), 0.0f00001f, 0.999999f));
		boolean[] obj = ObjectExtraction.objectFromImage(diffused, nx,ny,nz, threshold, ObjectExtraction.SUPERIOR);
		int[] vessels = ObjectLabeling.connected26Object3D(obj, nx,ny,nz);
		
		// outputs
		byte[][][] seg = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (obj[id]) seg[x][y][z] = 1;
			else seg[x][y][z] = 0;
		}
		obj = null;

		String outname = inputImg.getName()+"_"+filterParam.getValue().substring(0,3);
		if (maxiter>0) outname+="d";
		
		ImageDataUByte vesselData = new ImageDataUByte(seg);	
		vesselData.setHeader(inputImg.getHeader());
		vesselData.setName(outname+"_seg");
		vesselImage.setValue(vesselData);
		
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			if (maxiter==0) result[x][y][z] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0f), 0.0f, 1.0f);
			else result[x][y][z] = Numerics.bounded((float)(FastMath.exp(diffused[id])-1.0f)/(1.0f+2.0f*factor), 0.0f, 1.0f);
			//result[x][y][z] = diffused[id]/(1.0f+2.0f*factor);
		}
		diffused = null;
		
		ImageDataFloat resData = new ImageDataFloat(result);	
		resData.setHeader(inputImg.getHeader());
		resData.setName(outname+"_proba");
		probaImage.setValue(resData);
		
		int[][][] vasc = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int id = x + nx*y + nx*ny*z;
			vasc[x][y][z] = vessels[id];
		}
		vessels = null;
		
		ImageDataInt labelData = new ImageDataInt(vasc);	
		labelData.setHeader(inputImg.getHeader());
		labelData.setName(outname+"_label");
		labelImage.setValue(labelData);
		*/
		
		if (debug) {
			// debug
			float[][][] sh = new float[nx][ny][nz]; 
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				sh[x][y][z] = shape[id];
			}
			resData = new ImageDataFloat(sh);	
			resData.setHeader(inputImg.getHeader());
			resData.setName(inputImg.getName()+"_shape");
			shapeImage.setValue(resData);
				
			float[][][][] sim = new float[nx][ny][nz][2]; 
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				sim[x][y][z][0] = similarity1[id];
				sim[x][y][z][1] = similarity2[id];
			}
			resData = new ImageDataFloat(sim);	
			resData.setHeader(inputImg.getHeader());
			resData.setName(inputImg.getName()+"_sim");
			similarityImage.setValue(resData);
			
			byte[][][][] ngb = new byte[nx][ny][nz][2]; 
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				ngb[x][y][z][0] = neighbor1[id];
				ngb[x][y][z][1] = neighbor2[id];
			}
			ImageDataUByte bytData = new ImageDataUByte(ngb);	
			bytData.setHeader(inputImg.getHeader());
			bytData.setName(inputImg.getName()+"_ngb");
			neighborImage.setValue(bytData);
			
			sim = new float[nx][ny][nz][3]; 
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int id = x + nx*y + nx*ny*z;
				sim[x][y][z][X] = direction[X][id];
				sim[x][y][z][Y] = direction[Y][id];
				sim[x][y][z][Z] = direction[Z][id];
			}
			resData = new ImageDataFloat(sim);	
			resData.setHeader(inputImg.getHeader());
			resData.setName(inputImg.getName()+"_dir");
			dirImage.setValue(resData);
		}		
		return;
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
				if (similarityParam.getValue().equals("recomputed")) {
					w1 = (float)FastMath.sqrt(w1*ngb1);
					w2 = (float)FastMath.sqrt(w2*ngb2);
					// need to avoid saturation?
					w1 = Numerics.min(w1,1.0f);
					w2 = Numerics.min(w2,1.0f);
				}
				
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
	
	private final float[] integralPropagation(float[] proba, float[] similarity1, float[] similarity2, byte[] ngb1, byte[] ngb2, boolean[] mask, int nxyz, int iter, float maxdiff) {
		
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
			message1fg[xyz] = 0.0f;
			message2fg[xyz] = 0.0f;
			message1bg[xyz] = 0.0f;
			message2bg[xyz] = 0.0f;
			newmsg1fg[xyz] = 0.0f;
			newmsg2fg[xyz] = 0.0f;
			newmsg1bg[xyz] = 0.0f;
			newmsg2bg[xyz] = 0.0f;
		}
		
		// min value for probas, similarities
		float minp = 0.00f;
		float mins = 0.00f;
		
		float offset = 0.0f;
		
		// message passing
		for (int t=0;t<iter;t++) {
			BasicInfo.displayMessage("iteration "+(t+1)+": ");
			float diff = 0.0f;
			//float prev;
			int ngb;
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float simfgfg1 = similarity1[xyz];
				float simbgbg1 = similarity1[xyz];
				
				float simfgfg2 = similarity2[xyz];
				float simbgbg2 = similarity2[xyz];
				
				// first neighbor message, fg label
				ngb = neighborIndex(ngb1[xyz],xyz);
				float mfgfg1 = Numerics.max(mins,simfgfg1)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) mfgfg1 += message1fg[ngb]-offset;
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) mfgfg1 += message2fg[ngb]-offset;
				newmsg1fg[xyz] = mfgfg1;
				diff = Numerics.max(diff, Numerics.abs(newmsg1fg[xyz]-message1fg[xyz]));
				
				// first neighbor message, bg label
				ngb = neighborIndex(ngb1[xyz],xyz);
				float mbgbg1 = Numerics.max(mins,simbgbg1)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) mbgbg1 += message1bg[ngb]-offset;
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) mbgbg1 += message2bg[ngb]-offset;
				newmsg1bg[xyz] = mbgbg1;
				diff = Numerics.max(diff, Numerics.abs(newmsg1bg[xyz]-message1bg[xyz]));
				
				// second neighbor message, fg label
				ngb = neighborIndex(ngb2[xyz],xyz);
				float mfgfg2 = Numerics.max(mins,simfgfg2)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) mfgfg2 += message1fg[ngb]-offset;
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) mfgfg2 += message2fg[ngb]-offset;
				newmsg2fg[xyz] = mfgfg2;
				diff = Numerics.max(diff, Numerics.abs(newmsg2fg[xyz]-message2fg[xyz]));
				
				// second neighbor message, bg label
				ngb = neighborIndex(ngb2[xyz],xyz);
				float mbgbg2 = Numerics.max(mins,simbgbg2)*Numerics.max(minp,(1.0f-proba[ngb]));
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) mbgbg2 += message1bg[ngb]-offset;
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) mbgbg2 += message2bg[ngb]-offset;
				newmsg2bg[xyz] = mbgbg2;
				diff = Numerics.max(diff, Numerics.abs(newmsg2bg[xyz]-message2bg[xyz]));
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
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
			belief[xyz] = proba[xyz] + message1fg[xyz] + message2fg[xyz] - 2.0f*offset;
			bgbelief = 1.0f-proba[xyz] + message1bg[xyz] + message2bg[xyz] - 2.0f*offset;
			
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
		}
		return belief;
	}
	
	private final float[] basicBeliefPropagation(float[] proba, float[] similarity1, float[] similarity2, byte[] ngb1, byte[] ngb2, boolean[] mask, int nxyz, int iter, float maxdiff) {
		
		float[] message1 = new float[nxyz];
		float[] message2 = new float[nxyz];
		
		float[] newmsg1 = new float[nxyz];
		float[] newmsg2 = new float[nxyz];
		
		// init
		for (int xyz=0;xyz<nxyz;xyz++) {
			message1[xyz] = 1.0f;
			message2[xyz] = 1.0f;
			newmsg1[xyz] = 1.0f;
			newmsg2[xyz] = 1.0f;
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
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				float sim1 = similarity1[xyz];
				float sim2 = similarity2[xyz];
				
				// first neighbor message, fg label
				ngb = neighborIndex(ngb1[xyz],xyz);
				float msg1 = Numerics.max(mins,sim1)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) msg1 *= message1[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) msg1 *= message2[ngb];
				newmsg1[xyz] = msg1;
				diff = Numerics.max(diff, Numerics.abs(newmsg1[xyz]-message1[xyz]));
				
				// second neighbor message, fg label
				ngb = neighborIndex(ngb2[xyz],xyz);
				float msg2 = Numerics.max(mins,sim2)*Numerics.max(minp,proba[ngb]);
				if (neighborIndex(ngb1[ngb],ngb)!=xyz) msg2 *= message1[ngb];
				if (neighborIndex(ngb2[ngb],ngb)!=xyz) msg2 *= message2[ngb];
				newmsg2[xyz] = msg2;
				diff = Numerics.max(diff, Numerics.abs(newmsg2[xyz]-message2[xyz]));
			}
			// copy new messages onto older ones
			for (int xyz=0;xyz<nxyz;xyz++) if (mask[xyz]) {
				message1[xyz] = newmsg1[xyz];
				message2[xyz] = newmsg2[xyz];
			}
			BasicInfo.displayMessage("diff "+diff+"\n");
			if (diff < maxdiff) t = iter;
		}
		
		// compute final belief
		float[] belief = new float[nxyz];
		float 	bgbelief;
		
		for (int xyz=0;xyz<nxyz;xyz++) {
			belief[xyz] = proba[xyz]*message1[xyz]*message2[xyz];
			//bgbelief = ??
			
			// normalize
			//belief[xyz] = belief[xyz]/Numerics.max(1e-9f,belief[xyz]+bgbelief);
		}
		return belief;
	}
	
	private final void shapeDirectionFrangi(float[] img, boolean[] mask, float[] shape, float[][] direction, float scaling) {
        int x,y,z,n;
		double[][] hessian = new double[3][3];
		
		double imin=1e16, imax=-1e16;
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
        	int id = x + nx*y + nx*ny*z;
			if (mask[id]) {
				if (img[id]<imin) imin = img[id];
				if (img[id]>imax) imax = img[id];
		    }
		}
		double a2 = 2.0f*scaling*scaling;
		double b2 = 2.0f*scaling*scaling;
		double c2 = 2.0f*scaling*(imax-imin)*scaling*(imax-imin);
					
		
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
        	int id = x + nx*y + nx*ny*z;
			if (mask[id]) {
				ImageGeometry.computeHessianAt(img, id, hessian, nx, ny, nz);
				
				// compute eigenvalues
				double[] vals = Matrix3D.eigenvalues(hessian);
				double[][] dirs = Matrix3D.eigenvectors(hessian, vals);
				
				// shape: Frangi score
				if (vals[0]<0 && vals[1]<0) {
					double Rb = Numerics.abs(vals[2])/FastMath.sqrt(Numerics.abs(vals[0]*vals[1]));
					double Ra = Numerics.abs(vals[1])/Numerics.abs(vals[0]);
					double S = FastMath.sqrt(vals[0]*vals[0]+vals[1]*vals[1]+vals[2]*vals[2]);
					
					shape[id] = (float)( (1.0f-FastMath.exp(-Ra*Ra/a2))*FastMath.exp(-Rb*Rb/b2)*(1-FastMath.exp(-S*S/c2)) );
				}
				
				// main direction: lowest eigenvalue
				direction[X][id] = (float)dirs[X][2];
				direction[Y][id] = (float)dirs[Y][2];
				direction[Z][id] = (float)dirs[Z][2];
				
			}
		}

		return;
    }
        
	private final void shapeDirectionFromHessian(float[] img, boolean[] mask, float[] shape, float[][] direction) {
       double[][] hessian = new double[3][3];
		       
        for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int id = x + nx*y + nx*ny*z;
			if (mask[id]) {
				ImageGeometry.computeHessianAt(img, id, hessian, nx, ny, nz);
				
				// compute eigenvalues
				double[] vals = Matrix3D.eigenvalues(hessian);
				double[][] dirs = Matrix3D.eigenvectors(hessian, vals);
				
				// shape: basic score is not good enough (too much background noise)
				if (vals[0]!=0) shape[id] = (float)(-vals[0]*(Numerics.abs(vals[0])-Numerics.abs(vals[2]))/Numerics.abs(vals[0]));
				
				// main direction: lowest eigenvalue
				direction[X][id] = (float)dirs[X][2];
				direction[Y][id] = (float)dirs[Y][2];
				direction[Z][id] = (float)dirs[Z][2];
				
			}
		}

		// normalization
		double mean = 0.0f;
		double max = 0.0f;
		double den = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (shape[xyz]<0) shape[xyz] = 0.0f;
				else {
					// fit only to non-zero data: much nicer!
					mean += shape[xyz];
					den++;
					max = Numerics.max(Numerics.abs(shape[xyz]), max);
				}
			}
		}
		mean /= den;
		
		// exponential fit or gaussian fit ? gaussian seems more sensitive, better for simple tubular filter
		double normg = 2.0f/(mean*FastMath.PI);
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				//double pe = FastMath.exp( -0.5f*shape[xyz]/mean)/mean;
				//shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
				double pg = normg*FastMath.exp(-shape[xyz]*shape[xyz]/(FastMath.PI*mean*mean));
				shape[xyz] = (float)(1.0f/max/( 1.0f/max+pg));
			}
		}

		return;
    }
    
	private final void shapeDirectionFromTubularFilter(float[] img, boolean[] mask, float[] shape, float[][] direction) {
		
		double avg = 0.0f;
		double sig = 0.0f;
		double den = 0.0f;
		float[] filter = new float[nx*ny*nz];
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			filter[xyz] = 0.0f;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// check for zero-valued neighbors as well
				float best = 0.0f;
				int dmax = 13;
				int dbest = -1;
				for (int d=0;d<dmax;d++) {
					float val = tubularScore(img, x,y,z,d);
					if (val*val>best*best) {
						best = val;
						dbest = d;
					}
				}
				filter[xyz] = best;
				float[] vec = directionVector(dbest);
				direction[X][xyz] = vec[X];
				direction[Y][xyz] = vec[Y];
				direction[Z][xyz] = vec[Z];
				avg += best;
				den++;
			}
		}
		// normalization: best is the iterative robust exponential (others are removed)
		int nb = 0;
		double max = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (filter[xyz]<=0) filter[xyz] = 0.0f;
				else {
					// fit exp only to non-zero data
					nb++;
					max = Numerics.max(Numerics.abs(filter[xyz]), max);
				}
			}
		}
		
		// robust measures
		double[] response = new double[nb];
		int n=0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (filter[xyz]>0) {
					response[n] = filter[xyz];
					n++;
				}
			}
		}
		Percentile measure = new Percentile();
		double median = measure.evaluate(response, 50.0f);
		double beta = median/FastMath.log(2.0f);
	
		BasicInfo.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
		
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				double pe = FastMath.exp( -filter[xyz]/beta)/beta;
				shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
			}
		}
		/*
		// iterative?
		double betaprev;
		for (int t=0;t<20;t++) {
			BasicInfo.displayMessage("iteration "+(t+1)+"\n");
			betaprev = beta;
				
			nb = 0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					if (filter[xyz]>0 && shape[xyz]<0.5f) nb++;
				}
			}
			response = new double[nb];
			n=0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					if (filter[xyz]>0 && shape[xyz]<0.5f) {
						response[n] = filter[xyz];
						n++;
					}
				}
			}
			median = measure.evaluate(response, 50.0f);
			beta = median/FastMath.log(2.0f);
			BasicInfo.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
			
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					double pe = FastMath.exp( -filter[xyz]/beta)/beta;
					shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
				}
			}
			if (2.0f*Numerics.abs(beta-betaprev)/(beta+betaprev)<0.01f) t=1000;
		}
		*/
		return;
    }
    
	private final void shapeDirectionFromStrictTubularFilter(float[] img, boolean[] mask, float[] shape, float[][] direction) {
		
		// get the tubular filter response
		float[][][] planescore = new float[nx][ny][nz];
		byte[][][] planedir = new byte[nx][ny][nz];
		byte[][][] linedir = new byte[nx][ny][nz];
		float[][][] image = new float[nx][ny][nz];
		float[] filter = new float[nx*ny*nz];
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
					//filter[xyz] = Numerics.sign(linescore)*Numerics.min(planescore[x][y][z],linescore);
				} else {
					filter[xyz] = 0.0f;
				}
				float[] vec = directionVector(linedir[x][y][z]);
				direction[X][xyz] = vec[X];
				direction[Y][xyz] = vec[Y];
				direction[Z][xyz] = vec[Z];
			}
		}
		planescore = null;
		planedir = null;
		linedir = null;
		image = null;
		
		// normalization: best is the iterative robust exponential (others are removed)
		int nb = 0;
		double max = 0.0f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			//if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (filter[xyz]<=0) filter[xyz] = 0.0f;
				else {
					// fit exp only to non-zero data
					nb++;
					//max = Numerics.max(Numerics.abs(filter[xyz]), max);
					max = Numerics.max(filter[xyz], max);
				}
			//}
		}
		
		// robust measures
		double[] response = new double[nb];
		int n=0;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			//if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				// keep only the proper sign
				if (filter[xyz]>0) {
					response[n] = filter[xyz];
					n++;
				}
			//}
		}
		Percentile measure = new Percentile();
		double median = measure.evaluate(response, 50.0f);
		double beta = median/FastMath.log(2.0f);
	
		BasicInfo.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
		
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (filter[xyz]>0) {
			//if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				double pe = FastMath.exp( -filter[xyz]/beta)/beta;
				shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
			}
		}
		/*
		// iterative?
		double betaprev;
		for (int t=0;t<20;t++) {
			BasicInfo.displayMessage("iteration "+(t+1)+"\n");
			betaprev = beta;
				
			nb = 0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					if (filter[xyz]>0 && shape[xyz]<0.5f) nb++;
				}
			}
			response = new double[nb];
			n=0;
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					if (filter[xyz]>0 && shape[xyz]<0.5f) {
						response[n] = filter[xyz];
						n++;
					}
				}
			}
			median = measure.evaluate(response, 50.0f);
			beta = median/FastMath.log(2.0f);
			BasicInfo.displayMessage("parameter estimates: median "+median+", beta "+beta+",\n");
			
			for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
					double pe = FastMath.exp( -filter[xyz]/beta)/beta;
					shape[xyz] = (float)(1.0f/max/( 1.0f/max+pe));
				}
			}
			if (2.0f*Numerics.abs(beta-betaprev)/(beta+betaprev)<0.01f) t=1000;
		}
		*/
		return;
    }
    
	private final void shapeDirectionFromPrecomputedFilter(float[] img, boolean[] mask, float[] shape, float[][] direction, boolean recompute) {
		
		float max = -1e9f;
		float min = 1e9f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			shape[xyz] = 0.0f;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				if (recompute) {
					// check for zero-valued neighbors as well
					float best = 0.0f;
					int dmax = 13;
					int dbest = -1;
					for (int d=0;d<dmax;d++) {
						byte[] ngb = directionNeighbor(d);
						float val = img[xyz+ngb[X]+nx*ngb[Y]+nx*ny*ngb[Z]]+img[xyz-ngb[X]-nx*ngb[Y]-nx*ny*ngb[Z]];
						if (val>best) {
							best = val;
							dbest = d;
						}
					}
					float[] vec = directionVector(dbest);
					direction[X][xyz] = vec[X];
					direction[Y][xyz] = vec[Y];
					direction[Z][xyz] = vec[Z];
				}
				shape[xyz] = img[xyz];
				min = Numerics.min(min, img[xyz]);
				max = Numerics.max(max, img[xyz]);
			}
		}
		// proba: just normalize the filter results
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz] && !zeroNeighbor(img, mask, x,y,z,2)) {
				shape[xyz] = (shape[xyz]-min)/(max-min);
			}
		}
		return;
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
    
    // too complicated?
    private final void estimateDiffusionSimilarity(float[][] imdir, boolean[] mask, byte[] neighbor1, byte[] neighbor2, 
    													float[] similarity1, float[] similarity2) {
		float[]	ds1 = new float[NC2];
		float[]	ds2 = new float[NC2];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int id = x+nx*y+nx*ny*z;
			if (mask[id]) {
				for (byte d=0;d<NC2;d++) {
					// lower-bounded by zero (not assuming conflicts between vessels)
					ds1[d] = 0.0f;
					ds2[d] = 0.0f;
					
					int idn = neighborIndex(d,id);
					
					if (mask[idn]) {
						if (sameDirection(d, id, imdir)) {
							ds1[d] = signedAngleSimilarityGain(imdir, d,id,idn);
						} else {
							ds2[d] = signedAngleSimilarityGain(imdir, d,id,idn);
						}
					}
				}
				// find the best for each side
				neighbor1[id] = (byte)Numerics.argmax(ds1);
				neighbor2[id] = (byte)Numerics.argmax(ds2);
				similarity1[id] = Numerics.max(ds1[neighbor1[id]],0.0f);
				similarity2[id] = Numerics.max(ds2[neighbor2[id]],0.0f);
			}
		}
		
		return;
	} // 
	
	private final float signedAngleSimilarityGain(float[][] imdir, int dir, int id1, int id2) {
		float prod = (imdir[0][id1]*imdir[0][id2]+imdir[1][id1]*imdir[1][id2]+imdir[2][id1]*imdir[2][id2]);
		float dataAngle = (float)FastMath.acos(Numerics.min(1.0f,Numerics.abs(prod)))/PI2;
		
		float dv1=0, dv2=0;
			
			 if (dir==X) { dv1 = imdir[0][id1]; dv2 = imdir[0][id2]; }
		else if (dir==Y) { dv1 = imdir[1][id1]; dv2 = imdir[1][id2]; }
		else if (dir==Z) { dv1 = imdir[2][id1]; dv2 = imdir[2][id2]; }
		else if (dir==XpY) { dv1 = imdir[0][id1]/SQRT2+imdir[1][id1]/SQRT2; dv2 = imdir[0][id2]/SQRT2+imdir[1][id2]/SQRT2; }
		else if (dir==YpZ) { dv1 = imdir[1][id1]/SQRT2+imdir[2][id1]/SQRT2; dv2 = imdir[1][id2]/SQRT2+imdir[2][id2]/SQRT2; }
		else if (dir==ZpX) { dv1 = imdir[2][id1]/SQRT2+imdir[0][id1]/SQRT2; dv2 = imdir[2][id2]/SQRT2+imdir[0][id2]/SQRT2; }
		else if (dir==XmY) { dv1 = imdir[0][id1]/SQRT2-imdir[1][id1]/SQRT2; dv2 = imdir[0][id2]/SQRT2-imdir[1][id2]/SQRT2; }
		else if (dir==YmZ) { dv1 = imdir[1][id1]/SQRT2-imdir[2][id1]/SQRT2; dv2 = imdir[1][id2]/SQRT2-imdir[2][id2]/SQRT2; }
		else if (dir==ZmX) { dv1 = imdir[2][id1]/SQRT2-imdir[0][id1]/SQRT2; dv2 = imdir[2][id2]/SQRT2-imdir[0][id2]/SQRT2; }
		else if (dir==XpYpZ) { dv1 = imdir[0][id1]/SQRT3+imdir[1][id1]/SQRT3+imdir[2][id1]/SQRT3; dv2 = imdir[0][id2]/SQRT3+imdir[1][id2]/SQRT3+imdir[2][id2]/SQRT3; }
		else if (dir==XmYpZ) { dv1 = imdir[0][id1]/SQRT3-imdir[1][id1]/SQRT3+imdir[2][id1]/SQRT3; dv2 = imdir[0][id2]/SQRT3-imdir[1][id2]/SQRT3+imdir[2][id2]/SQRT3; }
		else if (dir==XpYmZ) { dv1 = imdir[0][id1]/SQRT3+imdir[1][id1]/SQRT3-imdir[2][id1]/SQRT3; dv2 = imdir[0][id2]/SQRT3+imdir[1][id2]/SQRT3-imdir[2][id2]/SQRT3; }
		else if (dir==XmYmZ) { dv1 = imdir[0][id1]/SQRT3-imdir[1][id1]/SQRT3-imdir[2][id1]/SQRT3; dv2 = imdir[0][id2]/SQRT3-imdir[1][id2]/SQRT3-imdir[2][id2]/SQRT3; }
		else if (dir==mX) { dv1 = -imdir[0][id1]; dv2 = -imdir[0][id2]; }
		else if (dir==mY) { dv1 = -imdir[1][id1]; dv2 = -imdir[1][id2]; }
		else if (dir==mZ) { dv1 = -imdir[2][id1]; dv2 = -imdir[2][id2]; }
		else if (dir==mXpY) { dv1 = -(imdir[0][id1]/SQRT2+imdir[1][id1]/SQRT2); dv2 = -(imdir[0][id2]/SQRT2+imdir[1][id2]/SQRT2); }
		else if (dir==mYpZ) { dv1 = -(imdir[1][id1]/SQRT2+imdir[2][id1]/SQRT2); dv2 = -(imdir[1][id2]/SQRT2+imdir[2][id2]/SQRT2); }
		else if (dir==mZpX) { dv1 = -(imdir[2][id1]/SQRT2+imdir[0][id1]/SQRT2); dv2 = -(imdir[2][id2]/SQRT2+imdir[0][id2]/SQRT2); }
		else if (dir==mXmY) { dv1 = -(imdir[0][id1]/SQRT2-imdir[1][id1]/SQRT2); dv2 = -(imdir[0][id2]/SQRT2-imdir[1][id2]/SQRT2); }
		else if (dir==mYmZ) { dv1 = -(imdir[1][id1]/SQRT2-imdir[2][id1]/SQRT2); dv2 = -(imdir[1][id2]/SQRT2-imdir[2][id2]/SQRT2); }
		else if (dir==mZmX) { dv1 = -(imdir[2][id1]/SQRT2-imdir[0][id1]/SQRT2); dv2 = -(imdir[2][id2]/SQRT2-imdir[0][id2]/SQRT2); }
		else if (dir==mXpYpZ) { dv1 = -(imdir[0][id1]/SQRT3+imdir[1][id1]/SQRT3+imdir[2][id1]/SQRT3); dv2 = -(imdir[0][id2]/SQRT3+imdir[1][id2]/SQRT3+imdir[2][id2]/SQRT3); }
		else if (dir==mXmYpZ) { dv1 = -(imdir[0][id1]/SQRT3-imdir[1][id1]/SQRT3+imdir[2][id1]/SQRT3); dv2 = -(imdir[0][id2]/SQRT3-imdir[1][id2]/SQRT3+imdir[2][id2]/SQRT3); }
		else if (dir==mXpYmZ) { dv1 = -(imdir[0][id1]/SQRT3+imdir[1][id1]/SQRT3-imdir[2][id1]/SQRT3); dv2 = -(imdir[0][id2]/SQRT3+imdir[1][id2]/SQRT3-imdir[2][id2]/SQRT3); }
		else if (dir==mXmYmZ) { dv1 = -(imdir[0][id1]/SQRT3-imdir[1][id1]/SQRT3-imdir[2][id1]/SQRT3); dv2 = -(imdir[0][id2]/SQRT3-imdir[1][id2]/SQRT3-imdir[2][id2]/SQRT3); }

		float dirAngle = (float)FastMath.acos(Numerics.min(1.0f,Numerics.max(Numerics.abs(dv1),Numerics.abs(dv2))))/PI2;
		
		return (1.0f-dirAngle)*(1.0f-2.0f*dataAngle);
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

	private final boolean sameDirection(int dir, int id, float[][] imdir) {
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

		if (dv>0) return true;
		else return false;
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
	
	/*
	private final float[][][] probaDiffusion(float[][][] orig, boolean[][][] mask, int iter, float scale, float factor) {
    	
		float[][][] map = new float[nx][ny][nz];
    	for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    		map[x][y][z] = orig[x][y][z];
    	}
				
    	// propagate the values : diffusion
    	float maxdiff = 1.0f;
		for (int t=0;t<iter && maxdiff>0.0f5f;t++) {
			BasicInfo.displayMessage(".");
			maxdiff = 0.0f;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				if (mask[x][y][z]) {
					//float den = 1.0f;
					float den = diffusionWeightFunction(map[x][y][z],orig[x][y][z],scale);
					float num = den*orig[x][y][z];
					float weight;
					float prev = map[x][y][z];
						
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x-1][y][z],scale);
					num += weight*map[x-1][y][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x+1][y][z],scale);
					num += weight*map[x+1][y][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y-1][z],scale);
					num += weight*map[x][y-1][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y+1][z],scale);
					num += weight*map[x][y+1][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y][z-1],scale);
					num += weight*map[x][y][z-1];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y][z+1],scale);
					num += weight*map[x][y][z+1];
					den += weight;
	
					map[x][y][z] = num/den;
					
					maxdiff = Numerics.max(maxdiff, Numerics.abs(map[x][y][z]-prev));
				} else {
					map[x][y][z] = orig[x][y][z];
				}
			}
			BasicInfo.displayMessage("max diff. "+maxdiff+"\n");
		}
		return map;
    }*/
    
    private final float semiProduct(float s, float g) {
    	if (s>0 || g>0) return s*g;
    	else return 0.0f;
    }
    
    private final float diffusionWeightFunction(float val, float ngb, float scale) {
    	
    	//return 1.0f;
    	//return 1.0f/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	//return 1.0f/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	// making it exact for the cases p=0 or 1 ?
    	return Numerics.square(Numerics.min(1.0f-val, val)/scale)/(1.0f+Numerics.square(Numerics.min(1.0f-val, val)/scale))
    				/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    }
    
    private final float probaWeight(float val) {    	
    	return 2.0f*Numerics.cube(Numerics.min(1.0f-val, val)/0.5f)/(1.0f+Numerics.cube(Numerics.min(1.0f-val, val)/0.5f));	
    }

    boolean zeroNeighbor(float[] image, boolean[] mask, int x, int y, int z, int d) {
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			//if (mask[x+i+nx*(y+j)+nx*ny*(z+l)]==false && i*i+j*j+l*l<=2*d*d) return true;
			if (image[x+i+nx*(y+j)+nx*ny*(z+l)]!=image[x+nx*y+nx*ny*z] && i*i+j*j+l*l<=2*d*d) return false;
		}
		return true;
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
}
