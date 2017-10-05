package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
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

import Jama.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistShapeLevelsetPCASegmentation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 			origImage;
	private ParamVolume 			smeanImage;
	private ParamVolume 			spcaImage;
	private	ParamNumberCollection	seigenValues;
	private ParamVolume 			imeanImage;
	private ParamVolume 			ipcaImage;
	private	ParamNumberCollection	ieigenValues;
	
	//private ParamDouble ratioParam;
	private ParamDouble distanceParam;
	//private ParamBoolean cropParam;
	//private ParamString nameParam;
	//private ParamDouble scaleParam;
	
	private ParamVolume 			probaImage;
	private ParamVolume 			shapeImage;
	private ParamVolume 			intensImage;
	
	private static final byte X=0, Y=1, Z=2;
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(origImage = new ParamVolume("Image to Segment"));
		inputParams.add(smeanImage = new ParamVolume("Mean Shape Image"));
		inputParams.add(spcaImage = new ParamVolume("Shape PCA Image"));
		inputParams.add(seigenValues = new ParamNumberCollection("Shape Eigenvalues"));
		inputParams.add(imeanImage = new ParamVolume("Mean Intensity Image"));
		inputParams.add(ipcaImage = new ParamVolume("Intensity PCA Image"));
		inputParams.add(ieigenValues = new ParamNumberCollection("Intensity Eigenvalues"));
		
		inputParams.add(distanceParam = new ParamDouble("Distance to levelset (mm)", 0.0f, 100.0f, 5.0f));
		/*
		inputParams.add(ratioParam = new ParamDouble("Ratio of first to last kept eigenvalue", 0.0f, 1.0f, 0.1f));
		inputParams.add(cropParam = new ParamBoolean("Output cropped atlas", false));
		inputParams.add(nameParam = new ParamString("Atlas name"));
		inputParams.add(scaleParam = new ParamDouble("Eigenmode visualization scaling (x stdev)", 0.0f, 10.0f, 1.0f));
		*/	
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Shape.devel");
		inputParams.setLabel("Levelset PCA Segmentation");
		inputParams.setName("LevelsetPCASegmentation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compares an active shape levelset model to segmentations.");
		
		info.setVersion("1.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(probaImage = new ParamVolume("Probability Image",VoxelType.FLOAT));
		outputParams.add(shapeImage = new ParamVolume("Shape Image",VoxelType.FLOAT));
		outputParams.add(intensImage = new ParamVolume("Intensity Image",VoxelType.FLOAT));
		
		outputParams.setName("pca seg images");
		outputParams.setLabel("pca seg images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		System.out.println("get all images");
		ImageDataFloat	origImg = new ImageDataFloat(origImage.getImageData());
		
		int nx = origImg.getRows();
		int ny = origImg.getCols();
		int nz = origImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = origImg.getHeader().getDimResolutions()[0];
		float ry = origImg.getHeader().getDimResolutions()[1];
		float rz = origImg.getHeader().getDimResolutions()[2];
		
		float[][][] intensity= origImg.toArray3d();
		
		ImageDataFloat	smeanImg = new ImageDataFloat(smeanImage.getImageData());
		float[][][] smean = smeanImg.toArray3d();
		
		ImageDataFloat	spcaImg = new ImageDataFloat(spcaImage.getImageData());
		float[][][][] spca = spcaImg.toArray4d();
		
		int nsp = spcaImg.getComponents();
		System.out.println("Shape PCA components: "+nsp);
		float[] seigenval = new float[nsp];
		for (int n=0;n<nsp;n++) {
			seigenval[n] = seigenValues.getValue(n).floatValue();
		}

		ImageDataFloat	imeanImg = new ImageDataFloat(imeanImage.getImageData());
		float[][][] imean = imeanImg.toArray3d();
		
		ImageDataFloat	ipcaImg = new ImageDataFloat(ipcaImage.getImageData());
		float[][][][] ipca = ipcaImg.toArray4d();
		
		int nip = ipcaImg.getComponents();
		System.out.println("Intensity PCA components: "+nip);
		float[] ieigenval = new float[nip];
		for (int n=0;n<nip;n++) {
			ieigenval[n] = ieigenValues.getValue(n).floatValue();
		}

		// generate mask far away from shape
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = false;
			if (smean[x][y][z]<20.0f) mask[x][y][z] = true;
		}
		
		float bdist = 0.5f*distanceParam.getValue().floatValue()/rx;
		// distance for projection inside: make sure it's not bigger than the structure
		float minval = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			if (smean[x][y][z]<minval) minval = smean[x][y][z];
		}
		minval = Numerics.min(bdist, -0.5f*minval);
		System.out.println("min distance inside (voxels): "+minval);
		
		// find mean intensity inside, outside
		float meanIin = 0.0f, invol = 0.0f;
		float meanIout = 0.0f, outvol = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			if (smean[x][y][z]<=0) {
				meanIin += imean[x][y][z];
				invol++;
			} else if (smean[x][y][z]<bdist) {
				meanIout += imean[x][y][z];
				outvol++;
			}
		}
		meanIin /= invol;
		meanIout /= outvol;
		System.out.println("mean intensity inside: "+meanIin);
		System.out.println("mean intensity outside: "+meanIout);
		
		// first step: mean shape and intensity init
		float[][][] proba = new float[nx][ny][nz];
		double[] dproba = new double[nsp];
		double sumproba = 0.0;
		for (int n=0;n<nsp;n++) dproba[n] = 0.0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			
			// compute the shape distance
			float pshape = 1.0f/(1.0f + (float)FastMath.exp( smean[x][y][z] / bdist ));
		
			// compute the intensity similarity
			float pintens = (float)FastMath.exp( -0.5f*Numerics.square( (intensity[x][y][z]-imean[x][y][z])/(meanIin-meanIout) ) );
			
			//proba[x][y][z] = (float)FastMath.sqrt(pshape*pintens);
			//proba[x][y][z] =  0.5f+pintens*(pshape-0.5f);
			//proba[x][y][z] = pintens;
			//if (pshape>0.5f) proba[x][y][z] = (2.0f*(pshape+pintens)-1.0f)/3.0f;
			//else proba[x][y][z] = 2.0f*(pshape-pintens+1.0f)/3.0f;
			//proba[x][y][z] = Numerics.min(pshape*pintens, 1.0f-(1.0f-pshape)*pintens);
			proba[x][y][z] = (float)FastMath.sqrt(pshape*pintens*(1.0f-(1.0f-pshape)*pintens));
			
			// compute the derivatives at the same time
			for (int n=0;n<nsp;n++) {
				dproba[n] += bdist*seigenval[n]*spca[x][y][z][n]*pshape*(1.0-pshape)*pintens*(pintens-1.0-2.0*pintens*pshape)
							+(intensity[x][y][z]-imean[x][y][z])/Numerics.square(meanIin-meanIout)*ieigenval[n]*ipca[x][y][z][n]
							*(1.0-2.0*pintens+4.0*pintens*pshape);
			}
			sumproba += proba[x][y][z];
		}
		System.out.println("proba score: "+sumproba);
			
		// iterate?
		float ratio = 0.1f;
		float[][][] shape = new float[nx][ny][nz];
		float[][][] intens = new float[nx][ny][nz];
		float[] alpha = new float[nsp];
		for (int p=0;p<nsp;p++) {
			System.out.println("PCA dimension: "+(p+1));
			for (int t=0;t<10;t++) {
				// 1. adjust the PCA shape, intensity to minimize distance
				alpha[p] = (float)(sumproba/Numerics.max(1e-6, dproba[p]));
				
				double prevproba = sumproba;
				for (int m=0;m<10 && prevproba>=sumproba;m++) {
					System.out.print("alpha: "+alpha[p]);
				
					// 2. recompute the shape
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
						shape[x][y][z] = smean[x][y][z];
						for (int q=0;q<=p;q++)
							shape[x][y][z] += alpha[q]*seigenval[q]*spca[x][y][z][q];
					}
					
					// 3. map the intensity
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
						float[] pt = projectOntoLevelset(shape, smean[x][y][z], x, y, z, nx, ny, nz);
						
						intens[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(imean, pt[X], pt[Y], pt[Z], nx, ny, nz);
						for (int q=0;q<=p;q++)
							intens[x][y][z] += alpha[q]*ieigenval[q]*ImageInterpolation.nearestNeighborClosestInterpolation(ipca, pt[X], pt[Y], pt[Z], q, nx, ny, nz, nsp);
					}
					
					// 4. compute the new proba, check for progress
					sumproba = 0.0;
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
						// compute the shape distance
						float pshape = 1.0f/(1.0f + (float)FastMath.exp( shape[x][y][z] / bdist ));
					
						// map the intensity
						// compute the intensity similarity
						float pintens = (float)FastMath.exp( -0.5f*Numerics.square( (intensity[x][y][z]-intens[x][y][z])/(meanIin-meanIout) ) );
					
						proba[x][y][z] = (float)FastMath.sqrt(pshape*pintens*(1.0f-(1.0f-pshape)*pintens));
	
						sumproba += proba[x][y][z];	
					}
					System.out.println(", proba score: "+sumproba);
					alpha[p] *= 0.5f;
				}
			
				// 5. compute the final proba derivatives
				for (int n=0;n<nsp;n++) dproba[n] = 0.0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
					// compute the shape distance
					float pshape = 1.0f/(1.0f + (float)FastMath.exp( shape[x][y][z] / bdist ));
				
					// map the intensity
					// compute the intensity similarity
					float pintens = (float)FastMath.exp( -0.5f*Numerics.square( (intensity[x][y][z]-intens[x][y][z])/(meanIin-meanIout) ) );
				
					// compute the derivatives at the same time
					for (int n=0;n<nsp;n++) {
						dproba[n] += bdist*seigenval[n]*spca[x][y][z][n]*pshape*(1.0-pshape)*pintens*(pintens-1.0-2.0*pintens*pshape)
									+(intensity[x][y][z]-imean[x][y][z])/Numerics.square(meanIin-meanIout)*ieigenval[n]*ipca[x][y][z][n]
									*(1.0-2.0*pintens+4.0*pintens*pshape);
					}
				}
				System.out.println("difference: "+(2.0*(sumproba-prevproba)/(sumproba+prevproba)));
				if (Numerics.abs(2.0*(sumproba-prevproba)/(sumproba+prevproba))<0.01) t=10;
			}
		}
		ImageDataFloat probaData = new ImageDataFloat(proba);		
		probaData.setHeader(origImage.getImageData().getHeader());
		probaData.setName(origImage.getImageData().getName()+"_proba");
		probaImage.setValue(probaData);
		probaData = null;
		proba = null;
		
		ImageDataFloat shapeData = new ImageDataFloat(shape);		
		shapeData.setHeader(origImage.getImageData().getHeader());
		shapeData.setName(origImage.getImageData().getName()+"_shape");
		shapeImage.setValue(shapeData);
		shapeData = null;
		shape = null;
		
		ImageDataFloat intensData = new ImageDataFloat(intens);		
		intensData.setHeader(origImage.getImageData().getHeader());
		intensData.setName(origImage.getImageData().getName()+"_intens");
		intensImage.setValue(intensData);
		intensData = null;
		intens = null;
	}

	private final float[] projectOntoLevelset(float[][][] levelset, float level, int x, int y, int z, int nx, int ny, int nz) {
		
		float[] pt = new float[3];
		
		float err = level - levelset[x][y][z];
		float ratio = 1.0f;
		
		pt[X] = x; pt[Y] = y; pt[Z] = z;
		
		for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
			float err0 = err;
			// get levelset gradient at init location
			float dX = ImageInterpolation.linearInterpolationXderivative(levelset, pt[X], pt[Y], pt[Z], nx, ny, nz);
			float dY = ImageInterpolation.linearInterpolationYderivative(levelset, pt[X], pt[Y], pt[Z], nx, ny, nz);
			float dZ = ImageInterpolation.linearInterpolationZderivative(levelset, pt[X], pt[Y], pt[Z], nx, ny, nz);
			float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
		
			pt[X] += ratio*err/norm2*dX;
			pt[Y] += ratio*err/norm2*dY;
			pt[Z] += ratio*err/norm2*dZ;
				
			err = level - ImageInterpolation.linearInterpolation(levelset,levelset[x][y][z], pt[X], pt[Y], pt[Z], nx, ny, nz);
			if (Numerics.abs(err)>Numerics.abs(err0)) {
				pt[X] -= ratio*err/norm2*dX;
				pt[Y] -= ratio*err/norm2*dY;
				pt[Z] -= ratio*err/norm2*dZ;
				ratio *= 0.67f;
				err = err0;
			} else {
				ratio = 1.0f;
			}
		}
		//if (err>0.01f) System.out.print("!");
		return pt;
	}
}
