package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
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

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.algorithms.AlgorithmWSinc;

/*
 * @author Pierre-Louis Bazin
 */
public class JistRegistrationApplyDeformations extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume sourceImage;
	private ParamVolume referenceImage;
	private ParamVolume deformation1Image;
	private ParamVolume deformation2Image;
	private ParamOption type1Option;
	private ParamOption type2Option;
	private ParamOption interpOption;
	private ParamOption padOption;
	
	private static final String[] types1 = {"none", "deformation(voxels)", "mapping(voxels)", "deformation(mm)", "mapping(mm)"};
	private static final String[] types2 = {"none", "deformation(voxels)", "mapping(voxels)", "deformation(mm)", "mapping(mm)"};
	private static final String[] interp = {"WSinc","linear", "NN"};
	private static final String[] pads = {"closest", "zero", "min", "max"};
	
	private ParamVolume mappedImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	private static final byte T = 3;

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(sourceImage = new ParamVolume("Image to deform",null,-1,-1,-1,-1));
		inputParams.add(referenceImage = new ParamVolume("Reference Image",null,-1,-1,-1,-1));
		inputParams.add(deformation1Image = new ParamVolume("Deformation 1",null,-1,-1,-1,-1));
		inputParams.add(type1Option = new ParamOption("Deformation type 1",types1));
		inputParams.add(deformation2Image = new ParamVolume("Deformation 2",null,-1,-1,-1,-1));
		inputParams.add(type2Option = new ParamOption("Deformation type 2",types2));
		inputParams.add(interpOption = new ParamOption("Interpolation type",interp));
		inputParams.add(padOption = new ParamOption("Image padding",pads));
		deformation2Image.setMandatory(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		inputParams.setLabel("Apply Deformations");
		inputParams.setName("ApplyDeformations");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Apply multiple deformations from non-linear transformations. Transformations can be a deformation field or a coordinate mapping.");
		
		info.setVersion("3.0.1");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mappedImage = new ParamVolume("Deformed Image",null,-1,-1,-1,-1));
		
		outputParams.setName("deformed images");
		outputParams.setLabel("deformed images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	sourceImg = new ImageDataFloat(sourceImage.getImageData());
		ImageDataFloat	referenceImg = new ImageDataFloat(referenceImage.getImageData());
		
		int nsx = sourceImg.getRows();
		int nsy = sourceImg.getCols();
		int nsz = sourceImg.getSlices();
		int nst = sourceImg.getComponents();
		int nsxyz = nsx*nsy*nsz;
		float rsx = sourceImg.getHeader().getDimResolutions()[X];
		float rsy = sourceImg.getHeader().getDimResolutions()[Y];
		float rsz = sourceImg.getHeader().getDimResolutions()[Z];
		float osx = sourceImg.getHeader().getOrigin()[X];
		float osy = sourceImg.getHeader().getOrigin()[Y];
		float osz = sourceImg.getHeader().getOrigin()[Z];
		
		int nrx = referenceImg.getRows();
		int nry = referenceImg.getCols();
		int nrz = referenceImg.getSlices();
		int nrt = referenceImg.getComponents();
		int nrxyz = nrx*nry*nrz;
		float rrx = referenceImg.getHeader().getDimResolutions()[X];
		float rry = referenceImg.getHeader().getDimResolutions()[Y];
		float rrz = referenceImg.getHeader().getDimResolutions()[Z];
		float orx = referenceImg.getHeader().getOrigin()[X];
		float ory = referenceImg.getHeader().getOrigin()[Y];
		float orz = referenceImg.getHeader().getOrigin()[Z];
		
		// different paths for 3d and 4d data
		if (nst==1) {
			// pull the image data (not needed for reference)
			float[][][] source = sourceImg.toArray3d();
			
			// deformation: in reference space
			System.out.println("load deformation 1");
			ImageDataFloat deformation1Img = new ImageDataFloat(deformation1Image.getImageData());
			int nd1x = deformation1Img.getRows();
			int nd1y = deformation1Img.getCols();
			int nd1z = deformation1Img.getSlices();
			float rd1x = referenceImg.getHeader().getDimResolutions()[X];
			float rd1y = referenceImg.getHeader().getDimResolutions()[Y];
			float rd1z = referenceImg.getHeader().getDimResolutions()[Z];
			float[][][][] deformation1 = deformation1Img.toArray4d();
			// scale to voxels if needed
			if (type1Option.getValue().endsWith("(mm)")) {
				System.out.println("normalize to resolution ("+rd1x+", "+rd1y+", "+rd1z+")");
				for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
					deformation1[x][y][z][X] /= rd1x; 	
					deformation1[x][y][z][Y] /= rd1y; 	
					deformation1[x][y][z][Z] /= rd1z; 	
				}
			}
			// turn into a mapping if needed
			if (type1Option.getValue().startsWith("deformation")) {
				for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
					deformation1[x][y][z][X] += x;
					deformation1[x][y][z][Y] += y;
					deformation1[x][y][z][Z] += z;
				}
			}
			
			// second deformation if given
			if (!type2Option.getValue().equals("none") && deformation2Image.getImageData()!=null) {
				System.out.println("load deformation 2");
				ImageDataFloat deformation2Img = new ImageDataFloat(deformation2Image.getImageData());
				int nd2x = deformation2Img.getRows();
				int nd2y = deformation2Img.getCols();
				int nd2z = deformation2Img.getSlices();
				float rd2x = referenceImg.getHeader().getDimResolutions()[X];
				float rd2y = referenceImg.getHeader().getDimResolutions()[Y];
				float rd2z = referenceImg.getHeader().getDimResolutions()[Z];
				float[][][][] deformation2 = deformation2Img.toArray4d();
				// scale to voxels if needed
				if (type2Option.getValue().endsWith("(mm)")) {
					System.out.println("normalize to resolution ("+rd2x+", "+rd2y+", "+rd2z+")");
					for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
						deformation2[x][y][z][X] /= rd2x; 	
						deformation2[x][y][z][Y] /= rd2y; 	
						deformation2[x][y][z][Z] /= rd2z; 	
					}
				}
				// turn into a mapping if needed
				if (type2Option.getValue().startsWith("deformation")) {
					for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
						deformation2[x][y][z][X] += x;
						deformation2[x][y][z][Y] += y;
						deformation2[x][y][z][Z] += z;
					}
				}
				// compose the deformations: X' = def1(def2(X))
				System.out.println("compose deformations");
				float[][][][] composed = new float[nrx][nry][nrz][3];
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
						composed[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], X, nd1x, nd1y, nd1z, 3);
						composed[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Y, nd1x, nd1y, nd1z, 3);
						composed[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Z, nd1x, nd1y, nd1z, 3);
					//}
				}
				deformation1 = composed;
				deformation2 = null;
			}
			// new image
			System.out.println("deform image");
			float min = 1e10f, max = -1e10f;
			if (padOption.getValue().equals("min") || padOption.getValue().equals("max")) {
				for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
					if (source[x][y][z]<min) min = source[x][y][z];
					if (source[x][y][z]>max) max = source[x][y][z];
				}
			}
			float[][][] deformed = new float[nrx][nry][nrz];
			if (interpOption.getValue().equals("NN")) {
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
						if (padOption.getValue().equals("closest"))
							deformed[x][y][z] = ImageInterpolation.nearestNeighborClosestInterpolation(source, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("zero"))
							deformed[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(source, 0.0f, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("min"))
							deformed[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(source, min, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("max"))
							deformed[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(source, max, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
					//}
				}				
			} else if (interpOption.getValue().equals("linear")) {
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
						if (padOption.getValue().equals("closest"))
							deformed[x][y][z] = ImageInterpolation.linearClosestInterpolation(source, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("zero"))
							deformed[x][y][z] = ImageInterpolation.linearInterpolation(source, 0.0f, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("min"))
							deformed[x][y][z] = ImageInterpolation.linearInterpolation(source, min, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
						else if (padOption.getValue().equals("max"))
							deformed[x][y][z] = ImageInterpolation.linearInterpolation(source, max, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], nsx, nsy, nsz);
					//}
				}
			} else if (interpOption.getValue().equals("WSinc")) {
				double[] imgBuf = new double[nsx*nsy*nsz];
				for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
					imgBuf[x+nsx*y+nsx*nsy*z] = source[x][y][z];
				}
				AlgorithmWSinc WSinc = new AlgorithmWSinc();
				final int[] inVolExtents = {nsx, nsy, nsz};
				WSinc.setup3DWSinc(imgBuf, inVolExtents, true);
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
						if (deformation1[x][y][z][X]>=0 && deformation1[x][y][z][X]<nsx
							&& deformation1[x][y][z][Y]>=0 && deformation1[x][y][z][Y]<nsy
							&& deformation1[x][y][z][Z]>=0 && deformation1[x][y][z][Z]<nsz) {
							deformed[x][y][z] = (float)WSinc.wSinc3D(deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z]);
						} else {
							if (padOption.getValue().equals("zero")) deformed[x][y][z] = 0.0f;
							else if (padOption.getValue().equals("min")) deformed[x][y][z] = min;
							else if (padOption.getValue().equals("max")) deformed[x][y][z] = max;
						}
					//}
				}
			}
			// output
			String imgname = sourceImage.getImageData().getName();
			
			ImageDataFloat mapData = new ImageDataFloat(deformed);		
			mapData.setHeader(referenceImage.getImageData().getHeader());
			mapData.setName(imgname+"_def");
			mappedImage.setValue(mapData);
			mapData = null;
			deformed = null;
		} else {
			// pull the image data (not needed for reference)
			float[][][][] source = sourceImg.toArray4d();
			
			// deformation: in reference space
			System.out.println("load deformation 1");
			ImageDataFloat deformation1Img = new ImageDataFloat(deformation1Image.getImageData());
			int nd1x = deformation1Img.getRows();
			int nd1y = deformation1Img.getCols();
			int nd1z = deformation1Img.getSlices();
			float rd1x = referenceImg.getHeader().getDimResolutions()[X];
			float rd1y = referenceImg.getHeader().getDimResolutions()[Y];
			float rd1z = referenceImg.getHeader().getDimResolutions()[Z];
			float[][][][] deformation1 = deformation1Img.toArray4d();
			// scale to voxels if needed
			if (type1Option.getValue().endsWith("(mm)")) {
				for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
					deformation1[x][y][z][X] /= rd1x; 	
					deformation1[x][y][z][Y] /= rd1y; 	
					deformation1[x][y][z][Z] /= rd1z; 	
				}
			}
			// turn into a mapping if needed
			if (type1Option.getValue().startsWith("deformation")) {
				for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) for (int z=0;z<nd1z;z++) {
					deformation1[x][y][z][X] += x;
					deformation1[x][y][z][Y] += y;
					deformation1[x][y][z][Z] += z;
				}
			}
			
			// second deformation if given
			if (!type2Option.getValue().equals("none") && deformation2Image.getImageData()!=null) {
				ImageDataFloat deformation2Img = new ImageDataFloat(deformation2Image.getImageData());
				int nd2x = deformation2Img.getRows();
				int nd2y = deformation2Img.getCols();
				int nd2z = deformation2Img.getSlices();
				float rd2x = referenceImg.getHeader().getDimResolutions()[X];
				float rd2y = referenceImg.getHeader().getDimResolutions()[Y];
				float rd2z = referenceImg.getHeader().getDimResolutions()[Z];
				float[][][][] deformation2 = deformation2Img.toArray4d();
				// scale to voxels if needed
				if (type2Option.getValue().endsWith("(mm)")) {
					for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
						deformation2[x][y][z][X] /= rd2x; 	
						deformation2[x][y][z][Y] /= rd2y; 	
						deformation2[x][y][z][Z] /= rd2z; 	
					}
				}
				// turn into a mapping if needed
				if (type2Option.getValue().startsWith("deformation")) {
					for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
						deformation2[x][y][z][X] += x;
						deformation2[x][y][z][Y] += y;
						deformation2[x][y][z][Z] += z;
					}
				}
				// compose the deformations: X' = def1(def2(X))
				float[][][][] composed = new float[nrx][nry][nrz][3];
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
						composed[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], X, nd1x, nd1y, nd1z, 3);
						composed[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Y, nd1x, nd1y, nd1z, 3);
						composed[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Z, nd1x, nd1y, nd1z, 3);
					//}
				}
				deformation1 = composed;
				deformation2 = null;
			}
			// new image
			System.out.println("deform image");
			float min = 1e10f, max = -1e10f;
			if (padOption.getValue().equals("min") || padOption.getValue().equals("max")) {
				for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int t=0;t<nst;t++) {
					if (source[x][y][z][t]<min) min = source[x][y][z][t];
					if (source[x][y][z][t]>max) max = source[x][y][z][t];
				}
			}
			float[][][][] deformed = new float[nrx][nry][nrz][nst];
			if (interpOption.getValue().equals("NN")) {
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
						for (int t=0;t<nst;t++) {
							if (padOption.getValue().equals("closest"))
								deformed[x][y][z][t] = ImageInterpolation.nearestNeighborClosestInterpolation(source, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("zero"))
								deformed[x][y][z][t] = ImageInterpolation.nearestNeighborInterpolation(source, 0.0f, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("min"))
								deformed[x][y][z][t] = ImageInterpolation.nearestNeighborInterpolation(source, min, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("max"))
								deformed[x][y][z][t] = ImageInterpolation.nearestNeighborInterpolation(source, max, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst); 
						}
					//}
				}
			} else if (interpOption.getValue().equals("linear")) {
				for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
					//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
						for (int t=0;t<nst;t++) {
							if (padOption.getValue().equals("closest"))
								deformed[x][y][z][t] = ImageInterpolation.linearClosestInterpolation(source, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("zero"))
								deformed[x][y][z][t] = ImageInterpolation.linearInterpolation(source, 0.0f, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("min"))
								deformed[x][y][z][t] = ImageInterpolation.linearInterpolation(source, min, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
							else if (padOption.getValue().equals("max"))
								deformed[x][y][z][t] = ImageInterpolation.linearInterpolation(source, max, deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z], t, nsx, nsy, nsz, nst);
						}
					//}					
				}
			} else if (interpOption.getValue().equals("WSinc")) {
				double[] imgBuf = new double[nsx*nsy*nsz];
				AlgorithmWSinc WSinc = new AlgorithmWSinc();
				final int[] inVolExtents = {nsx, nsy, nsz};
				for (int t=0;t<nst;t++) {
					for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
						imgBuf[x+nsx*y+nsx*nsy*z] = source[x][y][z][t];
					}
					WSinc.setup3DWSinc(imgBuf, inVolExtents, true);
					for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) for (int z=0;z<nrz;z++) {
						//if (deformation1[x][y][z][X]>=1 && deformation1[x][y][z][Y]>=1 && deformation1[x][y][z][Z]>=1) {
							if (deformation1[x][y][z][X]>=0 && deformation1[x][y][z][X]<nsx
								&& deformation1[x][y][z][Y]>=0 && deformation1[x][y][z][Y]<nsy
								&& deformation1[x][y][z][Z]>=0 && deformation1[x][y][z][Z]<nsz) {
								deformed[x][y][z][t] = (float)WSinc.wSinc3D(deformation1[x][y][z][X], deformation1[x][y][z][Y], deformation1[x][y][z][Z]);
							} else {
								if (padOption.getValue().equals("zero")) deformed[x][y][z][t] = 0.0f;
								else if (padOption.getValue().equals("min")) deformed[x][y][z][t] = min;
								else if (padOption.getValue().equals("max")) deformed[x][y][z][t] = max;
							}
						//}
					}
				}
			}
			// output
			String imgname = sourceImage.getImageData().getName();
			
			ImageDataFloat mapData = new ImageDataFloat(deformed);		
			mapData.setHeader(referenceImage.getImageData().getHeader());
			mapData.setName(imgname+"_def");
			mappedImage.setValue(mapData);
			mapData = null;
			deformed = null;
		}

		
	}
	
	public final String displayTransform(float[][] trans) {
		String info = "transform: \n";
		for (int m=0;m<3;m++) {
			info+="(";
			for (int n=0;n<3;n++) info += trans[m][n]+", ";
			info += trans[m][3]+")\n";
		}
		return info;
	}
	
	/** display the atlas data */
	public final String displayVector(float[] vect) {
		String info = "vector: (";
		for (int n=0;n<vect.length-1;n++) info += vect[n]+", ";
		info += vect[vect.length-1]+")\n";
		
		return info;
	}
	

}
