package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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

//import gov.nih.mipav.model.algorithms.AlgorithmWSinc;
import de.mpg.cbs.core.registration.RegistrationApplyDeformations;

/*
 * @author Pierre-Louis Bazin
 */
public class JistRegistrationApplyDeformations extends ProcessingAlgorithm {

    private RegistrationApplyDeformation algorithm;
    
	// jist containers
	private ParamVolume sourceImage;
	private ParamVolume referenceImage;
	private ParamVolume deformation1Image;
	private ParamVolume deformation2Image;
	private ParamVolume deformation3Image;
	private ParamOption type1Option;
	private ParamOption type2Option;
	private ParamOption type3Option;
	private ParamOption interpOption;
	private ParamOption padOption;
	
	//private static final String[] types = {"none", "deformation(voxels)", "mapping(voxels)", "deformation(mm)", "mapping(mm)"};
	//private static final String[] interp = {"NN", "linear", "WSinc"};
	//private static final String[] pads = {"closest", "zero", "min", "max"};
	
	private ParamVolume deformedImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	private static final byte T = 3;

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(sourceImage = new ParamVolume("Image to deform",null,-1,-1,-1,-1));
		inputParams.add(referenceImage = new ParamVolume("Reference Image",null,-1,-1,-1,-1));
		inputParams.add(deformation1Image = new ParamVolume("Deformation 1",null,-1,-1,-1,-1));
		inputParams.add(type1Option = new ParamOption("Deformation type 1",RegistrationApplyDeformations.types));
		type1Option.setValue("mapping(voxels)");
		inputParams.add(deformation2Image = new ParamVolume("Deformation 2",null,-1,-1,-1,-1));
		inputParams.add(type2Option = new ParamOption("Deformation type 2",RegistrationApplyDeformations.types));
		deformation2Image.setMandatory(false);
		type2Option.setValue("none");
		inputParams.add(deformation3Image = new ParamVolume("Deformation 3",null,-1,-1,-1,-1));
		inputParams.add(type3Option = new ParamOption("Deformation type 3",RegistrationApplyDeformations.types));
		deformation3Image.setMandatory(false);
		type3Option.setValue("none");
		inputParams.add(interpOption = new ParamOption("Interpolation type",RegistrationApplyDeformations.interp));
		inputParams.add(padOption = new ParamOption("Image padding",RegistrationApplyDeformations.pads));
		
		sourceImage.setLoadAndSaveOnValidate(false);
		referenceImage.setLoadAndSaveOnValidate(false);
		deformation1Image.setLoadAndSaveOnValidate(false);
		deformation2Image.setLoadAndSaveOnValidate(false);
		
		algorithm = new RegistrationApplyDeformations();
		
		inputParams.setPackage(algorithm.getPackage());
		inputParams.setCategory(algorithm.getCategory());
		inputParams.setLabel(algorithm.getLabel());
		inputParams.setName(algorithm.getName());

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(References.getAuthor(algorithm.getAlgorithmAuthors()[0]));
		info.setAffiliation(algorithm.getAffiliation());
		info.setDescription(algorithm.getDescription());
		
		info.setVersion(algorithm.getVersion());
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(deformedImage = new ParamVolume("Deformed Image",null,-1,-1,-1,-1));
		
		outputParams.setName("deformed images");
		outputParams.setLabel("deformed images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions4D(sourceImage);
		float[] res = Interface.getResolutions(sourceImage);
		String name = Interface.getName(sourceImage);
		ImageHeader header = Interface.getHeader(sourceImage);
		
		// main algorithm
		algorithm = new RegistrationApplyDeformations();
		
		if (dims[T]>1) algorithm.setSourceImage(Interface.getFloatImage4D(sourceImage));
		else algorithm.setSourceImage(Interface.getFloatImage3D(sourceImage));
		algorithm.setImageDimensions(dims);
		algorithm.setImageResolutions(res);
		
		// load the deformations
		algorithm.setDeformationMapping1(Interface.getFloatImage4D(deformation1Image));
		algorithm.setDeformationType1(type1Option.getValue()); 
		algorithm.setDeformation1Dimensions(Interface.getDimensions(deformation1Image));
		algorithm.setDeformation1Resolutions(nterface.getResolutions(deformation1Image));
		
		if (!type2Option.getValue().equals("none") && Interface.isValid(deformation2Image)) {
		    algorithm.setDeformationMapping2(Interface.getFloatImage4D(deformation2Image));
            algorithm.setDeformationType2(type2Option.getValue()); 
            algorithm.setDeformation2Dimensions(Interface.getDimensions(deformation2Image));
            algorithm.setDeformation2Resolutions(nterface.getResolutions(deformation2Image));
            
            if (!type3Option.getValue().equals("none") && Interface.isValid(deformation3Image)) {
                algorithm.setDeformationMapping3(Interface.getFloatImage4D(deformation3Image));
                algorithm.setDeformationType3(type3Option.getValue()); 
                algorithm.setDeformation3Dimensions(Interface.getDimensions(deformation3Image));
                algorithm.setDeformation3Resolutions(nterface.getResolutions(deformation3Image));
                
                
            }
               
		}
		 
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	sourceImg = new ImageDataFloat(sourceImage.getImageData());
		ImageDataFloat	referenceImg = new ImageDataFloat(referenceImage.getImageData());
		
		sourceImage.dispose();
		referenceImage.dispose();
		
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
		String imgname = sourceImage.getImageData().getName();
			
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
		ImageHeader refheader = referenceImg.getHeader();
		//referenceImg.finalize();
		referenceImg = null;
		
		// different paths for 3d and 4d data
		if (nst==1) {
			// pull the image data (not needed for reference)
			float[][][] source = sourceImg.toArray3d();
			//sourceImg.finalize();
			sourceImg = null;
			
			// deformation: in reference space
			System.out.println("load deformation 1");
			ImageDataFloat deformation1Img = new ImageDataFloat(deformation1Image.getImageData());
			deformation1Image.dispose();
			int nd1x = deformation1Img.getRows();
			int nd1y = deformation1Img.getCols();
			int nd1z = deformation1Img.getSlices();
			float rd1x = deformation1Img.getHeader().getDimResolutions()[X];
			float rd1y = deformation1Img.getHeader().getDimResolutions()[Y];
			float rd1z = deformation1Img.getHeader().getDimResolutions()[Z];
			float[][][][] deformation1 = deformation1Img.toArray4d();
			//deformation1Img.finalize();
			deformation1Img = null;
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
				deformation2Image.dispose();
				int nd2x = deformation2Img.getRows();
				int nd2y = deformation2Img.getCols();
				int nd2z = deformation2Img.getSlices();
				float rd2x = deformation2Img.getHeader().getDimResolutions()[X];
				float rd2y = deformation2Img.getHeader().getDimResolutions()[Y];
				float rd2z = deformation2Img.getHeader().getDimResolutions()[Z];
				float[][][][] deformation2 = deformation2Img.toArray4d();
				//deformation2Img.finalize();
				deformation2Img = null;
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
				float[][][][] composed12 = new float[nd2x][nd2y][nd2z][3];
				for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
					//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
						composed12[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], X, nd1x, nd1y, nd1z, 3);
						composed12[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Y, nd1x, nd1y, nd1z, 3);
						composed12[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Z, nd1x, nd1y, nd1z, 3);
					//}
				}
				deformation1 = composed12;
				deformation2 = null;

				// third deformation if given
				if (!type3Option.getValue().equals("none") && deformation3Image.getImageData()!=null) {
					System.out.println("load deformation 3");
					ImageDataFloat deformation3Img = new ImageDataFloat(deformation3Image.getImageData());
					deformation3Image.dispose();
					int nd3x = deformation3Img.getRows();
					int nd3y = deformation3Img.getCols();
					int nd3z = deformation3Img.getSlices();
					float rd3x = deformation3Img.getHeader().getDimResolutions()[X];
					float rd3y = deformation3Img.getHeader().getDimResolutions()[Y];
					float rd3z = deformation3Img.getHeader().getDimResolutions()[Z];
					float[][][][] deformation3 = deformation3Img.toArray4d();
					//deformation3Img.finalize();
					deformation3Img = null;
					// scale to voxels if needed
					if (type3Option.getValue().endsWith("(mm)")) {
						System.out.println("normalize to resolution ("+rd3x+", "+rd3y+", "+rd3z+")");
						for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
							deformation3[x][y][z][X] /= rd3x; 	
							deformation3[x][y][z][Y] /= rd3y; 	
							deformation3[x][y][z][Z] /= rd3z; 	
						}
					}
					// turn into a mapping if needed
					if (type3Option.getValue().startsWith("deformation")) {
						for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
							deformation3[x][y][z][X] += x;
							deformation3[x][y][z][Y] += y;
							deformation3[x][y][z][Z] += z;
						}
					}
					// compose the deformations: X' = def1(def2(X))
					System.out.println("compose deformations");
					float[][][][] composed23 = new float[nd3x][nd3y][nd3z][3];
					for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
						//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
							composed23[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], X, nd2x, nd2y, nd2z, 3);
							composed23[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], Y, nd2x, nd2y, nd2z, 3);
							composed23[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], Z, nd2x, nd2y, nd2z, 3);
						//}
					}
					deformation1 = composed23;
					deformation3 = null;
				}
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
			source = null;
			// output
			ImageDataFloat mapData = new ImageDataFloat(deformed);		
			mapData.setHeader(refheader);
			mapData.setName(imgname+"_def");
			mappedImage.setValue(mapData);
			mappedImage.writeAndFreeNow(this);
			mapData = null;
			deformed = null;
		} else {
			// pull the image data (not needed for reference)
			float[][][][] source = sourceImg.toArray4d();
			//sourceImg.finalize();
			sourceImg = null;
			
			// deformation: in reference space
			System.out.println("load deformation 1");
			ImageDataFloat deformation1Img = new ImageDataFloat(deformation1Image.getImageData());
			deformation1Image.dispose();
			int nd1x = deformation1Img.getRows();
			int nd1y = deformation1Img.getCols();
			int nd1z = deformation1Img.getSlices();
			float rd1x = deformation1Img.getHeader().getDimResolutions()[X];
			float rd1y = deformation1Img.getHeader().getDimResolutions()[Y];
			float rd1z = deformation1Img.getHeader().getDimResolutions()[Z];
			float[][][][] deformation1 = deformation1Img.toArray4d();
			//deformation1Img.finalize();
			deformation1Img = null;
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
				deformation2Image.dispose();
				int nd2x = deformation2Img.getRows();
				int nd2y = deformation2Img.getCols();
				int nd2z = deformation2Img.getSlices();
				float rd2x = deformation2Img.getHeader().getDimResolutions()[X];
				float rd2y = deformation2Img.getHeader().getDimResolutions()[Y];
				float rd2z = deformation2Img.getHeader().getDimResolutions()[Z];
				float[][][][] deformation2 = deformation2Img.toArray4d();
				//deformation2Img.finalize();
				deformation2Img = null;
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
				float[][][][] composed12 = new float[nrx][nry][nrz][3];
				for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) for (int z=0;z<nd2z;z++) {
					//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
						composed12[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], X, nd1x, nd1y, nd1z, 3);
						composed12[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Y, nd1x, nd1y, nd1z, 3);
						composed12[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation2[x][y][z][X], deformation2[x][y][z][Y], deformation2[x][y][z][Z], Z, nd1x, nd1y, nd1z, 3);
					//}
				}
				deformation1 = composed12;
				deformation2 = null;

				// third deformation if given
				if (!type3Option.getValue().equals("none") && deformation3Image.getImageData()!=null) {
					System.out.println("load deformation 3");
					ImageDataFloat deformation3Img = new ImageDataFloat(deformation3Image.getImageData());
					deformation3Image.dispose();
					int nd3x = deformation3Img.getRows();
					int nd3y = deformation3Img.getCols();
					int nd3z = deformation3Img.getSlices();
					float rd3x = deformation3Img.getHeader().getDimResolutions()[X];
					float rd3y = deformation3Img.getHeader().getDimResolutions()[Y];
					float rd3z = deformation3Img.getHeader().getDimResolutions()[Z];
					float[][][][] deformation3 = deformation3Img.toArray4d();
					//deformation2Img.finalize();
					deformation3Img = null;
					// scale to voxels if needed
					if (type3Option.getValue().endsWith("(mm)")) {
						System.out.println("normalize to resolution ("+rd3x+", "+rd3y+", "+rd3z+")");
						for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
							deformation3[x][y][z][X] /= rd3x; 	
							deformation3[x][y][z][Y] /= rd3y; 	
							deformation3[x][y][z][Z] /= rd3z; 	
						}
					}
					// turn into a mapping if needed
					if (type2Option.getValue().startsWith("deformation")) {
						for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
							deformation3[x][y][z][X] += x;
							deformation3[x][y][z][Y] += y;
							deformation3[x][y][z][Z] += z;
						}
					}
					// compose the deformations: X' = def1(def2(X))
					System.out.println("compose deformations");
					float[][][][] composed23 = new float[nd3x][nd3y][nd3z][3];
					for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) for (int z=0;z<nd3z;z++) {
						//if (deformation2[x][y][z][X]>=1 && deformation2[x][y][z][Y]>=1 && deformation2[x][y][z][Z]>=1) {
							composed23[x][y][z][X] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], X, nd2x, nd2y, nd2z, 3);
							composed23[x][y][z][Y] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], Y, nd2x, nd2y, nd2z, 3);
							composed23[x][y][z][Z] = ImageInterpolation.linearClosestInterpolation(deformation1, deformation3[x][y][z][X], deformation3[x][y][z][Y], deformation3[x][y][z][Z], Z, nd2x, nd2y, nd2z, 3);
						//}
					}
					deformation1 = composed23;
					deformation3 = null;
				}
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
			source = null;
			// output
			ImageDataFloat mapData = new ImageDataFloat(deformed);		
			mapData.setHeader(refheader);
			mapData.setName(imgname+"_def");
			mappedImage.setValue(mapData);
			mappedImage.writeAndFreeNow(this);
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
