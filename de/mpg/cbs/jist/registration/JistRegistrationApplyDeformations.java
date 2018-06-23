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

    private RegistrationApplyDeformations algorithm;
    
	// jist containers
	private ParamVolume sourceImage;
	//private ParamVolume referenceImage;
	private ParamVolume deformation1Image;
	private ParamVolume deformation2Image;
	private ParamVolume deformation3Image;
	private ParamVolume deformation4Image;
	private ParamOption type1Option;
	private ParamOption type2Option;
	private ParamOption type3Option;
	private ParamOption type4Option;
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
		//inputParams.add(referenceImage = new ParamVolume("Reference Image",null,-1,-1,-1,-1));
		//referenceImage.setMandatory(false);
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
		inputParams.add(deformation4Image = new ParamVolume("Deformation 4",null,-1,-1,-1,-1));
		inputParams.add(type4Option = new ParamOption("Deformation type 4",RegistrationApplyDeformations.types));
		deformation4Image.setMandatory(false);
		type4Option.setValue("none");
		inputParams.add(interpOption = new ParamOption("Interpolation type",RegistrationApplyDeformations.interp));
		inputParams.add(padOption = new ParamOption("Image padding",RegistrationApplyDeformations.pads));
		
		sourceImage.setLoadAndSaveOnValidate(false);
		//referenceImage.setLoadAndSaveOnValidate(false);
		deformation1Image.setLoadAndSaveOnValidate(false);
		deformation2Image.setLoadAndSaveOnValidate(false);
		deformation3Image.setLoadAndSaveOnValidate(false);
		deformation4Image.setLoadAndSaveOnValidate(false);
		
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
		
		if (dims[T]>1) algorithm.setImageToDeform(Interface.getFloatImage4D(sourceImage));
		else algorithm.setImageToDeform(Interface.getFloatImage3D(sourceImage));
		algorithm.setImageDimensions(dims);
		algorithm.setImageResolutions(res);
		
		// load the deformations
		algorithm.setDeformationMapping1(Interface.getFloatImage4D(deformation1Image));
		algorithm.setDeformationType1(type1Option.getValue()); 
		algorithm.setDeformation1Dimensions(Interface.getDimensions(deformation1Image));
		algorithm.setDeformation1Resolutions(Interface.getResolutions(deformation1Image));
		
		if (!type2Option.getValue().equals("none") && Interface.isValid(deformation2Image)) {
		    algorithm.setDeformationMapping2(Interface.getFloatImage4D(deformation2Image));
            algorithm.setDeformationType2(type2Option.getValue()); 
            algorithm.setDeformation2Dimensions(Interface.getDimensions(deformation2Image));
            algorithm.setDeformation2Resolutions(Interface.getResolutions(deformation2Image));
            
            if (!type3Option.getValue().equals("none") && Interface.isValid(deformation3Image)) {
                algorithm.setDeformationMapping3(Interface.getFloatImage4D(deformation3Image));
                algorithm.setDeformationType3(type3Option.getValue()); 
                algorithm.setDeformation3Dimensions(Interface.getDimensions(deformation3Image));
                algorithm.setDeformation3Resolutions(Interface.getResolutions(deformation3Image));
                
                if (!type4Option.getValue().equals("none") && Interface.isValid(deformation4Image)) {
                    algorithm.setDeformationMapping4(Interface.getFloatImage4D(deformation4Image));
                    algorithm.setDeformationType4(type4Option.getValue()); 
                    algorithm.setDeformation4Dimensions(Interface.getDimensions(deformation4Image));
                    algorithm.setDeformation4Resolutions(Interface.getResolutions(deformation4Image));                    
                }
            }
		}
		
		// parameters
		algorithm.setInterpolationType(interpOption.getValue());
		algorithm.setImagePadding(padOption.getValue());
		 
		algorithm.execute();
		
		if (dims[T]>1) Interface.setFloatImage4D(algorithm.getDeformedImage(), dims, dims[T], deformedImage, name+"_def_img", header);
		else Interface.setFloatImage3D(algorithm.getDeformedImage(), dims, deformedImage, name+"def_img", header);
	}
}
