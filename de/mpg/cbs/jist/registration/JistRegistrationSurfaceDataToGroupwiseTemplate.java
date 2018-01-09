package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
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
import de.mpg.cbs.core.registration.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistRegistrationSurfaceDataToGroupwiseTemplate extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume sourceContrastImage;
	private ParamVolume sourceMaskImage;
	private ParamVolume sourceMappingImage;
	private ParamVolume sourceLevelsetImage;
	private ParamVolume templateMappingImage;
	private ParamVolume templateLevelsetImage;
	
	private ParamOption interpParam;
	private ParamOption methodParam;
		
	private ParamVolume mappedDataImage;
	private ParamVolume mappedMaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private RegistrationSurfaceDataToGroupwiseTemplate algorithm;
		
	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(sourceContrastImage = new ParamVolume("Source Contrast Image",null,-1,-1,-1,-1));
		inputParams.add(sourceMaskImage = new ParamVolume("Source Contrast Mask (opt)",null,-1,-1,-1,-1));
		sourceMaskImage.setMandatory(false);
		inputParams.add(sourceLevelsetImage = new ParamVolume("Source Surface Image",null,-1,-1,-1,-1));
		inputParams.add(sourceMappingImage = new ParamVolume("Contrast To Template Mapping (opt)",null,-1,-1,-1,-1));
		sourceMappingImage.setMandatory(false);
		inputParams.add(templateLevelsetImage = new ParamVolume("Template Surface Image",null,-1,-1,-1,-1));
		inputParams.add(templateMappingImage = new ParamVolume("Template To Canonical Space Mapping (opt)",null,-1,-1,-1,-1));
		templateMappingImage.setMandatory(false);
		
		inputParams.add(methodParam = new ParamOption("Surface Mapping Method", RegistrationSurfaceDataToGroupwiseTemplate.mappingTypes));
		inputParams.add(interpParam = new ParamOption("Interpolation", RegistrationSurfaceDataToGroupwiseTemplate.interpTypes));
			
		algorithm = new RegistrationSurfaceDataToGroupwiseTemplate();
		
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
		outputParams.add(mappedDataImage = new ParamVolume("Cortex-mapped contrast image",null,-1,-1,-1,-1));
		outputParams.add(mappedMaskImage = new ParamVolume("Cortex-mapped contrast mask",null,-1,-1,-1,-1));
		
		outputParams.setName("registered images");
		outputParams.setLabel("registered images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(sourceContrastImage);
		ImageHeader header = null;
		if (Interface.isValid(templateMappingImage)) header = Interface.getHeader(templateMappingImage);
		else header = Interface.getHeader(templateLevelsetImage);
		
		int[] sdims = Interface.getDimensions4D(sourceContrastImage);
		float[] sres = Interface.getResolutions(sourceContrastImage);
		
		int[] tdims = Interface.getDimensions(templateLevelsetImage);
		float[] tres = Interface.getResolutions(templateLevelsetImage);
		
		int[] odims = null;
		float[] ores = null;
		if (Interface.isValid(templateMappingImage)) {
			odims = Interface.getDimensions(templateMappingImage);
			ores = Interface.getResolutions(templateMappingImage);
		} else {
			odims = Interface.getDimensions(templateLevelsetImage);
			ores = Interface.getResolutions(templateLevelsetImage);
		}		
		
		// main algorithm
		algorithm = new RegistrationSurfaceDataToGroupwiseTemplate();
		
		if (Interface.isImage4D(sourceContrastImage)) 
			algorithm.setSourceContrastImage(Interface.getFloatImage4D(sourceContrastImage));
		else
			algorithm.setSourceContrastImage(Interface.getFloatImage3D(sourceContrastImage));
		if (Interface.isValid(sourceMaskImage)) 
			algorithm.setSourceMaskImage(Interface.getIntegerImage3D(sourceMaskImage));
		
		algorithm.setSourceLevelsetImage(Interface.getFloatImage3D(sourceLevelsetImage));
		if (Interface.isValid(sourceMappingImage)) 
			algorithm.setSourceMappingImage(Interface.getFloatImage4D(sourceMappingImage));
		
		algorithm.setTemplateLevelsetImage(Interface.getFloatImage3D(templateLevelsetImage));
		if (Interface.isValid(templateMappingImage)) 
			algorithm.setTemplateMappingImage(Interface.getFloatImage4D(templateMappingImage));
		
		algorithm.setSourceDimensions(sdims);
		algorithm.setSourceResolutions(sres);

		algorithm.setTemplateDimensions(tdims);
		algorithm.setTemplateResolutions(tres);

		algorithm.setOutputDimensions(odims);
		algorithm.setOutputResolutions(ores);

		
		algorithm.setSurfaceMappingMethod(methodParam.getValue());
		algorithm.setInterpolation(interpParam.getValue());
		
		algorithm.execute();
		
		// output
		String imgname = sourceContrastImage.getImageData().getName();
		
		if (Interface.isImage4D(sourceContrastImage)) 
			Interface.setFloatImage4D(algorithm.getMappedData(), odims, sdims[3], mappedDataImage, name+"_map2grp_data", header);
		else
			Interface.setFloatImage3D(algorithm.getMappedData(), odims, mappedDataImage, name+"_map2grp_data", header);
			
		Interface.setIntegerImage3D(algorithm.getmappedDataMask(), odims, mappedMaskImage, name+"_map2grp_mask", header);
	}


}
