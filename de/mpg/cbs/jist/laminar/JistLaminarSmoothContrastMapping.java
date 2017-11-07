package de.mpg.cbs.jist.laminar;

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
import de.mpg.cbs.core.laminar.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarSmoothContrastMapping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume layersImage;
	private ParamVolume contrastImage;
	private ParamVolume mappingImage;
	private ParamVolume maskImage;
	
	private ParamOption interpParam;
	private ParamOption measureParam;
	private ParamFloat distwmParam;	
	private ParamFloat distcsfParam;	
	private ParamFloat fwhmParam;	
	private ParamBoolean smoothfirstParam;	
	
	
	private ParamVolume mappedContrastImage;
	private ParamVolume mappedMaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private LaminarSmoothContrastMapping algorithm;
		
	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(layersImage = new ParamVolume("Layer Surface Image",null,-1,-1,-1,-1));
		inputParams.add(contrastImage = new ParamVolume("Contrast Image",null,-1,-1,-1,-1));
		inputParams.add(mappingImage = new ParamVolume("Contrast To Layer Mapping (opt)",null,-1,-1,-1,-1));
		mappingImage.setMandatory(false);
		inputParams.add(maskImage = new ParamVolume("Contrast Mask (opt)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		
		inputParams.add(measureParam = new ParamOption("Measure along profile", LaminarSmoothContrastMapping.measureTypes));
		inputParams.add(distwmParam = new ParamFloat("Distance ratio to WM [0,1]", 0.0f, 1.0f, 0.0f));
		inputParams.add(distcsfParam = new ParamFloat("Distance ratio to CSF [0,1]", 0.0f, 1.0f, 0.0f));	
		inputParams.add(interpParam = new ParamOption("Interpolation", LaminarSmoothContrastMapping.interpTypes));
		inputParams.add(fwhmParam = new ParamFloat("Smoothing FWHM (mm)", 0.0f, 20.0f, 0.0f));
		inputParams.add(smoothfirstParam = new ParamBoolean("Smooth before measure", false));
			
		algorithm = new LaminarSmoothContrastMapping();
		
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
		outputParams.add(mappedContrastImage = new ParamVolume("Cortex-mapped contrast image",null,-1,-1,-1,-1));
		outputParams.add(mappedMaskImage = new ParamVolume("Cortex-mapped contrast mask",null,-1,-1,-1,-1));
		
		outputParams.setName("smooth layers images");
		outputParams.setLabel("smooth layers images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(contrastImage);
		int[] idims = Interface.getDimensions4D(contrastImage);
		float[] ires = Interface.getResolutions(contrastImage);
		
		ImageHeader header = Interface.getHeader(layersImage);
		int[] dims = Interface.getDimensions4D(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		
		// main algorithm
		algorithm = new LaminarSmoothContrastMapping();
		
		algorithm.setLayerSurfaceImage(Interface.getFloatImage4D(layersImage));
		algorithm.setIntensityImage(Interface.getFloatImage3D(contrastImage));
		if (Interface.isValid(maskImage)) 
			algorithm.setIntensityMaskImage(Interface.getIntegerImage3D(maskImage));
		if (Interface.isValid(mappingImage)) 
			algorithm.setContrastToLayersMapping(Interface.getFloatImage4D(mappingImage));
		
		algorithm.setLayersDimensions(dims);
		algorithm.setLayersResolutions(res);

		algorithm.setContrastDimensions(idims);
		algorithm.setContrastResolutions(ires);

		
		algorithm.setMeasureAlongProfile(measureParam.getValue());
		algorithm.setDistanceRatioToWM(distwmParam.getValue().floatValue());
		algorithm.setDistanceRatioToCSF(distcsfParam.getValue().floatValue());
		algorithm.setInterpolation(interpParam.getValue());
		algorithm.setSmoothingFWHMmm(fwhmParam.getValue().floatValue());
		algorithm.setSmoothBeforeMeasure(smoothfirstParam.getValue().booleanValue());
		
		algorithm.execute();
		
		// output
		String imgname = contrastImage.getImageData().getName();
		
		Interface.setFloatImage4D(algorithm.getCortexMappedImage(), dims, idims[3], mappedContrastImage, name+"_lscm_data", header);
		Interface.setIntegerImage3D(algorithm.getCortexMappedMask(), dims, mappedMaskImage, name+"_lscm_mask", header);
	}


}
