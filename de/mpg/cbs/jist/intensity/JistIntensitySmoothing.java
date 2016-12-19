package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.core.intensity.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensitySmoothing extends ProcessingAlgorithm {

	private IntensitySmoothing algorithm;
	
	// jist containers
	private ParamVolume inputImage;
	private ParamVolume maskImage;
	private ParamOption methodParam;
	private ParamFloat scaleParam;
	private ParamBoolean	skip0Param;
	
	private ParamVolume smoothedImage;
	
	// parameters
	private		static final String[]	methods = IntensitySmoothing.methods;
	private		String		method = "gaussian";
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("mask Image (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(methodParam = new ParamOption("Smoothing method", methods));
		methodParam.setValue(method);
		
		inputParams.add(scaleParam = new ParamFloat("Smoothing scale (mm)", 0.0f, 100.0f, 1.0f));
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", true));
		
		algorithm = new IntensitySmoothing();
		
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
		outputParams.add(smoothedImage = new ParamVolume("Smoothed Image",VoxelType.FLOAT));
		
		outputParams.setName("smoothed images");
		outputParams.setLabel("smoothed images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions(inputImage);
		float[] res = Interface.getResolutions(inputImage);
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		
		// main algorithm
		algorithm = new IntensitySmoothing();
		
		algorithm.setInputImage(Interface.getFloatImage3D(inputImage));
		if (Interface.isValid(maskImage)) algorithm.setInputMask(Interface.getUByteImage3D(maskImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setSmoothingScale_mm(scaleParam.getValue().floatValue());
		algorithm.setSmoothingMethod(methodParam.getValue());
		algorithm.setSkipZeroValues(skip0Param.getValue());
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getSmoothedImage(), dims, smoothedImage, name+"_smoothed", header);
	}


}
