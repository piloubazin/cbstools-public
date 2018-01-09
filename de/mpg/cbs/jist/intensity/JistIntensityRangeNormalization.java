package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
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

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.core.intensity.IntensityRangeNormalization;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityRangeNormalization extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	normtypes = {"linear","robust","robust-min","robust-max"};
	private		String		normtype = "robust";
	
	// jist containers
	private ParamVolume inImage;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private	ParamOption	 normParam;
	private ParamFloat	 ratioParam;
	private ParamBoolean ignoreNegParam;
	private ParamBoolean ignoreZeroParam;
	private ParamFloat	 scalingParam;
	
	private IntensityRangeNormalization algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("Mask Image"));
		maskImage.setMandatory(false);
		inputParams.add(normParam = new ParamOption("normalization type", normtypes));
		normParam.setDescription("linear: use the whole range [Imin,Imax]; robust: use a truncated version of [Imin,Imax]; robust-max: use a truncated version of [0,Imax]");
		normParam.setValue(normtype);
		inputParams.add(ratioParam = new ParamFloat("Robustness ratio", 0, 1, 0.01f));
		ratioParam.setDescription("ratio of discarded values below Imin and above Imax (in [0,1])");
		inputParams.add(ignoreNegParam = new ParamBoolean("set negative values to zero", true));
		inputParams.add(ignoreZeroParam = new ParamBoolean("ignore zero values", true));
		inputParams.add(scalingParam = new ParamFloat("Output scaling", 0, 1e10f, 1.0f));
		scalingParam.setDescription("scaling the output image into [0,S]");
		
		algorithm = new IntensityRangeNormalization();
		
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
		outputParams.add(resultImage = new ParamVolume("Normalized Image",VoxelType.FLOAT));
		outputParams.setName("normalize_image");
		outputParams.setLabel("Normalize Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(inImage);
		ImageHeader header = Interface.getHeader(inImage);
		int[] dims = Interface.getDimensions(inImage);
		float[] res = Interface.getResolutions(inImage);
		
		// main algorithm
		algorithm = new IntensityRangeNormalization();
		
		algorithm.setInputImage(Interface.getFloatImage3D(inImage));
		if (Interface.isValid(maskImage)) 
			algorithm.setMaskImage(Interface.getIntegerImage3D(maskImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		// parameters
		algorithm.setNormalization(normParam.getValue());
		algorithm.setRobustnessRatio(ratioParam.getValue().floatValue());
		algorithm.setNegativeValuesToZero(ignoreNegParam.getValue().booleanValue());
		algorithm.setIgnoreZeroValues(ignoreZeroParam.getValue().booleanValue());
		algorithm.setOutputScaling(scalingParam.getValue().floatValue());
		
		algorithm.execute();
		
		Interface.setFloatImage3D(algorithm.getNormalizedImage(), dims, resultImage, name+"_rnorm_img", header);
	}
	
}
