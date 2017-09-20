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

import de.mpg.cbs.core.intensity.IntensityPropagate;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityPropagate extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	normtypes = {"max","mean","min"};
	private		String		normtype = "max";
	private		static final String[]	targettypes = {"zero","mask","lower","higher"};
	private		String		targettype = "zero";
	
	// jist containers
	private ParamVolume inImage;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private	ParamOption	 normParam;
	private ParamFloat distParam;
	private	ParamOption	 targetParam;
	private ParamFloat scalingParam;
	//private ParamBoolean ignoreNegParam;
	//private ParamBoolean ignoreZeroParam;
	
	private IntensityPropagate algorithm;
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("Mask Image"));
		maskImage.setMandatory(false);
		inputParams.add(normParam = new ParamOption("combination method", normtypes));
		normParam.setDescription("use the mean, max or min data from neighboring voxels");
		normParam.setValue(normtype);
		inputParams.add(distParam = new ParamFloat("Propagation distance (mm)", 0, 50, 5));
		distParam.setDescription("distance for the propagation (note: this algorithm will be slow for large distances)");
		inputParams.add(targetParam = new ParamOption("target values", targettypes));
		targetParam.setDescription("propagate into zero, lower or higher neighboring voxels");
		targetParam.setValue(targettype);
		inputParams.add(scalingParam = new ParamFloat("Propagation scaling factor", 0, 1, 1.0f));
		scalingParam.setDescription("scaling factor for the propagated values");
		//inputParams.add(ignoreNegParam = new ParamBoolean("set negative values to zero", true));
		//inputParams.add(ignoreZeroParam = new ParamBoolean("ignore zero values", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("Intensity Propagation");
		inputParams.setName("IntensityPropagation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Propagates the values inside the mask (or non-zero) into the neighboring voxels");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Propagated Image",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.setName("propagate_image");
		outputParams.setLabel("Propagate Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		algorithm = new IntensityPropagate();
		
		String name = Interface.getName(inImage);
		ImageHeader header = Interface.getHeader(inImage);
		
		int[] dims = Interface.getDimensions(inImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(inImage));
		
		if (Interface.isImage4D(inImage)) algorithm.setInputImage(Interface.getFloatImage4D(inImage));
		else algorithm.setInputImage(Interface.getFloatImage3D(inImage));

		if (maskImage.getValue()!=null) algorithm.setMaskImage(Interface.getUByteImage3D(maskImage));
		
		algorithm.setTargetVoxels(targetParam.getValue());
		algorithm.setPropagationDistance(distParam.getValue().floatValue());
		algorithm.setCombinationMethod(normParam.getValue());
		algorithm.setPropogationScalingFactor(scalingParam.getValue().floatValue());
		
		algorithm.execute();
		
		// outputs
		if (Interface.isImage4D(inImage)) Interface.setFloatImage4D(algorithm.getResultImage(), dims, dims[3], resultImage, name+"_propag", header);
		else Interface.setFloatImage3D(algorithm.getResultImage(), dims, resultImage, name+"_propag", header);
		
	}
	
}
