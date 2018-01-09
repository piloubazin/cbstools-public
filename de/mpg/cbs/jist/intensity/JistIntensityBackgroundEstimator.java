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
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.intensity.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityBackgroundEstimator extends ProcessingAlgorithm {

	private IntensityBackgroundEstimator algorithm;
	
	// jist containers
	private ParamVolume inputImage;
	private ParamOption distribParam;
	private ParamFloat ratioParam;
	
	private ParamVolume maskedImage;
	private ParamVolume probaImage;
	private ParamVolume maskImage;
	
	// parameters
	private		static final String[]	distribs = {"exponential","half-normal"};
	private		String		distrib = "exponential";
	
	private		static final byte	HNORM = 1;
	private		static final byte	EXP = 2;
	
	private ParamBoolean	skip0Param;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		
		inputParams.add(distribParam = new ParamOption("Background distribution", distribs));
		distribParam.setDescription("Model distribution for background noise (default is exponential, half-normal is less stringent).");
		distribParam.setValue(distrib);
		
		inputParams.add(ratioParam = new ParamFloat("Robust min, max thresholding", 0.0f, 1.0f, 0.0f));
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", true));
		
		algorithm = new IntensityBackgroundEstimator();
		
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
		outputParams.add(maskedImage = new ParamVolume("Masked Normalized Image",VoxelType.FLOAT));
		outputParams.add(probaImage = new ParamVolume("Foreground Proba Image",VoxelType.FLOAT));
		outputParams.add(maskImage = new ParamVolume("Foreground Mask Image",VoxelType.UBYTE));
		
		outputParams.setName("background images");
		outputParams.setLabel("background images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions(inputImage);
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		
		// main algorithm
		algorithm = new IntensityBackgroundEstimator();
		
		algorithm.setInputImage(Interface.getFloatImage3D(inputImage));
		algorithm.setDimensions(dims);
		algorithm.setRobustMinMaxThresholding(ratioParam.getValue().floatValue());
		algorithm.setBackgroundDistribution(distribParam.getValue());
		algorithm.setSkipZeroValues(skip0Param.getValue().booleanValue());
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getMaskedImage(), dims, maskedImage, name+"_masked", header);
		Interface.setFloatImage3D(algorithm.getProbaImage(), dims, probaImage, name+"_proba", header);
		Interface.setUByteImage3D(algorithm.getMask(), dims, maskImage, name+"_mask", header);
	}


}
