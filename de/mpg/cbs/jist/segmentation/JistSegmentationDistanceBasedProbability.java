package de.mpg.cbs.jist.segmentation;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.segmentation.SegmentationDistanceBasedProbability;

/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationDistanceBasedProbability extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume segImage;
	private ParamFloat ratioParam;
	private ParamFloat bgscaleParam;
	private ParamBoolean bgDistParam;
	
	private ParamVolume 	probaImage;
		
	private SegmentationDistanceBasedProbability algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation"));
		inputParams.add(ratioParam = new ParamFloat("Distance ratio",0.0f,1.0f,0.5f));
		inputParams.add(bgscaleParam = new ParamFloat("Background Distance (mm)",0.0f,1.0e15f,3.0f));
		inputParams.add(bgDistParam = new ParamBoolean("Use fixed background distance",true));
		
		
		algorithm = new SegmentationDistanceBasedProbability();
		
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
		outputParams.add(probaImage = new ParamVolume("Distance-based Probability Image (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName(algorithm.getName());
		outputParams.setLabel(algorithm.getLabel());
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new SegmentationDistanceBasedProbability();
		
		// i/o variables
		String name = Interface.getName(segImage);
		ImageHeader header = Interface.getHeader(segImage);
		
		int[] dims = Interface.getDimensions(segImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		// pass the input parameters
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setDistanceRatio(ratioParam.getValue().floatValue());
		algorithm.setBackgroundDistance_mm(bgscaleParam.getValue().floatValue());
		algorithm.setUseBackgroundDistance(bgDistParam.getValue().booleanValue());
				
		algorithm.execute();
		
		// outputs
		Interface.setFloatImage4D(algorithm.getProbabilityImage(), dims, algorithm.getLabelNumber(), probaImage, name+"_prob", header);
		
		return;
	}
}
