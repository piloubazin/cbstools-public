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
	private ParamVolume priorImage;
	private ParamFloat ratioParam;
	private ParamFloat bgscaleParam;
	private ParamFloat bgprobaParam;
	private ParamBoolean bgincludedParam;
	private ParamOption mergeParam;
	
	private ParamVolume 	probaImage;
	private ParamVolume 	bgmaskImage;
	private ParamVolume 	mgdmImage;
		
	private SegmentationDistanceBasedProbability algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation (opt if priors)"));
		segImage.setMandatory(false);
		inputParams.add(priorImage = new ParamVolume("Prior probabilities (4D, opt if segmentation)"));
		priorImage.setMandatory(false);
		inputParams.add(ratioParam = new ParamFloat("Distance ratio",0.0f,1.0f,0.5f));
		inputParams.add(bgscaleParam = new ParamFloat("Background Distance (mm)",0.0f,1.0e15f,3.0f));
		inputParams.add(bgprobaParam = new ParamFloat("Background Probability",0.0f,1.0f,0.5f));
		inputParams.add(bgincludedParam = new ParamBoolean("Background included",true));
		inputParams.add(mergeParam = new ParamOption("Probability merging",algorithm.mergingTypes));
		mergeParam.setValue(algorithm.mergingType);
		
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
		outputParams.add(bgmaskImage = new ParamVolume("Background mask",VoxelType.INT));
		outputParams.add(mgdmImage = new ParamVolume("Distance-based MGDM Function",VoxelType.FLOAT));
		
		outputParams.setName(algorithm.getName());
		outputParams.setLabel(algorithm.getLabel());
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new SegmentationDistanceBasedProbability();
		
		// i/o variables
		String name;
		ImageHeader header;
		int[] dims;
		if (segImage.getValue()!=null) {
			name = Interface.getName(segImage);
			header = Interface.getHeader(segImage);
		
			dims = Interface.getDimensions(segImage);
			algorithm.setDimensions(dims);
			algorithm.setResolutions(Interface.getResolutions(segImage));
		} else {
			name = Interface.getName(priorImage);
			header = Interface.getHeader(priorImage);
		
			dims = Interface.getDimensions(priorImage);
			algorithm.setDimensions(dims);
			algorithm.setResolutions(Interface.getResolutions(priorImage));
		}
		// pass the input parameters
		if (segImage.getValue()!=null) algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		if (priorImage.getValue()!=null) {
			if (Interface.isImage4D(priorImage)) algorithm.setPriorProbabilityImage(Interface.getFloatImage4D(priorImage));
			else  algorithm.setPriorProbabilityImage(Interface.getFloatImage3D(priorImage));
		}
		algorithm.setDistanceRatio(ratioParam.getValue().floatValue());
		algorithm.setBackgroundDistance_mm(bgscaleParam.getValue().floatValue());
		algorithm.setBackgroundProbability(bgprobaParam.getValue().floatValue());
		algorithm.setBackgroundIncluded(bgincludedParam.getValue().booleanValue());
		algorithm.setProbabilityMerging(mergeParam.getValue());
				
		algorithm.execute();
		
		// outputs
		Interface.setFloatImage4D(algorithm.getProbabilityImage(), dims, algorithm.getLabelNumber(), probaImage, name+"_prob", header);
		Interface.setIntegerImage3D(algorithm.getBackgroundMaskImage(), dims, bgmaskImage, name+"_mask", header);
		Interface.setFloatImage3D(algorithm.getMgdmImage(), dims, mgdmImage, name+"_mgdm", header);
		
		return;
	}
}
