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
import de.mpg.cbs.core.segmentation.SegmentationSureSeg;

/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationSureSeg extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;
	
	private ParamVolume var1Image;
	private ParamVolume var2Image;
	private ParamVolume var3Image;
	
	private ParamFloat scale1Param;
	private ParamFloat scale2Param;
	private ParamFloat scale3Param;
	
	private ParamVolume 	probaImage;
	private ParamVolume 	labelImage;
	private ParamVolume 	maskImage;
	
	private ParamInteger 	nbestParam;
	private ParamBoolean	includeBgParam;
	private ParamBoolean	rescaleProbaParam;
	private ParamBoolean	computeDistributionParam;
	private ParamBoolean	computeNoiseParam;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 		imgscaleParam;
	private	ParamFloat 		certainscaleParam;
	private ParamFloat 		mincertaintyParam;
	private ParamFloat	 	scaleParam;
	private ParamInteger 	neighborParam;
	private ParamFloat	 	diffratioParam;
	
	private ParamVolume segmentImage;
	private ParamVolume maxprobaImage;
	private ParamVolume maxidImage;
	private ParamVolume debugImage;
	
	private SegmentationSureSeg algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Data");
		
		imageParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		imageParams.add(var1Image = new ParamVolume("Noise Image 1 (opt)"));
		var1Image.setMandatory(false);
		imageParams.add(scale1Param = new ParamFloat("Contrast Scale 1",0.0f,1e15f,1.0f));
		
		imageParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		imageParams.add(var2Image = new ParamVolume("Noise Image 2 (opt)"));
		var2Image.setMandatory(false);
		imageParams.add(scale2Param = new ParamFloat("Contrast Scale 2",0.0f,1e15f,1.0f));
		input2Image.setMandatory(false);
		
		imageParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		imageParams.add(var3Image = new ParamVolume("Noise Image 3 (opt)"));
		var3Image.setMandatory(false);
		imageParams.add(scale3Param = new ParamFloat("Contrast Scale 3",0.0f,1e15f,1.0f));
		input3Image.setMandatory(false);
		
		imageParams.add(probaImage = new ParamVolume("Label Probabilities (3D or 4D)"));
		imageParams.add(labelImage = new ParamVolume("Label Segmentation (only if 3D proba)"));
		labelImage.setMandatory(false);
		
		imageParams.add(maskImage = new ParamVolume("Image Mask (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		
		mainParams.add(nbestParam = new ParamInteger("Labeling depth",0,20,4));
		mainParams.add(includeBgParam = new ParamBoolean("Background included",true));
		mainParams.add(rescaleProbaParam = new ParamBoolean("Rescale individual probabilities",true));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(imgscaleParam = new ParamFloat("Image Scale", 0.0f, 100.0f, 0.5f));
		mainParams.add(certainscaleParam = new ParamFloat("Certainty Scale", 0.0f, 100.0f, 2.0f));
		mainParams.add(mincertaintyParam = new ParamFloat("Min Certainty", 0, 1, 0.5f));
		mainParams.add(neighborParam = new ParamInteger("Neighborhood size", 0, 26, 6));
		mainParams.add(diffratioParam = new ParamFloat("Min difference ratio", 0, 1, 0.01f));
		
		mainParams.add(computeDistributionParam = new ParamBoolean("Re-estimate intensity distributions",true));
		mainParams.add(computeNoiseParam = new ParamBoolean("Estimate image noise",true));
		
		inputParams.add(mainParams);
		
		algorithm = new SegmentationSureSeg();
		
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
		outputParams.add(segmentImage = new ParamVolume("Segmented Image",VoxelType.INT));
		outputParams.add(maxprobaImage = new ParamVolume("Maximum Probability Image (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(maxidImage = new ParamVolume("Segmentation Ids Image (4D)",VoxelType.UBYTE,-1,-1,-1,-1));
		outputParams.add(debugImage = new ParamVolume("Debug Image (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("sure segmentation images");
		outputParams.setLabel("sure segmentation images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new SegmentationSureSeg();
		
		// i/o variables
		String name = Interface.getName(input1Image);
		ImageHeader header = Interface.getHeader(input1Image);
		
		int[] dims = Interface.getDimensions(input1Image);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(input1Image));
		
		// pass the input parameters
		algorithm.setContrastImage1(Interface.getFloatImage3D(input1Image));
		if (Interface.isValid(input2Image)) algorithm.setContrastImage2(Interface.getFloatImage3D(input2Image));
		if (Interface.isValid(input3Image)) algorithm.setContrastImage3(Interface.getFloatImage3D(input3Image));
		
		algorithm.setContrastScale1(scale1Param.getValue().floatValue());
		algorithm.setContrastScale2(scale2Param.getValue().floatValue());
		algorithm.setContrastScale3(scale3Param.getValue().floatValue());
		
		if (Interface.isValid(var1Image)) algorithm.setNoiseImage1(Interface.getFloatImage3D(var1Image));
		if (Interface.isValid(var2Image)) algorithm.setNoiseImage2(Interface.getFloatImage3D(var2Image));
		if (Interface.isValid(var3Image)) algorithm.setNoiseImage3(Interface.getFloatImage3D(var3Image));
		
		if (Interface.isValid(labelImage)) {
			algorithm.setLabelSegmentation(Interface.getIntegerImage3D(labelImage));
			if (Interface.isImage4D(probaImage)) algorithm.setLabel4DProbabilities(Interface.getFloatImage4D(probaImage));
			else algorithm.setLabelMaxProbability(Interface.getFloatImage3D(probaImage));
		} else {
			algorithm.setLabel4DProbabilities(Interface.getFloatImage4D(probaImage));
		}
		if (Interface.isValid(maskImage)) algorithm.setImageMask(Interface.getUByteImage3D(maskImage));
		
		algorithm.setLabelDepth(nbestParam.getValue().byteValue());
		algorithm.setBackgroundIncluded(includeBgParam.getValue().booleanValue());
		algorithm.setRescaleIndividualProbabilities(rescaleProbaParam.getValue().booleanValue());
		
		algorithm.setMaxIterations(iterationParam.getValue().intValue());
		algorithm.setImageScale(imgscaleParam.getValue().floatValue());
		algorithm.setCertaintyScale(certainscaleParam.getValue().floatValue());
		algorithm.setMinCertainty(mincertaintyParam.getValue().floatValue());
		algorithm.setNeighborhoodSize(neighborParam.getValue().intValue());
		algorithm.setMinDifferenceRatio(diffratioParam.getValue().floatValue());
		
		algorithm.setReestimateIntensityDistributions(computeDistributionParam.getValue().booleanValue());
		algorithm.setEstimateNoise(computeNoiseParam.getValue().booleanValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setIntegerImage3D(algorithm.getSegmentationImage(), dims, segmentImage, name+"_seg", header);
		Interface.setFloatImage4D(algorithm.getMaxProbabilityImage(), dims, nbestParam.getValue().byteValue(), maxprobaImage, name+"_mems", header);
		Interface.setUByteImage4D(algorithm.getSegmentedIdsImage(), dims, nbestParam.getValue().byteValue(), maxidImage, name+"_lbls", header);
		
		Interface.setFloatImage4D(algorithm.getDebugImage(), dims, 2*neighborParam.getValue().byteValue()+1, debugImage, name+"_debug", header);
		
		return;
	}
}
