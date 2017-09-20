package de.mpg.cbs.jist.brain;

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
import de.mpg.cbs.core.brain.BrainMgdmMultiSegmentation2;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMgdmMultiSegmentation2 extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;
	private ParamVolume input4Image;
	
	private ParamOption type1Param;
	private ParamOption type2Param;
	private ParamOption type3Param;
	private ParamOption type4Param;
	private static final String[] inputTypes = BrainMgdmMultiSegmentation2.inputTypes;
	
	private static final String inputTypeInfo = BrainMgdmMultiSegmentation2.inputTypeInfo;									
													
	private ParamFile atlasParam;

	private ParamInteger 	iterationParam;
	private ParamInteger 	stepParam;
	private ParamFloat 	changeParam;
	private	ParamFloat 	forceParam;
	private ParamFloat 	curvParam;
	private ParamFloat 	scaleParam;
	
	private ParamBoolean	computePosteriors;
	private ParamBoolean	adjustIntensPriors;
	private ParamBoolean	diffuseProbabilities;
	private ParamFloat		diffuseParam;
	
	private ParamOption 	topologyParam;
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamOption 	outputParam;
	private static final String[] outputTypes = {"label_memberships","raw_memberships","segmentation","debug"};
	//private static final String[] outputTypes = {"segmentation","memberships","cortex"};
	private ParamBoolean	normalizeQuantitative;
	
	private ParamVolume segmentImage;
	private ParamVolume mgdmImage;
	private ParamVolume idImage;
	
	private ParamVolume membershipImage;
	private ParamVolume labelImage;
	
	private BrainMgdmMultiSegmentation2 algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Data");
		
		
		
		imageParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		imageParams.add(type1Param = new ParamOption("Contrast Type 1",inputTypes));
		type1Param.setValue("none");
		type1Param.setDescription(inputTypeInfo);
		
		imageParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		imageParams.add(type2Param = new ParamOption("Contrast Type 2",inputTypes));
		input2Image.setMandatory(false);
		type2Param.setValue("none");
		type2Param.setDescription(inputTypeInfo);
		
		imageParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		imageParams.add(type3Param = new ParamOption("Contrast Type 3",inputTypes));
		input3Image.setMandatory(false);
		type3Param.setValue("none");
		type3Param.setDescription(inputTypeInfo);
		
		imageParams.add(input4Image = new ParamVolume("Contrast Image 4 (opt)"));
		imageParams.add(type4Param = new ParamOption("Contrast Type 4",inputTypes));
		input4Image.setMandatory(false);
		type4Param.setValue("none");
		type4Param.setDescription(inputTypeInfo);
		
		imageParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		
		mainParams.add(forceParam = new ParamFloat("Data weight", -1E10f, 1E10f, 0.1f));
		mainParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.4f));
		mainParams.add(scaleParam = new ParamFloat("Posterior scale (mm)", -1E10f, 1E10f, 5.0f));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(changeParam = new ParamFloat("Min change", 0f, 1f, 0.001f));
		mainParams.add(stepParam = new ParamInteger("Steps", 0, 100000, 5));
		
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(computePosteriors = new ParamBoolean("Compute posteriors", true));
		mainParams.add(adjustIntensPriors = new ParamBoolean("Adjust intensity priors", true));
		mainParams.add(diffuseProbabilities = new ParamBoolean("Diffuse probabilities", true));
		mainParams.add(diffuseParam = new ParamFloat("Diffusion factor", 0f, 10f, 0.5f));
		
		mainParams.add(outputParam = new ParamOption("Output images", outputTypes));
		outputParam.setValue("label_memberships");

		mainParams.add(normalizeQuantitative = new ParamBoolean("Normalize quantitative maps", true));
		
		inputParams.add(mainParams);
		
		algorithm = new BrainMgdmMultiSegmentation2();
		
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
		outputParams.add(segmentImage = new ParamVolume("Segmented Brain Image",VoxelType.BYTE));
		outputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image",VoxelType.BYTE));
		outputParams.add(idImage = new ParamVolume("Segmentation Ids Image",VoxelType.BYTE));
		
		outputParams.add(membershipImage = new ParamVolume("Posterior Maximum Memberships (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(labelImage = new ParamVolume("Posterior Maximum Labels (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		membershipImage.setMandatory(false);
		labelImage.setMandatory(false);

		outputParams.setName("brain segmentation images");
		outputParams.setLabel("brain segmentation images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new BrainMgdmMultiSegmentation2();
		
		// i/o variables
		String name = Interface.getName(input1Image);
		ImageHeader header = Interface.getHeader(input1Image);
		
		int[] dims = Interface.getDimensions(input1Image);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(input1Image));
		algorithm.setOrientations(Interface.getOrientations(input1Image));
	
		// pass the input parameters
		algorithm.setContrastImage1(Interface.getFloatImage3D(input1Image));
		if (Interface.isValid(input2Image)) algorithm.setContrastImage2(Interface.getFloatImage3D(input2Image));
		if (Interface.isValid(input3Image)) algorithm.setContrastImage3(Interface.getFloatImage3D(input3Image));
		if (Interface.isValid(input4Image)) algorithm.setContrastImage4(Interface.getFloatImage3D(input4Image));
		
		algorithm.setContrastType1(type1Param.getValue());
		algorithm.setContrastType2(type2Param.getValue());
		algorithm.setContrastType3(type3Param.getValue());
		algorithm.setContrastType4(type4Param.getValue());
			
		algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
		
		algorithm.setDataWeight(forceParam.getValue().floatValue());
		algorithm.setCurvatureWeight(curvParam.getValue().floatValue());
		algorithm.setPosteriorScale_mm(scaleParam.getValue().floatValue());
		
		algorithm.setMaxIterations(iterationParam.getValue().intValue());
		algorithm.setMinChange(changeParam.getValue().floatValue());
		algorithm.setSteps(stepParam.getValue().intValue());
		
		algorithm.setTopology(topologyParam.getValue());
		
		algorithm.setComputePosterior(computePosteriors.getValue().booleanValue());
		algorithm.setAdjustIntensityPriors(adjustIntensPriors.getValue().booleanValue());
		algorithm.setDiffuseProbabilities(diffuseProbabilities.getValue().booleanValue());
		algorithm.setDiffusionFactor(diffuseParam.getValue().floatValue());
		
		algorithm.setOutputImages(outputParam.getValue());
		
		algorithm.setNormalizeQuantitativeMaps(normalizeQuantitative.getValue().booleanValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setIntegerImage3D(algorithm.getSegmentedBrainImage(), dims, segmentImage, name+"_seg", header);
		Interface.setFloatImage3D(algorithm.getLevelsetBoundaryImage(), dims, mgdmImage, name+"_mgdm", header);
		Interface.setUByteImage3D(algorithm.getSegmentedIdsImage(), dims, idImage, name+"_segid", header);
		
		Interface.setUByteImage4D(algorithm.getPosteriorMaximumLabels4D(), dims, algorithm.getOutput4Dlength(), labelImage, name+"_lbls", header);
		Interface.setFloatImage4D(algorithm.getPosteriorMaximumMemberships4D(), dims, algorithm.getOutput4Dlength(), membershipImage, name+"_mems", header);
	
		return;
	}
}
