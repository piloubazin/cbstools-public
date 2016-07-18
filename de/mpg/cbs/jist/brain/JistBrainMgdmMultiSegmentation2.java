package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
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
	private static final String[] inputTypes = {"-- 3T --", "MPRAGE3T", "T1MAP3T", "MP2RAGE3T", 
													"HCPT1w", "HCPT2w", "NormMPRAGE", "FLAIR3T",
													"DWIFA3T", "DWIMD3T",
													"-- 7T --", "T1MAP7T", "MP2RAGE7T", "T2SW7T", "QSM7T", 
													"-- 9.4T --", "T1MAP9T", "MP2RAGE9T",
													"-- MPM --", "mpmPD3T", "mpmMT3T", "mpmR13T", "mpmR2s3T",
													"-- misc --", "Filters", "WMLesions", "PVDURA", "Labeling","none"};
	
	private static final String inputTypeInfo = "Currently available contrasts:\n"
			+"T1MAP7T: a 7T quantitative T1 map, \n"
			+"MP2RAGE7T: a T1-weighted image from 7T MP2RAGE (UNI), \n"
			+"T2SW7T: a 7T T2*-weighted image (devel), \n"
			+"SQSM7T: a 7T quantitative susceptibility map (devel), \n"
			+"T1MAP9T: a 9.4T quantitative T1 map (devel), \n"
			+"MP2RAGE9T: a T1-weighted image from 9.4T MP2RAGE (UNI) (devel), \n"
			+"MPRAGE3T: a 3T T1-weighted MPRAGE image, \n"
			+"T1MAP3T: a 3T quantitative T1 map, \n"
			+"MP2RAGE3T: a T1-weighted image from MP2RAGE (UNI), \n"
			+"HCPT1w: a T1-weighted image using the HCP sequence, \n"
			+"HCPT2w: a T2-weighted image using the HCP sequence, \n"
			+"NormMPRAGE: a 3T MPRAGE normalised for B1 shading, \n"
			+"FLAIR3T: a 3T FLAIR image, \n"
			+"dwiFA: a DWI fractional anisotropy image, \n"
			+"dwiMD: a DWI mean diffusivity image, \n"
			+"MPM: maps from the multi parameter mapping (MPM) sequence, \n"
			+"Filters: a composite image of outputs from dura, pv and arteries pre-processing, \n"
			+"WMLesions: a (probabilistic or binary) mask of detected white matter lesions, \n"
			+"PVDURA: a composite image of outputs from dura and pv (obsolete).\n";												
													
	private ParamFile atlasParam;

	private ParamInteger 	iterationParam;
	private ParamInteger 	stepParam;
	private ParamDouble 	changeParam;
	private	ParamDouble 	forceParam;
	private ParamDouble 	curvParam;
	private ParamDouble 	scaleParam;
	
	private ParamBoolean	computePosteriors;
	private ParamBoolean	adjustIntensPriors;
	private ParamBoolean	diffuseProbabilities;
	private ParamDouble		diffuseParam;
	
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
		
		mainParams.add(forceParam = new ParamDouble("Data weight", -1E10, 1E10, 0.2));
		mainParams.add(curvParam = new ParamDouble("Curvature weight", -1E10, 1E10, 0.8));
		mainParams.add(scaleParam = new ParamDouble("Posterior scale (mm)", -1E10, 1E10, 5.0));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(changeParam = new ParamDouble("Min change", 0, 1, 0.001));
		mainParams.add(stepParam = new ParamInteger("Steps", 0, 100000, 5));
		
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(computePosteriors = new ParamBoolean("Compute posteriors", true));
		mainParams.add(adjustIntensPriors = new ParamBoolean("Adjust intensity priors", true));
		mainParams.add(diffuseProbabilities = new ParamBoolean("Diffuse probabilities", true));
		mainParams.add(diffuseParam = new ParamDouble("Diffusion factor", 0, 10, 0.5));
		
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
