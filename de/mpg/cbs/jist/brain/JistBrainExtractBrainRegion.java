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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.core.brain.BrainExtractBrainRegion;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainExtractBrainRegion extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume segImage;
	private ParamVolume mgdmImage;
	private ParamVolume labelImage;
	private ParamVolume functionImage;
	private ParamFile atlasParam;
	private ParamOption 	regionParam;
	private static final String[] regionTypes = {"left_cerebrum", "right_cerebrum", "cerebrum", "cerebellum", 
													"cerebellum_brainstem", "subcortex", "tissues(anat)", "tissues(func)",
													"brain_mask"};
	private ParamBoolean	normalizeParam;
	private ParamBoolean	densityParam;
	private ParamFloat		scalingParam;
		
	private ParamVolume segStructureImage;
	private ParamVolume segInsideImage;
	private ParamVolume segBackgroundImage;
	private ParamVolume lvlStructureImage;
	private ParamVolume lvlInsideImage;
	private ParamVolume lvlBackgroundImage;
	private ParamVolume probaStructureImage;
	private ParamVolume probaInsideImage;
	private ParamVolume probaBackgroundImage;
		
	private BrainExtractBrainRegion algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		inputParams.add(functionImage = new ParamVolume("Maximum Membership Image (4D)"));
		inputParams.add(labelImage = new ParamVolume("Maximum Label Image (4D)"));
		
		segImage.setLoadAndSaveOnValidate(false);
		mgdmImage.setLoadAndSaveOnValidate(false);
		functionImage.setLoadAndSaveOnValidate(false);
		labelImage.setLoadAndSaveOnValidate(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		/*
		String filename = "/scr/watt1/bazin/Datasets/shape-atlases/whole-brain-8mm/separate-brain-atlas-nov11.txt";
        try {
			//ClassLoader cl = Thread.currentThread().getContextClassLoader();
			//atlasParam.setValue(cl.getResource(filename).getFile());
			atlasParam.setValue(filename);
        } catch (Exception e) {
        	System.out.print("Error: Unable to set default atlas\n");
        }
        */
		 
        inputParams.add(regionParam = new ParamOption("Extracted Region", regionTypes));
		inputParams.add(normalizeParam = new ParamBoolean("Normalize probabilities", true));
		inputParams.add(densityParam = new ParamBoolean("Estimate tissue densities", false));
		inputParams.add(scalingParam = new ParamFloat("Partial voluming distance (voxels)", 0.0f, 10.0f, 1.0f));
		
		algorithm = new BrainExtractBrainRegion();
		
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
		outputParams.add(segInsideImage = new ParamVolume("Inside (WM) Mask",VoxelType.BYTE));
		outputParams.add(segStructureImage = new ParamVolume("Structure (GM) Mask",VoxelType.BYTE));
		outputParams.add(segBackgroundImage = new ParamVolume("Background (CSF) Mask",VoxelType.BYTE));
		outputParams.add(lvlInsideImage = new ParamVolume("Inside (WM) Leveset",VoxelType.FLOAT));
		outputParams.add(lvlStructureImage = new ParamVolume("Structure (GM) Levelset",VoxelType.FLOAT));
		outputParams.add(lvlBackgroundImage = new ParamVolume("Background (CSF) Levelset",VoxelType.FLOAT));
		outputParams.add(probaInsideImage = new ParamVolume("Inside (WM) Probability",VoxelType.FLOAT));
		outputParams.add(probaStructureImage = new ParamVolume("Structure (GM) Probability",VoxelType.FLOAT));
		outputParams.add(probaBackgroundImage = new ParamVolume("Background (CSF) Probability",VoxelType.FLOAT));
		
		/*
		segInsideImage.setLoadAndSaveOnValidate(false);
		segStructureImage.setLoadAndSaveOnValidate(false);
		segBackgroundImage.setLoadAndSaveOnValidate(false);
		lvlInsideImage.setLoadAndSaveOnValidate(false);
		lvlStructureImage.setLoadAndSaveOnValidate(false);
		lvlBackgroundImage.setLoadAndSaveOnValidate(false);
		probaInsideImage.setLoadAndSaveOnValidate(false);
		probaStructureImage.setLoadAndSaveOnValidate(false);
		probaBackgroundImage.setLoadAndSaveOnValidate(false);
		*/
		
		outputParams.setName("extract atlas images");
		outputParams.setLabel("extract atlas images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		algorithm = new BrainExtractBrainRegion(); 
		
		// i/o variables
		String functName = Interface.getName(functionImage);
		String segName = Interface.getName(segImage);
		String mgdmName = Interface.getName(mgdmImage);
		
		ImageHeader header = Interface.getHeader(functionImage);
		int[] dims = Interface.getDimensions(segImage);
		int comps = Interface.getComponents(functionImage);
		
		algorithm.setDimensions(dims);
		algorithm.setComponents(comps);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setLevelsetBoundaryImage(Interface.getFloatImage3D(mgdmImage));
		algorithm.setMaximumMembershipImage(Interface.getFloatImage4D(functionImage));
		algorithm.setMaximumLabelImage(Interface.getIntegerImage4D(labelImage));
		
		algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
		algorithm.setExtractedRegion(regionParam.getValue());
		algorithm.setNormalizeProbabilities(normalizeParam.getValue().booleanValue());
		algorithm.setEstimateTissueDensities(densityParam.getValue().booleanValue());
		algorithm.setPartialVolumingDistance(scalingParam.getValue().floatValue());
		
		algorithm.execute();
		
		Interface.setUByteImage3D(algorithm.getInsideWMmask(), dims, segInsideImage, segName+"_"+algorithm.getInsideName(), header);
		Interface.setUByteImage3D(algorithm.getStructureGMmask(), dims, segStructureImage, segName+"_"+algorithm.getStructureName(), header);
		Interface.setUByteImage3D(algorithm.getBackgroundCSFmask(), dims, segBackgroundImage, segName+"_"+algorithm.getBackgroundName(), header);
		
		Interface.setFloatImage3D(algorithm.getInsideWMprobability(), dims, probaInsideImage, functName+"_"+algorithm.getInsideName(), header);
		Interface.setFloatImage3D(algorithm.getStructureGMprobability(), dims, probaStructureImage, functName+"_"+algorithm.getStructureName(), header);
		Interface.setFloatImage3D(algorithm.getBackgroundCSFprobability(), dims, probaBackgroundImage, functName+"_"+algorithm.getBackgroundName(), header);
		
		Interface.setFloatImage3D(algorithm.getInsideWMlevelset(), dims, lvlInsideImage, mgdmName+"_"+algorithm.getInsideName(), header);
		Interface.setFloatImage3D(algorithm.getStructureGMlevelset(), dims, lvlStructureImage, mgdmName+"_"+algorithm.getStructureName(), header);
		Interface.setFloatImage3D(algorithm.getBackgroundCSFlevelset(), dims, lvlBackgroundImage, mgdmName+"_"+algorithm.getBackgroundName(), header);
		
		return;
	}
	
}
