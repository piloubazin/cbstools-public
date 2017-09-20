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
import de.mpg.cbs.core.brain.BrainEnhanceRegionContrast;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainEnhanceRegionContrast extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume intensImage;
	private ParamVolume segImage;
	private ParamVolume mgdmImage;
	
	private ParamFile atlasParam;
	private ParamOption 	regionParam;
	private ParamOption 	backgroundParam;
	private ParamFloat		distanceParam;
	
	private ParamVolume regionImage;
	private ParamVolume backgroundImage;
	private ParamVolume probaRegionImage;
	private ParamVolume probaBackgroundImage;
	private ParamVolume pvRegionImage;
	private ParamVolume pvBackgroundImage;
		
	private BrainEnhanceRegionContrast algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensImage = new ParamVolume("Intensity Image"));
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		
		inputParams.add(regionParam = new ParamOption("Enhanced Region", BrainEnhanceRegionContrast.regionTypes));
		inputParams.add(backgroundParam = new ParamOption("Contrast Background", BrainEnhanceRegionContrast.backgroundTypes));
		inputParams.add(distanceParam = new ParamFloat("Partial voluming distance (voxels)", 0.0f, 10.0f, 1.0f));
		
		algorithm = new BrainEnhanceRegionContrast();
		
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
		outputParams.add(regionImage = new ParamVolume("Region Mask",VoxelType.BYTE));
		outputParams.add(backgroundImage = new ParamVolume("Background Mask",VoxelType.BYTE));
		outputParams.add(probaRegionImage = new ParamVolume("Region Probability",VoxelType.FLOAT));
		outputParams.add(probaBackgroundImage = new ParamVolume("Background Probability",VoxelType.FLOAT));
		outputParams.add(pvRegionImage = new ParamVolume("Region partial volume",VoxelType.FLOAT));
		outputParams.add(pvBackgroundImage = new ParamVolume("Background partial volume",VoxelType.FLOAT));
		
		outputParams.setName("enhanced images");
		outputParams.setLabel("enhanced images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		algorithm = new BrainEnhanceRegionContrast(); 
		
		// i/o variables
		String intensName = Interface.getName(intensImage);
		String segName = Interface.getName(segImage);
		String mgdmName = Interface.getName(mgdmImage);
		
		ImageHeader header = Interface.getHeader(intensImage);
		int[] dims = Interface.getDimensions(segImage);
		int comps = Interface.getComponents(intensImage);
		
		algorithm.setDimensions(dims);
		algorithm.setComponents(comps);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setLevelsetBoundaryImage(Interface.getFloatImage3D(mgdmImage));
		if (comps>1) algorithm.setIntensityImage(Interface.getFloatImage4D(intensImage));
		else algorithm.setIntensityImage(Interface.getFloatImage3D(intensImage));
		
		algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
		algorithm.setEnhancedRegion(regionParam.getValue());
		algorithm.setContrastBackground(backgroundParam.getValue());
		algorithm.setPartialVolumingDistance(distanceParam.getValue().floatValue());
		
		algorithm.execute();
		
		Interface.setUByteImage3D(algorithm.getRegionMask(), dims, regionImage, segName+algorithm.getRegionName(), header);
		Interface.setUByteImage3D(algorithm.getBackgroundMask(), dims, backgroundImage, segName+algorithm.getBackgroundName(), header);
		
		Interface.setFloatImage3D(algorithm.getRegionProbability(), dims, probaRegionImage, intensName+algorithm.getRegionName(), header);
		Interface.setFloatImage3D(algorithm.getBackgroundProbability(), dims, probaBackgroundImage, intensName+algorithm.getBackgroundName(), header);
		
		Interface.setFloatImage3D(algorithm.getRegionPartialVolume(), dims, pvRegionImage, mgdmName+algorithm.getRegionName(), header);
		Interface.setFloatImage3D(algorithm.getBackgroundPartialVolume(), dims, pvBackgroundImage, mgdmName+algorithm.getBackgroundName(), header);
		
		return;
	}
	
}
