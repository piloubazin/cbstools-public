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
import de.mpg.cbs.core.segmentation.SegmentationLesionExtraction;

/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationLesionExtraction extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume probaImage;
	private ParamVolume segImage;
	private ParamVolume mgdmImage;
	private ParamVolume priorImage;
	
	private ParamFile atlasParam;
	private ParamFloat		gmdistanceParam;
	private ParamFloat		csfdistanceParam;
	private ParamFloat		lesdistanceParam;
	private ParamFloat		minprobaThresholdParam;
	private ParamFloat		maxprobaThresholdParam;
	private ParamFloat		minSizeParam;
	
	private ParamVolume scoreImage;
	private ParamVolume pvRegionImage;
	private ParamVolume lesionsizeImage;
	private ParamVolume lesionprobaImage;
	private ParamVolume boundarypvImage;
	private ParamVolume labelImage;
		
	private SegmentationLesionExtraction algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(probaImage = new ParamVolume("Probability Image"));
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		inputParams.add(priorImage = new ParamVolume("Location Prior Image (opt)"));
		priorImage.setMandatory(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		inputParams.add(gmdistanceParam = new ParamFloat("GM boundary partial voluming distance (voxels)",0.0f,100.0f,1.0f));
		inputParams.add(csfdistanceParam = new ParamFloat("CSF boundary partial voluming distance (voxels)",0.0f,100.0f,2.0f));
		inputParams.add(lesdistanceParam = new ParamFloat("Lesion clustering distance (voxels)",0.0f,100.0f,2.0f));
		inputParams.add(minprobaThresholdParam = new ParamFloat("Probability minimum threshold (main threshold)",0.0f,1.0f,0.84f));
		inputParams.add(maxprobaThresholdParam = new ParamFloat("Probability maximum threshold (to exclude dirty WM)",0.0f,1.0f,0.97f));
		inputParams.add(minSizeParam = new ParamFloat("Small lesion size (voxels)",0.0f,10000.0f,4.0f));

		algorithm = new SegmentationLesionExtraction();
		
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
		outputParams.add(scoreImage = new ParamVolume("Lesions score",VoxelType.FLOAT));
		outputParams.add(pvRegionImage = new ParamVolume("Region partial volume",VoxelType.FLOAT));
		outputParams.add(lesionsizeImage = new ParamVolume("Lesions size",VoxelType.FLOAT));
		outputParams.add(lesionprobaImage = new ParamVolume("Lesions probabilities",VoxelType.FLOAT));
		outputParams.add(boundarypvImage = new ParamVolume("Boundari partial volume",VoxelType.FLOAT));
		outputParams.add(labelImage = new ParamVolume("Lesion labels",VoxelType.INT));
		
		outputParams.setName("lesion images");
		outputParams.setLabel("lesion images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		algorithm = new SegmentationLesionExtraction(); 
		
		// i/o variables
		String probaName = Interface.getName(probaImage);
		String segName = Interface.getName(segImage);
		String mgdmName = Interface.getName(mgdmImage);
		
		ImageHeader header = Interface.getHeader(probaImage);
		int[] dims = Interface.getDimensions(segImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setLevelsetBoundaryImage(Interface.getFloatImage3D(mgdmImage));
		algorithm.setProbaImage(Interface.getFloatImage3D(probaImage));
		if (Interface.isValid(priorImage)) algorithm.setLocationPriorImage(Interface.getFloatImage3D(priorImage));
		
		algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
		
		algorithm.setGMPartialVolumingDistance(gmdistanceParam.getValue().floatValue());
		algorithm.setCSFPartialVolumingDistance(csfdistanceParam.getValue().floatValue());
		algorithm.setLesionClusteringDistance(lesdistanceParam.getValue().floatValue());
		algorithm.setMinProbabilityThreshold(minprobaThresholdParam.getValue().floatValue());
		algorithm.setMaxProbabilityThreshold(maxprobaThresholdParam.getValue().floatValue());
		algorithm.setMinimumSize(minSizeParam.getValue().floatValue());
		
		algorithm.execute();
		
		Interface.setFloatImage3D(algorithm.getLesionScore(), dims, scoreImage, probaName+"_lesion_score", header);
		Interface.setFloatImage3D(algorithm.getRegionPrior(), dims, pvRegionImage, probaName+"_lesion_prior", header);
		Interface.setFloatImage3D(algorithm.getLesionSize(), dims, lesionsizeImage, probaName+"_lesion_size", header);
		Interface.setFloatImage3D(algorithm.getLesionProba(), dims, lesionprobaImage, probaName+"_lesion_proba", header);
		Interface.setFloatImage3D(algorithm.getBoundaryPartialVolume(), dims, boundarypvImage, probaName+"_lesion_pv", header);
		Interface.setIntegerImage3D(algorithm.getLesionLabels(), dims, labelImage, probaName+"_lesion_labels", header);
		
		return;
	}
	
}
