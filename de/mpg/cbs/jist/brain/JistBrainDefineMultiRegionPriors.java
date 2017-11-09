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
import de.mpg.cbs.core.brain.BrainDefineMultiRegionPriors;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainDefineMultiRegionPriors extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume segImage;
	private ParamVolume mgdmImage;
	
	private ParamFile atlasParam;
	//private ParamOption 	regionParam;
	private ParamOption 	methodParam;
	private ParamFloat		distanceParam;
	
	private ParamVolume interImage;
	private ParamVolume hornsImage;
	private ParamVolume icapsImage;
		
	private BrainDefineMultiRegionPriors algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		
		//inputParams.add(regionParam = new ParamOption("Defined Region", BrainDefineMultiRegionPriors.regionTypes));
		//inputParams.add(methodParam = new ParamOption("Definition Method", BrainDefineMultiRegionPriors.methodTypes));
		inputParams.add(distanceParam = new ParamFloat("Distance Offset (voxels)", 0.0f, 10.0f, 1.0f));
		
		algorithm = new BrainDefineMultiRegionPriors();
		
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
		outputParams.add(interImage = new ParamVolume("Inter-ventricular PV",VoxelType.FLOAT));
		outputParams.add(hornsImage = new ParamVolume("Ventricular horns PV",VoxelType.FLOAT));
		outputParams.add(icapsImage = new ParamVolume("Internal capsule PV",VoxelType.FLOAT));
		
		outputParams.setName("defined images");
		outputParams.setLabel("defined images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		algorithm = new BrainDefineMultiRegionPriors(); 
		
		// i/o variables
		String segName = Interface.getName(segImage);
		String mgdmName = Interface.getName(mgdmImage);
		
		ImageHeader header = Interface.getHeader(segImage);
		int[] dims = Interface.getDimensions(segImage);
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setLevelsetBoundaryImage(Interface.getFloatImage3D(mgdmImage));
		
		algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
		//algorithm.setDefinedRegion(regionParam.getValue());
		//algorithm.setDefinitionMethod(methodParam.getValue());
		algorithm.setDistanceOffset(distanceParam.getValue().floatValue());
		
		algorithm.execute();
		
		Interface.setFloatImage3D(algorithm.getInterVentricularPV(), dims, interImage, mgdmName+"_iv", header);
		Interface.setFloatImage3D(algorithm.getVentricularHornsPV(), dims, hornsImage, mgdmName+"_vh", header);
		Interface.setFloatImage3D(algorithm.getInternalCapsulePV(), dims, icapsImage, mgdmName+"_ic", header);
		
		return;
	}
	
}
