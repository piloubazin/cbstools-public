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
import de.mpg.cbs.core.segmentation.SegmentationCellMgdm;

/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationCellMgdm extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;

	private ParamOption type1Param;
	private ParamOption type2Param;
	private ParamOption type3Param;
	private static final String[] inputTypes = SegmentationCellMgdm.inputTypes;
														
	private ParamOption 	dimParam;
	private static final String[] dimTypes = SegmentationCellMgdm.dimTypes;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 	changeParam;
	private	ParamFloat 	forceParam;
	private ParamFloat 	curvParam;
	private ParamFloat 	cellthresholdParam;
	
	private ParamOption 	topologyParam;
	private static final String[] topoTypes = SegmentationCellMgdm.topoTypes;
	
	private ParamVolume segmentImage;
	private ParamVolume mgdmImage;

	private SegmentationCellMgdm algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Data");
		
		inputParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		inputParams.add(type1Param = new ParamOption("Contrast Type 1",inputTypes));
		
		inputParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		inputParams.add(type2Param = new ParamOption("Contrast Type 2",inputTypes));
		input2Image.setMandatory(false);
		
		inputParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		inputParams.add(type3Param = new ParamOption("Contrast Type 3",inputTypes));
		input3Image.setMandatory(false);
		
		inputParams.add(dimParam = new ParamOption("Data dimensions", dimTypes));
		
		inputParams.add(forceParam = new ParamFloat("Data weight", -1E10f, 1E10f, 0.1f));
		inputParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.4f));
		inputParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		inputParams.add(changeParam = new ParamFloat("Min change", 0f, 1f, 0.001f));
		inputParams.add(cellthresholdParam = new ParamFloat("Cell threshold", 0f, 1f, 0.8f));
		
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		algorithm = new SegmentationCellMgdm();
		
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

		outputParams.setName("segmentation images");
		outputParams.setLabel("segmentation images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new SegmentationCellMgdm();
		
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
		
		algorithm.setContrastType1(type1Param.getValue());
		algorithm.setContrastType2(type2Param.getValue());
		algorithm.setContrastType3(type3Param.getValue());
		
		algorithm.setDataStackDimension(dimParam.getValue());
			
		algorithm.setDataWeight(forceParam.getValue().floatValue());
		algorithm.setCurvatureWeight(curvParam.getValue().floatValue());
		algorithm.setMaxIterations(iterationParam.getValue().intValue());
		algorithm.setMinChange(changeParam.getValue().floatValue());
		algorithm.setCellThreshold(cellthresholdParam.getValue().floatValue());
		
		algorithm.setTopology(topologyParam.getValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setIntegerImage3D(algorithm.getSegmentedImage(), dims, segmentImage, name+"_seg", header);
		Interface.setFloatImage3D(algorithm.getLevelsetBoundaryImage(), dims, mgdmImage, name+"_mgdm", header);
		
		return;
	}
}
