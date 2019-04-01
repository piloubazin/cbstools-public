package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.shape.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistShapeTopologyCorrection extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	private ParamVolume startImage;
	
	private ParamOption inputParam;
	private ParamOption connectParam;
	private ParamOption propagParam;
	
    private ParamFloat highestParam;
	private ParamFloat lowestParam;
	private ParamBoolean normalizeParam;
	private ParamFloat mindistParam;
		
	private ParamVolume correctImage;
	private ParamVolume correctobjImage;
	
	private ShapeTopologyCorrection2 algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Shape image"));
		inputParams.add(startImage = new ParamVolume("Starting object image (opt)"));
		startImage.setMandatory(false);
		
		inputParams.add(inputParam = new ParamOption("Shape image type", ShapeTopologyCorrection2.inputTypes));
		inputParams.add(propagParam = new ParamOption("Correction direction", ShapeTopologyCorrection2.propagTypes));
		inputParams.add(connectParam = new ParamOption("Topology", ShapeTopologyCorrection2.connectTypes));
		
		//inputParams.add(highestParam = new ParamFloat("Highest (absolute or relative) value", -1e16f, 1e16f, 1.0f));
		//inputParams.add(lowestParam = new ParamFloat("Lowest (absolute or relative) value", -1e16f, 1e16f, 0.0f));
		//inputParams.add(normalizeParam = new ParamBoolean("Normalize intensity range", true));
		inputParams.add(mindistParam = new ParamFloat("Minimum distance", 0.0f, 10.0f, 0.0001f));
		
		algorithm = new ShapeTopologyCorrection2();
		
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
		outputParams.add(correctImage = new ParamVolume("Corrected Image",null));
		outputParams.add(correctobjImage = new ParamVolume("Corrected Object",null));
		
		outputParams.setName("skeleton images");
		outputParams.setLabel("skeleton images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new ShapeTopologyCorrection2();
		
		// i/o variables
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		
		int[] dims = Interface.getDimensions(inputImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(inputImage));
		
		algorithm.setShapeImage(Interface.getFloatImage3D(inputImage));
		algorithm.setShapeImageType(inputParam.getValue());
		
		if (Interface.isValid(startImage))
            algorithm.setStartingObjectImage(Interface.getIntegerImage3D(startImage));
		
        algorithm.setPropagationDirection(propagParam.getValue());
        algorithm.setTopology(connectParam.getValue());
        
		//algorithm.setLowerThreshold(lowestParam.getValue().floatValue());
		//algorithm.setHigherThreshold(highestParam.getValue().floatValue());
		//algorithm.setNormalizeIntensityRange(normalizeParam.getValue().booleanValue());
		algorithm.setMinimumDistance(mindistParam.getValue().floatValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setFloatImage3D(algorithm.getCorrectedImage(), dims, correctImage, name+"_topo_corr", header);
		Interface.setIntegerImage3D(algorithm.getCorrectedObjectImage(), dims, correctobjImage, name+"_topo_obj", header);
	}
}
