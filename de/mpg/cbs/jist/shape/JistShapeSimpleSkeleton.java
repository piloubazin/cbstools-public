package de.mpg.cbs.jist.shape;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
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
public class JistShapeSimpleSkeleton extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	
	private ParamFloat boundParam;
	private ParamFloat distParam;
	private ParamOption featureParam;
	private static final String[] featureTypes = {"signed_distance","probability_map"}; 
	
	private ParamVolume medialImage;
	private ParamVolume skelImage;
	
	private ShapeSimpleSkeleton algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Shape image"));
		inputParams.add(featureParam = new ParamOption("Shape image type", featureTypes));
		
		inputParams.add(boundParam = new ParamFloat("Boundary threshold (>0: inside, <0: outside)", -100.0f, 100.0f, 0.0f));
		inputParams.add(distParam = new ParamFloat("Skeleton threshold (>0: inside, <0: outside)", -100.0f, 100.0f, 2.0f));
		
		algorithm = new ShapeSimpleSkeleton();
		
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
		outputParams.add(medialImage = new ParamVolume("Medial Surface Image",VoxelType.UBYTE));
		outputParams.add(skelImage = new ParamVolume("Medial Curve Image",VoxelType.UBYTE));

		outputParams.setName("skeleton images");
		outputParams.setLabel("skeleton images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new ShapeSimpleSkeleton();
		
		// i/o variables
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		
		int[] dims = Interface.getDimensions(inputImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(inputImage));
		
		algorithm.setShapeImage(Interface.getFloatImage3D(inputImage));
		algorithm.setShapeImageType(featureParam.getValue());
		
		algorithm.setBoundaryThreshold(boundParam.getValue().floatValue());
		algorithm.setSkeletonThreshold(distParam.getValue().floatValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setUByteImage3D(algorithm.getMedialSurfaceImage(), dims, medialImage, name+"_medial", header);
		Interface.setUByteImage3D(algorithm.getMedialCurveImage(), dims, skelImage, name+"_skel", header);
	}
}
