package de.mpg.cbs.jist.surface;

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
import de.mpg.cbs.core.surface.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceLevelsetBoundary extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume lvlImage;
	
	private ParamFloat levelParam;
	
	private ParamVolume boundaryImage;
	
	private SurfaceLevelsetBoundary algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(lvlImage = new ParamVolume("Level Set Image"));
		inputParams.add(levelParam = new ParamFloat("Level Set Boundary Value", -100.0f, 100.0f, 0.0f));
			
		algorithm = new SurfaceLevelsetBoundary();
		
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
		outputParams.add(boundaryImage = new ParamVolume("Boundary Image",VoxelType.INT));

		outputParams.setName("boundary image");
		outputParams.setLabel("boundary image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(lvlImage);
		ImageHeader header = Interface.getHeader(lvlImage);
		int[] dims = Interface.getDimensions(lvlImage);
		float[] res = Interface.getResolutions(lvlImage);
		
		// main algorithm
		algorithm = new SurfaceLevelsetBoundary();
		
		algorithm.setLevelSetImage(Interface.getFloatImage3D(lvlImage));
		algorithm.setLevelSetBoundaryValue(levelParam.getValue().floatValue());
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.execute();
		
		// output
		String imgname = lvlImage.getImageData().getName();
		
		Interface.setIntegerImage3D(algorithm.getBoundaryImage(), dims, boundaryImage, name+"_bound", header);
	}


}
