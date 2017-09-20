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
public class JistSurfaceProbabilityToLevelset extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume lvlImage;
	
	private ParamFloat scaleParam;
	
	private ParamVolume probaImage;
	
	private SurfaceProbabilityToLevelset algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(probaImage = new ParamVolume("Probability Image"));
		inputParams.add(scaleParam = new ParamFloat("Scale (mm)", 0.0f, 100.0f, 5.0f));
			
		algorithm = new SurfaceProbabilityToLevelset();
		
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
		outputParams.add(lvlImage = new ParamVolume("Level Set Image",VoxelType.FLOAT));

		outputParams.setName("levelset image");
		outputParams.setLabel("levelset image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(probaImage);
		ImageHeader header = Interface.getHeader(probaImage);
		int[] dims = Interface.getDimensions(probaImage);
		float[] res = Interface.getResolutions(probaImage);
		
		// main algorithm
		algorithm = new SurfaceProbabilityToLevelset();
		
		algorithm.setProbabilityImage(Interface.getFloatImage3D(probaImage));
		algorithm.setScale_mm(scaleParam.getValue().floatValue());
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.execute();
		
		// output
		String imgname = probaImage.getImageData().getName();
		
		Interface.setFloatImage3D(algorithm.getLevelSetImage(), dims, lvlImage, name+"_P2L", header);
	}


}
