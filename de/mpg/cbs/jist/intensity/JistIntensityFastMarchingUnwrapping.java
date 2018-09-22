package de.mpg.cbs.jist.intensity;

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
import de.mpg.cbs.core.intensity.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityFastMarchingUnwrapping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	//private ParamVolume magImage;
	private ParamVolume maskImage;
	private ParamFloat ratioParam;
	private ParamInteger quadrantParam;
	private ParamBoolean posttvParam;
		
	private ParamVolume correctImage;
	private ParamVolume countImage;
	private ParamVolume scoreImage;
	
	private IntensityFastMarchingUnwrapping algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Phase image"));
		inputParams.add(maskImage = new ParamVolume("Mask image (opt)"));
		maskImage.setMandatory(false);
		inputParams.add(quadrantParam = new ParamInteger("Seeding factor", 1, 10, 3));
		inputParams.add(posttvParam = new ParamBoolean("Background removal", false));
		inputParams.add(ratioParam = new ParamFloat("Background scaling ratio", 0.0f, 10.0f, 0.33f));
				
		algorithm = new IntensityFastMarchingUnwrapping();
		
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
		outputParams.add(correctImage = new ParamVolume("Unwrapped Phase Image",null));
		outputParams.add(countImage = new ParamVolume("Wrap Count Image",null));
		outputParams.add(scoreImage = new ParamVolume("Phase score Image",null));
		
		outputParams.setName("unwrapped images");
		outputParams.setLabel("unwrapped images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new IntensityFastMarchingUnwrapping();
		
		// i/o variables
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		
		int[] dims = Interface.getDimensions(inputImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(inputImage));
		
		algorithm.setPhaseImage(Interface.getFloatImage3D(inputImage));
		if (Interface.isValid(maskImage))
            algorithm.setMaskImage(Interface.getIntegerImage3D(maskImage));
		
        algorithm.setQuadrantNumber(quadrantParam.getValue().intValue());
        
        if (posttvParam.getValue().booleanValue()) algorithm.setTVPostProcessing("TV-residuals");
        else algorithm.setTVPostProcessing("none");
        
        algorithm.setTVScale(ratioParam.getValue().floatValue());
        
 		algorithm.execute();
		
		// outputs
		Interface.setFloatImage3D(algorithm.getCorrectedImage(), dims, correctImage, name+"_uwp", header);
		Interface.setIntegerImage3D(algorithm.getCountImage(), dims, countImage, name+"_uwc", header);
		Interface.setFloatImage3D(algorithm.getScoreImage(), dims, scoreImage, name+"_uws", header);
	}
}
