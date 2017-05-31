package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.registration.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistRegistrationTargetBasedReorientation extends ProcessingAlgorithm {

	private RegistrationTargetBasedReorientation algorithm;
	
	// jist containers
	private ParamVolume headImage;
	private ParamVolume extraImage;
	private ParamVolume targetImage;
	private ParamFloat distanceParam;
	private ParamOption neckdirParam;
	
	private ParamVolume reorientImage;
	private ParamVolume locatorImage;
	
	private		static final String[]	dirtypes = {"none","-X","+X","-Y","+Y","-Z","+Z"};
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(headImage = new ParamVolume("Head Probability Image"));
		inputParams.add(extraImage = new ParamVolume("Image to Reorient"));
		inputParams.add(targetImage = new ParamVolume("Target Region Image"));
		
		inputParams.add(distanceParam = new ParamFloat("Distance to target (mm)", 0.0f, 100.0f, 3.5f));
		inputParams.add(neckdirParam = new ParamOption("Neck extension direction", dirtypes));
		
		algorithm = new RegistrationTargetBasedReorientation();
		
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
		outputParams.add(reorientImage = new ParamVolume("ReorientedImage",VoxelType.FLOAT));
		outputParams.add(locatorImage = new ParamVolume("Target Locator Image",VoxelType.UBYTE));
		
		outputParams.setName("target images");
		outputParams.setLabel("target images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions(headImage);
		float[] res = Interface.getResolutions(headImage);
		String name = Interface.getName(headImage);
		ImageHeader header = Interface.getHeader(headImage);
		
		// main algorithm
		algorithm = new RegistrationTargetBasedReorientation();
		
		algorithm.setHeadProbabilityImage(Interface.getFloatImage3D(headImage));
		algorithm.setExtraImage(Interface.getFloatImage3D(extraImage));
		algorithm.setTargetImage(Interface.getUByteImage3D(targetImage));
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);
		algorithm.setDistance_mm(distanceParam.getValue().floatValue());
		
		     if (neckdirParam.getValue().equals("-X")) algorithm.setNeckDirection(algorithm.mX);
		else if (neckdirParam.getValue().equals("+X")) algorithm.setNeckDirection(algorithm.pX);
		else if (neckdirParam.getValue().equals("-Y")) algorithm.setNeckDirection(algorithm.mY);
		else if (neckdirParam.getValue().equals("+Y")) algorithm.setNeckDirection(algorithm.pY);
		else if (neckdirParam.getValue().equals("-Z")) algorithm.setNeckDirection(algorithm.mZ);
		else if (neckdirParam.getValue().equals("+Z")) algorithm.setNeckDirection(algorithm.pZ);
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getReorientedImage(), dims, reorientImage, name+"_rot", header);
		Interface.setUByteImage3D(algorithm.getLocatorImage(), dims, locatorImage, name+"_loc", header);
	}


}
