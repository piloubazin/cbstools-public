package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.core.brain.BrainMp2rageDuraEstimation;


/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMp2rageDuraEstimation extends ProcessingAlgorithm {

    BrainMp2rageDuraEstimation algorithm;
    
	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private ParamOption outParam;
	private ParamFloat distParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) image"));
		inputParams.add(maskImage = new ParamVolume("Skull stripping mask"));
		inputParams.add(distParam = new ParamFloat("Distance to background (mm)", -1E10f, 1E10f, 5.0f));
		inputParams.add(outParam = new ParamOption("Output type", BrainMp2rageDuraEstimation.outTypes));
		outParam.setDescription("Outputs an estimate of the dura / CSF boundary or an estimate of the entire dura region.");
		
		algorithm = new BrainMp2rageDuraEstimation();
		
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
		outputParams.add(resultImage = new ParamVolume("Dura Image",VoxelType.FLOAT));
		outputParams.setName("dura_image");
		outputParams.setLabel("Dura Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(inv2Image);
		ImageHeader header = Interface.getHeader(inv2Image);
		int[] dims = Interface.getDimensions(inv2Image);
		float[] res = Interface.getResolutions(inv2Image);
		
		// main algorithm
		algorithm = new BrainMp2rageDuraEstimation();

		algorithm.setSecondInversionImage(Interface.getFloatImage3D(inv2Image));
		algorithm.setSkullStrippingMask(Interface.getIntegerImage3D(maskImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);
		
		algorithm.setDistanceToBackground_mm(distParam.getValue().floatValue());
		algorithm.setOutputType(outParam.getValue());
		
		algorithm.execute();
	
		Interface.setFloatImage3D(algorithm.getDuraImage(), dims, resultImage, name+"_dura_img", header);
		
	}
}
