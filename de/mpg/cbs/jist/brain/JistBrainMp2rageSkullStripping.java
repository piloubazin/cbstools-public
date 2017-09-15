package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.core.brain.BrainMp2rageSkullStripping;
import de.mpg.cbs.utilities.Interface;
import de.mpg.cbs.utilities.References;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMp2rageSkullStripping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume t1mapImage;
	private ParamVolume isoImage;
	/*
	private ParamOption	 bgParam;
	private String		bgType = "exponential";
	private 	static final String[]	bgTypes = {"exponential","half-normal"};
	*/
	private		static final byte	HNORM = 101;
	private		static final byte	EXP = 102;
	
	private ParamBoolean	skip0Param;
	private ParamVolume brainmaskImage;
	private ParamVolume t1maskImage;
	private ParamVolume isomaskImage;
	private ParamVolume inv2maskImage;
	
	private BrainMp2rageSkullStripping algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) Image"));
		
		inputParams.add(t1mapImage = new ParamVolume("T1 Map (T1_Images) Image (opt)"));
		inputParams.add(isoImage = new ParamVolume("T1-weighted (UNI) Image (opt)"));
		t1mapImage.setMandatory(false);
		isoImage.setMandatory(false);
		/*
		inputParams.add(bgParam = new ParamOption("Background noise model", bgTypes));
		bgParam.setValue(bgType);
		*/
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", false));

		algorithm = new BrainMp2rageSkullStripping();
		
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
		outputParams.add(brainmaskImage = new ParamVolume("Brain Mask Image",VoxelType.UBYTE));
		outputParams.add(t1maskImage = new ParamVolume("Masked T1 Map Image",VoxelType.FLOAT));
		outputParams.add(isomaskImage = new ParamVolume("Masked T1-weighted Image",VoxelType.FLOAT));
		outputParams.add(inv2maskImage = new ParamVolume("Masked Inv2 Image",VoxelType.FLOAT));
		
		t1maskImage.setMandatory(false);
		isomaskImage.setMandatory(false);
		
		outputParams.setName("skull stripped images");
		outputParams.setLabel("skull stripped images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String inv2name = Interface.getName(inv2Image);
		ImageHeader header = Interface.getHeader(inv2Image);
		int[] dims = Interface.getDimensions(inv2Image);
		float[] res = Interface.getResolutions(inv2Image);
		
		String t1name = null;
		if (Interface.isValid(t1mapImage)) t1name = Interface.getName(t1mapImage);
		String isoname = null;
		if (Interface.isValid(isoImage)) isoname = Interface.getName(isoImage);
		
		// import the image data into 1D arrays
		algorithm = new BrainMp2rageSkullStripping();
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setSecondInversionImage(Interface.getFloatImage3D(inv2Image));
		if (Interface.isValid(t1mapImage)) algorithm.setT1MapImage(Interface.getFloatImage3D(t1mapImage));
		if (Interface.isValid(isoImage)) algorithm.setT1weightedImage(Interface.getFloatImage3D(isoImage));
		
		algorithm.setSkipZeroValues(skip0Param.getValue().booleanValue());

		algorithm.execute();
		
		Interface.setUByteImage3D(algorithm.getBrainMaskImage(), dims, brainmaskImage, inv2name+"_stripmask", header);
		
		Interface.setFloatImage3D(algorithm.getMaskedSecondInversionImage(), dims, inv2maskImage, inv2name+"_strip", header);
		if (t1mapImage.getImageData() != null)
			Interface.setFloatImage3D(algorithm.getMaskedT1MapImage(), dims, t1maskImage, t1name+"_strip", header);
		if (isoImage.getImageData() != null)
			Interface.setFloatImage3D(algorithm.getMaskedT1weightedImage(), dims, isomaskImage, isoname+"_strip", header);
	}


}
