package de.mpg.cbs.jist.intensity;

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
import de.mpg.cbs.core.intensity.IntensityMp2rageMasking;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityMp2rageMasking extends ProcessingAlgorithm {

    private IntensityMp2rageMasking algorithm;
    
	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume t1mapImage;
	private ParamVolume isoImage;
	private ParamOption distribParam;
	private ParamOption maskingParam;
	
	private ParamVolume probaImage;
	private ParamVolume maskImage;
	private ParamVolume t1maskImage;
	private ParamVolume isomaskImage;
	
	// parameters
	private ParamBoolean	skip0Param;
	private ParamBoolean	noiterParam;
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) Image"));
		
		inputParams.add(t1mapImage = new ParamVolume("Quantitative T1 Map (T1_Images) Image"));
		inputParams.add(isoImage = new ParamVolume("T1-weighted (UNI) Image"));
		t1mapImage.setMandatory(false);
		isoImage.setMandatory(false);
		
		inputParams.add(distribParam = new ParamOption("Background distribution", IntensityMp2rageMasking.distribs));
		distribParam.setDescription("Model distribution for background noise (default is half-normal, exponential is more stringent).");
		distribParam.setValue(IntensityMp2rageMasking.distrib);
		
		inputParams.add(skip0Param = new ParamBoolean("Skip zero values", true));
		inputParams.add(noiterParam = new ParamBoolean("non-iterative estimate", false));
		
		inputParams.add(maskingParam = new ParamOption("Masking method", IntensityMp2rageMasking.maskings));
		maskingParam.setDescription("Whether to use a binary threshold or a weighted average based on the probability.");
		maskingParam.setValue(IntensityMp2rageMasking.masking);
		
		algorithm = new IntensityMp2rageMasking();
		
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
		outputParams.add(probaImage = new ParamVolume("Signal Proba Image",VoxelType.FLOAT));
		outputParams.add(maskImage = new ParamVolume("Signal Mask Image",VoxelType.UBYTE));
		outputParams.add(t1maskImage = new ParamVolume("Masked T1_Map Image",VoxelType.FLOAT));
		outputParams.add(isomaskImage = new ParamVolume("Masked T1-weighted Image",VoxelType.FLOAT));
		
		t1maskImage.setMandatory(false);
		isomaskImage.setMandatory(false);
		
		outputParams.setName("masking images");
		outputParams.setLabel("masking images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		int[] dims = Interface.getDimensions(inv2Image);
		String name = Interface.getName(inv2Image);
		ImageHeader header = Interface.getHeader(inv2Image);

		// main algorithm
		algorithm = new IntensityMp2rageMasking();
		
		algorithm.setSecondInversionImage(Interface.getFloatImage3D(inv2Image));
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(inv2Image));
		
		if (Interface.isValid(t1mapImage))
		    algorithm.setQuantitativeT1MapImage(Interface.getFloatImage3D(t1mapImage));
		if (Interface.isValid(isoImage))
		    algorithm.setT1weightedImage(Interface.getFloatImage3D(isoImage));
		
		algorithm.setBackgroundDistribution(distribParam.getValue());
		algorithm.setSkipZeroValues(skip0Param.getValue().booleanValue());
		algorithm.setNonIterativeEstimate(noiterParam.getValue().booleanValue());
		algorithm.setMaskingMethod(maskingParam.getValue());
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getSignalProbaImage(), dims, probaImage, name+"_mp2msk_proba", header);
		Interface.setIntegerImage3D(algorithm.getSignalMaskImage(), dims, maskImage, name+"_mp2msk_mask", header);
		if (Interface.isValid(t1mapImage)) {
		    String t1name = Interface.getName(t1mapImage);
		    Interface.setFloatImage3D(algorithm.getMaskedT1MapImage(), dims, t1maskImage, t1name+"_mp2msk_masked", header);
		}
		if (Interface.isValid(isoImage)) {
		    String isoname = Interface.getName(isoImage);
		    Interface.setFloatImage3D(algorithm.getMaskedT1weightedImage(), dims, isomaskImage, isoname+"_mp2msk_masked", header);
		}
		return;
	}
}
