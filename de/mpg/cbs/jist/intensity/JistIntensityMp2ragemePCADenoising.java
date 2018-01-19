package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
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
import de.mpg.cbs.core.intensity.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityMp2ragemePCADenoising extends ProcessingAlgorithm {

	private IntensityMp2ragemePCADenoising algorithm;
	
	// jist containers
	private ParamVolume 	inv1Image;
	private ParamVolume 	inv2e1Image;
	private ParamVolume 	inv2e2Image;
	private ParamVolume 	inv2e3Image;
	private ParamVolume 	inv2e4Image;
	
	private ParamFloat 		cutoffParam;
	private ParamInteger    mindimParam;
	
	private ParamVolume 	denoisedInv1Image;
	private ParamVolume 	denoisedInv2e1Image;
	private ParamVolume 	denoisedInv2e2Image;
	private ParamVolume 	denoisedInv2e3Image;
	private ParamVolume 	denoisedInv2e4Image;
	private ParamVolume 	localdimImage;
	private ParamVolume 	noisemagImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inv1Image = new ParamVolume("First Inversion Image (4D:mag+phs)"));
		inputParams.add(inv2e1Image = new ParamVolume("Second Inversion Image Echo1 (4D:mag+phs)"));
		inputParams.add(inv2e2Image = new ParamVolume("Second Inversion Image Echo2 (4D:mag+phs)"));
		inputParams.add(inv2e3Image = new ParamVolume("Second Inversion Image Echo3 (4D:mag+phs)"));
		inputParams.add(inv2e4Image = new ParamVolume("Second Inversion Image Echo4 (4D:mag+phs)"));
		
		inputParams.add(cutoffParam = new ParamFloat("Stdev cutoff", 0.0f, 100.0f, 2.3f));
		inputParams.add(mindimParam = new ParamInteger("Minimum dimension", 0, 100, 3));

		algorithm = new IntensityMp2ragemePCADenoising();
		
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
		outputParams.add(denoisedInv1Image = new ParamVolume("Denoised Inv1 Image",null,-1,-1,-1,-1));
		outputParams.add(denoisedInv2e1Image = new ParamVolume("Denoised Inv2 E1 Image",null,-1,-1,-1,-1));
		outputParams.add(denoisedInv2e2Image = new ParamVolume("Denoised Inv2 E2 Image",null,-1,-1,-1,-1));
		outputParams.add(denoisedInv2e3Image = new ParamVolume("Denoised Inv2 E3 Image",null,-1,-1,-1,-1));
		outputParams.add(denoisedInv2e4Image = new ParamVolume("Denoised Inv2 E4 Image",null,-1,-1,-1,-1));
        outputParams.add(localdimImage = new ParamVolume("PCA Dimension Image",null,-1,-1,-1,-1));
        outputParams.add(noisemagImage = new ParamVolume("Noise Magnitude Image",null,-1,-1,-1,-1));

		outputParams.setName("denoised images");
		outputParams.setLabel("denoised images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions(inv1Image);
		float[] res = Interface.getResolutions(inv1Image);
		String name = Interface.getName(inv1Image);
		ImageHeader header = Interface.getHeader(inv1Image);
		
		// main algorithm
		algorithm = new IntensityMp2ragemePCADenoising();
		
		// important: set image size before loading images
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		System.out.println("load images");
		algorithm.setFirstInversionImage(Interface.getFloatImage4D(inv1Image));
		algorithm.setSecondInversionEcho1Image(Interface.getFloatImage4D(inv2e1Image));
		algorithm.setSecondInversionEcho2Image(Interface.getFloatImage4D(inv2e2Image));
		algorithm.setSecondInversionEcho3Image(Interface.getFloatImage4D(inv2e3Image));
		algorithm.setSecondInversionEcho4Image(Interface.getFloatImage4D(inv2e4Image));
		
		algorithm.setStdevCutoff(cutoffParam.getValue().floatValue());
		algorithm.setMinimumDimension(mindimParam.getValue().intValue());
		
		System.out.println("run the algorithm");
		algorithm.execute();

		System.out.println("save denoised images");
		Interface.setFloatImage4D(algorithm.getFirstInversionImage(), dims, 2, denoisedInv1Image, Interface.getName(inv1Image)+"_lpca", Interface.getHeader(inv1Image));
		Interface.setFloatImage4D(algorithm.getSecondInversionEcho1Image(), dims, 2, denoisedInv2e1Image, Interface.getName(inv2e1Image)+"_lpca", Interface.getHeader(inv2e1Image));
		Interface.setFloatImage4D(algorithm.getSecondInversionEcho2Image(), dims, 2, denoisedInv2e2Image, Interface.getName(inv2e2Image)+"_lpca", Interface.getHeader(inv2e2Image));
		Interface.setFloatImage4D(algorithm.getSecondInversionEcho3Image(), dims, 2, denoisedInv2e3Image, Interface.getName(inv2e3Image)+"_lpca", Interface.getHeader(inv2e3Image));
		Interface.setFloatImage4D(algorithm.getSecondInversionEcho4Image(), dims, 2, denoisedInv2e4Image, Interface.getName(inv2e4Image)+"_lpca", Interface.getHeader(inv2e4Image));
		Interface.setFloatImage3D(algorithm.getLocalDimensionImage(), dims, localdimImage, Interface.getName(inv1Image)+"_ldim", Interface.getHeader(inv1Image));
		Interface.setFloatImage3D(algorithm.getNoiseMagnitudeImage(), dims, noisemagImage, Interface.getName(inv1Image)+"_lerr", Interface.getHeader(inv1Image));
	}


}
