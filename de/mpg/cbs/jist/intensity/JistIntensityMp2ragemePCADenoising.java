package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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
	private ParamVolume 	inv2Image;
	private ParamVolume 	inv2e1Image;
	private ParamVolume 	inv2e2Image;
	private ParamVolume 	inv2e3Image;
	private ParamVolume 	inv2e4Image;
	
	private ParamFloat 		cutoffParam;
	private ParamInteger    mindimParam;
	private ParamInteger    maxdimParam;
	private ParamInteger    sizeParam;
	//private ParamBoolean    separateParam;
	//private ParamBoolean    tvphsParam;
	//private ParamBoolean    tvmagParam;
	
	private ParamVolume 	denoisedInv1Image;
	private ParamVolume 	denoisedInv2Image;
	private ParamVolume 	denoisedInv2e1Image;
	private ParamVolume 	denoisedInv2e2Image;
	private ParamVolume 	denoisedInv2e3Image;
	private ParamVolume 	denoisedInv2e4Image;
	private ParamVolume 	localdimImage;
	private ParamVolume 	noisemagImage;
	private ParamVolume 	rawcomplexImage;
	private ParamVolume 	dencomplexImage;
	private ParamVolume 	eigvalImage;
	private ParamVolume 	eigvecImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inv1Image = new ParamVolume("First Inversion Image (4D:mag+phs)"));
		inputParams.add(inv2Image = new ParamVolume("Second Inversion Image (4D: N mag+ N phs)"));
		inv2Image.setMandatory(false);
		inputParams.add(inv2e1Image = new ParamVolume("Second Inversion Image Echo1 (4D:mag+phs)"));
		inv2e1Image.setMandatory(false);
		inputParams.add(inv2e2Image = new ParamVolume("Second Inversion Image Echo2 (4D:mag+phs)"));
		inv2e2Image.setMandatory(false);
		inputParams.add(inv2e3Image = new ParamVolume("Second Inversion Image Echo3 (4D:mag+phs)"));
		inv2e3Image.setMandatory(false);
		inputParams.add(inv2e4Image = new ParamVolume("Second Inversion Image Echo4 (4D:mag+phs)"));
		inv2e4Image.setMandatory(false);
		
		inputParams.add(cutoffParam = new ParamFloat("Stdev cutoff", 0.0f, 100.0f, 1.0f));
		inputParams.add(mindimParam = new ParamInteger("Minimum dimension", 0, 100, 2));
        inputParams.add(maxdimParam = new ParamInteger("Maximum dimension", -1, 100, -1));
        inputParams.add(sizeParam = new ParamInteger("Patch size", 2, 20, 5));
        //inputParams.add(separateParam = new ParamBoolean("Separate mag/phs", false));
        //inputParams.add(tvmagParam = new ParamBoolean("Magnitude TV subtraction", false));
        //inputParams.add(tvphsParam = new ParamBoolean("Phase TV subtraction", false));

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
		outputParams.add(denoisedInv2Image = new ParamVolume("Denoised Inv2 Image",null,-1,-1,-1,-1));
		denoisedInv2Image.setMandatory(false);
		outputParams.add(denoisedInv2e1Image = new ParamVolume("Denoised Inv2 E1 Image",null,-1,-1,-1,-1));
		denoisedInv2e1Image.setMandatory(false);
		outputParams.add(denoisedInv2e2Image = new ParamVolume("Denoised Inv2 E2 Image",null,-1,-1,-1,-1));
		denoisedInv2e2Image.setMandatory(false);
		outputParams.add(denoisedInv2e3Image = new ParamVolume("Denoised Inv2 E3 Image",null,-1,-1,-1,-1));
		denoisedInv2e3Image.setMandatory(false);
		outputParams.add(denoisedInv2e4Image = new ParamVolume("Denoised Inv2 E4 Image",null,-1,-1,-1,-1));
        denoisedInv2e4Image.setMandatory(false);
		outputParams.add(localdimImage = new ParamVolume("PCA Dimension Image",null,-1,-1,-1,-1));
        outputParams.add(noisemagImage = new ParamVolume("Noise Magnitude Image",null,-1,-1,-1,-1));
        outputParams.add(rawcomplexImage = new ParamVolume("Complex Image",null,-1,-1,-1,-1));
        outputParams.add(dencomplexImage = new ParamVolume("Denoised Complex Image",null,-1,-1,-1,-1));
        outputParams.add(eigvecImage = new ParamVolume("Eigenvectors Image",null,-1,-1,-1,-1));
        outputParams.add(eigvalImage = new ParamVolume("Eigenvalues Image",null,-1,-1,-1,-1));

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
		if (Interface.isValid(inv2Image)) {
		    algorithm.setSecondInversionImage(Interface.getFloatImage4D(inv2Image));    
		} else {
            algorithm.setSecondInversionEcho1Image(Interface.getFloatImage4D(inv2e1Image));
            algorithm.setSecondInversionEcho2Image(Interface.getFloatImage4D(inv2e2Image));
            algorithm.setSecondInversionEcho3Image(Interface.getFloatImage4D(inv2e3Image));
            algorithm.setSecondInversionEcho4Image(Interface.getFloatImage4D(inv2e4Image));
		}
		
		algorithm.setStdevCutoff(cutoffParam.getValue().floatValue());
		algorithm.setMinimumDimension(mindimParam.getValue().intValue());
		algorithm.setMaximumDimension(maxdimParam.getValue().intValue());
		algorithm.setPatchSize(sizeParam.getValue().intValue());
		//algorithm.setProcessSeparately(separateParam.getValue().booleanValue());
		//algorithm.setMagnitudeTVSubtraction(tvmagParam.getValue().booleanValue());
		//algorithm.setPhaseTVSubtraction(tvphsParam.getValue().booleanValue());
		
		System.out.println("run the algorithm");
		algorithm.execute();

		System.out.println("save denoised images");
		Interface.setFloatImage4D(algorithm.getFirstInversionImage(), dims, 2, denoisedInv1Image, Interface.getName(inv1Image)+"_lpca", Interface.getHeader(inv1Image));
		if (Interface.isValid(inv2Image)) {
		    Interface.setFloatImage4D(algorithm.getSecondInversionImage(), dims, 8, denoisedInv2Image, Interface.getName(inv2Image)+"_lpca", Interface.getHeader(inv2Image));
		} else {
            Interface.setFloatImage4D(algorithm.getSecondInversionEcho1Image(), dims, 2, denoisedInv2e1Image, Interface.getName(inv2e1Image)+"_lpca", Interface.getHeader(inv2e1Image));
            Interface.setFloatImage4D(algorithm.getSecondInversionEcho2Image(), dims, 2, denoisedInv2e2Image, Interface.getName(inv2e2Image)+"_lpca", Interface.getHeader(inv2e2Image));
            Interface.setFloatImage4D(algorithm.getSecondInversionEcho3Image(), dims, 2, denoisedInv2e3Image, Interface.getName(inv2e3Image)+"_lpca", Interface.getHeader(inv2e3Image));
            Interface.setFloatImage4D(algorithm.getSecondInversionEcho4Image(), dims, 2, denoisedInv2e4Image, Interface.getName(inv2e4Image)+"_lpca", Interface.getHeader(inv2e4Image));
        }
		Interface.setFloatImage3D(algorithm.getLocalDimensionImage(), dims, localdimImage, Interface.getName(inv1Image)+"_ldim", Interface.getHeader(inv1Image));
		Interface.setFloatImage3D(algorithm.getNoiseMagnitudeImage(), dims, noisemagImage, Interface.getName(inv1Image)+"_lerr", Interface.getHeader(inv1Image));
		Interface.setFloatImage4D(algorithm.getRawComplexImage(), dims, 10, rawcomplexImage, Interface.getName(inv1Image)+"_lrcx", Interface.getHeader(inv1Image));
		Interface.setFloatImage4D(algorithm.getDenComplexImage(), dims, 10, dencomplexImage, Interface.getName(inv1Image)+"_ldcx", Interface.getHeader(inv1Image));
		Interface.setFloatImage4D(algorithm.getEigenvectorImage(), dims, 10, eigvecImage, Interface.getName(inv1Image)+"_lvec", Interface.getHeader(inv1Image));
		Interface.setFloatImage4D(algorithm.getEigenvalueImage(), dims, 10, eigvalImage, Interface.getName(inv1Image)+"_lval", Interface.getHeader(inv1Image));
	}


}
