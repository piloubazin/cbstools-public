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
public class JistIntensityComplexPCADenoising extends ProcessingAlgorithm {

	private IntensityComplexPCADenoising algorithm;
	
	// jist containers
	private ParamVolume 	magImage;
	private ParamVolume 	phsImage;
	private ParamVolume 	combiImage;
	
	private ParamFloat 		cutoffParam;
	private ParamInteger    mindimParam;
	private ParamInteger    maxdimParam;
	private ParamInteger    sizeParam;
	private ParamInteger    windowParam;
	//private ParamBoolean    separateParam;
	//private ParamBoolean    tvphsParam;
	//private ParamBoolean    tvmagParam;
	
	private ParamVolume 	denoisedmagImage;
	private ParamVolume 	denoisedphsImage;
	private ParamVolume 	denoisedcombiImage;
	private ParamVolume 	localdimImage;
	private ParamVolume 	noisefitImage;
//	private ParamVolume 	rawcomplexImage;
//	private ParamVolume 	dencomplexImage;
//	private ParamVolume 	eigvalImage;
//	private ParamVolume 	eigvecImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(magImage = new ParamVolume("Magnitude Image (4D)"));
		magImage.setMandatory(false);
		inputParams.add(phsImage = new ParamVolume("Phase Image (4D)"));
		phsImage.setMandatory(false);
		inputParams.add(combiImage = new ParamVolume("Combined Image (4D: N mag+ N phs)"));
		combiImage.setMandatory(false);
		
		inputParams.add(cutoffParam = new ParamFloat("Stdev cutoff", 0.0f, 100.0f, 1.0f));
		inputParams.add(mindimParam = new ParamInteger("Minimum dimension", 0, 100, 2));
        inputParams.add(maxdimParam = new ParamInteger("Maximum dimension", -1, 100, -1));
        inputParams.add(sizeParam = new ParamInteger("Patch size", 2, 20, 5));
        inputParams.add(windowParam = new ParamInteger("Window size", -1, 1000, 5));
        //inputParams.add(separateParam = new ParamBoolean("Separate mag/phs", false));
        //inputParams.add(tvmagParam = new ParamBoolean("Magnitude TV subtraction", false));
        //inputParams.add(tvphsParam = new ParamBoolean("Phase TV subtraction", false));

		algorithm = new IntensityComplexPCADenoising();
		
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
		outputParams.add(denoisedmagImage = new ParamVolume("Denoised Magnitude Image",null,-1,-1,-1,-1));
		denoisedmagImage.setMandatory(false);
		outputParams.add(denoisedphsImage = new ParamVolume("Denoised Phase Image",null,-1,-1,-1,-1));
		denoisedphsImage.setMandatory(false);
		outputParams.add(denoisedcombiImage = new ParamVolume("Denoised Combined Image",null,-1,-1,-1,-1));
		denoisedcombiImage.setMandatory(false);
		
		outputParams.add(localdimImage = new ParamVolume("PCA Dimension Image",null,-1,-1,-1,-1));
        outputParams.add(noisefitImage = new ParamVolume("Noise Fit Image",null,-1,-1,-1,-1));
        //outputParams.add(rawcomplexImage = new ParamVolume("Complex Image",null,-1,-1,-1,-1));
        //outputParams.add(dencomplexImage = new ParamVolume("Denoised Complex Image",null,-1,-1,-1,-1));
        //outputParams.add(eigvecImage = new ParamVolume("Eigenvectors Image",null,-1,-1,-1,-1));
        //outputParams.add(eigvalImage = new ParamVolume("Eigenvalues Image",null,-1,-1,-1,-1));

		outputParams.setName("denoised images");
		outputParams.setLabel("denoised images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = null;
		int nimg = 0;
		float[] res = null;
		String name = null;
		ImageHeader header = null;
		if (Interface.isValid(magImage)) {
		    dims = Interface.getDimensions(magImage);
		    nimg = Interface.getComponents(magImage);
		    res = Interface.getResolutions(magImage);
		    name = Interface.getName(magImage);
		    header = Interface.getHeader(magImage);
		} else {
		    dims = Interface.getDimensions(combiImage);
		    nimg = Interface.getComponents(combiImage)/2;
		    res = Interface.getResolutions(combiImage);
		    name = Interface.getName(combiImage);
		    header = Interface.getHeader(combiImage);
		}
		
		// main algorithm
		algorithm = new IntensityComplexPCADenoising();
		
		// important: set image size before loading images
		algorithm.setDimensions(dims);
		algorithm.setImageNumber(nimg);
		algorithm.setResolutions(res);

		System.out.println("load images");
		if (Interface.isValid(magImage)) {
		    algorithm.setMagnitudeImages(Interface.getFloatImage4D(magImage));
		    if (Interface.isValid(phsImage)) algorithm.setPhaseImages(Interface.getFloatImage4D(phsImage));
		} else {
		    algorithm.setMagnitudeAndPhaseImage(Interface.getFloatImage4D(combiImage));    
		}
		
		algorithm.setStdevCutoff(cutoffParam.getValue().floatValue());
		algorithm.setMinimumDimension(mindimParam.getValue().intValue());
		algorithm.setMaximumDimension(maxdimParam.getValue().intValue());
		algorithm.setPatchSize(sizeParam.getValue().intValue());
		algorithm.setWindowSize(windowParam.getValue().intValue());
		//algorithm.setProcessSeparately(separateParam.getValue().booleanValue());
		//algorithm.setMagnitudeTVSubtraction(tvmagParam.getValue().booleanValue());
		//algorithm.setPhaseTVSubtraction(tvphsParam.getValue().booleanValue());
		
		int nsample = Numerics.ceil(nimg/Numerics.floor(windowParam.getValue().intValue()/2));
		
		System.out.println("run the algorithm");
		algorithm.execute();

		System.out.println("save denoised images");
		if (Interface.isValid(magImage)) {
            Interface.setFloatImage4D(algorithm.getDenoisedMagnitudeImages(), dims, nimg, denoisedmagImage, Interface.getName(magImage)+"_lpca", Interface.getHeader(magImage));
           if (Interface.isValid(phsImage)) Interface.setFloatImage4D(algorithm.getDenoisedPhaseImages(), dims, nimg, denoisedphsImage, Interface.getName(phsImage)+"_lpca", Interface.getHeader(phsImage));
		} else {
            Interface.setFloatImage4D(algorithm.getDenoisedMagnitudeAndPhaseImage(), dims, 2*nimg, denoisedcombiImage, Interface.getName(combiImage)+"_lpca", Interface.getHeader(combiImage));
        }
		Interface.setFloatImage4D(algorithm.getLocalDimensionImage(), dims, nimg, localdimImage, name+"_ldim", header);
		Interface.setFloatImage4D(algorithm.getNoiseFitImage(), dims, nimg, noisefitImage, name+"_lerr", header);
		//Interface.setFloatImage4D(algorithm.getRawComplexImage(), dims, 10, rawcomplexImage, Interface.getName(inv1Image)+"_lrcx", Interface.getHeader(inv1Image));
		//Interface.setFloatImage4D(algorithm.getDenComplexImage(), dims, 10, dencomplexImage, Interface.getName(inv1Image)+"_ldcx", Interface.getHeader(inv1Image));
		//Interface.setFloatImage4D(algorithm.getEigenvectorImage(), dims, 10, eigvecImage, Interface.getName(inv1Image)+"_lvec", Interface.getHeader(inv1Image));
		//Interface.setFloatImage4D(algorithm.getEigenvalueImage(), dims, 10, eigvalImage, Interface.getName(inv1Image)+"_lval", Interface.getHeader(inv1Image));
	}


}
