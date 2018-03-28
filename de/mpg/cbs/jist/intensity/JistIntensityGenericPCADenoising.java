package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
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
public class JistIntensityGenericPCADenoising extends ProcessingAlgorithm {

	private IntensityGenericPCADenoising algorithm;
	
	// jist containers
	private ParamVolumeCollection 	inputImages;
	
	private ParamFloat 		cutoffParam;
	private ParamInteger    mindimParam;
	private ParamInteger    maxdimParam;
	private ParamInteger    sizeParam;
	
	private ParamVolumeCollection 	denoisedImages;
	private ParamVolume 	localdimImage;
	private ParamVolume 	noisemagImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImages = new ParamVolumeCollection("Input Images (3D or 4D)"));
		
		inputParams.add(cutoffParam = new ParamFloat("Stdev cutoff", 0.0f, 100.0f, 2.3f));
		inputParams.add(mindimParam = new ParamInteger("Minimum dimension", 0, 100, 2));
        inputParams.add(maxdimParam = new ParamInteger("Maximum dimension", -1, 100, -1));
        inputParams.add(sizeParam = new ParamInteger("Patch size", 3, 20, 4));

		algorithm = new IntensityGenericPCADenoising();
		
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
		outputParams.add(denoisedImages = new ParamVolumeCollection("Denoised Images",null,-1,-1,-1,-1));
		outputParams.add(localdimImage = new ParamVolume("PCA Dimension Image",null,-1,-1,-1,-1));
        outputParams.add(noisemagImage = new ParamVolume("Noise Magnitude Image",null,-1,-1,-1,-1));

		outputParams.setName("denoised images");
		outputParams.setLabel("denoised images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
	    // retrieve images
	    int nimg = Interface.getVolumeCollectionSize(inputImages);
	    int nvol = 0;
	    for (int n=0;n<nimg;n++) {
	        nvol += Interface.getComponents(Interface.getVolumeFromCollection(n, inputImages));
	    }
	    
		// i/o variables
		int[] dims = Interface.getDimensions(Interface.getVolumeFromCollection(0, inputImages));
		float[] res = Interface.getResolutions(Interface.getVolumeFromCollection(0, inputImages));
		String[] names = new String[nimg];
		for (int n=0;n<nimg;n++) names[n] = Interface.getName(Interface.getVolumeFromCollection(n, inputImages));
		ImageHeader[] headers = new ImageHeader[nimg];
		for (int n=0;n<nimg;n++) headers[n] = Interface.getHeader(Interface.getVolumeFromCollection(n, inputImages));
		
		// main algorithm
		algorithm = new IntensityGenericPCADenoising();
		
		// important: set image size before loading images
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		System.out.println("load images");
		algorithm.setImageNumber(nvol);
		int v=0;
		for (int n=0;n<nimg;n++) {
		    if (Interface.isImage4D(Interface.getVolumeFromCollection(n, inputImages))) {
		        for (int c=0;c<Interface.getComponents(Interface.getVolumeFromCollection(n, inputImages));c++) {
		            algorithm.setImageAt(v, Interface.getFloatImageComponent4D(Interface.getVolumeFromCollection(n, inputImages), c));
		            v++;
		        }
		    } else {
                algorithm.setImageAt(v, Interface.getFloatImage3D(Interface.getVolumeFromCollection(n, inputImages)));
                v++;
            }
        }
		
		algorithm.setStdevCutoff(cutoffParam.getValue().floatValue());
		algorithm.setMinimumDimension(mindimParam.getValue().intValue());
		algorithm.setMaximumDimension(maxdimParam.getValue().intValue());
		algorithm.setPatchSize(sizeParam.getValue().intValue());
		
		System.out.println("run the algorithm");
		algorithm.execute();

		System.out.println("save denoised images");
		v = 0;
		for (int n=0;n<nimg;n++) {
		    if (Interface.isImage4D(Interface.getVolumeFromCollection(n, inputImages))) {
		        Interface.addFloatImage4D(algorithm.getDenoisedImageAt(v,v+Interface.getComponents(Interface.getVolumeFromCollection(n, inputImages))), dims, Interface.getComponents(Interface.getVolumeFromCollection(n, inputImages)), denoisedImages, names[n]+"_lpca", headers[n]);
		        v+=Interface.getComponents(Interface.getVolumeFromCollection(n, inputImages));
		    } else {
		        Interface.addFloatImage3D(algorithm.getDenoisedImageAt(v,v+1), dims, denoisedImages, names[n]+"_lpca", headers[n]);
                v++;
            }
        }
		Interface.setFloatImage3D(algorithm.getLocalDimensionImage(), dims, localdimImage, names[0]+"_ldim", headers[0]);
		Interface.setFloatImage3D(algorithm.getNoiseMagnitudeImage(), dims, noisemagImage, names[0]+"_lerr", headers[0]);
		//Interface.setFloatImage4D(algorithm.getNoiseMagnitudeImage(), dims, nimg+2, noisemagImage, names[0]+"_lerr", headers[0]);
	}


}
