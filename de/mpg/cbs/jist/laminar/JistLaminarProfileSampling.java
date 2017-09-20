package de.mpg.cbs.jist.laminar;

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
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.laminar.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileSampling extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume intensityImage;
	private ParamVolume maskImage;
	private ParamOption interpParam;
	private static final String[] interpTypes = {"linear", "nearest"};
		
	private ParamVolume mappedImage;
	private ParamVolume mappedmaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private LaminarProfileSampling algorithm;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		imageParams.add(intensityImage = new ParamVolume("Intensity Image",null,-1,-1,-1,-1));
		imageParams.add(maskImage = new ParamVolume("Cortex Mask (opt)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		imageParams.add(interpParam = new ParamOption("Interpolation", interpTypes));
		
		inputParams.add(imageParams);
			
		algorithm = new LaminarProfileSampling();
		
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
		outputParams.add(mappedImage = new ParamVolume("Profile-mapped Intensity Image",null,-1,-1,-1,-1));
		
		outputParams.add(mappedmaskImage = new ParamVolume("Profile 4D Mask",null,-1,-1,-1,-1));
		
		outputParams.setName("layers images");
		outputParams.setLabel("layers images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.getName(intensityImage);
		ImageHeader header = Interface.getHeader(layersImage);
		int[] dims = Interface.getDimensions4D(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		
		// main algorithm
		algorithm = new LaminarProfileSampling();
		
		algorithm.setProfileSurfaceImage(Interface.getFloatImage4D(layersImage));
		algorithm.setIntensityImage(Interface.getFloatImage3D(intensityImage));
		algorithm.setCortexMask(Interface.getUByteImage3D(maskImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setInterpolation(interpParam.getValue());
		algorithm.execute();
		
		// output
		String imgname = intensityImage.getImageData().getName();
		
		Interface.setFloatImage4D(algorithm.getProfileMappedIntensityImage(), dims, dims[3], mappedImage, name+"_profiles", header);
		
		Interface.setUByteImage4D(algorithm.getProfile4Dmask(), dims, dims[3], mappedmaskImage, name+"_4dmask", header);
	}


}
