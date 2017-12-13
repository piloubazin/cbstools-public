package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
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

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;


import gov.nih.mipav.model.structures.jama.*;
import org.apache.commons.math3.util.FastMath;


public class JistLaminarIterativeSmoothing extends ProcessingAlgorithm{
	private ParamVolume intensityImage;
	private ParamVolume layersImage;
	private ParamVolume maskImage;
	
	private ParamFloat fwhmParam;
	
	private ParamVolume smoothIntensityImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final float HASQ3 = (float)(FastMath.sqrt(3.0f)/2.0f);

	private LaminarIterativeSmoothing algorithm;

	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(intensityImage = new ParamVolume("Intensity Image",null,-1,-1,-1,-1));
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("ROI Mask (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(fwhmParam = new ParamFloat("FWHM (mm)", 0.0f, 50.0f, 5.0f));
		
		algorithm = new LaminarIterativeSmoothing();
		
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
		outputParams.add(smoothIntensityImage = new ParamVolume("Smoothed Data",null,-1,-1,-1,-1));
		
		outputParams.setName("LaminarSmoothedIntensity");
		outputParams.setLabel("Laminar Smoothed Intensity");
	}

	
	@Override
	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		
		// import the image data
		System.out.println("Import data");
		
		String name = Interface.getName(intensityImage);
		ImageHeader header = Interface.getHeader(intensityImage);
		
		int[] dims = Interface.getDimensions(intensityImage);
		float[] res = Interface.getResolutions(intensityImage);
		
		int nlayers = Interface.getComponents(layersImage)-1;
		int nt = 1;
		if (Interface.isImage4D(intensityImage)) nt = Interface.getDimensions4D(intensityImage)[3];
		
		// main algorithm
		algorithm = new LaminarIterativeSmoothing();
		
		algorithm.setProfileSurfaceImage(Interface.getFloatImage4D(layersImage));
		if (Interface.isImage4D(intensityImage)) {
		    algorithm.setIntensityImage(Interface.getFloatImage4D(intensityImage));
		} else {
            algorithm.setIntensityImage(Interface.getFloatImage3D(intensityImage));
		} 
		if (Interface.isValid(maskImage)) {
		    algorithm.setROIMask(Interface.getIntegerImage3D(maskImage));
		}
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);
		algorithm.setLayers(nlayers);
		algorithm.set4thDimension(nt);
		
		algorithm.setFWHMmm(fwhmParam.getValue().floatValue());

		algorithm.execute();
		
		// output
		if (Interface.isImage4D(intensityImage)) {
		    Interface.setFloatImage4D(algorithm.getSmoothedIntensityImage(), dims, dims[3], smoothIntensityImage, name+"_lis_smoothed", header);
		} else {
		    Interface.setFloatImage3D(algorithm.getSmoothedIntensityImage(), dims, smoothIntensityImage, name+"_lis_smoothed", header);
		}
	}
}
