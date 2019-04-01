package de.mpg.cbs.jist.filter;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataDouble;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.core.filter.FilterRecursiveRidgeDiffusion;

import java.io.*;

public class JistFilterRecursiveRidgeDiffusion extends ProcessingAlgorithm {
	
	// jist containers
	private ParamVolume inputImage;
	private ParamOption brightParam;
	
	private ParamVolume surfaceImage;
	private ParamOption	orientationParam;
	private ParamFloat	angleParam;
	
	private ParamVolume locationImage;
	
	private ParamInteger minscaleParam;
	private ParamInteger maxscaleParam;
	
	//private ParamFloat thresholdParam;
	
	private ParamOption filterParam;
	
	private ParamOption propagationParam;
	private ParamFloat difffactorParam;
	private ParamFloat simscaleParam;
	private ParamInteger ngbParam;
	private ParamInteger iterParam;
	private ParamFloat maxdiffParam;
	
	private ParamVolume pvImage;
	private ParamVolume filterImage;
	private ParamVolume probaImage;
	private ParamVolume propagImage;
	private ParamVolume scaleImage;
	private ParamVolume directionImage;
	private ParamVolume correctImage;
	private ParamVolume sizeImage;
	
	private FilterRecursiveRidgeDiffusion algorithm;
			
	@Override
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(brightParam = new ParamOption("Ridge intensities", FilterRecursiveRidgeDiffusion.brightTypes));
		inputParams.add(filterParam = new ParamOption("Ridge filter", FilterRecursiveRidgeDiffusion.filterTypes));
		inputParams.add(surfaceImage = new ParamVolume("Surface(s) level sets"));
		surfaceImage.setMandatory(false);
		inputParams.add(orientationParam = new ParamOption("Orientation relative to surface", FilterRecursiveRidgeDiffusion.orientationTypes));
		inputParams.add(angleParam = new ParamFloat("Angular factor", 0.1f, 10.0f, 1.0f));
		
		inputParams.add(locationImage = new ParamVolume("Location prior (opt)"));
		locationImage.setMandatory(false);
		
		//inputParams.add(thresholdParam = new ParamFloat("Probability threshold", 0.0, 1.0, 0.5));
		//inputParams.add(angleParam = new ParamFloat("Scale factor ", 0.25f, 2.0f, 1.0f));
		inputParams.add(minscaleParam = new ParamInteger("Minimum scale ", 0, 100, 0));
		inputParams.add(maxscaleParam = new ParamInteger("Maximum scale ", 0, 100, 3));

		inputParams.add(propagationParam = new ParamOption("Propagation model", FilterRecursiveRidgeDiffusion.propagationTypes));
		inputParams.add(difffactorParam = new ParamFloat("Diffusion factor", 0.0f, 100.0f, 1.0f));
		inputParams.add(simscaleParam = new ParamFloat("Similarity scale", 0.0f, 100.0f, 0.1f));
		inputParams.add(ngbParam = new ParamInteger("Neighborhood size", 0, 26, 4));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000, 100));
		inputParams.add(maxdiffParam = new ParamFloat("Max difference", 0.0f, 1.0f, 0.001f));
		
		algorithm = new FilterRecursiveRidgeDiffusion();
		
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
		outputParams.add(pvImage = new ParamVolume("Ridge partial volume",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(filterImage = new ParamVolume("Filter response",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(probaImage = new ParamVolume("Probability response",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(propagImage = new ParamVolume("Propagated response",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(scaleImage = new ParamVolume("Detection scale",VoxelType.BYTE,-1,-1,-1,-1));
		outputParams.add(directionImage = new ParamVolume("Ridge direction",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(correctImage = new ParamVolume("Directional correction",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(sizeImage = new ParamVolume("Ridge size",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("Recursive Ridge Diffusion");
		outputParams.setLabel("Recursive Ridge Diffusion");
	}
	
	@Override
	protected void execute(CalculationMonitor monitor){
	    
		// i/o variables
		String name = Interface.getName(inputImage);
		ImageHeader header = Interface.getHeader(inputImage);
		int[] dims = Interface.getDimensions(inputImage);
		float[] res = Interface.getResolutions(inputImage);
		
		// main algorithm
		algorithm = new FilterRecursiveRidgeDiffusion();
	    
		algorithm.setInputImage(Interface.getFloatImage3D(inputImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setRidgeIntensities(brightParam.getValue());
		algorithm.setRidgeFilter(filterParam.getValue());
		if (Interface.isValid(surfaceImage)) algorithm.setSurfaceLevelSet(Interface.getFloatImage3D(surfaceImage));
		algorithm.setOrientationToSurface(orientationParam.getValue());
		algorithm.setAngularFactor(angleParam.getValue().floatValue());
		
		if (Interface.isValid(locationImage)) algorithm.setLocationPrior(Interface.getFloatImage3D(locationImage));
		
		algorithm.setMinimumScale(minscaleParam.getValue().intValue());
		algorithm.setMaximumScale(maxscaleParam.getValue().intValue());
		
		algorithm.setPropagationModel(propagationParam.getValue());
		algorithm.setDiffusionFactor(difffactorParam.getValue().floatValue());
		algorithm.setSimilarityScale(simscaleParam.getValue().floatValue());
		algorithm.setNeighborhoodSize(ngbParam.getValue().intValue());
		algorithm.setMaxIterations(iterParam.getValue().intValue());
		algorithm.setMaxDifference(maxdiffParam.getValue().floatValue());
	
		algorithm.execute();
		
		Interface.setFloatImage3D(algorithm.getRidgePartialVolumeImage(), dims, pvImage, name+"_rrf_pv", header);
		Interface.setFloatImage3D(algorithm.getFilterResponseImage(), dims, filterImage, name+"_rrf_fr", header);
		Interface.setFloatImage3D(algorithm.getProbabilityResponseImage(), dims, probaImage, name+"_rrf_proba", header);
		Interface.setFloatImage3D(algorithm.getPropagatedResponseImage(), dims, propagImage, name+"_rrf_propag", header);
		Interface.setIntegerImage3D(algorithm.getDetectionScaleImage(), dims, scaleImage, name+"_rrf_scale", header);
		Interface.setFloatImage4D(algorithm.getRidgeDirectionImage(), dims, 3, directionImage, name+"_rrf_dir", header);
		Interface.setFloatImage3D(algorithm.getDirectionalCorrectionImage(), dims, correctImage, name+"_rrf_corr", header);
		Interface.setFloatImage3D(algorithm.getRidgeSizeImage(), dims, sizeImage, name+"_rrf_size", header);
		
	}
	
}