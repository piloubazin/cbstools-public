package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
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

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.laminar.LaminarVolumetricLayering;

import javax.vecmath.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarVolumetricLayering extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume gwImage;
	private ParamVolume cgImage;
	//private ParamVolume segImage;
	
	private ParamInteger 	layersParam;
	private ParamInteger 	iterationParamNarrowBand;
	private ParamFloat 	minimumParamNarrowBand;
	private ParamOption 	algoParam;
	private ParamOption 	dirParam;
	private ParamOption 	topologyParam;
	private ParamInteger	kernelParam;
	private ParamFloat 	ratioKernelParam;
	//private ParamFloat 	curvscaleParam;
	private	ParamBoolean	presmoothParam;
	
	private static final String[] algoTypes = {"distance-preserving", "volume-preserving"};
	private static final String[] dirTypes = {"outward", "inward"};
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume layeringImage;
	private ParamVolume labelsImage;
	private ParamVolume surfImage;
	private ParamVolume midsurfImage;
	private ParamVolume debugImage;
	
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	
	private static final boolean debugOutput=false;
	
	private LaminarVolumetricLayering algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(gwImage = new ParamVolume("Inner Distance Image (GM/WM boundary)"));
		imageParams.add(cgImage = new ParamVolume("Outer Distance Image (CSF/GM boundary)"));
		//imageParams.add(segImage = new ParamVolume("Cortex and WM Mask (opt, GM=1, WM=2)"));
		//segImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		mainParams.add(layersParam = new ParamInteger("Number of layers", 1, 50, 10)); // min=2, max=100, default=2
		mainParams.add(iterationParamNarrowBand = new ParamInteger("Max iterations for narrow band evolution", 0, 100000, 500)); 
		mainParams.add(minimumParamNarrowBand = new ParamFloat("Min change ratio for narrow band evolution", 0, 0.05f, 0.0005f)); 
		
		mainParams.add(algoParam = new ParamOption("Layering method", algoTypes));
		algoParam.setValue("volume-preserving");

		mainParams.add(dirParam = new ParamOption("Layering direction", dirTypes));
		dirParam.setValue("outward");

		mainParams.add(kernelParam=new ParamInteger("curvature approximation scale (voxels)", 2, 15, 3));
		mainParams.add(ratioKernelParam=new ParamFloat("ratio smoothing kernel size (voxels)", 0.0f, 6.0f, 1.0f));
		mainParams.add(presmoothParam=new ParamBoolean("pre-smooth cortical surfaces", false));
		
		//mainParams.add(iterationParamBisekt = new ParamInteger("Max. iterations for bisection method", 0, 100000, 500));
		//mainParams.add(minDiffParamBisekt = new ParamFloat("epsilon of bisection method", 0.0f, 1e-3f, 1e-6f));
		
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("no");

		//mainParams.add(curvscaleParam=new ParamFloat("Curvature scale", 0.0f, 2.0f, 1.0f));
		
		inputParams.add(mainParams);
		
		algorithm = new LaminarVolumetricLayering();
		
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
		outputParams.add(layeringImage = new ParamVolume("Continuous depth measurement",VoxelType.FLOAT));
		outputParams.add(labelsImage = new ParamVolume("Discrete sampled layers",VoxelType.UBYTE));
		outputParams.add(surfImage = new ParamVolume("Layer boundary surfaces",null,-1,-1,-1,-1));
		outputParams.add(midsurfImage = new ParamVolume("Layer-centered surfaces",null,-1,-1,-1,-1));
		if (debugOutput) outputParams.add(debugImage = new ParamVolume("Debug data",null,-1,-1,-1,-1));
		
		outputParams.setName("layering images");
		outputParams.setLabel("layering images");
		
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		String name = Interface.commonNameBase(Interface.getName(gwImage),Interface.getName(cgImage));
		ImageHeader header = Interface.getHeader(gwImage);
		int[] dims = Interface.getDimensions(gwImage);
		float[] res = Interface.getResolutions(gwImage);
		
		// main algorithm
		algorithm = new LaminarVolumetricLayering();
		
		algorithm.setInnerDistanceImage(Interface.getFloatImage3D(gwImage));
		algorithm.setOuterDistanceImage(Interface.getFloatImage3D(cgImage));
		
		algorithm.setNumberOfLayers(layersParam.getValue().intValue());
		algorithm.setMaxNarrowBandIterations(iterationParamNarrowBand.getValue().intValue()); 
		algorithm.setMinNarrowBandChange(minimumParamNarrowBand.getValue().intValue()); 
		algorithm.setLayeringMethod(algoParam.getValue());
		algorithm.setLayeringDirection(dirParam.getValue());
		algorithm.setCurvatureApproximationScale(kernelParam.getValue().intValue());
		algorithm.setRatioSmoothingKernelSize(ratioKernelParam.getValue().floatValue());
		algorithm.setPresmoothCorticalSurfaces(presmoothParam.getValue().booleanValue());
		algorithm.setTopology(topologyParam.getValue());
			
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.execute();

		// output
		BasicInfo.displayMessage("...depth\n");
		
		Interface.setFloatImage3D(algorithm.getContinuousDepthMeasurement(), dims, layeringImage, name+"_depth", header);
		
		BasicInfo.displayMessage("...labels\n");
		
		Interface.setUByteImage3D(algorithm.getDiscreteSampledLayers(), dims, labelsImage, name+"_labels", header);
		
		BasicInfo.displayMessage("...surfaces\n");
		Interface.setFloatImage4D(algorithm.getLayerBoundarySurfaces(), dims, algorithm.getLayerBoundarySurfacesLength(), surfImage, name+"_layer_boundaries", header);
	
		Interface.setFloatImage4D(algorithm.getLayerCenteredSurfaces(), dims, algorithm.getLayerCenteredSurfacesLength(), midsurfImage, name+"_layer_centers", header);
	
		BasicInfo.displayMessage("done\n");		
	}


}
