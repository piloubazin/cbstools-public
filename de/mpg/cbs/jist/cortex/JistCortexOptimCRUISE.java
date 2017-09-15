package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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

import de.mpg.cbs.core.cortex.CortexOptimCRUISE;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistCortexOptimCRUISE extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume initImage;
	private ParamVolume wmImage;
	private ParamVolume gmImage;
	private ParamVolume csfImage;
	private ParamVolume vdImage;
	
	//private ParamFloat 	edgeParam;
	private	ParamFloat 	balloonParam;
	private ParamFloat 	curvParam;
	private ParamInteger 	iterationParam;
	//private ParamOption 	gwImageParam;
	//private ParamOption 	cgImageParam;
	private ParamOption 	topologyParam;
	private ParamBoolean	normalizeParam;
	private ParamBoolean	pvwmParam;
	//private ParamFloat 	offsetParam;
	private ParamFloat 		wmdistParam;
	
	//private static final String[] opts = {"Iso image", "T1 map", "Probabilities", "Iso+Proba", "T1map+Proba", "All"};
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume cortexImage;
	private ParamVolume gwbImage;
	private ParamVolume cgbImage;
	private ParamVolume pgmImage;
	private ParamVolume avgImage;
	private ParamVolume pwmImage;
	private ParamVolume thickImage;
	private ParamVolume pcsfImage;
	
	private CortexOptimCRUISE algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(initImage = new ParamVolume("Initial WM Segmentation Image"));
		imageParams.add(wmImage = new ParamVolume("Filled WM Probability Image"));
		imageParams.add(gmImage = new ParamVolume("GM Probability Image"));
		imageParams.add(csfImage = new ParamVolume("CSF and BG Probability Image"));
		imageParams.add(vdImage = new ParamVolume("Veins and Dura Probability Image (opt)"));
		vdImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		//mainParams.add(edgeParam = new ParamFloat("Edge weight", -1E10f, 1E10f, 0.1f25));
		mainParams.add(balloonParam = new ParamFloat("Data weight (balloon force)", -1E10f, 1E10f, 0.9f));
		mainParams.add(curvParam = new ParamFloat("Regularization weight (curvature)", -1E10f, 1E10f, 0.1f));
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
			
		//mainParams.add(gwImageParam = new ParamOption("Image used for WM / GM", opts));
		//gwImageParam.setValue("Iso+Proba");

		//mainParams.add(cgImageParam = new ParamOption("Image used for GM / CSF", opts));
		//cgImageParam.setValue("Iso+Proba");

		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(normalizeParam = new ParamBoolean("Normalize probabilities", true));
		mainParams.add(pvwmParam = new ParamBoolean("Correct for WM-GM partial voluming", true));
		//mainParams.add(offsetParam = new ParamFloat("WM/GM offset", -1, 1, 0.0f));
		mainParams.add(wmdistParam = new ParamFloat("WM drop-off distance", 0.1f, 100.0f, 1.0f));
		
		inputParams.add(mainParams);
			
		algorithm = new CortexOptimCRUISE();
		
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
		outputParams.add(cortexImage = new ParamVolume("Cortex Mask",VoxelType.UBYTE));
		outputParams.add(gwbImage = new ParamVolume("WM-GM Levelset",VoxelType.FLOAT));
		outputParams.add(cgbImage = new ParamVolume("GM-CSF Levelset",VoxelType.FLOAT));
		outputParams.add(avgImage = new ParamVolume("Central Levelset",VoxelType.FLOAT));
		outputParams.add(thickImage = new ParamVolume("Cortical thickness",VoxelType.FLOAT));
		outputParams.add(pwmImage = new ParamVolume("Cerebral WM probability",VoxelType.FLOAT));
		outputParams.add(pgmImage = new ParamVolume("Cortical GM probability",VoxelType.FLOAT));
		outputParams.add(pcsfImage = new ParamVolume("Sulcal CSF probability",VoxelType.FLOAT));
		
		outputParams.setName("cortical images");
		outputParams.setLabel("cortical images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// i/o variables
		// import the image data into 1D arrays
		String name = Interface.getName(gmImage);
		ImageHeader header = Interface.getHeader(gmImage);
		int[] dims = Interface.getDimensions(gmImage);
		float[] res = Interface.getResolutions(gmImage);
		
		// main algorithm
		algorithm = new CortexOptimCRUISE();
		
		algorithm.setInitialWMSegmentationImage(Interface.getUByteImage3D(initImage));
		algorithm.setFilledWMProbabilityImage(Interface.getFloatImage3D(wmImage));
		algorithm.setGMProbabilityImage(Interface.getFloatImage3D(gmImage));
		algorithm.setCSFandBGProbabilityImage(Interface.getFloatImage3D(csfImage));
		if (Interface.isValid(vdImage)) algorithm.setVeinsAndDuraProbabilityImage(Interface.getFloatImage3D(vdImage));
		
		algorithm.setDataWeight(balloonParam.getValue().floatValue());
		algorithm.setRegularizationWeight(curvParam.getValue().floatValue()); 
		algorithm.setMaxIterations(iterationParam.getValue().intValue()); 
		algorithm.setNormalizeProbabilities(normalizeParam.getValue().booleanValue());
		algorithm.setCorrectForWMGMpartialVoluming(pvwmParam.getValue().booleanValue());
		algorithm.setWMdropoffDistance(wmdistParam.getValue().floatValue());
		algorithm.setTopology(topologyParam.getValue());
			
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.execute();

		// output
		Interface.setUByteImage3D(algorithm.getCortexMask(), dims, cortexImage, name+"_cruise_cortex", header);
		Interface.setFloatImage3D(algorithm.getWMGMLevelset(), dims, gwbImage, name+"_cruise_gwb", header);
		Interface.setFloatImage3D(algorithm.getGMCSFLevelset(), dims, cgbImage, name+"_cruise_cgb", header);
		Interface.setFloatImage3D(algorithm.getCentralLevelset(), dims, avgImage, name+"_cruise_avg", header);
		Interface.setFloatImage3D(algorithm.getCorticalThickness(), dims, thickImage, name+"_cruise_thick", header);
		Interface.setFloatImage3D(algorithm.getCerebralWMprobability(), dims, pwmImage, name+"_cruise_pwm", header);
		Interface.setFloatImage3D(algorithm.getCorticalGMprobability(), dims, pgmImage, name+"_cruise_pgm", header);
		Interface.setFloatImage3D(algorithm.getSulcalCSFprobability(), dims, pcsfImage, name+"_cruise_pcsf", header);
	}
}
