package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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
public class JistIntensityMp2rageT1Fitting extends ProcessingAlgorithm {

	private IntensityMp2rageT1Fitting algorithm;
	
	// jist containers
	private ParamVolume inv1Image;
	private ParamVolume inv2Image;
	private ParamFloat repetitionTimeParam;
	private ParamFloat inversionTime1Param;
	private ParamFloat inversionTime2Param;
	private ParamFloat flipAngle1Param;
	private ParamFloat flipAngle2Param;
	//private ParamOption methodParam;
	
	private ParamVolume uniformImage;
	private ParamVolume t1mapImage;
	private ParamVolume relativeSnrImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inv1Image = new ParamVolume("First Inversion Image"));
		inputParams.add(inv2Image = new ParamVolume("Second Inversion Image"));
		
		inputParams.add(repetitionTimeParam = new ParamFloat("Repetition time (TR)", 0.0f, 10000.0f, 1.0f));
		inputParams.add(inversionTime1Param = new ParamFloat("First inversion time (TI1)", 0.0f, 10000.0f, 1.0f));
		inputParams.add(inversionTime2Param = new ParamFloat("Second inversion time (TI2)", 0.0f, 10000.0f, 1.0f));
		inputParams.add(flipAngle1Param = new ParamFloat("First flip angle (deg)", 0.0f, 180.0f, 1.0f));
		inputParams.add(flipAngle2Param = new ParamFloat("Second flip angle (deg)", 0.0f, 180.0f, 1.0f));
		
		algorithm = new IntensityMp2ragerT1Fitting();
		
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
		outputParams.add(uniformImage = new ParamVolume("Uniform T1-weighted Image",VoxelType.FLOAT));
		outputParams.add(t1mapImage = new ParamVolume("Quantitative T1 map Image",VoxelType.FLOAT));
		outputParams.add(relativeSnrImage = new ParamVolume("SNR Ratio Image",VoxelType.FLOAT));
		
		outputParams.setName("reconstructed images");
		outputParams.setLabel("reconstructed images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		// i/o variables
		int[] dims = Interface.getDimensions(inv1Image);
		float[] res = Interface.getResolutions(inv1Image);
		String name = Interface.getName(inv1Image);
		ImageHeader header = Interface.getHeader(inv1Image);
		
		// main algorithm
		algorithm = new IntensityMp2rageT1Fitting();
		
		algorithm.setFirstInversionImage(Interface.getFloatImage4D(inv1Image));
		algorithm.setSecondInversionImage(Interface.getFloatImage4D(inv2Image));
		if (Interface.isValid(maskImage)) algorithm.setInputMask(Interface.getUByteImage3D(maskImage));
		
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setRepetitionTime(repetitionTime1Param.getValue().floatValue());
		algorithm.setFirstInversionTime(inversionTime1Param.getValue().floatValue());
		algorithm.setSecondInversionTime(inversionTime2Param.getValue().floatValue());
		algorithm.setFirstFlipAngle(flipAngle1Param.getValue().floatValue());
		algorithm.setSecondFlipAngle(flipAngle2Param.getValue().floatValue());
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getUniformT1weightedImage(), dims, uniformImage, name+"_uni", header);
		Interface.setFloatImage3D(algorithm.getQuantitativeT1mapImage(), dims, t1mapImage, name+"_qt1", header);
		Interface.setFloatImage3D(algorithm.getRelativeSnrImage(), dims, relativeSnrImage, name+"_snr", header);
	}


}
