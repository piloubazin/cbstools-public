package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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
	private ParamCollection	imgParam;
	private ParamCollection	mrParam;
	private ParamCollection scaleParam;
	
	private ParamVolume 	inv1Image;
	private ParamVolume 	inv2Image;
	private ParamVolume 	b1Image;
	private ParamFloat 		invertRepTimeParam;
	private ParamFloat 		excite1RepTimeParam;
	private ParamFloat 		excite2RepTimeParam;
	private ParamInteger 	exciteNumberParam;
	private ParamFloat 		inversionTime1Param;
	private ParamFloat 		inversionTime2Param;
	private ParamFloat 		flipAngle1Param;
	private ParamFloat 		flipAngle2Param;
	private ParamFloat 		invertEffParam;
	//private ParamOption methodParam;
	private ParamBoolean	useB1Param;
	private ParamFloat		b1ScalingParam;
	
	private ParamVolume uniformImage;
	private ParamVolume t1mapImage;
	private ParamVolume r1mapImage;
	private ParamVolume relativeSnrImage;
	private ParamVolume lutImage;
	private ParamVolume invlutImage;
	
	// parameters
	//private		static final String[]	methods = IntensityMp2rageT1Fitting.methods;
	//private		String		method = "normal";
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imgParam = new ParamCollection("images");
		imgParam.add(inv1Image = new ParamVolume("First Inversion Image (4D:mag+phs)"));
		imgParam.add(inv2Image = new ParamVolume("Second Inversion Image (4D:mag+phs)"));
		imgParam.add(b1Image = new ParamVolume("B1 map (calculated & co-registered)"));
		b1Image.setMandatory(false);
		inputParams.add(imgParam);
		
		mrParam = new ParamCollection("MR parameters");
		mrParam.add(inversionTime1Param = new ParamFloat("First inversion time (sec)", 0.0f, 10000.0f, 0.8f));
		mrParam.add(inversionTime2Param = new ParamFloat("Second inversion time (sec)", 0.0f, 10000.0f, 2.7f));
		mrParam.add(flipAngle1Param = new ParamFloat("First flip angle (deg)", 0.0f, 180.0f, 7.0f));
		mrParam.add(flipAngle2Param = new ParamFloat("Second flip angle (deg)", 0.0f, 180.0f, 5.0f));
		mrParam.add(invertEffParam = new ParamFloat("Inversion efficiency", 0.0f, 1.0f, 0.96f));
		
		mrParam.add(invertRepTimeParam = new ParamFloat("Inversion repetition time (sec)", 0.0f, 10000.0f, 5.5f));
		mrParam.add(excite1RepTimeParam = new ParamFloat("First Excitation repetition time (sec)", 0.0f, 10000.0f, 0.0062f));
		mrParam.add(excite2RepTimeParam = new ParamFloat("Second Excitation repetition time (sec)", 0.0f, 10000.0f, 0.0062f));
		mrParam.add(exciteNumberParam = new ParamInteger("Number of excitations", 0, 10000, 160));
		
		mrParam.add(useB1Param = new ParamBoolean("correct B1 inhomogeneities", true));
		mrParam.add(b1ScalingParam = new ParamFloat("B1 map scaling", 0.0f, 10000000000.0f, 1000000.0f));
		inputParams.add(mrParam);
		
		algorithm = new IntensityMp2rageT1Fitting();
		
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
		outputParams.add(r1mapImage = new ParamVolume("Quantitative R1 map Image",VoxelType.FLOAT));
		outputParams.add(relativeSnrImage = new ParamVolume("SNR Ratio Image",VoxelType.FLOAT));
		
		//outputParams.add(lutImage = new ParamVolume("T1 LUT Image",VoxelType.FLOAT));
		//outputParams.add(invlutImage = new ParamVolume("UNI LUT Image",VoxelType.FLOAT));
		
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
		if (Interface.isValid(b1Image)) algorithm.setB1mapImage(Interface.getFloatImage3D(b1Image));
			
		algorithm.setDimensions(dims);
		algorithm.setResolutions(res);

		algorithm.setFirstInversionTime(inversionTime1Param.getValue().floatValue());
		algorithm.setSecondInversionTime(inversionTime2Param.getValue().floatValue());
		algorithm.setFirstFlipAngle(flipAngle1Param.getValue().floatValue());
		algorithm.setSecondFlipAngle(flipAngle2Param.getValue().floatValue());
		algorithm.setInversionEfficiency(invertEffParam.getValue().floatValue());
		
		algorithm.setInversionRepetitionTime(invertRepTimeParam.getValue().floatValue());
		algorithm.setFirstExcitationRepetitionTime(excite1RepTimeParam.getValue().floatValue());
		algorithm.setSecondExcitationRepetitionTime(excite2RepTimeParam.getValue().floatValue());
		algorithm.setNumberExcitations(exciteNumberParam.getValue().intValue());
		
		algorithm.setCorrectB1inhomogeneities(useB1Param.getValue().booleanValue());
		algorithm.setB1mapScaling(b1ScalingParam.getValue().floatValue());
		
		algorithm.execute();

		Interface.setFloatImage3D(algorithm.getUniformT1weightedImage(), dims, uniformImage, name+"_uni", header);
		Interface.setFloatImage3D(algorithm.getQuantitativeT1mapImage(), dims, t1mapImage, name+"_qt1", header);
		Interface.setFloatImage3D(algorithm.getQuantitativeR1mapImage(), dims, r1mapImage, name+"_qr1", header);
		Interface.setFloatImage3D(algorithm.getRelativeSnrImage(), dims, relativeSnrImage, name+"_snr", header);
		
		//Interface.setFloatImage2D(algorithm.generateT1LookupImage(), new int[]{200,200}, lutImage, name+"_lut", header);
		//Interface.setFloatImage2D(algorithm.generateUniLookupImage(), new int[]{200,200}, invlutImage, name+"_invlut", header);
		
	}


}
