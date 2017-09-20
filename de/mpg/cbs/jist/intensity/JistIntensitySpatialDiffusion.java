package de.mpg.cbs.jist.intensity;

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
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;
import de.mpg.cbs.core.intensity.IntensitySpatialDiffusion;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensitySpatialDiffusion extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;
	
	private ParamVolume var1Image;
	private ParamVolume var2Image;
	private ParamVolume var3Image;
	
	private ParamVolume 	probaImage;
	private ParamVolume 	maskImage;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 		imgscaleParam;
	private	ParamFloat 		certainscaleParam;
	private ParamFloat 		mincertaintyParam;
	private ParamFloat	 	scaleParam;
	private ParamInteger 	neighborParam;
	private ParamFloat	 	diffratioParam;
	
	private ParamVolume output1Image;
	private ParamVolume output2Image;
	private ParamVolume output3Image;
	private ParamVolume outprobaImage;
	
	private IntensitySpatialDiffusion algorithm;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Data");
		
		imageParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		imageParams.add(var1Image = new ParamVolume("Noise Image 1"));
		
		imageParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		imageParams.add(var2Image = new ParamVolume("Noise Image 2 (opt)"));
		var2Image.setMandatory(false);
		input2Image.setMandatory(false);
		
		imageParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		imageParams.add(var3Image = new ParamVolume("Noise Image 3 (opt)"));
		var3Image.setMandatory(false);
		input3Image.setMandatory(false);
		
		imageParams.add(probaImage = new ParamVolume("Signal Probability Map"));
		imageParams.add(maskImage = new ParamVolume("Image Mask (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(imgscaleParam = new ParamFloat("Image Scale", 0.0f, 100.0f, 0.5f));
		mainParams.add(certainscaleParam = new ParamFloat("Certainty Scale", 0.0f, 100.0f, 2.0f));
		mainParams.add(mincertaintyParam = new ParamFloat("Min Certainty", 0, 1, 0.5f));
		mainParams.add(neighborParam = new ParamInteger("Neighborhood size", 0, 26, 6));
		mainParams.add(diffratioParam = new ParamFloat("Min difference ratio", 0, 1, 0.01f));
		
		inputParams.add(mainParams);
		
		algorithm = new IntensitySpatialDiffusion();
		
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
		outputParams.add(output1Image = new ParamVolume("Diffused Image 1",VoxelType.FLOAT));
		outputParams.add(output2Image = new ParamVolume("Diffused Image 2",VoxelType.FLOAT));
		output2Image.setMandatory(false);
		outputParams.add(output3Image = new ParamVolume("Diffused Image 3",VoxelType.FLOAT));
		output3Image.setMandatory(false);
		outputParams.add(outprobaImage = new ParamVolume("Diffused Probability Image",VoxelType.FLOAT));
		
		outputParams.setName("sid images");
		outputParams.setLabel("sid images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// main algorithm
		algorithm = new IntensitySpatialDiffusion();
		
		// i/o variables
		String name1 = Interface.getName(input1Image);
		ImageHeader header1 = Interface.getHeader(input1Image);
		String namep = Interface.getName(probaImage);
		ImageHeader headerp = Interface.getHeader(probaImage);
		
		String name2 =null, name3 = null;
		ImageHeader header2 = null, header3 = null;
		if (Interface.isValid(input2Image)) { name2 = Interface.getName(input2Image); header2 = Interface.getHeader(input2Image); }
		if (Interface.isValid(input3Image)) { name3 = Interface.getName(input3Image); header3 = Interface.getHeader(input3Image); }
		
		int[] dims = Interface.getDimensions(input1Image);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(input1Image));
		
		// pass the input parameters
		algorithm.setContrastImage1(Interface.getFloatImage3D(input1Image));
		if (Interface.isValid(input2Image)) algorithm.setContrastImage2(Interface.getFloatImage3D(input2Image));
		if (Interface.isValid(input3Image)) algorithm.setContrastImage3(Interface.getFloatImage3D(input3Image));
		
		algorithm.setNoiseImage1(Interface.getFloatImage3D(var1Image));
		if (Interface.isValid(var2Image)) algorithm.setNoiseImage2(Interface.getFloatImage3D(var2Image));
		if (Interface.isValid(var3Image)) algorithm.setNoiseImage3(Interface.getFloatImage3D(var3Image));
		
		algorithm.setSignalProbabilityImage(Interface.getFloatImage3D(probaImage));
		
		algorithm.setMaxIterations(iterationParam.getValue().intValue());
		algorithm.setImageScale(imgscaleParam.getValue().floatValue());
		algorithm.setCertaintyScale(certainscaleParam.getValue().floatValue());
		algorithm.setMinCertainty(mincertaintyParam.getValue().floatValue());
		algorithm.setNeighborhoodSize(neighborParam.getValue().intValue());
		algorithm.setMinDifferenceRatio(diffratioParam.getValue().floatValue());
		
		algorithm.execute();
		
		// outputs
		Interface.setFloatImage3D(algorithm.getDiffusedImage1(), dims, output1Image, name1+"_df", header1);
		if (Interface.isValid(input2Image)) Interface.setFloatImage3D(algorithm.getDiffusedImage2(), dims, output2Image, name2+"_df", header2);
		if (Interface.isValid(input3Image)) Interface.setFloatImage3D(algorithm.getDiffusedImage3(), dims, output3Image, name3+"_df", header3);
		
		Interface.setFloatImage3D(algorithm.getDiffusedProbabilityImage(), dims, outprobaImage, namep+"_df", headerp);
		
		return;
	}
}
