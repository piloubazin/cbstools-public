package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityRangeNormalization extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	normtypes = {"linear","robust","robust-min","robust-max"};
	private		String		normtype = "robust";
	
	// jist containers
	private ParamVolume inImage;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private	ParamOption	 normParam;
	private ParamFloat	 ratioParam;
	private ParamBoolean ignoreNegParam;
	private ParamBoolean ignoreZeroParam;
	private ParamFloat	 scalingParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("Mask Image"));
		maskImage.setMandatory(false);
		inputParams.add(normParam = new ParamOption("normalization type", normtypes));
		normParam.setDescription("linear: use the whole range [Imin,Imax]; robust: use a truncated version of [Imin,Imax]; robust-max: use a truncated version of [0,Imax]");
		normParam.setValue(normtype);
		inputParams.add(ratioParam = new ParamFloat("Robustness ratio", 0, 1, 0.01f));
		ratioParam.setDescription("ratio of discarded values below Imin and above Imax (in [0,1])");
		inputParams.add(ignoreNegParam = new ParamBoolean("set negative values to zero", true));
		inputParams.add(ignoreZeroParam = new ParamBoolean("ignore zero values", true));
		inputParams.add(scalingParam = new ParamFloat("Output scaling", 0, 1e10f, 1.0f));
		scalingParam.setDescription("scaling the output image into [0,S]");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("Intensity Range Normalization");
		inputParams.setName("IntensityRangeNormalization");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Outputs a normalized version of the image in [0,1].");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Normalized Image",VoxelType.FLOAT));
		outputParams.setName("normalize_image");
		outputParams.setLabel("Normalize Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat input = new ImageDataFloat(inImage.getImageData());
		float[][][] image = input.toArray3d();
		int nx = input.getRows();
		int ny= input.getCols();
		int nz = input.getSlices();
		float rx = input.getHeader().getDimResolutions()[0];
		float ry = input.getHeader().getDimResolutions()[1];
		float rz = input.getHeader().getDimResolutions()[2];
		
		// use a mask by default
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = true;
		}
		// input mask
		if (maskImage.getImageData()!=null) {
			ImageDataUByte maskI = new ImageDataUByte(maskImage.getImageData());
			byte[][][] tmp = maskI.toArray3d();
			mask = new boolean[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				mask[x][y][z] = (tmp[x][y][z]!=0);
			}
		}
		// negative to zero
		if (ignoreNegParam.getValue().booleanValue()) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (image[x][y][z]<0) image[x][y][z] = 0.0f;
			}
		}
		// ignore zero values: masked
		if (ignoreZeroParam.getValue().booleanValue()) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (image[x][y][z]==0) mask[x][y][z] = false;
			}
		}
		
		
		// main algorithm
		
		// 1. estimate the min, max
		float Imin = 0, Imax = 0;
		if (normParam.getValue().equals("linear")) {
			Imin = ImageStatistics.minimum(image, mask, nx, ny, nz); 
			Imax = ImageStatistics.maximum(image, mask, nx, ny, nz);				
		} else if (normParam.getValue().equals("robust")) {
			Imin = ImageStatistics.robustMinimum(image, mask, ratioParam.getValue().floatValue(), 4, nx, ny, nz); 
			Imax = ImageStatistics.robustMaximum(image, mask, ratioParam.getValue().floatValue(), 4, nx, ny, nz);
		} else if (normParam.getValue().equals("robust-min")) {
			Imin = ImageStatistics.robustMinimum(image, mask, ratioParam.getValue().floatValue(), 4, nx, ny, nz); 
			Imax = ImageStatistics.maximum(image, nx, ny, nz);
		} else if (normParam.getValue().equals("robust-max")) {
			Imin = 0.0f;
			Imax = ImageStatistics.robustMaximum(image, mask, ratioParam.getValue().floatValue(), 4, nx, ny, nz);
		}
		BasicInfo.displayMessage("image min, max: "+Imin+", "+Imax+"\n");
			
		// 2. scale the data
		float scaling = scalingParam.getValue().floatValue();
		float[][][] result = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = scaling*Numerics.bounded( (image[x][y][z]-Imin)/(Imax-Imin), 0.0f, 1.0f);
		}

		ImageDataFloat resultData = new ImageDataFloat(result);		
		resultData.setHeader(input.getHeader());
		resultData.setName(input.getName()+"_norm");
		resultImage.setValue(resultData);
		resultData = null;
		result = null;
		
	}
	
}
