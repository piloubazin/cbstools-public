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
public class JistIntensityPropagate extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	normtypes = {"max","mean","min"};
	private		String		normtype = "max";
	
	// jist containers
	private ParamVolume inImage;
	private ParamVolume maskImage;
	private ParamVolume resultImage;
	private	ParamOption	 normParam;
	private ParamFloat distParam;
	//private ParamBoolean ignoreNegParam;
	//private ParamBoolean ignoreZeroParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("Mask Image"));
		maskImage.setMandatory(false);
		inputParams.add(normParam = new ParamOption("combination method", normtypes));
		normParam.setDescription("use the mean, max or min data from neighboring voxels");
		normParam.setValue(normtype);
		inputParams.add(distParam = new ParamFloat("Propagation distance (mm)", 0, 50, 5));
		distParam.setDescription("distance for the propagation (note: this algorithm will be slow for large distances)");
		//inputParams.add(ignoreNegParam = new ParamBoolean("set negative values to zero", true));
		//inputParams.add(ignoreZeroParam = new ParamBoolean("ignore zero values", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("Intensity Propagation");
		inputParams.setName("IntensityPropagation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Propagates the values inside the mask (or non-zero) into the neighboring voxels");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Propagated Image",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.setName("propagate_image");
		outputParams.setLabel("Propagate Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat input = new ImageDataFloat(inImage.getImageData());
		float[][][] image3 = null;
		float[][][][] image4 = null;
		int nx = input.getRows();
		int ny= input.getCols();
		int nz = input.getSlices();
		int nc = input.getComponents();
		if (nc==1) image3 = input.toArray3d();
		else image4 = input.toArray4d();
		
		float rx = input.getHeader().getDimResolutions()[0];
		float ry = input.getHeader().getDimResolutions()[1];
		float rz = input.getHeader().getDimResolutions()[2];
		
		// use a mask by default
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			mask[x][y][z] = true;
		}
		// input mask or mask zero values
		if (maskImage.getImageData()!=null) {
			System.out.println("input mask");
			ImageDataUByte maskI = new ImageDataUByte(maskImage.getImageData());
			byte[][][] tmp = maskI.toArray3d();
			mask = new boolean[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				mask[x][y][z] = (tmp[x][y][z]!=0);
			}
		} else {
			System.out.println("mask from zero values");
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				if (nc==1 && image3[x][y][z]==0) mask[x][y][z] = false;
				if (nc>1) {
					mask[x][y][z] = false;
					for (int c=0;c<nc;c++) if (image4[x][y][z][c]!=0) mask[x][y][z] = true;
				}
			}
		}
		
		// main algorithm
		int nd = Numerics.ceil(distParam.getValue().floatValue()/Numerics.min(rx,ry,rz));
		
		byte MIN = 1, MAX = 2, MEAN = 3;
		byte merge = MAX;
		if (normParam.getValue().equals("mean")) merge = MEAN;
		if (normParam.getValue().equals("min")) merge = MIN;
		
		float[][][] result3 = null;
		float[][][][] result4 = null;
		if (nc==1) result3 = new float[nx][ny][nz];
		else result4 = new float[nx][ny][nz][nc];
		byte[][][] count = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (nc==1) result3[x][y][z] = image3[x][y][z];
			else for (int c=0;c<nc;c++) result4[x][y][z][c] = image4[x][y][z][c];
		}
		for (int n=0;n<nd;n++) {
			System.out.println("step "+(n+1));
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				count[x][y][z] = 0;
			}
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			
				for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) if (!mask[x+dx][y+dy][z+dz]) {
					//System.out.print(".");
					if (nc==1) {
						if (merge==MIN) {
							if (count[x+dx][y+dy][z+dz]==0) result3[x+dx][y+dy][z+dz] = result3[x][y][z];
							else result3[x+dx][y+dy][z+dz] = Numerics.min(result3[x+dx][y+dy][z+dz], result3[x][y][z]);
						} else if (merge==MAX) {
							if (count[x+dx][y+dy][z+dz]==0) result3[x+dx][y+dy][z+dz] = result3[x][y][z];
							else result3[x+dx][y+dy][z+dz] = Numerics.max(result3[x+dx][y+dy][z+dz], result3[x][y][z]);
						} else if (merge==MEAN) {
							result3[x+dx][y+dy][z+dz] += result3[x][y][z];
						}
					} else {
						for (int c=0;c<nc;c++) {
							if (merge==MIN) {
								if (count[x+dx][y+dy][z+dz]==0) result4[x+dx][y+dy][z+dz][c] = result4[x][y][z][c];
								else result4[x+dx][y+dy][z+dz][c] = Numerics.min(result4[x+dx][y+dy][z+dz][c], result4[x][y][z][c]);
							} else if (merge==MAX) {
								if (count[x+dx][y+dy][z+dz]==0) result4[x+dx][y+dy][z+dz][c] = result4[x][y][z][c];
								else result4[x+dx][y+dy][z+dz][c] = Numerics.max(result4[x+dx][y+dy][z+dz][c], result4[x][y][z][c]);
							} else if (merge==MEAN) {
								result4[x+dx][y+dy][z+dz][c] += result4[x][y][z][c];
							}
						}
					}					
					count[x+dx][y+dy][z+dz]++;
				}
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (count[x][y][z]>0) {
				if (merge==MEAN) {
					if (nc==1) result3[x][y][z] /= (float)count[x][y][z];
					else for (int c=0;c<nc;c++) result4[x][y][z][c] /= (float)count[x][y][z];
				}
				mask[x][y][z] = true;
			}
		}

		ImageDataFloat resultData;
		if (nc==1) resultData = new ImageDataFloat(result3);	
		else resultData = new ImageDataFloat(result4);	
		resultData.setHeader(input.getHeader());
		resultData.setName(input.getName()+"_propag");
		resultImage.setValue(resultData);
		resultData = null;
		result3 = null;
		result4 = null;
		
	}
	
}
