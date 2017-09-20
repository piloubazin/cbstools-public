package de.mpg.cbs.jist.surface
;

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
public class JistSurfaceIntensityIntegration extends ProcessingAlgorithm {

	private		static final String[]	normtypes = {"max","mean","min"};
	private		String		normtype = "max";
	
	// jist containers
	private ParamVolume lvlImage;
	private ParamVolume valImage;
	
	private ParamFloat startParam;
	private ParamFloat stopParam;
	private	ParamOption	 normParam;
	
	private ParamVolume integralImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(lvlImage = new ParamVolume("Surface Level Set Image"));
		inputParams.add(valImage = new ParamVolume("Image to integrate"));
		
		inputParams.add(startParam = new ParamFloat("Starting signed distance (voxels)", -1000.0f, 1000.0f, 0.0f));
		inputParams.add(stopParam = new ParamFloat("Stopping signed distance (voxels)", -1000.0f, 1000.0f, 0.0f));
			
		inputParams.add(normParam = new ParamOption("combination method", normtypes));
		normParam.setDescription("use the mean, max or min data from neighboring voxels");
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setLabel("Intensity Integration");
		inputParams.setName("IntensityIntegration");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Integrates an image along the normal to the level set surface.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(integralImage = new ParamVolume("Integrated Image",VoxelType.FLOAT));

		outputParams.setName("integral image");
		outputParams.setLabel("integral image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		System.out.println("levelset input");
		ImageDataFloat	lvlImg = new ImageDataFloat(lvlImage.getImageData());
		
		int nx = lvlImg.getRows();
		int ny = lvlImg.getCols();
		int nz = lvlImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = lvlImg.getHeader().getDimResolutions()[0];
		float ry = lvlImg.getHeader().getDimResolutions()[1];
		float rz = lvlImg.getHeader().getDimResolutions()[2];
		
		float[][][] lvl = lvlImg.toArray3d();
		
		ImageDataFloat	valImg = new ImageDataFloat(valImage.getImageData());
		float[][][] image = valImg.toArray3d();
		
		// integrate along the path: fast marching from start to stop
		float start = startParam.getValue().floatValue();
		float stop = stopParam.getValue().floatValue();
		
		byte MIN = 1, MAX = 2, MEAN = 3;
		byte merge = MAX;
		if (normParam.getValue().equals("mean")) merge = MEAN;
		if (normParam.getValue().equals("min")) merge = MIN;
		System.out.println("merging method "+merge);
		
		// main algorithm
		int nd = Numerics.ceil(Numerics.abs(start-stop));
			
		boolean[][][] mask = new boolean[nx][ny][nz];
		float[][][] result = new float[nx][ny][nz];
		float[][][] weight = new float[nx][ny][nz];
		
		float factor = 1.0f;
		if (start>stop) factor = -1.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (factor*lvl[x][y][z]<factor*start || factor*lvl[x][y][z]>factor*stop) result[x][y][z] = 0.0f;
			else if (merge==MIN || merge==MAX) result[x][y][z] = image[x][y][z];
			else if (merge==MEAN) result[x][y][z] = 0.0f;
			if (factor*lvl[x][y][z]<factor*start) mask[x][y][z] = true;
			else mask[x][y][z] = false;
		}
		System.out.println("growth direction "+factor);
		for (int n=0;n<nd;n++) {
			System.out.println("step "+(n+1));
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				weight[x][y][z] = 0;
			}
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			
				for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) if (!mask[x+dx][y+dy][z+dz]) {
					//System.out.print(".");
					float w = Numerics.max(0.0f, factor*(lvl[x+dx][y+dy][z+dz]-lvl[x][y][z]));
					if (merge==MEAN) {
						result[x+dx][y+dy][z+dz] += w*result[x][y][z];
						weight[x+dx][y+dy][z+dz]+= w;
					} else
					if (merge==MAX) {
						if (w>weight[x+dx][y+dy][z+dz]) {
							if (result[x][y][z]>result[x+dx][y+dy][z+dz]) result[x+dx][y+dy][z+dz] = result[x][y][z];
							weight[x+dx][y+dy][z+dz] = w;
						}
					} else
					if (merge==MIN) {
						if (w>weight[x+dx][y+dy][z+dz]) {
							if (result[x][y][z]<result[x+dx][y+dy][z+dz]) result[x+dx][y+dy][z+dz] = result[x][y][z];
							weight[x+dx][y+dy][z+dz] = w;
						}
					} 
				}
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (weight[x][y][z]>0) {
				if (merge==MEAN) result[x][y][z] = image[x][y][z]+result[x][y][z]/weight[x][y][z];
				mask[x][y][z] = true;
			}
		}

		ImageDataFloat resultData = new ImageDataFloat(result);	
		resultData.setHeader(valImg.getHeader());
		resultData.setName(valImg.getName()+"_sint");
		integralImage.setValue(resultData);
		resultData = null;
		result = null;
	}
	
}
