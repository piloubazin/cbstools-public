package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
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

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfacePartialVolumeToLevelset extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume lvlImage;
	
	private ParamFloat scaleParam;
	
	private ParamVolume pvImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(pvImage = new ParamVolume("Partial Volume Image"));
		inputParams.add(scaleParam = new ParamFloat("Scale (mm)", 0.0f, 100.0f, 5.0f));
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Partial Volume To Level Set");
		inputParams.setName("PartialVolumeToLevelSet");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Convert a partial voluming map into a level set surface (incl. re-initialization).");
		
		info.setVersion("3.1f");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(lvlImage = new ParamVolume("Level Set Image",VoxelType.FLOAT));

		outputParams.setName("levelset image");
		outputParams.setLabel("levelset image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		System.out.println("partial voluming to levelset");
		ImageDataFloat	pvImg = new ImageDataFloat(pvImage.getImageData());
		
		int nx = pvImg.getRows();
		int ny = pvImg.getCols();
		int nz = pvImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = pvImg.getHeader().getDimResolutions()[0];
		float ry = pvImg.getHeader().getDimResolutions()[1];
		float rz = pvImg.getHeader().getDimResolutions()[2];
		
		float[][][] pv = pvImg.toArray3d();
		
		// transform the values
		float[] levelset = new float[nx*ny*nz];
		//float scale = scaleParam.getValue().floatValue();
		//System.out.println("transform data (scale: "+scale);
		
		// buidl levelset on boundary
		float epsilon = 1e-9f;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (pv[x][y][z]>epsilon && (pv[x+1][y][z]<=epsilon || pv[x-1][y][z]<=epsilon
									   || pv[x][y+1][z]<=epsilon|| pv[x][y-1][z]<=epsilon
									   || pv[x][y][z+1]<=epsilon|| pv[x][y][z-1]<=epsilon)) 
				levelset[xyz] = -0.5f/(2.0f/pv[x][y][z]-1.0f);
			else if (pv[x][y][z]<=epsilon&& (pv[x+1][y][z]>epsilon || pv[x-1][y][z]>epsilon
										   || pv[x][y+1][z]>epsilon || pv[x][y-1][z]>epsilon
									   	   || pv[x][y][z+1]>epsilon || pv[x][y][z-1]>epsilon))
				levelset[xyz] = 1.0f;
			else if (pv[x][y][z]>epsilon) levelset[xyz] = -1.0f;
			else levelset[xyz] = +1.0f;
		}
		// no masking
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = true;
		}
		// expand boundary? yes!
		InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		gdm.evolveNarrowBand(0, 1.0f);
		
		System.out.println("\n done");
		
		ImageDataFloat lvlData = new ImageDataFloat(gdm.exportLevelset());		
		lvlData.setHeader(pvImage.getImageData().getHeader());
		lvlData.setName(pvImage.getImageData().getName()+"_PV2L");
		lvlImage.setValue(lvlData);
		lvlData = null;
		levelset = null;
		
	}


}
