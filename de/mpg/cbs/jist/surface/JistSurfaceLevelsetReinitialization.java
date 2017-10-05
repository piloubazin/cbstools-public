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


/*
 * @authors Pierre-Louis Bazin, Christine Lucas Tardif
 */
public class JistSurfaceLevelsetReinitialization extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume levelsetImage;
	
	//private ParamOption topologyParam;
	//private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6"};
	
	private ParamVolume reinitLevelsetImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImage = new ParamVolume("Surface Level Set Image"));
		
		//inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		//topologyParam.setValue("no");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Levelset Reinitialization");
		inputParams.setName("LevelsetReinitialization");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Reinitialises the level set distance function.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(reinitLevelsetImage = new ParamVolume("Reinitialized level set",VoxelType.FLOAT));
					
		outputParams.setName("newlevelset");
		outputParams.setLabel("new levelset");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		
		ImageDataFloat	levelsetImg = new ImageDataFloat(levelsetImage.getImageData());
		
		int nx = levelsetImg.getRows();
		int ny = levelsetImg.getCols();
		int nz = levelsetImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = levelsetImg.getHeader().getDimResolutions()[0];
		float ry = levelsetImg.getHeader().getDimResolutions()[1];
		float rz = levelsetImg.getHeader().getDimResolutions()[2];
		
		float[] levelset = new float[nxyz];
		float[][][] buffer = levelsetImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			levelset[xyz] = buffer[x][y][z];
		}
		buffer = null;
		levelsetImg = null;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = true; //(levelset[xyz]>=-5.0f && levelset[xyz]<=5.0f);
		}
		// important: increase the data range (useful for better smoothing)
		//int delta = 30;
		//ObjectMorphology.fastDilateObject(bgmask, nx, ny, nz, delta);
				
		InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
				
		gdm.evolveNarrowBand(0, 1.0f);
		//gdm.fastMarchingReinitialization(false);
		
		// output
		String imgname = levelsetImage.getImageData().getName();
		
		ImageDataFloat reinitLevelsetImg = new ImageDataFloat(gdm.exportLevelset());		
		reinitLevelsetImg.setHeader(levelsetImage.getImageData().getHeader());
		reinitLevelsetImg.setName(imgname+"_reinit");
		reinitLevelsetImage.setValue(reinitLevelsetImg);
		reinitLevelsetImg = null;		
		
	}


}
