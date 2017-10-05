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
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceLevelsetInflation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume levelsetImage;
	
	private ParamFloat 	smoothingParam;
	private ParamFloat 	precisionParam;
	private ParamInteger 	iterationParam;
	private ParamInteger 	smoothiterParam;
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume inflatedImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImage = new ParamVolume("Levelset Image"));
				
		inputParams.add(smoothingParam = new ParamFloat("Smoothing scale (voxels)", 0.0f, 50.0f, 1.0f));
		inputParams.add(smoothiterParam = new ParamInteger("Nb of smoothing iterations", 0, 100000, 20));
		inputParams.add(iterationParam = new ParamInteger("Max. iterations", 0, 100000, 500));
		inputParams.add(precisionParam = new ParamFloat("Precision", 0, 1E9f, 0.001f));
			
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Levelset Inflation");
		inputParams.setName("LevelsetInflation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Inflates a cortical surface using a levelset approach.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(inflatedImage = new ParamVolume("Inflated levelset",VoxelType.FLOAT));
			
		outputParams.setName("inflated images");
		outputParams.setLabel("inflated images");
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
		
		BasicInfo.displayMessage("mask creation...\n");

		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = (levelset[xyz]>=-5.0f && levelset[xyz]<=5.0f);
		}
		// important: increase the data range (useful for better smoothing)
		int delta = 30;
		ObjectMorphology.fastDilateObject(bgmask, nx, ny, nz, delta);
		
		BasicInfo.displayMessage("re-build levelset...\n");

		// main algorithm
		InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, 
										bgmask, 0.4f, 
										0.4f, topologyParam.getValue(),null);
		
		/* seems to work for the inflation part (still a few pbs) */
		double basis = 1.0f;
		double scale = smoothingParam.getValue().floatValue();
		if (smoothiterParam.getValue().floatValue()>1) {
			basis = Math.pow(2.0f*smoothingParam.getValue().floatValue(), 1.0f/(smoothiterParam.getValue().floatValue()-1.0f));
			scale = 0.5f;
		}
		for (int t=0;t<smoothiterParam.getValue().intValue();t++) {
			BasicInfo.displayMessage("step "+(t+1)+"\n");
			
			BasicInfo.displayMessage("input smoothing ("+scale+")...\n");
			gdm.smoothLevelset((float)scale);
    		
			BasicInfo.displayMessage("inflated surface estimation...\n");
			gdm.evolveNarrowBand(iterationParam.getValue().intValue(), precisionParam.getValue().floatValue());
			
			BasicInfo.displayMessage("copy target...\n");
			gdm.updateTarget();
			scale *= basis;
		}
		
		// output
		String imgname = levelsetImage.getImageData().getName();
		
		ImageDataFloat inflatedData = new ImageDataFloat(gdm.exportLevelset());		
		inflatedData.setHeader(levelsetImage.getImageData().getHeader());
		inflatedData.setName(imgname+"_inf");
		inflatedImage.setValue(inflatedData);
		inflatedData = null;		
	}


}
