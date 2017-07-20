package de.mpg.cbs.jist.cortex;

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
public class JistCortexCentralSurface extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume lv1Image;
	private ParamVolume lv2Image;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 	regularizationParam;
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume avgImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(lv1Image = new ParamVolume("GM-WM Boundary Levelset Image"));
		imageParams.add(lv2Image = new ParamVolume("CSF-GM Boundary Levelset Image"));
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		
		mainParams.add(regularizationParam = new ParamFloat("Regularization ratio (in [0,1])", 0.0f, 1.0f, 0.1f));
			
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing.prev");
		inputParams.setLabel("Central Surface");
		inputParams.setName("CentralSurface");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimate a central levelset surface from the inner and outer cortical boundaries with consistent topology.");
		
		info.setVersion("2.0f");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(avgImage = new ParamVolume("Central Levelset Image",VoxelType.FLOAT));

		outputParams.setName("central levelset image");
		outputParams.setLabel("central levelset image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		
		ImageDataFloat	lv1Img = new ImageDataFloat(lv1Image.getImageData());
		ImageDataFloat	lv2Img = new ImageDataFloat(lv2Image.getImageData());
		
		int nx = lv1Img.getRows();
		int ny = lv1Img.getCols();
		int nz = lv1Img.getSlices();
		int nxyz = nx*ny*nz;
		float rx = lv1Img.getHeader().getDimResolutions()[0];
		float ry = lv1Img.getHeader().getDimResolutions()[1];
		float rz = lv1Img.getHeader().getDimResolutions()[2];
		
		float[] lv1 = new float[nxyz];
		float[][][] buffer = lv1Img.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lv1[xyz] = buffer[x][y][z];
		}
		buffer = null;
		lv1Img = null;
		
		float[] lv2 = new float[nxyz];
		buffer = lv2Img.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lv2[xyz] = buffer[x][y][z];
		}
		buffer = null;
		lv2Img = null;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = (lv1[xyz]>=-5.0f && lv2[xyz]<=5.0f);
		}
		float balloon = 0.8f*(1.0f-regularizationParam.getValue().floatValue());
		float curv = 0.8f*(regularizationParam.getValue().floatValue());
		
		float[] ctr = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctr[xyz] = 0.5f*(lv1[xyz]+lv2[xyz]);
		}
		
		// main algorithm
		SmoothGdm gdm = new SmoothGdm(lv1, ctr, nx, ny, nz, rx, ry, rz,
											bgmask, balloon, curv, 	topologyParam.getValue(), null);
		//CorticalLayersGdm gdm = new CorticalLayersGdm(lv1, lv2, 0.5f, nx, ny, nz, rx, ry, rz,
		//												bgmask, balloon, curv, topologyParam.getValue());
			

		BasicInfo.displayMessage("average estimation...\n");
		//gdm.setFraction(0.5f);
		gdm.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		float[] average = gdm.exportLevelset();
		
		// output
		String imgname = lv1Image.getImageData().getName();
		
		float[][][] avg = new float[nx][ny][nz];
		float max = 0.0f;
		for (int xyz=0;xyz<nxyz;xyz++) 
			if (bgmask[xyz] && average[xyz]>max) max = average[xyz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (bgmask[xyz]) avg[x][y][z] = average[xyz];
			else if (lv1[xyz]<=0) avg[x][y][z] = -max;
			else if (lv2[xyz]>=0) avg[x][y][z] = +max;
		}
		average = null;
		
		ImageDataFloat avgData = new ImageDataFloat(avg);		
		avgData.setHeader(lv1Image.getImageData().getHeader());
		avgData.setName(imgname+"_ctr");
		avgImage.setValue(avgData);
		avgData = null;
		avg = null;
		
	}


}
