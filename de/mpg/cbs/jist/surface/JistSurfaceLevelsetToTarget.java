package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceLevelsetToTarget extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume initGdmImage;
	private ParamVolume targetGdmImage;
	private ParamVolume resultGdmImage;
	
	private	ParamFloat 	balloonParam;
	private ParamFloat 	curvParam;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 	changeParam;
	
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(initGdmImage = new ParamVolume("Initial Level Set Surface"));
		inputParams.add(targetGdmImage = new ParamVolume("Target Level Set Surface"));
		
		inputParams.add(balloonParam = new ParamFloat("Data weight", -1E10f, 1E10f, 0.5f));
		inputParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.2f));
		
		inputParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		inputParams.add(changeParam = new ParamFloat("Min change", 0, 1, 0.001f));
			
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Level Set To Target");
		inputParams.setName("LevelSetToTarget");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Evolves a level set surface toward a target boundary, while preserving smoothness and topology.");
		
		info.setVersion("3.0.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultGdmImage = new ParamVolume("Final Level Set Surface",VoxelType.FLOAT));
		
		outputParams.setLabel("GDM To Target");
		outputParams.setName("GdmToTarget");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat initgdmImg = new ImageDataFloat(initGdmImage.getImageData());
		int nx = initgdmImg.getRows();
		int ny = initgdmImg.getCols();
		int nz = initgdmImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = initgdmImg.getHeader().getDimResolutions()[0];
		float ry = initgdmImg.getHeader().getDimResolutions()[1];
		float rz = initgdmImg.getHeader().getDimResolutions()[2];
		float[][][] buffer = initgdmImg.toArray3d();
		float[] initgdm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			initgdm[xyz] = buffer[x][y][z];
		}
		buffer = null;
		
		ImageDataFloat targetGdmImg = new ImageDataFloat(targetGdmImage.getImageData());
		buffer = targetGdmImg.toArray3d();
		float[] targetgdm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			targetgdm[xyz] = buffer[x][y][z];
		}
		buffer = null;
		
		// data mask : exclude anything too far away
		boolean[] mask = new boolean[nxyz];
		float maxdist = 10.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (Numerics.abs(initgdm[xyz])>maxdist && Numerics.abs(targetgdm[xyz])>maxdist) 
				mask[xyz] = false;
			else 
				mask[xyz] = true;
		}
			
		// MGDM evolution		
		SmoothGdm evolve = new SmoothGdm(initgdm, targetgdm, nx, ny, nz, rx, ry, rz,
											mask, balloonParam.getValue().floatValue(), 
											curvParam.getValue().floatValue(), 
											topologyParam.getValue(),null);
		
		evolve.fastMarchingReinitialization(true);
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			evolve.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());
		}
		
		ImageDataFloat lvlData = new ImageDataFloat(evolve.exportLevelset3d());	
		lvlData.setHeader(initgdmImg.getHeader());
		lvlData.setName(initgdmImg.getName()+"_2trg");
		resultGdmImage.setValue(lvlData);
		lvlData = null;
	}

}
