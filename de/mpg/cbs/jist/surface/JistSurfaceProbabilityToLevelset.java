package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamDouble;
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
public class JistSurfaceProbabilityToLevelset extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume lvlImage;
	
	private ParamDouble scaleParam;
	
	private ParamVolume probaImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(probaImage = new ParamVolume("Probability Image"));
		inputParams.add(scaleParam = new ParamDouble("Scale (mm)", 0.0f, 100.0f, 5.0f));
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Probability To Level Set");
		inputParams.setName("ProbabilityToLevelSet");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Convert a probability map thresholded at 0.5 into a level set surface (incl. re-initialization).");
		
		info.setVersion("3.0");
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
		System.out.println("probability to levelset");
		ImageDataFloat	probaImg = new ImageDataFloat(probaImage.getImageData());
		
		int nx = probaImg.getRows();
		int ny = probaImg.getCols();
		int nz = probaImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = probaImg.getHeader().getDimResolutions()[0];
		float ry = probaImg.getHeader().getDimResolutions()[1];
		float rz = probaImg.getHeader().getDimResolutions()[2];
		
		float[][][] proba = probaImg.toArray3d();
		
		// transform the values
		float[] levelset = new float[nx*ny*nz];
		//float scale = scaleParam.getValue().floatValue();
		//System.out.println("transform data (scale: "+scale);
		
		// buidl levelset on boundary
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (proba[x][y][z]>=0.5f && (proba[x+1][y][z]<0.5f || proba[x-1][y][z]<0.5f
									   || proba[x][y+1][z]<0.5f || proba[x][y-1][z]<0.5f
									   || proba[x][y][z+1]<0.5f || proba[x][y][z-1]<0.5f)) 
				levelset[xyz] = 0.5f-proba[x][y][z];
			else if (proba[x][y][z]<0.5f && (proba[x+1][y][z]>=0.5f || proba[x-1][y][z]>=0.5f
										   || proba[x][y+1][z]>=0.5f || proba[x][y-1][z]>=0.5f
									   	   || proba[x][y][z+1]>=0.5f || proba[x][y][z-1]>=0.5f))
				levelset[xyz] = 0.5f-proba[x][y][z];
			else if (proba[x][y][z]>=0.5f) levelset[xyz] = -1.0f;
			else levelset[xyz] = +1.0f;
		}
		// no masking
		boolean[] bgmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			bgmask[xyz] = true;
		}
		// expand boundary? yes!
		InflateGdm gdm = new InflateGdm(levelset, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no");
		gdm.evolveNarrowBand(0, 1.0f);
		
		System.out.println("\n done");
		
		ImageDataFloat lvlData = new ImageDataFloat(gdm.exportLevelset());		
		lvlData.setHeader(probaImage.getImageData().getHeader());
		lvlData.setName(probaImage.getImageData().getName()+"_P2L");
		lvlImage.setValue(lvlData);
		lvlData = null;
		levelset = null;
		
	}


}
