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
public class JistSurfaceLevelsetToProbability extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume lvlImage;
	
	private ParamFloat scaleParam;
	
	private ParamVolume probaImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(lvlImage = new ParamVolume("Surface Level Set Image"));
		
		inputParams.add(scaleParam = new ParamFloat("Scale (mm)", 0.0f, 100.0f, 5.0f));
			
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Level Set To Probability");
		inputParams.setName("LevelSetToProbability");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Convert a level set surface representation to a smooth probability map (0 outside, 1 inside).");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(probaImage = new ParamVolume("Probability Image",VoxelType.FLOAT));

		outputParams.setName("probability image");
		outputParams.setLabel("probability image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		System.out.println("levelset to probability");
		ImageDataFloat	lvlImg = new ImageDataFloat(lvlImage.getImageData());
		
		int nx = lvlImg.getRows();
		int ny = lvlImg.getCols();
		int nz = lvlImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = lvlImg.getHeader().getDimResolutions()[0];
		float ry = lvlImg.getHeader().getDimResolutions()[1];
		float rz = lvlImg.getHeader().getDimResolutions()[2];
		
		float[][][] lvl = lvlImg.toArray3d();
		
		// transform the values
		float[][][] proba = new float[nx][ny][nz];
		float scale = scaleParam.getValue().floatValue();
		System.out.println("transform data (scale: "+scale);
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			proba[x][y][z] = (float)(1.0f/(1.0f+FastMath.exp(lvl[x][y][z]*rx/scale)));
		}
		System.out.println("\n done");
		
		ImageDataFloat probaData = new ImageDataFloat(proba);		
		probaData.setHeader(lvlImage.getImageData().getHeader());
		probaData.setName(lvlImage.getImageData().getName()+"_L2P");
		probaImage.setValue(probaData);
		probaData = null;
		proba = null;
		
	}


}
