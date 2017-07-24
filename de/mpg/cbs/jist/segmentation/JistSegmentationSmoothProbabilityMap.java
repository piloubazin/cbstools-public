package de.mpg.cbs.jist.segmentation;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationSmoothProbabilityMap extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume labelImage;
	private ParamVolume probaImage;
	private ParamFile atlasParam;
	private ParamFloat scaleParam;
	private ParamFloat factorParam;
	
	private ParamVolume labelSmoothImage;
	private ParamVolume probaSmoothImage;
	
	private static final byte EMPTY=-1;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(probaImage = new ParamVolume("Probability Function Image"));
		inputParams.add(labelImage = new ParamVolume("Probability Label Image"));
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		inputParams.add(scaleParam = new ParamFloat("Scale", 0.0f, 100.0f, 0.5f));
		inputParams.add(factorParam = new ParamFloat("Factor", 0.0f, 100.0f, 0.5f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Smooth Max Probability Map");
		inputParams.setName("SmoothMaxProbabilityMap");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooth a multi-object maximum probability map with a MRF diffusion technique.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(labelSmoothImage = new ParamVolume("Smoothed Labels",VoxelType.BYTE,-1,-1,-1,-1));
		outputParams.add(probaSmoothImage = new ParamVolume("Smoothed Probability",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("smooth probability images");
		outputParams.setLabel("smooth probability images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat probaImg = new ImageDataFloat(probaImage.getImageData());
		int nx = probaImg.getRows();
		int ny = probaImg.getCols();
		int nz = probaImg.getSlices();
		int nc = probaImg.getComponents();
		int nxyz = nx*ny*nz;
		BasicInfo.displayMessage("Image dims: "+nx+", "+ny+", "+nz+"; "+nc+"\n");
		float rx = probaImg.getHeader().getDimResolutions()[0];
		float ry = probaImg.getHeader().getDimResolutions()[1];
		float rz = probaImg.getHeader().getDimResolutions()[2];
		
		int orient = probaImg.getHeader().getImageOrientation().ordinal();
		int orx = probaImg.getHeader().getAxisOrientation()[0].ordinal();
		int ory = probaImg.getHeader().getAxisOrientation()[1].ordinal();
		int orz = probaImg.getHeader().getAxisOrientation()[2].ordinal();
		
		ImageDataByte	labelImg = new ImageDataByte(labelImage.getImageData());
		
		// load mask and build boolean signature for each region
		BasicInfo.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam.getValue().getAbsolutePath());

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		
		// new multi-proba: empty
		float[][][][] smoothproba = new float[nx][ny][nz][nc];
		byte[][][][] smoothlabel = new byte[nx][ny][nz][nc];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int c=0;c<nc;c++) {
			smoothlabel[x][y][z][c] = EMPTY;
		}
		
		// proceed object by object
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
		
			BasicInfo.displayMessage("Object label "+atlas.getLabels()[nobj]+"\n");
	
			// reconstruct the probabilities
			float[][][] proba = new float[nx][ny][nz];
			boolean[][][] mask = new boolean[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				proba[x][y][z] = 0.0f;
				for (int c=0;c<nc;c++) {
					if (labelImg.getByte(x,y,z,c)==atlas.getLabels()[nobj]) {
						proba[x][y][z] = probaImg.getFloat(x,y,z,c);
					}
				}
				// mask out values with absolute certainty (0 and 1) - for speed -
				if (proba[x][y][z]!=0.0f && proba[x][y][z]!=1.0f) mask[x][y][z] = true;
				else mask[x][y][z] = false;
			}
			
			// smooth things
			BasicInfo.displayMessage("diffusion smoothing ("+scaleParam.getValue().floatValue()+", "+factorParam.getValue().floatValue()+")\n");
			proba = probaDiffusion(proba, mask, nx, ny, nz, 20, scaleParam.getValue().floatValue(), factorParam.getValue().floatValue());
			
			// store to multi-label image
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				for (int c=0;c<nc;c++) {
					if (smoothlabel[x][y][z][c]==EMPTY) {
						smoothproba[x][y][z][c] = proba[x][y][z];
						smoothlabel[x][y][z][c] = atlas.getLabels()[nobj];
						c = nc;
					} else if (smoothproba[x][y][z][c]<proba[x][y][z]) {
						for (int d=nc-1;d>c;d--) {
							smoothlabel[x][y][z][d] = smoothlabel[x][y][z][d-1];
							smoothproba[x][y][z][d] = smoothproba[x][y][z][d-1];
						}
						smoothproba[x][y][z][c] = proba[x][y][z];
						smoothlabel[x][y][z][c] = atlas.getLabels()[nobj];
						c = nc;
					}
				}
			}
					
			
		}	
		
		ImageDataFloat resData = new ImageDataFloat(smoothproba);	
		smoothproba = null;
		resData.setHeader(probaImg.getHeader());
		resData.setName(probaImg.getName()+"_smooth");
		probaSmoothImage.setValue(resData);
		
		ImageDataByte lbData = new ImageDataByte(smoothlabel);	
		smoothlabel = null;
		lbData.setHeader(labelImg.getHeader());
		lbData.setName(labelImg.getName()+"_smooth");
		labelSmoothImage.setValue(lbData);
				
		return;
	}
	
	private final float[][][] probaDiffusion(float[][][] orig, boolean[][][] mask, int nx, int ny, int nz, int iter, float scale, float factor) {
    	
		float[][][] map = new float[nx][ny][nz];
    	for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
    		map[x][y][z] = orig[x][y][z];
    	}
				
    	// propagate the values : diffusion
    	float maxdiff = 1.0f;
		for (int t=0;t<iter && maxdiff>0.05f;t++) {
			BasicInfo.displayMessage(".");
			maxdiff = 0.0f;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				if (mask[x][y][z]) {
					//float den = 1.0f;
					float den = diffusionWeightFunction(map[x][y][z],orig[x][y][z],scale);
					float num = den*orig[x][y][z];
					float weight;
					float prev = map[x][y][z];
						
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x-1][y][z],scale);
					num += weight*map[x-1][y][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x+1][y][z],scale);
					num += weight*map[x+1][y][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y-1][z],scale);
					num += weight*map[x][y-1][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y+1][z],scale);
					num += weight*map[x][y+1][z];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y][z-1],scale);
					num += weight*map[x][y][z-1];
					den += weight;
					
					weight = factor/6.0f*diffusionWeightFunction(map[x][y][z],map[x][y][z+1],scale);
					num += weight*map[x][y][z+1];
					den += weight;
	
					map[x][y][z] = num/den;
					
					maxdiff = Numerics.max(maxdiff, Numerics.abs(map[x][y][z]-prev));
				} else {
					map[x][y][z] = orig[x][y][z];
				}
			}
			BasicInfo.displayMessage("max diff. "+maxdiff+"\n");
		}
		return map;
    }
    
    private final float diffusionWeightFunction(float val, float ngb, float scale) {
    	
    	//return 1.0f;
    	//return 1.0f/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	//return 1.0f/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	// making it exact for the cases p=0 or 1 ?
    	return Numerics.square(Numerics.min(1.0f-val, val)/scale)/(1.0f+Numerics.square(Numerics.min(1.0f-val, val)/scale))
    				/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    }

}
