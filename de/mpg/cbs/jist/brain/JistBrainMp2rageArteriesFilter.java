package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
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
import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMp2rageArteriesFilter extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inv2Image;
	private ParamVolume t1mImage;
	private ParamVolume t1wImage;
	private ParamVolume resultImage;
	
	private static final boolean verbose=true;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(t1mImage = new ParamVolume("T1 map (T1_Images)Image"));
		inputParams.add(t1wImage = new ParamVolume("T1-weighted (UNI) Image"));
		inputParams.add(inv2Image = new ParamVolume("Second inversion (Inv2) Image"));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("MP2RAGE Arteries Filter");
		inputParams.setName("Mp2rageArteriesFilter");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Filters a MP2RAGE brain image to obtain a probability map of arteries.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Arteries Image",VoxelType.FLOAT));
		outputParams.setName("Arteries_Image");
		outputParams.setLabel("Arteries Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat t1m = new ImageDataFloat(t1mImage.getImageData());
		ImageDataFloat t1w = new ImageDataFloat(t1wImage.getImageData());
		ImageDataFloat inv2 = new ImageDataFloat(inv2Image.getImageData());
		float[][][] t1mimg = t1m.toArray3d();
		float[][][] t1wimg = t1w.toArray3d();
		float[][][] inv2img = inv2.toArray3d();
		int nx = inv2.getRows();
		int ny= inv2.getCols();
		int nz = inv2.getSlices();
		float rx = inv2.getHeader().getDimResolutions()[0];
		float ry = inv2.getHeader().getDimResolutions()[1];
		float rz = inv2.getHeader().getDimResolutions()[2];
		
		// main algorithm
		
		// Estimate regions of potential tubular aspect in each image
		byte T1M=0;
		byte T1W=1;
		byte INV=2;
		
		float[][][][] map = new float[nx][ny][nz][3];
		float[] avg = new float[3];
		float[] sig = new float[3];
		float[] den = new float[3];
		avg[T1M] = 0.0f; sig[T1M] = 0.0f; den[T1M] = 0.0f;
		avg[T1W] = 0.0f; sig[T1W] = 0.0f; den[T1W] = 0.0f;
		avg[INV] = 0.0f; sig[INV] = 0.0f; den[INV] = 0.0f;

		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (!zeroNeighbor(t1mimg,x,y,z,2) || !zeroNeighbor(t1wimg,x,y,z,2) || !zeroNeighbor(inv2img,x,y,z,2)) {
				mask[x][y][z] = true;
			} else {
				mask[x][y][z] = false;
			}
		}
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			if (mask[x][y][z]) {
				// check for T1 map
				map[x][y][z][T1M] = 0.0f;
				float best = 0.0f;
				float dmax = 13;
				for (int d=0;d<dmax;d++) {
					float val = -tubularScore(t1mimg, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				map[x][y][z][T1M] = best;
				avg[T1M] += best;
				den[T1M]++;
				// check for T1-weighted
				map[x][y][z][T1W] = 0.0f;
				best = 0.0f;
				for (int d=0;d<dmax;d++) {
					float val = tubularScore(t1wimg, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				map[x][y][z][T1W] = best;
				avg[T1W] += best;
				den[T1W]++;
				// check for INV2
				map[x][y][z][INV] = 0.0f;
				best = 0.0f;
				for (int d=0;d<dmax;d++) {
					float val = tubularScore(inv2img, x,y,z,d);
					if (val*val>best*best) best = val;
				}
				map[x][y][z][INV] = best;
				avg[INV] += best;
				den[INV]++;
			}
		}
		for (int n=0;n<3;n++) {
			avg[n] /= den[n];
		}
		if (verbose) {
			System.out.println("Average T1M response: "+avg[T1M]);
			System.out.println("Average T1W response: "+avg[T1W]);
			System.out.println("Average INV response: "+avg[INV]);
		}
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
			sig[T1M] += Numerics.square(map[x][y][z][T1M]-avg[T1M]);
			sig[T1W] += Numerics.square(map[x][y][z][T1W]-avg[T1W]);
			sig[INV] += Numerics.square(map[x][y][z][INV]-avg[INV]);
		}
		for (int n=0;n<3;n++) {
			sig[n] /= (den[n]-1.0f);
		}
		if (verbose) {
			System.out.println("Stdev T1M response: "+sig[T1M]);
			System.out.println("Stdev T1W response: "+sig[T1W]);
			System.out.println("Stdev INV response: "+sig[INV]);
		}
		// probability
		float[][][] result = new float[nx][ny][nz];
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
			
			// null the negative filter responses (i.e. to zero probability)
			if (map[x][y][z][T1M]<0 || map[x][y][z][T1W]<0 || map[x][y][z][INV]<0) { 
				result[x][y][z] = 0.0f;
			} else { 
				result[x][y][z] = (1.0f - (float)FastMath.exp( -0.5f*(map[x][y][z][T1M]-avg[T1M])*(map[x][y][z][T1M]-avg[T1M])/sig[T1M] ) )
								  *(1.0f - (float)FastMath.exp( -0.5f*(map[x][y][z][T1W]-avg[T1W])*(map[x][y][z][T1W]-avg[T1W])/sig[T1W] ) )
								  *(1.0f - (float)FastMath.exp( -0.5f*(map[x][y][z][INV]-avg[INV])*(map[x][y][z][INV]-avg[INV])/sig[INV] ) );
			}
		}
		float[][][] result2 = new float[nx][ny][nz];
		float ratio = 0.5f;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) if (mask[x][y][z]) {
			result2[x][y][z] = result[x][y][z];
			if (ratio*result[x+1][y][z]>result2[x][y][z]) result2[x][y][z] = ratio*result[x+1][y][z];
			if (ratio*result[x-1][y][z]>result2[x][y][z]) result2[x][y][z] = ratio*result[x-1][y][z];
			if (ratio*result[x][y+1][z]>result2[x][y][z]) result2[x][y][z] = ratio*result[x][y+1][z];
			if (ratio*result[x][y-1][z]>result2[x][y][z]) result2[x][y][z] = ratio*result[x][y-1][z];
			if (ratio*result[x][y][z+1]>result2[x][y][z]) result2[x][y][z] = ratio*result[x][y][z+1];
			if (ratio*result[x][y][z-1]>result2[x][y][z]) result2[x][y][z] = ratio*result[x][y][z-1];
		}
		ImageDataFloat resultData = new ImageDataFloat(result2);		
		resultData.setHeader(t1m.getHeader());
		resultData.setName(t1m.getName()+"_art");
		resultImage.setValue(resultData);
		resultData = null;
		result = null;
		
	}

	boolean zeroNeighbor(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		for (int i=-d;i<=d;i++) for (int j=-d;j<=d;j++) for (int l=-d;l<=d;l++) {
			if (image[x+i][y+j][z+l]==val && i*i+j*j+l*l<=2*d*d) return true;
		}
		return false;
	}
	
	float planeScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {		
			val =(2.0f*image[x][y][z]-image[x-1][y][z]-image[x+1][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z]-image[x+1][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y][z-1]-image[x+1][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y][z+1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z-1]-image[x-1][y-1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x][y-1][z+1]-image[x-1][y-1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x-1][y+1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x][y+1][z+1]-image[x-1][y+1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==1) {
			val =(2.0f*image[x][y][z]-image[x][y-1][z]-image[x][y+1][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z]-image[x-1][y+1][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]-image[x][y-1][z-1]-image[x][y+1][z-1]
				 +2.0f*image[x][y][z+1]-image[x][y-1][z+1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y][z-1]-image[x-1][y-1][z-1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y][z+1]-image[x-1][y-1][z+1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x+1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y][z+1]-image[x+1][y-1][z+1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==2) { 			
			val =(2.0f*image[x][y][z]-image[x][y][z-1]-image[x][y][z+1]
				 +2.0f*image[x-1][y][z]-image[x-1][y][z-1]-image[x-1][y][z+1]
				 +2.0f*image[x+1][y][z]-image[x+1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z]-image[x][y-1][z-1]-image[x][y-1][z+1]
				 +2.0f*image[x][y+1][z]-image[x][y+1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x-1][y-1][z]-image[x-1][y-1][z-1]-image[x-1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z]-image[x-1][y+1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z]-image[x+1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y+1][z]-image[x+1][y+1][z-1]-image[x+1][y+1][z+1])/18.0f;
		} else if (d==3) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z]-image[x+1][y+1][z]
				 +2.0f*image[x-1][y+1][z]-image[x-2][y][z]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y-1][z-1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]-image[x][y+2][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y-1][z+1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]-image[x][y+2][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==4) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z-1]-image[x][y+1][z+1]
				 +2.0f*image[x][y+1][z-1]-image[x][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x][y-1][z+1]-image[x][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]-image[x-1][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]-image[x-1][y][z+2]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]-image[x+1][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]-image[x+1][y][z+2])/18.0f;
		} else if (d==5) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z-1]-image[x+1][y][z+1]
				 +2.0f*image[x+1][y][z-1]-image[x][y][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y][z+1]-image[x-2][y][z]-image[x][y][z+2]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]-image[x+2][y-1][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]-image[x][y-1][z+2]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]-image[x+2][y+1][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]-image[x][y+1][z+2])/18.0f;
		} else if (d==6) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z]-image[x+1][y-1][z]
				 +2.0f*image[x-1][y-1][z]-image[x-2][y][z]-image[x][y-2][z]
				 +2.0f*image[x+1][y+1][z]-image[x][y-2][z]-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]-image[x-1][y+1][z-1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]-image[x][y-2][z-1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]-image[x-1][y+1][z+1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]-image[x][y-2][z+1]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/18.0f;
		} else if (d==7) { 			
			val =(2.0f*image[x][y][z]-image[x][y-1][z+1]-image[x][y+1][z-1]
				 +2.0f*image[x][y-1][z-1]-image[x][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x][y+1][z+1]-image[x][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x-1][y][z]-image[x-1][y-1][z+1]-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]-image[x-1][y][z-2]
				 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]-image[x-1][y+2][z]
				 +2.0f*image[x+1][y][z]-image[x+1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]-image[x+1][y][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]-image[x+1][y+2][z])/18.0f;
		} else if (d==8) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y][z+1]-image[x+1][y][z-1]
				 +2.0f*image[x-1][y][z-1]-image[x-2][y][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y][z+1]-image[x][y][z+2]-image[x+2][y][z]
				 +2.0f*image[x][y-1][z]-image[x-1][y-1][z+1]-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]-image[x][y-1][z-2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]-image[x+2][y-1][z]
				 +2.0f*image[x][y+1][z]-image[x-1][y+1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]-image[x][y+1][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]-image[x+2][y+1][z])/18.0f;
		} else if (d==9) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z-1]-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]-image[x+2][y][z]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]-image[x+2][y+2][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]-image[x][y+2][z+2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==10) { 			
			val =(2.0f*image[x][y][z]-image[x+1][y-1][z-1]-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]-image[x][y+2][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]-image[x-2][y][z]
				 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]-image[x-2][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]-image[x-2][y][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]-image[x][y+2][z+2])/14.0f;
		} else if (d==11) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y+1][z-1]-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]-image[x][y-2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]-image[x+2][y-2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]-image[x][y-2][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]-image[x+2][y][z+2])/14.0f;
		} else if (d==12) { 			
			val =(2.0f*image[x][y][z]-image[x-1][y-1][z+1]-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]-image[x][y][z-2]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]-image[x+2][y][z-2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]-image[x][y+2][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]-image[x+2][y+2][z])/14.0f;
		}
		return 0.5f*val;
	}

	public static float tubularScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z])
						-image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z]
						-image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z]
						-image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1]
						-image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1]
						-image[x][y-1][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1]
						-image[x][y-1][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1]
						-image[x][y+1][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1]
						-image[x][y+1][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1])/24.0f;		
		} else if (d==1) {
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z])
						-image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z]
						-image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z]
						-image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1]
						-image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1]
						-image[x-1][y][z-1]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1]
						-image[x-1][y][z+1]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1]
						-image[x+1][y][z-1]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1]
						-image[x+1][y][z+1]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1])/24.0f;
		} else if (d==2) {
			val =(8.0f*(image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1])
						-image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1]
						-image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1]
						-image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1]
						-image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1]
						-image[x-1][y-1][z]		-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1]
						-image[x-1][y+1][z]		-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1]
						-image[x+1][y-1][z]		-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1]
						-image[x+1][y+1][z]		-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1])/24.0f;
		} else if (d==3) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z])
						-image[x-1][y+1][z]		-image[x-2][y][z]		-image[x][y+2][z]
						-image[x+1][y-1][z]		-image[x][y-2][z]		-image[x+2][y][z]
						-image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-1]		-image[x][y+2][z-1]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
						-image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1]
						-image[x-1][y+1][z+1]	-image[x-2][y][z+1]		-image[x][y+2][z+1]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
		} else if (d==4) {		
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1])
						-image[x][y+1][z-1]		-image[x][y][z-2]		-image[x][y+2][z]
						-image[x][y-1][z+1]		-image[x][y-2][z]		-image[x][y][z+2]
						-image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1]
						-image[x-1][y+1][z-1]	-image[x-1][y][z-2]		-image[x-1][y-2][z]
						-image[x-1][y-1][z+1]	-image[x-1][y-2][z]		-image[x-1][y][z+2]
						-image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x+1][y][z-2]		-image[x+1][y][z+2]
						-image[x+1][y-1][z+1]	-image[x+1][y-2][z]		-image[x+1][y][z+2])/24.0f;
		} else if (d==5) {	
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1])
						-image[x+1][y][z-1]		-image[x][y][z-2]		-image[x+2][y][z]
						-image[x-1][y][z+1]		-image[x-2][y][z]		-image[x][y][z+2]
						-image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1]
						-image[x+1][y-1][z-1]	-image[x][y-1][z-2]		-image[x+2][y-1][z]
						-image[x-1][y-1][z+1]	-image[x-2][y-1][z]		-image[x][y-1][z+2]
						-image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1]
						-image[x+1][y+1][z-1]	-image[x][y+1][z-2]		-image[x+2][y+1][z]
						-image[x-1][y+1][z+1]	-image[x-2][y+1][z]		-image[x][y+1][z+2])/24.0f;
		} else if (d==6) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z])
						-image[x-1][y-1][z]		-image[x-2][y][z]		-image[x][y-2][z]
						-image[x+1][y+1][z]		-image[x][y-2][z]		-image[x+2][y][z]
						-image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y][z-1]		-image[x][y-2][z-1]
						-image[x+1][y+1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
						-image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1]
						-image[x-1][y-1][z+1]	-image[x-2][y][z+1]		-image[x][y-2][z+1]
						-image[x+1][y+1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
		} else if (d==7) {
			val =(8.0f*(image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1])
						-image[x][y-1][z-1]		-image[x][y-2][z]		-image[x][y][z-2]
						-image[x][y+1][z+1]		-image[x][y][z+2]		-image[x][y+2][z]
						-image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1]
						-image[x-1][y-1][z-1]	-image[x-1][y-2][z]		-image[x-1][y][z-2]
						-image[x-1][y+1][z+1]	-image[x-1][y][z+2]		-image[x-1][y+2][z]
						-image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1]
						-image[x+1][y-1][z-1]	-image[x+1][y-2][z]		-image[x+1][y][z-2]
						-image[x+1][y+1][z+1]	-image[x+1][y][z+2]		-image[x+1][y+2][z])/24.0f;
		} else if (d==8) {
			val =(8.0f*(image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1])
						-image[x-1][y][z-1]		-image[x-2][y][z]		-image[x][y][z-2]
						-image[x+1][y][z+1]		-image[x][y][z+2]		-image[x+2][y][z]
						-image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1]
						-image[x-1][y-1][z-1]	-image[x-2][y-1][z]		-image[x][y-1][z-2]
						-image[x+1][y-1][z+1]	-image[x][y-1][z+2]		-image[x+2][y-1][z]
						-image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1]
						-image[x-1][y+1][z-1]	-image[x-2][y+1][z]		-image[x][y+1][z-2]
						-image[x+1][y+1][z+1]	-image[x][y+1][z+2]		-image[x+2][y+1][z])/24.0f;
		} else if (d==9) {
			val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1])
						-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2]
						-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z]
						-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z]
						-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z]
						-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2]
						-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2])/18.0f;
		} else if (d==10) {
			val =(6.0f*(image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1])
						-image[x+1][y-1][z+1]	-image[x+2][y-2][z]		-image[x][y][z+2]
						-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z]
						-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z]
						-image[x-1][y+1][z-1]	-image[x][y][z-2]		-image[x-2][y+2][z]
						-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2]
						-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2])/18.0f;
		} else if (d==11) {				
			val =(6.0f*(image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1])
						 -image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2]
						 -image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z]
						 -image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z]
						 -image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z]
						 -image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2]
						 -image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2])/18.0f;
		} else if (d==12) {				
			val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1])
						-image[x-1][y+1][z+1]	-image[x-2][y][z-2]		-image[x][y+2][z]
						-image[x+1][y-1][z+1]	-image[x][y-2][z+2]		-image[x+2][y][z]
						-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2]
						-image[x+1][y-1][z-1]	-image[x][y-2][z]		-image[x+2][y][z-2]
						-image[x-1][y+1][z-1]	-image[x-2][y][z]		-image[x][y+2][z-2]
						-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z])/18.0f;
		}
		return 0.5f*val;
	}

}
