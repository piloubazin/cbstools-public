package de.mpg.cbs.jist.laminar;

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


/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarProfileDepthMapping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	intensityImage;
	private ParamVolume 	maskImage;
	//private ParamOption		calcParam;
	private ParamFloat		mindepthParam;
	private ParamFloat		maxdepthParam;
	
	//private static final String[] calcTypes = {"mean", "stdev", "skewness", "kurtosis"};
	
	private ParamVolume calcImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensityImage = new ParamVolume("Intensity Profile Image",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("Mask Image (opt, 3D or 4D)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		
		//inputParams.add(calcParam = new ParamOption("computed statistic", calcTypes));
		
		inputParams.add(mindepthParam=new ParamFloat("Minimum relative depth", 0, 1, 0.0));
		mindepthParam.setDescription("minimum depth to be mapped (in [0,1] from CSF boundary to WM boundary)");
		inputParams.add(maxdepthParam=new ParamFloat("Maximum relative depth", 0, 1, 0.0));
		maxdepthParam.setDescription("maximum depth to be mapped (in [0,1] from CSF boundary to WM boundary)");
		//inputParams.add(scaleParam=new ParamDouble("Scaling (for Gaussian)", 0, 5, 0.3));
		//scaleParam.setDescription("scale parameter for Gaussian filtering (profiles normalized in [-1;+1])");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis");
		inputParams.setLabel("Profile Depth Mapping");
		inputParams.setName("ProfileDepthMapping");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Marcel Weiss", "weiss@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute an average profile within a range of cortical depth.");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(calcImage = new ParamVolume("Result", VoxelType.FLOAT));
		
		outputParams.setName("profile depth mapping");
		outputParams.setLabel("profile depth mapping");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays : TO DO
		
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		
		int nx = intensImg.getRows();
		int ny = intensImg.getCols();
		int nz = intensImg.getSlices();
		int nlayers = intensImg.getComponents()-1;
		int nxyz = nx*ny*nz;
		float rx = intensImg.getHeader().getDimResolutions()[0];
		float ry = intensImg.getHeader().getDimResolutions()[1];
		float rz = intensImg.getHeader().getDimResolutions()[2];
		
		float[][] intensity = new float[nlayers+1][nxyz];
		float[][][][] buffer4 = intensImg.toArray4d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			intensity[l][xyz] = buffer4[x][y][z][l];
		}
		buffer4 = null;
		intensImg = null;
		
		boolean[] ctxmask = new boolean[nxyz];
		boolean[][] layermask = null;
		if (maskImage.getImageData()==null) {
			// create a mask for regions of value 0 on center image (kinda arbitrary..)
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				ctxmask[xyz] = (intensity[nlayers/2][xyz]!=0.0);
			}
		} else {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			if (maskImg.getComponents()==1) {
				// 3D mask : replace ctxmask
				byte[][][] bytebuffer = maskImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					ctxmask[xyz] = (bytebuffer[x][y][z]!=0);
				}
				bytebuffer = null;
				maskImg = null;
			} else {
				// 4D mask : create layer mask and replace ctxmask
				layermask = new boolean[nlayers+1][nxyz];
				byte[][][][] bytebuffer = maskImg.toArray4d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					ctxmask[xyz] = false;
					for (int l=0;l<=nlayers;l++) {
						layermask[l][xyz] = (bytebuffer[x][y][z][l]!=0);
						if (layermask[l][xyz]) ctxmask[xyz] = true;
					}
				}
				bytebuffer = null;
				maskImg = null;
			}
		}
				
		// main algorithm
		String imgname = intensityImage.getImageData().getName();
				
		float[][][] calc = new float[nx][ny][nz];
		float mindepth = mindepthParam.getValue().floatValue();
		float maxdepth = maxdepthParam.getValue().floatValue();
				
		// compute statistics
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				if (layermask==null) {
					float mean = 0.0f;
					for (int l=0;l<=nlayers;l++) mean += intensity[l][xyz];
					calc[x][y][z] = mean/(nlayers+1.0f);
				} else {
					float mean = 0.0f, den = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						mean += intensity[l][xyz];
						den++;
					}
					calc[x][y][z] = mean/den;
				}		
			}
		}
		ImageDataFloat calcData = new ImageDataFloat(calc);		
		calcData.setHeader(intensityImage.getImageData().getHeader());
		calcData.setName(imgname+"_"+calcParam.getValue());
		calcImage.setValue(calcData);
		calcData = null;
		calc = null;
	}


}
