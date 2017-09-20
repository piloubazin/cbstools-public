package de.mpg.cbs.jist.laminar;

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
public class JistLaminarProfileCalculator extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	intensityImage;
	private ParamVolume 	maskImage;
	private ParamOption		calcParam;
	private ParamFloat		minParam;
	private ParamFloat		maxParam;
	
	private static final String[] calcTypes = {"mean", "stdev", "skewness", "kurtosis"};
	
	private ParamVolume calcImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensityImage = new ParamVolume("Intensity Profile Image",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("Mask Image (opt, 3D or 4D)",null,-1,-1,-1,-1));
		maskImage.setMandatory(false);
		
		inputParams.add(calcParam = new ParamOption("computed statistic", calcTypes));
		
		inputParams.add(minParam=new ParamFloat("Minimum depth", 0, 1, 0.0f));
		minParam.setDescription("minimum depth included in measurement");
		inputParams.add(maxParam=new ParamFloat("Maximum depth", 0, 1, 1.0f));
		maxParam.setDescription("maximum depth included in measurement");
		
		//inputParams.add(centerParam=new ParamFloat("Center (for Gaussian)", -2, 2, 0.0f));
		//centerParam.setDescription("center parameter for Gaussian filtering (profiles normalized in [-1;+1])");
		//inputParams.add(scaleParam=new ParamFloat("Scaling (for Gaussian)", 0, 5, 0.3f));
		//scaleParam.setDescription("scale parameter for Gaussian filtering (profiles normalized in [-1;+1])");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis");
		inputParams.setLabel("Profile Calculator");
		inputParams.setName("ProfileCalculator");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various moments for intensities mapped along a cortical profile.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(calcImage = new ParamVolume("Result", VoxelType.FLOAT));
		
		outputParams.setName("profile calculator");
		outputParams.setLabel("profile calculator");
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
		
		// 4D mask : create layer mask and replace ctxmask
		boolean[][] layermask = new boolean[nlayers+1][nxyz];
		int minlayer = Numerics.floor(minParam.getValue().floatValue()*nlayers);
		int maxlayer = Numerics.ceil(maxParam.getValue().floatValue()*nlayers);
		if (maskImage.getImageData()==null) {
			// create a mask for regions of value 0 on center image (kinda arbitrary..)
			for (int xyz=0;xyz<nx*ny*nz;xyz++) for (int l=minlayer;l<=maxlayer;l++) {
				layermask[l][xyz] = (intensity[nlayers/2][xyz]!=0.0f);
			}
		} else {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			if (maskImg.getComponents()==1) {
				// 3D mask : replace ctxmask
				byte[][][] bytebuffer = maskImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=minlayer;l<=maxlayer;l++) {
					int xyz = x+nx*y+nx*ny*z;
					layermask[l][xyz] = (bytebuffer[x][y][z]!=0);
				}
				bytebuffer = null;
				maskImg = null;
			} else {
				byte[][][][] bytebuffer = maskImg.toArray4d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					for (int l=minlayer;l<=maxlayer;l++) {
						layermask[l][xyz] = (bytebuffer[x][y][z][l]!=0);
					}
				}
				bytebuffer = null;
				maskImg = null;
			}
		}
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) for (int l=minlayer;l<=maxlayer;l++) {
			if (layermask[l][xyz]) ctxmask[xyz] = true;
		}
		
		// main algorithm
		String imgname = intensityImage.getImageData().getName();
				
		float[][][] calc = new float[nx][ny][nz];
		
		/*
		// pre-compute filter
		float center = centerParam.getValue().floatValue();
		float scale = scaleParam.getValue().floatValue();
		float[] filter = new float[nlayers+1];
		
		if (calcParam.getValue().equals("Gaussian")) {
			for (int l=0;l<=nlayers;l++) {
				float x = (l - nlayers/2.0f)/(nlayers/2.0f);
				
				filter[l] = (float)(Math.exp( -(x-center)*(x-center)/(2.0f*scale*scale) )
									/(Math.sqrt(2.0f*Math.PI)*scale) );
			}
		} else if (calcParam.getValue().equals("1st_Gauss_deriv")) {
			for (int l=0;l<=nlayers;l++) {
				float x = (l - nlayers/2.0f)/(nlayers/2.0f);
				
				filter[l] = (float)(Math.exp( -(x-center)*(x-center)/(2.0f*scale*scale) )
									*(x-center)/(Math.sqrt(2.0f*Math.PI)*scale*scale*scale) );
			}
		} else if (calcParam.getValue().equals("2nd_Gauss_deriv")) {
			for (int l=0;l<=nlayers;l++) {
				float x = (l - nlayers/2.0f)/(nlayers/2.0f);
				
				filter[l] = (float)(Math.exp( -(x-center)*(x-center)/(2.0f*scale*scale) )
									*( (x-center)*(x-center)/(scale*scale) - 1.0f)
									/(Math.sqrt(2.0f*Math.PI)*scale*scale*scale) );
			}
		}
		*/
		
		// compute statistics
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				if (calcParam.getValue().equals("mean")) {
					float mean = 0.0f, den = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						mean += intensity[l][xyz];
						den++;
					}
					calc[x][y][z] = mean/den;
				} else if (calcParam.getValue().equals("stdev")) {
					float mean = 0.0f, den = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						mean += intensity[l][xyz];
						den++;
					}
					mean /= den;
					float std = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						std += (intensity[l][xyz]-mean)*(intensity[l][xyz]-mean);
					}
					if (den>1) calc[x][y][z] = (float)Math.sqrt(std/(den-1.0f));
				} else if (calcParam.getValue().equals("skewness")) {
					float mean = 0.0f, den = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						mean += intensity[l][xyz];
						den++;
					}	
					mean /= den;
					float std = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						std += (intensity[l][xyz]-mean)*(intensity[l][xyz]-mean);
					}
					if (den>1) {
						std = (float)Math.sqrt(std/(den-1.0f));
						float skew = 0.0f;
						for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
							skew += (intensity[l][xyz]-mean)*(intensity[l][xyz]-mean)*(intensity[l][xyz]-mean)/(std*std*std);
						}
						calc[x][y][z] = skew/den;
					}
				} else if (calcParam.getValue().equals("kurtosis")) {
					float mean = 0.0f, den = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						mean += intensity[l][xyz];
						den++;
					}
					mean /= den;
					float var = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						var += (intensity[l][xyz]-mean)*(intensity[l][xyz]-mean);
					}
					var /= den;
					float kurt = 0.0f;
					for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
						kurt += (intensity[l][xyz]-mean)*(intensity[l][xyz]-mean)
								*(intensity[l][xyz]-mean)*(intensity[l][xyz]-mean)/(var*var);
					}
					calc[x][y][z] = kurt/den-3.0f;
				} /*else if (calcParam.getValue().contains("Gauss")) {
					if (layermask==null) {
						float data = 0.0f;
						for (int l=0;l<=nlayers;l++) {
							data += filter[l]*intensity[l][xyz];
						}
						calc[x][y][z] = data;	
					} else {
						float data = 0.0f;
						for (int l=0;l<=nlayers;l++) if (layermask[l][xyz]) {
							data += filter[l]*intensity[l][xyz];
						}
						calc[x][y][z] = data;	
					}
				}*/
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
