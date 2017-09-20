package de.mpg.cbs.core.intensity;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensitySpatialDiffusion {

	private float[] input1Image = null;
	private float[] input2Image = null;
	private float[] input3Image = null;
	
	private float[] var1Image = null;
	private float[] var2Image = null;
	private float[] var3Image = null;
	
	
	private byte[] 	maskImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private float[] probaImage = null;
	
	private int 	iterationParam	=	100;
	private float 	imgscaleParam		=	1.0f;
	private	float 	certainscaleParam		= 	1.0f;
	private float 	mincertaintyParam		=	1.0f;
	private int 	neighborParam	=	4;
	private float 	diffratioParam		=	0.01f;
	
		// outputs
	private float[] output1Image = null;
	private float[] output2Image = null;
	private float[] output3Image = null;
	private float[] outprobaImage = null;
	
	// create inputs
	public final void setContrastImage1(float[] val) { input1Image = val; }
	public final void setContrastImage2(float[] val) { input2Image = val; }
	public final void setContrastImage3(float[] val) { input3Image = val; }
	
	public final void setNoiseImage1(float[] val) { var1Image = val; }
	public final void setNoiseImage2(float[] val) { var2Image = val; }
	public final void setNoiseImage3(float[] val) { var3Image = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setImageMask(byte[] val) { maskImage = val; }
	
	public final void setSignalProbabilityImage(float[] val) { probaImage = val; }
	
	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setImageScale(float val) { imgscaleParam = val; }
	public final void setCertaintyScale(float val) { certainscaleParam = val; }
	public final void setMinCertainty(float val) { mincertaintyParam = val; }
	public final void setNeighborhoodSize(int val) { neighborParam = val; }
	public final void setMinDifferenceRatio(float val) { diffratioParam = val; }
	
	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Intensity.devel"; }
	public static final String getLabel() { return "Spatial Intensity Diffusion"; }
	public static final String getName() { return "SID"; }

	public static final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public static final String getAffiliation() { return "Netherlands Institute for Neuroscience | Spinoza Centre for Neuroimaging | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public static final String getDescription() { return "Diffuse image intensities from high to low probability regions."; }
	public static final String getLongDescription() { return getDescription(); }
		
	public static final String getVersion() { return "3.1.1"; };

	// create outputs
	public final float[] getDiffusedProbabilityImage() { return outprobaImage; }
	public final float[] getDiffusedImage1() { return output1Image; }
	public final float[] getDiffusedImage2() { return output2Image; }
	public final float[] getDiffusedImage3() { return output3Image; }
	
	public void execute(){
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image != null) nimg++;
		if (input3Image != null) nimg++;
		
		float[][] image = new float[nimg][nxyz];
		float[][] noise = new float[nimg][nxyz];
		int n = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[n][xyz] = input1Image[xyz];
			noise[n][xyz] = var1Image[xyz];
		}
		input1Image = null;
		var1Image = null;
		if (input2Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input2Image[xyz];
				noise[n][xyz] = var2Image[xyz];
			}
		}	
		input2Image = null;		
		var2Image = null;		
		if (input3Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input3Image[xyz];
				noise[n][xyz] = var3Image[xyz];
			}
		}			
		input3Image = null;		
		var3Image = null;		
				
		// main algorithm
		
		SpatialImageDiffusion sid = new SpatialImageDiffusion(image, noise, probaImage, maskImage, nimg, nx, ny, nz);
		
		BasicInfo.displayMessage("main certainty diffusion...\n");
		sid.diffuseCertainty(iterationParam, imgscaleParam, certainscaleParam, neighborParam, mincertaintyParam, diffratioParam);
				
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
		
		// map to 1D arrays
		output1Image = sid.getImage()[0];
		if (nimg>1) output2Image = sid.getImage()[1];
		if (nimg>2) output3Image = sid.getImage()[2];
		
		outprobaImage = sid.getProba();
		
		return;
	}

}
