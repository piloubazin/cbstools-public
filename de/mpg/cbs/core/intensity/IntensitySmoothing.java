package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensitySmoothing {

	// input parameters
	private		float[] 	image;
	private		byte[] 	mask = null;
	private		int		nx, ny, nz, nxyz;
	private 		float 	rx, ry, rz;

	private		float		scale = 1.0f;
	private		String		method = "gaussian";
	private 	boolean		skip0 = true;
	
	// output parameters
	private		float[] smoothed = null;
	
	public		static final String[]	methods = {"gaussian"};
	
	// set inputs
	public final void setInputImage(float[] in) { image = in; }
	public final void setInputMask(byte[] in) { mask = in; }
	public final void setSmoothingScale_mm(float s) { scale = s; }
	public final void setSmoothingMethod(String m) { method = m; }
	public final void setSkipZeroValues(boolean sk) { skip0 = sk; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; } 
	public final String getLabel() { return "Smoothing"; }
	public final String getName() { return "Smoothing"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Basic intensity smoothing."; }
		
	public final String getVersion() { return "3.1"; }

	// get outputs
	public float[] getSmoothedImage() { return smoothed; }
	
	
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// build a boolean mask
		boolean[] bmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (mask!=null) bmask[xyz] = (mask[xyz]>0);
			else bmask[xyz] = true;
			if (skip0 && image[xyz]==0) bmask[xyz] = false;
		}
		
		if (method.equals("gaussian")) {
			float[][] kernel = ImageFilters.separableGaussianKernel(scale/rx,scale/ry,scale/rz);
			 smoothed =  ImageFilters.separableMaskedConvolution(image, bmask, nx, ny, nz, kernel);
		}
		
		return;
	}


}
