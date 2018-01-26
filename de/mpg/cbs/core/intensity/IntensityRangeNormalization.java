package de.mpg.cbs.core.intensity;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
public class IntensityRangeNormalization {

	// parameters
	public		static final String[]	normtypes = {"linear","robust","robust-min","robust-max"};
	private		String		normtype = "robust";
	
	// jist containers
	private float[] inImage;
	private int[] maskImage = null;
	
	private float[] resultImage;

	private int nx, ny, nz, nxyz, nt;
	private float rx, ry, rz;
	
	private	String	 normParam = normtype;
	private float	 ratioParam = 0.01f;
	private boolean ignoreNegParam = true;
	private boolean ignoreZeroParam = true;
	private float	 scalingParam = 1.0f;
	
	// create inputs
	public final void setInputImage(float[] val) { inImage = val; }
	public final void setMaskImage(int[] val) { maskImage = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setNormalization(String val) { normParam = val; }
	public final void setRobustnessRatio(float val) { ratioParam = val; }
	public final void setNegativeValuesToZero(boolean val) { ignoreNegParam = val; }
	public final void setIgnoreZeroValues(boolean val) { ignoreZeroParam = val; }
	public final void setOutputScaling(float val) { scalingParam = val; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; }
	public final String getLabel() { return "Intensity Range Normalization"; }
	public final String getName() { return "IntensityRangeNormalization"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for Neuroimaging, Netherlands Institute for Neuroscience, Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Outputs a normalized version of the image in [0,S]."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
	
	// create outputs
	public final float[] getNormalizedImage() { return resultImage; }
	
	public final void execute(){
				
		// use a mask by default
		boolean[] mask = new boolean[nxyz];
		if (maskImage!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				mask[xyz] = (maskImage[xyz]!=0);
			}
		} else {
			for (int xyz=0;xyz<nxyz;xyz++) {
				mask[xyz]  =true;
			}
		}
		// negative to zero
		if (ignoreNegParam) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (inImage[xyz]<0) inImage[xyz] = 0.0f;
			}
		}
		// ignore zero values: masked
		if (ignoreZeroParam) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (inImage[xyz]==0) mask[xyz] = false;
			}
		}
		    
		// main algorithm
		
		// 1. estimate the min, max
		float Imin = 0, Imax = 0;
		BasicInfo.displayMessage("normalization method: "+normParam+"\n");
		
		if (normParam.equals("linear")) {
		    Imin = ImageStatistics.minimum(inImage, mask, nx, ny, nz); 
			Imax = ImageStatistics.maximum(inImage, mask, nx, ny, nz);				
		} else if (normParam.equals("robust")) {
			Imin = ImageStatistics.robustMinimum(inImage, mask, ratioParam, 4, nx, ny, nz); 
			Imax = ImageStatistics.robustMaximum(inImage, mask, ratioParam, 4, nx, ny, nz);
		} else if (normParam.equals("robust-min")) {
			Imin = ImageStatistics.robustMinimum(inImage, mask, ratioParam, 4, nx, ny, nz); 
			Imax = ImageStatistics.maximum(inImage, mask, nx, ny, nz);
		} else if (normParam.equals("robust-max")) {
			Imin = 0.0f;
			Imax = ImageStatistics.robustMaximum(inImage, mask, ratioParam, 4, nx, ny, nz);
		}
		BasicInfo.displayMessage("image min, max: "+Imin+", "+Imax+"\n");
			
		// 2. scale the data
		float[] result = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			result[xyz] = scalingParam*Numerics.bounded( (inImage[xyz]-Imin)/(Imax-Imin), 0.0f, 1.0f);
		}
		resultImage = result;
	}
	
}
