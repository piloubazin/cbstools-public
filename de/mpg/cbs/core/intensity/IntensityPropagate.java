package de.mpg.cbs.core.intensity;

import java.net.URL;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis Bazin
 */
 
public class IntensityPropagate {

	// parameters
	public		static final String[]	normtypes = {"max","mean","min"};
	public		static final String[]	targettypes = {"zero","mask","lower","higher"};
	
	// variables
	private float[] inputImage = null;
	private byte[] maskImage = null;
	private float[] resultImage;
	private	String	 normParam = "max";
	private float distParam = 5;
	private		String		targetParam = "zero";
	private float scalingParam = 1.0f;
	
	private int nx, ny, nz, nc, nxyz;
	private float rx, ry, rz;

	// set inputs
	public final void setInputImage(float[] val) { inputImage = val; }
	public final void setMaskImage(byte[] val) { maskImage = val; }
	public final void setCombinationMethod(String val) { normParam = val; }
	public final void setPropagationDistance(float val) { distParam = val; }
	public final void setTargetVoxels(String val) { targetParam = val; }
	public final void setPropogationScalingFactor(float val) { scalingParam = val; }
	
	// set generic inputs	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nc=1; nxyz=nx*ny*nz; }
	public final void setDimensions(int x, int y, int z, int c) { nx=x; ny=y; nz=z; nc=c; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; if (dim.length>3) nc=dim[3]; else nc=1; nxyz=nx*ny*nz; }
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	//set JIST definitions
	//to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; }
	public final String getLabel() { return "Intensity Propogation"; }
	public final String getName() { return "IntensityPropogation"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Propagates the values inside the mask (or non-zero) into the neighboring voxels"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	//set outputs
	public float[] getResultImage() { return resultImage;}
	
	public void execute() {
	
		// different settings
		byte ZERO = 1, MASK = 2, LOWER = 3, HIGHER = 4;
		byte target = ZERO;
		if (targetParam.equals("mask")) target = MASK;
		if (targetParam.equals("lower")) target = LOWER;
		if (targetParam.equals("higher")) target = HIGHER;
		
		// main algorithm
		int nd = Numerics.ceil(distParam/Numerics.min(rx,ry,rz));
		
		byte MIN = 1, MAX = 2, MEAN = 3;
		byte merge = MAX;
		if (normParam.equals("mean")) merge = MEAN;
		if (normParam.equals("min")) merge = MIN;

		// main algorithm
		resultImage = new float[nxyz*nc];
		byte[] count = new byte[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (target==MASK && maskImage[xyz]==0) {
				for (int c=0;c<nc;c++) resultImage[xyz+nxyz*c] = 0.0f;
			} else {
				for (int c=0;c<nc;c++) resultImage[xyz+nxyz*c] = inputImage[xyz+nxyz*c];
			}
		}
		for (int n=0;n<nd;n++) {
			System.out.println("step "+(n+1));
			for (int xyz=0;xyz<nxyz;xyz++) {
				count[xyz] = 0;
			}
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				for (int c=0;c<nc;c++) {								
					// only for places with actual values inside..
					if ( (target==ZERO && inputImage[xyz+nxyz*c]!=0)
						|| (target==MASK && maskImage[xyz]!=0)
						|| target==LOWER || target==HIGHER) {
						// propagate to neighbors
						for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) {
							int xyzd = xyz + dx+nx*dy+nx*ny*dz;
							if ( (target==ZERO && inputImage[xyzd + nxyz*c]==0)
								|| (target==MASK && maskImage[xyzd]==0)
								|| (target==LOWER && inputImage[xyzd + nxyz*c] < scalingParam*inputImage[xyz + nxyz*c])
								|| (target==HIGHER && inputImage[xyzd + nxyz*c] > scalingParam*inputImage[xyz + nxyz*c]) ) {
								//System.out.print(".");
								if (merge==MIN) {
									if (count[xyzd]==0) resultImage[xyzd+nxyz*c] = scalingParam*inputImage[xyz+nxyz*c];
									else resultImage[xyzd+nxyz*c] = Numerics.min(resultImage[xyzd+nxyz*c], scalingParam*inputImage[xyz+nxyz*c]);
								} else if (merge==MAX) {
									if (count[xyzd]==0) resultImage[xyzd+nxyz*c] = scalingParam*inputImage[xyz+nxyz*c];
									else resultImage[xyzd+nxyz*c] = Numerics.max(resultImage[xyzd+nxyz*c], scalingParam*inputImage[xyz+nxyz*c]);
								} else if (merge==MEAN) {
									resultImage[xyzd+nxyz*c] += scalingParam*inputImage[xyz+nxyz*c];
								}
								if (c==nc-1) count[xyzd]++;
							}
						}
					}
				}
			}
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (merge==MEAN && (target==ZERO || target==MASK) && count[xyz]>0) {
					for (int c=0;c<nc;c++) resultImage[xyz+nxyz*c] /= (count[xyz]);
				} else if (merge==MEAN) {
					for (int c=0;c<nc;c++) resultImage[xyz+nxyz*c] /= (1.0f+count[xyz]);
				}
				for (int c=0;c<nc;c++) inputImage[xyz+nxyz*c] = resultImage[xyz+nxyz*c];
				if (target==MASK && inputImage[xyz]!=0) maskImage[xyz] = (byte)1;
			}
		}
	}
}
