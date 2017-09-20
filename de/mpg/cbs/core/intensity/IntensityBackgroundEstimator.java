package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensityBackgroundEstimator {

	// input parameters
	private		float[] 	image;
	private		int			nx, ny, nz, nxyz;
	private		float		ratio = 0.001f;
	private		String		distribution = "exponential";
	private 	boolean		skip0 = true;
	
	// output parameters
	private		byte[] mask;
	private		float[] masked;
	private		float[] proba;
	
	public		static final String[]	distributions = {"exponential","half-normal"};
	private		static final byte	HNORM = 1;
	private		static final byte	EXP = 2;
	
	
	
	// set inputs
	public final void setInputImage(float[] in) { image = in; }
	public final void setRobustMinMaxThresholding(float r) { ratio = r; }
	public final void setBackgroundDistribution(String dist) { distribution = dist; }
	public final void setSkipZeroValues(boolean sk) { skip0 = sk; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; } 
	public final String getLabel() { return "Background Estimator"; }
	public final String getName() { return "BackgroundEstimator"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimate the background data region."; }
		
	public final String getVersion() { return "3.1"; }

		// get outputs
	public float[] getMaskedImage() { return masked; }
	public float[] getProbaImage() { return proba; }
	public byte[] getMask() { return mask; }
	
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm

		// assume  an exponential + outlier distributions for the image
		
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 50;
		
		// min max; optionally a robust measurement 
		float min, max;
		if (ratio>0) {
			min = ImageStatistics.robustMinimum(image, ratio, 5, nx, ny, nz);
			max = ImageStatistics.robustMaximum(image, ratio, 5, nx, ny, nz);
		} else {
			min = ImageStatistics.minimum(image, nx, ny, nz);
			max = ImageStatistics.maximum(image, nx, ny, nz);
		}
		// re-normalize
		for (int xyz=0;xyz<nxyz;xyz++) {
			image[xyz] = Numerics.bounded( (image[xyz]-min)/(max-min), 0.0f, 1.0f);
		}
		System.out.println("image range: "+min+", "+max);

		double mean = 0.0f;
		double den = 0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[xyz]!=0) {
			mean += image[xyz];
			den++;
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		double pi = Math.PI;
		
		byte model = HNORM;
		if (distribution.equals("exponential")) model = EXP;
		
		// re-normalized probability map
		proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*FastMath.exp(-image[xyz]*image[xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(FastMath.exp(-image[xyz]/mean)/mean);
			// bg vs. outlier test : proba is p(xyz is background)
			proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
		}
		// loop
		double diff = 1.0;
		double prev;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			prev = mean;
			mean = 0.0;
			den = 0.0;
			for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[xyz]!=0) {
				mean += proba[xyz]*image[xyz];
				den += proba[xyz];
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(prev-mean)/(0.5f*(prev+mean));
			System.out.println("diff parameters: "+diff);
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (model==HNORM)
					// half-normal model
					proba[xyz] = (float)(2.0/(mean*pi)*FastMath.exp(-image[xyz]*image[xyz]/(pi*mean*mean)));
				else if (model==EXP)
					// exponential model
					proba[xyz] = (float)(FastMath.exp(-image[xyz]/mean)/mean);
				// bg vs. outlier test : proba is p(xyz is background)
				proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
			}
		}
		
		mask = new byte[nxyz];
		masked = new float[nxyz];
			
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (proba[xyz]<=0.5f) {
				mask[xyz] = 1;
				masked[xyz] = image[xyz];
			} else {
				mask[xyz] = 0;
				masked[xyz] = 0.0f;
			}
			proba[xyz] = 1.0f-proba[xyz];
		}
		
		return;
	}


}
