package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class IntensityMp2rageMasking {

	// containers
	private float[] inv2Image;
	private float[] t1mapImage = null;
	private float[] isoImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private float[] probaImage;
	private int[] maskImage;
	private float[] t1maskImage;
	private float[] isomaskImage;
	
	// parameters
	public		static final String[]	distribs = {"exponential","half-normal"};
	public		static final String		distrib = "half-normal";
	
	public		static final String[]	maskings = {"binary","dilated","proba"};
	public		static final String		masking = "binary";
	
	private		static final byte	HNORM = 1;
	private		static final byte	EXP = 2;
	
	private String distribParam = distrib;
	private String maskingParam = masking;
	private boolean	skip0Param = true;
	private boolean	noiterParam = false;
	
	
	// create inputs
	public final void setSecondInversionImage(float[] val) { inv2Image = val; }
	public final void setQuantitativeT1MapImage(float[] val) { t1mapImage = val; }
	public final void setT1weightedImage(float[] val) { isoImage = val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setBackgroundDistribution(String val) {  distribParam = val; }
	public final void setSkipZeroValues(boolean val) { skip0Param = val; }
	public final void setNonIterativeEstimate(boolean val) { noiterParam = val; }
	public final void setMaskingMethod(String val) { maskingParam = val; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; }
	public final String getLabel() { return "MP2RAGE Background Masking"; }
	public final String getName() { return "MP2RAGEBackgroundMasking"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for Neuroimaging, Netherlands Institute for Neuroscience, Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimate a background signal mask for a MP2RAGE dataset."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
	
	// create outputs
	public final float[] getSignalProbaImage() { return probaImage; }
	public final int[] getSignalMaskImage() { return maskImage; }
	public final float[] getMaskedT1MapImage() { return t1maskImage; }
	public final float[] getMaskedT1weightedImage() { return isomaskImage; }

	public final void execute(){
		
		// import the image data into 1D arrays
		int inv2 = 0;
		int nimg = 1;
		
		int t1map = -1, iso = -1;
		if (t1mapImage != null) { t1map = nimg; nimg++; }
		if (isoImage != null) { iso = nimg; nimg++; }
		
		float[][] image = new float[nimg][nxyz];
		image[inv2] = inv2Image;
		if (t1map>-1) image[t1map] = t1mapImage;
		if (iso>-1) image[iso] = isoImage;
		
		// main algorithm

		// assume  an exponential + outlier distributions for the inv2 image
		
		// use default parameters: 
		float maxdiff = 0.0001f;
		int itermax = 20;
		
		// min max for t1 map fixing
		float max = -1e10f;
		float min = 1e10f;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (image[inv2][xyz]>max) max = image[inv2][xyz];
			if (image[inv2][xyz]<min) min = image[inv2][xyz];
		}
		// re-normalize
		for (int xyz=0;xyz<nxyz;xyz++) {
			image[inv2][xyz] = (image[inv2][xyz]-min)/(max-min);
		}
		System.out.println("image range: "+min+", "+max);

		double mean = 0.0f;
		double den = 0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0Param || image[inv2][xyz]!=0) {
			mean += image[inv2][xyz];
			den++;
		}
		mean /= den;
		System.out.println("mean parameters: "+mean);
		
		double pi = FastMath.PI;
		
		byte model = HNORM;
		if (distribParam.equals("exponential")) model = EXP;
		
		// re-normalized probability map
		float[] proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*FastMath.exp(-image[inv2][xyz]*image[inv2][xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(FastMath.exp(-image[inv2][xyz]/mean)/mean);
			// bg vs. outlier test : proba is p(xyz is background)
			proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
		}
		// loop
		if (noiterParam) itermax = 0;
		double diff = 1.0;
		double prev;
		for (int t=0;t<itermax && diff>maxdiff;t++) {
			System.out.println("iteration "+(t+1));
			prev = mean;
			mean = 0.0;
			den = 0.0;
			for (int xyz=0;xyz<nxyz;xyz++) if (!skip0Param || image[inv2][xyz]!=0) {
				mean += proba[xyz]*image[inv2][xyz];
				den += proba[xyz];
			}
			mean /= den;
			System.out.println("mean parameters: "+mean);
			diff = Numerics.abs(prev-mean)/(0.5f*(prev+mean));
			System.out.println("diff parameters: "+diff);
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (model==HNORM)
					// half-normal model
					proba[xyz] = (float)(2.0/(mean*pi)*FastMath.exp(-image[inv2][xyz]*image[inv2][xyz]/(pi*mean*mean)));
				else if (model==EXP)
					// exponential model
					proba[xyz] = (float)(FastMath.exp(-image[inv2][xyz]/mean)/mean);
				// bg vs. outlier test : proba is p(xyz is background)
				proba[xyz] = proba[xyz]/(1.0f+proba[xyz]);
			}
		}
		
		int[] mask = new int[nxyz];
		float[] prob = new float[nxyz];
			
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (proba[xyz]<=0.5f) mask[xyz] = 1;
			else mask[xyz] = 0;
			prob[xyz] = 1.0f-proba[xyz];
		}
		proba = null;

		if (maskingParam.equals("dilated"))
			mask = ObjectMorphology.dilateObject(mask, nx, ny, nz, 2, 2, 2);
		
		probaImage = prob;
		maskImage = mask;
		
		if (t1map>-1) {
			if (maskingParam.equals("proba")) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					image[t1map][xyz] *= prob[xyz];
				}
			} else {
				for (int xyz=0;xyz<nxyz;xyz++) {
					image[t1map][xyz] *= mask[xyz];
				}
			}
			t1maskImage = image[t1map];
		}
		
		if (iso>-1) {
			if (maskingParam.equals("proba")) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					image[iso][xyz] *= prob[xyz];
				}
			} else {
				for (int xyz=0;xyz<nxyz;xyz++) {
					image[iso][xyz] *= mask[xyz];
				}
			}
			isomaskImage = image[iso];
		}
		return;
	}

}
