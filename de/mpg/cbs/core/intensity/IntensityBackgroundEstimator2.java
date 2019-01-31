package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensityBackgroundEstimator2 {

	// input parameters
	private		float[] 	image;
	private		int			nx, ny, nz, nxyz;
	private		float		ratio = 0.001f;
	private		String		distribution = "exponential";
	private 	boolean		skip0 = true;
	private 	boolean		iterate = true;
	private     int         dilate = 0;
	private     float       threshold = 0.5f;
	
	// output parameters
	private		byte[] mask;
	private		float[] masked;
	private		float[] proba;
	
	public		static final String[]	distributions = {"exponential","half-normal","exp+log-normal","half+log-normal"};
	private		static final byte	HNORM = 101;
	private		static final byte	EXP = 102;
	private		static final byte	UNI = 21;
	private		static final byte	LOGN = 22;
	
	
	
	// set inputs
	public final void setInputImage(float[] in) { image = in; }
	public final void setRobustMinMaxThresholding(float r) { ratio = r; }
	public final void setBackgroundDistribution(String dist) { distribution = dist; }
	public final void setSkipZeroValues(boolean sk) { skip0 = sk; }
	public final void setIterative(boolean it) { iterate = it; }
	public final void setDilateMask(int dl) {dilate = dl; }
	public final void setMaskThreshold(float th) {threshold = th; }
	
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

		int ndata=0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[xyz]!=0) ndata++;
		
		double[] data = new double[ndata];
		int n=0;
		double max=0.0;
		for (int xyz=0;xyz<nxyz;xyz++) if (!skip0 || image[xyz]!=0) {
			data[n] = image[xyz];
			if (data[n]>max) max = data[n];
			n++;
		}
		Percentile measure = new Percentile();
		measure.setData(data, 0, ndata);
		max = measure.evaluate(100.0*(1.0-ratio));
		
		BasicInfo.displayMessage("image robust max: "+max+"\n");
		
		double mean = 0.0;		
		double pi = Math.PI;
		
		byte model = EXP;
		if (distribution.startsWith("half")) model = HNORM;
		byte fgmodel = UNI;
		if (distribution.endsWith("log-normal")) fgmodel = LOGN;
		
		if (model==EXP) mean = ImageStatistics.robustExponentialFit(data, iterate, ndata);
		else if (model==HNORM) mean = ImageStatistics.robustHalfGaussianFit(data, iterate, ndata);
		
		BasicInfo.displayMessage("background mean: "+mean+"\n");
		
		double fmean = 0.0;
		double fdev = 0.0;
		if (fgmodel==LOGN) {
			// model the filter response as something more interesting, e.g. log-normal (removing the bg samples)
			int nb=0;
			for (int xyz=0;xyz<nxyz;xyz++) if (image[xyz]>0) {
				// only count positive values
				nb++;
			}
			double[] response = new double[nb];
			n=0;
			for (int xyz=0;xyz<nxyz;xyz++) if (image[xyz]>0) {
				response[n] = image[xyz];
				n++;
			}
			double[] weights = new double[nb];
			for (int b=0;b<nb;b++) { 
				if (model==EXP) weights[b] = (1.0-FastMath.exp( -response[b]/mean));
				else if (model==HNORM)  weights[b] = (1.0-Math.exp(-response[b]*response[b]/(pi*mean*mean)));
				else weights[b] = 1.0;
				response[b] = FastMath.log(response[b]);
			}
		
			fmean = ImageStatistics.weightedPercentile(response,weights,50.0,nb);
		
			// stdev: 50% +/- 1/2*erf(1/sqrt(2)) (~0.341344746..)
			double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
			fdev = 0.5*(ImageStatistics.weightedPercentile(response,weights,50.0+dev,nb) - ImageStatistics.weightedPercentile(response,weights,50.0-dev,nb));
		
			BasicInfo.displayMessage("Log-normal parameter estimates: mean = "+FastMath.exp(fmean)+", stdev = "+FastMath.exp(fdev)+",\n");
		}
		
		// re-normalized probability map
		proba = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (model==HNORM)
				// half-normal model
				proba[xyz] = (float)(2.0/(mean*pi)*Math.exp(-image[xyz]*image[xyz]/(pi*mean*mean)));
			else if (model==EXP)
				// exponential model
				proba[xyz] = (float)(Math.exp(-image[xyz]/mean)/mean);
				
				
			// bg vs. outlier test : proba is p(xyz is background)
			if (fgmodel==UNI) {
			    float uni = 1.0f;
			    if (model==HNORM) uni = (float)((max/(max-mean*pi/2.0))*(1.0 - Math.exp(-image[xyz]*image[xyz]/(pi*mean*mean))));
			    else if (model==EXP) uni = (float)((max/(max-mean))*(1.0 - Math.exp(-image[xyz]/mean)));
			    proba[xyz] = (float)max*proba[xyz]/(uni+(float)max*proba[xyz]);
			    
			} else if (fgmodel==LOGN) {
				double plg = FastMath.exp(-Numerics.square(FastMath.log(image[xyz])-fmean)/(2.0*fdev*fdev))/FastMath.sqrt(2.0*FastMath.PI*fdev*fdev);
				proba[xyz] = (float)(proba[xyz]/(plg+proba[xyz]));
			}
		}
		mask = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (proba[xyz]<=1.0-threshold) mask[xyz] = 1;
			else mask[xyz] = 0;
		}
		if (dilate>0) {
		    if (nz>1) mask = Morphology.dilateObject(mask, nx,ny,nz, dilate,dilate,dilate);
		    else mask = Morphology.dilateObject(mask, nx,ny,nz, dilate,dilate,0);
		} else if (dilate<0) {
			if (nz>1) mask = Morphology.erodeObject(mask, nx,ny,nz, -dilate,-dilate,-dilate);
			else mask = Morphology.erodeObject(mask, nx,ny,nz, -dilate,-dilate,0);
		}

		masked = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz]>0) masked[xyz] = image[xyz];
			else masked[xyz] = 0.0f;
			proba[xyz] = 1.0f-proba[xyz];
		}
		
		return;
	}


}
