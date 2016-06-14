package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.random.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm estimates image variance based on neighbor variations
 *
 *	Several options are possible: histogram thresholding, random sampling, outlier estimation, etc.
 *
 *	@version    Sept 2014
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class RobustVarianceEstimation {
		
	// data buffers
	private 	float[]		image;  			// original data (3D)
	private 	int			nx,ny,nz;   		// image dimensions
	private 	float		rx,ry,rz;   		// image resolutions
	private		boolean[]	mask;
	private		float		Imin, Imax;
	//private		static 	boolean		simplified;
	//private		String[]	modes = {"profiles", "binary", "scalar", "noise"};
	
	private static final double PI2 = Math.PI/2.0;
	private static final double SQPI = Math.sqrt(Math.PI);
	private static final float INF = 1e12f;
	
	private static int[] ngbx;
    private static int[] ngby;
    private static int[] ngbz;
    private static int		connect;

	static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	/**
	 *  constructor
	 */
	public RobustVarianceEstimation(float[] img_, boolean[] msk_,
										int nx_, int ny_, int nz_, 
										float rx_, float ry_, float rz_,
										int connect_) {
		image = img_;
		mask = msk_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;		
		
		connect = connect_;
		
		// 6/18/26-neighborhood: pre-compute the index offsets
		int ox = 1;
		int oy = nx;
		int oz = nx*ny;
		
		ngbx = new int[]{ ox,  0,  0, -ox,  0,  0,  ox,  0, ox,  ox, 0, -ox, -ox,   0,  ox, -ox,   0, -ox,  ox,  ox,  ox,  ox, -ox, -ox, -ox, -ox};
		ngby = new int[]{ 0,  oy,  0,  0, -oy,  0,  oy, oy,  0, -oy, oy,  0,  oy, -oy,   0, -oy, -oy,   0,  oy, -oy, -oy,  oy,  oy, -oy, -oy,  oy};
		ngbz = new int[]{ 0,  0,  oz,  0,  0, -oz,  0,  oz, oz,  0, -oz, oz,   0,  oz, -oz,   0, -oz, -oz,  oz, -oz,  oz, -oz,  oz, -oz,  oz, -oz};
	
		// image min, max
		Imin = INF; 
		Imax = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}	
	
		return;		
	}

	final public void finalize() {
		image = null;
		System.gc();
	}
	
	// variance estimates
	public final float spatialMeanStdev() {
		
		if (debug) System.out.println("-- basic spatial variance estimate --");
		
		double var = 0.0;
		int nb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				// 6-C
				if (x<nx-1 && mask[xyz+1])	{ var += Numerics.square(image[xyz]-image[xyz+1]); nb++; }
				if (y<ny-1 && mask[xyz+nx]) { var += Numerics.square(image[xyz]-image[xyz+nx]); nb++; }
				if (z<nz-1 && mask[xyz+nx*ny])	{ var += Numerics.square(image[xyz]-image[xyz+nx*ny]); nb++; }
				// 18-C
				if (connect>6) {
					if (x>0 && y<ny-1 && mask[xyz-1+nx]) { var += Numerics.square(image[xyz]-image[xyz-1+nx]); nb++; }
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx]) { var += Numerics.square(image[xyz]-image[xyz+1+nx]); nb++; }
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz-nx+nx*ny]); nb++; }
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz+nx+nx*ny]); nb++; }
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1]) { var += Numerics.square(image[xyz]-image[xyz-nx*ny+1]); nb++; }
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1]) { var += Numerics.square(image[xyz]-image[xyz+nx*ny+1]); nb++; }
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz-1-nx+nx*ny]); nb++; }
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz-1+nx+nx*ny]); nb++; }
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz+1-nx+nx*ny]); nb++; }
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny]) { var += Numerics.square(image[xyz]-image[xyz+1+nx+nx*ny]); nb++; }
				}
			}
		}
		
		// simple mean
		float meanstd = (float)FastMath.sqrt(var/(nb-1.0));
		
		return meanstd;
	}

	// variance estimates
	public final float histogramKIThresholdStdev(int nbins) {
		
		if (debug) System.out.println("-- basic histogram-based variance estimate --");
		
		int ngb = 1; // 6C neighborhood
		if (connect>6) ngb = 2;
		else if (connect>18) ngb = 3;
		
		System.out.println("build local difference histogram");
		
		Histogram hist = new Histogram(nbins);
		hist.buildFromDifferences(image, mask, ngb, nx, ny, nz);
		
		int threshold = hist.computeKIThresholdExhaustive();
		// get noise model from threshold
		float stdev = hist.getPartialMean(0,threshold);
		
		return stdev;
	}
	
	// variance estimates
	public final float medianHistogramStdev(int nbins) {
		
		if (debug) System.out.println("-- basic histogram-based variance estimate --");
		
		int ngb = 1; // 6C neighborhood
		if (connect>6) ngb = 2;
		else if (connect>18) ngb = 3;
		
		System.out.println("build local difference histogram");
		
		Histogram hist = new Histogram(nbins);
		hist.buildFromDifferences(image, mask, ngb, nx, ny, nz);
		// simply get the median
		float stdev = hist.percentage(0.50f);
		
		return stdev;
	}
	
	// variance estimates
	public final float medianSamplingStdev(int nsample, int niter, float mindist) {
		return medianSamplingStdev(nsample, niter, mindist, 19.15); 
	}
	
	// variance estimates
	public final float medianSamplingStdev(int nsample, int niter, float mindist, double halfvol) {
		
		double[] range = new double[niter];
		Well19937c rand = new Well19937c();
		
		// sample uniformly location and difference direction
		float dist = 1e9f;
		double medrng = 0.0;
		int maxtrial = (nx*ny+ny*nz+nz*nx)*connect;
		//float mindist = 0.01f;
		boolean stop = false;
		int nstop = 0;
		for (int t=0;t<niter && !stop;t++) {
			// build a sample array
			double[] sample = new double[nsample];
			int count=0;
			int trial=0;
			while (count<nsample && trial<maxtrial) {
				trial++;
				// generate pseudo-random locations
				int x = 1+rand.nextInt(nx-1);		
				int y = 1+rand.nextInt(ny-1);		
				int z = 1+rand.nextInt(nz-1);	
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					// pseudo-random direction
					int d = rand.nextInt(connect);
					if (mask[xyz+ngbx[d]+ngby[d]+ngbz[d]]) {
						// centered difference (shape is closer to Gaussian distribution)
						sample[count] = image[xyz]-image[xyz+ngbx[d]+ngby[d]+ngbz[d]];
						count++;
					}
				}
			}
			if (trial>=maxtrial) System.out.println("too many failed trials!! (succeeded: "+count+")");
			
			// estimate mean, variance robustly
			Percentile measure = new Percentile();
			measure.setData(sample);
			// compute the range mean +/- 1/2 stdev
			// 0.5 +/- 0.1915
			range[t] = measure.evaluate(50.0 + halfvol) - measure.evaluate(50.0 - halfvol);
			
			// check for convergence with the mean (rather than median that has non-zero proba to be zero randomly)??
			if (t>FastMath.sqrt(niter)) {
				
				// find the median of the range itself
				Percentile median = new Percentile();
				double newrng = median.evaluate(range, 0, t+1, 50.0);
				
				dist = (float)Numerics.abs(medrng - newrng);
				if (dist<mindist*(Imax-Imin)) {
					nstop++;
					if (nstop>=5) stop = true;
				} else {
					nstop=0;
				}
				medrng = newrng;
			}
			//System.out.println((t+1)+". sample (median, range): "+median[t]+", "+range[t]+" --> estimate: "+(meanmed/nsample)+", "+(meanrng/nsample));
		}
		System.out.println("median range: "+medrng);
		
		return (float)medrng;
	}
	
	// variance estimates
	public final float recursiveMedianSamplingStdev(int nsample, int niter, float mindist) {

		// step 1: assume no outliers
		double medrng = medianSamplingStdev(nsample, niter, mindist);
		
		// estimate the outliers
		double outratio = 0.0;
		double sum = 0.0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				for (int d=0;d<connect;d++) {
					if (mask[xyz+ngbx[d]+ngby[d]+ngbz[d]]) {
						float diff = image[xyz]-image[xyz+ngbx[d]+ngby[d]+ngbz[d]];
						float proba = (float)(FastMath.exp( -0.5*Numerics.square(diff/medrng))
											  /FastMath.sqrt(2.0*FastMath.PI*medrng*medrng) );
						float ratio = 1.0f/(Imax-Imin)/(1.0f/(Imax-Imin)+proba);
						outratio += ratio;
						sum++;
					}
				}
			}
		}
		outratio /= sum;
		System.out.println("outlier ratio: "+outratio);
		
		for (int n=0;n<50;n++) {
			double prevratio = outratio;
			double halfvol = (0.1915*(1-outratio) + outratio*medrng/(Imax-Imin))*100.0;
			
			System.out.println("half volume under the [mu-sigma,mu+sigma] curves (%): "+halfvol);
		
			medrng = medianSamplingStdev(nsample, niter, mindist, halfvol);
		
			outratio = 0.0f;
			sum = 0.0f;
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					for (int d=0;d<connect;d++) {
						if (mask[xyz+ngbx[d]+ngby[d]+ngbz[d]]) {
							float diff = image[xyz]-image[xyz+ngbx[d]+ngby[d]+ngbz[d]];
							float proba = (float)(FastMath.exp( -0.5*Numerics.square(diff/medrng))
												  /FastMath.sqrt(2.0*FastMath.PI*medrng*medrng) );
							float ratio = 1.0f/(Imax-Imin)/(1.0f/(Imax-Imin)+proba);
							outratio += ratio;
							sum++;
						}
					}
				}
				outratio /= sum;
				
				System.out.println("outlier ratio: "+outratio);
				
				if (Numerics.abs(prevratio-outratio)<0.0001) n = 10000;
			}
		}

		return (float)medrng;
	}
}
