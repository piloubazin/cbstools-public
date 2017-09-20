package de.mpg.cbs.structures;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class makes various histogram manipulations
 *	
 *	@version    May 2012
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Histogram {
	
	double[] hist;
	float	min, max;
	int 	bins;
	
	private static final float INF = 1e9f;

	private static final boolean debug = false;
	private static final boolean verbose = false;
	
	/**
     *    1D histogram building
     */
    public Histogram(int bins) {
    	this.bins = bins;
    	this.min = 0;
    	this.max = 1;
    	
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;			
	}

	/**
     *    1D histogram building
     */
    public Histogram(float[][][] image, float min, float max, int bins, int nx, int ny, int nz) {
    	this.bins = bins;
    	this.min = min;
    	this.max = max;
    	
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
			
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			// compute histogram within min, max (rest is ignored)
			if (  (image[x][y][z] >= min )
					&& (image[x][y][z] <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
			
			for (int n=1;n<bins;n++) {
				if (  (image[x][y][z] >  min + (float)n/(float)bins*(max-min) )
					&& (image[x][y][z] <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
			}
		}
	}

	/**
     *    1D histogram building
     */
    public Histogram(float[][][] image, int bins, int nx, int ny, int nz) {
    	this.bins = bins;
    	
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
			
		min = INF;
		max = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z] > max) max = image[x][y][z];
			if (image[x][y][z] < min) min = image[x][y][z];
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			// compute histogram within min, max (rest is ignored)
			if (  (image[x][y][z] >= min )
					&& (image[x][y][z] <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
			
			for (int n=1;n<bins;n++) {
				if (  (image[x][y][z] >  min + (float)n/(float)bins*(max-min) )
					&& (image[x][y][z] <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
			}
		}
	}

	/**
     *    1D histogram building
     */
    public Histogram(float[][][] image, boolean[][][] mask, int bins, int nx, int ny, int nz) {
    	this.bins = bins;
    	
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
			
		min = INF;
		max = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			if (image[x][y][z] > max) max = image[x][y][z];
			if (image[x][y][z] < min) min = image[x][y][z];
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			// compute histogram within min, max (rest is ignored)
			if (  (image[x][y][z] >= min )
					&& (image[x][y][z] <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
			
			for (int n=1;n<bins;n++) {
				if (  (image[x][y][z] >  min + (float)n/(float)bins*(max-min) )
					&& (image[x][y][z] <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
			}
		}
	}

	/**
     *    1D histogram building
     */
    public Histogram(float[] data, int bins, int size) {
    	this.bins = bins;
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
		
		min = INF;
		max = -INF;
		for (int s=0;s<size;s++) {
			if (data[s] > max) max = data[s];
			if (data[s] < min) min = data[s];
		}
		
		for (int s=0;s<size;s++) {
			// compute histogram within min, max (rest is ignored)
			if (  (data[s] >= min )
					&& (data[s] <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
			
			for (int n=1;n<bins;n++) {
				if (  (data[s] >  min + (float)n/(float)bins*(max-min) )
					&& (data[s] <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
			}
		}
	}

	/**
     *    1D histogram building
     */
    public Histogram(double[] data, int bins, int size) {
    	this.bins = bins;
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
		
		min = INF;
		max = -INF;
		for (int s=0;s<size;s++) {
			if (data[s] > max) max = (float)data[s];
			if (data[s] < min) min = (float)data[s];
		}
		
		for (int s=0;s<size;s++) {
			// compute histogram within min, max (rest is ignored)
			if (  (data[s] >= min )
					&& (data[s] <= min + 1.0f/(double)bins*(max-min) ) ) hist[0]++;
			
			for (int n=1;n<bins;n++) {
				if (  (data[s] >  min + (double)n/(double)bins*(max-min) )
					&& (data[s] <= min + (double)(n+1)/(double)bins*(max-min) ) ) hist[n]++;	
			}
		}
	}

	/**
     *    1D histogram building with unequal weights
     */
    public Histogram(float[] data, int[] weight, int bins, int size) {
    	this.bins = bins;
		hist = new double[bins];
		
		for (int n=0;n<bins;n++) hist[n] = 0;
		
		min = INF;
		max = -INF;
		for (int s=0;s<size;s++) {
			if (data[s] > max) max = data[s];
			if (data[s] < min) min = data[s];
		}
		
		for (int s=0;s<size;s++) {
			// compute histogram within min, max (rest is ignored)
			if (  (data[s] >= min )
					&& (data[s] <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]+=weight[s];
			
			for (int n=1;n<bins;n++) {
				if (  (data[s] >  min + (float)n/(float)bins*(max-min) )
					&& (data[s] <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]+=weight[s];	
			}
		}
	}
	/**
     *    1D histogram building
     */
    public void buildFromDifferences(float[][][] image, boolean[][][] mask, int ngb, int nx, int ny, int nz) {
    	
    	// ngb: 1 for 6C, 2 for 18C, 3 for 26C
    	
    	// estimate max difference
		min = 0.0f;
		max = -INF;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l>0 && i*i+j*j+l*l<=ngb && mask[x+i][y+j][z+l]) {
					float val = Numerics.abs(image[x][y][z]-image[x+i][y+j][z+l]);
					if (val > max) max = val;
				}
			}
		}
		if (debug) System.out.println("max: "+max);
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l>0 && i*i+j*j+l*l<=ngb && mask[x+i][y+j][z+l]) {
					float val = Numerics.abs(image[x][y][z]-image[x+i][y+j][z+l]);
					// compute histogram within min, max (rest is ignored)
					if (  (val >= min )
							&& (val <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
					
					for (int n=1;n<bins;n++) {
						if (  (val >  min + (float)n/(float)bins*(max-min) )
							&& (val <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
					}
				}
			}
		}
	}
	/**
     *    1D histogram building
     */
    public void buildFromDifferences(float[] image, boolean[] mask, int ngb, int nx, int ny, int nz) {
    	
    	// ngb: 1 for 6C, 2 for 18C, 3 for 26C
    	
    	// estimate max difference
		min = 0.0f;
		max = -INF;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz]) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					int ijl = i+nx*j+nx*ny*l;
					if (i*i+j*j+l*l>0 && i*i+j*j+l*l<=ngb && mask[xyz+ijl]) {
						float val = Numerics.abs(image[xyz]-image[xyz+ijl]);
						if (val > max) max = val;
					}
				}
			}
		}
		if (debug) System.out.println("max: "+max);
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz]) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					int ijl = i+nx*j+nx*ny*l;
					if (i*i+j*j+l*l>0 && i*i+j*j+l*l<=ngb && mask[xyz+ijl]) {
						float val = Numerics.abs(image[xyz]-image[xyz+ijl]);
						// compute histogram within min, max (rest is ignored)
						if (  (val >= min )
								&& (val <= min + 1.0f/(float)bins*(max-min) ) ) hist[0]++;
						
						for (int n=1;n<bins;n++) {
							if (  (val >  min + (float)n/(float)bins*(max-min) )
								&& (val <= min + (float)(n+1)/(float)bins*(max-min) ) ) hist[n]++;	
						}
					}
				}
			}
		}
	}

	public final float percentage(float percent) {
		if (percent<=0) return min;
		else if (percent>=1) return max;
		
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
		
		double partial = 0.0f;
		float value = 0.0f;
		for (int n=0;n<bins && partial<percent*total;n++) {
			partial += hist[n];
			if (partial>=percent*total) {
				value = (float)( (min + (float)n/(float)bins*(max-min))*(partial - percent*total)/hist[n]
								+(min + (float)(n+1)/(float)bins*(max-min))*(percent*total - partial + hist[n])/hist[n] );
			}
		}
		return value;
	}
	
	public final float mean() {
		
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
		
		double mean = 0.0f;
		for (int n=0;n<bins;n++) {
			mean += hist[n]*(min + (float)(n+0.5f)/(float)bins*(max-min));
		}
		return (float)(mean/total);
	}
	
	public final float stdev() {
		
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
		
		double mean = 0.0f;
		for (int n=0;n<bins;n++) {
			mean += hist[n]*(min + (float)(n+0.5f)/(float)bins*(max-min));
		}
		mean /= total;
		
		double var = 0.0f;
		for (int n=0;n<bins;n++) {
			var += hist[n]*Numerics.square( ( (min + (float)(n+0.5f)/(float)bins*(max-min)) - mean) );
		}
		
		return (float)Math.sqrt(var/(total-1));
	}
	
	public final float argmax() {
		
		int nmax = 0;
		for (int n=0;n<bins;n++) if (hist[n]>hist[nmax]) nmax = n;
		
		return (min + (float)(nmax+0.5f)/(float)bins*(max-min));
	}
	
	public final float rayleighParameter() {
		
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
		
		double sum = 0.0f;
		for (int n=0;n<bins;n++) {
			sum += hist[n]*Numerics.square( (min + (float)(n+0.5f)/(float)bins*(max-min)) );
		}
		return (float)Math.sqrt( sum/(2.0*total) );
	}
	
	public final float min() { return min; }
	
	public final float max() { return max; }

	public final void normalizeAsPdf() {
		
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
		
		for (int n=0;n<bins;n++) hist[n] /= total;
		
		return;
	}
	
	public final void normalizeToMaximum() {
		
		double hmax = 0.0f;
		for (int n=0;n<bins;n++) if (hist[n]>hmax) hmax = hist[n];
		
		for (int n=0;n<bins;n++) hist[n] /= hmax;
		
		return;
	}
	
	public final void normalizeToNonZeroMaximum() {
		
		double hmax = 0.0f;
		for (int n=1;n<bins;n++) if (hist[n]>hmax) hmax = hist[n];
		
		hist[0] = Numerics.min(hist[0]/hmax, 1.0);
		for (int n=1;n<bins;n++) hist[n] /= hmax;
		
		return;
	}
	
	public final void inverseCumulativeHistogram() {
		// always start at zero?
		double total = 0.0;
		double[] cumul = new double[bins];
		for (int n=0;n<bins;n++) {
			cumul[n] = total;
			total += hist[n];
		}
		for (int n=0;n<bins;n++) hist[n] = 1.0 - cumul[n]/total;
	}
	
	public final float getHistogramCount(float value) {
		int level = Numerics.bounded(Numerics.round(bins*(value-min)/(max-min)),0, bins-1);
		
		return (float)hist[level];
	}
	
	public final String printHistogram() {
		String output = "Histogram (min: "+min+", max: "+max+", nbins: "+bins+"):\n";
		for (int n=0;n<bins;n++) output += "	"+hist[n];
		return output;
	}
	
	public final byte[][] plotLogHistogram() {
		byte[][] histimg = new byte[bins][bins];
		// 1. find max
		double hmax = 0.0;
		for (int n=0;n<bins;n++) if (hist[n]>hmax) hmax = hist[n];
		hmax = FastMath.log(1+hmax);
		// 2. plot log histogram
		for (int n=0;n<bins;n++) {
			int level = Numerics.round(bins*FastMath.log(1+hist[n])/hmax);
			for (int m=bins-1;m>bins-1-level;m--) histimg[n][m] = 1;
		}
		return histimg;
	}
	
	public final byte[][] plotLogHistogram(int threshold) {
		byte[][] histimg = new byte[bins][bins];
		// 1. find max
		double hmax = 0.0;
		for (int n=0;n<bins;n++) if (hist[n]>hmax) hmax = hist[n];
		hmax = FastMath.log(1+hmax);
		// 2. plot log histogram
		for (int n=0;n<bins;n++) {
			byte val = 1;
			if (n>=threshold) val = 2;
			int level = Numerics.round(bins*FastMath.log(1+hist[n])/hmax);
			for (int m=bins-1;m>bins-1-level;m--) histimg[n][m] = val;
		}
		return histimg;
	}
	
	public final float fitHalfGaussian() {
		double sqrt = FastMath.sqrt(Math.PI/2.0)*mean();
		
		normalizeAsPdf();
		double res = 0.0f;
		for (int n=0;n<bins;n++) {
			double x = (min + (float)(n+0.5f)/(float)bins*(max-min));
			res += Numerics.abs(hist[n]-1.0/FastMath.sqrt(Math.PI/2.0)/sqrt*FastMath.exp(- Numerics.square(x/2.0/sqrt)));
		}
		return (float)res;
	}
	
	public final float fitBeta() {
		double expectation = (mean()-min)/(max-min);
		double variance = Numerics.square(stdev()/(max-min));
		
		double alpha = expectation*(expectation*(1.0-expectation)/variance-1.0);
		double beta = (1.0-expectation)*(expectation*(1.0-expectation)/variance-1.0);
		
		normalizeAsPdf();
		double res = 0.0f;
		for (int n=0;n<bins;n++) {
			double x = (float)(n+0.5f)/(float)bins;
			//res += Numerics.abs(hist[n]-1.0/FastMath.sqrt(Math.PI/2.0)/sqrt*FastMath.exp(- Numerics.square(x/2.0/sqrt)));
		}
		return (float)res;
	}
	
	/** histogram thresholding with the Kittler Illingworth algorithm */
	public final int computeKIThresholdExhaustive() {
		
		if (verbose) System.out.println("KI Thresholding values: ");
		
		// need normalized quantities
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
				
		int best = -1;
		double Jbest = 1e10;
		double zero = 0.0001/total;
		for (int t=2;t<bins-1;t++) { // we skip the first and last bins
			// compute prior, mean and stdev
			
			double P1 = 0.0, P2 = 0.0;
			for (int n=0;n<t;n++) P1 += hist[n]/total;
			for (int n=t;n<bins;n++) P2 += hist[n]/total;
			
			double mu1 = 0.0, mu2 = 0.0;
			for (int n=0;n<t;n++) mu1 += hist[n]/total*n;
			if (P1>0) mu1 /= P1;
			for (int n=t;n<bins;n++) mu2 += hist[n]/total*n;
			if (P2>0) mu2 /= P2;
			
			double sig1 = 0.0, sig2 = 0.0;
			for (int n=0;n<t;n++) sig1 += hist[n]/total*(n-mu1)*(n-mu1);
			if (P1>0) sig1 /= P1;
			for (int n=t;n<bins;n++) sig2 += hist[n]/total*(n-mu2)*(n-mu2);
			if (P2>0) sig2 /= P2;
			
			double Jt = 1.0 + 2.0*(P1*FastMath.log(Numerics.max(sig1,zero))+P2*FastMath.log(Numerics.max(sig2,zero))) 
							 - 2.0*(P1*FastMath.log(Numerics.max(P1,zero))+P2*FastMath.log(Numerics.max(P2,zero)));
			if (Jt<Jbest) {
				Jbest = Jt;
				best = t;
			}
			if (debug) System.out.print(Jt+", ");
		}
		
		if (verbose) System.out.println("\n -> best: "+Jbest+" (at "+best+")");
		
		/*
		if (best!=-1) return (min + (float)(best+0.5f)/(float)bins*(max-min));
		else return -1.0f;
		*/
		return best;
	}
	
	/** histogram thresholding with the Kittler Illingworth algorithm */
	public final int computeExpNormalKIThresholdExhaustive() {
		
		if (verbose) System.out.println("KI Thresholding values (exponential / normal distributions): ");
		
		// need normalized quantities
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
				
		int best = -1;
		double Jbest = 1e10;
		double zero = 0.0001/total;
		for (int t=2;t<bins-1;t++) { // we skip the first and last bins
			// compute prior, mean and stdev
			
			double P1 = 0.0, P2 = 0.0;
			for (int n=0;n<t;n++) P1 += hist[n]/total;
			for (int n=t;n<bins;n++) P2 += hist[n]/total;
			
			double mu1 = 0.0, mu2 = 0.0;
			for (int n=0;n<t;n++) mu1 += hist[n]/total*n;
			if (P1>0) mu1 /= P1;
			for (int n=t;n<bins;n++) mu2 += hist[n]/total*n;
			if (P2>0) mu2 /= P2;
			
			double sig1 = 0.0, sig2 = 0.0;
			for (int n=0;n<t;n++) sig1 += hist[n]/total*(n-mu1)*(n-mu1);
			if (P1>0) sig1 /= P1;
			for (int n=t;n<bins;n++) sig2 += hist[n]/total*(n-mu2)*(n-mu2);
			if (P2>0) sig2 /= P2;
			
			double Jt = 1.0 + 2.0*(P1*FastMath.log(Numerics.max(0.5*(mu1+sig1),zero))+P2*FastMath.log(Numerics.max(sig2,zero))) 
							 - 2.0*(P1*FastMath.log(Numerics.max(P1,zero))+P2*FastMath.log(Numerics.max(P2,zero)));
			if (Jt<Jbest) {
				Jbest = Jt;
				best = t;
			}
			if (debug) System.out.print(Jt+", ");
		}
		
		if (verbose) System.out.println("\n -> best: "+Jbest+" (at "+best+")");
		
		/*
		if (best!=-1) return (min + (float)(best+0.5f)/(float)bins*(max-min));
		else return -1.0f;
		*/
		return best;
	}
	
	/** histogram thresholding with the Kittler Illingworth algorithm */
	public final int computeNonnegKIThresholdExhaustive() {
		
		if (verbose) System.out.println("Non-negative KI Thresholding values: ");
		
		// need normalized quantities
		double total = 0.0f;
		for (int n=0;n<bins;n++) total += hist[n];
				
		int best = -1;
		double Jbest = 1e10;
		double zero = 0.0001/total;
		for (int t=2;t<bins-1;t++) { // we skip the first and last bins
			// compute prior, mean and stdev
			
			double P1 = 0.0, P2 = 0.0;
			for (int n=0;n<t;n++) P1 += hist[n]/total;
			for (int n=t;n<bins;n++) P2 += hist[n]/total;
			
			double mu1 = 0.0, mu2 = 0.0;
			// mu1 is always zero by hypothesis (half-normal vs. normal)
			/*
			for (int n=0;n<t;n++) mu1 += hist[n]/total*n;
			if (P1>0) mu1 /= P1;
			*/
			for (int n=t;n<bins;n++) mu2 += hist[n]/total*n;
			if (P2>0) mu2 /= P2;
			
			double sig1 = 0.0, sig2 = 0.0;
			for (int n=0;n<t;n++) sig1 += hist[n]/total*(n-mu1)*(n-mu1);
			if (P1>0) sig1 /= P1;
			for (int n=t;n<bins;n++) sig2 += hist[n]/total*(n-mu2)*(n-mu2);
			if (P2>0) sig2 /= P2;
			
			double Jt = 1.0 + 2.0*(P1*FastMath.log(Numerics.max(sig1,zero))+P2*FastMath.log(Numerics.max(sig2,zero))) 
							 - 2.0*(P1*FastMath.log(Numerics.max(P1,zero))+P2*FastMath.log(Numerics.max(P2,zero)));
			if (Jt<Jbest) {
				Jbest = Jt;
				best = t;
			}
			if (debug) System.out.print(Jt+", ");
		}
		
		if (verbose) System.out.println("\n -> best: "+Jbest+" (at "+best+")");
		
		/*
		if (best!=-1) return (min + (float)(best+0.5f)/(float)bins*(max-min));
		else return -1.0f;
		*/
		return best;
	}
	
	public final float getPartialMean(int from, int to) {
		// always start at zero?
		double mean = 0.0;
		double count = 0.0;
		for (int n=from;n<to;n++) {
			mean += hist[n]*(min + (float)(n+0.5f)/(float)bins*(max-min));
			count += hist[n];
		}
		if (count>0) mean /= count;
		return (float)mean;
	}

	public final float getIntensityAt(int bin) {
		if (bin<0) return min;
		else if (bin>=bins) return max;
		else return (min + (float)(bin+0.5f)/(float)bins*(max-min));
	}
}
