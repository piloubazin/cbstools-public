package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm builds tools for analyzing the intrinsic scale
 *	of (large) images based on various models
 *
 *	@version    September 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class RandomWalkSegmentation {
		
	// data buffers
	private 	float[]		image;  			// original images
	private 	int			nx,ny,nz;   		// image dimensions
	private 	float		rx,ry,rz;   		// image resolutions
	private		float		Imin, Imax;

	
	private		byte[]		labels;
	private		int		nlb;
	private		float		scale;				// difference scale used here
	
	private		float[][]	weight;
	private		float[][]	intens;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
			
	public	static	final	byte	pX = 0;
	public	static	final	byte	pY = 1;
	public	static	final	byte	pZ = 2;
	public	static	final	byte	mX = 3;
	public	static	final	byte	mY = 4;
	public	static	final	byte	mZ = 5;
	public	static	final	byte	S = 6;

	/**
	 *  constructor
	 */
	public RandomWalkSegmentation(float[] img_, byte[] lbs_,
										int nx_, int ny_, int nz_,
										float rx_, float ry_, float rz_,
										int nlb_, float scale_) {
		image = img_;
		labels = lbs_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		nlb = nlb_;
		
		scale = scale_;
		
		Imin = ImageStatistics.robustMinimum(image, 0.01f, 2, nx, ny, nz, 4);
		Imax = ImageStatistics.robustMaximum(image, 0.01f, 2, nx, ny, nz, 4);

		if (debug) System.out.println("image: ["+Imin+", "+Imax+"]");
	}

	final public void finalize() {
		image = null;
		labels = null;
		System.gc();
	}
	
	public final float[][] getIntensity() { return intens; }
	
	public final void computeEdgeWeights() {
		weight = new float[7][nx*ny*nz];
		
		float r0 = Numerics.min(rx, ry, rz);
		
		for (int d=0;d<6;d++) for (int xyz=0;xyz<nx*ny*nz;xyz++) weight[d][xyz] = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (x<nx-1)	weight[pX][xyz] = 1.0f/(1.0f + Numerics.square( rx/r0*(image[xyz]-image[xyz+1])/(scale*(Imax-Imin)) ) );
			if (y<ny-1)	weight[pY][xyz] = 1.0f/(1.0f + Numerics.square( ry/r0*(image[xyz]-image[xyz+nx])/(scale*(Imax-Imin)) ) );
			if (z<nz-1)	weight[pZ][xyz] = 1.0f/(1.0f + Numerics.square( rz/r0*(image[xyz]-image[xyz+nx*ny])/(scale*(Imax-Imin)) ) );
			if (x>0)	weight[mX][xyz] = 1.0f/(1.0f + Numerics.square( rx/r0*(image[xyz]-image[xyz-1])/(scale*(Imax-Imin)) ) );
			if (y>0)	weight[mY][xyz] = 1.0f/(1.0f + Numerics.square( ry/r0*(image[xyz]-image[xyz-nx])/(scale*(Imax-Imin)) ) );
			if (z>0)	weight[mZ][xyz] = 1.0f/(1.0f + Numerics.square( rz/r0*(image[xyz]-image[xyz-nx*ny])/(scale*(Imax-Imin)) ) );
			weight[S][xyz] = weight[pX][xyz] + weight[pY][xyz] + weight[pZ][xyz]
							+ weight[mX][xyz] + weight[mY][xyz] + weight[mZ][xyz];
		}
	}

	public final void compute2DEdgeWeights() {
		weight = new float[7][nx*ny*nz];
		
		for (int d=0;d<6;d++) for (int xyz=0;xyz<nx*ny*nz;xyz++) weight[d][xyz] = 0.0f;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			weight[pX][xyz] = 1.0f/(1.0f + Numerics.square( (image[xyz]-image[xyz+1])/(scale*(Imax-Imin)) ) );
			weight[pY][xyz] = 1.0f/(1.0f + Numerics.square( (image[xyz]-image[xyz+nx])/(scale*(Imax-Imin)) ) );
			weight[mX][xyz] = 1.0f/(1.0f + Numerics.square( (image[xyz]-image[xyz-1])/(scale*(Imax-Imin)) ) );
			weight[mY][xyz] = 1.0f/(1.0f + Numerics.square( (image[xyz]-image[xyz-nx])/(scale*(Imax-Imin)) ) );
			weight[S][xyz] = weight[pX][xyz] + weight[pY][xyz]
							+ weight[mX][xyz] + weight[mY][xyz];
		}
	}

	public final float[] exportMinWeight() {
		float[] minw = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
		
			minw[xyz] = 1.0f;
			if(weight[mX][xyz]>0 && weight[mX][xyz]<minw[xyz]) minw[xyz] = weight[mX][xyz];
			if(weight[pX][xyz]>0 && weight[pX][xyz]<minw[xyz]) minw[xyz] = weight[pX][xyz];
			if(weight[mY][xyz]>0 && weight[mY][xyz]<minw[xyz]) minw[xyz] = weight[mY][xyz];
			if(weight[pY][xyz]>0 && weight[pY][xyz]<minw[xyz]) minw[xyz] = weight[pY][xyz];
			if(weight[mZ][xyz]>0 && weight[mZ][xyz]<minw[xyz]) minw[xyz] = weight[mZ][xyz];
			if(weight[pZ][xyz]>0 && weight[pZ][xyz]<minw[xyz]) minw[xyz] = weight[pZ][xyz];
		}
		return minw;
	}
	
	public final float[] exportMaxWeight() {
		float[] maxw = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
		
			maxw[xyz] = 0.0f;
			if(weight[mX][xyz]>0 && weight[mX][xyz]>maxw[xyz]) maxw[xyz] = weight[mX][xyz];
			if(weight[pX][xyz]>0 && weight[pX][xyz]>maxw[xyz]) maxw[xyz] = weight[pX][xyz];
			if(weight[mY][xyz]>0 && weight[mY][xyz]>maxw[xyz]) maxw[xyz] = weight[mY][xyz];
			if(weight[pY][xyz]>0 && weight[pY][xyz]>maxw[xyz]) maxw[xyz] = weight[pY][xyz];
			if(weight[mZ][xyz]>0 && weight[mZ][xyz]>maxw[xyz]) maxw[xyz] = weight[mZ][xyz];
			if(weight[pZ][xyz]>0 && weight[pZ][xyz]>maxw[xyz]) maxw[xyz] = weight[pZ][xyz];
		}
		return maxw;
	}
	
	public final void computeIntensityWeights() {
		// build histograms
		int nbins = 100;
		int sum;
		float[][] hist = new float[nlb][nbins];
		for (int lb=1;lb<=nlb;lb++) {
			sum=0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==lb) {
					int bin = Numerics.bounded(Numerics.floor(image[xyz]/(Imax-Imin)*nbins), 0, nbins-1);
					hist[lb-1][bin]++;
					sum++;
				}
			}
			for (int n=0;n<nbins;n++) hist[lb-1][n] /= sum;
		}
		// smooth them
		for (int t=0;t<nbins;t++) {
			for (int lb=0;lb<nlb;lb++) {
				for (int n=0;n<nbins-1;n++) {
					hist[lb][n] = (hist[lb][n] + 0.5f*hist[lb][n+1])/1.5f;
				}
				for (int n=nbins-1;n>0;n--) {
					hist[lb][n] = (hist[lb][n] + 0.5f*hist[lb][n-1])/1.5f;
				}
			}
		}
		/*
		// normalize to [0, 1]
		for (int lb=0;lb<nlb;lb++) {
			float max = 0;
			for (int n=0;n<nbins;n++) {
				if (hist[lb][n]>max) max = hist[lb][n];
			}
			for (int n=0;n<nbins;n++) {
				hist[lb][n] /= max;
			}
		}
		*/
		// normalize over the labels
		for (int n=0;n<nbins;n++) {
			float hsum = 0;
			for (int lb=0;lb<nlb;lb++) hsum += hist[lb][n];
			hsum = Numerics.max(hsum, 1e-6f);
			for (int lb=0;lb<nlb;lb++) {
				hist[lb][n] /= hsum;
			}
		}
		// estimate the probability for each point
		intens = new float[nlb][nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			int bin = Numerics.bounded(Numerics.floor(image[xyz]/(Imax-Imin)*nbins), 0, nbins-1);
			for (int lb=0;lb<nlb;lb++) {
				intens[lb][xyz] = hist[lb][bin];
			}
		}
	}
	
	public final float[] compute3DLabelProbability(int lb, int iter) {
		
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			p[xyz] = 0.0f;
		}
		
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
					 if (labels[xyz]==lb) p[xyz] = (1.0f 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz+1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz+nx]
												+ weight[pZ][xyz]*p[xyz+nx*ny] + weight[mZ][xyz]*p[xyz-nx*ny])
												/(1.0f+weight[S][xyz]);
				else if (labels[xyz]==0) p[xyz] = (0.5f 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz+1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz+nx]
												+ weight[pZ][xyz]*p[xyz+nx*ny] + weight[mZ][xyz]*p[xyz-nx*ny])
												/(1.0f+weight[S][xyz]);
				else 					 p[xyz] = (0.0f 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz+1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz+nx]
												+ weight[pZ][xyz]*p[xyz+nx*ny] + weight[mZ][xyz]*p[xyz-nx*ny])
												/(1.0f+weight[S][xyz]);
			}
		}
		return p;
	}
	
	public final float[] compute2DLabelProbability(int lb, int iter) {
		
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			p[xyz] = 0.0f;
		}
		
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
					 if (labels[xyz]==lb) p[xyz] = (1.0f 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz-1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz-nx])
												/(1.0f+weight[S][xyz]);
				else if (labels[xyz]==0) p[xyz] = (p[xyz] 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz-1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz-nx])
												/(1.0f+weight[S][xyz]);
				else 					 p[xyz] = (0.0f 
												+ weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz-1]
												+ weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz-nx])
												/(1.0f+weight[S][xyz]);
			}
		}
		return p;
	}
	
	public final float[] compute2DLabelGain(int lb, int iter) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		update2DLabelGain(p, labels, lb, iter);
		
		return p;
	}
		
	public final void update2DLabelGain(float[] p, byte[] lbs, int lb, int iter) {
		labels = lbs;
		
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) {
					p[xyz] = (p[xyz] + weight[pX][xyz]*p[xyz+1] + weight[mX][xyz]*p[xyz-1] 
									 + weight[pY][xyz]*p[xyz+nx] + weight[mY][xyz]*p[xyz-nx])
									/(1.0f + weight[S][xyz]);
				}
			}
		}
		return;
	}
	
	public final float[] computeFast2DLabelGain(int lb, int iter) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum;
		float ratio;
		for (int t=0;t<iter;t++) {
			ratio = (iter-1-t)/(iter-1);
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) {
					val = p[xyz];
					sum = 1.0f;
					if (p[xyz+1]*p[xyz+1]>ratio*p[xyz]*p[xyz]) {
						val += weight[pX][xyz]*p[xyz+1];
						sum += weight[pX][xyz];
					}
					if (p[xyz-1]*p[xyz-1]>ratio*p[xyz]*p[xyz]) {
						val += weight[mX][xyz]*p[xyz-1];
						sum += weight[mX][xyz];
					}
					if (p[xyz+nx]*p[xyz+nx]>ratio*p[xyz]*p[xyz]) {
						val += weight[pY][xyz]*p[xyz+nx];
						sum += weight[pY][xyz];
					}
					if (p[xyz-nx]*p[xyz-nx]>ratio*p[xyz]*p[xyz]) {
						val += weight[mY][xyz]*p[xyz-nx];
						sum += weight[mY][xyz];
					}
					p[xyz] = val/sum;
				}
			}
		}
		return p;
	}
	
	public final float[] computeFast1DLabelGain(int lb, int iter, float ratio) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum;
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[mX][xyz]*p[xyz-1])/(1.0f + ratio*weight[mX][xyz]);
			}
			for (int x=nx-2;x>=0;x--) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[pX][xyz]*p[xyz+1])/(1.0f + ratio*weight[pX][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=1;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[mY][xyz]*p[xyz-nx])/(1.0f + ratio*weight[mY][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=ny-2;y>=0;y--) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[pY][xyz]*p[xyz+nx])/(1.0f + ratio*weight[pY][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=1;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[mZ][xyz]*p[xyz-nx*ny])/(1.0f + ratio*weight[mZ][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=nz-2;z>=0;z--) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*weight[pZ][xyz]*p[xyz+nx*ny])/(1.0f + ratio*weight[pZ][xyz]);
			}
		}
		return p;
	}
	
	public final float[] computeFast2DLabelGain(int lb, int iter, float ratio) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum;
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) {
					val = p[xyz];
					sum = 1.0f;
					if (p[xyz+1]*p[xyz+1]>ratio*p[xyz]*p[xyz]) {
						val += weight[pX][xyz]*p[xyz+1];
						sum += weight[pX][xyz];
					}
					if (p[xyz-1]*p[xyz-1]>ratio*p[xyz]*p[xyz]) {
						val += weight[mX][xyz]*p[xyz-1];
						sum += weight[mX][xyz];
					}
					if (p[xyz+nx]*p[xyz+nx]>ratio*p[xyz]*p[xyz]) {
						val += weight[pY][xyz]*p[xyz+nx];
						sum += weight[pY][xyz];
					}
					if (p[xyz-nx]*p[xyz-nx]>ratio*p[xyz]*p[xyz]) {
						val += weight[mY][xyz]*p[xyz-nx];
						sum += weight[mY][xyz];
					}
					p[xyz] = val/sum;
				}
			}
		}
		return p;
	}
	
	public final float[] computeFastIntens2DLabelGain(int lb, int iter, float ratio) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum, gI;
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) {
					val = p[xyz];
					sum = 1.0f;
					if (p[xyz+1]*p[xyz+1]>ratio*p[xyz]*p[xyz]) {
						val += intens[lb-1][xyz]*weight[pX][xyz]*p[xyz+1];
						sum += intens[lb-1][xyz]*weight[pX][xyz];
					}
					if (p[xyz-1]*p[xyz-1]>ratio*p[xyz]*p[xyz]) {
						val += intens[lb-1][xyz]*weight[mX][xyz]*p[xyz-1];
						sum += intens[lb-1][xyz]*weight[mX][xyz];
					}
					if (p[xyz+nx]*p[xyz+nx]>ratio*p[xyz]*p[xyz]) {
						val += intens[lb-1][xyz]*weight[pY][xyz]*p[xyz+nx];
						sum += intens[lb-1][xyz]*weight[pY][xyz];
					}
					if (p[xyz-nx]*p[xyz-nx]>ratio*p[xyz]*p[xyz]) {
						val += intens[lb-1][xyz]*weight[mY][xyz]*p[xyz-nx];
						sum += intens[lb-1][xyz]*weight[mY][xyz];
					}
					p[xyz] = val/sum;
				}
			}
		}
		return p;
	}
	
	public final float[] computeFastIntensLabelGain(int lb, int iter, float ratio) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum, gI;
		for (int t=0;t<iter;t++) {
			// do the borders as well
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) {
					val = p[xyz];
					sum = 1.0f;
					if (x<nx-1 && p[xyz+1]*p[xyz+1]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+1]*weight[pX][xyz]*p[xyz+1];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+1]*weight[pX][xyz];
					}
					if (x>0 && p[xyz-1]*p[xyz-1]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-1]*weight[mX][xyz]*p[xyz-1];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-1]*weight[mX][xyz];
					}
					if (y<ny-1 && p[xyz+nx]*p[xyz+nx]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx]*weight[pY][xyz]*p[xyz+nx];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx]*weight[pY][xyz];
					}
					if (y>0 && p[xyz-nx]*p[xyz-nx]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx]*weight[mY][xyz]*p[xyz-nx];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx]*weight[mY][xyz];
					}
					if (z<nz-1 && p[xyz+nx*ny]*p[xyz+nx*ny]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx*ny]*weight[pY][xyz]*p[xyz+nx*ny];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx*ny]*weight[pY][xyz];
					}
					if (z>0 && p[xyz-nx*ny]*p[xyz-nx*ny]>ratio*p[xyz]*p[xyz]) {
						val += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx*ny]*weight[mY][xyz]*p[xyz-nx*ny];
						sum += 8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx*ny]*weight[mY][xyz];
					}
					p[xyz] = val/sum;
				}
			}
		}
		return p;
	}
	
	public final float[] computeFast1DLabelIntensityGain(int lb, int iter, float ratio) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		float val, sum;
		for (int t=0;t<iter;t++) {
			for (int x=1;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-1]*weight[mX][xyz]*p[xyz-1])/(1.0f + ratio*weight[mX][xyz]);
			}
			for (int x=nx-2;x>=0;x--) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+1]*weight[pX][xyz]*p[xyz+1])/(1.0f + ratio*weight[pX][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=1;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx]*weight[mY][xyz]*p[xyz-nx])/(1.0f + ratio*weight[mY][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=ny-2;y>=0;y--) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx]*weight[pY][xyz]*p[xyz+nx])/(1.0f + ratio*weight[pY][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=1;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz-nx*ny]*weight[mZ][xyz]*p[xyz-nx*ny])/(1.0f + ratio*weight[mZ][xyz]);
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=nz-2;z>=0;z--) {
				int xyz = x+nx*y+nx*ny*z;
				if (labels[xyz]==0) p[xyz] = (p[xyz] + ratio*8.0f*intens[lb-1][xyz]*intens[lb-1][xyz+nx*ny]*weight[pZ][xyz]*p[xyz+nx*ny])/(1.0f + ratio*weight[pZ][xyz]);
			}
		}
		return p;
	}
	
	public final float[] computeFastMarchingGain(byte lb, float geom) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		BinaryHeap2D heap = new BinaryHeap2D(Numerics.ceil(1.25f*nx*ny), Numerics.ceil(0.1f*nx*ny), BinaryHeap2D.MAXTREE);
    	 
		heap.reset();
		float speed, dist;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (p[xyz]==1.0f) {
				dist = 1.0f;
				if (x>0 && p[xyz-1]==0.0f) {
					// add to the heap
					speed = dist*weight[mX][xyz]/(geom+weight[mX][xyz]);
					heap.addValue(speed, xyz-1, lb);
				}
				if (x<nx-1 && p[xyz+1]==0.0f) {
					// add to the heap
					speed = dist*weight[pX][xyz]/(geom+weight[pX][xyz]);
					heap.addValue(speed, xyz+1, lb);
				}
				if (y>0 && p[xyz-nx]==0.0f) {
					// add to the heap
					speed = dist*weight[mY][xyz]/(geom+weight[mY][xyz]);
					heap.addValue(speed, xyz-nx, lb);
				}
				if (y<ny-1 && p[xyz+nx]==0.0f) {
					// add to the heap
					speed = dist*weight[pY][xyz]/(geom+weight[pY][xyz]);
					heap.addValue(speed, xyz+nx, lb);
				}
				if (z>0 && p[xyz-nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[mZ][xyz]/(geom+weight[mZ][xyz]);
					heap.addValue(speed, xyz-nx*ny, lb);
				}
				if (z<nz-1 && p[xyz+nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[pZ][xyz]/(geom+weight[pZ][xyz]);
					heap.addValue(speed, xyz+nx*ny, lb);
				}
			}
		}
			
		// grow the labels and functions	
		int t = 0;		
		while (heap.isNotEmpty() && t<10000000) {
			//t++;
			
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			heap.removeFirst();

			// if the location has been processed already, we stop
			if (p[xyz]!=0.0f)  continue;
			
			// update the value
			p[xyz] = dist;
			//System.out.print(".");
			
			int z = xyz/(nx*ny);
			int y = (xyz-nx*ny*z)/nx;
			int x = xyz-nx*ny*z-nx*y;

			// find new neighbors
			if (x>0 && p[xyz-1]==0.0f) {
				// add to the heap
				speed = dist*weight[mX][xyz]/(geom+weight[mX][xyz]);
				heap.addValue(speed, xyz-1, lb);
			}
			if (x<nx-1 && p[xyz+1]==0.0f) {
				// add to the heap
				speed = dist*weight[pX][xyz]/(geom+weight[pX][xyz]);
				heap.addValue(speed, xyz+1, lb);
			}
			if (y>0 && p[xyz-nx]==0.0f) {
				// add to the heap
				speed = dist*weight[mY][xyz]/(geom+weight[mY][xyz]);
				heap.addValue(speed, xyz-nx, lb);
			}
			if (y<ny-1 && p[xyz+nx]==0.0f) {
				// add to the heap
				speed = dist*weight[pY][xyz]/(geom+weight[pY][xyz]);
				heap.addValue(speed, xyz+nx, lb);
			}
			if (z>0 && p[xyz-nx*ny]==0.0f) {
				// add to the heap
				speed = dist*weight[mZ][xyz]/(geom+weight[mZ][xyz]);
				heap.addValue(speed, xyz-nx*ny, lb);
			}
			if (z<nz-1 && p[xyz+nx*ny]==0.0f) {
				// add to the heap
				speed = dist*weight[pZ][xyz]/(geom+weight[pZ][xyz]);
				heap.addValue(speed, xyz+nx*ny, lb);
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]!=lb && labels[xyz]!=0) p[xyz] = 0.0f;
		}
		
		return p;
	}
	
		public final float[] computeDualFastMarchingGain(byte lb, float geom) {
		float[] p = new float[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]==lb) p[xyz] = 1.0f; 
			else if (labels[xyz]!=0) p[xyz] = -1.0f;
			else p[xyz] = 0.0f;
		}
		
		BinaryHeap2D heap = new BinaryHeap2D(Numerics.ceil(1.25f*nx*ny), Numerics.ceil(0.1f*nx*ny), BinaryHeap2D.MAXTREE);
    	 
		heap.reset();
		float speed, dist;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (p[xyz]==1.0f) {
				dist = 1.0f;
				if (x>0 && p[xyz-1]==0.0f) {
					// add to the heap
					speed = dist*weight[mX][xyz]/(geom+weight[mX][xyz]);
					heap.addValue(speed, xyz-1, lb);
				}
				if (x<nx-1 && p[xyz+1]==0.0f) {
					// add to the heap
					speed = dist*weight[pX][xyz]/(geom+weight[pX][xyz]);
					heap.addValue(speed, xyz+1, lb);
				}
				if (y>0 && p[xyz-nx]==0.0f) {
					// add to the heap
					speed = dist*weight[mY][xyz]/(geom+weight[mY][xyz]);
					heap.addValue(speed, xyz-nx, lb);
				}
				if (y<ny-1 && p[xyz+nx]==0.0f) {
					// add to the heap
					speed = dist*weight[pY][xyz]/(geom+weight[pY][xyz]);
					heap.addValue(speed, xyz+nx, lb);
				}
				if (z>0 && p[xyz-nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[mZ][xyz]/(geom+weight[mZ][xyz]);
					heap.addValue(speed, xyz-nx*ny, lb);
				}
				if (z<nz-1 && p[xyz+nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[pZ][xyz]/(geom+weight[pZ][xyz]);
					heap.addValue(speed, xyz+nx*ny, lb);
				}
			} else if (p[xyz]==-1.0f) {
				dist = 1.0f;
				if (x>0 && p[xyz-1]==0.0f) {
					// add to the heap
					speed = dist*weight[mX][xyz]/(geom+weight[mX][xyz]);
					heap.addValue(speed, xyz-1, (byte)(-lb));
				}
				if (x<nx-1 && p[xyz+1]==0.0f) {
					// add to the heap
					speed = dist*weight[pX][xyz]/(geom+weight[pX][xyz]);
					heap.addValue(speed, xyz+1, (byte)(-lb));
				}
				if (y>0 && p[xyz-nx]==0.0f) {
					// add to the heap
					speed = dist*weight[mY][xyz]/(geom+weight[mY][xyz]);
					heap.addValue(speed, xyz-nx, (byte)(-lb));
				}
				if (y<ny-1 && p[xyz+nx]==0.0f) {
					// add to the heap
					speed = dist*weight[pY][xyz]/(geom+weight[pY][xyz]);
					heap.addValue(speed, xyz+nx, (byte)(-lb));
				}
				if (z>0 && p[xyz-nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[mZ][xyz]/(geom+weight[mZ][xyz]);
					heap.addValue(speed, xyz-nx*ny, (byte)(-lb));
				}
				if (z<nz-1 && p[xyz+nx*ny]==0.0f) {
					// add to the heap
					speed = dist*weight[pZ][xyz]/(geom+weight[pZ][xyz]);
					heap.addValue(speed, xyz+nx*ny, (byte)(-lb));
				}
			}

		}
			
		// grow the labels and functions	
		int t = 0;		
		while (heap.isNotEmpty() && t<10000000) {
			//t++;
			
			// extract point with minimum distance
			dist = heap.getFirst();
			int xyz = heap.getFirstId();
			lb = heap.getFirstState();
			heap.removeFirst();

			// if the location has been processed already, we stop
			if (p[xyz]!=0.0f)  continue;
			
			// update the value
			if (lb>0) p[xyz] = dist;
			else p[xyz] = -dist;
			//System.out.print(".");
			
			int z = xyz/(nx*ny);
			int y = (xyz-nx*ny*z)/nx;
			int x = xyz-nx*ny*z-nx*y;

			// find new neighbors
			if (x>0 && p[xyz-1]==0.0f) {
				// add to the heap
				speed = dist*weight[mX][xyz]/(geom+weight[mX][xyz]);
				heap.addValue(speed, xyz-1, lb);
			}
			if (x<nx-1 && p[xyz+1]==0.0f) {
				// add to the heap
				speed = dist*weight[pX][xyz]/(geom+weight[pX][xyz]);
				heap.addValue(speed, xyz+1, lb);
			}
			if (y>0 && p[xyz-nx]==0.0f) {
				// add to the heap
				speed = dist*weight[mY][xyz]/(geom+weight[mY][xyz]);
				heap.addValue(speed, xyz-nx, lb);
			}
			if (y<ny-1 && p[xyz+nx]==0.0f) {
				// add to the heap
				speed = dist*weight[pY][xyz]/(geom+weight[pY][xyz]);
				heap.addValue(speed, xyz+nx, lb);
			}
			if (z>0 && p[xyz-nx*ny]==0.0f) {
				// add to the heap
				speed = dist*weight[mZ][xyz]/(geom+weight[mZ][xyz]);
				heap.addValue(speed, xyz-nx*ny, lb);
			}
			if (z<nz-1 && p[xyz+nx*ny]==0.0f) {
				// add to the heap
				speed = dist*weight[pZ][xyz]/(geom+weight[pZ][xyz]);
				heap.addValue(speed, xyz+nx*ny, lb);
			}
		}
		/*
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labels[xyz]!=lb && labels[xyz]!=0) p[xyz] = 0.0f;
		}
		*/
		return p;
	}
}
