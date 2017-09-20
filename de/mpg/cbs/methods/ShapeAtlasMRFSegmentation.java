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
 *  This algorithm performs a simple alignment of the image with a topology template
 *	<p>
 *	The algorithm handles all the little things needed for image cropping, conversion
 * 	to an image buffer, padding, etc.
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class ShapeAtlasMRFSegmentation {
		
	// data buffers
	private 	float[][]		images;  			// original images
	private		int				nc;					// number of images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String[]		modality;			// image modality / contrast
	
	// atlas parameters
	private		SimpleShapeAtlas	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	int				nsx,nsy,nsz;   			// shape dimensions
	private 	float			rsx,rsy,rsz;   			// shape resolutions
	private		float[][]		intensity;		// atlas intensities
	
	// output results
	private		float[][]		avg, rng;
	private		float[][]		similarity;
	private		byte[][]		neighbor;
	private		int				nsim;					// number of similar neighbors
	private		float[][]		gain;
	private		float			sigmaI, sigmaS, sigmaR, sigmaN;
	private		float			factor;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
		
	private	static	final	byte	EMPTY = -1;
			
	public	static	final	byte	pX = 0;
	public	static	final	byte	pY = 1;
	public	static	final	byte	pZ = 2;
	public	static 	final 	byte 	pXpY = 3;
	public	static 	final 	byte 	pYpZ = 4;
	public	static 	final 	byte 	pZpX = 5;
	public	static 	final 	byte 	pXmY = 6;
	public	static 	final 	byte 	pYmZ = 7;
	public	static 	final 	byte 	pZmX = 8;
	public	static 	final 	byte 	pXpYpZ = 9;
	public	static 	final 	byte 	pXmYmZ = 10;
	public	static 	final 	byte 	pXmYpZ = 11;
	public	static 	final 	byte 	pXpYmZ = 12;
	public	static	final	byte	mX = 13;
	public	static	final	byte	mY = 14;
	public	static	final	byte	mZ = 15;
	public	static 	final 	byte 	mXpY = 16;
	public	static 	final 	byte 	mYpZ = 17;
	public	static 	final 	byte 	mZpX = 18;
	public	static 	final 	byte 	mXmY = 19;
	public	static 	final 	byte 	mYmZ = 20;
	public	static 	final 	byte 	mZmX = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	public	static 	final 	byte 	NGB = 26;

	
	/**
	 *  constructor
	 */
	public ShapeAtlasMRFSegmentation(float[][] img_, String[] mod_, int nc_,
										int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_,
										SimpleShapeAtlas atlas_,
										int nsim_,
										float p1_, float p2_, float p3_, float p4_) {
		images = img_;
		modality = mod_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		atlas = atlas_;
		nobj = atlas.getNumber();
		nsx = atlas.getShapeDim()[0];
		nsy = atlas.getShapeDim()[1];
		nsz = atlas.getShapeDim()[2];
		rsx = atlas.getShapeRes()[0];
		rsy = atlas.getShapeRes()[1];
		rsz = atlas.getShapeRes()[2];		
		intensity = atlas.getIntensityPriors(modality, nc);
		
		nsim = nsim_;
		
		//sigma1 = ZERO;
		//sigma2 = 1.0f;
		sigmaI = p1_;
		sigmaS = p2_;
		sigmaN = p3_;
		sigmaR = p4_;
	}

	final public void finalize() {
		images = null;
		System.gc();
	}

	public final float[][] getGain() { return gain; }
	
	/** 
	 *  compute the average and std of image intensity at each atlas voxel
	 */
    final public void computeImageMinMaxValues() {
    	avg = new float[nc][nsx*nsy*nsz];
    	rng = new float[nc][nsx*nsy*nsz];
    	
     	for (int c=0;c<nc;c++) for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
    		avg[c][xyz] = INF;
    		rng[c][xyz] = -INF;
    	}
    	
		float[] Xs = new float[3];
		int xyzs, xyzi;
		// scans the entire image only once (note: some atlas locations may not have values)
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
        	xyzi = x + nix*y + nix*niy*z;
        	
        	atlas.imageToShapeCoordinates(Xs, x, y, z);
        	
        	// maps to closest neighbor
        	xyzs = Numerics.round(Numerics.bounded(Xs[0],0,nsx-1)) 
        		  + nsx*Numerics.round(Numerics.bounded(Xs[1],0,nsy-1)) 
        		  + nsx*nsy*Numerics.round(Numerics.bounded(Xs[2],0,nsz-1));
        	
        	for (int c=0;c<nc;c++) {
        		avg[c][xyzs] = Numerics.min(avg[c][xyzs], images[c][xyzi]);
        		rng[c][xyzs] = Numerics.max(rng[c][xyzs], images[c][xyzi]);
        	}
		}
		// compute exact image range for normalization in [0,1]
		float Imin, Imax;
		for (int c=0;c<nc;c++) {
			Imin = INF;
			Imax =-INF;
			for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
				Imin = Numerics.min(Imin, avg[c][xyz]);
				Imax = Numerics.max(Imax, rng[c][xyz]);
			}
			// then normalize everything (not robust, though)
			for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
				avg[c][xyz] = (avg[c][xyz] - Imin)/(Imax - Imin);
				rng[c][xyz] = (rng[c][xyz] - Imin)/(Imax - Imin);
			}
		}
		// replace the values with avg, range
		float average, range;
		for (int c=0;c<nc;c++) for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			if (avg[c][xyz]>1.0f || rng[c][xyz]<0.0f) {
				avg[c][xyz] = 0.0f;
				rng[c][xyz] = 0.0f;
			} else {
				average = 0.5f*(rng[c][xyz]+avg[c][xyz]);
				range  =  0.5f*(rng[c][xyz]-avg[c][xyz]);
				avg[c][xyz] = average;
				rng[c][xyz] = range;
			}
		}
        return;
    }// image avg rng

	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeInitialGainFunction() {
    	gain = new float[nobj][nsx*nsy*nsz];
    	
    	// compute the local priors
    	float intens, shape, complexity;
    	float val, diff;
    	float shapesum, shapemax;
    	int xyzs;
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	
        	for (int n=0;n<nobj;n++) {
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					//diff += (avg[c][xyzs]-intensity[c][n])/(sigmaI*Numerics.max(ZERO,rng[c][xyzs]));
					diff += (avg[c][xyzs]-intensity[c][n])/sigmaI;
				}
				diff /= nc;
				//intens = 1.0f/( (1.0f+diff*diff)*(1.0f+diff*diff) );
				intens = (1.0f-diff*diff)/(1.0f+diff*diff);
				
				val = atlas.getShapes()[n][xyzs];
				//shape = val;
				shape = (val - sigmaS)/(val + sigmaS - 2.0f*sigmaS*val);
				
				if (shape>0 && intens>0) {
					gain[n][xyzs] = shape*intens;
				} else if (shape<0 && intens<0) {
					gain[n][xyzs] = -shape*intens;
				} else {
					gain[n][xyzs] = 0;	
				}
			}
		}    	
    }
    
	/**
	 *	compute the intensity priors given the atlas
	 */
	/* 
    public final void propagateGainFunction() {
    	
    	// compute the local priors
    	float intens, shape, complexity;
    	float val, diff;
    	float shapesum, shapemax;
    	int xyzs,xyzn;
    	
    	float[][] prev = new float[nobj][nsx*nsy*nsz];
    	for (int n=0;n<nobj;n++)
    		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) 
    			prev[n][xyz] = gain[n][xyz];
    	
        for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	
        	// local gain	
        	shapesum = 0.0f;
        	shapemax = 0.0f;
        	for (int n=0;n<nobj;n++) {
        		shapesum += atlas.getShapes()[n][xyzs];
        		shapemax = Numerics.max(shapemax, atlas.getShapes()[n][xyzs]);
        	}
  			
        	val = 0.0f;
        	for (int c=0;c<nc;c++) {
        		val += rng[c][xyzs]/sigma4;
        	}
        	val /= nc;
        	
        	complexity = 1.0f/(1.0f+val*val);
        	
        	for (int n=0;n<nobj;n++) if (atlas.getShapes()[n][xyzs]>0) {
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					diff += (avg[c][xyzs]-intensity[c][n])/sigmaI;
				}
				diff /= nc;
				intens = 1.0f/( (1.0f+diff*diff)*(1.0f+diff*diff) );
				
				shape = atlas.getShapes()[n][xyzs];
				
				gain[n][xyzs] = shape*intens*complexity;
			}
			
			// propagation from neighbors
			for (int s=0;s<nsim;s++) {
				xyzn = neighborIndex(neighbor[s][xyzs],xyzs);
				for (int n=0;n<nobj;n++) {
					gain[n][xyzs] += factor/nsim*similarity[s][xyzs]*prev[n][xyzn];
				}
			}
        }    	
    }
    */
    /*
    public final void propagateAvgGainFunction() {
    	
    	float intens, shape, complexity;
    	float val, diff;
    	float shapesum, shapemax;
    	int xyzs,xyzn;
    	float sim;
    	
    	float[][] prev = new float[nobj][nsx*nsy*nsz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) { 
    			prev[n][xyz] = gain[n][xyz];
    		}
    	}
    	float change = 0.0f;
		float maxchange = 0.0f;
		
		float[] sum = new float[nobj];
        for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	
        	for (int n=0;n<nobj;n++) sum[n] = 1.0f;
        	
			for (byte d=0;d<NGB;d++) {
				xyzn = neighborIndex(d,xyzs);
				
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					diff += (avg[c][xyzs]-avg[c][xyzn])/sigmaN;
				}
				diff /= nc;
				sim = 1.0f/(1.0f+diff*diff);
				
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					diff += (rng[c][xyzs]-rng[c][xyzn])/sigmaR;
				}
				diff /= nc;
				complexity = 1.0f/(1.0f+diff*diff);
				
				for (int n=0;n<nobj;n++) {
					val = atlas.getShapes()[n][xyzs]/sigmaS;
					shape = val*val/(1.0f + val*val);
					
					gain[n][xyzs] += sim*complexity*shape*prev[n][xyzn];
					sum[n] += sim*complexity*shape;
				}
			}
			
			for (int n=0;n<nobj;n++) {
				gain[n][xyzs] /= sum[n];
				
				change += Numerics.abs(gain[n][xyzs]-prev[n][xyzs]);
				maxchange = Numerics.max(maxchange, Numerics.abs(gain[n][xyzs]-prev[n][xyzs]) );
			}
        }
        change /= nsx*nsy*nsz*nobj;
        
        BasicInfo.displayMessage("avg change: "+change+", max change: "+maxchange+"\n");
    }
    */
    
    public final void propagateMaxGainFunction() {
    	
    	float intens, shape, complexity;
    	float val, diff;
    	float shapesum, shapemax;
    	int xyzs,xyzn;
    	float sim;
    	
    	float[][] prev = new float[nobj][nsx*nsy*nsz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) { 
    			prev[n][xyz] = gain[n][xyz];
    		}
    	}
    	float change = 0.0f;
		float maxchange = 0.0f;
        for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	
			for (byte d=0;d<NGB;d++) {
				xyzn = neighborIndex(d,xyzs);
				
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					diff += (avg[c][xyzs]-avg[c][xyzn])/sigmaN;
				}
				diff /= nc;
				sim = 1.0f/(1.0f+diff*diff);
				
				diff = 0.0f;
				for (int c=0;c<nc;c++) {
					diff += (rng[c][xyzs]-rng[c][xyzn])/sigmaR;
				}
				diff /= nc;
				complexity = 1.0f/(1.0f+diff*diff);
				
				for (int n=0;n<nobj;n++) {
					val = (atlas.getShapes()[n][xyzs]*atlas.getShapes()[n][xyzn])/(sigmaS*sigmaS);
					shape = val/(1.0f + val);
					
					gain[n][xyzs] = Numerics.maxmag(gain[n][xyzs], sim*complexity*shape*prev[n][xyzn]);
				}
			}
			
			for (int n=0;n<nobj;n++) {
				change += Numerics.abs(gain[n][xyzs]-prev[n][xyzs]);
				maxchange = Numerics.max(maxchange, Numerics.abs(gain[n][xyzs]-prev[n][xyzs]) );
			}
        }
        change /= nsx*nsy*nsz*nobj;
        
        BasicInfo.displayMessage("avg change: "+change+", max change: "+maxchange+"\n");
    }
    
	private final int neighborIndex(byte d, int id) {
		switch (d) {
			case pX		: 	return id+1; 		
			case mX		:	return id-1;
			case pY		:	return id+nsx;
			case mY		:	return id-nsx;
			case pZ		:	return id+nsx*nsy;
			case mZ		:	return id-nsx*nsy;
			case pXpY	:	return id+1+nsx;
			case mXpY	:	return id-1+nsx;
			case pYpZ	:	return id+nsx+nsx*nsx;
			case mYpZ	:	return id-nsx+nsx*nsy;
			case pZpX	:	return id+nsx*nsy+1;	
			case mZpX	:	return id-nsx*nsy+1;
			case pXmY	:	return id+1-nsx;	
			case mXmY	:	return id-1-nsx;
			case pYmZ	:	return id+nsx-nsx*nsy;
			case mYmZ	:	return id-nsx-nsx*nsy;
			case pZmX	:	return id+nsx*nsy-1;
			case mZmX	:	return id-nsx*nsy-1;
			case pXpYpZ	:	return id+1+nsx+nsx*nsy;
			case mXpYpZ	:	return id-1+nsx+nsx*nsy;
			case pXmYmZ	:	return id+1-nsx-nsx*nsy; 
			case mXmYmZ	:	return id-1-nsx-nsx*nsy;
			case pXmYpZ	:	return id+1-nsx+nsx*nsy;
			case mXmYpZ	:	return id-1-nsx+nsx*nsy;
			case pXpYmZ	:	return id+1+nsx-nsx*nsy; 
			case mXpYmZ	:	return id-1+nsx-nsx*nsy;
			default		:	return id;
		}
	}
    
    public final float[] getBinaryGainValue(int id) {
    	float[] seg = new float[nsx*nsy*nsz];
    	
    	int best;
    	float max;
    	int xyz;
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	xyz = x + nsx*y + nsx*nsy*z;
        	
        	best = 0;
        	max = gain[0][xyz];
        	for (int n=1;n<nobj;n++) {
				if (gain[n][xyz]>0 && gain[n][xyz]>max) best = n;
			}
			if (best==id) seg[xyz] = 1.0f;
			else seg[xyz] = 0.0f;
		}    	
		return seg;
    }
    
     public final float[] exportAvg(int c) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(avg[c], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportRng(int c) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(rng[c], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportSimilarity(int s) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(similarity[s], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportNeighbor(int s) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(neighbor[s], (byte)(-1), Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportGain(int n) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(gain[n], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportGainSegmentation() {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
    	int id;
    	float val,best;
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		id = 0;
    		best = ImageInterpolation.linearClosestInterpolation(gain[0], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
    		for (int n=1;n<nobj;n++) {
    			val = ImageInterpolation.linearClosestInterpolation(gain[n], Xs[0], Xs[1], Xs[2], nsx, nsy, nsz); 
    			if (val>best) { best = val; id = n; }
    		}
    		result[x+nix*y+nix*niy*z] = atlas.getLabels()[id];
		}
		return result;
    }
    
}
