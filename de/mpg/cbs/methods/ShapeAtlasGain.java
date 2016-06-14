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
 
public class ShapeAtlasGain {
		
	// data buffers
	private 	float[][]		images;  			// original images
	private		int				nc;					// number of images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String[]		modality;			// image modality / contrast
	
	// atlas parameters
	private		SimpleShapeAtlas	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	int				nax,nay,naz;   			// shape dimensions
	private 	float			rax,ray,raz;   			// shape resolutions
	private		float[][]		intensity;		// atlas intensities
	
	// output results
	private		float[][]		min, max;
	private		float[][]		gain;
	private		float			sigmaI, sigmaS;
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
	 public ShapeAtlasGain(float[][] img_, String[] mod_, int nc_, int nix_, int niy_, int niz_, 
							float rix_, float riy_, float riz_, SimpleShapeAtlas atlas_, float p1_, float p2_) {
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
		nax = atlas.getShapeDim()[0];
		nay = atlas.getShapeDim()[1];
		naz = atlas.getShapeDim()[2];
		rax = atlas.getShapeRes()[0];
		ray = atlas.getShapeRes()[1];
		raz = atlas.getShapeRes()[2];		
		intensity = atlas.getIntensityPriors(modality, nc);
		for (int c=0;c<nc;c++) {
			BasicInfo.displayMessage("channel "+c+" intensity priors: "+atlas.displayVector(intensity[c])+"\n");
		}
		
		sigmaI = p1_;
		sigmaS = p2_;
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
    	min = new float[nc][nax*nay*naz];
    	max = new float[nc][nax*nay*naz];
    	
     	for (int c=0;c<nc;c++) for (int xyz=0;xyz<nax*nay*naz;xyz++) {
    		min[c][xyz] = INF;
    		max[c][xyz] = -INF;
    	}
    	
		float[] Xs = new float[3];
		int xyzs, xyzi;
		// scans the entire image only once (note: some atlas locations may not have values)
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
        	xyzi = x + nix*y + nix*niy*z;
        	
        	atlas.imageToShapeCoordinates(Xs, x, y, z);
        	
        	// maps to closest neighbor
        	xyzs = Numerics.round(Numerics.bounded(Xs[0],0,nax-1)) 
        		  + nax*Numerics.round(Numerics.bounded(Xs[1],0,nay-1)) 
        		  + nax*nay*Numerics.round(Numerics.bounded(Xs[2],0,naz-1));
        	
        	for (int c=0;c<nc;c++) {
        		min[c][xyzs] = Numerics.min(min[c][xyzs], images[c][xyzi]);
        		max[c][xyzs] = Numerics.max(max[c][xyzs], images[c][xyzi]);
        	}
		}
		// compute exact image range for normalization in [0,1]
		float Imin, Imax;
		for (int c=0;c<nc;c++) {
			Imin = INF;
			Imax =-INF;
			for (int xyz=0;xyz<nax*nay*naz;xyz++) {
				Imin = Numerics.min(Imin, min[c][xyz]);
				Imax = Numerics.max(Imax, max[c][xyz]);
			}
			// then normalize everything (not robust, though)
			for (int xyz=0;xyz<nax*nay*naz;xyz++) {
				min[c][xyz] = (min[c][xyz] - Imin)/(Imax - Imin);
				max[c][xyz] = (max[c][xyz] - Imin)/(Imax - Imin);
			}
		}
        return;
    }// image min max

	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeAtlasGainFunction() {
    	gain = new float[nobj][nax*nay*naz];
    	
    	// compute the local priors
    	float shape;
    	float[] intens = new float[nc];
    	float val, diff, mindiff, maxdiff;
    	float shapesum, shapemax;
    	int xyzs;
        for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
        	xyzs = x + nax*y + nax*nay*z;
        	
        	for (int n=0;n<nobj;n++) {
				for (int c=0;c<nc;c++) {
					mindiff = (min[c][xyzs]-intensity[c][n])/sigmaI;
					maxdiff = (max[c][xyzs]-intensity[c][n])/sigmaI;
					
					// more sensitive to scale issues
					diff = Numerics.min(mindiff*mindiff,  maxdiff*maxdiff);
					// oblivious of scale
					//val = 0.5f*(mindiff + maxdiff);

					//diff = val*val;
					//if (c==0) diff = val*val;
					//else diff = Numerics.min(diff, val*val);
					intens[c] = (1.0f - diff)/(1.0f + diff);
				}
				
				val = atlas.getShapes()[n][xyzs];
				shape = (val - sigmaS)/(val + sigmaS - 2.0f*sigmaS*val);
				
				float product = intens[0];
				for (int c=1;c<nc;c++) {
					if (intens[c]>0 && product>0) product = intens[c]*product;
					else if (intens[c]<0 && product<0) product = -intens[c]*product;
					else product = 0;
				}
				if (shape>0 && product>0) {
					gain[n][xyzs] = shape*product;
				} else if (shape<0 && product<0) {
					gain[n][xyzs] = -shape*product;
				} else {
					gain[n][xyzs] = 0;	
				}
			}
		}    	
		// quantify the results
		float[] positive = new float[nobj];
		float[] negative = new float[nobj];
		float[]	zero = new float[nobj];
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
        	for (int n=0;n<nobj;n++) {
        		if (gain[n][xyz]>0) positive[n] += gain[n][xyz];
        		else if (gain[n][xyz]<0) negative[n] += gain[n][xyz];
        		else zero[n]++;
        	}
        }
        System.out.println("positive gain: "+atlas.displayVector(positive));
        System.out.println("negative gain: "+atlas.displayVector(negative));
        System.out.println(" zero    gain: "+atlas.displayVector(zero));
    }
    
	private final int neighborIndex(byte d, int id) {
		switch (d) {
			case pX		: 	return id+1; 		
			case mX		:	return id-1;
			case pY		:	return id+nax;
			case mY		:	return id-nax;
			case pZ		:	return id+nax*nay;
			case mZ		:	return id-nax*nay;
			case pXpY	:	return id+1+nax;
			case mXpY	:	return id-1+nax;
			case pYpZ	:	return id+nax+nax*nax;
			case mYpZ	:	return id-nax+nax*nay;
			case pZpX	:	return id+nax*nay+1;	
			case mZpX	:	return id-nax*nay+1;
			case pXmY	:	return id+1-nax;	
			case mXmY	:	return id-1-nax;
			case pYmZ	:	return id+nax-nax*nay;
			case mYmZ	:	return id-nax-nax*nay;
			case pZmX	:	return id+nax*nay-1;
			case mZmX	:	return id-nax*nay-1;
			case pXpYpZ	:	return id+1+nax+nax*nay;
			case mXpYpZ	:	return id-1+nax+nax*nay;
			case pXmYmZ	:	return id+1-nax-nax*nay; 
			case mXmYmZ	:	return id-1-nax-nax*nay;
			case pXmYpZ	:	return id+1-nax+nax*nay;
			case mXmYpZ	:	return id-1-nax+nax*nay;
			case pXpYmZ	:	return id+1+nax-nax*nay; 
			case mXpYmZ	:	return id-1+nax-nax*nay;
			default		:	return id;
		}
	}
    
    public final float[] getBinaryGainValue(int id) {
    	float[] seg = new float[nax*nay*naz];
    	
    	int best;
    	float max;
    	int xyz;
        for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
        	xyz = x + nax*y + nax*nay*z;
        	
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
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(min[c], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
		}
		return result;
    }
    
     public final float[] exportRng(int c) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(max[c], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
		}
		return result;
    }
        
    public final float[] exportGain(int n) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(gain[n], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
		}
		return result;
    }
    
    public final int[] generateGainSegmentation() {
    	int[] result = new int[nax*nay*naz];
    	
    	float[] Xs = new float[3];
    	int id;
    	float val,best;
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
    		id = 0;
    		best = gain[0][xyz];
    		for (int n=1;n<nobj;n++) {
    			val = gain[n][xyz];
    			if (val>best) { best = val; id = n; }
    		}
    		result[xyz] = atlas.getLabels()[id];
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
    		best = ImageInterpolation.linearClosestInterpolation(gain[0], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
    		for (int n=1;n<nobj;n++) {
    			val = ImageInterpolation.linearClosestInterpolation(gain[n], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
    			if (val>best) { best = val; id = n; }
    		}
    		result[x+nix*y+nix*niy*z] = atlas.getLabels()[id];
		}
		return result;
    }
    
}
