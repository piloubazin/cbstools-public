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
 
public class FastMarchingSkullStripping {
		
	// data buffers
	private 	float[]			image;  			// original images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String			modality;			// image modality / contrast
	
	// atlas parameters
	private		SimpleShapeAtlas	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	int				nax,nay,naz;   			// shape dimensions
	private 	float			rax,ray,raz;   			// shape resolutions
	private		float[]			intensity;		// atlas intensities

	// segmentation parameters
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private 	BinaryHeap2D		tree;		   		// the binary tree used for the fast marching
	private		float[]				score;
	private		byte[]				segmentation;
	private		boolean[]			region;
	private		byte				from, in, to;
	
	// output results
	private		int				nsx,nsy,nsz;
	private 	float			rsx,rsy,rsz;   			// shape resolutions
	private		float			scaling;
	private		float[]			min, max;
	private		float[][]		gain;
	private		float			sigmaI, sigmaS, sigmaN;
	private		float			factor;
	private		float			threshold1;
	private		float			threshold2;
	private		float			slope;
	private		int				ngbsize;
	
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
	public	static	final	byte	mX = 3;
	public	static	final	byte	mY = 4;
	public	static	final	byte	mZ = 5;
	public	static 	final 	byte 	pXpY = 6;
	public	static 	final 	byte 	pYpZ = 7;
	public	static 	final 	byte 	pZpX = 8;
	public	static 	final 	byte 	pXmY = 9;
	public	static 	final 	byte 	pYmZ = 10;
	public	static 	final 	byte 	pZmX = 11;
	public	static 	final 	byte 	mXpY = 12;
	public	static 	final 	byte 	mYpZ = 13;
	public	static 	final 	byte 	mZpX = 14;
	public	static 	final 	byte 	mXmY = 15;
	public	static 	final 	byte 	mYmZ = 16;
	public	static 	final 	byte 	mZmX = 17;
	public	static 	final 	byte 	pXpYpZ = 18;
	public	static 	final 	byte 	pXmYmZ = 19;
	public	static 	final 	byte 	pXmYpZ = 20;
	public	static 	final 	byte 	pXpYmZ = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	public	static 	final 	byte 	NGB = 26;

	private static final byte[] ngbx = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	private static final byte[] ngby = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	private static final byte[] ngbz = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	
	/**
	 *  constructor
	 */
	public FastMarchingSkullStripping(float[] img_, String mod_, 
										int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_,
										SimpleShapeAtlas atlas_,
										float scaling_,
										float p1_, float p2_, float p3_,
										float p4_, float p5_, float p6_,
										int n1_) {
		image = img_;
		modality = mod_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		scaling = scaling_;
		nsx = Numerics.ceil(scaling*nix);
		nsy = Numerics.ceil(scaling*niy);
		nsz = Numerics.ceil(scaling*niz);
		
		rsx = rix/scaling;
		rsy = riy/scaling;
		rsz = riz/scaling;
		
		atlas = atlas_;
		nobj = atlas.getNumber();
		nax = atlas.getShapeDim()[0];
		nay = atlas.getShapeDim()[1];
		naz = atlas.getShapeDim()[2];
		rax = atlas.getShapeRes()[0];
		ray = atlas.getShapeRes()[1];
		raz = atlas.getShapeRes()[2];		
		intensity = atlas.getIntensityPriors(new String[] {modality}, 1)[0];
		BasicInfo.displayMessage("intensity priors: "+atlas.displayVector(intensity)+"\n");
		
		sigmaI = p1_;
		sigmaN = p2_;
		
		sigmaS = p3_;
		
		threshold1 = p4_;
		threshold2 = p5_;
		
		slope	= p6_;
		
		ngbsize = n1_;
	}

	final public void finalize() {
		image = null;
		System.gc();
	}

	public final float[][] getGain() { return gain; }
	
	/** 
	 *  compute the average and std of image intensity at each atlas voxel
	 */
    final public void computeImageMinMaxValues() {
    	min = new float[nsx*nsy*nsz];
    	max = new float[nsx*nsy*nsz];
    	
     	for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
    		min[xyz] = INF;
    		max[xyz] = -INF;
    	}
    	
		float[] Xs = new float[3];
		int xyzs, xyzi;
		// scans the entire image only once (note: some atlas locations may not have values)
       for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
       //for (int x=ngbsize;x<nix-ngbsize;x++) for (int y=ngbsize;y<niy-ngbsize;y++) for (int z=ngbsize;z<niz-ngbsize;z++) {
	        xyzi = x + nix*y + nix*niy*z;
        	xyzs = Numerics.floor(x*scaling) + nsx*Numerics.floor(y*scaling) + nsx*nsy*Numerics.floor(z*scaling);
        	
        	min[xyzs] = Numerics.min(min[xyzs], image[xyzi]);
        	max[xyzs] = Numerics.max(max[xyzs], image[xyzi]);
        	// use a region
        	//for (int i=-ngbsize;i<=ngbsize;i++) for (int j=-ngbsize;j<=ngbsize;j++) for (int l=-ngbsize;l<=ngbsize;l++) {
        	for (int n=0;n<ngbsize;n++) { 
				min[xyzs] = Numerics.min(min[xyzs], image[xyzi+ngbx[n]+nix*ngby[n]+nix*niy*ngbz[n]]);
				max[xyzs] = Numerics.max(max[xyzs], image[xyzi+ngbx[n]+nix*ngby[n]+nix*niy*ngbz[n]]);
			}
        }
		// compute exact image range for normalization in [0,1]
		float Imin, Imax;

		Imin = INF;
		Imax =-INF;
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			Imin = Numerics.min(Imin, min[xyz]);
			Imax = Numerics.max(Imax, max[xyz]);
		}
		// then normalize everything (not robust, though)
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			min[xyz] = (min[xyz] - Imin)/(Imax - Imin);
			max[xyz] = (max[xyz] - Imin)/(Imax - Imin);
		}
		/*
		// replace the values with min, range
		float average, range;
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			if (min[xyz]>1.0f || max[xyz]<0.0f) {
				min[xyz] = 0.0f;
				max[xyz] = 0.0f;
			} else {
				average = 0.5f*(max[xyz]+min[xyz]);
				range  =  0.5f*(max[xyz]-min[xyz]);
				min[xyz] = average;
				max[xyz] = range;
			}
		}
		*/
        return;
    }// image min max
   
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeAtlasGainFunction() {
    	gain = new float[nobj][nsx*nsy*nsz];
    	
    	// compute the local priors
    	float intens, shape, complexity;
    	float val, diff, mindiff, maxdiff;
    	float shapesum, shapemax;
    	int xyzs;
    	float[] XA = new float[3];
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	atlas.imageToShapeCoordinates(XA, (x+0.5f)/scaling, (y+0.5f)/scaling, (z+0.5f)/scaling);
        	
        	for (int n=0;n<nobj;n++) {
				mindiff = (min[xyzs]-intensity[n])/sigmaI;
				maxdiff = (max[xyzs]-intensity[n])/sigmaI;
				diff = Numerics.max(mindiff*mindiff,  maxdiff*maxdiff);
				intens = (1.0f - diff)/(1.0f + diff);
				
				val = ImageInterpolation.linearInterpolation(atlas.getShape(n),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
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
	 *	morphological post-processing
	 */
    public final void closeMask(int d) {
    	
    	region = ObjectExtraction.objectFromLabelImage(segmentation, nsx, nsy, nsz, (byte)0, ObjectExtraction.SUPERIOR);
    	ObjectMorphology.fastDilateObject(region, nsx, nsy, nsz, d);
    	region = ObjectLabeling.removeHoles(region, nsx, nsy, nsz, 6);
    	ObjectMorphology.fastErodeObject(region, nsx, nsy, nsz, d); 
    	
    }
    
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void propagateGainSegmentation(String start, String through, String stop) {
    	
    	from = (byte)atlas.getObject(start);
    	in = (byte)atlas.getObject(through);
    	to = (byte)atlas.getObject(stop);
    	
    	// start from the largest region with positive gain
    	boolean[] obj = new boolean[nsx*nsy*nsz];
    	boolean[] boundary = new boolean[nsx*nsy*nsz];
    	segmentation = new byte[nsx*nsy*nsz];
    	score = new float[nsx*nsy*nsz];
    	boolean[] mask = new boolean[nsx*nsy*nsz];
    	int xyzs,ngb;
     	for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
    		if (x<2 || x>nsx-3 || y<2 || y>nsy-3 || z<2 || z>nsz-3) {
        		mask[xyzs] = false;
        		obj[xyzs] = false;
        	} else {
        		mask[xyzs] = true;
        		if (gain[from][xyzs]>0) {
        			obj[xyzs] = true;
				} else {
					obj[xyzs] = false;
				}
			}
    	}
    	obj = ObjectLabeling.largestObject(obj, nsx, nsy, nsz, 6);
    	
     	for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	xyzs = x + nsx*y + nsx*nsy*z;
        	if (!obj[xyzs]) {
				segmentation[xyzs] = 0;
				score[xyzs] = 0.0f;
			} else {
				segmentation[xyzs] = from;
				score[xyzs] = 1.0f;
			}
			boundary[xyzs] = false;
		}
		
    	// build a binary tree at the boundary
    	tree = new BinaryHeap2D(nix+niy+niz, BinaryHeap2D.MAXTREE);
    	float diff, sim;
    	float mindiff,maxdiff;
    	float valf, vali, valt;
    	float val, prev;
    	byte lb;
    	float sgn = intensity[to]-intensity[from];
    	float proba;
    	float[] XA = new float[3];
    	for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x + nsx*y + nsx*nsy*z;
			if (obj[xyzs]) {
				// on a boundary?
				for (byte n=0;n<26;n++) {
					ngb = neighborIndex(n, xyzs);
					if (!obj[ngb] && mask[ngb]) {
						if (!boundary[ngb]) {
							if (sgn>0) addIncreasingNeighbor(xyzs, ngb, 1.0f, from, boundary);
							else addDecreasingNeighbor(xyzs, ngb, 1.0f, from, boundary);
						}
					}
				}
			}
    	}
    	
    	// propagate out following regions of lowest variation
    	while (tree.isNotEmpty()) {
    		
    		// get the next value
			prev = tree.getFirst();
			xyzs = tree.getFirstId();
			lb = tree.getFirstState();
			tree.removeFirst();
	
    		if (!obj[xyzs]) {
    			
    			// add to object (topology free)
    			obj[xyzs]=true;
    			segmentation[xyzs] = lb;
    			score[xyzs] = prev;
    			boundary[xyzs] = false;
    			
				// find the neighbors
				for (byte n=0;n<26;n++) {
					ngb = neighborIndex(n, xyzs);
					if (!obj[ngb] && mask[ngb]) {
						if (!boundary[xyzs]) {
							if (sgn>0) addIncreasingNeighbor(xyzs, ngb, prev, lb, boundary);
							else addDecreasingNeighbor(xyzs, ngb, prev, lb, boundary);
						}
					}
				}
    		}
    	}
    	return;
    }
    
    private final void addIncreasingNeighbor(int xyz, int ngb, float prev, byte lb0, boolean[] used) {
		float diff, rng, sim, val;
		byte lb;
		
		// new label
		if (max[ngb]>intensity[to]) lb = to;
		else if (min[ngb]>intensity[from]) lb = in;
		else lb = from;
		
		// similarity function : depends on type?
		//mindiff = -1.0f;
		//if (lb==to) mindiff = Numerics.min(max[ngb]-max[xyz], min[ngb]-min[xyz]);
		//if (lb==in) mindiff = 0.5f*(max[ngb]-max[xyz] + min[ngb]-min[xyz]);
		//if (lb==from) mindiff = Numerics.max(max[ngb]-max[xyz], min[ngb]-min[xyz]);
		
		// trying diffrent diff values:
		diff = Numerics.min(max[ngb]-max[xyz], min[ngb]-min[xyz]);
		//diff = 0.5f*(max[ngb]-max[xyz]+min[ngb]-min[xyz]);
		//rng =  0.5f*(max[ngb]-max[xyz]-min[ngb]+min[xyz]);
		//maxdiff = Numerics.max(max[ngb]-max[xyz], min[ngb]-min[xyz]);
		//diff = min[ngb]-min[xyz];
		//diff = max[ngb]-max[xyz];
		
		//diff = ( (1.0f-min[ngb])*(min[ngb]-min[xyz]) + (max[ngb])*(max[ngb]-max[xyz]) ) / (1.0f-min[ngb]+max[ngb]);
		
		sim = 1.0f/( 1.0f + (diff*diff)/(sigmaN*sigmaN) );
				
		// label differently with cases
		val = sim*slope*prev;
			
		// stop at target? not good at growing the csf
		//if (lb0==to) return;
		//if (lb0==to) if (min[ngb]<intensity[in]) return;
		
		// no going back!
		if (lb0==in && lb==from) return;
		if (lb0==to && (lb==in || lb==from) ) return;
		
		// stop if decreasing, definitely in wm and gm
		if (diff<0) {
			//if (rng>0) {
				if (lb==from && max[ngb]>intensity[from]) return;
				if (lb==in) return;
				if (lb==to && min[ngb]<intensity[to]) return; 
			//}
		}
		
		tree.addValue(val, ngb, lb);
		used[ngb] = true;
    }
        
    private final void addDecreasingNeighbor(int xyz, int ngb, float prev, byte lb0, boolean[] used) {
		float diff;
		
		//diff = (max[xyz]-max[ngb])/sigmaN;
		diff = (min[xyz]-min[ngb])/sigmaN;
		
		float sim;
		if (diff>0) {
			sim = 1.0f/(1.0f+diff*diff);
			//sim = diff*diff/(1.0f+diff*diff);
		} else {
			sim = 0.0f;
		}
		
		// estimate whether or not we've arrived to the target region
		float valt,vali,valf;
		
		diff = (max[ngb]-intensity[from])/sigmaI;
		valf = (1.0f - diff*diff)/(1.0f + diff*diff);
		diff = (0.5f*(min[ngb]+max[ngb])-intensity[in])/sigmaI;
		vali = (1.0f - diff*diff)/(1.0f + diff*diff);
		diff = (min[ngb]-intensity[to])/sigmaI;
		valt = (1.0f - diff*diff)/(1.0f + diff*diff);

		// label differently with cases
		byte lb;
		if (valf>vali && valf>valt) {
			lb = from;
		} else if (vali>valf && vali>valt) {
			lb = in;
		} else {
			lb = to;
		}
		// label differently with cases
		float val = sim*prev;
			
		tree.addValue(val, ngb, lb);
    	used[ngb] = true;
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
	
	private final float getAtlasPrior(byte lb, int xyz) {
		int x,y,z;
		z = xyz/(nsx*nsy);
		y = (xyz-nsx*nsy*z)/nsx;
		x = xyz-nsx*nsy*z-nsx*y;
		float[] XA = new float[3];
		atlas.imageToShapeCoordinates(XA, (x+0.5f)/scaling, (y+0.5f)/scaling, (z+0.5f)/scaling);
 		
		return ImageInterpolation.linearInterpolation(atlas.getShape(lb),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
	}
    
    public final float[] exportAvg() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
        	result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(min, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
     public final float[] exportRng() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(max, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportGain(int n) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(gain[n], x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    
    public final float[] exportProba(byte n) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		int xyz = Numerics.floor(x*scaling) + nsx*Numerics.floor(y*scaling) + nsx*nsy*Numerics.floor(z*scaling);
    		result[x+nix*y+nix*niy*z] = getAtlasPrior(n, xyz);
		}
		return result;
    }
    
    public final float[] exportSegmentation() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.nearestNeighborInterpolation(segmentation, (byte)0, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
 
    public final float[] exportRegion() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		if (ImageInterpolation.nearestNeighborInterpolation(region, false, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz))
    			result[x+nix*y+nix*niy*z] = 1.0f;
    		else
    			result[x+nix*y+nix*niy*z] = 0.0f;
		}
		return result;
    }
 
    public final float[] exportScore() {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(score, x*scaling, y*scaling, z*scaling, nsx, nsy, nsz); 
		}
		return result;
    }
    

}
