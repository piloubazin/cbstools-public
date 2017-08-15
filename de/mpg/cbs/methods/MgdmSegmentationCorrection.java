package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import org.apache.commons.math3.util.FastMath;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm uses the MGDM framework to evolve a labeling
 *	according to internal (curvature) and external (vector field) forces
 *	note: for now, we are assuming isotropic voxels 	
 *
 *	@version    August 2012
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class MgdmSegmentationCorrection {
	
	
	// image data
	private 	byte[]			corrected;  		// corrected label image
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		float[][]		bestgain;			// best gain function
	private		byte[][]		bestlabel;			// corresponding labels
	private		int				ngain;
	
	// segmentation parameters
	private 	int 			nobj;    			// number of shapes
	private 	int[]			objLabel;			// label values in the original image
	    
	// data and membership buffers
	private 	float[][] 		mgdmfunctions;  	// MGDM's pseudo level set mgdmfunctions
	private 	byte[][] 		mgdmlabels;   		// MGDM's label maps
	private 	byte[] 			segmentation;   	// MGDM's segmentation
	private		short[]			counter;
	private		boolean[]		mask;				// masking regions not used in computations
	private static	byte 	  	nmgdm;					// total number of MGDM mgdmlabels and mgdmfunctions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// parameters
	private	double		smoothweight, forceweight;
	private double		gaindist2 = 25.0;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 4.0f;
	private	float		narrowBandDist = landmineDist+1.8f;
	private	float		extraDist = narrowBandDist+1.0f;
	private	short		maxcount = 5;
	
	// computation variables to avoid re-allocating
	
	// for levesetForces
	double[] phi;
	double[] Dmx, Dmy, Dmz, Dpx, Dpy, Dpz;
	double D0x,D0y,D0z;
	double SD0x, SD0y, SD0z, GPhi;
	double Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
	double K, G, tmp;
	double Div;
	byte bestlb, gainlb;
	double bestval, gainval;
	float[] smoothfactor;
	boolean done;
	double[] distval;
	
	// for homeomorphicLabeling
	boolean[][][] obj = new boolean[3][3][3];
	int Nconfiguration;
	short[] lbs = new short[26];
	boolean found;
		
	// for minimumMarchingDistance
	double s, s2; 
    int count;
    double dist;
        
	// useful constants & flags
	private	static	final	byte	EMPTY = -1;
	
	// neighborhood flags (ordering is important!!!)
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
	public	static 	final 	byte 	CTR = 26;

	private static final byte[] ngbx = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	private static final byte[] ngby = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	private static final byte[] ngbz = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

    public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;
	public	static	final	byte	N = 3;
	
	// numerical quantities
	private static final	float   INF=1e15f;
	private static final	float   ZERO=1e-15f;
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private final static float SQR2 = (float) Math.sqrt(2.0f);
    private final static float SQR3 = (float) Math.sqrt(3.0f);
    private final static float diagdist = 1/(2*SQR2);
    private final static float cubedist = 1/(2*SQR3);
	private static final	float	UNKNOWN = -1.0f;

	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public MgdmSegmentationCorrection(byte[] corr_,
									int nix_, int niy_, int niz_, float rix_, float riy_, float riz_,
									byte[] init_, float[] bound_, int nmgdm_, 
									byte[][] lb_, float[][] fn_, int ngain_,
									float fw_, float sw_,
									String connectivityType_) {
	
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		nmgdm = (byte)nmgdm_;
		
		bestgain = fn_;
		bestlabel = lb_;
		ngain = ngain_;
						
		forceweight = fw_;
		smoothweight = sw_;
		
		if (verbose) System.out.print("MGDM parameters: "+fw_+" (force), "+sw_+" (smoothing)\n");
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nix, -nix, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nix*niy, -nix*niy};

		// init all the arrays
		try {
			segmentation = new byte[nix*niy*niz];	
			corrected = new byte[nix*niy*niz];	
			counter = new short[nix*niy*niz];	
			mgdmfunctions = new float[nmgdm][nix*niy*niz];
			mgdmlabels = new byte[nmgdm+1][nix*niy*niz];	
			phi = new double[27];
			
			Dmx = new double[nmgdm+1];
			Dmy = new double[nmgdm+1];
			Dmz = new double[nmgdm+1];
			Dpx = new double[nmgdm+1];
			Dpy = new double[nmgdm+1];
			Dpz = new double[nmgdm+1];
			
			distval = new double[nmgdm+1];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nix*niy+niy*niz+niz*nix, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT("critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT("criticalWCLUT.raw.gz",200);
				checkComposed=false;
			}
			else if (connectivityType_.equals("wco")) {
				lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
				checkComposed=true;
			}
			else if (connectivityType_.equals("no")) {
				lut = null;
				checkTopology=false;
			}
			else {
				lut = null;
				checkTopology=false;
			}
			if (checkTopology) {
				if (!lut.loadCompressedPattern()) {
					finalize();
					System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
					BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
				} else {
					//if (verbose) System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		
		// init segmentation
		objLabel = ObjectLabeling.listLabels(init_, nix, niy, niz);
		nobj = objLabel.length;
		int[] inverseLabel = ObjectLabeling.invertList(objLabel);
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			segmentation[xyz] = (byte)inverseLabel[init_[xyz]];
			corrected[xyz] = (byte)inverseLabel[corr_[xyz]];
			for (int b=0;b<ngain;b++) {
				if (bestlabel[b][xyz]!= EMPTY) {
					bestlabel[b][xyz] = (byte)inverseLabel[bestlabel[b][xyz]];
				}
			}
			/*
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==init_[xyz]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
			
			nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==corr_[xyz]) {
					nlb = n;
					continue;
				}
			}
			corrected[xyz] = nlb;
			
			for (int b=0;b<ngain;b++) {
				nlb = EMPTY;
				for (byte n=0; n<nobj; n++) {
					if (objLabel[n]==bestlabel[b][xyz]) {
						nlb = n;
						continue;
					}
				}
				bestlabel[b][xyz] = nlb;
			}
			*/
			
			counter[xyz] = 0;
		}
		computeCorrectionMask(corr_, init_);
		
		// init MGDM boundaries
		mgdmfunctions[0] = bound_;

		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		mgdmfunctions = null;
		mgdmlabels = null;
		segmentation = null;
		lut = null;
		heap = null;
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
		mgdmfunctions = null;
		mgdmlabels = null;
		heap.finalize();
		heap = null;
		System.gc();
	}

	public final float[][] getFunctions() { return mgdmfunctions; }
	
	public final byte[][] getLabels() { return mgdmlabels; }
    
	public final byte[] getSegmentation() { return segmentation; }
	
   public final boolean[] getMask() { return mask; }
	
    
    /** precompute the edge map and its gradient based on a smoothed version of the image
     *	(smoothing scale is in mm) */
    public final void computeCorrectionMask(byte[] lb1, byte[] lb2) {
    	
    	mask = new boolean[nix*niy*niz];
    	// 1. find all the changes
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		if (lb1[xyz]!=lb2[xyz]) mask[xyz] = true;
    		else mask[xyz] = false;
 		}
 		
 		// extend a mask of size equal to the narrow band
 		ObjectMorphology.fastDilateObject(mask, nix, niy, niz, Numerics.ceil(narrowBandDist));

		// remove two layers off the images (for avoiding limits)
		for (int x=0; x<nix; x++) for (int y=0; y<niy; y++) for (int z = 0; z<niz; z++) {
			if (x<2 && x>nix-3 && y<2 && y>niy-3 && z<2 && z>niz-3) mask[x+nix*y+nix*niy*z] = false;
		}
		
    	return;
    }
    
	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedLevelSets() {
    	float[][] levelsets = new float[nobj][nix*niy*niz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    			if (mgdmlabels[0][xyz]==n) levelsets[n][xyz] = -mgdmfunctions[0][xyz];
    			else  levelsets[n][xyz] = 0.0f;
    			
    			for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=n;l++) {
    				levelsets[n][xyz] += mgdmfunctions[l][xyz];
    			}
    		}
    	}
    	return levelsets;
    }
    /**
     *	reconstruct the levelset only where it is guaranteed to be exact 
     */
    public final float[][] reconstructedExactLevelSets() {
    	float[][] levelsets = new float[nobj][nix*niy*niz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    			if (mgdmlabels[0][xyz]==n) levelsets[n][xyz] = -mgdmfunctions[0][xyz];
    			else  levelsets[n][xyz] = 0.0f;
    			
    			int max=0;
    			for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=n;l++) {
    				levelsets[n][xyz] += mgdmfunctions[l][xyz];
    				max++;
    			}
    			if (max==nmgdm) levelsets[n][xyz] = UNKNOWN;
    		}
    	}
    	return levelsets;
    }

	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[] reconstructedLevelsetForce(byte n) {
    	float[] forces = new float[nix*niy*niz];
    	double[] force = new double[nmgdm+1];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
    		levelsetForces(xyz, force);
    		boolean found=false;
    		int lb = -1;
    		for (int l=0;l<nmgdm;l++) if (mgdmlabels[l][xyz]==n && mgdmfunctions[l][xyz]!=UNKNOWN) {
    			found = true;
    			lb = l;
    		}
    		if (found) {
				forces[xyz] = (float)force[lb];
			} else {
				forces[xyz] = 0;
			}
    	}
    	return forces;
    }
 	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final int[][] reconstructedLabels() {
    	int[][] labels = new int[nmgdm][nix*niy*niz];
    	for (int n=0;n<nmgdm;n++) {
    		for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    			if (mgdmlabels[n][xyz]>-1) {
					labels[n][xyz] = objLabel[mgdmlabels[n][xyz]];
				}
			}
    	}
    	return labels;
    }
	/**
	 *	get a final segmentation
	 */
    public final float[] exportSegmentation() {
    	float[] seg = new float[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		//if (mgdmlabels[0][xyz]>-1) {
			if (segmentation[xyz]>-1) {
				//seg[xyz] = objLabel[mgdmlabels[0][xyz]];
				seg[xyz] = objLabel[segmentation[xyz]];
			}
    	}
    	return seg;
    }
    public final float[] exportCounter() {
    	float[] ct = new float[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		ct[xyz] = counter[xyz];
    	}
    	return ct;
    }
    public final float[] exportLabel(int n) {
    	float[] lbs = new float[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		if (mgdmlabels[n][xyz]>-1) {
				lbs[xyz] = objLabel[mgdmlabels[n][xyz]];
			}
    	}
    	return lbs;
    }
    public final float[][][][] exportForces() {
    	float[][][][] result = new float[nix][niy][niz][nmgdm+1];
    	
    	double[] vals = new double[nmgdm+1];
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) if (mask[x+nix*y+nix*niy*z]) {
    		levelsetForces(x+nix*y+nix*niy*z, vals);
    		for (int n=0;n<=nmgdm;n++) result[x][y][z][n] = (float)vals[n];
		}
		return result;
    }
    public final float[] exportMap(float[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = map[x+nix*y+nix*niy*z]; 
		}
		return result;
    }
    
     public final float[] exportMap(byte[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = map[x+nix*y+nix*niy*z]; 
		}
		return result;
    }
    
     public final float[] exportMap(int[] map) {
    	float[] result = new float[nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		result[x+nix*y+nix*niy*z] = map[x+nix*y+nix*niy*z]; 
		}
		return result;
    }
    
    public final float[][] exportBestGainFunctions() {
    	float[][] result = new float[ngain][nix*niy*niz];
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) for (byte n=0;n<ngain;n++) {
     		result[n][x+nix*y+nix*niy*z] = bestgain[n][x+nix*y+nix*niy*z];
		}
		return result;
    }
    
    public final float[] exportBestGainLabels(int n) {
    	float[] seg = new float[nix*niy*niz];
    	int id;
    	for (int x=0; x<nix; x++) for (int y=0; y<niy; y++) for (int z=0; z<niz; z++) {
    		id = bestlabel[n][x+nix*y+nix*niy*z];
    		if (id>-1) seg[x+nix*y+nix*niy*z] = objLabel[id];
    		else seg[x+nix*y+nix*niy*z] = -1;
    	}
		return seg;
    }

	/**
	 *	get the processing flags
	 */
    public final float[] exportProcessingType() {
    	float[] type = new float[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		type[xyz] = 0;
			if (mask[xyz]) type[xyz] = -1;
			if (mgdmfunctions[0][xyz]<narrowBandDist) type[xyz] = 1;
			if (mgdmfunctions[0][xyz]>=landmineDist) type[xyz] = 2;
			if (mgdmfunctions[0][xyz]<1.0f) type[xyz] = 3;
    		if (mgdmfunctions[0][xyz]==UNKNOWN) type[xyz] = -2;
    	}
    	return type;
    }
    
	/**
	 *	get a final segmentation
	 */
    public final int[] labelSegmentation() {
    	int[] seg = new int[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) {
    		//if (mgdmlabels[0][xyz]>-1) {
			if (segmentation[xyz]>-1) {
				//seg[xyz] = objLabel[mgdmlabels[0][xyz]];
				seg[xyz] = objLabel[segmentation[xyz]];
			}
    	}
    	return seg;
    }
   
	/**
	 *	get a final segmentation
	 */
    public final int[][][] exportLabelSegmentation() {
    	int[][][] seg = new int[nix][niy][niz];
    	for (int x=0; x<nix; x++) for (int y=0; y<niy; y++) for (int z=0; z<niz; z++) {
    		int xyz = x + nix*y + nix*niy*z;
    		if (mgdmlabels[0][xyz]>-1) {
				seg[x][y][z] = objLabel[mgdmlabels[0][xyz]];
			}
    	}
    	return seg;
    }

    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	boolean fullMarching = false;
    	
		// init decomposition
		fastMarchingInitialization(false, fullMarching);
			
    	// first estimate the narrow band size
    	int size = 0;
		int boundarysize=0;
		
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) size++;
			if (Numerics.abs(mgdmfunctions[0][xyz])<1.0 && mgdmfunctions[0][xyz]!=UNKNOWN) boundarysize++;
		}
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) {
				narrowband.addPoint(xyz, mgdmlabels, mgdmfunctions);
				// in addition, if close to the narrow band boundariy, set a landmine
				if (mgdmfunctions[0][xyz]>=landmineDist) {
					landmines.set(xyz,true);
				}
			}
		}
		int[] nswap = new int[nmgdm];
			
		double[] forces = new double[nmgdm+1];
		int newlb;
		boolean reinitLM, reinitOL;
		int ncounted;
		
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && (t<5 || diff>mindiff);t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        	if (verbose) System.out.print(t+": ");
        			
			reinitLM = false;
			reinitOL = false;
			
			for (int lb=0;lb<nmgdm;lb++) nswap[lb] = 0;
			
			ncounted = 0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				//if (debug) System.out.print(".");
				// skip points that have seen little change
				if (counter[xyz]<maxcount) {
					// evolve the MGDM functions
					
					// compute the forces from current levelset values, update the narrow band from it
					levelsetForces(xyz, forces);
					//if (debug) System.out.print(":");
				
					for (int lb=nmgdm-1;lb>=0;lb--) if (mgdmlabels[lb][xyz]!=EMPTY) {
						
						// update the narrow band values, not the original data
						narrowband.functions[lb][n] += Numerics.bounded(forces[lb] - forces[lb+1], -0.9f, 0.9f);
						
						// change of sign ?
						if (narrowband.functions[lb][n]<0) {
							//if (debug) System.out.print(""+lb);
							
							// try all possible labels
							newlb = EMPTY;
							if (mgdmlabels[lb+1][xyz]!=EMPTY && (lb>0 || homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz])))
								newlb = lb+1;
							
							if (newlb!=EMPTY) {
								// for all levels
								nswap[lb]++;
								narrowband.labels[lb][n] = mgdmlabels[newlb][xyz];
								narrowband.functions[lb][n] = -narrowband.functions[lb][n];
								// never switch the last label
								if (newlb<nmgdm) narrowband.labels[newlb][n] = mgdmlabels[lb][xyz];
								// update the segmentation with first label
								if (lb==0) {
									segmentation[xyz] = mgdmlabels[newlb][xyz];
									// check for boundary changes in the landmines : force reinitialization
									if (landmines.get(xyz)) reinitLM = true;
									// check for far labels getting mixed in: time to re-initialize
									if (narrowband.labels[0][n]==mgdmlabels[nmgdm][xyz]) reinitOL = true;
								}
							} else {
								// reset to low value?
								narrowband.functions[lb][n] = lowlevel;
								// or to original value ?
								//narrowband.functions[lb][n] -= Numerics.bounded(forces[lb] - forces[lb+1], -0.9f, 0.9f);
							}
							// reset the counter; also reset all the neighbors
							if (lb==0) {
								counter[xyz]=0;
								for (int b=0;b<NGB;b++) {
									int xyzn = xyz + ngbx[b] + ngby[b]*nix + ngbz[b]*nix*niy;
									if (counter[xyzn]>maxcount) counter[xyzn] = (short)(maxcount-1);
									else if (counter[xyzn]>0) counter[xyzn]--;
								}
							}
						} else {
							if (lb==0) counter[xyz]++;
						}
					}
				} else {
					counter[xyz]++;
					ncounted++;
				}
			}
			diff = (nswap[0]/(float)boundarysize);
			if (debug) for (int lb=0;lb<nmgdm;lb++) System.out.print("changed labels ("+lb+"): "+nswap[lb]+" ("+(nswap[lb]/(float)boundarysize*100.0f)+" % of boundary)\n");
			if (debug) System.out.print("frozen points: "+(ncounted/(float)narrowband.currentsize*100.0f)+" % of narrow band\n");
			if (verbose) System.out.print("changed "+(nswap[0]/(float)boundarysize*100.0f)+" % of boundary "
											+"(frozen points: "+(ncounted/(float)narrowband.currentsize*100.0f)+" % of narrow band)\n");
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				if (counter[xyz]<maxcount) {
					for (int lb=0;lb<nmgdm;lb++) {
						mgdmlabels[lb][xyz] = narrowband.labels[lb][n];
						mgdmfunctions[lb][xyz] = narrowband.functions[lb][n];
					}
				}
			}
			// important to check for changes in labels (can get stuck otherwise)
			if (t<iter-1 && (t<5 || diff>mindiff) && (reinitLM || reinitOL) ) {
			//if (t<iter-1 && reinitLM) {
				if (debug) System.out.print("re-initialization (LM: "+reinitLM+" | OL: "+reinitOL+" )\n");
				if (verbose) System.out.print("(*)");
        		
				//resetIsosurfaceNarrowBand(narrowband);
				resetIsosurfaceBoundary();
				fastMarchingReinitialization(false, fullMarching, true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
					// the criterion for being in the narrow band is to have a shortdistance to closest boundaries
					if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) {
						narrowband.addPoint(xyz, mgdmlabels, mgdmfunctions);
						if (mgdmfunctions[0][xyz]>=landmineDist) {
							landmines.set(xyz,true);
						}
						if (Numerics.abs(mgdmfunctions[0][xyz])<1.0) {
							boundarysize++;
						}
					}
				}
			}	
     	}
     	
		
		// end of the evolution: recompute the level sets (disabled for debugging)
		resetIsosurfaceBoundary();
		fastMarchingReinitialization(false, fullMarching, false);
		
        return;
    }
    
  	/** specific forces applied to the level sets (application dependent) */
 	private final void levelsetForces(int xyz, double[] forces) {
    	
		// simple option: rebuild each level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		
		for (int n=0;n<=nmgdm;n++) {
			
			// label
			byte lb = mgdmlabels[n][xyz];
			
			// do the center point first
			if (n==0) phi[CTR] = -mgdmfunctions[0][xyz];
			else  phi[CTR] = 0.0f;
			for (int l=0;l<n;l++) {
				phi[CTR] += mgdmfunctions[l][xyz];
			}
				
			// neighbors
			for (int b=0;b<NGB;b++) {
				int xyzn = xyz + ngbx[b] + ngby[b]*nix + ngbz[b]*nix*niy;

				if (mask[xyzn] && mgdmlabels[0][xyzn]!=EMPTY && mgdmfunctions[0][xyzn]!=UNKNOWN) {
					if (mgdmlabels[0][xyzn]==lb) phi[b] = -mgdmfunctions[0][xyzn];
					else  phi[b] = 0.0f;
					
					for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb && mgdmfunctions[l][xyzn]!=UNKNOWN;l++) {
						phi[b] += mgdmfunctions[l][xyzn];
					}
				} else {
					// filling in values outside the mask?? center value
					phi[b] = phi[CTR];
				}
			}
			
			// first derivatives
			
			Dmx[n] = phi[CTR] - phi[mX];
			Dmy[n] = phi[CTR] - phi[mY];
			Dmz[n] = phi[CTR] - phi[mZ];
			
			Dpx[n] = phi[pX] - phi[CTR];
			Dpy[n] = phi[pY] - phi[CTR];
			Dpz[n] = phi[pZ] - phi[CTR];
			
			D0x = (phi[pX] - phi[mX])/2.0;
			D0y = (phi[pY] - phi[mY])/2.0;
			D0z = (phi[pZ] - phi[mZ])/2.0;
		
			// second derivatives
			Dxx = phi[mX] + phi[pX] - 2.0*phi[CTR];
			Dyy = phi[mY] + phi[pY] - 2.0*phi[CTR];
			Dzz = phi[mZ] + phi[pZ] - 2.0*phi[CTR];
			
			Dxy = (phi[mXmY] + phi[pXpY] - phi[mXpY] - phi[pXmY])/4.0;
			Dyz = (phi[mYmZ] + phi[pYpZ] - phi[mYpZ] - phi[pYmZ])/4.0;
			Dzx = (phi[mZmX] + phi[pZpX] - phi[mZpX] - phi[pZmX])/4.0;
			
			// gradient norm
			SD0x = D0x * D0x;
			SD0y = D0y * D0y;
			SD0z = D0z * D0z;
			GPhi = Math.sqrt(SD0x + SD0y + SD0z);
			
			// mean curvature
			K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
						- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
				
			// gaussian curvature
			G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
						+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));

			// curvature smoothing force:
			if(GPhi > 0.0000001){
				tmp = GPhi*GPhi;
				K = K/(GPhi*tmp);
				G = G/(tmp*tmp);
				tmp = K*K - 2*G;
				if(tmp > 0 ) K = K*G/tmp;
			} else {
				K = 0;
			}
			// regularization
			forces[n] = -smoothweight*stepsize*K*GPhi;
			
			// distance-based atenuation
			if (phi[CTR]<0) distval[n] = 1.0;
			else distval[n] = 1.0/(1.0+phi[CTR]*phi[CTR]/gaindist2);
		}
		
		// find the best target in the neighborhood to build the corresponding force
		// (product of gain and proximity => must check all)
		bestlb = EMPTY;
		bestval = -1.0f;
		done=false;
		for (int n=0;n<ngain && !done;n++) {
			gainlb = bestlabel[n][xyz];
			gainval = 0.5f+0.5f*bestgain[n][xyz];
			if (gainval<=0) done = true;
			for (byte l=0;l<=nmgdm && !done;l++) {
				if (mgdmlabels[l][xyz]==gainlb  && distval[l]*gainval>bestval) {
					bestlb = l;
					bestval = distval[l]*gainval;
				}
			}
		}
		/*
		// (product of gain and proximity => must check all)
		bestlb = EMPTY;
		bestval = -1.0f;
		done=false;
		gainlb = corrected[xyz];
		for (byte l=0;l<=nmgdm && !done;l++) {
			if (mgdmlabels[l][xyz]==gainlb) {
				bestlb = l;
				bestval = distval[l];
				done=true;
			}
		}
		*/
		// second pass for the data force
		if (bestlb!=EMPTY) {
		//if (bestlb!=EMPTY && gainval>0) {	// use only positive forces??
			for (int n=0;n<=nmgdm;n++) {
				byte lb = mgdmlabels[n][xyz];
				if (lb!=EMPTY) {
					
					// upwind scheme
					forces[n] += forceweight*stepsize*bestval
								*(Numerics.max(Dmx[bestlb]+Dpx[bestlb],0)*Dmx[n] + Numerics.min(Dmx[bestlb]+Dpx[bestlb],0)*Dpx[n] 
								 +Numerics.max(Dmy[bestlb]+Dpy[bestlb],0)*Dmy[n] + Numerics.min(Dmy[bestlb]+Dpy[bestlb],0)*Dpy[n] 
								 +Numerics.max(Dmz[bestlb]+Dpz[bestlb],0)*Dmz[n] + Numerics.min(Dmz[bestlb]+Dpz[bestlb],0)*Dpz[n]);
					
					// add a balloon term if prod ~ 0: not useful
				}
			}
		}
		return;
    }

    /**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean homeomorphicLabeling(int xyz, byte lb) {
    	// if we don't check, just exit
    	if (!checkTopology) return true;
    	
		// is the new segmentation homeomorphic ? 
		
		// inside the original object ?
		if (segmentation[xyz]==lb) return true;
		
		// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[xyz+i+j*nix+l*nix*niy]==lb) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		// does it change the topology of the object it modifies ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[xyz+i+j*nix+l*nix*niy]==segmentation[xyz]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = false;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;

		// does it change the topology of a relation between the modified object and its neighbors ?
		Nconfiguration = 0;
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[xyz+i+j*nix+l*nix*niy]!=lb) 
				&& (segmentation[xyz+i+j*nix+l*nix*niy]!=segmentation[xyz]) ) {
				found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[n]) { found = true; break; }
				
				if (!found) {
					lbs[Nconfiguration] = segmentation[xyz+i+j*nix+l*nix*niy];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xyz+i+j*nix+l*nix*niy]==segmentation[xyz])
					|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = false;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}
		for (int n=0;n<Nconfiguration;n++) {
			// in relation with new object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xyz+i+j*nix+l*nix*niy]==lb)
					|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = true;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}

		// triplets
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with previous object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[xyz+i+j*nix+l*nix*niy]==segmentation[xyz])
						|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[n])
						|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = false;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with new object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[xyz+i+j*nix+l*nix*niy]==lb)
						|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[n]) 
						|| (segmentation[xyz+i+j*nix+l*nix*niy]==lbs[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = true;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}

		// else, it works
		return true;
    }
	
	public final void fastMarchingInitialization(boolean narrowBandOnly, boolean almostEverywhere) {
         // initialize the quantities
         for (int xyz = 0; xyz<nix*niy*niz; xyz++) {
         	 // mgdm functions
         	 mgdmlabels[0][xyz] = EMPTY;
         	 for (int n = 1; n<nmgdm; n++) {
            	mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
            mgdmlabels[nmgdm][xyz] = EMPTY;
        }
		
		// computation variables
        byte[] processed = new byte[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;
		boolean done, isprocessed;
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
        	processed[xyz] = 0;
        	// search for boundaries
        	for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				if (segmentation[xyzn]!=segmentation[xyz]) if (mask[xyzn]) {
					// add to the heap
					heap.addValue(mgdmfunctions[0][xyzn],xyzn,segmentation[xyz]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	curdist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			done = false;
			for (int n=0; n<processed[xyz]; n++)
				if (mgdmlabels[n][xyz]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xyz]][xyz] = curdist;
			mgdmlabels[processed[xyz]][xyz] = lb;
			processed[xyz]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					isprocessed = false;
					if (segmentation[xyzn]==lb) isprocessed = true;
					else {
						for (int n=0; n<processed[xyzn]; n++)
							if (mgdmlabels[n][xyzn]==lb) isprocessed = true;
					}
					
					if (!isprocessed) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							for (int n=0; n<processed[xyznb]; n++) if (mask[xyznb]) if (mgdmlabels[n][xyznb]==lb) {
								nbdist[l] = mgdmfunctions[n][xyznb];
								nbflag[l] = true;
							}			
						}
						newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if ( (!narrowBandOnly && !almostEverywhere)
							|| (narrowBandOnly && newdist<=narrowBandDist+extraDist)
							|| (almostEverywhere && (segmentation[xyzn]!=0 || newdist<=narrowBandDist+extraDist) ) ) {
							// add to the heap
							heap.addValue(newdist,xyzn,lb);
						}
					}
				}
			}
		}
		// re-buidld a better mask region
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;

			if (x==0 || x==nix-1 || y==0 || y==niy-1 || z==0 || z==niz-1) mask[xyz] = false;
			else if (segmentation[xyz]==EMPTY || (processed[xyz]==0 && segmentation[xyz]==0)) mask[xyz] = false;
			else mask[xyz] = true;
		}
		
		// to create the MGDM functions, we need to copy the segmentation, forget the last labels
		// and compute differences between distance functions
		if (debug) BasicInfo.displayMessage("transform into MGDM functions\n");		
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			// label permutation
			for (int n=nmgdm;n>0;n--) {
				mgdmlabels[n][xyz] = mgdmlabels[n-1][xyz];
			}
			mgdmlabels[0][xyz] = segmentation[xyz];
			
			// distance function difference
        	for (int n = nmgdm-1; n>0; n--) {
        		mgdmfunctions[n][xyz] = Numerics.max(UNKNOWN, mgdmfunctions[n][xyz]
        														-mgdmfunctions[n-1][xyz]);
			}
        }
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }
     
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingReinitialization(boolean narrowBandOnly, boolean almostEverywhere, boolean stopCounter) {
        // computation variables
        byte[] processed = new byte[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;	
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		

		long start_time = System.currentTimeMillis(); 

        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
        	// mgdm functions : reinit everiywhere
			for (int n = 0; n<nmgdm; n++) {
            	if (n>0) mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
            mgdmlabels[nmgdm][xyz] = EMPTY;
            processed[xyz] = 0;
        	// not needed, should be kept the same
        	//segmentation[xyz] = mgdmlabels[0][xyz];
        	// search for boundaries
        	for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				if (segmentation[xyzn]!=segmentation[xyz]) if (mask[xyzn]) {
					
					// add to the heap with previous value
					heap.addValue(mgdmfunctions[0][xyzn],xyzn,segmentation[xyz]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	curdist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			done = false;
			for (int n=0; n<processed[xyz]; n++)
				if (mgdmlabels[n][xyz]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xyz]][xyz] = curdist;
			mgdmlabels[processed[xyz]][xyz] = lb;
			processed[xyz]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn] && (!stopCounter || counter[xyzn]<2*maxcount) ) {
					// must be in outside the object or its processed neighborhood
					isprocessed = false;
					if (segmentation[xyzn]==lb) isprocessed = true;
					else {
						for (int n=0; n<processed[xyzn]; n++)
							if (mgdmlabels[n][xyzn]==lb) isprocessed = true;
					}
					
					if (!isprocessed) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							for (int n=0; n<processed[xyznb]; n++) if (mask[xyznb]) if (mgdmlabels[n][xyznb]==lb) {
								nbdist[l] = mgdmfunctions[n][xyznb];
								nbflag[l] = true;
							}			
						}
						newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if ( (!narrowBandOnly && !almostEverywhere)
							|| (narrowBandOnly && newdist<=narrowBandDist+extraDist)
							|| (almostEverywhere && (segmentation[xyzn]!=0 || newdist<=narrowBandDist+extraDist) ) ) {
							// add to the heap
							heap.addValue(newdist,xyzn,lb);
						}
					}
				}
			}
		}
		// to create the MGDM functions, we need to copy the segmentation, forget the last labels
		// and compute differences between distance functions
		if (debug) BasicInfo.displayMessage("transform into MGDM functions\n");		
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			// label permutation
			for (int n=nmgdm;n>0;n--) {
				mgdmlabels[n][xyz] = mgdmlabels[n-1][xyz];
			}
			mgdmlabels[0][xyz] = segmentation[xyz];
			
			// distance function difference
        	for (int n = nmgdm-1; n>0; n--) {
        		mgdmfunctions[n][xyz] = Numerics.max(UNKNOWN, mgdmfunctions[n][xyz]
        														-mgdmfunctions[n-1][xyz]);
        	}
        }
		if (debug) BasicInfo.displayMessage("done (time: " + (System.currentTimeMillis()-start_time)+")\n"); 

       return;
    }

	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * 
     */
    public final float minimumMarchingDistance(float[] val, boolean[] flag) {

        // s = a + b +c; s2 = a*a + b*b +c*c
        s = 0;
        s2 = 0;
        count = 0;

        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.min(val[n], val[n+1]); // Take the smaller one if both are processed
				s += tmp;
				s2 += tmp*tmp;
				count++;
			} else if (flag[n]) {
				s += val[n]; // Else, take the processed one
				s2 += val[n]*val[n];
				count++;
			} else if (flag[n+1]) {
				s += val[n+1];
				s2 += val[n+1]*val[n+1];
				count++;
			}
		}
         // count must be greater than zero since there must be at least one processed pt in the neighbors
        
        tmp = (s+Math.sqrt( (s*s-count*(s2-1.0f))))/count;

        // The larger root
        return (float)tmp;
    }
	/**
     * the isosurface distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * (the input values are all positive, the flags are true only if the isosurface crosses)
     */
    public final float isoSurfaceDistance(double cur, float[] val, boolean[] flag) {
    	
    	if (cur==0) return 0;
    	
        s = 0;
        dist = 0;
        
        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.max(val[n], val[n+1]); // Take the largest distance (aka closest to current point) if both are across the boundariy
				s = cur/(cur+tmp);
				dist += 1.0/(s*s);
			} else if (flag[n]) {
				s = cur/(cur+val[n]); // Else, take the boundariy point
				dist += 1.0/(s*s);
			} else if (flag[n+1]) {
				s = cur/(cur+val[n+1]);
				dist += 1.0/(s*s);
			}
		}
		// triangular (tetrahedral?) relationship of height in right triangles gives correct distance
        tmp = Math.sqrt(1.0/dist);

        // The larger root
        return (float)tmp;
    }
    /** 
    *   isosurface distance re-initialization at the boundariy
    */
    private final void resetIsosurfaceBoundary() {
    	if (debug) System.out.print("fast marching evolution: iso-surface reinit\n");

    	float[] nbdist = new float[6];
    	boolean[] nbflag = new boolean[6];
    	boolean boundary;
    	float[] tmp = new float[nix*niy*niz];
    	boolean[] processed = new boolean[nix*niy*niz];
    	
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			
			boundary = false;
			for (int l=0; l<6; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				int xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
				if (segmentation[xyznb]!=segmentation[xyz] && mask[xyznb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(mgdmfunctions[0][xyznb]);
					nbflag[l] = true;
					boundary = true;
				}
			}
			if (boundary) {
				tmp[xyz] = isoSurfaceDistance(mgdmfunctions[0][xyz], nbdist, nbflag);
				processed[xyz] = true;
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) {
			if (processed[xyz]) mgdmfunctions[0][xyz] = tmp[xyz];
			else mgdmfunctions[0][xyz] = UNKNOWN;
		}
			
        return;
    }
    
	private static class NarrowBand {
		public int[] id;
		public byte[][] labels;
		public float[][] functions;
		public int currentsize;
		public int capacity;
		public int update;
		
		/** init: create an array of a given size for efficient storage */
		public NarrowBand(int size, int increase) {
			capacity = size;
			update = increase;
			currentsize = 0;
			
			id = new int[capacity];
			labels = new byte[nmgdm+1][capacity];
			functions = new float[nmgdm][capacity];
		}
		
		public void finalize() {
			capacity = -1;
			id = null;
			labels = null;
			functions = null;
		}
		
		public final void addPoint(int xyz, byte[][] mgdmlabels, float[][] mgdmfn) {
			// check for size
			if (currentsize>=capacity-1) {
				capacity += update;
				
				int[] oldid = id;
				byte[][] oldlabels = labels;
				float[][] oldfunctions = functions;
				
				id = new int[capacity];
				labels = new byte[nmgdm+1][capacity];
				functions = new float[nmgdm][capacity];
				
				for (int n=0;n<currentsize;n++) {
					id[n] = oldid[n];
					for (int l=0;l<nmgdm;l++) {
						labels[l][n] = oldlabels[l][n];
						functions[l][n] = oldfunctions[l][n];
					}
					labels[nmgdm][n] = oldlabels[nmgdm][n];
				}
				
				oldid = null; oldlabels = null; oldfunctions = null;
			}
			
			// add the new point (use the MGDM variables)
			id[currentsize] = xyz;
			for (int l=0;l<nmgdm;l++) {
				labels[l][currentsize] = mgdmlabels[l][xyz];
				functions[l][currentsize] = mgdmfn[l][xyz];
			}
			labels[nmgdm][currentsize] = mgdmlabels[nmgdm][xyz];
			currentsize++;
		}
		
		public void reset() {
			currentsize = 0;
		}
	}

}

