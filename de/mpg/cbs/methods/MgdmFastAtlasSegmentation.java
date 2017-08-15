package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm uses the MGDM framework to evolve a labeling
 *	according to internal (curvature) and external (vector field) forces
 *	note: for now, we are assuming isotropic voxels 	
 *
 *	@version    August 2011
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class MgdmFastAtlasSegmentation {
	
	
	// image data
	private 	float[][]		image;  			// original images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String[]		modality;			// image modality / contrast
	private		float[]			imrange;			// image intensity range (robust estimate)
	private		int				nc;					// number of channels

	// atlas parameters
	private		SimpleShapeAtlas	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	byte[]			objLabel;			// label values in the original image
	private 	int				nax,nay,naz;   			// atlas dimensions
	private 	float			rax,ray,raz;   			// atlas resolutions
	    
	// gain functions
	private		float[][]		bestgain;			// best gain function
	private		byte[][]		bestlabel;			// corresponding labels
	private		float[][]		bestgainHD;			// best gain function (high-res version)
	private		byte[][]		bestlabelHD;		// corresponding labels (high-res version)
	private static	byte 	  	ngain;				// total number of gain functions
	private		float			pvsize;				// window size
	private		float			sigmaI, sigmaS;
	private		float[] 		sigmaN;
	private		float			factor;
	private		float			threshold1;
	private		float			threshold2;
	private		float			slope;
	private		int				ngbsize;
	private		float			alpha;
	
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
	private	double		smoothweight, forceweight, divweight;
	private	double		K0;
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
	//double[] Dmx, Dmy, Dmz, Dpx, Dpy, Dpz;
	double[] Dmx,Dmy,Dmz, Dpx, Dpy, Dpz;
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
	
	// convenience labels
	public static final byte	X=0;
	public static final byte	Y=1;
	public static final byte	Z=2;
	
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
	public MgdmFastAtlasSegmentation(float[][] img_, String[] mod_, float[] rng_, int nc_,
									int nix_, int niy_, int niz_, float rix_, float riy_, float riz_,
									SimpleShapeAtlas atlas_,
									int nmgdm_, int ngain_,
									float fw_, float sw_, float dw_, float k0_, float gd_, 
									String connectivityType_) {
	
		image = img_;
		imrange = rng_;
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
		
		sigmaI = 0.2f;
		sigmaN = new float[nc];
		for (int c=0;c<nc;c++) sigmaN[c] = 0.02f;
		sigmaS = 0.1f;
		
		threshold1 = 0.0f;
		threshold2 = 0.0f;
		slope	= 0.9f;
		ngbsize = 0;
		alpha = 0.1f;
		pvsize = 0.5f;
		
		nmgdm = (byte)nmgdm_;
		ngain = (byte)ngain_;
						
		forceweight = fw_;
		smoothweight = sw_;
		divweight = dw_;
		
		K0 = k0_;
		
		gaindist2 = gd_*gd_;
		
		if (verbose) System.out.print("MGDMA parameters: "+fw_+" (force), "+sw_+" (smoothing), "+dw_+" (div), "+k0_+" (K0), "+(rax/rix)+" (scale), "+gd_+" (distance)\n");
		
		smoothfactor = atlas.getRegularizationFactor();
		
		objLabel = atlas.getLabels();

		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nax, -nax, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nax*nay, -nax*nay};

		// init all the arrays in atlas space
		try {
			segmentation = new byte[nax*nay*naz];	
			mask = new boolean[nax*nay*naz];	
			counter = new short[nax*nay*naz];	
			mgdmfunctions = new float[nmgdm][nax*nay*naz];
			mgdmlabels = new byte[nmgdm+1][nax*nay*naz];	
			phi = new double[27];
			/*
			Dx = new double[nmgdm+1];
			Dy = new double[nmgdm+1];
			Dz = new double[nmgdm+1];
			*/
			Dmx = new double[nmgdm+1];
			Dmy = new double[nmgdm+1];
			Dmz = new double[nmgdm+1];
			Dpx = new double[nmgdm+1];
			Dpy = new double[nmgdm+1];
			Dpz = new double[nmgdm+1];
			
			distval = new double[nmgdm+1];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nax*nay+nay*naz+naz*nax, BinaryHeap2D.MINTREE);
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
					//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nax; x++) for (int y=0; y<nay; y++) for (int z = 0; z<naz; z++) {
			if (x>1 && x<nax-2 && y>1 && y<nay-2 && z>1 && z<naz-2) mask[x+nax*y+nax*nay*z] = true;
			else mask[x+nax*y+nax*nay*z] = false;
		}
		// init segmentation from atlas
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			int xyz = x + nax*y + nax*nay*z;
			
			int init = atlas.getTemplate()[xyz];
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (atlas.getLabels()[n]==init) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
			
			counter[xyz] = 0;
		}
		
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
		lut.finalize();
		lut = null;
		System.gc();
	}

	public final float[][] getFunctions() { return mgdmfunctions; }
	
	public final byte[][] getLabels() { return mgdmlabels; }
	
	public final byte[] getSegmentation() { return segmentation; }
    
	public final float[][] getBestGainFunctionHD() { return bestgainHD; }
	
	public final byte[][] getBestGainLabelHD() { return bestlabelHD; }
    
	public final void setWeights(float fw_, float sw_) {
		forceweight = fw_;
		smoothweight = sw_;
		
		if (verbose) System.out.print("MGDMA forces: "+fw_+" (force), "+sw_+" (smoothing)\n");
	}
	
	/**
	 *	get a final segmentation
	 */
    public final byte[] exportSegmentationByte() {
    	byte[] seg = new byte[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			int lbl = ImageInterpolation.nearestNeighborInterpolation(mgdmlabels[0],(byte)-1,XA[0],XA[1],XA[2],nax,nay,naz);
    		if (lbl>-1) {
				seg[xyz] = (byte)objLabel[lbl];
			}
    	}
    	return seg;
    }
    public final byte[] getLabeledSegmentation() {
    	byte[] seg = new byte[nax*nay*naz];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			int xyz = x + nax*y + nax*nay*z;
			
			int lbl = segmentation[xyz];
    		if (lbl>-1) {
				seg[xyz] = (byte)objLabel[lbl];
			} else {
				seg[xyz] = 0;
			}
    	}
    	return seg;
    }
    public final float[] exportAtlasSegmentation() {
    	float[] seg = new float[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			int xa = Numerics.round(XA[X]);
			int ya = Numerics.round(XA[Y]);
			int za = Numerics.round(XA[Z]);
			if (xa>=0 && xa<nax && ya>=0 && ya<nay && za>=0 && za<naz) {
				int xyza = xa + nax*ya + nax*nay*za;
			
				int lbl = mgdmlabels[0][xyza];
				if (lbl>-1) {
					seg[xyz] = (float)objLabel[lbl];
				}
			} else {
				seg[xyz] = 0;
			}
    	}
    	return seg;
    }
    public final float[] exportAtlasBestGain() {
    	float[] seg = new float[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			int lbl = ImageInterpolation.nearestNeighborInterpolation(bestlabel[0],(byte)-1,XA[0],XA[1],XA[2],nax,nay,naz);
    		if (lbl>-1) {
				seg[xyz] = (float)objLabel[lbl];
			}
    	}
    	return seg;
    }
    public final int[][][] exportImageBestGain() {
    	int[][][] seg = new int[nix][niy][niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			int lbl = bestlabelHD[0][xyz];
    		if (lbl>-1) {
				seg[x][y][z] = objLabel[lbl];
			}
    	}
    	return seg;
    }
    public final float[] exportAtlasGainLabels(int n) {
    	float[] seg = new float[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			int lbl = ImageInterpolation.nearestNeighborInterpolation(bestlabel[n],(byte)-1,XA[0],XA[1],XA[2],nax,nay,naz);
    		if (lbl>-1) {
				seg[xyz] = (float)objLabel[lbl];
			}
    	}
    	return seg;
    }
    public final float[] exportAtlasGainFunctions(int n) {
    	float[] seg = new float[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			seg[xyz] = ImageInterpolation.nearestNeighborInterpolation(bestgain[n],(byte)-1,XA[0],XA[1],XA[2],nax,nay,naz);
    	}
    	return seg;
    }
    public final float[] exportAtlasObjectGain(int n) {
    	float[] seg = new float[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			for (int g=0;g<=ngain;g++) {
				int lb = ImageInterpolation.nearestNeighborInterpolation(bestlabel[g],(byte)-1,XA[0],XA[1],XA[2],nax,nay,naz);
				if (lb==n) {
					seg[xyz] = ImageInterpolation.nearestNeighborInterpolation(bestgain[g],(byte)0,XA[0],XA[1],XA[2],nax,nay,naz);
				}
			}
    	}
    	return seg;
    }
    public final short[] exportFrozenPointCounter() {
    	short[] seg = new short[nix*niy*niz];
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x + nix*y + nix*niy*z;
			
			float[] XA = new float[3];
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			seg[xyz] = ImageInterpolation.nearestNeighborInterpolation(counter,(short)0,XA[0],XA[1],XA[2],nax,nay,naz);
    	}
    	return seg;
    }
   
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final void computeAtlasBestGainFunction() {
    	bestgainHD = new float[ngain+1][nix*niy*niz];
    	bestlabelHD = new byte[ngain+1][nix*niy*niz];
    	
    	byte mp2rage7T=-1;
    	byte t1map7T = -1;
    	byte mp2rage3T=-1;
    	byte t1map3T = -1;
    	byte mprage3T=-1;
    	byte t2sw7T=-1;
    	byte qsm7T=-1;
    	byte pv=-1;
    	byte filters=-1;
    	byte hcpt1w=-1;
    	byte hcpt2w=-1;
    	byte normmprage=-1;
    	for (byte n=0;n<nc;n++) {
			if (modality[n].equals("Iso") || modality[n].equals("FlatImg") || modality[n].equals("Mp2rage7T") || modality[n].equals("MP2RAGE7T")) mp2rage7T=n;
			else if (modality[n].equals("T1map") || modality[n].equals("T1map7T") || modality[n].equals("T1MAP7T")) t1map7T=n;
			else if (modality[n].equals("T1w") || modality[n].equals("T1_SPGR") || modality[n].equals("T1_MPRAGE") || modality[n].equals("Mprage3T") || modality[n].equals("MPRAGE3T")) mprage3T=n;
			else if (modality[n].equals("T2SW7T") || modality[n].equals("Flash7T")) t2sw7T=n;
			else if (modality[n].equals("QSM7T") || modality[n].equals("qsm7T")) qsm7T=n;
			else if (modality[n].equals("Mp2rage3T") || modality[n].equals("MP2RAGE3T")) mp2rage3T=n;
			else if (modality[n].equals("T1map3T") || modality[n].equals("T1MAP3T")) t1map3T=n;
			else if (modality[n].equals("PV") || modality[n].equals("PVDURA")) pv=n;
			else if (modality[n].equals("Filters") ||  modality[n].equals("FILTER")) filters=n;
			else if (modality[n].equals("HCPT1w") ||  modality[n].equals("HCPT1W")) hcpt1w=n;
			else if (modality[n].equals("HCPT2w") ||  modality[n].equals("HCPT2W")) hcpt2w=n;
			else if (modality[n].equals("NORMMPRAGE") ||  modality[n].equals("NormMPRAGE")) normmprage=n;
		}
		if (verbose) {
			System.out.println("Modalities used in segmentation: ");
			//System.out.print("(");
			//for (byte n=0;n<nc-1;n++) System.out.print(modality[n]+", ");
			//System.out.println(modality[nc-1]+")->");
			if (mp2rage7T!=-1) System.out.println("MP2RAGE7T ("+mp2rage7T+")");
			if (t1map7T!=-1) System.out.println("T1MAP7T ("+t1map7T+")");
			if (mp2rage3T!=-1) System.out.println("MP2RAGE3T ("+mp2rage3T+")");
			if (t1map3T!=-1) System.out.println("T1MAP3T ("+t1map3T+")");
			if (mprage3T!=-1) System.out.println("MPRAGE3T ("+mprage3T+")");
			if (t2sw7T!=-1) System.out.println("T2SW7T ("+t2sw7T+")");
			if (qsm7T!=-1) System.out.println("QSM7T ("+qsm7T+")");
			if (pv!=-1) System.out.println("PVDURA ("+pv+")");
			if (filters!=-1) System.out.println("FILTERS ("+filters+")");
			if (hcpt1w!=-1) System.out.println("HCPT1W ("+hcpt1w+")");
			if (hcpt2w!=-1) System.out.println("HCPT2W ("+hcpt2w+")");
			if (normmprage!=-1) System.out.println("NORMMPRAGE ("+normmprage+")");
		}

    	// compute the local priors
    	float shape;
    	float[] intens = new float[nobj];
    	//float[] pvscore = new float[nobj];
    	float[] contrast = new float[3];
    	float val, diff, mindiff, maxdiff, proba;
    	float pvproba, pvval;
    	float filterproba, filterval;
    	float shapesum, shapemax;
    	float img;
    	int xyza, xyzi, xyzb;
    	int xa, ya, za;
    	float[] XA = new float[3];

    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		xyzi = x + nix*y + nix*niy*z;
        	atlas.imageToShapeCoordinates(XA, x, y, z);
        	xa = Numerics.round(XA[X]);
        	ya = Numerics.round(XA[Y]);
        	za = Numerics.round(XA[Z]);
        	
			for (int n=0;n<nobj;n++) {
				intens[n] = -1.0f;
			}
			
			// just look at the given scale
					
			// first: check for mask
			boolean isMasked=true;
			if (ImageInterpolation.nearestNeighborInterpolation(segmentation,(byte)0,XA[X],XA[Y],XA[Z],nax,nay,naz)!=0) isMasked = false;
			if (ImageInterpolation.linearInterpolation(atlas.getShape(0),1.0f,XA[X],XA[Y],XA[Z],nax,nay,naz)<0.9f) isMasked = false;
			if (mp2rage7T>-1 && image[mp2rage7T][xyzi]!=0) isMasked = false;
			if (t1map7T>-1 && image[t1map7T][xyzi]!=0) isMasked = false;
			if (mp2rage3T>-1 && image[mp2rage3T][xyzi]!=0) isMasked = false;
			if (t1map3T>-1 && image[t1map3T][xyzi]!=0) isMasked = false;
			if (mprage3T>-1 && image[mprage3T][xyzi]!=0) isMasked = false;
			if (t2sw7T>-1 && image[t2sw7T][xyzi]!=0) isMasked = false;
			if (qsm7T>-1 && image[qsm7T][xyzi]!=0) isMasked = false;
			if (pv>-1 && image[pv][xyzi]!=0) isMasked = false;
			if (filters>-1 && image[filters][xyzi]!=0) isMasked = false;
			if (hcpt1w>-1 && image[hcpt1w][xyzi]!=0) isMasked = false;
			if (hcpt2w>-1 && image[hcpt2w][xyzi]!=0) isMasked = false;
			if (normmprage>-1 && image[normmprage][xyzi]!=0) isMasked = false;
			// disable this?
			//isMasked=false;
			
			if (isMasked) {
				bestlabelHD[0][xyzi] = 0;
				bestgainHD[0][xyzi] = 1.0f;
				for (int n=1;n<ngain+1;n++) {
					bestlabelHD[n][xyzi] = EMPTY;
					bestgainHD[n][xyzi] = -1.0f;
				}				
			} else {
				for (int n=0;n<nobj;n++) {
					// compute the exisiting contrasts
					if (mp2rage7T>-1) {
						contrast[mp2rage7T] = -1.0f;
						if  (image[mp2rage7T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.MP2RAGE7T)[n].length;t+=3) {
								diff = Numerics.abs(image[mp2rage7T][xyzi]/imrange[mp2rage7T]-atlas.getMap(atlas.MP2RAGE7T)[n][t])/atlas.getMap(atlas.MP2RAGE7T)[n][t+1];
								//contrast[iso] = Numerics.max(contrast[iso], atlas.getMap(atlas.ISO)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.MP2RAGE7T)[n][t+2]>0)
									contrast[mp2rage7T] = Numerics.max(contrast[mp2rage7T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[mp2rage7T] = 1.0f;
					}
					if (t1map7T>-1) {
						contrast[t1map7T] = -1.0f;
						if (image[t1map7T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.T1MAP7T)[n].length;t+=3) {
								diff = Numerics.abs(image[t1map7T][xyzi]/imrange[t1map7T]-atlas.getMap(atlas.T1MAP7T)[n][t])/atlas.getMap(atlas.T1MAP7T)[n][t+1];
								//contrast[t1map] = Numerics.max(contrast[t1map], atlas.getMap(atlas.T1M)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.T1MAP7T)[n][t+2]>0)
									contrast[t1map7T] = Numerics.max(contrast[t1map7T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[t1map7T] = 1.0f;
					}
					if (mp2rage3T>-1) {
						contrast[mp2rage3T] = -1.0f;
						if (image[mp2rage3T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.MP2RAGE3T)[n].length;t+=3) {
								diff = Numerics.abs(image[mp2rage3T][xyzi]/imrange[mp2rage3T]-atlas.getMap(atlas.MP2RAGE3T)[n][t])/atlas.getMap(atlas.MP2RAGE3T)[n][t+1];
								//contrast[iso] = Numerics.max(contrast[iso], atlas.getMap(atlas.ISO)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.MP2RAGE3T)[n][t+2]>0)
									contrast[mp2rage3T] = Numerics.max(contrast[mp2rage3T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[mp2rage3T] = 1.0f;
					}
					if (t1map3T>-1) {
						contrast[t1map3T] = -1.0f;
						if (image[t1map3T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.T1MAP3T)[n].length;t+=3) {
								diff = Numerics.abs(image[t1map3T][xyzi]/imrange[t1map3T]-atlas.getMap(atlas.T1MAP3T)[n][t])/atlas.getMap(atlas.T1MAP3T)[n][t+1];
								//contrast[t1map] = Numerics.max(contrast[t1map], atlas.getMap(atlas.T1M)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.T1MAP3T)[n][t+2]>0)
									contrast[t1map3T] = Numerics.max(contrast[t1map3T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[t1map3T] = 1.0f;
					}
					if (mprage3T>-1) {
						contrast[mprage3T] = -1.0f;
						if (image[mprage3T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.MPRAGE3T)[n].length;t+=3) {
								diff = Numerics.abs(image[mprage3T][xyzi]/imrange[mprage3T]-atlas.getMap(atlas.MPRAGE3T)[n][t])/atlas.getMap(atlas.MPRAGE3T)[n][t+1];
								//contrast[t1w] = Numerics.max(contrast[t1w], atlas.getMap(atlas.T1w)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.MPRAGE3T)[n][t+2]>0)
									contrast[mprage3T] = Numerics.max(contrast[mprage3T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[mprage3T] = 1.0f;
					}
					if (t2sw7T>-1) {
						contrast[t2sw7T] = -1.0f;
						if (image[t2sw7T][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.T2SW7T)[n].length;t+=3) {
								diff = Numerics.abs(image[t2sw7T][xyzi]/imrange[t2sw7T]-atlas.getMap(atlas.T2SW7T)[n][t])/atlas.getMap(atlas.T2SW7T)[n][t+1];
								//contrast[flair] = Numerics.max(contrast[flair], atlas.getMap(atlas.FLAIR)[n][t+2]*(1.0f - diff*diff)/(1.0f + diff*diff) );
								if (atlas.getMap(atlas.T2SW7T)[n][t+2]>0)
									contrast[t2sw7T] = Numerics.max(contrast[t2sw7T], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[t2sw7T] = 1.0f;
					}
					if (hcpt1w>-1) {
						contrast[hcpt1w] = -1.0f;
						if (image[hcpt1w][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.HCPT1W)[n].length;t+=3) {
								diff = Numerics.abs(image[hcpt1w][xyzi]/imrange[hcpt1w]-atlas.getMap(atlas.HCPT1W)[n][t])/atlas.getMap(atlas.HCPT1W)[n][t+1];
								if (atlas.getMap(atlas.HCPT1W)[n][t+2]>0)
									contrast[hcpt1w] = Numerics.max(contrast[hcpt1w], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[hcpt1w] = 1.0f;
					}
					if (hcpt2w>-1) {
						contrast[hcpt2w] = -1.0f;
						if (image[hcpt2w][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.HCPT2W)[n].length;t+=3) {
								diff = Numerics.abs(image[hcpt2w][xyzi]/imrange[hcpt2w]-atlas.getMap(atlas.HCPT2W)[n][t])/atlas.getMap(atlas.HCPT2W)[n][t+1];
								if (atlas.getMap(atlas.HCPT2W)[n][t+2]>0)
									contrast[hcpt2w] = Numerics.max(contrast[hcpt2w], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[hcpt2w] = 1.0f;
					}
					if (normmprage>-1) {
						contrast[normmprage] = -1.0f;
						if (image[normmprage][xyzi]!=0) {
							for (int t=0;t<atlas.getMap(atlas.NORMMPRAGE)[n].length;t+=3) {
								diff = Numerics.abs(image[normmprage][xyzi]/imrange[normmprage]-atlas.getMap(atlas.NORMMPRAGE)[n][t])/atlas.getMap(atlas.NORMMPRAGE)[n][t+1];
								if (atlas.getMap(atlas.NORMMPRAGE)[n][t+2]>0)
									contrast[normmprage] = Numerics.max(contrast[normmprage], (1.0f - diff*diff)/(1.0f + diff*diff) );
							}
						} else if (n==0) contrast[normmprage] = 1.0f;
					}
					/*
					if (t1pv>-1) {
						proba = Numerics.min(1.0f,Numerics.abs(image[t1pv][xyzi]/imrange[t1pv]));
						contrast[t1pv] = Numerics.max(-1.0f, (proba-atlas.getMap(atlas.pT1)[n][1])
															/(proba + atlas.getMap(atlas.pT1)[n][1]*(1.0f-2.0f*proba) ) );
					}
					*/
					pvproba = 0.0f;
					pvval = 0.0f;
					if (pv>-1) {
						if (image[pv][xyzi]<0) {
							pvproba = Numerics.min(1.0f, -image[pv][xyzi]/imrange[pv]);
							pvval = atlas.getMap(atlas.PVDURA)[n][0];
						} else {
							pvproba = Numerics.min(1.0f, image[pv][xyzi]/imrange[pv]);
							pvval = atlas.getMap(atlas.PVDURA)[n][1];
						}
					}
					filterproba = 0.0f;
					filterval = 0.0f;
					if (filters>-1) {
						if (image[filters][xyzi]>=0 && image[filters][xyzi]<=1) {
							// first level
							filterproba = image[filters][xyzi];
							filterval = atlas.getMap(atlas.FILTER)[n][0];
						} else if (image[filters][xyzi]>=2 && image[filters][xyzi]<=3) {
							// second level
							filterproba = image[filters][xyzi]-2.0f;
							filterval = atlas.getMap(atlas.FILTER)[n][1];
						} else if (image[filters][xyzi]>=4 && image[filters][xyzi]<=5) {
							// third level
							filterproba = image[filters][xyzi]-4.0f;
							filterval = atlas.getMap(atlas.FILTER)[n][2];
						}
					}
					
					// combine the intensities and probabilities
					intens[n] = 1.0f;
					
					if (mp2rage7T>-1) intens[n] = Numerics.min(intens[n], contrast[mp2rage7T]);
					if (t1map7T>-1) intens[n] = Numerics.min(intens[n], contrast[t1map7T]);
					if (mp2rage3T>-1) intens[n] = Numerics.min(intens[n], contrast[mp2rage3T]);
					if (t1map3T>-1) intens[n] = Numerics.min(intens[n], contrast[t1map3T]);
					if (mprage3T>-1) intens[n] = Numerics.min(intens[n], contrast[mprage3T]);
					if (t2sw7T>-1) intens[n] = Numerics.min(intens[n], contrast[t2sw7T]);
					if (qsm7T>-1) intens[n] = Numerics.min(intens[n], contrast[qsm7T]);
					if (hcpt1w>-1) intens[n] = Numerics.min(intens[n], contrast[hcpt1w]);
					if (hcpt2w>-1) intens[n] = Numerics.min(intens[n], contrast[hcpt2w]);
					if (normmprage>-1) intens[n] = Numerics.min(intens[n], contrast[normmprage]);
					
					if (pv>-1) if (pvval!=0) intens[n] = pvproba*pvval + (1.0f-pvproba)*intens[n];
					if (filters>-1) if (filterval!=0) intens[n] = filterproba*filterval + (1.0f-filterproba)*intens[n];
				}
				for (int n=0;n<nobj;n++) {			
					if (n==0) val = ImageInterpolation.linearInterpolation(atlas.getShape(n),1.0f,XA[0],XA[1],XA[2],nax,nay,naz);
					else val = ImageInterpolation.linearInterpolation(atlas.getShape(n),0.0f,XA[0],XA[1],XA[2],nax,nay,naz);
					shape = (val - sigmaS)/(val + sigmaS - 2.0f*sigmaS*val);
					
					// use gain composition
					float den = (1.01f+intens[n]) + (1.01f+shape);
					float ws = Numerics.bounded((1.01f+intens[n])/den,0,1);
					float wi = Numerics.bounded((1.01f+shape)/den,0,1);
					intens[n] = wi*intens[n] + ws*shape;
				}
				// keep only the ngain best
				for (int n=0;n<ngain+1;n++) {
					byte nmax=0;
					// careful here: if (almost) all values are -1 to start, you have duplicate values
					for (byte m=1;m<nobj;m++) if (intens[m]>intens[nmax]) {
						nmax = m;
					}
					bestlabelHD[n][xyzi] = nmax;
					bestgainHD[n][xyzi] = intens[nmax];
					intens[nmax] = -2.0f;
				}
			} 
		}
		
		// second pass: compute all the functions, then get best
		float[][] allgains = new float[nobj][nax*nay*naz];
		float[] warpcount = new float[nax*nay*naz];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyza = x + nax*y + nax*nay*z;
    		for (int n=0;n<nobj;n++) allgains[n][xyza] = 0.0f;
    		warpcount[xyza] = 0;
    	}	
    	
    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			xyzi = x + nix*y + nix*niy*z;
			atlas.imageToShapeCoordinates(XA, x, y, z);
			xa = Numerics.round(XA[X]);
			ya = Numerics.round(XA[Y]);
			za = Numerics.round(XA[Z]);
			if (xa>=0 && xa<nax && ya>=0 && ya<nay && za>=0 && za<naz) {
				xyza = xa + nax*ya + nax*nay*za;
				
				warpcount[xyza]++;
				for (int l=0;l<ngain+1;l++) {
					if (bestlabelHD[l][xyzi]!=UNKNOWN) {
						allgains[bestlabelHD[l][xyzi]][xyza] += 0.5f*(1.0f+bestgainHD[l][xyzi]);
					}
				}
			}
		}
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyza = x + nax*y + nax*nay*z;
    		if (warpcount[xyza]>0) {
    			for (int n=0;n<nobj;n++) allgains[n][xyza] /= warpcount[xyza];
    		}
    	}
    	bestgain = new float[ngain+1][nax*nay*naz];
    	bestlabel = new byte[ngain+1][nax*nay*naz];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyza = x + nax*y + nax*nay*z;
    		if (warpcount[xyza]>0) {
				for (int n=0;n<ngain+1;n++) {
					byte nmax=0;
					
					for (byte m=1;m<nobj;m++) if (allgains[m][xyza]>allgains[nmax][xyza]) {
						nmax = m;
					}
					bestlabel[n][xyza] = nmax;
					bestgain[n][xyza] = 2.0f*allgains[nmax][xyza]-1.0f;
					allgains[nmax][xyza] = -1.0f;
				}
			} else {
				for (int n=0;n<ngain+1;n++) {
					bestgain[n][xyza] = -1.0f;
					bestlabel[n][xyza] = -1;
				}
			}
     	}
		allgains = null;
		warpcount = null;
		
		/* method above is better
		// second pass for the atlas version
		bestgain = new float[ngain+1][nax*nay*naz];
    	bestlabel = new byte[ngain+1][nax*nay*naz];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyza = x + nax*y + nax*nay*z;
    		for (int n=0;n<ngain+1;n++) {
    			bestgain[n][xyza] = -1.0f;
    			bestlabel[n][xyza] = -1;
    		}
     	}
    	
    	for (int n=0;n<ngain+1;n++) {
			// from the best to the last
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				xyzi = x + nix*y + nix*niy*z;
				atlas.imageToShapeCoordinates(XA, x, y, z);
				xa = Numerics.round(XA[X]);
				ya = Numerics.round(XA[Y]);
				za = Numerics.round(XA[Z]);
				if (xa>=0 && xa<nax && ya>=0 && ya<nay && za>=0 && za<naz) {
					xyza = xa + nax*ya + nax*nay*za;
					
					for (int l=0;l<ngain+1;l++) {
						if (bestgainHD[l][xyzi]>bestgain[n][xyza]) {
							// check if this label is already set
							boolean found=false;
							for (int m=0;m<n;m++) if (bestlabel[m][xyza]==bestlabelHD[l][xyzi]) {
								found=true;
							}
							if (!found) {
								// replace by the new one
								bestlabel[n][xyza] = bestlabelHD[l][xyzi];
								bestgain[n][xyza] = bestgainHD[l][xyzi];
							}
						}
					}
					
					// also set the mask for the background
					if (segmentation[xyza]==0 && bestlabel[0][xyza]==0 && bestgain[0][xyza]==1 && bestgain[1][xyza]==-1) {
						mask[xyza] = false;
					}
				}
			}
		} 
		*/
    }
    
    // probability diffusion
    public final void diffuseBestGainFunctions(int iter, float scale, float factor) {
    	
    	//mix with the neighbors?
    	float[] map = new float[nix*niy*niz];
    	float[] orig = new float[nix*niy*niz];
    	float[][] newgainHD = new float[ngain+1][nix*niy*niz];
    	byte[][] newlabelHD = new byte[ngain+1][nix*niy*niz];
    	boolean[] diffmask = new boolean[nix*niy*niz];
		for (int m=0;m<ngain+1;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			newlabelHD[m][xyzi] = EMPTY;
		}
    	
    	for (byte n=0;n<nobj;n++) {	
    		BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			// get the gain ; normalize
			for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
				float sum = 0.0f;
				map[xyzi] = 0.0f;
				for (int m=0;m<ngain+1;m++) {
					if (bestlabelHD[m][xyzi]==n) map[xyzi] = 0.5f+0.5f*bestgainHD[m][xyzi];
					//sum += 0.5f+0.5f*bestgainHD[m][xyzi];
				}
				
				// normalize over sum proba?
				//map[xyzi] /= sum;
				
				orig[xyzi] = map[xyzi];
				if (orig[xyzi]==0 || orig[xyzi]==1) diffmask[xyzi] = false;
				else diffmask[xyzi] = true;
			}
			// propagate the values : diffusion
			float maxdiff = 1.0f;
			for (int t=0;t<iter && maxdiff>0.05f;t++) {
				BasicInfo.displayMessage(".");
				maxdiff = 0.0f;
				for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
					int xyzi = x+nix*y+nix*niy*z;
					if (diffmask[xyzi]) {
						float den = diffusionWeightFunction(map[xyzi],orig[xyzi],scale);
						float num = den*orig[xyzi];
						float weight;
						float prev = map[xyzi];
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi-1],scale);
						num += weight*map[xyzi-1];
						den += weight;
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi+1],scale);
						num += weight*map[xyzi+1];
						den += weight;
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi-nix],scale);
						num += weight*map[xyzi-nix];
						den += weight;
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi+nix],scale);
						num += weight*map[xyzi+nix];
						den += weight;
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi-nix*niy],scale);
						num += weight*map[xyzi-nix*niy];
						den += weight;
					
						weight = factor/6.0f*diffusionWeightFunction(map[xyzi],map[xyzi+nix*niy],scale);
						num += weight*map[xyzi+nix*niy];
						den += weight;
						
						map[xyzi] = num/den;
						
						maxdiff = Numerics.max(maxdiff, Numerics.abs(map[xyzi]-prev));
					} else {
						map[xyzi] = orig[xyzi];
					}
				}
				BasicInfo.displayMessage("max diff. "+maxdiff+"\n");
			}
			// store in the gain if large enough
			for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
				for (int m=0;m<ngain+1;m++) {
					if (newlabelHD[m][xyzi]==EMPTY) {
						newgainHD[m][xyzi] = 2.0f*map[xyzi]-1.0f;
						newlabelHD[m][xyzi] = n;
						m = ngain+1;
					} else if (2.0f*map[xyzi]-1.0f>newgainHD[m][xyzi]) {
						for (int p=ngain;p>m;p--) {
							newgainHD[p][xyzi] = newgainHD[p-1][xyzi];
							newlabelHD[p][xyzi] = newlabelHD[p-1][xyzi];
						}
						newgainHD[m][xyzi] = 2.0f*map[xyzi]-1.0f;
						newlabelHD[m][xyzi] = n;
						m=ngain+1;
					}
				}
			}
		}
		// make a hard copy
		for (int m=0;m<ngain+1;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			bestgainHD[m][xyzi] = newgainHD[m][xyzi];
			bestlabelHD[m][xyzi] = newlabelHD[m][xyzi];
		}
		newgainHD = null;
		newlabelHD = null;
    }
    /*
    // probability diffusion
    public final void diffuseBestGainFunctions2(int iter, float scale) {
    	
    	//mix with the neighbors?
    	float[][][] map = new float[nix][niy][niz];
    	float[][][] orig = new float[nix][niy][niz];
    	float[][] newgainHD = new float[ngain+1][nix*niy*niz];
    	byte[][] newlabelHD = new byte[ngain+1][nix*niy*niz];
    	boolean[][][] mask = new boolean[nix][niy][niz];
		for (int m=0;m<ngain+1;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			newlabelHD[m][xyzi] = EMPTY;
		}
    	
    	for (byte n=0;n<nobj;n++) {	
    		BasicInfo.displayMessage("propagate gain for label "+n+"\n");
			// get the gain ; normalize
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				int xyzi = x+nix*y+nix*niy*z;
				float sum = 0.0f;
				map[x][y][z] = 0.0f;
				for (int m=0;m<ngain+1;m++) {
					if (bestlabelHD[m][xyzi]==n) map[x][y][z] = 0.5f+0.5f*bestgainHD[m][xyzi];
					sum += 0.5f+0.5f*bestgainHD[m][xyzi];
				}
				
				// normalize over sum proba?
				//map[x][y][z] /= sum;
				
				orig[x][y][z] = map[x][y][z];
				if (orig[x][y][z]==0 || orig[x][y][z]==1) mask[x][y][z] = false;
				else mask[x][y][z] = true;
			}
			// propagate the values : diffusion
			float maxdiff = 1.0f;
			for (int t=0;t<iter && maxdiff>0.05f;t++) {
				BasicInfo.displayMessage(".");
				maxdiff = 0.0f;
				for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
					int xyzi = x+nix*y+nix*niy*z;
					if (mask[x][y][z]) {
						float den = diffusionWeightFunction(map[x][y][z],orig[x][y][z],scale);
						float num = den*orig[x][y][z];
						float weight;
						float prev = map[x][y][z];
						
						weight = diffusionWeightFunction(map[x][y][z],map[x-1][y][z],scale);
						num += weight*map[x-1][y][z];
						den += weight;
						
						weight = diffusionWeightFunction(map[x][y][z],map[x+1][y][z],scale);
						num += weight*map[x+1][y][z];
						den += weight;
						
						weight = diffusionWeightFunction(map[x][y][z],map[x][y-1][z],scale);
						num += weight*map[x][y-1][z];
						den += weight;
						
						weight = diffusionWeightFunction(map[x][y][z],map[x][y+1][z],scale);
						num += weight*map[x][y+1][z];
						den += weight;
						
						weight = diffusionWeightFunction(map[x][y][z],map[x][y][z-1],scale);
						num += weight*map[x][y][z-1];
						den += weight;
						
						weight = diffusionWeightFunction(map[x][y][z],map[x][y][z+1],scale);
						num += weight*map[x][y][z+1];
						den += weight;
		
						map[x][y][z] = num/den;
						
						maxdiff = Numerics.max(maxdiff, Numerics.abs(map[x][y][z]-prev));
					} else {
						map[x][y][z] = orig[x][y][z];
					}
				}
				BasicInfo.displayMessage("max diff. "+maxdiff+"\n");
			}
			// store in the gain if large enough
			for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
				int xyzi = x+nix*y+nix*niy*z;
				for (int m=0;m<ngain+1;m++) {
					if (newlabelHD[m][xyzi]==EMPTY) {
						newgainHD[m][xyzi] = 2.0f*map[x][y][z]-1.0f;
						newlabelHD[m][xyzi] = n;
						m = ngain+1;
					} else if (2.0f*map[x][y][z]-1.0f>newgainHD[m][xyzi]) {
						for (int p=ngain;p>m;p--) {
							newgainHD[p][xyzi] = newgainHD[p-1][xyzi];
							newlabelHD[p][xyzi] = newlabelHD[p-1][xyzi];
						}
						newgainHD[m][xyzi] = 2.0f*map[x][y][z]-1.0f;
						newlabelHD[m][xyzi] = n;
						m=ngain+1;
					}
				}
			}
		}
		// make a hard copy
		for (int m=0;m<ngain+1;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			bestgainHD[m][xyzi] = newgainHD[m][xyzi];
			bestlabelHD[m][xyzi] = newlabelHD[m][xyzi];
		}
		newgainHD = null;
		newlabelHD = null;
    }
    */
    
    private final float diffusionWeightFunction(float val, float ngb, float scale) {
    	
    	//return 1.0f/(1.0f+Numerics.square( (val-ngb)/scale ));
    	//return 1.0f/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	//return 1.0f/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    	return Numerics.square(Numerics.min(1.0f-val, val)/scale)/(1.0f+Numerics.square(Numerics.min(1.0f-val, val)/scale))
    				/(1.0f+Numerics.square( (val-ngb)/scale ))/(1.0f+Numerics.square(Numerics.min(1.0f-ngb, ngb)/scale));	
    }
    
    public final void propagateBestGainFunctions(int iter, float scale) {
    	/* use the gain itself instead?*/
    	byte iso=-1;
    	byte t1map = -1;
    	for (byte n=0;n<nc;n++) {
			if (modality[n].equals("Iso") || modality[n].equals("FlatImg")) iso=n;
			else if (modality[n].equals("T1map")) t1map=n;
		}
    	
		//mix with the neighbors?
    	float[] map = new float[nax*nay*naz];
    	float[] orig = new float[nax*nay*naz];
    	float[][] newgain = new float[ngain+1][nax*nay*naz];
    	byte[][] newlabel = new byte[ngain+1][nax*nay*naz];
    	float[][] weight = new float[6][nax*nay*naz];
		
    	boolean useBoth = true;
    	
    	// init: create memberships
		for (int xyza=0;xyza<nax*nay*naz;xyza++) {
			float sum = 0.0f;
			for (int m=0;m<ngain+1;m++) {
				sum += 1.0f+bestgain[m][xyza];
			}
			for (int m=0;m<ngain+1;m++) {
				bestgain[m][xyza] = (1.0f+bestgain[m][xyza])/sum;
			}
			weight[pX][xyza] = 2.0f;
			weight[mX][xyza] = 2.0f;
			weight[pY][xyza] = 2.0f;
			weight[mY][xyza] = 2.0f;
			weight[pZ][xyza] = 2.0f;
			weight[mZ][xyza] = 2.0f;
		}

		// init: create a weight map
		float[] XA = new float[3];
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x + nix*y + nix*niy*z;
			atlas.imageToShapeCoordinates(XA, x, y, z);
			int xa = Numerics.round(XA[X]);
			int ya = Numerics.round(XA[Y]);
			int za = Numerics.round(XA[Z]);
			if (xa>=0 && xa<nax && ya>=0 && ya<nay && za>=0 && za<naz) {
				int xyza = xa + nax*ya + nax*nay*za;
				float wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi+1])/(scale*imrange[t1map]) ));
				if (wij<weight[pX][xyza]) weight[pX][xyza] = wij;
				
				wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi-1])/(scale*imrange[t1map]) ));
				if (wij<weight[mX][xyza]) weight[mX][xyza] = wij;
				
				wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi+nix])/(scale*imrange[t1map]) ));
				if (wij<weight[pY][xyza]) weight[pY][xyza] = wij;
				
				wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi-nix])/(scale*imrange[t1map]) ));
				if (wij<weight[mY][xyza]) weight[mY][xyza] = wij;
				
				wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi+nix*niy])/(scale*imrange[t1map]) ));
				if (wij<weight[pZ][xyza]) weight[pZ][xyza] = wij;
				
				wij = 1.0f/(1.0f + Numerics.square( (image[t1map][xyzi]-image[t1map][xyzi-nix*niy])/(scale*imrange[t1map]) ));
				if (wij<weight[mZ][xyza]) weight[mZ][xyza] = wij;
			}
		}				
		for (int xyza=0;xyza<nax*nay*naz;xyza++) for (int d=0;d<6;d++) {
			if (weight[d][xyza]==2) weight[d][xyza]=0;
		}
    	for (byte n=0;n<nobj;n++) {				
     		BasicInfo.displayMessage("propagate gain for label"+n+"\n");
			// get the gain ; normalize
			for (int xyza=0;xyza<nax*nay*naz;xyza++) {
				for (int m=0;m<ngain+1;m++) if (bestlabel[m][xyza]==n) {
					map[xyza] = bestgain[m][xyza];
					orig[xyza] = map[xyza];
				}
			}
			// propagate the values : diffusion
			for (int t=0;t<iter;t++) {
				BasicInfo.displayMessage(".");
				for (int x=1;x<nax-1;x++) for (int y=1;y<nay-1;y++) for (int z=1;z<naz-1;z++) {
					int xyza = x+nax*y+nax*nay*z;
					float num = orig[xyza];
					float den = 1.0f;
					
					if (t1map>-1) {
						if (useBoth || map[xyza-1]>map[xyza]) {
							num += weight[mX][xyza]*(map[xyza-1]-map[xyza]);
							den += weight[mX][xyza];
						}
						if (useBoth || map[xyza+1]>map[xyza]) {
							num += weight[pX][xyza]*(map[xyza+1]-map[xyza]);
							den += weight[pX][xyza];
						}
						if (useBoth || map[xyza-nax]>map[xyza]) {
							num += weight[mY][xyza]*(map[xyza-nax]-map[xyza]);
							den += weight[mY][xyza];
						}
						if (useBoth || map[xyza+nax]>map[xyza]) {
							num += weight[pY][xyza]*(map[xyza+nax]-map[xyza]);
							den += weight[pY][xyza];
						}
						if (useBoth || map[xyza-nax*nay]>map[xyza]) {
							num += weight[mZ][xyza]*(map[xyza-nax*nay]-map[xyza]);
							den += weight[mZ][xyza];
						}
						if (useBoth || map[xyza+nax*nay]>map[xyza]) {
							num += weight[pZ][xyza]*(map[xyza+nax*nay]-map[xyza]);
							den += weight[pZ][xyza];
						}
					}
					map[xyza] = num/den;
				}
			}
			// store in the gain if large enough
			for (int xyza=0;xyza<nax*nay*naz;xyza++) {
				for (int m=0;m<ngain+1;m++) {
					if (map[xyza]>newgain[m][xyza]) {
						for (int p=ngain;p>m;p--) {
							newgain[p][xyza] = newgain[p-1][xyza];
							newlabel[p][xyza] = newlabel[p-1][xyza];
						}
						newgain[m][xyza] = map[xyza];
						newlabel[m][xyza] = n;
						m=ngain+1;
					}
				}
			}
		}
		bestgain = newgain;
		bestlabel = newlabel;
    }
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// option flags (for testing)
    	boolean fullMarching=false;
    	// not needed at all, it seems (removed from code)
    	//boolean recomputeCritical=false;
    	
		// init decomposition
		fastMarchingInitializationFromSegmentation(false, fullMarching);
			
    	// first estimate the narrow band size
    	int size = 0;
		int boundarysize=0;
		
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
			if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) size++;
			if (Numerics.abs(mgdmfunctions[0][xyz])<1.0 && mgdmfunctions[0][xyz]!=UNKNOWN) boundarysize++;
		}
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
     	if (debug) System.out.print("init ("+size+")\n");
        
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
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
				//if (true) {
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
							
							// try all possible labels: no, next label is the closest by definition
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
								// reset to low value ?
								narrowband.functions[lb][n] = lowlevel;
								// or to original value ?
								//narrowband.functions[lb][n] -= Numerics.bounded(forces[lb] - forces[lb+1], -0.9f, 0.9f);
							}
							// reset the counter; also reset all the neighbors
							if (lb==0) {
								counter[xyz]=0;
								for (int b=0;b<NGB;b++) {
									int xyzn = xyz + ngbx[b] + ngby[b]*nax + ngbz[b]*nax*nay;
									// reset all to zero?
									counter[xyzn] = 0;
									/*
									if (counter[xyzn]>maxcount) counter[xyzn] = (short)(maxcount-1);
									else if (counter[xyzn]>0) counter[xyzn]--;
									*/
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
				for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
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
				int xyzn = xyz + ngbx[b] + ngby[b]*nax + ngbz[b]*nax*nay;

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
			/*
			Dmx[n] = phi[CTR] - phi[mX];
			Dmy[n] = phi[CTR] - phi[mY];
			Dmz[n] = phi[CTR] - phi[mZ];
			
			Dpx[n] = phi[pX] - phi[CTR];
			Dpy[n] = phi[pY] - phi[CTR];
			Dpz[n] = phi[pZ] - phi[CTR];
			
			D0x[n] = (phi[pX] - phi[mX])/2.0;
			D0y[n] = (phi[pY] - phi[mY])/2.0;
			D0z[n] = (phi[pZ] - phi[mZ])/2.0;
			*/
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
			
			forces[n] = -smoothweight*stepsize*(K-K0)*GPhi;
			
			// divergence-smoothing force:
			forces[n] = -divweight*stepsize*(phi[CTR]*phi[CTR]/(1.0+phi[CTR]*phi[CTR]))*(Dxx+Dyy+Dzz);
			
			// distance-based atenuation
			/* one-sided? */
			if (phi[CTR]<0) distval[n] = 1.0;
			else distval[n] = 1.0/(1.0+phi[CTR]*phi[CTR]/gaindist2);
			/* or two sided? not much difference it seems
			distval[n] = 1.0/(1.0+phi[CTR]*phi[CTR]/gaindist2);
			*/
		}
		
		// find the best target in the neighborhood to build the corresponding force
		/*
		bestlb = EMPTY;
		gainval = 0.5f;
		for (int n=0;n<=ngain && bestlb==EMPTY;n++) {
			gainlb = bestlabel[n][xyz];
			gainval = bestgain[n][xyz];
			for (byte l=0;l<=nmgdm && bestlb==EMPTY;l++) {
				if (mgdmlabels[l][xyz]==gainlb) {
					bestlb = l;
				}
			}
		}
		*/
		// (product of gain and proximity => must check all)
		bestlb = EMPTY;
		bestval = -1.0f;
		done=false;
		for (int n=0;n<=ngain && !done;n++) {
			gainlb = bestlabel[n][xyz];
			gainval = 0.5f + 0.5f*bestgain[n][xyz];
			//if (gainval<=0) done = true;
			for (byte l=0;l<=nmgdm && !done;l++) {
				if (mgdmlabels[l][xyz]==gainlb  && distval[l]*gainval>bestval) {
					bestlb = l;
					bestval = distval[l]*gainval;
				}
			}
		}
		
		/*
		// use the difference with next best??
		if (bestlb<ngain) gainval -= bestgain[bestlb+1][xyz];
		else gainval += 1.0;
		*/
		/*
		// only the first label?
		gainlb = bestlabel[0][xyz];
		for (byte l=0;l<=nmgdm && bestlb==EMPTY;l++) {
			if (mgdmlabels[l][xyz]==gainlb) bestlb = l;
		}
		*/
		
		// second pass for the data force
		if (bestlb!=EMPTY) {
		//if (bestlb!=EMPTY && gainval>0) {	// use only positive forces??
			for (int n=0;n<=nmgdm;n++) {
				byte lb = mgdmlabels[n][xyz];
				if (lb!=EMPTY) {
					/*
					// central differences?
					forces[n] += forceweight/smoothfactor[lb]*stepsize*bestval*(D0x[bestlb]*D0x[n] + D0y[bestlb]*D0y[n] + D0z[bestlb]*D0z[n]);
					*/
					// upwind scheme?
					forces[n] += forceweight/smoothfactor[lb]*stepsize*bestval
								*(Numerics.max(Dpx[bestlb]+Dmx[bestlb],0)*Dmx[n] + Numerics.min(Dpx[bestlb]+Dmx[bestlb],0)*Dpx[n] 
								 +Numerics.max(Dpy[bestlb]+Dmy[bestlb],0)*Dmy[n] + Numerics.min(Dpy[bestlb]+Dmy[bestlb],0)*Dpy[n] 
								 +Numerics.max(Dpz[bestlb]+Dmz[bestlb],0)*Dmz[n] + Numerics.min(Dpz[bestlb]+Dmz[bestlb],0)*Dpz[n]);
					
					// add a balloon term if prod ~ 0: not useful
				} else {
					// regularize more ? do nothing ?
					//forces[n] += -forceweight*stepsize*0.5f;
				}
			}
		} else {
			// if best is not a neighbor, just shrink the label
			// or should we just set things to zero then (and let the smoothness drive everything)? 
			// ->not good, creates static points
			
			// increase smoothness in that location instead? do nothing?
			/*
			for (int n=0;n<=nmgdm;n++) {
				forces[n] += -forceweight*stepsize*0.5f;
			}
			*/
		}
		
		/*
		//just balloon forces??
		for (int n=0;n<=nmgdm;n++) {
			byte lb = mgdmlabels[n][xyz];
			if (lb!=EMPTY) {
				boolean found=false;
				for (byte l=0;l<=nmgdm && !found;l++) if (lb==bestlabel[l][xyz]) {
					found=true;
					forces[n] += forceweight/smoothfactor[lb]*stepsize*bestgain[l][xyz];
				}
				if (!found) forces[n] += -1.0*forceweight/smoothfactor[lb]*stepsize;
			}
		}
		*/
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
			if (segmentation[xyz+i+j*nax+l*nax*nay]==lb) {
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
			if (segmentation[xyz+i+j*nax+l*nax*nay]==segmentation[xyz]) {
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
				&& (segmentation[xyz+i+j*nax+l*nax*nay]!=lb) 
				&& (segmentation[xyz+i+j*nax+l*nax*nay]!=segmentation[xyz]) ) {
				found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[n]) { found = true; break; }
				
				if (!found) {
					lbs[Nconfiguration] = segmentation[xyz+i+j*nax+l*nax*nay];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xyz+i+j*nax+l*nax*nay]==segmentation[xyz])
					|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[n]) ) {
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
				if ( (segmentation[xyz+i+j*nax+l*nax*nay]==lb)
					|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[n]) ) {
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
					if ( (segmentation[xyz+i+j*nax+l*nax*nay]==segmentation[xyz])
						|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[n])
						|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[m]) ) {
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
					if ( (segmentation[xyz+i+j*nax+l*nax*nay]==lb)
						|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[n]) 
						|| (segmentation[xyz+i+j*nax+l*nax*nay]==lbs[m]) ) {
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
	
	public final void fastMarchingInitializationFromSegmentation(boolean narrowBandOnly, boolean almostEverywhere) {
         // initialize the quantities
         for (int xyz = 0; xyz<nax*nay*naz; xyz++) {
         	 // mgdm functions
			for (int n = 0; n<nmgdm; n++) {
            	mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
            mgdmlabels[nmgdm][xyz] = EMPTY;
        }
		
        // basic mask region: boundaries
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			int xyz = x+nax*y+nax*nay*z;
			if (x==0 || x==nax-1 || y==0 || y==nay-1 || z==0 || z==naz-1) mask[xyz] = false;
			else mask[xyz] = true;
		}
		// computation variables
        byte[] processed = new byte[nax*nay*naz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;
		boolean done, isprocessed;
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
        	processed[xyz] = 0;
        	// search for boundaries
        	for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				if (segmentation[xyzn]!=segmentation[xyz]) if (mask[xyzn]) {
					// add to the heap
					heap.addValue(0.5f,xyzn,segmentation[xyz]);
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
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			int xyz = x+nax*y+nax*nay*z;

			if (x==0 || x==nax-1 || y==0 || y==nay-1 || z==0 || z==naz-1) mask[xyz] = false;
			else if (segmentation[xyz]==EMPTY || (processed[xyz]==0 && segmentation[xyz]==0)) mask[xyz] = false;
			else mask[xyz] = true;
		}
		
		// to create the MGDM functions, we need to copy the segmentation, forget the last labels
		// and compute differences between distance functions
		if (debug) BasicInfo.displayMessage("transform into MGDM functions\n");		
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
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
        byte[] processed = new byte[nax*nay*naz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;	
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		

		long start_time = System.currentTimeMillis(); 

        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
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
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
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
    	float[] tmp = new float[nax*nay*naz];
    	boolean[] processed = new boolean[nax*nay*naz];
    	
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) if (mask[xyz]) {
			
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
		for (int xyz = 0; xyz<nax*nay*naz; xyz++) {
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

