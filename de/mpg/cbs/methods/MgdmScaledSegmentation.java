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
 *	@version    May 2011
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class MgdmScaledSegmentation {
	
	// object types

	private	static	final	byte	EMPTY = -1;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    private final static byte W = 3;
    	
	private final static byte NORTH 	= 0;
    private final static byte SOUTH 	= 1;
	private final static byte EAST 		= 2;
    private final static byte WEST 		= 3;
	private final static byte FRONT 	= 4;
    private final static byte BACK 		= 5;
    
	public final static byte ROUNDED 	= 1;
    public final static byte POINTED 	= 2;
	public final static byte SMOOTH	= 3;
	public final static byte NONE		= 4;
	
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
	
	
    // numerical quantities
	private static final	float   INF=1e15f;
	private static final	float   ZERO=1e-15f;
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private final static float SQR2 = (float) Math.sqrt(2.0f);
    private final static float SQR3 = (float) Math.sqrt(3.0f);
    private final static float diagdist = 1/(2*SQR2);
    private final static float cubedist = 1/(2*SQR3);
	private static final	float	UNKNOWN = -1.0f;
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;


	// data and membership buffers
	private 	float[][] 		mgdmfunctions;  	// MGDM's pseudo level set mgdmfunctions
	private 	byte[][] 		mgdmlabels;   		// MGDM's label maps
	private 	byte[] 			segmentation;   	// MGDM's segmentation
	private 	float[][] 		fieldforce;  		// original image forces, indep. object (e.g. boundaries)
	private 	byte[] 			fieldtarget;  		// original image forces, indep. object (e.g. boundaries)
	private 	float[][] 		balloonforces;  	// original image forces, along object normals (e.g. from memberships)
	private 	float[] 		curvatureforce;  	// weighting factor for the curvature force
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nix,niy,niz;   		// images dimensions
	private static	float 		rix,riy,riz;   		// images resolutions
	private static	int 		nbx,nby,nbz;   		// balloon forces dimensions
	private	static	float 		rbx,rby,rbz;   		// balloon forces resolutions
	private	static	boolean		scaledBalloon;		// whether to use specific dim,res for balloon forces
	private static	int 		ncx,ncy,ncz;   		// balloon forces dimensions
	private	static	float 		rcx,rcy,rcz;   		// balloon forces resolutions
	private	static	boolean		scaledCurv;		// whether to use specific dim,res for balloon forces
	private static	int 		nfx,nfy,nfz;   		// field forces dimensions
	private	static	float 		rfx,rfy,rfz;   		// field forces resolutions
	private static 	boolean		scaledField;		// whether to use specific dim,res for field forces
	//private		byte[]		curvaturedir;
	
	
	private static	int 	  	nobj;					// total number of objects to represent (including background)
	private static	int 	  	nmgdm;					// total number of MGDM mgdmlabels and mgdmfunctions
	private 	int[]			objLabel;			// label values in the original image
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// parameters
	private	double		smoothweight, fieldweight, balloonweight;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 5.0f;
	//private float		landmineDist = 2.0f;
	private	float		narrowBandDist = landmineDist+1.8f;
	private	float		extraDist = narrowBandDist+1.0f;
	
	// computation variables to avoid re-allocating
	float[][][] phi = new float[3][3][3];
	float[][] phii;
	int xyz, xyzn, xyznb;
	int x,y,z;
	byte lb,olb;
	double Dmx, Dmy, Dmz, Dpx, Dpy, Dpz, D0x, D0y, D0z;
	double SD0x, SD0y, SD0z, GPhi;
	double Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
	double K, G, tmp;
	double DeltaP, DeltaM;
	double smooth, balloon, field;
	double pos, neg, curr, next, val, vc, vn;
	double Dmx0,Dmy0,Dmz0,Dpx0,Dpy0,Dpz0;
	int npos, nneg;
	float dist, maxdist, newdist, avgdiff;
	boolean changed;
	double boundaryfactor;
	double	sigmaB = 0.25*narrowBandDist*narrowBandDist;
	double s, s2; 
    int count;
        
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
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
		
		public final void addPoint(int xyz, byte[][] mgdmlb, float[][] mgdmfn) {
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
				labels[l][currentsize] = mgdmlb[l][xyz];
				functions[l][currentsize] = mgdmfn[l][xyz];
			}
			labels[nmgdm][currentsize] = mgdmlb[nmgdm][xyz];
			currentsize++;
		}
		
		public void reset() {
			currentsize = 0;
		}
	}


	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public MgdmScaledSegmentation(int[] init_, int nix_, int niy_, int niz_, float rix_, float riy_, float riz_,
								int nobj_, int nmgdm_,
								float[][] field_, byte[] trg_, int nfx_, int nfy_, int nfz_, float rfx_, float rfy_, float rfz_,
								float[][] balloon_, int nbx_, int nby_, int nbz_, float rbx_, float rby_, float rbz_,
								float[] curv_, int ncx_, int ncy_, int ncz_, float rcx_, float rcy_, float rcz_,
								float bw_, float fw_, float sw_, 
								String connectivityType_) {
	
		nix = nix_;
		niy = niy_;
		niz = niz_;
		
		nobj = nobj_;
		nmgdm = nmgdm_;
		
		rix = rix_;
		riy = riy_;
		riz = riz_;

		fieldforce = field_;
		fieldtarget = trg_;
		nfx = nfx_;
		nfy = nfy_;
		nfz = nfz_;
		rfx = rfx_;
		rfy = rfy_;
		rfz = rfz_;
		if (nfx==nix && nfy==niy && nfz==niz) scaledField = false;
		else scaledField = true;

		balloonforces = balloon_;
		nbx = nbx_;
		nby = nby_;
		nbz = nbz_;
		rbx = rbx_;
		rby = rby_;
		rbz = rbz_;
		if (nbx==nix && nby==niy && nbz==niz) scaledBalloon = false;
		else scaledBalloon = true;
		
		curvatureforce = curv_;
		ncx = ncx_;
		ncy = ncy_;
		ncz = ncz_;
		rcx = rcx_;
		rcy = rcy_;
		rcz = rcz_;
		if (ncx==nix && ncy==niy && ncz==niz) scaledCurv = false;
		else scaledCurv = true;
				
		balloonweight = bw_;
		fieldweight = fw_;
		smoothweight = sw_;
		
		if (debug) System.out.print("MGDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
		
		objLabel = ObjectLabeling.listOrderedLabels(init_, nix, niy, niz);
		// note: we do expect that there are nb objects (not checked)
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nix, -nix, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nix*niy, -nix*niy};
		
		// init all the arrays
		try {
			mgdmfunctions = new float[nmgdm][nix*niy*niz];
			mgdmlabels = new byte[nmgdm+1][nix*niy*niz];	
			segmentation = new byte[nix*niy*niz];	
			mask = new boolean[nix*niy*niz];
			phii = new float[nmgdm+1][19];
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
					//System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nix; x++) for (int y=0; y<niy; y++) for (int z = 0; z<niz; z++) {
			if (x>1 && x<nix-2 && y>1 && y<niy-2 && z>1 && z<niz-2) mask[x+nix*y+nix*niy*z] = true;
			else mask[x+nix*y+nix*niy*z] = false;
		}
		// init decomposition
		fastMarchingInitializationFromSegmentation(init_);
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		mgdmfunctions = null;
		mgdmlabels = null;
		segmentation = null;
		fieldforce = null;
		balloonforces = null;
		heap.finalize();
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
    
	public final void setWeights(float fw_, float bw_, float sw_) {
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		
		if (debug) System.out.print("MGDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing)\n");
	}
	
	public final void reduceMGDMsize(int nred) {
		
		float[][] redfunctions = new float[nred][nix*niy*niz];
		byte[][] redlabels = new byte[nred+1][nix*niy*niz];	
				
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) {
			for (int n = 0; n<nred; n++) {
         		redfunctions[n][xyz] = mgdmfunctions[n][xyz];	
         		redlabels[n][xyz] = mgdmlabels[n][xyz];	
         	}
         	redlabels[nred][xyz] = mgdmlabels[nred][xyz];
        }
        
        // replace the maps
        nmgdm = nred;
        mgdmfunctions = redfunctions;
        mgdmlabels = redlabels;
	}
	
	public final void fastMarchingInitializationFromSegmentation(int[] init) {
         // initialize the quantities
         for (int xyz = 0; xyz<nix*niy*niz; xyz++) {
         	 // mgdm functions
			for (int n = 0; n<nmgdm; n++) {
            	mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
            mgdmlabels[nmgdm][xyz] = EMPTY;
            
            // segmentation
            byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==init[xyz]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
        }
		
        // computation variables
        byte[] processed = new byte[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
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
					heap.addValue(0.5f,xyzn,segmentation[xyz]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			boolean done = false;
			for (int n=0; n<processed[xyz]; n++)
				if (mgdmlabels[n][xyz]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xyz]][xyz] = dist;
			mgdmlabels[processed[xyz]][xyz] = lb;
			processed[xyz]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					boolean isprocessed = false;
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
						
						if (newdist<=narrowBandDist+extraDist) {
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
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }
     
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingReinitialization(boolean narrowBandOnly) {
        // computation variables
        byte[] processed = new byte[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
			
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
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
				xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
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
        	dist = heap.getFirst();
        	xyz = heap.getFirstId();
        	lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			done = false;
			for (int n=0; n<processed[xyz]; n++)
				if (mgdmlabels[n][xyz]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xyz]][xyz] = dist;
			mgdmlabels[processed[xyz]][xyz] = lb;
			processed[xyz]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
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
							xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							for (int n=0; n<processed[xyznb]; n++) if (mask[xyznb]) if (mgdmlabels[n][xyznb]==lb) {
								nbdist[l] = mgdmfunctions[n][xyznb];
								nbflag[l] = true;
							}			
						}
						newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if (!narrowBandOnly || newdist<=narrowBandDist+extraDist) {
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
		if (debug) BasicInfo.displayMessage("done\n");		
		
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
    	double[][] bforce = new double[nmgdm+1][4];
    	double[] bval = new double[nmgdm+1];
    	double[] best = new double[nmgdm+1];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
    		basicBalloonForces(xyz, bval, best);
    		boolean found=false;
    		int lb = -1;
    		for (int l=0;l<nmgdm;l++) if (mgdmlabels[l][xyz]==n && mgdmfunctions[l][xyz]!=UNKNOWN) {
    			found = true;
    			lb = l;
    		}
    		if (found) {
				forces[xyz] = (float)levelsetForces(xyz, n, bval[lb]);
			} else {
				forces[xyz] = 0;
			}
    	}
    	return forces;
    }
    public final int[] reconstructedLevelsetType() {
    	int[] type = new int[nix*niy*niz];
    	double[][] bforce = new double[nmgdm+1][4];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			centralBalloonForces(xyz, bforce);
			if (bforce[0][W]>0 && bforce[1][W]<0) {
				type[xyz] = 1;
			} else if (bforce[0][W]<0 && bforce[1][W]>0) {
				type[xyz] = 2;
			} else if (bforce[0][W]>0 && bforce[1][W]>0) {
				type[xyz] = 3;
			} else if (bforce[0][W]<0 && bforce[1][W]<0) {
				type[xyz] = 4;
			} else {
				type[xyz] = 0;
			}
    	}
    	return type;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[] reconstructedFastMarchingForces(byte n) {
    	float[] forces = new float[nix*niy*niz];
    	for (int xyz=0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)fastMarchingForces(xyz,n);
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
    		if (mgdmlabels[0][xyz]>-1) {
				seg[xyz] = objLabel[mgdmlabels[0][xyz]];
			}
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
    		if (mgdmlabels[0][xyz]>-1) {
				seg[xyz] = objLabel[mgdmlabels[0][xyz]];
			}
    	}
    	return seg;
    }

    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBand(int iter) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// init
    	
    	// first estimate the narrow band size
    	int size = 0;
		int boundarysize=0;
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) size++;
			if (Numerics.abs(mgdmfunctions[0][xyz])<1.0) boundarysize++;
		}
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
		for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		double[][] bforce = new double[nmgdm+1][4];
		int lbmax, lbsec;
		double curr, next;
		boolean reinitLM, reinitOL;

		// evolve until a landmine is closer than minDist of the boundaries
		for (int t=0;t<iter;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			reinitLM = false;
			reinitOL = false;
			
			for (int lb=0;lb<nmgdm;lb++) nswap[lb] = 0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// select best label from previous initialization
				olb = mgdmlabels[nmgdm][xyz];
				
				// evolve the MGDM functions
				
				// compute the forces from current levelset values, update the narrow band from it
				if (olb!=EMPTY) forces[nmgdm] = levelsetForces(xyz, olb);
				else forces[nmgdm] = 0.0f;
				lbmax = nmgdm;
				lbsec = -1;
				for (int lb=nmgdm-1;lb>=0;lb--) {
					if (mgdmlabels[lb][xyz]!=EMPTY) {
						forces[lb] = levelsetForces(xyz, mgdmlabels[lb][xyz]);
						if (forces[lb]>forces[lbmax]) {
							lbsec = lbmax;
							lbmax = lb;
						} else if (lbsec==-1) {
							lbsec = lb;
						} else if (forces[lb]>forces[lbsec]) {
							lbsec = lb;
						}
					} else forces[lb] = 0.0f;					
				}
				
				for (int lb=nmgdm-1;lb>=0;lb--) if (mgdmlabels[lb][xyz]!=EMPTY) {
					
					// update the narrow band values, not the original data
					narrowband.functions[lb][n] += Numerics.bounded(forces[lb] - forces[lb+1], -0.9f, 0.9f);

					// change of sign ?
					if (narrowband.functions[lb][n]<0) {
						//if (debug) System.out.print(""+lb);
						nswap[lb]++;
						
						if (lb==nmgdm-1) {
							if (olb!=EMPTY) {
								//if (homeomorphicLabeling(xyz, olb)) {
								narrowband.labels[lb][n] = olb;
								narrowband.functions[lb][n] = -narrowband.functions[lb][n];
							} else {
								// reset to low value
								narrowband.functions[lb][n] = lowlevel;
							}
						} else if (mgdmlabels[lb+1][xyz]!=EMPTY) {
							if (lb==0) {
								// check for topology here (optional)
								if (homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz])) {
									narrowband.labels[lb+1][n] = mgdmlabels[lb][xyz];
									narrowband.labels[lb][n] = mgdmlabels[lb+1][xyz];
									narrowband.functions[lb][n] = -narrowband.functions[lb][n];
									segmentation[xyz] = mgdmlabels[lb+1][xyz];
									// check for boundariy changes in the landmines : force reinitialization
									if (landmines.get(xyz)) reinitLM = true;
									// check for far labels getting mixed in: time to re-initialize
									if (narrowband.labels[0][n]==mgdmlabels[nmgdm][xyz]) reinitOL = true;
								} else {
									/* optional?
									// first, check for other possible labels (nearby)
									dist = mgdmfunctions[lb][xyz];
									changed=false;
									maxdist = 2.0f;
									for (int k=1;lb+k<nmgdm && dist<maxdist && !changed;k++) {
										dist += mgdmfunctions[lb+k][xyz];
										if (dist<maxdist) {
											if (lb+k+1<nmgdm && homeomorphicLabeling(xyz, mgdmlabels[lb+k+1][xyz])) {
												//System.out.print(".");
												narrowband.labels[lb+k+1][n] = mgdmlabels[lb][xyz];
												narrowband.labels[lb][n] = mgdmlabels[lb+k+1][xyz];
												narrowband.functions[lb][n] = -narrowband.functions[lb][n];
												segmentation[xyz] = mgdmlabels[lb+k+1][xyz];
												// check for boundary changes in the landmines : force reinitialization
												if (landmines.get(xyz)) reinitLM = true;
												// check for far labels getting mixed in: time to re-initialize
												if (narrowband.labels[0][n]==otherlabels[xyz]) reinitOL = true;
												changed=true;
											}
										}
									}
									if (!changed) {
										// reset to low value
										narrowband.functions[lb][n] = lowlevel;
									}
									*/
									
									// reset to low value
									narrowband.functions[lb][n] = lowlevel;
								}
							} else {
								narrowband.labels[lb+1][n] = mgdmlabels[lb][xyz];
								narrowband.labels[lb][n] = mgdmlabels[lb+1][xyz];
								narrowband.functions[lb][n] = -narrowband.functions[lb][n];
							}
						} else {
							// reset to low value
							narrowband.functions[lb][n] = lowlevel;
						}
					}
				}
			}
			if (debug) for (int lb=0;lb<nmgdm;lb++) System.out.print("changed labels ("+lb+"): "+nswap[lb]+" ("+(nswap[lb]/(float)boundarysize*100.0f)+" % of boundary)\n");
			
			// once all the new values are computed, copy into original MGDM functions
			avgdiff = 0.0f;
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				// measure the changes at the base level
				if (mgdmlabels[0][xyz]==narrowband.labels[0][n]) avgdiff += Numerics.abs(mgdmfunctions[0][xyz]-narrowband.functions[0][n]);
				else avgdiff += Numerics.abs(mgdmfunctions[0][xyz]+narrowband.functions[0][n]);
				
				for (int lb=0;lb<nmgdm;lb++) {
					mgdmlabels[lb][xyz] = narrowband.labels[lb][n];
					mgdmfunctions[lb][xyz] = narrowband.functions[lb][n];
				}
			}
			if (debug) System.out.print("mean distance function change: "+(avgdiff/narrowband.currentsize)+"\n");
	
			//if (t<iter-1 && (reinitLM || reinitOL) ) {
			if (t<iter-1 && reinitLM) {
				if (debug) System.out.print("re-initialization (LM: "+reinitLM+" | OL: "+reinitOL+" )\n");
        		
				//resetIsosurfaceNarrowBand(narrowband);
				resetIsosurfaceBoundary();
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		fastMarchingReinitialization(true);
		
        return;
    }
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveStrictNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// init
    	
    	// first estimate the narrow band size
    	int size = 0;
		int boundarysize=0;
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
			if (mgdmfunctions[0][xyz]<narrowBandDist && mgdmfunctions[0][xyz]!=UNKNOWN) size++;
			if (Numerics.abs(mgdmfunctions[0][xyz])<1.0) boundarysize++;
		}
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
		for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		double[] bestforces = new double[nmgdm+1];
		byte[] bestid = new byte[nmgdm+1];
		double[][] bforce = new double[nmgdm+1][4];
		int lbmax, lbsec;
		int curlb, newlb;
		double curr, next, sum;
		boolean reinitLM, reinitOL;

		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && (t<5 || diff>mindiff);t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			reinitLM = false;
			reinitOL = false;
			
			for (int lb=0;lb<nmgdm;lb++) nswap[lb] = 0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// evolve the MGDM functions
				
				// compute the forces from current levelset values, update the narrow band from it
				/*
				for (int lb=nmgdm;lb>=0;lb--) {
					if (mgdmlabels[lb][xyz]!=EMPTY) forces[lb] = levelsetForces(xyz, mgdmlabels[lb][xyz]);
					else forces[lb] = 0.0f;
				}
				*/
				if (balloonweight>0) basicBalloonForces(xyz, forces, bestforces);
				
				for (int lb=nmgdm-1;lb>=0;lb--) if (mgdmlabels[lb][xyz]!=EMPTY) {
					
					// update the narrow band values, not the original data
					narrowband.functions[lb][n] += Numerics.bounded(levelsetForces(xyz, mgdmlabels[lb][xyz], forces[lb])
																	- levelsetForces(xyz, mgdmlabels[lb+1][xyz], bestforces[lb]),
																	 -0.9f, 0.9f);
					
					// change of sign ?
					if (narrowband.functions[lb][n]<0) {
						//if (debug) System.out.print(""+lb);
						
						/*
						// try all possible labels
						newlb = EMPTY;
						curlb = lb+1;
						while (newlb==EMPTY && curlb<=nmgdm) {
							if (mgdmlabels[curlb][xyz]!=EMPTY && homeomorphicLabeling(xyz, mgdmlabels[curlb][xyz]) )
								newlb = curlb;
							else
								curlb++;
						}
						/*
						newlb = EMPTY;
						if (mgdmlabels[lb+1][xyz]!=EMPTY && homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz]) ) newlb = lb+1;
						*/
						// allow any change of topology at second and third order labels
						if (mgdmlabels[lb+1][xyz]!=EMPTY && (lb>0 || homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz]))) {
							newlb = lb+1;
						} else {
							newlb = EMPTY;
						}
						// restrict changes of labeling at all levels
						//if (mgdmlabels[lb+1][xyz]!=EMPTY && homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz]) ) {
						if (newlb!=EMPTY) {
							// for all levels
							nswap[lb]++;
							narrowband.labels[lb][n] = mgdmlabels[newlb][xyz];
							narrowband.functions[lb][n] = -narrowband.functions[lb][n];
							// never switch the last label
							if (newlb<nmgdm) narrowband.labels[newlb][n] = mgdmlabels[lb][xyz];
							// update the segmentation with first label
							if (lb==0) segmentation[xyz] = mgdmlabels[newlb][xyz];
							// check for boundary changes in the landmines : force reinitialization
							if (lb==0 && landmines.get(xyz)) reinitLM = true;
							// check for far labels getting mixed in: time to re-initialize
							if (lb==0 && narrowband.labels[0][n]==mgdmlabels[nmgdm][xyz]) reinitOL = true;
						} else {
							// reset to low value
							narrowband.functions[lb][n] = lowlevel;
						}
					}
				}
			}
			diff = (nswap[0]/(float)boundarysize);
			if (debug) for (int lb=0;lb<nmgdm;lb++) System.out.print("changed labels ("+lb+"): "+nswap[lb]+" ("+(nswap[lb]/(float)boundarysize*100.0f)+" % of boundary)\n");
			
			// once all the new values are computed, copy into original MGDM functions
			avgdiff = 0.0f;
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				// measure the changes at the base level
				if (mgdmlabels[0][xyz]==narrowband.labels[0][n]) avgdiff += Numerics.abs(mgdmfunctions[0][xyz]-narrowband.functions[0][n]);
				else avgdiff += Numerics.abs(mgdmfunctions[0][xyz]+narrowband.functions[0][n]);
				
				for (int lb=0;lb<nmgdm;lb++) {
					mgdmlabels[lb][xyz] = narrowband.labels[lb][n];
					mgdmfunctions[lb][xyz] = narrowband.functions[lb][n];
				}
			}
			if (debug) System.out.print("mean distance function change: "+(avgdiff/narrowband.currentsize)+"\n");
	
			//if (t<iter-1 && (reinitLM || reinitOL) ) {
			if (t<iter-1 && reinitLM) {
				if (debug) System.out.print("re-initialization (LM: "+reinitLM+" | OL: "+reinitOL+" )\n");
        		
				//resetIsosurfaceNarrowBand(narrowband);
				resetIsosurfaceBoundary();
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		fastMarchingReinitialization(true);
		
        return;
    }
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    /*
    public final void evolveFastNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// init
    	
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
        
		for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		double[] bestforces = new double[nmgdm+1];
		byte[] bestid = new byte[nmgdm+1];
		double[][] bforce = new double[nmgdm+1][4];
		int lbmax, lbsec;
		int curlb, newlb;
		double curr, next, sum;
		boolean reinitLM, reinitOL;

		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && (t<5 || diff>mindiff);t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			reinitLM = false;
			reinitOL = false;
			
			for (int lb=0;lb<nmgdm;lb++) nswap[lb] = 0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// evolve the MGDM functions
				
				// compute the forces from current levelset values, update the narrow band from it
				allLevelsetForces(xyz, forces, bestid);
				
				for (int lb=nmgdm-1;lb>=0;lb--) if (mgdmlabels[lb][xyz]!=EMPTY) {
					
					// update the narrow band values, not the original data
					narrowband.functions[lb][n] += Numerics.bounded(levelsetForces(xyz, mgdmlabels[lb][xyz], forces[lb])
																	- levelsetForces(xyz, mgdmlabels[lb+1][xyz], bestforces[lb]),
																	 -0.9f, 0.9f);
					
					// change of sign ?
					if (narrowband.functions[lb][n]<0) {
						//if (debug) System.out.print(""+lb);
						
						// try all possible labels
						newlb = EMPTY;
						curlb = lb+1;
						while (newlb==EMPTY && curlb<=nmgdm) {
							if (mgdmlabels[curlb][xyz]!=EMPTY && homeomorphicLabeling(xyz, mgdmlabels[curlb][xyz]) )
								newlb = curlb;
							else
								curlb++;
						}
						// allow any change of topology at second and third order labels
						//if (mgdmlabels[lb+1][xyz]!=EMPTY && (lb>0 || homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz]))) {
						
						// restrict changes of labeling at all levels
						//if (mgdmlabels[lb+1][xyz]!=EMPTY && homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz]) ) {
						if (newlb!=EMPTY) {
							// for all levels
							nswap[lb]++;
							narrowband.labels[lb][n] = mgdmlabels[newlb][xyz];
							narrowband.functions[lb][n] = -narrowband.functions[lb][n];
							// never switch the last label
							if (newlb<nmgdm) narrowband.labels[newlb][n] = mgdmlabels[lb][xyz];
							// update the segmentation with first label
							if (lb==0) segmentation[xyz] = mgdmlabels[newlb][xyz];
							// check for boundary changes in the landmines : force reinitialization
							if (lb==0 && landmines.get(xyz)) reinitLM = true;
							// check for far labels getting mixed in: time to re-initialize
							if (lb==0 && narrowband.labels[0][n]==mgdmlabels[nmgdm][xyz]) reinitOL = true;
						} else {
							// reset to low value
							narrowband.functions[lb][n] = lowlevel;
						}
					}
				}
			}
			diff = (nswap[0]/(float)boundarysize);
			if (debug) for (int lb=0;lb<nmgdm;lb++) System.out.print("changed labels ("+lb+"): "+nswap[lb]+" ("+(nswap[lb]/(float)boundarysize*100.0f)+" % of boundary)\n");
			
			// once all the new values are computed, copy into original MGDM functions
			avgdiff = 0.0f;
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				// measure the changes at the base level
				if (mgdmlabels[0][xyz]==narrowband.labels[0][n]) avgdiff += Numerics.abs(mgdmfunctions[0][xyz]-narrowband.functions[0][n]);
				else avgdiff += Numerics.abs(mgdmfunctions[0][xyz]+narrowband.functions[0][n]);
				
				for (int lb=0;lb<nmgdm;lb++) {
					mgdmlabels[lb][xyz] = narrowband.labels[lb][n];
					mgdmfunctions[lb][xyz] = narrowband.functions[lb][n];
				}
			}
			if (debug) System.out.print("mean distance function change: "+(avgdiff/narrowband.currentsize)+"\n");
	
			//if (t<iter-1 && (reinitLM || reinitOL) ) {
			if (t<iter-1 && reinitLM) {
				if (debug) System.out.print("re-initialization (LM: "+reinitLM+" | OL: "+reinitOL+" )\n");
        		
				//resetIsosurfaceNarrowBand(narrowband);
				resetIsosurfaceBoundary();
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
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
		fastMarchingReinitialization(true);
		
        return;
    }
    
 	/** specific forces applied to the level sets (application dependent) */
 	private final double levelsetForces(int xyz, int lb, double balloonval) {
    	
		// simple option: rebuild level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		
		
		if (mgdmfunctions[0][xyz]==UNKNOWN) return 0.0;
		
		// do the center point first
		if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
		else  phi[1][1][1] = 0.0f;
		for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb && mgdmfunctions[l][xyz]!=UNKNOWN;l++) {
			phi[1][1][1] += mgdmfunctions[l][xyz];
		}
		// neighbors
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
    		int xyzn = xyz + i + j*nix + k*nix*niy;

    		if (mask[xyzn] && mgdmlabels[0][xyzn]!=EMPTY && mgdmfunctions[0][xyzn]!=UNKNOWN) {
				if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
				else  phi[i+1][j+1][k+1] = 0.0f;
					
				for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb && mgdmfunctions[l][xyzn]!=UNKNOWN;l++) {
					phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
				}
			} else {
				// filling in values outside the mask?? center value
				phi[i+1][j+1][k+1] = phi[1][1][1];
			}
    	}
    	// first derivatives
    	Dmx = phi[1][1][1] - phi[0][1][1];
    	Dmy = phi[1][1][1] - phi[1][0][1];
    	Dmz = phi[1][1][1] - phi[1][1][0];
    	
    	Dpx = phi[2][1][1] - phi[1][1][1];
    	Dpy = phi[1][2][1] - phi[1][1][1];
    	Dpz = phi[1][1][2] - phi[1][1][1];
    	
    	D0x = (phi[2][1][1] - phi[0][1][1])/2.0;
    	D0y = (phi[1][2][1] - phi[1][0][1])/2.0;
    	D0z = (phi[1][1][2] - phi[1][1][0])/2.0;
    	
    	// if using smoothing forces:
    	smooth = 0.0;
    	if (Numerics.abs(smoothweight)>0) {
    	
			// second derivatives
			Dxx = phi[0][1][1] + phi[2][1][1] - 2.0*phi[1][1][1];
			Dyy = phi[1][0][1] + phi[1][2][1] - 2.0*phi[1][1][1];
			Dzz = phi[1][1][0] + phi[1][1][2] - 2.0*phi[1][1][1];
			
			Dxy = (phi[0][0][1] + phi[2][2][1] - phi[0][2][1] - phi[2][0][1])/4.0;
			Dyz = (phi[1][0][0] + phi[1][2][2] - phi[1][0][2] - phi[1][2][0])/4.0;
			Dzx = (phi[0][1][0] + phi[2][1][2] - phi[2][1][0] - phi[0][1][2])/4.0;
			
			// gradient norm
			SD0x = D0x * D0x;
			SD0y = D0y * D0y;
			SD0z = D0z * D0z;
			GPhi = Math.sqrt(SD0x + SD0y + SD0z);
			
			/*
			if (curvaturedir[lb]==NONE) {
				K = 0.0;
			} else {
			*/
			// mean curvature
			K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
						- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
				
			// gaussian curvature
			G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
						+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));
	
			if(GPhi > 0.0000001){
				tmp = GPhi*GPhi;
				K = K/(GPhi*tmp);
				G = G/(tmp*tmp);
				tmp = K*K - 2*G;
				if(tmp > 0 ) K = K*G/tmp;
			} else {
				K = 0;
			}
			/*
				if (K>0 && curvaturedir[lb]==POINTED) K = 0.1*K;
				else if (K<0 && curvaturedir[lb]==ROUNDED) K = 0.1*K;
			}
			*/
			
			/*
			// curvature smoothing
			if (scaledCurv) {
				z = xyz/(nix*niy);
				y = (xyz-nix*niy*z)/nix;
				x = xyz-nix*niy*z-nix*y;
				
				val = ImageInterpolation.linearInterpolation(curvatureforce, 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
			} else {
				val = curvatureforce[xyz];
			}
			smooth = smoothweight*val*stepsize*K*GPhi;
			*/
			smooth = smoothweight*stepsize*K*GPhi;
		}
		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		balloon = 0.0;
		if (Numerics.abs(balloonweight)>0) {
			DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
									 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
									 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
			
			DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
									 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
									 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
			/*
			if (scaledBalloon) {
				z = xyz/(nix*niy);
				y = (xyz-nix*niy*z)/nix;
				x = xyz-nix*niy*z-nix*y;
				
				val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
			} else {
				val = balloonforces[lb][xyz];
			}
			balloon = balloonweight*stepsize*(Numerics.max(val, 0.0)*DeltaP + Numerics.min(val, 0.0)*DeltaM);
			*/
			balloon = balloonweight*stepsize*(Numerics.max(balloonval, 0.0)*DeltaP + Numerics.min(balloonval, 0.0)*DeltaM);
		}
		
		// external force field from target
		double field = 0.0;
		if (Numerics.abs(fieldweight)>0) {
			if (fieldtarget[xyz]!=EMPTY) {
				int tlb = fieldtarget[xyz];
				// check if in label list
				boolean process=false;
				for (int n=0;n<=nmgdm&&!process;n++) if (mgdmlabels[n][xyz]==tlb) process=true;
				
				if (process) {
					// do the center point first
					if (mgdmlabels[0][xyz]==tlb) phi[1][1][1] = -mgdmfunctions[0][xyz];
					else  phi[1][1][1] = 0.0f;
					for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=tlb && mgdmfunctions[l][xyz]!=UNKNOWN;l++) {
						phi[1][1][1] += mgdmfunctions[l][xyz];
					}
					// neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
						int xyzn = xyz + i + j*nix + k*nix*niy;
			
						if (mask[xyzn] && mgdmlabels[0][xyzn]!=EMPTY && mgdmfunctions[0][xyzn]!=UNKNOWN) {
							if (mgdmlabels[0][xyzn]==tlb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
							else  phi[i+1][j+1][k+1] = 0.0f;
								
							for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=tlb && mgdmfunctions[l][xyzn]!=UNKNOWN;l++) {
								phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
							}
						} else {
							// filling in values outside the mask?? center value
							phi[i+1][j+1][k+1] = phi[1][1][1];
						}
					}
					// first derivatives
					Dmx0 = phi[2][1][1] - phi[0][1][1];
					Dmy0 = phi[1][2][1] - phi[1][0][1];
					Dmz0 = phi[1][1][2] - phi[1][1][0];
					
					field = fieldweight*stepsize*(Dmx0*D0x + Dmy0*D0y + Dmz0*D0z);
				} else {
					// if best is not a neighbor, just shrink the label
					field = -fieldweight*stepsize*0.5f;
				}
			}
		}

		
		return - (smooth - balloon - field);
    }

 	/** specific forces applied to the level sets (application dependent) */
 	/*
 	private final void allLevelsetForces(int xyz, double[] force, byte[] bestid) {
    	
		// simple option: rebuild each level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		
		for (byte n=0;n<=nmgdm+1;n++) {
			int lb = mgdmlabels[0][xyz];
			// do the center point first
			if (n==0) phii[n][0] = -mgdmfunctions[0][xyz];
			else  phii[n][0] = 0.0f;
			for (int l=0;l<n;l++) {
				phii[n][0] += mgdmfunctions[l][xyz];
			}
			// neighbors
			for (int b=0;b<18;b++) {
				int xyzn = xyz + ngbx[b] + ngby[b]*nix + ngbz[b]*nix*niy;

				if (mask[xyzn] && mgdmlabels[0][xyzn]!=EMPTY && mgdmfunctions[0][xyzn]!=UNKNOWN) {
					if (mgdmlabels[0][xyzn]==) phii[n][b+1] = -mgdmfunctions[0][xyzn];
					else  phii[n][b+1] = 0.0f;
					
					for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb && mgdmfunctions[l][xyzn]!=UNKNOWN;l++) {
						phii[n][b+1] += mgdmfunctions[l][xyzn];
					}
				} else {
					// filling in values outside the mask?? center value
					phii[n][b+1] = phii[n][0];
				}
			}
			// first derivatives
			Dmx = phi[n][0] - phi[0][1][1];
			Dmy = phi[1][1][1] - phi[1][0][1];
			Dmz = phi[1][1][1] - phi[1][1][0];
			
			Dpx = phi[2][1][1] - phi[1][1][1];
			Dpy = phi[1][2][1] - phi[1][1][1];
			Dpz = phi[1][1][2] - phi[1][1][1];
			
			Dx[n] = (phi[2][1][1] - phi[0][1][1])/2.0;
			Dy[n] = (phi[1][2][1] - phi[1][0][1])/2.0;
			Dz[n] = (phi[1][1][2] - phi[1][1][0])/2.0;
    	
    	// if using smoothing forces:
    	smooth = 0.0;
    	if (Numerics.abs(smoothweight)>0) {
    	
			// second derivatives
			Dxx = phi[0][1][1] + phi[2][1][1] - 2.0*phi[1][1][1];
			Dyy = phi[1][0][1] + phi[1][2][1] - 2.0*phi[1][1][1];
			Dzz = phi[1][1][0] + phi[1][1][2] - 2.0*phi[1][1][1];
			
			Dxy = (phi[0][0][1] + phi[2][2][1] - phi[0][2][1] - phi[2][0][1])/4.0;
			Dyz = (phi[1][0][0] + phi[1][2][2] - phi[1][0][2] - phi[1][2][0])/4.0;
			Dzx = (phi[0][1][0] + phi[2][1][2] - phi[2][1][0] - phi[0][1][2])/4.0;
			
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
	
			if(GPhi > 0.0000001){
				tmp = GPhi*GPhi;
				K = K/(GPhi*tmp);
				G = G/(tmp*tmp);
				tmp = K*K - 2*G;
				if(tmp > 0 ) K = K*G/tmp;
			} else {
				K = 0;
			}
			
			smooth = smoothweight*stepsize*K*GPhi;
		}
		
		// external force field from target
		double field = 0.0;
		if (Numerics.abs(fieldweight)>0) {
			if (fieldtarget[xyz]!=EMPTY) {
				int tlb = fieldtarget[xyz];
				// check if in label list
				boolean process=false;
				for (int n=0;n<=nmgdm&&!process;n++) if (mgdmlabels[n][xyz]==tlb) process=true;
				
				if (process) {
					// do the center point first
					if (mgdmlabels[0][xyz]==tlb) phi[1][1][1] = -mgdmfunctions[0][xyz];
					else  phi[1][1][1] = 0.0f;
					for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=tlb && mgdmfunctions[l][xyz]!=UNKNOWN;l++) {
						phi[1][1][1] += mgdmfunctions[l][xyz];
					}
					// neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
						int xyzn = xyz + i + j*nix + k*nix*niy;
			
						if (mask[xyzn] && mgdmlabels[0][xyzn]!=EMPTY && mgdmfunctions[0][xyzn]!=UNKNOWN) {
							if (mgdmlabels[0][xyzn]==tlb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
							else  phi[i+1][j+1][k+1] = 0.0f;
								
							for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=tlb && mgdmfunctions[l][xyzn]!=UNKNOWN;l++) {
								phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
							}
						} else {
							// filling in values outside the mask?? center value
							phi[i+1][j+1][k+1] = phi[1][1][1];
						}
					}
					// first derivatives
					Dmx0 = phi[2][1][1] - phi[0][1][1];
					Dmy0 = phi[1][2][1] - phi[1][0][1];
					Dmz0 = phi[1][1][2] - phi[1][1][0];
					
					field = fieldweight*stepsize*(Dmx0*D0x + Dmy0*D0y + Dmz0*D0z);
				} else {
					// if best is not a neighbor, just shrink the label
					field = -fieldweight*stepsize*0.5f;
				}
			}
		}

		
		return - (smooth - field);
    }

	/** specific forces applied to the level sets (application dependent) */
    private final double levelsetForces(int xyz, int lb) {
    	
		// simple option: rebuild level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		float[][][] phi = new float[3][3][3];
		
		// do the center point first
		if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
		else  phi[1][1][1] = 0.0f;
		for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb;l++) {
			phi[1][1][1] += mgdmfunctions[l][xyz];
		}
		// neighbors
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
    		int xyzn = xyz + i + j*nix + k*nix*niy;

			if (mask[xyzn]) {
				if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
				else  phi[i+1][j+1][k+1] = 0.0f;
					
				for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb;l++) {
					phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
				}
			} else {
				// filling in values outside the mask?? center value
				phi[i+1][j+1][k+1] = phi[1][1][1];
			}
    	}
    	// first derivatives
    	double Dmx = phi[1][1][1] - phi[0][1][1];
    	double Dmy = phi[1][1][1] - phi[1][0][1];
    	double Dmz = phi[1][1][1] - phi[1][1][0];
    	
    	double Dpx = phi[2][1][1] - phi[1][1][1];
    	double Dpy = phi[1][2][1] - phi[1][1][1];
    	double Dpz = phi[1][1][2] - phi[1][1][1];
    	
    	double D0x = (phi[2][1][1] - phi[0][1][1])/2.0;
    	double D0y = (phi[1][2][1] - phi[1][0][1])/2.0;
    	double D0z = (phi[1][1][2] - phi[1][1][0])/2.0;
    	
    	// if using smoothing forces:
    	double smooth = 0.0;
    	if (Numerics.abs(smoothweight)>0) {
    	
			// second derivatives
			double Dxx = phi[0][1][1] + phi[2][1][1] - 2.0*phi[1][1][1];
			double Dyy = phi[1][0][1] + phi[1][2][1] - 2.0*phi[1][1][1];
			double Dzz = phi[1][1][0] + phi[1][1][2] - 2.0*phi[1][1][1];
			
			double Dxy = (phi[0][0][1] + phi[2][2][1] - phi[0][2][1] - phi[2][0][1])/4.0;
			double Dyz = (phi[1][0][0] + phi[1][2][2] - phi[1][0][2] - phi[1][2][0])/4.0;
			double Dzx = (phi[0][1][0] + phi[2][1][2] - phi[2][1][0] - phi[0][1][2])/4.0;
			
			// gradient norm
			double SD0x = D0x * D0x;
			double SD0y = D0y * D0y;
			double SD0z = D0z * D0z;
			double GPhi = Math.sqrt(SD0x + SD0y + SD0z);
			
			// mean curvature
			double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
						- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
				
			// gaussian curvature
			double G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
						+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));
	
			if(GPhi > 0.0000001){
				double tmp = GPhi*GPhi;
				K = K/(GPhi*tmp);
				G = G/(tmp*tmp);
				tmp = K*K - 2*G;
				if(tmp > 0 ) K = K*G/tmp;
			} else {
				K = 0;
			}
			
			// curvature smoothing
			smooth = smoothweight*stepsize*K*GPhi;
		}
		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		double balloon = 0.0;
		if (Numerics.abs(balloonweight)>0) {
			double DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
									 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
									 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
			
			double DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
									 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
									 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
			
			balloon = balloonweight*stepsize*(Numerics.max(balloonforces[lb][xyz], 0.0)*DeltaP + Numerics.min(balloonforces[lb][xyz], 0.0)*DeltaM);
		}
				
		return - (smooth - balloon);
    }
    
    private final double balloonForce(int lb, int xyz) {
    	if (scaledBalloon) {
			z = xyz/(nix*niy);
			y = (xyz-nix*niy*z)/nix;
			x = xyz-nix*niy*z-nix*y;
				
			val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
		} else {
			val = balloonforces[lb][xyz];
		}
		return val;
	}
    
 	/** specific forces applied to the level sets (application dependent) */
 	/*
    private final double levelsetSpeed(int xyz, int lb, double[][] forces, int lbmax, int lbsec) {
    	
		// simple option: rebuild level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		//float[][][] phi = new float[3][3][3];
		
		// do the center point first
		if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
		else  phi[1][1][1] = 0.0f;
		for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb && mgdmlabels[l][xyz]!=EMPTY;l++) {
			phi[1][1][1] += mgdmfunctions[l][xyz];
		}
		// neighbors
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
    		xyzn = xyz + i + j*nix + k*nix*niy;

			if (mask[xyzn]) {
				if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
				else  phi[i+1][j+1][k+1] = 0.0f;
					
				for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb && mgdmlabels[l][xyzn]!=EMPTY;l++) {
					phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
				}
			} else {
				// filling in values outside the mask?? center value
				phi[i+1][j+1][k+1] = phi[1][1][1];
			}
    	}
    	// first derivatives
    	Dmx = phi[1][1][1] - phi[0][1][1];
    	Dmy = phi[1][1][1] - phi[1][0][1];
    	Dmz = phi[1][1][1] - phi[1][1][0];
    	
    	Dpx = phi[2][1][1] - phi[1][1][1];
    	Dpy = phi[1][2][1] - phi[1][1][1];
    	Dpz = phi[1][1][2] - phi[1][1][1];
    	
    	D0x = (phi[2][1][1] - phi[0][1][1])/2.0;
    	D0y = (phi[1][2][1] - phi[1][0][1])/2.0;
    	D0z = (phi[1][1][2] - phi[1][1][0])/2.0;
    	
		// gradient norm
		SD0x = D0x * D0x;
		SD0y = D0y * D0y;
		SD0z = D0z * D0z;
		GPhi = Math.sqrt(SD0x + SD0y + SD0z);
		
		// second derivatives
		Dxx = phi[0][1][1] + phi[2][1][1] - 2.0*phi[1][1][1];
		Dyy = phi[1][0][1] + phi[1][2][1] - 2.0*phi[1][1][1];
		Dzz = phi[1][1][0] + phi[1][1][2] - 2.0*phi[1][1][1];
		
		Dxy = (phi[0][0][1] + phi[2][2][1] - phi[0][2][1] - phi[2][0][1])/4.0;
		Dyz = (phi[1][0][0] + phi[1][2][2] - phi[1][0][2] - phi[1][2][0])/4.0;
		Dzx = (phi[0][1][0] + phi[2][1][2] - phi[2][1][0] - phi[0][1][2])/4.0;

		// mean curvature
		K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
					- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
			
		// gaussian curvature
		G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
					+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));

		if(GPhi > 0.0000001){
			tmp = GPhi*GPhi;
			K = K/(GPhi*tmp);
			G = G/(tmp*tmp);
			tmp = K*K - 2*G;
			if(tmp > 0 ) K = K*G/tmp;
		} else {
			K = 0;
		}
		
		// curvature smoothing
		smooth = smoothweight*stepsize*K*GPhi;
		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		balloon = 0.0;
		if (Numerics.abs(balloonweight)+Numerics.abs(fieldweight)>0) {
			pos = 0.0; neg = 0.0;
			npos = 0; nneg = 0;
			
			for (int n=0;n<nmgdm;n++) if (mgdmlabels[n][xyz]==lb) {
				// current
				curr = forces[n][W];
				// next
				next = forces[n+1][W];
				
				/*
				if (curr>0 && next<0) {
					//balloon = balloonweight*stepsize;
					balloon = 1;
				} else if (curr<0 && next>0) {
					//balloon = -balloonweight*stepsize;
					balloon = 2;
				} else if (curr>0 && next>0) {
					//balloon = balloonweight*stepsize*(curr-next);
					balloon = 3;
				} else if (curr<0 && next<0) {
					/*
					// use stronger force
					if (lbmax!=EMPTY && forces[lbmax][W]>0) {
						val = forces[lbmax][X]*D0x + forces[lbmax][Y]*D0y + forces[lbmax][Z]*D0z;
						if (val>0) balloon = balloonweight*stepsize;
						else balloon = -balloonweight*stepsize;
					}	
					*/
					/*
					balloon = 4;
				}
				*/
				/*
				// for testing: remove main force...
				balloon = balloonweight*stepsize*(curr-next);
				//balloon = balloonweight*stepsize*(curr-next);
				// use stronger force
				if (lbmax!=EMPTY && forces[lbmax][W]>0) {
				//if (lbmax!=EMPTY) {
					vc = forces[lbmax][X]*forces[n][X] + forces[lbmax][Y]*forces[n][Y] + forces[lbmax][Z]*forces[n][Z];
					vn = forces[lbmax][X]*forces[n+1][X] + forces[lbmax][Y]*forces[n+1][Y] + forces[lbmax][Z]*forces[n+1][Z];
					balloon += fieldweight*stepsize*forces[lbmax][W]*(vc-vn);
				}
				
			}

			/*
			for (int n=0;n<nmgdm;n++) if (mgdmlabels[n][xyz]==lb) {
				// current
				val = forces[n][X]*D0x + forces[n][Y]*D0y + forces[n][Z]*D0z;
				if (GPhi>1) val = val/GPhi;
				balloon += 0.5*balloonweight*stepsize*forces[n][W]*val;
				// next
				val = forces[n+1][X]*D0x + forces[n+1][Y]*D0y + forces[n+1][Z]*D0z;
				if (GPhi>1) val = val/GPhi;
				balloon += 0.5*balloonweight*stepsize*forces[n+1][W]*val;
			}
			*/
			
			/* not bad at all, but slow
			for (int n=0;n<nmgdm;n++) if (mgdmlabels[n][xyz]!=EMPTY) {
				val = forces[n][X]*D0x + forces[n][Y]*D0y + forces[n][Z]*D0z;
				if (GPhi>1) val = val/GPhi;
				if (val>0) {
					pos += forces[n][W]*val;
					npos++;
				} else if (val<0) {
					neg += forces[n][W]*val;
					nneg++;
				}
			}
			if (otherlabels[xyz]!=EMPTY) {
				val = forces[nmgdm][X]*D0x + forces[nmgdm][Y]*D0y + forces[nmgdm][Z]*D0z;
				if (GPhi>1) val = val/GPhi;
				if (val>0) {
					pos += forces[nmgdm][W]*val;
					npos++;
				} else if (val<0) {
					neg += forces[nmgdm][W]*val;
					nneg++;
				}
			}
			
			if (npos>0) balloon += 0.5*balloonweight*stepsize*pos/npos;
			if (nneg>0) balloon += 0.5*balloonweight*stepsize*neg/nneg;
			
			//if (npos==0) System.out.print("o");
			*/
			
			/*
			for (int n=0;n<nmgdm;n++) if (mgdmlabels[n][xyz]!=EMPTY) {
				/*
				if (n==0) weight = 1.0;
				else if (n==1) weight = 1.0;
				else if (n==2) weight = 0.5;
				else if (n==3) weight = 0.25;
				else weight = 0.0;
				*/
				/*
				weight = 1.0;
				// weighting??
				/*
				balloon += weight*balloonweight*stepsize*(Numerics.max(forces[n][X],0.0)*Dmx + Numerics.min(forces[n][X],0.0)*Dpx
													+Numerics.max(forces[n][Y],0.0)*Dmy + Numerics.min(forces[n][Y],0.0)*Dpy
											   		+Numerics.max(forces[n][Z],0.0)*Dmz + Numerics.min(forces[n][Z],0.0)*Dpz);
				*/
				/*
				balloon += weight*balloonweight*stepsize*(forces[n][X]*D0x + forces[n][Y]*D0y + forces[n][Z]*D0z);
			}
			*/
			// balloon force for this object
			/*
			double DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
									 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
									 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
			
			double DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
									 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
									 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
			double val;
			if (scaledBalloon) {
				int z = xyz/(nix*niy);
				int y = (xyz-nix*niy*z)/nix;
				int x = xyz-nix*niy*z-nix*y;
				
				val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
			} else {
				val = balloonforces[lb][xyz];
			}
			balloon = balloonweight*stepsize*(Numerics.max(val, 0.0)*DeltaP + Numerics.min(val, 0.0)*DeltaM);
			*/
			
			/*
			for (int n=0;n<nmgdm;n++) if (mgdmlabels[n][xyz]==lb) {
				balloon = balloonweight*stepsize*(Numerics.max(forces[n][X],0.0)*Dmx + Numerics.min(forces[n][X],0.0)*Dpx
												  +Numerics.max(forces[n][Y],0.0)*Dmy + Numerics.min(forces[n][Y],0.0)*Dpy
												  +Numerics.max(forces[n][Z],0.0)*Dmz + Numerics.min(forces[n][Z],0.0)*Dpz);
			}			
			// balloon force for the largest object
			if (lbmax!=EMPTY && lbmax<nmgdm && lbsec!=EMPTY && mgdmlabels[lbmax][xyz]==lb) {
				balloon += balloonweight*stepsize*(Numerics.max(forces[lbsec][X],0.0)*Dmx + Numerics.min(forces[lbsec][X],0.0)*Dpx
													+Numerics.max(forces[lbsec][Y],0.0)*Dmy + Numerics.min(forces[lbsec][Y],0.0)*Dpy
											   		+Numerics.max(forces[lbsec][Z],0.0)*Dmz + Numerics.min(forces[lbsec][Z],0.0)*Dpz);
			} else if (lbmax!=EMPTY) {
				balloon += balloonweight*stepsize*(Numerics.max(forces[lbmax][X],0.0)*Dmx + Numerics.min(forces[lbmax][X],0.0)*Dpx
													+Numerics.max(forces[lbmax][Y],0.0)*Dmy + Numerics.min(forces[lbmax][Y],0.0)*Dpy
											   		+Numerics.max(forces[lbmax][Z],0.0)*Dmz + Numerics.min(forces[lbmax][Z],0.0)*Dpz);				
			}
			*/
			/*
			if (lbmax!=EMPTY) {
				balloon = balloonweight*stepsize*(Numerics.max(forces[lbmax][X],0.0)*Dmx + Numerics.min(forces[lbmax][X],0.0)*Dpx
												  +Numerics.max(forces[lbmax][Y],0.0)*Dmy + Numerics.min(forces[lbmax][Y],0.0)*Dpy
											   	  +Numerics.max(forces[lbmax][Z],0.0)*Dmz + Numerics.min(forces[lbmax][Z],0.0)*Dpz);				
			}
			*/
			/*
			if (lbmax!=EMPTY) {
				balloon = balloonweight*stepsize*(forces[lbmax][X]*D0x + forces[lbmax][Y]*D0y + forces[lbmax][Z]*D0z);				
			}
			*/
			/*
		}
		
		// external force field
		field = 0.0;
		/*
		if (Numerics.abs(fieldweight)>0) {
			double vx,vy,vz;
			if (scaledField) {
				int z = xyz/(nix*niy);
				int y = (xyz-nix*niy*z)/nix;
				int x = xyz-nix*niy*z-nix*y;
				
				vx = ImageInterpolation.linearInterpolation(fieldforce[X], 0.0f, x*rix/rfx, y*riy/rfy, z*riz/rfz, nfx, nfy, nfz);
				vy = ImageInterpolation.linearInterpolation(fieldforce[Y], 0.0f, x*rix/rfx, y*riy/rfy, z*riz/rfz, nfx, nfy, nfz);
				vz = ImageInterpolation.linearInterpolation(fieldforce[Z], 0.0f, x*rix/rfx, y*riy/rfy, z*riz/rfz, nfx, nfy, nfz);
			} else {
				vx = fieldforce[X][xyz];
				vy = fieldforce[Y][xyz];
				vz = fieldforce[Z][xyz];
			}
			
			field =  fieldweight*stepsize*(Numerics.max(vx,0.0)*Dmx + Numerics.min(vx,0.0)*Dpx
										   +Numerics.max(vy,0.0)*Dmy + Numerics.min(vy,0.0)*Dpy
										   +Numerics.max(vz,0.0)*Dmz + Numerics.min(vz,0.0)*Dpz);
		}
		*/
		/*
		// lower forces far from the boundary
		boundaryfactor = 1.0/(1.0+phi[1][1][1]*phi[1][1][1]/sigmaB);
		//boundaryfactor = 1.0;
		//if (phi[1][1][1]*phi[1][1][1]>=landmineDist) boundaryfactor = 0.0;
		
		return boundaryfactor*(balloon + field) - smooth;
    }


	/** specific forces applied to the level sets (application dependent) */
	/*
    private final void directedBalloonForces(int xyz, double[][] forces) {
    	// first, the active labels
    	for (int n=0;n<nmgdm;n++) {
    		lb = mgdmlabels[n][xyz];
    		if (lb!=EMPTY) {
 				// simple option: rebuild level set locally
				// note: we go back to the convention of usual level sets with negative value inside, positive value outside
				
				// do the center point first
				if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
				else  phi[1][1][1] = 0.0f;
				for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb;l++) {
					phi[1][1][1] += mgdmfunctions[l][xyz];
				}
				// neighbors
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
					xyzn = xyz + i + j*nix + k*nix*niy;
		
					if (mask[xyzn]) {
						if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
						else  phi[i+1][j+1][k+1] = 0.0f;
							
						for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb;l++) {
							phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
						}
					} else {
						// filling in values outside the mask?? center value
						phi[i+1][j+1][k+1] = phi[1][1][1];
					}
				}
				// first derivatives
				Dmx = phi[1][1][1] - phi[0][1][1];
				Dmy = phi[1][1][1] - phi[1][0][1];
				Dmz = phi[1][1][1] - phi[1][1][0];
				
				Dpx = phi[2][1][1] - phi[1][1][1];
				Dpy = phi[1][2][1] - phi[1][1][1];
				Dpz = phi[1][1][2] - phi[1][1][1];
				
				// external object-specific balloon forces
				// we assume b > 0 inside objects, b<0 outside, as in membership functions
	
				DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
										 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
										 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
				
				DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
										 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
										 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
				if (scaledBalloon) {
					z = xyz/(nix*niy);
					y = (xyz-nix*niy*z)/nix;
					x = xyz-nix*niy*z-nix*y;
					
					val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
				} else {
					val = balloonforces[lb][xyz];
				}
				
				forces[n][X] = Numerics.max(val,0.0)*(Numerics.max(Dmx, 0.0) + Numerics.min(Dpx, 0.0))/DeltaP
							  + Numerics.min(val,0.0)*(Numerics.max(Dpx, 0.0) + Numerics.min(Dmx, 0.0))/DeltaM;
				forces[n][Y] = Numerics.max(val,0.0)*(Numerics.max(Dmy, 0.0) + Numerics.min(Dpy, 0.0))/DeltaP
							  + Numerics.min(val,0.0)*(Numerics.max(Dpy, 0.0) + Numerics.min(Dmy, 0.0))/DeltaM;
				forces[n][Z] = Numerics.max(val,0.0)*(Numerics.max(Dmz, 0.0) + Numerics.min(Dpz, 0.0))/DeltaP
							  + Numerics.min(val,0.0)*(Numerics.max(Dpz, 0.0) + Numerics.min(Dmz, 0.0))/DeltaM;
							  
				forces[n][W] = Numerics.max(val, 0.0)*DeltaP + Numerics.min(val, 0.0)*DeltaM;
			}
		}
		// repeat for the outside label
		lb = otherlabels[xyz];
		if (lb!=EMPTY) {
			// simple option: rebuild level set locally
			// note: we go back to the convention of usual level sets with negative value inside, positive value outside
			
			// do the center point first
			if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
			else  phi[1][1][1] = 0.0f;
			for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb;l++) {
				phi[1][1][1] += mgdmfunctions[l][xyz];
			}
			// neighbors
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
				xyzn = xyz + i + j*nix + k*nix*niy;
	
				if (mask[xyzn]) {
					if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
					else  phi[i+1][j+1][k+1] = 0.0f;
						
					for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb;l++) {
						phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
					}
				} else {
					// filling in values outside the mask?? center value
					phi[i+1][j+1][k+1] = phi[1][1][1];
				}
			}
			// first derivatives
			Dmx = phi[1][1][1] - phi[0][1][1];
			Dmy = phi[1][1][1] - phi[1][0][1];
			Dmz = phi[1][1][1] - phi[1][1][0];
			
			Dpx = phi[2][1][1] - phi[1][1][1];
			Dpy = phi[1][2][1] - phi[1][1][1];
			Dpz = phi[1][1][2] - phi[1][1][1];
			
			// external object-specific balloon forces
			// we assume b > 0 inside objects, b<0 outside, as in membership functions

			DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
									 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
									 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
			
			DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
									 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
									 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
			if (scaledBalloon) {
				z = xyz/(nix*niy);
				y = (xyz-nix*niy*z)/nix;
				x = xyz-nix*niy*z-nix*y;
				
				val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
			} else {
				val = balloonforces[lb][xyz];
			}
			
			forces[nmgdm][X] = Numerics.max(val,0.0)*(Numerics.max(Dmx, 0.0) + Numerics.min(Dpx, 0.0))/DeltaP
						  	  + Numerics.min(val,0.0)*(Numerics.max(Dpx, 0.0) + Numerics.min(Dmx, 0.0))/DeltaM;
			forces[nmgdm][Y] = Numerics.max(val,0.0)*(Numerics.max(Dmy, 0.0) + Numerics.min(Dpy, 0.0))/DeltaP
						  	  + Numerics.min(val,0.0)*(Numerics.max(Dpy, 0.0) + Numerics.min(Dmy, 0.0))/DeltaM;
			forces[nmgdm][Z] = Numerics.max(val,0.0)*(Numerics.max(Dmz, 0.0) + Numerics.min(Dpz, 0.0))/DeltaP
						  	  + Numerics.min(val,0.0)*(Numerics.max(Dpz, 0.0) + Numerics.min(Dmz, 0.0))/DeltaM;
						  
			forces[nmgdm][W] = Numerics.max(val, 0.0)*DeltaP + Numerics.min(val, 0.0)*DeltaM;
		}

		return;
    }
    */
    
	/** specific forces applied to the level sets (application dependent) */
    private final void centralBalloonForces(int xyz, double[][] forces) {
    	// first, the active labels
    	for (int n=0;n<=nmgdm;n++) {
    		int lb = mgdmlabels[n][xyz];
    		if (lb!=EMPTY) {
 				// simple option: rebuild level set locally
				// note: we go back to the convention of usual level sets with negative value inside, positive value outside
				
				// do the center point first
				if (mgdmlabels[0][xyz]==lb) phi[1][1][1] = -mgdmfunctions[0][xyz];
				else  phi[1][1][1] = 0.0f;
				for (int l=0;l<nmgdm && mgdmlabels[l][xyz]!=lb && mgdmlabels[l][xyz]!=EMPTY;l++) {
					phi[1][1][1] += mgdmfunctions[l][xyz];
				}
				// neighbors
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i>0 || j*j>0 || k*k>0) {
					int xyzn = xyz + i + j*nix + k*nix*niy;
		
					if (mask[xyzn]) {
						if (mgdmlabels[0][xyzn]==lb) phi[i+1][j+1][k+1] = -mgdmfunctions[0][xyzn];
						else  phi[i+1][j+1][k+1] = 0.0f;
							
						for (int l=0;l<nmgdm && mgdmlabels[l][xyzn]!=lb && mgdmlabels[l][xyzn]!=EMPTY;l++) {
							phi[i+1][j+1][k+1] += mgdmfunctions[l][xyzn];
						}
					} else {
						// filling in values outside the mask?? center value
						phi[i+1][j+1][k+1] = phi[1][1][1];
					}
				}
				// first derivatives
				D0x = (phi[2][1][1] - phi[0][1][1])/2.0;
				D0y = (phi[1][2][1] - phi[1][0][1])/2.0;
				D0z = (phi[1][1][2] - phi[1][1][0])/2.0;
    	
				// gradient norm
				SD0x = D0x * D0x;
				SD0y = D0y * D0y;
				SD0z = D0z * D0z;
				GPhi = Numerics.max(1e-16, Math.sqrt(SD0x + SD0y + SD0z));
    		
				if (scaledBalloon) {
					z = xyz/(nix*niy);
					y = (xyz-nix*niy*z)/nix;
					x = xyz-nix*niy*z-nix*y;
					
					val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
				} else {
					val = balloonforces[lb][xyz];
				}
				/*
				forces[n][X] = val*D0x/GPhi;
				forces[n][Y] = val*D0y/GPhi;
				forces[n][Z] = val*D0z/GPhi;
				*/
				forces[n][X] = D0x/GPhi;
				forces[n][Y] = D0y/GPhi;
				forces[n][Z] = D0z/GPhi;
				forces[n][W] = val;
			} else {
				forces[n][X] = 0.0;
				forces[n][Y] = 0.0;
				forces[n][Z] = 0.0;
				forces[n][W] = -1.0;
			}
		}
		return;
    }

	/** specific forces applied to the level sets (application dependent) */
    private final void basicBalloonForces(int xyz, double[] forces, double[] best) {
    	// first, the active labels
    	for (int n=0;n<=nmgdm;n++) {
    		int lb = mgdmlabels[n][xyz];
    		if (lb!=EMPTY) {
				if (scaledBalloon) {
					z = xyz/(nix*niy);
					y = (xyz-nix*niy*z)/nix;
					x = xyz-nix*niy*z-nix*y;
					
					val = ImageInterpolation.linearInterpolation(balloonforces[lb], 0.0f, x*rix/rbx, y*riy/rby, z*riz/rbz, nbx, nby, nbz);
				} else {
					val = balloonforces[lb][xyz];
				}
				forces[n] = val;
			} else {
				forces[n] = 0.0;
			}
		}
    	for (int n=0;n<=nmgdm;n++) {
    		best[n] = -1.0;
    		for (int m=n+1;m<=nmgdm;m++) if (forces[m]>best[n]) {
    			best[n] = forces[m];
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
		
		boolean [][][] obj = new boolean[3][3][3];
		
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
		int  Nconfiguration = 0;
		short[] lbs = new short[26];
		for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[xyz+i+j*nix+l*nix*niy]!=lb) 
				&& (segmentation[xyz+i+j*nix+l*nix*niy]!=segmentation[xyz]) ) {
				boolean found = false;
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
	
    /** 
    *   isosurface distance re-initialization inside the narrow band 
    *	(assumes a correct narrow band structure enclosing the boundariy)
    */
    private final void resetIsosurfaceNarrowBand(NarrowBand narrowband) {
    	if (debug) System.out.print("level set evolution: iso-surface reinit\n");

    	float[] nbdist = new float[6];
    	boolean[] nbflag = new boolean[6];
    	boolean boundariy;
    	
		for (int n=0; n<narrowband.currentsize;n++) {
			xyz = narrowband.id[n];
			
			boundariy = false;
			for (int l=0; l<6; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
				if (mgdmlabels[0][xyznb]!=mgdmlabels[0][xyz] && mask[xyznb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(mgdmfunctions[0][xyznb]);
					nbflag[l] = true;
					boundariy = true;
				}
			}
			if (boundariy) {
				narrowband.functions[0][n] = isoSurfaceDistance(mgdmfunctions[0][xyz], nbdist, nbflag);
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int n=0; n<narrowband.currentsize;n++) {
			xyz = narrowband.id[n];
			mgdmfunctions[0][xyz] = narrowband.functions[0][n];
		}
			
        return;
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
				
				xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
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
    
    /** 
    *  	Evolution using a joint fast marching scheme
    */
    /*
    public final void evolveFastMarching(int iter, int sign, float offset) {
    	
    	if (debug) System.out.print("level set evolution: fast marching\n");

    	// computation variables
        boolean[] processed = new boolean[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
    	
        // init
    	
    	// first estimate the binariy tree size
    	int size = 0;
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) if (mgdmfunctions[0][xyz]<=0.5f) size++;
		// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	
		for (int t=0;t<iter;t++) {
			if (debug) System.out.print("iteration "+t+"\n");

			
			if (debug) System.out.print("init ("+size+")\n");
        
			heap.reset();
			double curr, next;
			float speed = 0;
			int nbound = 0;
			for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
				// the criterion for being in the narrow band is to have a short distance to closest boundaries
				if (mgdmfunctions[0][xyz]<1.0f && mgdmlabels[0][xyz]!=EMPTY && mgdmlabels[1][xyz]!=EMPTY) {
					// add to the heap
					curr = balloonweight*balloonforces[mgdmlabels[0][xyz]][xyz];
					next = balloonweight*balloonforces[mgdmlabels[1][xyz]][xyz];
					if (next>curr+offset) {
						speed =  (float)(1.0/smoothweight);
							 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
						else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
						else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
						
						heap.addValue(speed, xyz, mgdmlabels[1][xyz]);
					}
				
					nbound++;
				}
				processed[xyz] = false;
			}
			System.out.println("initial boundariy points: "+nbound);
				
			// grow the labels and functions			
			float dist = 0.0f;
			while (heap.isNotEmpty()) {
				// extract point with minimum distance
				dist = heap.getFirst();
				int xyz = heap.getFirstId();
				byte lb = heap.getFirstState();
				heap.removeFirst();
	
				// if more than nmgdm labels have been found already, this is done
				if (processed[xyz])  continue;
				
				// check topology
				if (!homeomorphicLabeling(xyz, lb)) continue; 
				
				// update the segmentation and distance function (?) at the current level
				mgdmfunctions[0][xyz] = dist;
				mgdmlabels[0][xyz] = lb;
				processed[xyz] = true;
				segmentation[xyz] = lb;
				//System.out.print(".");
				
				// find new neighbors
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					
					if (mask[xyzn]) {
						// must be in outside the object or its processed neighborhood
						if (!processed[xyzn]) {
							if (mgdmlabels[0][xyzn]!=lb && mgdmlabels[0][xyzn]!=EMPTY) {
								curr = balloonweight*balloonforces[mgdmlabels[0][xyzn]][xyzn];
								next = balloonweight*balloonforces[lb][xyzn];
								if (next>curr+offset) {
									// add to the heap
									speed =  (float)(1.0/smoothweight);
										 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
									else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
									else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
									
									//System.out.print("+");
									heap.addValue(dist+speed, xyzn, lb);
								}
							}
						}
					}
				}
			}
			for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (!processed[xyz]) {
				mgdmfunctions[0][xyz] = dist;
			}
			// re-init
			resetIsosurfaceBoundary();
			fastMarchingReinitialization(true);
		}
        return;
    }
    */
    
    /** 
    *  	Evolution using a joint fast marching scheme
    */
    /*
    public final void evolveFastMarchingProbas(int iter, float p0, float pmin) {
    	
    	if (debug) System.out.print("level set evolution: fast marching probas\n");

    	// computation variables
        boolean[] processed = new boolean[nix*niy*niz]; // note: using a byte instead of boolean for the second pass
    	
        // init
    	
    	// first estimate the binariy tree size
    	int size = 0;
		for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) if (mgdmfunctions[0][xyz]<=0.5f) size++;
		// create the tree with initial estimates of size
    	heap = new BinaryHeap2D(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MAXTREE);
    	
		for (int t=0;t<iter;t++) {
			if (debug) System.out.print("iteration "+t+"\n");

			
			if (debug) System.out.print("init ("+size+")\n");
        
			heap.reset();
			double curr, next;
			float speed = 0;
			int nbound = 0;
			for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (mask[xyz]) {
				// the criterion for being in the narrow band is to have a short distance to closest boundaries
				if (mgdmfunctions[0][xyz]<1.0f && mgdmlabels[0][xyz]!=EMPTY && mgdmlabels[1][xyz]!=EMPTY) {
					// add to the heap
					curr = balloonweight*balloonforces[mgdmlabels[0][xyz]][xyz];
					next = balloonweight*balloonforces[mgdmlabels[1][xyz]][xyz];
					
					/*
					if (next>curr) {
						speed =  (float)(1.0/smoothweight);
							 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
						else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
						else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
						
						heap.addValue(speed, xyz, mgdmlabels[1][xyz]);
					}
					*/
					/*
					speed = (float)(p0*(2.0+next-curr)/(2.0-(next-curr)*(1-2.0*p0) ) );
					
					if (speed>pmin) heap.addValue(speed, xyz, mgdmlabels[1][xyz]);
					nbound++;
				}
				processed[xyz] = false;
			}
			System.out.println("initial boundariy points: "+nbound);
				
			// grow the labels and functions			
			float dist = 0.0f;
			while (heap.isNotEmpty()) {
				// extract point with minimum distance
				dist = heap.getFirst();
				int xyz = heap.getFirstId();
				byte lb = heap.getFirstState();
				heap.removeFirst();
	
				// if more than nmgdm labels have been found already, this is done
				if (processed[xyz])  continue;
				
				// check topology
				if (!homeomorphicLabeling(xyz, lb)) continue; 
				
				// update the segmentation and distance function (?) at the current level
				mgdmfunctions[0][xyz] = dist;
				mgdmlabels[0][xyz] = lb;
				processed[xyz] = true;
				segmentation[xyz] = lb;
				//System.out.print(".");
				
				// find new neighbors
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					
					if (mask[xyzn]) {
						// must be in outside the object or its processed neighborhood
						if (!processed[xyzn]) {
							if (mgdmlabels[0][xyzn]!=lb && mgdmlabels[0][xyzn]!=EMPTY) {
								curr = balloonweight*balloonforces[mgdmlabels[0][xyzn]][xyzn];
								next = balloonweight*balloonforces[lb][xyzn];
								/*
								if (next>curr) {
									// add to the heap
									speed =  (float)(1.0/smoothweight);
										 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
									else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
									else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
									
									//System.out.print("+");
									heap.addValue(dist+speed, xyzn, lb);
								}
								*/
								/*
								speed = (float)(p0*(2.0+next-curr)/(2.0-(next-curr)*(1-2.0*p0) ) );
					
								if (dist*speed>pmin) heap.addValue(dist*speed, xyzn, lb);
							}
						}
					}
				}
			}
			for (int xyz = 0; xyz<nix*niy*niz; xyz++) if (!processed[xyz]) {
				mgdmfunctions[0][xyz] = dist;
			}
			// re-init
			//resetIsosurfaceBoundariy();
			//fastMarchingReinitialization(false);
		}
        return;
    }
    */
    private final float fastMarchingForces(int xyz, byte lb) {
    	if (balloonforces[lb][xyz]>0) return 1.0f/(float)(smoothweight+balloonweight*balloonforces[lb][xyz]);
    	else return 1.0f/(float)smoothweight;
    }
}

