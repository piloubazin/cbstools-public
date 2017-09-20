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
 *  This algorithm uses the classical GDM framework to evolve a single label
 *	according to internal (curvature) and external (vector field) forces
 *	note: for now, we are assuming isotropic voxels, and boundary conditions are not enforced
 *	(not a problem unless the narrowband reaches a boundary of the image, 
 *	mirroring conditions would solve this problem if needed)
 *
 *	@version    Jul 2011
 *	@author     Marcel Weiss 
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class Mp2rageCortexHgdm {
	
	// object types

	private	static	final	byte	EMPTY = -1;
	private	static	final	byte	OBJ = 1;
	private	static	final	byte	BG = 0;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    	
	private final static byte NORTH 	= 0;
    private final static byte SOUTH 	= 1;
	private final static byte EAST 		= 2;
    private final static byte WEST 		= 3;
	private final static byte FRONT 	= 4;
    private final static byte BACK 		= 5;
    
    // numerical quantities
	private static final	float   INF=1e15f;
	private static final	float   ZERO=1e-15f;
	private static final	double	PIx2 = Math.PI*2.0;
	private final static double SQR2 = Math.sqrt(2.0);
    private final static double SQR3 = Math.sqrt(3.0);
    private final static double diagdist = 1/(2*SQR2);
    private final static double cubedist = 1/(2*SQR3);
	private final	float	UNKNOWN;
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private 	float[] 		levelset;  			// level set functions
	private 	byte[] 			segmentation;   	// segmentation
	private 	float[][] 		images;  			// original images
	private 	float[] 		wmproba;  			// probability map for inside
	private 	float[] 		gmproba;  			// probability map for outside
	private 	float[] 		csfproba;  		// probability map for other regions
	private		float[]			probafield;		
	//private 	float[] 		bgproba;  		// probability map for other regions
	private		boolean[]		mask;				// masking regions not used in computations (true if computed)
	private		boolean[]		included;			// regions to be included by default (true if included)
	private		boolean[]		critical;
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
		
	// parameters
	private	double		curvweight, divweight, edgeweight, balloonweight;
	private double 		fieldscale;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 5.0f;
	private	float		narrowBandDist = landmineDist+1.8f;
	private	boolean		gwb = false;
	private	boolean		central = true;
	private double		ratio;
	private	double 		mindist, maxdist;
	private	int			npos,nneg,nboth,nnone;
	private	float		offset = 1.0f;
	
	// secondary functions
	private float[]		edgemap;
	private float[][]		edgegrad;
	private	double		Iwm, Igm, Icsf;
	private double		Swm, Sgm, Scsf;	
	private	float[]		Imin, Imax;
	private byte[]		labeling;
	private final static byte REGULAR 	= 0;
    private final static byte FORWARD 	= 1;
	private final static byte BACKWARD	= 2;
    private final static byte BOTH		= 3;
	private final static byte NONE		= 4;
	private float[]		deltamap;
	
	// internal variables for computations
	int	xyz, xyzmx, xyzmy, xyzmz, xyzpx, xyzpy, xyzpz;
	int xyzpxpy, xyzpypz, xyzpzpx, xyzmxmy, xyzmymz, xyzmzmx;
	int xyzmxpy, xyzmypz, xyzmzpx, xyzpxmy, xyzpymz, xyzpzmx;
	int	xyzmx2, xyzmy2, xyzmz2, xyzpx2, xyzpy2, xyzpz2;
	int xyzpxpypz, xyzmxpypz, xyzpxmypz, xyzmxmypz;
	int xyzpxpymz, xyzmxpymz, xyzpxmymz, xyzmxmymz;
	double Dmx, Dmy, Dmz, Dpx, Dpy, Dpz, D0x, D0y, D0z;
	double D0xm, D0ym, D0zm, D0xp, D0yp, D0zp;
	double GPhixp, GPhiyp, GPhizp, GPhixm, GPhiym, GPhizm;
	double SD0x, SD0y, SD0z, GPhi;
	double Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
	double K, G, tmp;
	double DeltaP, DeltaM;
	double delta, forces, balloonforce, pv, div, pfg, pbg;

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

	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=false;
	
	private static class NarrowBand {
		public int[] id;
		public byte[] seg;
		public float[] levels;
		public int currentsize;
		public int capacity;
		public int update;
		
		/** init: create an array of a given size for efficient storage */
		public NarrowBand(int size, int increase) {
			capacity = size;
			update = increase;
			currentsize = 0;
			
			id = new int[capacity];
			seg = new byte[capacity];
			levels = new float[capacity];
		}
		
		public void finalize() {
			capacity = -1;
			id = null;
			seg = null;
			levels = null;
		}
		
		public final void addPoint(int xyz, byte[] lb, float[] lvl) {
			// check for size
			if (currentsize>=capacity-1) {
				capacity += update;
				
				int[] oldid = id;
				byte[] oldseg = seg;
				float[] oldlevels = levels;
				
				id = new int[capacity];
				seg = new byte[capacity];
				levels = new float[capacity];
				
				for (int n=0;n<currentsize;n++) {
					id[n] = oldid[n];
					seg[n] = oldseg[n];
					levels[n] = oldlevels[n];
				}
				
				oldid = null; oldseg = null; oldlevels = null;
			}
			
			// add the new point (use the MGDM variables)
			id[currentsize] = xyz;
			seg[currentsize] = lb[xyz];
			levels[currentsize] = lvl[xyz];
			
			currentsize++;
		}
		
		public void reset() {
			currentsize = 0;
		}
	}


	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public Mp2rageCortexHgdm(byte[] init_, float[][] img_, int nx_, int ny_, int nz_,
							float rx_, float ry_, float rz_,
							float[] wm_, float[] gm_, float[] csf_, float[] bg_,
							boolean[] mask_, boolean[] include_, 
							float ew_, float bw_, float cw_, float dw_,
							String connectivityType_, boolean gwb_, boolean central_) {
		
		images = img_;
		segmentation = init_;
		mask = mask_;
		included = include_;
		wmproba = wm_;
		gmproba = gm_;
		csfproba = csf_;
		//bgproba = bg_;
		
		edgeweight = 0.0f;
		balloonweight = bw_;
		curvweight = cw_;
		divweight = 0.0f;
		fieldscale = 1.0f;
		ratio = dw_;
		offset = ew_;
		
		gwb = gwb_;
		central = central_;
		
		if (debug) System.out.print("GDM forces: "+ew_+" (edges), "+bw_+" (balloon), "+cw_+" (curvature), "+dw_+" (divergence)\n");
		if (verbose) System.out.print("GDM forces: "+ew_+" (edges), "+bw_+" (balloon), "+cw_+" (curvature), "+dw_+" (divergence)\n");
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		UNKNOWN = Numerics.max(nx,ny,nz);
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
				
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		// init all the arrays
		try {
			levelset = new float[nx*ny*nz];
			probafield = new float[nx*ny*nz];
			critical = new boolean[nx*ny*nz];
			labeling = new byte[nx*ny*nz];
			deltamap = new float[nx*ny*nz];
			
			if (mask==null) {
				mask = new boolean[nx*ny*nz];
				for (int xyz=0; xyz<nx*ny*nz; xyz++) {
					mask[xyz] = true;
				}
			}
			
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT("critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT("critical66LUT.raw.gz",200);
				checkComposed=true;
			}
			else if (connectivityType_.equals("wco")) {
				lut = null;
				checkTopology=false;
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
		if (debug) BasicInfo.displayMessage("initial GDM decomposition\n");		
		
		// update mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[x+nx*y+nx*ny*z] = false;
		}
		// init decomposition
		fastMarchingInitializationFromSegmentation(false);
				
		// init proba gdm
		fastMarchingProbabilityInitialization(false);
		
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		levelset = null;
		segmentation = null;
		images = null;
		mask = null;
		included = null;
		heap.finalize();
		heap = null;
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
		
		heap.finalize();
		heap = null;
		System.gc();
	}

	public final float[] getLevelSet() { return levelset; }
	
	public final float[] getProbaField() { return probafield; }
	
	public final byte[] getSegmentation() { return segmentation; }
    
	public final float[] exportSegmentation() {
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) res[xyz] = segmentation[xyz];
		return res;
	}
	
	public final float[] generateProbaSeg() {
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			if (probafield[xyz]<0) 	res[xyz] = OBJ;
			else res[xyz] = BG;
		}
		return res;
	}
	
	public final float[] generateCriticalMap() {
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			if (critical[xyz]) 	res[xyz] = OBJ;
			else res[xyz] = BG;
		}
		return res;
	}
	
	public final float[] generateLabelingMap() {
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			res[xyz] = labeling[xyz];
		}
		return res;
	}
	
	public final float[] getDeltaMap() { return deltamap; }
	
	public final float[] getEdgeMap() { return edgemap; }
	
	public final void fastMarchingInitializationFromSegmentation(boolean narrowBandOnly) {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			levelset[xyz] = UNKNOWN;                            
           
			// segmentation
            if (segmentation[xyz]>0) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
         }
		
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz];
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (mask[xyz]) {
        		processed[xyz] = false;
				// search for boundaries
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (segmentation[xyzn]!=segmentation[xyz]) if (mask[xyzn]) {
						
						// add to the heap
						heap.addValue(0.5f,xyzn,segmentation[xyzn]);
					}
				}
			}
        }
		if (debug) BasicInfo.displayMessage("init\n");		
		
		
        // grow the labels and functions
		float maxdist = 0.0f;
        while ( heap.isNotEmpty() && maxdist<=narrowBandDist+SQR2) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==OBJ) levelset[xyz] = -dist;
			else levelset[xyz] = +dist;
			processed[xyz]=true; // update the current level
 			if (narrowBandOnly) maxdist = dist;	// keep track of distance if stopping at the narrow band
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (segmentation[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							if (processed[xyznb]) if (mask[xyznb]) if (segmentation[xyznb]==lb) {
								nbdist[l] = Numerics.abs(levelset[xyznb]);
								nbflag[l] = true;
							}			
						}
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						// add to the heap
						heap.addValue(newdist,xyzn,lb);
					}
				}
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
        boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
        	processed[xyz] = false;
        	// search for boundaries
        	for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				if (segmentation[xyzn]!=segmentation[xyz]) if (mask[xyzn]) {
					// we assume the levelset value has been recomputed on the boundary
					
					// add to the heap with previous value
					heap.addValue(Numerics.abs(levelset[xyzn]),xyzn,segmentation[xyzn]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float maxdist = 0.0f;
        while (heap.isNotEmpty() && maxdist<=narrowBandDist+SQR2) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==OBJ) levelset[xyz] = -dist;
			else levelset[xyz] = dist;
			processed[xyz]=true; // update the current level
 			if (narrowBandOnly) maxdist = dist;	// keep track of distance if stopping at the narrow band
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (segmentation[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							if (processed[xyznb]) if (mask[xyznb]) if (segmentation[xyznb]==lb) {
								nbdist[l] = Numerics.abs(levelset[xyznb]);
								nbflag[l] = true;
							}			
						}
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						// add to the heap
						heap.addValue(newdist,xyzn,lb);
					}
				}
			}
		}
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }

	public final void fastMarchingProbabilityInitialization(boolean narrowBandOnly) {
		// computation variables
		boolean[] processed = new boolean[nx*ny*nz];
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		byte[] probaseg = new byte[nx*ny*nz];
								
		// initialize the quantities
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			 // gdm functions
			probafield[xyz] = UNKNOWN;
			
			// segmentation
			float phi = csfproba[xyz]-wmproba[xyz];
			if (gwb) phi += gmproba[xyz];
			else phi -= gmproba[xyz];
			if (phi<=0) probaseg[xyz] = OBJ;
			else probaseg[xyz] = BG;
		}
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
		if (debug) BasicInfo.displayMessage("fast marching\n");		
		heap.reset();
		// initialize the heap from boundaries
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (mask[xyz]) {
        		processed[xyz] = false;
        		// search for boundaries
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (probaseg[xyzn]!=probaseg[xyz]) if (mask[xyzn]) {
						float phin = csfproba[xyzn]-wmproba[xyzn];
						if (gwb) phin += gmproba[xyzn];
						else phin -= gmproba[xyzn];
						
						// add to the heap
						if (probaseg[xyzn]==OBJ) heap.addValue(-phin,xyzn,OBJ);
						else heap.addValue(+phin,xyzn,BG);
					}
				}
			}
        }
		if (debug) BasicInfo.displayMessage("init\n");		
		
		
        // grow the labels and functions
		float maxdist = 0.0f;
        while ( heap.isNotEmpty() && maxdist<=narrowBandDist+SQR2) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==OBJ) probafield[xyz] = -dist;
			else probafield[xyz] = +dist;
			processed[xyz]=true; // update the current level
 			if (narrowBandOnly) maxdist = dist;	// keep track of distance if stopping at the narrow band
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (probaseg[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							if (processed[xyznb]) if (mask[xyznb]) if (probaseg[xyznb]==lb) {
								nbdist[l] = Numerics.abs(probafield[xyznb]);
								nbflag[l] = true;
							}			
						}
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						// add to the heap
						heap.addValue(newdist,xyzn,lb);
					}
				}
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
    public static final float minimumMarchingDistance(float[] val, boolean[] flag) {

        float s, s2; // s = a + b +c; s2 = a*a + b*b +c*c
        float tmp;
        int count;
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
        if (count==0) System.err.print("!");
        
        tmp = (s+(float) Math.sqrt((double) (s*s-count*(s2-1.0f))))/count;

        // The larger root
        return tmp;
    }
	/**
     * the isosurface distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * (the input values are all positive, the flags are true only if the isosurface crosses)
     */
    public static final float isoSurfaceDistance(float cur, float[] val, boolean[] flag) {
    	
    	if (cur==0) return 0;
    	
        float s; 
        double dist;
        float tmp;
        s = 0;
        dist = 0;
        
        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.max(val[n], val[n+1]); // Take the largest distance (aka closest to current point) if both are across the boundary
				s = cur/(cur+tmp);
				dist += 1.0/(s*s);
			} else if (flag[n]) {
				s = cur/(cur+val[n]); // Else, take the boundary point
				dist += 1.0/(s*s);
			} else if (flag[n+1]) {
				s = cur/(cur+val[n+1]);
				dist += 1.0/(s*s);
			}
		}
		// triangular (tetrahedral?) relationship of height in right triangles gives correct distance
        tmp = (float)Math.sqrt(1.0/dist);

        // The larger root
        return tmp;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[] generateForces() {
    	float[] forces = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	return forces;
    }
    public final float[] generateCurvatureForce() {
    	float[] forces = new float[nx*ny*nz];
    	
    	double ew = edgeweight;
    	double bw = balloonweight;
    	double cw = curvweight;
    	double dw = divweight;
    	
    	edgeweight = 0.0f; 
    	balloonweight = 0.0f; 
    	curvweight = 1.0f; 
    	divweight = 0.0f;
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	
    	edgeweight = ew;
    	balloonweight = bw;
    	curvweight = cw;
    	divweight = dw;
    	
    	return forces;
    }
    public final float[] generateBalloonForce() {
    	float[] forces = new float[nx*ny*nz];
    	
    	double ew = edgeweight;
    	double bw = balloonweight;
    	double cw = curvweight;
    	double dw = divweight;
    	
    	edgeweight = 0.0f; 
    	balloonweight = 1.0f; 
    	curvweight = 0.0f; 
    	divweight = 0.0f;
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	
    	edgeweight = ew;
    	balloonweight = bw;
    	curvweight = cw;
    	divweight = dw;
    	
    	return forces;
    }
    public final float[] generateDivergenceForce() {
    	float[] forces = new float[nx*ny*nz];
    	
    	double ew = edgeweight;
    	double bw = balloonweight;
    	double cw = curvweight;
    	double dw = divweight;
    	
    	edgeweight = 0.0f; 
    	balloonweight = 0.0f; 
    	curvweight = 0.0f; 
    	divweight = 1.0f;
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	
    	edgeweight = ew;
    	balloonweight = bw;
    	curvweight = cw;
    	divweight = dw;
    	
    	return forces;
    }
    public final float[] generateEdgeForce() {
    	float[] forces = new float[nx*ny*nz];
    	
    	double ew = edgeweight;
    	double bw = balloonweight;
    	double cw = curvweight;
    	double dw = divweight;
    	
    	edgeweight = 1.0f; 
    	balloonweight = 0.0f; 
    	curvweight = 0.0f; 
    	divweight = 0.0f;
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	
    	edgeweight = ew;
    	balloonweight = bw;
    	curvweight = cw;
    	divweight = dw;
    	
    	return forces;
    }
    public final float[] generateNarrowBandLevelset() {
    	float[] band = new float[nx*ny*nz];
    	
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) {
    		if (levelset[xyz]<-narrowBandDist) band[xyz] = -narrowBandDist;
    		else if (levelset[xyz]>narrowBandDist) band[xyz] = narrowBandDist;
    		else band[xyz] = levelset[xyz];
    	}
     	
    	return band;
    }
   
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// init
    	
    	// first estimate the narrow band size
    	int size = 0;
    	int boundarysize=0;
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
			if (Numerics.abs(levelset[xyz])<narrowBandDist && levelset[xyz]!=UNKNOWN) size++;
			if (Numerics.abs(levelset[xyz])<1.0 && levelset[xyz]!=UNKNOWN) boundarysize++;
		}
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
    	mindist = 1e6f;
    	maxdist = 0.0f;
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (Numerics.abs(levelset[xyz])<narrowBandDist) {
				narrowband.addPoint(xyz, segmentation, levelset);
				// in addition, if close to the narrow band boundary, set a landmine
				if (Numerics.abs(levelset[xyz])>=landmineDist) {
					landmines.set(xyz,true);
				}
				// find min, max diff on the narrow band
				if ( Numerics.abs(levelset[xyz]-probafield[xyz]) > maxdist) maxdist = Numerics.abs(levelset[xyz]-probafield[xyz]);
				if ( Numerics.abs(levelset[xyz]-probafield[xyz]) < mindist) mindist = Numerics.abs(levelset[xyz]-probafield[xyz]);
			}
			critical[xyz] = false;
		}
		
		// pre-compute intensity clusters and all
		
			
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && (diff>mindiff || t<5);t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        	if (verbose) System.out.print(t+": ");
        			
			boolean reinit = false;
			int nswap=0;
			
			npos=0; nneg=0; nboth=0; nnone=0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				double force = Numerics.bounded(levelsetForces(xyz), -0.9, 0.9);

				// update the narrow band values, not the original data
				double prev = narrowband.levels[n];
				narrowband.levels[n] -= force;
				
				// reset the critical points all the time?
				//critical[xyz] = false;
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						nswap++;
					
						// update the segmentation
						segmentation[xyz] = OBJ;
						// only when switching labels?
						critical[xyz] = false;
					} else {
						narrowband.levels[n] = lowlevel;
						critical[xyz] = true;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						nswap++;
					
						// update the segmentation
						segmentation[xyz] = BG;
						// only when switching labels?
						critical[xyz] = false;
					} else {
						narrowband.levels[n] = -lowlevel;
						critical[xyz] = true;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}				
			}
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed "+nswap+" ("+(diff*100.0f)+" % of boundary)"
				+", line search (+:"+npos+", -:"+nneg+", x:"+nboth+", o:"+nnone+")\n");
			if (verbose) System.out.print("changed "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
						
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
				
				// recompute min, max
				if ( Numerics.abs(levelset[xyz]-probafield[xyz]) > maxdist) maxdist = Numerics.abs(levelset[xyz]-probafield[xyz]);
				if ( Numerics.abs(levelset[xyz]-probafield[xyz]) < mindist) mindist = Numerics.abs(levelset[xyz]-probafield[xyz]);				
			}
			
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		if (verbose) System.out.print("(*)");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
					// the criterion for being in the narrow band is to have a shortdistance to closest boundaries
					if (Numerics.abs(levelset[xyz])<narrowBandDist && levelset[xyz]!=UNKNOWN) {
						narrowband.addPoint(xyz, segmentation, levelset);
						if (Numerics.abs(levelset[xyz])>=landmineDist) {
							landmines.set(xyz,true);
						}
					}
					if (Numerics.abs(levelset[xyz])<1.0 && levelset[xyz]!=UNKNOWN) {
						boundarysize++;
					}
					critical[xyz] = false;
				}
			}	
     	}
		
		// end of the evolution: recompute the level sets (disable for debugging)
		//resetIsosurfaceNarrowBand(narrowband);
		//fastMarchingReinitialization(false);
		
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
			if (segmentation[xyz+i+j*nx+l*nx*ny]==lb) {
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
			if (segmentation[xyz+i+j*nx+l*nx*ny]==segmentation[xyz]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = false;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		
		// else, it works
		return true;
    }
	
    /** 
    *   isosurface distance re-initialization inside the narrow band 
    *	(assumes a correct narrow band structure enclosing the boundary)
    */
    private final void resetIsosurfaceNarrowBand(NarrowBand narrowband) {
    	if (debug) System.out.print("level set evolution: iso-surface reinit\n");

    	float[] nbdist = new float[6];
    	boolean[] nbflag = new boolean[6];
    	boolean boundary;
    	
		for (int n=0; n<narrowband.currentsize;n++) {
			int xyz = narrowband.id[n];
			
			boundary = false;
			for (int l=0; l<6; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				int xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
				if (segmentation[xyznb]!=segmentation[xyz] && mask[xyznb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(levelset[xyznb]);
					nbflag[l] = true;
					boundary = true;
				}
			}
			if (boundary) {
				narrowband.levels[n] = isoSurfaceDistance(Numerics.abs(levelset[xyz]), nbdist, nbflag);
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int n=0; n<narrowband.currentsize;n++) {
			int xyz = narrowband.id[n];
			levelset[xyz] = narrowband.levels[n];
		}
			
        return;
    }
    
    public final void estimateIntensityClusters() {
    	
    	Imin = new float[1];
    	Imax = new float[1];
    	
    	for (int c=0;c<1;c++) {
    		Imin[c] = ImageStatistics.minimum(images[c],nx,ny,nz);
    		Imax[c] = ImageStatistics.maximum(images[c],nx,ny,nz);
    	}
    	Iwm = 0.0; Igm = 0.0; Icsf = 0.0;
    	Swm = 0.0; Sgm = 0.0; Scsf = 0.0;
    	double Nwm = 0.0;
    	double Ngm = 0.0;
		double Ncsf = 0.0;
		for (int xyz=0; xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
			Iwm += wmproba[xyz]*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
			Nwm += wmproba[xyz]/(nx*ny*nz);
			
			Igm += gmproba[xyz]*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
			Ngm += gmproba[xyz]/(nx*ny*nz);
			
			Icsf += csfproba[xyz]*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
			Ncsf += csfproba[xyz]/(nx*ny*nz);			
		}
		Iwm *= (Imax[0]-Imin[0])/Nwm;
		Igm *= (Imax[0]-Imin[0])/Ngm;
		Icsf *= (Imax[0]-Imin[0])/Ncsf;
    	
		for (int xyz=0; xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
			Swm += wmproba[xyz]*(images[0][xyz]-Iwm)/(Imax[0]-Imin[0])*(images[0][xyz]-Iwm)/(Imax[0]-Imin[0])/(nx*ny*nz);
			Sgm += gmproba[xyz]*(images[0][xyz]-Igm)/(Imax[0]-Imin[0])*(images[0][xyz]-Igm)/(Imax[0]-Imin[0])/(nx*ny*nz);
			Scsf += csfproba[xyz]*(images[0][xyz]-Icsf)/(Imax[0]-Imin[0])*(images[0][xyz]-Icsf)/(Imax[0]-Imin[0])/(nx*ny*nz);
		}
		Swm *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Nwm;
		Sgm *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Ngm;
		Scsf *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Ncsf;

		if (debug) {
			System.out.println("raw values: Iwm  = "+Iwm+", Swm  = "+Swm+", Nwm  = "+Nwm+"\n");
			System.out.println("raw values: Igm  = "+Igm+", Sgm  = "+Sgm+", Ngm  = "+Ngm+"\n");
			System.out.println("raw values: Icsf = "+Icsf+", Scsf = "+Scsf+", Ncsf = "+Ncsf+"\n");
		}			
			
		Swm =  Math.sqrt(Swm);
		Sgm =  Math.sqrt(Sgm);
		Scsf = Math.sqrt(Scsf);

		if (debug) {
			BasicInfo.displayMessage("Intensity clusters \n");
			BasicInfo.displayMessage(" wm:  "+Iwm+" +/- "+Swm+" \n");
			BasicInfo.displayMessage(" gm:  "+Igm+" +/- "+Sgm+" \n");
			BasicInfo.displayMessage(" csf: "+Icsf+" +/- "+Scsf+" \n");
		}
    }
    
    public final void robustIntensityClusters(int itermax, float diffmax) {
    	int t=0;
    	double diff = Imax[0]-Imin[0];
    	
    	while (t<itermax && diff>diffmax) {
    		t++;
    		
			double Nwm = 0.0;
			double Ngm = 0.0;
			double Ncsf = 0.0;
			double Irwm = 0.0;
			double Irgm = 0.0;
			double Ircsf = 0.0;
			double Srwm = 0.0;
			double Srgm = 0.0;
			double Srcsf = 0.0;
			double w;
			for (int xyz=0; xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
				w = wmproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Iwm)/Swm) );
				
				Irwm += w*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
				Nwm += w/(nx*ny*nz);

				w = gmproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Igm)/Sgm) );
				
				Irgm += w*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
				Ngm += w/(nx*ny*nz);

				w = csfproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Icsf)/Scsf) );
				
				Ircsf += w*images[0][xyz]/(Imax[0]-Imin[0])/(nx*ny*nz);
				Ncsf += w/(nx*ny*nz);

			}
			Irwm *= (Imax[0]-Imin[0])/Nwm;
			Irgm *= (Imax[0]-Imin[0])/Ngm;
			Ircsf *= (Imax[0]-Imin[0])/Ncsf;
			
			for (int xyz=0; xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
				w = wmproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Iwm)/Swm) );
				
				Srwm += w*(images[0][xyz]-Irwm)/(Imax[0]-Imin[0])*(images[0][xyz]-Irwm)/(Imax[0]-Imin[0])/(nx*ny*nz);
				
				w = gmproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Igm)/Sgm) );
				
				Srgm += w*(images[0][xyz]-Irgm)/(Imax[0]-Imin[0])*(images[0][xyz]-Irgm)/(Imax[0]-Imin[0])/(nx*ny*nz);
				
				w = csfproba[xyz]/(1.0+Numerics.square( (images[0][xyz]-Icsf)/Scsf) );
				
				Srcsf += w*(images[0][xyz]-Ircsf)/(Imax[0]-Imin[0])*(images[0][xyz]-Ircsf)/(Imax[0]-Imin[0])/(nx*ny*nz);
			}
			Srwm *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Nwm;
			Srgm *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Ngm;
			Srcsf *= (Imax[0]-Imin[0])*(Imax[0]-Imin[0])/Ncsf;
	
			Srwm =  Math.sqrt(Srwm);
			Srgm =  Math.sqrt(Srgm);
			Srcsf = Math.sqrt(Srcsf);
		
			diff = Numerics.max( Numerics.abs(Iwm-Irwm)/(Imax[0]-Imin[0]), 
								  Numerics.abs(Igm-Irgm)/(Imax[0]-Imin[0]), 
								  Numerics.abs(Icsf-Ircsf)/(Imax[0]-Imin[0]),
								  Numerics.abs(Swm-Srwm)/(Imax[0]-Imin[0]),
								  Numerics.abs(Sgm-Srgm)/(Imax[0]-Imin[0]),
								  Numerics.abs(Scsf-Srcsf)/(Imax[0]-Imin[0]) );
								  
			Iwm = Irwm;
			Igm = Irgm;
			Icsf = Ircsf;
			Swm = Srwm;
			Sgm = Srgm;
			Scsf = Srcsf;
			
			if (debug) {
				BasicInfo.displayMessage("Intensity clusters (iter "+t+"\n");
				BasicInfo.displayMessage(" wm:  "+Iwm+" +/- "+Swm+" \n");
				BasicInfo.displayMessage(" gm:  "+Igm+" +/- "+Sgm+" \n");
				BasicInfo.displayMessage(" csf: "+Icsf+" +/- "+Scsf+" \n");
			}
		}
    }
    
    public final void adjustIntensityBoundary(boolean gwb, boolean t1) {
    	if (gwb) {
    		if (t1) Icsf -= 2*Scsf;
    		else Icsf += 2*Scsf;
    	} else {
    		if (t1) Iwm += 2*Swm;
    		else	Iwm -= 2*Swm;
    	}
    }
    
    /** precompute the edge map and its gradient based on a smoothed version of the image
     *	(smoothing scale is in mm) */
    public final void computeEdgeMap(float smoothingScale, float gradientScale) {
    	float[][] kernel = ImageFilters.separableGaussianKernel(smoothingScale/rx,smoothingScale/ry,smoothingScale/rz);
    	float[] smoothed = ImageFilters.separableConvolution(images[0],nx,ny,nz,kernel);
    	
    	edgemap = new float[nx*ny*nz];
    	float dsx, dsy, dsz;
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) {
    		if (mask[xyz]) {
				dsx = 0.5f*(smoothed[xyz+1]-smoothed[xyz-1])/gradientScale;
				dsy = 0.5f*(smoothed[xyz+nx]-smoothed[xyz-nx])/gradientScale;
				dsz = 0.5f*(smoothed[xyz+nx*ny]-smoothed[xyz-nx*ny])/gradientScale;
				
				edgemap[xyz] = 1.0f/(1.0f + dsx*dsx + dsy*dsy + dsz*dsz);
			} else {
				edgemap[xyz] = 1.0f;
			}
		}
    	edgegrad = new float[3][nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		edgegrad[X][xyz] = 0.5f*(edgemap[xyz+1]-edgemap[xyz-1]);
    		edgegrad[Y][xyz] = 0.5f*(edgemap[xyz+nx]-edgemap[xyz-nx]);
    		edgegrad[Z][xyz] = 0.5f*(edgemap[xyz+nx*ny]-edgemap[xyz-nx*ny]);
    	}
    }
    
	/** specific forces applied to the level sets 
	 	(!! we assume forces > 0 expand the boundary, i.e. dphi/dt + f = 0, phi negative inside the object) */
    private final double levelsetForces(int xyz) {
    	forces = 0.0f;
    	
    	// get indices (adapts to the mask)
    	xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
    	xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
    	xyzpz = xyz+nx*ny;	if (!mask[xyzpz]) xyzpz = xyz;
    	xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
    	xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
    	xyzmz = xyz-nx*ny;	if (!mask[xyzmz]) xyzmz = xyz;
    	
    	xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
    	xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
    	xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
    	xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
    	
    	xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
    	xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
    	xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
    	xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
    	
    	xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
    	xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
    	xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
    	xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;    	    	
    	
    	xyzpxpypz = xyz+1+nx+nx*ny;		if (!mask[xyzpxpypz]) xyzpxpypz = xyz;
    	xyzmxpypz = xyz-1+nx+nx*ny;		if (!mask[xyzmxpypz]) xyzmxpypz = xyz;
    	xyzpxmypz = xyz+1-nx+nx*ny;		if (!mask[xyzpxmypz]) xyzpxmypz = xyz;
    	xyzmxmypz = xyz-1-nx+nx*ny;		if (!mask[xyzmxmypz]) xyzmxmypz = xyz;
    	xyzpxpymz = xyz+1+nx-nx*ny;		if (!mask[xyzpxpymz]) xyzpxpymz = xyz;
    	xyzmxpymz = xyz-1+nx-nx*ny;		if (!mask[xyzmxpymz]) xyzmxpymz = xyz;
    	xyzpxmymz = xyz+1-nx-nx*ny;		if (!mask[xyzpxmymz]) xyzpxmymz = xyz;
    	xyzmxmymz = xyz-1-nx-nx*ny;		if (!mask[xyzmxmymz]) xyzmxmymz = xyz;
    	
    	// first derivatives
    	Dmx = levelset[xyz] - levelset[xyzmx];
    	Dmy = levelset[xyz] - levelset[xyzmy];
    	Dmz = levelset[xyz] - levelset[xyzmz];
    	
    	Dpx = levelset[xyzpx] - levelset[xyz];
    	Dpy = levelset[xyzpy] - levelset[xyz];
    	Dpz = levelset[xyzpz] - levelset[xyz];
    	
    	D0x = (levelset[xyzpx] - levelset[xyzmx])/2.0;
    	D0y = (levelset[xyzpy] - levelset[xyzmy])/2.0;
    	D0z = (levelset[xyzpz] - levelset[xyzmz])/2.0;
    	
    	// second derivatives
		Dxx = levelset[xyzmx] + levelset[xyzpx] - 2.0*levelset[xyz];
		Dyy = levelset[xyzmy] + levelset[xyzpy] - 2.0*levelset[xyz];
		Dzz = levelset[xyzmz] + levelset[xyzpz] - 2.0*levelset[xyz];
		
		Dxy = (levelset[xyzmxmy] + levelset[xyzpxpy] - levelset[xyzmxpy] - levelset[xyzpxmy])/4.0;
		Dyz = (levelset[xyzmymz] + levelset[xyzpypz] - levelset[xyzmypz] - levelset[xyzpymz])/4.0;
		Dzx = (levelset[xyzmzmx] + levelset[xyzpzpx] - levelset[xyzmzpx] - levelset[xyzpzmx])/4.0;
			
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
		//if (K>1.0 || K<-1.0) K=0;
			
		// curvature smoothing
		forces += -curvweight*stepsize*K*GPhi;

		// divergence regularizing force: problem = requires a second order neighborhood :(
		xyzpx2 = xyz+2;			if (!mask[xyzpx2]) xyzpx2 = xyzpx;
		xyzpy2 = xyz+2*nx;		if (!mask[xyzpy2]) xyzpy2 = xyzpy;
		xyzpz2 = xyz+2*nx*ny;	if (!mask[xyzpz2]) xyzpz2 = xyzpz;
		
		xyzmx2 = xyz-2;			if (!mask[xyzmx2]) xyzmx2 = xyzmx;
		xyzmy2 = xyz-2*nx;		if (!mask[xyzmy2]) xyzmy2 = xyzmy;
		xyzmz2 = xyz-2*nx*ny;	if (!mask[xyzmz2]) xyzmz2 = xyzmz;
		
		D0xp = 0.5*(levelset[xyzpx2]-levelset[xyz]);
		D0yp = 0.5*(levelset[xyzpy2]-levelset[xyz]);
		D0zp = 0.5*(levelset[xyzpz2]-levelset[xyz]);
		
		D0xm = 0.5*(levelset[xyz]-levelset[xyzmx2]);
		D0ym = 0.5*(levelset[xyz]-levelset[xyzmy2]);
		D0zm = 0.5*(levelset[xyz]-levelset[xyzmz2]);
		
		GPhixp = Math.sqrt(D0xp*D0xp + Numerics.square(0.5*(levelset[xyzpxpy]-levelset[xyzpxmy])) + Numerics.square(0.5*(levelset[xyzpzpx]-levelset[xyzmzpx])));
		GPhiyp = Math.sqrt(Numerics.square(0.5*(levelset[xyzpxpy]-levelset[xyzmxpy])) + D0yp*D0yp + Numerics.square(0.5*(levelset[xyzpypz]-levelset[xyzpymz])));
		GPhizp = Math.sqrt(Numerics.square(0.5*(levelset[xyzpzpx]-levelset[xyzpzmx])) + Numerics.square(0.5*(levelset[xyzpypz]-levelset[xyzmypz])) + D0zp*D0zp);
		
		GPhixm = Math.sqrt(D0xm*D0xm + Numerics.square(0.5*(levelset[xyzmxpy]-levelset[xyzmxmy])) + Numerics.square(0.5*(levelset[xyzpzmx]-levelset[xyzmzmx])));
		GPhiym = Math.sqrt(Numerics.square(0.5*(levelset[xyzpxmy]-levelset[xyzmxmy])) + D0ym*D0ym + Numerics.square(0.5*(levelset[xyzmypz]-levelset[xyzmymz])));
		GPhizm = Math.sqrt(Numerics.square(0.5*(levelset[xyzmzpx]-levelset[xyzmzmx])) + Numerics.square(0.5*(levelset[xyzpymz]-levelset[xyzmymz])) + D0zm*D0zm);
		
		GPhixm = divergenceMeasure(GPhixm);
		GPhiym = divergenceMeasure(GPhiym);
		GPhizm = divergenceMeasure(GPhizm);
		GPhixp = divergenceMeasure(GPhixp);
		GPhiyp = divergenceMeasure(GPhiyp);
		GPhizp = divergenceMeasure(GPhizp);
		
		forces += -divweight*stepsize*( 0.5*(GPhixp*D0xp - GPhixm*D0xm) + 0.5*(GPhiyp*D0yp - GPhiym*D0ym) + 0.5*(GPhizp*D0zp - GPhizm*D0zm) );		
		/*
		// new force paradigm : check for cases
		double phi = csfproba[xyz]-wmproba[xyz];
		if (gwb) phi += gmproba[xyz];
		else phi -= gmproba[xyz];
		
		//balloonforce = Numerics.abs(phi);
		*/
		
		// almost working, but not for critical points
		// => try fully explicit search
		/*
		//balloonforce = Numerics.abs(probafield[xyz] - levelset[xyz])/maxdist;
		//balloonforce = 0.5 + 0.5*(1.0+Numerics.abs(probafield[xyz] - levelset[xyz])-mindist)/(1.0+maxdist-mindist);
		balloonforce = 1.0;
		
		//double sgn = -Numerics.sign(probafield[xyz] - levelset[xyz]);
		double Dfx, Dfy,Dfz;
		double prod;
		Dfx = - (probafield[xyzpx] - probafield[xyzmx])/2.0;
		Dfy = - (probafield[xyzpy] - probafield[xyzmy])/2.0;
		Dfz = - (probafield[xyzpz] - probafield[xyzmz])/2.0;
		prod = Dfx*D0x+Dfy*D0y+Dfz*D0z;
		
		if (probafield[xyz]>0) {
			// outside target : do nothing
		} else {
			// inside target
			Dfx = D0x;
			Dfy = D0y;
			Dfz = D0z;
		}
		*/
		
		balloonforce = 1.0;
		//balloonforce = Numerics.min(Numerics.abs(probafield[xyz] - levelset[xyz])/offset, 1.0f);
		//balloonforce = Numerics.bounded( (probafield[xyz] - levelset[xyz])/offset, -1.0f, 1.0f);
		
		double sign = 0.0;
		double Dfx, Dfy, Dfz;
		
		if (probafield[xyz]<=0) {
			Dfx = - (probafield[xyzpx] - probafield[xyzmx])/2.0;
			Dfy = - (probafield[xyzpy] - probafield[xyzmy])/2.0;
			Dfz = - (probafield[xyzpz] - probafield[xyzmz])/2.0;
			sign = 1.0;
		} else {			
			Dfx = - (probafield[xyzpx] - probafield[xyzmx])/2.0;
			Dfy = - (probafield[xyzpy] - probafield[xyzmy])/2.0;
			Dfz = - (probafield[xyzpz] - probafield[xyzmz])/2.0;
			sign = Dfx*D0x+Dfy*D0y+Dfz*D0z;			
		}

		// check if neighbor of critical point?
		boolean nonCritical = true;
		labeling[xyz] = REGULAR;
		deltamap[xyz] = 0.0f;
		
		if (critical[xyz] || 
			critical[xyzpx] || critical[xyzpy] || critical[xyzpz] || critical[xyzmx] || critical[xyzmy] || critical[xyzmz] || 
			critical[xyzpxpy] || critical[xyzpxmy] || critical[xyzmxpy] || critical[xyzmxmy] || 
			critical[xyzpypz] || critical[xyzpymz] || critical[xyzmypz] || critical[xyzmymz] || 
			critical[xyzpzpx] || critical[xyzpzmx] || critical[xyzmzpx] || critical[xyzmzmx] ||
			critical[xyzpxpypz] || critical[xyzmxpypz] || critical[xyzpxmypz] || critical[xyzmxmypz] ||
			critical[xyzpxpymz] || critical[xyzmxpymz] || critical[xyzpxmymz] || critical[xyzmxmymz] ) {
			nonCritical = false;
		}		
		
		if (!nonCritical) {
			// fully explicit search: move toward the closest object in the normal direction D0
			boolean found=false;
			if (GPhi==0) {
				found = true;
				//sign = +1.0f*Numerics.sign(probafield[xyz]); 
			}
			float delta = 1.0f;
			float highp = -1.0f;
			float highm = -1.0f;
			int d;
			int xyzp, xyzm;
			int id;
			int i,j,l;
			
			double xp,yp,zp;
			double xm,ym,zm;
			double vpx = D0x/GPhi;
			double vpy = D0y/GPhi;
			double vpz = D0z/GPhi;
			double vmx = -D0x/GPhi;
			double vmy = -D0y/GPhi;
			double vmz = -D0z/GPhi;
			double nvp, nvm;
			xyzp = xyz;
			xyzm = xyz;
			
			zp = Numerics.floor(xyz/(nx*ny));
			yp = Numerics.floor((xyz-nx*ny*zp)/nx);
			xp = Numerics.floor(xyz-nx*ny*zp-nx*yp);
			xm = xp;
			ym = yp;
			zm = zp;
			while (!found) {
				// linear increment (slow, but simple)
				
				/*
				xyzp = xyz + Numerics.round(delta*vpx) 
						    + Numerics.round(delta*vpy)*nx 
						    + Numerics.round(delta*vpz)*nx*ny;	
				xyzm = 2*xyz - xyzp;	
				if (vpx*Dfx+vpy*Dfy+vpz*Dfz>0) {
					vpx = (1.0-ratio)*vpx + ratio*Dfx;
					vpy = (1.0-ratio)*vpy + ratio*Dfy;
					vpz = (1.0-ratio)*vpz + ratio*Dfz;
				} else {
					vpx = (1.0-ratio)*vpx + ratio*Dfx;
					vpy = (1.0-ratio)*vpy + ratio*Dfy;
					vpz = (1.0-ratio)*vpz + ratio*Dfz;
				}
				nvp = Math.sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
				vpx /= nvp;
				vpy /= nvp;
				vpz /= nvp;
				*/
				
				xp += vpx;
				yp += vpy;
				zp += vpz;
				xyzp = Numerics.round(xp) + Numerics.round(yp)*nx + Numerics.round(zp)*nx*ny;
				xm += vmx;
				ym += vmy;
				zm += vmz;
				xyzm = Numerics.round(xm) + Numerics.round(ym)*nx + Numerics.round(zm)*nx*ny;
				
				if (xyzp>nx*ny && xyzp<nx*ny*(nz-1)) {
					// slowly change the vector to align with target directions : two different ones?
					vpx = (1.0-ratio)*vpx - ratio*0.5*(probafield[xyzp+1]-probafield[xyzp-1]);
					vpy = (1.0-ratio)*vpy - ratio*0.5*(probafield[xyzp+nx]-probafield[xyzp-nx]);
					vpz = (1.0-ratio)*vpz - ratio*0.5*(probafield[xyzp+nx*ny]-probafield[xyzp-nx*ny]);
					nvp = Math.sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
					vpx /= nvp;
					vpy /= nvp;
					vpz /= nvp;
				}
				
				if (xyzm>nx*ny && xyzm<nx*ny*(nz-1)) {
					vmx = (1.0-ratio)*vmx - ratio*0.5*(probafield[xyzm+1]-probafield[xyzm-1]);
					vmy = (1.0-ratio)*vmy - ratio*0.5*(probafield[xyzm+nx]-probafield[xyzm-nx]);
					vmz = (1.0-ratio)*vmz - ratio*0.5*(probafield[xyzm+nx*ny]-probafield[xyzm-nx*ny]);
					nvm = Math.sqrt(vmx*vmx+vmy*vmy+vmz*vmz);
					vmx /= nvm;
					vmy /= nvm;
					vmz /= nvm;
				}
				
				if (xyzp>=0 && xyzp<nx*ny*nz) if (probafield[xyz]*probafield[xyzp]<0) {
					highp = delta;
				}
				if (xyzm>=0 && xyzm<nx*ny*nz) if (probafield[xyz]*probafield[xyzm]<0) {
					highm = delta;
				}
				if (highp>-1 && highm==-1) { 
					sign = +1.0f; 
					//sign = +Numerics.sign(probafield[xyz]); 
					found = true; 
					npos++;
					labeling[xyz] = FORWARD;
				} else if (highm>-1 && highp==-1) { 
					sign = -1.0f; 
					//sign = -Numerics.sign(probafield[xyz]); 
					found = true; 
					nneg++;
					labeling[xyz] = BACKWARD;
				} else if (highp>-1 && highm>-1) { 
					found = true;
					//sign = -0.5f*Numerics.sign(probafield[xyz]);
					//sign = +0.5f*Numerics.sign(probafield[xyz]); 
					nboth++;
					labeling[xyz] = BOTH;
				} else if ( ( (xyzp<=0 || xyzp>=nx*ny*nz) && (xyzm<=0 || xyzm>=nx*ny*nz) ) 
							||  delta>5*nx || delta>5*ny || delta>5*nz) {	
					found  = true;
					//sign = -0.5f*Numerics.sign(probafield[xyz]);
					//sign = +0.5f*Numerics.sign(probafield[xyz]); 
					nnone++;
					labeling[xyz] = NONE;
				} else delta++;			
			}
			deltamap[xyz] = delta;
		}
		Dfx = sign*D0x;
		Dfy = sign*D0y;
		Dfz = sign*D0z;
		
		/*
		double alpha = 1.0;
		double beta = 0.0;
		
		if (levelset[xyz]<=0 && probafield[xyz]>0) {
			// inside levelset but outside target
			if (prod < 0) {
				Dfx = alpha*Dfx - beta*D0x;
				Dfy = alpha*Dfy - beta*D0y;
				Dfz = alpha*Dfz - beta*D0z;
			} else {
				Dfx = alpha*Dfx + beta*D0x;
				Dfy = alpha*Dfy + beta*D0y;
				Dfz = alpha*Dfz + beta*D0z;
			}
		} else if (levelset[xyz]>0 && probafield[xyz]<=0) {
			// inside target but outside levelset
			if (prod < 0) {
				Dfx = - beta*Dfx + alpha*D0x;
				Dfy = - beta*Dfy + alpha*D0y;
				Dfz = - beta*Dfz + alpha*D0z;					
			} else {
				Dfx = beta*Dfx + alpha*D0x;
				Dfy = beta*Dfy + alpha*D0y;
				Dfz = beta*Dfz + alpha*D0z;
			}
		} else if (levelset[xyz]<=0 && probafield[xyz]<=0) {
			// inside both
			if (prod < 0) {
				Dfx = - beta*Dfx + alpha*D0x;
				Dfy = - beta*Dfy + alpha*D0y;
				Dfz = - beta*Dfz + alpha*D0z;					
			} else {
				Dfx = beta*Dfx + alpha*D0x;
				Dfy = beta*Dfy + alpha*D0y;
				Dfz = beta*Dfz + alpha*D0z;
			}
		} else {
			// outside both
			if (prod < 0) {
				Dfx = alpha*Dfx - beta*D0x;
				Dfy = alpha*Dfy - beta*D0y;
				Dfz = alpha*Dfz - beta*D0z;					
			} else {
				Dfx = alpha*Dfx + beta*D0x;
				Dfy = alpha*Dfy + beta*D0y;
				Dfz = alpha*Dfz + beta*D0z;
			}			
		}
		*/
		forces += balloonweight*stepsize*balloonforce*(Numerics.max(Dfx,0)*Dmx + Numerics.min(Dfx,0)*Dpx 
														+Numerics.max(Dfy,0)*Dmy + Numerics.min(Dfy,0)*Dpy 
														+Numerics.max(Dfz,0)*Dmz + Numerics.min(Dfz,0)*Dpz);

		/*
		// add some expanding ballon force to counteract shrinking?
		forces += balloonweight*stepsize*ratio*(Numerics.max(D0x,0)*Dmx + Numerics.min(D0x,0)*Dpx 
												+Numerics.max(D0y,0)*Dmy + Numerics.min(D0y,0)*Dpy 
												+Numerics.max(D0z,0)*Dmz + Numerics.min(D0z,0)*Dpz);
		*/
					
		/*
		double Dm2 =	(Numerics.max(Dfx,0)*Dmx + Numerics.min(Dfx,0)*Dpx 
						+Numerics.max(Dfy,0)*Dmy + Numerics.min(Dfy,0)*Dpy 
						+Numerics.max(Dfz,0)*Dmz + Numerics.min(Dfz,0)*Dpz);
		
		double Dp2 =	(Numerics.max(Dfx,0)*Dpx + Numerics.min(Dfx,0)*Dmx 
						+Numerics.max(Dfy,0)*Dpy + Numerics.min(Dfy,0)*Dmy 
						+Numerics.max(Dfz,0)*Dpz + Numerics.min(Dfz,0)*Dmz);
		
		forces += balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*Dm2 + Numerics.min(balloonforce, 0.0)*Dp2);
		*/
		/*		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		balloonforce = -Numerics.bounded( (probafield[xyz] - levelset[xyz])/fieldscale, -1.0, 1.0);
		
		boolean variant = ( critical[xyzpx] || critical[xyzmx] || critical[xyzpy] || critical[xyzmy] || critical[xyzpz] || critical[xyzmz] 
							|| critical[xyzmxmy] || critical[xyzmxpy] || critical);
		
		if (!variant) {
			if (central) { 
				DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
										 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
										 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
				
				DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
										 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
										 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
	
				forces += balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*DeltaP + Numerics.min(balloonforce, 0.0)*DeltaM);
			} else {
				double Dm2 =	(Numerics.max(D0x,0)*Dmx + Numerics.min(D0x,0)*Dpx 
								+Numerics.max(D0y,0)*Dmy + Numerics.min(D0y,0)*Dpy 
								+Numerics.max(D0z,0)*Dmz + Numerics.min(D0z,0)*Dpz);
				
				double Dp2 =	(Numerics.max(D0x,0)*Dpx + Numerics.min(D0x,0)*Dmx 
								+Numerics.max(D0y,0)*Dpy + Numerics.min(D0y,0)*Dmy 
								+Numerics.max(D0z,0)*Dpz + Numerics.min(D0z,0)*Dmz);
				
				forces += balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*Dm2 + Numerics.min(balloonforce, 0.0)*Dp2);
			}
		} else {
			// field direction
			double Dfx = (probafield[xyzpx] - probafield[xyzmx])/2.0;
			double Dfy = (probafield[xyzpy] - probafield[xyzmy])/2.0;
			double Dfz = (probafield[xyzpz] - probafield[xyzmz])/2.0;
			
			if (central) {
				forces += balloonweight*stepsize*balloonforce*(Dfx*D0x + Dfy*D0y + Dfz*D0z);
			} else {
				double Dm2 =	(Numerics.max(Dfx,0)*Dmx + Numerics.min(Dfx,0)*Dpx 
								+Numerics.max(Dfy,0)*Dmy + Numerics.min(Dfy,0)*Dpy 
								+Numerics.max(Dfz,0)*Dmz + Numerics.min(Dfz,0)*Dpz);
				
				double Dp2 =	(Numerics.max(Dfx,0)*Dpx + Numerics.min(Dfx,0)*Dmx 
								+Numerics.max(Dfy,0)*Dpy + Numerics.min(Dfy,0)*Dmy 
								+Numerics.max(Dfz,0)*Dpz + Numerics.min(Dfz,0)*Dmz);
				
				forces += balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*Dm2 + Numerics.min(balloonforce, 0.0)*Dp2);
			}
		}
		*/
		return forces;
    }
    
    private final double divergenceMeasure(double grad) {
    	if(grad < 0.0000001) div = 1.0;
		else if (grad<1.0) div = Math.sin(PIx2*grad)/(PIx2*grad);
		else div = (grad - 1.0)/grad;
		return div;
	}

}

