package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

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
 *	@version    Jul 2010
 *	@author     Pierre-Louis Bazin 
 *	@author     John Bogovic
 * 	@author 	Hanlin Wan
 *		
 *
 */
 
public class Gdm3d {
	
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
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private final static float SQR2 = (float) Math.sqrt(2.0f);
    private final static float SQR3 = (float) Math.sqrt(3.0f);
    private final static float diagdist = 1/(2*SQR2);
    private final static float cubedist = 1/(2*SQR3);
	private final	float	UNKNOWN;
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private 	float[] 		levelset;  			// level set functions
	private 	byte[] 			segmentation;   	// segmentation
	private 	float[][] 		fieldforce;  		// original image forces, indep. object (e.g. boundaries)
	private 	float[] 		balloonforce;  		// original image forces, along object normals (e.g. from memberships)
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// parameters
	private	double		smoothweight, fieldweight, balloonweight;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private	float		narrowBandDist = 6.5f;
	private float		landmineDist = 5.0f;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
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
	public Gdm3d(int[] init_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_) {
		this(init_, nx_, ny_, nz_, rx_, ry_, rz_, field_, balloon_, fw_, bw_, sw_, connectivityType_, null);
	}
	public Gdm3d(int[] init_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_, String connectivityPath_) {
		
		fieldforce = field_;
		balloonforce = balloon_;
		
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		
		if (debug) System.out.print("GDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
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
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT(connectivityPath_, "critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT(connectivityPath_, "critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT(connectivityPath_, "critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT(connectivityPath_, "critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT(connectivityPath_, "critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT(connectivityPath_,"critical66LUT.raw.gz",200);
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
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
		}
		// init decomposition
		fastMarchingInitializationFromSegmentation(init_,false);
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		

	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public Gdm3d(float[] lvl_, float maxdist_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_) {
		this(lvl_, maxdist_, nx_, ny_, nz_, rx_, ry_, rz_, field_, balloon_, fw_, bw_, sw_,	connectivityType_, null);
	}	
	public Gdm3d(float[] lvl_, float maxdist_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_, String connectivityPath_) {
		
		levelset = lvl_;
	
		fieldforce = field_;
		balloonforce = balloon_;
		
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		
		if (debug) System.out.print("GDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
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
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT(connectivityPath_, "critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT(connectivityPath_, "critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT(connectivityPath_, "critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT(connectivityPath_, "critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT(connectivityPath_, "critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT(connectivityPath_,"criticalWCLUT.raw.gz",200);
				checkComposed=false;
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
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
		}
		// init decomposition
		initializationFromLevelset(true, maxdist_);
			
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		

	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public Gdm3d(boolean[] init_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_) {
		this(init_, nx_, ny_, nz_, rx_, ry_, rz_, field_, balloon_, fw_, bw_, sw_, connectivityType_, null);
	}
	public Gdm3d(boolean[] init_, int nx_, int ny_, int nz_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[] balloon_, 
						float fw_, float bw_, float sw_,
						String connectivityType_, String connectivityPath_) {
		
		fieldforce = field_;
		balloonforce = balloon_;
		
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		
		if (debug) System.out.print("GDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
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
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT(connectivityPath_, "critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT(connectivityPath_, "critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT(connectivityPath_, "critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT(connectivityPath_, "critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT(connectivityPath_, "critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT(connectivityPath_,"criticalWCLUT.raw.gz",200);
				checkComposed=false;
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
					//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial GDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
		}
		// init decomposition
		fastMarchingInitializationFromSegmentation(init_,true);
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		levelset = null;
		segmentation = null;
		fieldforce = null;
		balloonforce = null;
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
	
	public final byte[] getSegmentation() { return segmentation; }
    
	public final float[] exportSegmentation() {
    	float[] seg = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		seg[xyz] = (float)segmentation[xyz];
    	}
    	return seg;
    }
	public final float[] exportLevelset() {
    	float[] lvl = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) {
    		lvl[xyz] = (float)levelset[xyz];
    	}
    	return lvl;
    }
	public final void initializationFromLevelset(boolean restrictedRegion, float dist) {
		// initialize the quantities
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			// segmentation
			if (levelset[xyz]<=0) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
			// mask
			if (restrictedRegion && Numerics.abs(levelset[xyz])>dist) mask[xyz] = false;
		}
		return;
	}
     
	public final void fastMarchingInitializationFromSegmentation(int[] init, boolean narrowBandOnly) {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			levelset[xyz] = UNKNOWN;                            
           
			// segmentation
            if (init[xyz]>0) {
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
     
	public final void fastMarchingInitializationFromSegmentation(boolean[] init, boolean narrowBandOnly) {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			levelset[xyz] = UNKNOWN;                            
           
			// segmentation
            if (init[xyz]) {
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
    public final float[] levelsetForces() {
    	float[] forces = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	return forces;
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
        
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (Numerics.abs(levelset[xyz])<narrowBandDist) {
				narrowband.addPoint(xyz, segmentation, levelset);
				// in addition, if close to the narrow band boundary, set a landmine
				if (Numerics.abs(levelset[xyz])>=landmineDist) {
					landmines.set(xyz,true);
				}
			}
		}
			
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && (diff>mindiff || t<5);t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinit = false;
			int nswap=0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				double force = levelsetForces(xyz);

				// update the narrow band values, not the original data
				double prev = narrowband.levels[n];
				narrowband.levels[n] -= force;
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						nswap++;
					
						// update the segmentation
						segmentation[xyz] = OBJ;
					} else {
						narrowband.levels[n] = lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						nswap++;
					
						// do nothing
						segmentation[xyz] = BG;
					} else {
						narrowband.levels[n] = -lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
			}
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
			}
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines = new BitSet(Numerics.ceil(0.2f*size));
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
				}
			}	
     	}
		
		// end of the evolution: recompute the level sets (disabled for debugging)
		resetIsosurfaceNarrowBand(narrowband);
		fastMarchingReinitialization(false);
		
        return;
    }
    
    
	/** specific forces applied to the level sets (application dependent) */
    private final double levelsetForces(int xyz) {
    	
    	// first derivatives
    	double Dmx = levelset[xyz] - levelset[xyz-1];
    	double Dmy = levelset[xyz] - levelset[xyz-nx];
    	double Dmz = levelset[xyz] - levelset[xyz-nx*ny];
    	
    	double Dpx = levelset[xyz+1] - levelset[xyz];
    	double Dpy = levelset[xyz+nx] - levelset[xyz];
    	double Dpz = levelset[xyz+nx*ny] - levelset[xyz];
    	
    	double D0x = (levelset[xyz+1] - levelset[xyz-1])/2.0;
    	double D0y = (levelset[xyz+nx] - levelset[xyz-nx])/2.0;
    	double D0z = (levelset[xyz+nx*ny] - levelset[xyz-nx*ny])/2.0;
    	
    	// if using smoothing forces:
    	double smooth = 0.0;
    	if (Numerics.abs(smoothweight)>0) {
    	
			// second derivatives
			double Dxx = levelset[xyz-1] + levelset[xyz+1] - 2.0*levelset[xyz];
			double Dyy = levelset[xyz-nx] + levelset[xyz+nx] - 2.0*levelset[xyz];
			double Dzz = levelset[xyz-nx*ny] + levelset[xyz+nx*ny] - 2.0*levelset[xyz];
			
			double Dxy = (levelset[xyz-1-nx] + levelset[xyz+1+nx] - levelset[xyz-1+nx] - levelset[xyz+1-nx])/4.0;
			double Dyz = (levelset[xyz-nx-nx*ny] + levelset[xyz+nx+nx*ny] - levelset[xyz-nx+nx*ny] - levelset[xyz+nx-nx*ny])/4.0;
			double Dzx = (levelset[xyz-nx*ny-1] + levelset[xyz+nx*ny+1] - levelset[xyz-nx*ny+1] - levelset[xyz+nx*ny-1])/4.0;
			
			// gradient norm
			double SD0x = D0x * D0x;
			double SD0y = D0y * D0y;
			double SD0z = D0z * D0z;
			double GPhi = Math.sqrt(SD0x + SD0y + SD0z);
			
			// mean curvature: K = k1 + k2
			double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
						- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
				
			// gaussian curvature: G = k1*k2
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
			
			balloon = balloonweight*stepsize*(Numerics.max(balloonforce[xyz], 0.0)*DeltaP + Numerics.min(balloonforce[xyz], 0.0)*DeltaM);
		}
		
		// external force field
		double field = 0.0;
		if (Numerics.abs(fieldweight)>0) {
			field =  fieldweight*stepsize*(Numerics.max(fieldforce[X][xyz],0.0)*Dmx + Numerics.min(fieldforce[X][xyz],0.0)*Dpx
										  +Numerics.max(fieldforce[Y][xyz],0.0)*Dmy + Numerics.min(fieldforce[Y][xyz],0.0)*Dpy
										  +Numerics.max(fieldforce[Z][xyz],0.0)*Dmz + Numerics.min(fieldforce[Z][xyz],0.0)*Dpz);
		}
		
		return - (smooth - balloon - field);
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

		/*
		// does it change the topology of a relation between the modified object and its neighbors ?
		int  Nconfiguration = 0;
		short[] lbs = new short[26];
		for (short i=-1;i<=1;i++) for (short j=-1;j<=1;j++) for (short l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[xyz+i+j*nx+l*nx*ny]!=lb) 
				&& (segmentation[xyz+i+j*nx+l*nx*ny]!=mgdmlabels[0][xyz]) ) {
				boolean found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[n]) { found = true; break; }
				
				if (!found) {
					lbs[Nconfiguration] = segmentation[xyz+i+j*nx+l*nx*ny];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==mgdmlabels[0][xyz])
					|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[n]) ) {
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
				if ( (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lb)
					|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[n]) ) {
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
					if ( (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==mgdmlabels[0][xyz])
						|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[n])
						|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[m]) ) {
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
					if ( (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lb)
						|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[n]) 
						|| (mgdmlabels[0][xyz+i+j*nx+l*nx*ny]==lbs[m]) ) {
						obj[1+i][1+j][1+l] = true;
					} else {
						obj[1+i][1+j][1+l] = false;
					}
				}
				obj[1][1][1] = true;
				if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			}
		}
		*/
		
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
    

}

