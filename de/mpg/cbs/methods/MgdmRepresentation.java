package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm handles basic manipulations within the MGDM framework	
 *
 *	@version    December 2012
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class MgdmRepresentation {
	
	
	// data and membership buffers
	private 	float[][] 		mgdmfunctions;  	// MGDM's pseudo level set mgdmfunctions
	private 	byte[][] 		mgdmlabels;   		// MGDM's label maps
	private static	byte 	  	nmgdm;					// total number of MGDM mgdmlabels and mgdmfunctions
	private 	int				nx,ny,nz;   		// image dimensions
	private 	float			rx,ry,rz;   		// image resolutions
	
	// atlas parameters
	private 	int 				nobj;    	// number of shapes
	private 	byte[]			objLabel;			// label values in the original image
	    
	// data and membership buffers
	private 	byte[] 			segmentation;   	// MGDM's segmentation
	private		short[]			counter;
	private		boolean[]		mask;				// masking regions not used in computations
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	
	// parameters
	private	float		maxDist = 1e9f;
	private	short		maxcount = 5;
	private	boolean		stopDist = false;
	private	boolean		skipZero = false;
	
	// computation variables to avoid re-allocating
	
	// for minimumMarchingDistance
	double s, s2, tmp; 
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
	 *  constructors 
	 */
	public MgdmRepresentation(byte[] seg_, float[] bound_, 
									int nx_, int ny_, int nz_, float rx_, float ry_, float rz_,
									byte[] labels_, int nlb_,
									int nmgdm_, boolean skip_, float dist_) {
	
		nx = nx_;
		ny = ny_;
		nz = nz_;
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		nobj = nlb_;
		objLabel = labels_;
		
		nmgdm = (byte)nmgdm_;
		skipZero = skip_;
		stopDist = (dist_==-1);
		if (!stopDist) maxDist = dist_;
		
		maxcount = Numerics.max(maxcount, nmgdm);
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};

		// init all the arrays
		try {
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];	
			counter = new short[nx*ny*nz];	
			mgdmfunctions = new float[nmgdm][nx*ny*nz];
			mgdmlabels = new byte[nmgdm+1][nx*ny*nz];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (skipZero && seg_[x+nx*y+nx*ny*z]==0) mask[x+nx*y+nx*ny*z] = false;
		}
		// init segmentation (from atlas if null)
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==seg_[xyz]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
			counter[xyz] = 0;
			
			// init the MGDM model
			mgdmlabels[0][xyz] = nlb;
			if (bound_!=null) mgdmfunctions[0][xyz] = bound_[xyz];
			else mgdmfunctions[0][xyz] = 0.5f;
		}
		// build the full model
		fastMarchingReinitialization(stopDist, false, true);
		
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	/**
	 *  constructors 
	 */
	public MgdmRepresentation(int[] seg_, float[] bound_, 
									int nx_, int ny_, int nz_, float rx_, float ry_, float rz_,
									byte[] labels_, int nlb_,
									int nmgdm_, boolean skip_, float dist_) {
	
		nx = nx_;
		ny = ny_;
		nz = nz_;
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		nobj = nlb_;
		objLabel = labels_;
		
		nmgdm = (byte)nmgdm_;
		skipZero = skip_;
		stopDist = (dist_==-1);
		if (!stopDist) maxDist = dist_;
		
		maxcount = Numerics.max(maxcount, nmgdm);
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};

		// init all the arrays
		try {
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];	
			counter = new short[nx*ny*nz];	
			mgdmfunctions = new float[nmgdm][nx*ny*nz];
			mgdmlabels = new byte[nmgdm+1][nx*ny*nz];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (skipZero && seg_[x+nx*y+nx*ny*z]==0) mask[x+nx*y+nx*ny*z] = false;
		}
		// init segmentation (from atlas if null)
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==seg_[xyz]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
			counter[xyz] = 0;
			
			// init the MGDM model
			mgdmlabels[0][xyz] = nlb;
			if (bound_!=null) mgdmfunctions[0][xyz] = bound_[xyz];
			else mgdmfunctions[0][xyz] = 0.5f;
		}
		// build the full model
		fastMarchingReinitialization(stopDist, false, true);
		
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	/**
	 *  constructors 
	 */
	public MgdmRepresentation(byte[][][] seg_, float[][][] bound_, 
									int nx_, int ny_, int nz_, float rx_, float ry_, float rz_,
									byte[] labels_, int nlb_,
									int nmgdm_, boolean skip_, float dist_) {
	
		nx = nx_;
		ny = ny_;
		nz = nz_;
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		nobj = nlb_;
		objLabel = labels_;
		
		nmgdm = (byte)nmgdm_;
		skipZero = skip_;
		stopDist = (dist_==-1);
		if (!stopDist) maxDist = dist_;
		
		maxcount = Numerics.max(maxcount, nmgdm);
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};

		// init all the arrays
		try {
			segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];	
			counter = new short[nx*ny*nz];	
			mgdmfunctions = new float[nmgdm][nx*ny*nz];
			mgdmlabels = new byte[nmgdm+1][nx*ny*nz];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial MGDM decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (skipZero && seg_[x][y][z]==0) mask[x+nx*y+nx*ny*z] = false;
		}
		// init segmentation (from atlas if null)
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			
			byte nlb = EMPTY;
			for (byte n=0; n<nobj; n++) {
				if (objLabel[n]==seg_[x][y][z]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
			counter[xyz] = 0;
			
			// init the MGDM model
			mgdmlabels[0][xyz] = nlb;
			if (bound_!=null) mgdmfunctions[0][xyz] = bound_[x][y][z];
			else mgdmfunctions[0][xyz] = 0.5f;
		}
		// build the full model
		fastMarchingReinitialization(stopDist, false, true);
		
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		mgdmfunctions = null;
		mgdmlabels = null;
		segmentation = null;
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
	
	public final void reduceMGDMsize(int nred) {
		
		float[][] redfunctions = new float[nred][nx*ny*nz];
		byte[][] redlabels = new byte[nred+1][nx*ny*nz];	
				
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			for (byte n = 0; n<nred; n++) {
         		redfunctions[n][xyz] = mgdmfunctions[n][xyz];	
         		redlabels[n][xyz] = mgdmlabels[n][xyz];	
         	}
         	redlabels[nred][xyz] = mgdmlabels[nred][xyz];
        }
        
        // replace the maps
        nmgdm = (byte)nred;
        mgdmfunctions = redfunctions;
        mgdmlabels = redlabels;
	}
	
	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][][] exportBinaryLevelSet(boolean[] inside) {
    	float[][][] levelset = new float[nx][ny][nz];
    	
    	float maximumDist = maxDist;
    	
   		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
   			int xyz = x+nx*y+nx*ny*z;
   			
   			levelset[x][y][z] = 0.0f;
   			if (mgdmlabels[0][xyz]>-1) {
				if (inside[mgdmlabels[0][xyz]]) {
					// search for the next outside value, use constant if none
					int nout = -1;
					for (int n=1;n<nmgdm && nout==-1;n++) {
						if (mgdmlabels[n][xyz]>-1 && !inside[mgdmlabels[n][xyz]]) nout = n;
					}
					if (nout>-1) {
						for (int n=0;n<nout;n++) {
							levelset[x][y][z] -= mgdmfunctions[n][xyz];
						}
					} else {
						levelset[x][y][z] = -maximumDist;
					}
				} else {
					// search for the next inside value, use constant if none
					int nin = -1;
					for (int n=1;n<nmgdm && nin==-1;n++) {
						if (mgdmlabels[n][xyz]>-1 && inside[mgdmlabels[n][xyz]]) nin = n;
					}
					if (nin>-1) {
						for (int n=0;n<nin;n++) {
							levelset[x][y][z] += mgdmfunctions[n][xyz];
						}
					} else {
						levelset[x][y][z] = maximumDist;
					}
				}    
			} else {
				// outside by default
				levelset[x][y][z] = maximumDist;	
			}
    	}
    	return levelset;
    }
	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedLevelSets() {
    	float[][] levelsets = new float[nobj][nx*ny*nz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nx*ny*nz; xyz++) {
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
    	float[][] levelsets = new float[nobj][nx*ny*nz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nx*ny*nz; xyz++) {
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
    public final float reconstructedLevelSetAt(int xyz, byte rawlb) {
    	float levelset;
    	
		if (mgdmlabels[0][xyz]>-1 && mgdmlabels[0][xyz]==rawlb) levelset = -mgdmfunctions[0][xyz];
		else  levelset = 0.0f;
    			
		for (int l=0;l<nmgdm && mgdmlabels[l][xyz]>-1 && mgdmlabels[l][xyz]!=rawlb;l++) {
			levelset += mgdmfunctions[l][xyz];
		}
     	return levelset;
    }

 	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final int[][] reconstructedLabels() {
    	int[][] labels = new int[nmgdm][nx*ny*nz];
    	for (int n=0;n<nmgdm;n++) {
    		for (int xyz=0; xyz<nx*ny*nz; xyz++) {
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
    	float[] seg = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) {
    		//if (mgdmlabels[0][xyz]>-1) {
			if (segmentation[xyz]>-1) {
				//seg[xyz] = objLabel[mgdmlabels[0][xyz]];
				seg[xyz] = objLabel[segmentation[xyz]];
			}
    	}
    	return seg;
    }
    public final float[] exportCounter() {
    	float[] ct = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) {
    		ct[xyz] = counter[xyz];
    	}
    	return ct;
    }
    public final byte[][][][] exportLabels() {
    	byte[][][][] seg = new byte[nx][ny][nz][nmgdm+1];
    	for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) for (int n=0;n<=nmgdm;n++) {
    		int xyz = x + nx*y + nx*ny*z;
    		if (mgdmlabels[n][xyz]>-1) {
				seg[x][y][z][n] = objLabel[mgdmlabels[n][xyz]];
			} else {
				seg[x][y][z][n] = -1;
			}
    	}
    	return seg;
    }
    public final float[][][][] exportFunctions() {
    	float[][][][] fun = new float[nx][ny][nz][nmgdm];
    	for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) for (int n=0;n<nmgdm;n++) {
    		int xyz = x + nx*y + nx*ny*z;
    		fun[x][y][z][n] = mgdmfunctions[n][xyz];
    	}
    	return fun;
    }
    
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingReinitialization(boolean narrowBandOnly, boolean almostEverywhere, boolean stopCounter) {
        // computation variables
        byte[] processed = new byte[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;	
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		

		long start_time = System.currentTimeMillis(); 

        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
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
							|| (narrowBandOnly && newdist<=maxDist)
							|| (almostEverywhere && (segmentation[xyzn]!=0 || newdist<=maxDist) ) ) {
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
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
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
    	float[] tmp = new float[nx*ny*nz];
    	boolean[] processed = new boolean[nx*ny*nz];
    	
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
			
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
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
			if (processed[xyz]) mgdmfunctions[0][xyz] = tmp[xyz];
			else mgdmfunctions[0][xyz] = UNKNOWN;
		}
			
        return;
    }

}

