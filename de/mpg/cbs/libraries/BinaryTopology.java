package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm handles topology correction and other constraints
 *	for binary objects
 *
 *	@version    Jul 2010
 *	@author     Pierre-Louis Bazin 
 *	@author     John Bogovic
 * 	@author 	Hanlin Wan
 *		
 *
 */
 
public class BinaryTopology {
	
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
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private 	int[] 		object;   	// segmentation
	private 	byte[] 			segmentation;   	// segmentation
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)

	// computation variables
	boolean[][][] obj = new boolean[3][3][3];
			
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	

	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public BinaryTopology(int[] init_, int nx_, int ny_, int nz_,
							float rx_, float ry_, float rz_,
							String connectivityType_) {
	
		this(init_, nx_, ny_, nz_, rx_, ry_, rz_, connectivityType_, null);
	}
		
	public BinaryTopology(int[] init_, int nx_, int ny_, int nz_,
							float rx_, float ry_, float rz_,
							String connectivityType_, String connectivityPath_) {
		
		object = init_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
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
				lut = new CriticalPointLUT(connectivityPath_, "criticalWCLUT.raw.gz",200);
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
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
		}

		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		

	public void finalize() {
		segmentation = null;
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

	public final byte[] getSegmentation() { return segmentation; }
    
	public final float[] exportFloatSegmentation() {
    	float[] seg = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		seg[xyz] = (float)segmentation[xyz];
    	}
    	return seg;
    }
	public final int[] exportIntSegmentation() {
    	int[] seg = new int[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		seg[xyz] = (int)segmentation[xyz];
    	}
    	return seg;
    }
	public final void insideSphericalTopology() {
		
		// from the boundary, find the middle (approx)
		if (debug) BasicInfo.displayMessage("fast marching init\n");		
        heap.reset();
        boolean[] processed = new boolean[nx*ny*nz];
        float[] sqdist = new float[nx*ny*nz];
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	sqdist[xyz] = 0;
        	if (mask[xyz]) {
        		processed[xyz] = false;
        		if (object[xyz]==OBJ) {
					// search for boundaries
					for (int k = 0; k<6; k++) {
						int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
						if (object[xyzn]!=object[xyz]) if (mask[xyzn]) {
							
							// add to the heap
							heap.addValue(1.0f,xyzn,OBJ);
						}
					}
				}
			} else {
				processed[xyz] = true;
			}
        }
		if (debug) BasicInfo.displayMessage("init\n");		
		
        // grow the distance to boundary
		float maxdist = 0.0f;
        int xyzstart = -1;
		while (heap.isNotEmpty()) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			sqdist[xyz] = dist;
			processed[xyz]=true; // update the current level
			
			if (dist>maxdist) {
				maxdist = dist;
				xyzstart = xyz;
			}
			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (object[xyzn]==OBJ && !processed[xyzn]) {
						// add to the heap
						heap.addValue(dist+1.0f,xyzn,OBJ);
					}
				}
			}
		}
		
		// redo with topology constraints, from the center
		if (debug) BasicInfo.displayMessage("fast marching regrowth\n");		
        heap.reset();
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	sqdist[xyz] = 0;
        	if (mask[xyz]) {
        		processed[xyz] = false;
			} else {
				processed[xyz] = true;
			}
			segmentation[xyz] = BG;
        }
        processed[xyzstart] = true;
        sqdist[xyzstart] = 0;
        segmentation[xyzstart] = OBJ;
        
        for (int k = 0; k<6; k++) {
			int xyzn = xyzstart + xoff[k] + yoff[k] + zoff[k];
			if (mask[xyzn]) {
				// add to the heap
				heap.addValue(1.0f,xyzn,OBJ);
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
        // grow back to the boundary
		while (heap.isNotEmpty()) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// check for topology
			if (isSimplePoint(xyz, OBJ)) {
				// update the distance functions at the current level
				sqdist[xyz] = dist;
				processed[xyz]=true; // update the current level
				segmentation[xyz] = OBJ;
				
				// find new neighbors
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					
					if (mask[xyzn]) {
						// must be in outside the object or its processed neighborhood
						if (object[xyzn]==OBJ && !processed[xyzn]) {
							// add to the heap
							heap.addValue(dist+1.0f,xyzn,OBJ);
						}
					}
				}
			}
		}
		
		
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }
     
	public final void outsideSphericalTopology() {
		
		// from the image boundary, to the object
		if (debug) BasicInfo.displayMessage("fast marching init\n");		
        heap.reset();
        boolean[] processed = new boolean[nx*ny*nz];
        float[] sqdist = new float[nx*ny*nz];
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	sqdist[xyz] = 0;
        	if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) {
        		segmentation[xyz] = OBJ;
        		processed[xyz] = false;
        	} else if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) { 
				// boundary
        		segmentation[xyz] = OBJ;
        		processed[xyz] = false;
        		// add to the heap
				heap.addValue(1.0f,xyz,BG);
			} else {
				segmentation[xyz] = BG;
        		processed[xyz] = true;
			}
        }
		if (debug) BasicInfo.displayMessage("init\n");		
				
        // grow to the object boundary
		while (heap.isNotEmpty()) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// check for topology
			if (isSimplePoint(xyz, BG)) {
				// update the distance functions at the current level
				sqdist[xyz] = dist;
				processed[xyz]=true; // update the current level
				segmentation[xyz] = BG;
				
				// find new neighbors
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					
					// must be in outside the object or its processed neighborhood
					if (object[xyzn]==BG && !processed[xyzn]) {
						// add to the heap
						heap.addValue(dist+1.0f,xyzn,BG);
					}
				}
			}
		}
		
		
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }
     
    
    /**
	 *  simple point criterion
	 */
    private final boolean isSimplePoint(int xyz, byte lb) {
    	// if we don't check, just exit
    	if (!checkTopology) return true;
    	
		// is the new segmentation homeomorphic ? 
		
		// inside the original object ?
		if (segmentation[xyz]==lb) return true;
		
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
			
		return true;
    }

}

