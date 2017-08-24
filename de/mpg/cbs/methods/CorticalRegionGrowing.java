package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

/**
 *
 *  This algorithm grows geodesic regions of interest for cortical sheets
 *
 *	@version    Mar 2012
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class CorticalRegionGrowing {
	
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
	private 	float[] 		inner;  			// level set functions
	private 	float[] 		outer;  			// level set functions
	private 	float[] 		distance;  		// starting levelset
	private		boolean[]		region;
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// internal variables for computations
	int	xyz, xyzmx, xyzmy, xyzmz, xyzpx, xyzpy, xyzpz;
	int xyzpxpy, xyzpypz, xyzpzpx, xyzmxmy, xyzmymz, xyzmzmx;
	int xyzmxpy, xyzmypz, xyzmzpx, xyzpxmy, xyzpymz, xyzpzmx;
	double Dmx, Dmy, Dmz, Dpx, Dpy, Dpz, D0x, D0y, D0z;
	double SD0x, SD0y, SD0z, GPhi;
	double Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
	double K, G, tmp, K0;
	double DeltaP, DeltaM;
	double smooth, balloon;
	double balloonforce;

	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public CorticalRegionGrowing(float[] lvlsetin_, float[] lvlsetout_,
								int nx_, int ny_, int nz_,
								float rx_, float ry_, float rz_,
								boolean[] mask_,
								String connectivityType_) {
		
		inner = lvlsetin_;
		outer = lvlsetin_;
		
		
		mask = mask_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		UNKNOWN = Numerics.max(nx,ny,nz)+1;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
				
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		// init all the arrays
		try {
			distance = new float[nx*ny*nz];
			region = new boolean[nx*ny*nz];	
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
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[x+nx*y+nx*ny*z] = false;
		}
		// pre-processing?
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		inner = null;
		outer = null;
		region = null;
		distance = null;
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

	public final float[] getDistance() { return distance; }
	
	public final boolean[] getRegion() { return region; }
    
	public final byte[][][] exportRegion() {
		byte[][][] res = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (region[xyz]) res[x][y][z] = OBJ;
			else res[x][y][z] = BG;
		}
		return res;
	}

	public final float[][][] exportDistance() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = distance[xyz];
		}
		return res;
	}
    
	public final void growFromPoint(int xyz0, float maxdist) {
		
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz];
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float dist = 0.0f, newdist;
        int xyz, xyzn, xyznb;
        byte k,l;
        			        		
		// initialize the quantities
        for (xyz = 0; xyz<nx*ny*nz; xyz++) {
        	distance[xyz] = -1;
        	region[xyz] = false;
        	processed[xyz] = false;
        }
        
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();

        if (mask[xyz0]) {
			// we assume the levelset value has been recomputed on the boundary
					
			// add to the heap with previous value
			heap.addValue(0.0f,xyz0,OBJ);
		}
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty() && dist<=maxdist) {
        	// extract point with minimum distance
        	dist = heap.getFirst();
        	xyz = heap.getFirstId();
        	heap.removeFirst();

			// if already processed, move on
        	if (region[xyz])  continue;
			if (!homeomorphicLabeling(xyz)) continue;
			
			// update the distance functions at the current level
			distance[xyz] = dist;
			region[xyz]=true; // update the current level
 			
			// find new neighbors
			for (k = 0; k<6; k++) {
				xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be inside the cortical region 
					if (inner[xyzn]>0 && outer[xyzn]<0 && !region[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							// note that there is at most one value used here
							if (region[xyznb] && mask[xyznb]) {
								nbdist[l] = distance[xyznb];
								nbflag[l] = true;
							}			
						}
						newdist = minimumMarchingDistance(nbdist, nbflag);
						
						// add to the heap
						heap.addValue(newdist,xyzn,OBJ);
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
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean homeomorphicLabeling(int xyz) {
    	// if we don't check, just exit
    	if (!checkTopology) return true;
    	
		// is the new segmentation homeomorphic ? 
		if (checkComposed) if (!ObjectStatistics.isWellComposed(region,xyz,1,nx,nx*ny)) return false;		
		if (!lut.get(lut.keyFromPattern(region,xyz,1,nx,nx*ny))) return false;

		// else, it works
		return true;
    }
	
}

