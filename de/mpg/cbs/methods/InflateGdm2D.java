package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;


/**
 *
 *  This algorithm uses the classical GDM framework to inflate a single object
 *	to decrease curvature.
 *	note: for now, we are assuming isotropic voxels, and boundary conditions are not enforced
 *	(not a problem unless the narrowband reaches a boundary of the image, 
 *	mirroring conditions would solve this problem if needed)
 *
 *	@version    Nov 2012
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class InflateGdm2D {
	
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

	// data and membership buffers
	private 	float[] 		levelset;  			// level set functions
	private 	byte[] 			segmentation;   	// segmentation
	private 	float[] 		inlevelset;  		// starting levelset
	private 	float[] 		trglevelset;  		// target (smoothed) levelset
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
		
	private	double		smoothweight, balloonweight;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 5.0f;
	private	float		narrowBandDist = landmineDist+1.8f;	
	
	// internal variables for computations
	int	xyz, xyzmx, xyzmy, xyzmz, xyzpx, xyzpy, xyzpz;
	int xyzpxpy, xyzpypz, xyzpzpx, xyzmxmy, xyzmymz, xyzmzmx;
	int xyzmxpy, xyzmypz, xyzmzpx, xyzpxmy, xyzpymz, xyzpzmx;
	double Dmx, Dmy, Dmz, Dpx, Dpy, Dpz, D0x, D0y, D0z;
	double SD0x, SD0y, SD0z, GPhi;
	double Dxx, Dyy, Dzz, Dxy, Dyz, Dzx;
	double K, G, tmp;
	double DeltaP, DeltaM;
	double smooth, balloon;
	double balloonforce;

	// for debug and display
	private static final boolean		debug=false;
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
	public InflateGdm2D(float[] lvlsetin_, int nx_, int ny_, int nz_, float rx_, float ry_, float rz_,
						boolean[] mask_, float bw_, float sw_, String connectivityType_, String lutdir_) {
			
		inlevelset = lvlsetin_;
		trglevelset = lvlsetin_;
		
		balloonweight = bw_;
		smoothweight = sw_;
		
		mask = mask_;
		
		if (debug) System.out.print("GDM forces: "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
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
		
		// init all the arrays
		try {
			levelset = new float[nx*ny*nz];
			segmentation = new byte[nx*ny*nz];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT(lutdir_,"critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT(lutdir_,"critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT(lutdir_,"critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT(lutdir_,"critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT(lutdir_,"critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT(lutdir_,"critical66LUT.raw.gz",200);
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
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2) mask[x+nx*y+nx*ny*z] = false;
		}
		
		// init levelset evolution quantities
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			levelset[xyz] = inlevelset[xyz];                            
           
			// segmentation
            if (inlevelset[xyz]<0) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
        }
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		levelset = null;
		segmentation = null;
		inlevelset = null;
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
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) res[xyz] = segmentation[xyz];
		return res;
	}

	public final float[][][] exportLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			res[x][y][z] = levelset[x+y*nx+z*ny*nx];
		}
		return res;
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
    
    public final void smoothLevelset(float scale) {
    	trglevelset = ImageFilters.separableMaskedConvolution(trglevelset, mask, nx, ny, nz, 
    															ImageFilters.separableGaussianKernel(scale,scale,0));
    		
    	return;
    }
    
    public final void updateTarget() {
    	trglevelset = levelset;
    }
    
    public final void recursiveLevelsetInflation(int steps, float scale, int iter, float mindiff) {
    	for (int t=0;t<steps;t++) {
    		smoothLevelset(scale);
    		evolveNarrowBand(iter, mindiff);
    		updateTarget();
    	}
    	return;
    }
        
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingReinitialization(boolean narrowBandOnly) {
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[4];
		boolean[] nbflag = new boolean[4];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
        	processed[xyz] = false;
        	// search for boundaries
        	for (int k = 0; k<4; k++) {
				int xyzn = xyz + xoff[k] + yoff[k];
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
        while (heap.isNotEmpty() && (!narrowBandOnly || maxdist<=narrowBandDist+SQR2)) {
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
 			maxdist = dist;	// keep track of distance if stopping at the narrow band
			
			// find new neighbors
			for (int k = 0; k<4; k++) {
				int xyzn = xyz + xoff[k] + yoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (segmentation[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<4; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l];
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
		// make sure the outside is set to max distance
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (!mask[xyz]) {
			levelset[xyz] = maxdist;
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

        for (int n=0; n<4; n+=2) {
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
        
        for (int n=0; n<4; n+=2) {
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
    	
    	// create a second narrow band for critical points
    	NarrowBand criticalband = new NarrowBand(Numerics.ceil(0.1f*size), Numerics.ceil(0.1f*size));
    	
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
		for (int t=0;t<iter && diff>mindiff;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinit = false;
			int nswap=0;
			double prev;
			
			// reset the critical band every time
			criticalband.reset();
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				prev = narrowband.levels[n];
				narrowband.levels[n] -= Numerics.bounded(levelsetForces(xyz),-0.9,0.9);
				
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
						// add to the critical points
						criticalband.addPoint(xyz, segmentation, levelset);
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
					} else {
						narrowband.levels[n] = -lowlevel;
						// add to the critical points
						criticalband.addPoint(xyz, segmentation, levelset);
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
			}
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
		
			// loop over the critical points
			if (debug) System.out.print("critical points ("+criticalband.currentsize+")\n");
			
			int ncrit=criticalband.currentsize;
			int tc=0;
			while (ncrit>0 && tc<10) {
				ncrit=0;
				tc++;
				for (int n=0; n<criticalband.currentsize;n++) {
					xyz = criticalband.id[n];
					
					// compute the forces from current levelset values, update the narrow band from it
					prev = criticalband.levels[n];
					criticalband.levels[n] -= Numerics.bounded(levelsetForces(xyz),-0.9,0.9);
					
					// change of sign ? bg -> obj
					if (criticalband.levels[n]<0 && prev>0) {
						//if (debug) System.out.print(""+lb);
						// switch labels	
						if (homeomorphicLabeling(xyz, OBJ)) {
							ncrit++;
						
							// update the segmentation
							segmentation[xyz] = OBJ;
						} else {
							criticalband.levels[n] = lowlevel;
						}
					} else if (criticalband.levels[n]>0 && prev<0) {
						//if (debug) System.out.print(""+lb);
						// switch labels	
						if (homeomorphicLabeling(xyz, BG)) {
							ncrit++;
						
							// update the segmentation
							segmentation[xyz] = BG;
						} else {
							criticalband.levels[n] = -lowlevel;
						}
					}
				}
				if (debug) System.out.print("changed labels: "+ncrit+" ("+(ncrit*100.0f/criticalband.currentsize)+" % of boundary)\n");
			}
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
			}
			// and then change the critical ones
			for (int n=0; n<criticalband.currentsize;n++) {
				xyz = criticalband.id[n];
				levelset[xyz] = criticalband.levels[n];
			}
			
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				boundarysize = 0;
				for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
					// the criterion for being in the narrow band is to have a shortdistance to closest boundaries
					if (Numerics.abs(levelset[xyz])<narrowBandDist) {
						narrowband.addPoint(xyz, segmentation, levelset);
						if (Numerics.abs(levelset[xyz])>=landmineDist) {
							landmines.set(xyz,true);
						} 
					}
					if (Numerics.abs(levelset[xyz])<1.0) {
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
    	
    	xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
    	xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
    	xyzpz = xyz;
    	xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
    	xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
    	xyzmz = xyz;
    	
    	xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
    	xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
    	xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
    	xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
    	
    	xyzpypz = xyz;
    	xyzpymz = xyz;
    	xyzmypz = xyz;
    	xyzmymz = xyz;
    	
    	xyzpzpx = xyz;
    	xyzpzmx = xyz;
    	xyzmzpx = xyz;
    	xyzmzmx = xyz;

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
		
		// curvature smoothing
		smooth = smoothweight*stepsize*K*GPhi;
		
		// balloon forces to stay close to the original surface
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
								 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
								 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
		
		DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
								 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
								 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
		
		balloonforce = -Numerics.bounded( trglevelset[xyz] - levelset[xyz], -1.0, 1.0);
		
		balloon = balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*DeltaP + Numerics.min(balloonforce, 0.0)*DeltaM);
		
		return balloon - smooth;
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
			if (segmentation[xyz+i+j*nx]==lb) {
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
			if (segmentation[xyz+i+j*nx]==segmentation[xyz]) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = false;
		
		if (checkComposed) if (!ObjectStatistics.isWellComposed(obj,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		
		return true;
    }
	
    /** 
    *   isosurface distance re-initialization inside the narrow band 
    *	(assumes a correct narrow band structure enclosing the boundary)
    */
    private final void resetIsosurfaceNarrowBand(NarrowBand narrowband) {
    	if (debug) System.out.print("level set evolution: iso-surface reinit\n");

    	float[] nbdist = new float[4];
    	boolean[] nbflag = new boolean[4];
    	boolean boundary;
    	
		for (int n=0; n<narrowband.currentsize;n++) {
			int xyz = narrowband.id[n];
			
			boundary = false;
			for (int l=0; l<4; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				int xyznb = xyz + xoff[l] + yoff[l];
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

