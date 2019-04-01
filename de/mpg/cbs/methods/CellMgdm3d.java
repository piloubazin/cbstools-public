package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

/**
 *
 *  This algorithm uses the MGDM framework to evolve a labeling
 *	according to internal (curvature) and external (balloon, vector field) forces
 *	note: for now, we are assuming isotropic voxels 	
 *
 *	@version    Jul 2010
 *	@author     Pierre-Louis Bazin 
 *	@author     John Bogovic
 * 	@author 	Hanlin Wan
 *		
 *
 */
 
public class CellMgdm3d {
	
	// object types

	private	static	final	int	EMPTY = -1;
	
	// fast marching flags
	private final static int X = 0;
    private final static int Y = 1;
    private final static int Z = 2;
    	
	private final static int NORTH 	= 0;
    private final static int SOUTH 	= 1;
	private final static int EAST 		= 2;
    private final static int WEST 		= 3;
	private final static int FRONT 	= 4;
    private final static int BACK 		= 5;
    
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
	private 	int[][] 		mgdmlabels;   		// MGDM's label maps
	private 	int[] 			segmentation;   	// MGDM's segmentation
	private 	float[][] 		fieldforce;  		// original image forces, indep. object (e.g. boundaries)
	private 	float[][] 		balloonforces;  	// original image forces, along object normals (e.g. from memberships)
	private     float[]         intmax;
	private     float           intbg;
	private     int             bglb;
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private static	int 	  	nobj;					// total number of objects to represent (including background)
	private static	int 	  	nmgdm;					// total number of MGDM mgdmlabels and mgdmfunctions
	private 	int[]			objLabel;			// label values in the original image
	private		BinaryHeapPair	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private String	lutdir = null;
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	private		int[]			otherlabels;		// labels for the non-evolved region? 
	
	// parameters
	private	double		smoothweight, fieldweight, balloonweight, pressureweight;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 2.5f;
	private	float		narrowBandDist = landmineDist+1.8f;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
	private static class NarrowBand {
		public int[] id;
		public int[][] labels;
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
			labels = new int[nmgdm][capacity];
			functions = new float[nmgdm][capacity];
		}
		
		public void finalize() {
			capacity = -1;
			id = null;
			labels = null;
			functions = null;
		}
		
		public final void addPoint(int xyz, int[][] mgdmlb, float[][] mgdmfn) {
			// check for size
			if (currentsize>=capacity-1) {
				capacity += update;
				
				int[] oldid = id;
				int[][] oldlabels = labels;
				float[][] oldfunctions = functions;
				
				id = new int[capacity];
				labels = new int[nmgdm][capacity];
				functions = new float[nmgdm][capacity];
				
				for (int n=0;n<currentsize;n++) {
					id[n] = oldid[n];
					for (int l=0;l<nmgdm;l++) {
						labels[l][n] = oldlabels[l][n];
						functions[l][n] = oldfunctions[l][n];
					}
				}
				
				oldid = null; oldlabels = null; oldfunctions = null;
			}
			
			// add the new point (use the MGDM variables)
			id[currentsize] = xyz;
			for (int l=0;l<nmgdm;l++) {
				labels[l][currentsize] = mgdmlb[l][xyz];
				functions[l][currentsize] = mgdmfn[l][xyz];
			}
			currentsize++;
		}
		
		public void reset() {
			currentsize = 0;
		}
	}


	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public CellMgdm3d(int[] init_, int nx_, int ny_, int nz_, int nobj_, int nmgdm_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[][] balloon_, 
						float[] intmax_, float intbg_, int bglb_,
						float fw_, float bw_, float sw_, float pw_,
						String connectivityType_, String lutdir_) {
	
		init(init_, nx_, ny_, nz_, nobj_, nmgdm_, rx_, ry_, rz_, field_, balloon_,  intmax_, intbg_, bglb_, fw_, bw_, sw_, pw_, connectivityType_, lutdir_);
	}
	
	private void init(int[] init_, int nx_, int ny_, int nz_, int nobj_, int nmgdm_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[][] balloon_, 
						float[] intmax_, float intbg_, int bglb_,
						float fw_, float bw_, float sw_, float pw_,
						String connectivityType_, String lutdir_) {
		fieldforce = field_;
		balloonforces = balloon_;
		
		intmax = intmax_;
		intbg = intbg_;
		bglb = bglb_;
		
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		pressureweight = pw_;
		
		if (debug) System.out.print("MGDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing), "+pw_+" (pressure)\n");
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		nobj = nobj_;
		nmgdm = nmgdm_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		lutdir = lutdir_;
		
		objLabel = ObjectLabeling.listOrderedLabels(init_, nx, ny, nz);
		// note: we do expect that there are nb objects (not checked)
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		// init all the arrays
		try {
			mgdmfunctions = new float[nmgdm][nx*ny*nz];
			mgdmlabels = new int[nmgdm][nx*ny*nz];	
			segmentation = new int[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];
			otherlabels = new int[nx*ny*nz];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeapPair(nx*ny+ny*nz+nz*nx, BinaryHeapPair.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6")) lut = new CriticalPointLUT(lutdir, "critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26")) lut = new CriticalPointLUT(lutdir, "critical626LUT.raw.gz",200);
			else if (connectivityType_.equals("18/6")) lut = new CriticalPointLUT(lutdir, "critical186LUT.raw.gz",200);
			else if (connectivityType_.equals("6/18")) lut = new CriticalPointLUT(lutdir, "critical618LUT.raw.gz",200);
			else if (connectivityType_.equals("6/6")) lut = new CriticalPointLUT(lutdir, "critical66LUT.raw.gz",200);
			else if (connectivityType_.equals("wcs")) {
				lut = new CriticalPointLUT(lutdir, "criticalWCLUT.raw.gz",200);
				checkComposed=false;
			}
			else if (connectivityType_.equals("wco")) {
				lut = new CriticalPointLUT(lutdir, "critical66LUT.raw.gz",200);
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
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
		}
		// init decomposition
		fastMarchingInitializationFromSegmentation(init_);
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public CellMgdm3d(byte[] init_, int nx_, int ny_, int nz_, int nobj_, int nmgdm_,
						float rx_, float ry_, float rz_,
						float[][] field_, float[][] balloon_, 
						float[] intmax_, float intbg_, int bglb_,
						float fw_, float bw_, float sw_, float pw_,
						String connectivityType_, String lutdir_) {
		
		int[] tmp = new int[nx_*ny_*nz_];
		for (int n=0;n<nx_*ny_*nz_;n++) tmp[n] = init_[n];
		
		init(tmp, nx_, ny_, nz_, nobj_, nmgdm_, rx_, ry_, rz_, field_, balloon_,  intmax_, intbg_, bglb_, fw_, bw_, sw_, pw_, connectivityType_, lutdir_);
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
	
	public final int[][] getLabels() { return mgdmlabels; }
    
	public final void fastMarchingInitializationFromSegmentation(int[] init) {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // mgdm functions
			for (int n = 0; n<nmgdm; n++) {
            	mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
            // segmentation
            int nlb = EMPTY;
			for (int n=0; n<nobj; n++) {
				if (objLabel[n]==init[xyz]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xyz] = nlb;
        }
		
        // computation variables
        int[] processed = new int[nx*ny*nz]; // note: using a int instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
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
        	int xyz = heap.getFirstId1();
        	int lb = heap.getFirstId2();
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
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if (newdist<=narrowBandDist) {
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
			otherlabels[xyz] = mgdmlabels[nmgdm-1][xyz];
			for (int n=nmgdm-1;n>0;n--) {
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
     public final void fastMarchingReinitialization() {
        // computation variables
        int[] processed = new int[nx*ny*nz]; // note: using a int instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
        	// mgdm functions : reinit everywhere
			for (int n = 0; n<nmgdm; n++) {
            	if (n>0) mgdmfunctions[n][xyz] = UNKNOWN;                            
            	mgdmlabels[n][xyz] = EMPTY;
            }
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
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId1();
        	int lb = heap.getFirstId2();
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
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if (newdist<=narrowBandDist) {
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
			otherlabels[xyz] = mgdmlabels[nmgdm-1][xyz];
			for (int n=nmgdm-1;n>0;n--) {
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
    public final float[][] reconstructedLevelsetForces() {
    	float[][] forces = new float[nobj][nx*ny*nz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    			forces[n][xyz] = (float)levelsetForces(xyz,n);
    		}
    	}
    	return forces;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[] reconstructedLevelsetForceDifference(int lv) {
    	float[] forces = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		if (mgdmlabels[lv][xyz]!=EMPTY) forces[xyz] = (float)levelsetForces(xyz,mgdmlabels[lv][xyz]);
    		else forces[xyz] = 0.0f;
    		if (mgdmlabels[lv+1][xyz]!=EMPTY) forces[xyz]-= (float)levelsetForces(xyz,mgdmlabels[lv+1][xyz]);
    		
    	}
    	return forces;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedFastMarchingForces() {
    	float[][] forces = new float[nobj][nx*ny*nz];
    	for (int n=0;n<nobj;n++) {
    		for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    			forces[n][xyz] = (float)fastMarchingForces(xyz,n);
    		}
    	}
    	return forces;
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
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBand(int iter, float mindist) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	// init
    	
    	// first estimate the narrow band size
    	int size = 0;
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) if (mgdmfunctions[0][xyz]<narrowBandDist) size++;
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
			// the criterion for being in the narrow band is to have a int distance to closest boundaries
			if (mgdmfunctions[0][xyz]<narrowBandDist) {
				narrowband.addPoint(xyz, mgdmlabels, mgdmfunctions);
				// in addition, if close to the narrow band boundary, set a landmine
				if (mgdmfunctions[0][xyz]>=landmineDist) {
					landmines.set(xyz,true);
				}
			}
		}
		int[] nswap = new int[nmgdm];
			
		double[] forces = new double[nmgdm+1];
		int lbmax, lbsec;
		double curr, next;
		double curdist = mindist+1;
		// evolve until a landmine is closer than minDist of the boundaries
		for (int t=0;t<iter && curdist>mindist;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinitLM = false;
			boolean reinitOL = false;
			
			for (int lb=0;lb<nmgdm;lb++) nswap[lb] = 0;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				
				/* probably better; disabled for testing
				// select best label somehow (use the balloon forces with highest value (e.g. coming from memberships)
				int olb = lastNeighborCandidate(xyz);
				*/
				
				// select best label from previous initialization
				int olb = otherlabels[xyz];
				
				// evolve the MGDM functions
				
				// compute the forces from current levelset values, update the narrow band from it
				/*
				double curr, next;
				if (olb!=EMPTY) next = levelsetForces(xyz, olb);
				else next = levelsetForces(xyz, mgdmlabels[nmgdm-1][xyz]);
				*/
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
					/*curr = levelsetForces(xyz, mgdmlabels[lb][xyz]);*/
					curr = forces[lb];
					if (lb==lbmax) next = forces[lbsec];
					else next = forces[lbmax];
					
					// update the narrow band values, not the original data
					narrowband.functions[lb][n] += curr - next;
					
					/*next = curr;*/
	
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
									// check for boundary changes in the landmines : force reinitialization
									if (landmines.get(xyz)) reinitLM = true;
									// check for far labels getting mixed in: time to re-initialize
									if (narrowband.labels[0][n]==otherlabels[xyz]) reinitOL = true;
								} else {
									// reset to low value
									narrowband.functions[lb][n] = lowlevel;
								}
							} else {
								//if (homeomorphicLabeling(xyz, mgdmlabels[lb+1][xyz])) {
									narrowband.labels[lb+1][n] = mgdmlabels[lb][xyz];
									narrowband.labels[lb][n] = mgdmlabels[lb+1][xyz];
									narrowband.functions[lb][n] = -narrowband.functions[lb][n];
								//} else {
									// reset to low value
									//narrowband.functions[lb][n] = lowlevel;
								//}
							}
						} else {
							// reset to low value
							narrowband.functions[lb][n] = lowlevel;
						}
					}
				}
			}
			if (debug) for (int lb=0;lb<nmgdm;lb++) System.out.print("changed labels ("+lb+"): "+nswap[lb]+"\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			float avgdiff = 0.0f;
			for (int n=0; n<narrowband.currentsize;n++) {
				int xyz = narrowband.id[n];
				// measure the changes at the base level
				if (mgdmlabels[0][xyz]==narrowband.labels[0][n]) avgdiff += Numerics.abs(mgdmfunctions[0][xyz]-narrowband.functions[0][n]);
				else avgdiff += Numerics.abs(mgdmfunctions[0][xyz]+narrowband.functions[0][n]);
				
				for (int lb=0;lb<nmgdm;lb++) {
					mgdmlabels[lb][xyz] = narrowband.labels[lb][n];
					mgdmfunctions[lb][xyz] = narrowband.functions[lb][n];
				}
			}
			curdist = (avgdiff/narrowband.currentsize);
			if (debug) System.out.print("mean distance function change: "+curdist+"\n");
	
			if (reinitLM) {
				if (debug) System.out.print("re-initialization (LM: "+reinitLM+" | OL: "+reinitOL+" )\n");
        		
				//resetIsosurfaceNarrowBand(narrowband);
				resetIsosurfaceBoundary();
				fastMarchingReinitialization();
				
				// rebuild narrow band
				narrowband.reset();
				landmines.clear();
				for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
					// the criterion for being in the narrow band is to have a intdistance to closest boundaries
					if (mgdmfunctions[0][xyz]<narrowBandDist) {
						narrowband.addPoint(xyz, mgdmlabels, mgdmfunctions);
						if (mgdmfunctions[0][xyz]>=landmineDist) {
							landmines.set(xyz,true);
						}
					}
				}
			}	
     	}
		
		// end of the evolution: recompute the level sets (disabled for debugging)
		//resetIsosurfaceNarrowBand(narrowband);
		resetIsosurfaceBoundary();
		fastMarchingReinitialization();
		
        return;
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
    		int xyzn = xyz + i + j*nx + k*nx*ny;

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
			
			// build the balloon forces for different cells
			double bforce = 0.0;
			if (lb==bglb) {
			    bforce = Numerics.bounded( (2.0*intbg - balloonforces[0][xyz])/intbg, -1.0f, 1.0f);
			} else {
			    bforce = Numerics.bounded( (intmax[lb] - Numerics.abs(intmax[lb]-balloonforces[0][xyz]))/intmax[lb], -1.0f, 1.0f);
            }
			balloon = balloonweight*stepsize*(Numerics.max(bforce, 0.0)*DeltaP + Numerics.min(bforce, 0.0)*DeltaM);
		}
		
		// external force field
		double field = 0.0;
		if (Numerics.abs(fieldweight)>0) {
			field =  fieldweight*stepsize*(Numerics.max(fieldforce[X][xyz],0.0)*Dmx + Numerics.min(fieldforce[X][xyz],0.0)*Dpx
										  +Numerics.max(fieldforce[Y][xyz],0.0)*Dmy + Numerics.min(fieldforce[Y][xyz],0.0)*Dpy
										  +Numerics.max(fieldforce[Z][xyz],0.0)*Dmz + Numerics.min(fieldforce[Z][xyz],0.0)*Dpz);
		}
		
		// pressure forces of structures 1 voxel thick
		double pressure = 0.0;
		if (phi[1][1][1]<0 && Numerics.abs(pressureweight)>0) {
			for (int n=0;n<6;n++) {
				int xyzn = xyz+xoff[n]+yoff[n]+zoff[n];
				// neighbor-based pressure force
				pressure = Numerics.max(pressure, Numerics.min(Numerics.max(1.0-mgdmfunctions[1][xyzn],0.0),
																 Numerics.max(1.0-mgdmfunctions[0][xyzn],0.0)));
			}
			pressure = -pressureweight*stepsize*pressure;
		}
		
		return - (smooth - balloon - field - pressure);
    }

    /**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean homeomorphicLabeling(int xyz, int lb) {
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

		// does it change the topology of a relation between the modified object and its neighbors ?
		int  Nconfiguration = 0;
		int[] lbs = new int[26];
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if ( (i*i+j*j+l*l>0) 
				&& (segmentation[xyz+i+j*nx+l*nx*ny]!=lb) 
				&& (segmentation[xyz+i+j*nx+l*nx*ny]!=segmentation[xyz]) ) {
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
				if ( (segmentation[xyz+i+j*nx+l*nx*ny]==segmentation[xyz])
					|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[n]) ) {
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
				if ( (segmentation[xyz+i+j*nx+l*nx*ny]==lb)
					|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[n]) ) {
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
					if ( (segmentation[xyz+i+j*nx+l*nx*ny]==segmentation[xyz])
						|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[n])
						|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[m]) ) {
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
					if ( (segmentation[xyz+i+j*nx+l*nx*ny]==lb)
						|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[n]) 
						|| (segmentation[xyz+i+j*nx+l*nx*ny]==lbs[m]) ) {
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
				if (mgdmlabels[0][xyznb]!=mgdmlabels[0][xyz] && mask[xyznb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(mgdmfunctions[0][xyznb]);
					nbflag[l] = true;
					boundary = true;
				}
			}
			if (boundary) {
				narrowband.functions[0][n] = isoSurfaceDistance(mgdmfunctions[0][xyz], nbdist, nbflag);
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int n=0; n<narrowband.currentsize;n++) {
			int xyz = narrowband.id[n];
			mgdmfunctions[0][xyz] = narrowband.functions[0][n];
		}
			
        return;
    }

    /** 
    *   isosurface distance re-initialization at the boundary
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
    
    /** 
    *  	Evolution using a joint fast marching scheme
    */
    public final void evolveFastMarching(int iter, int sign) {
    	
    	if (debug) System.out.print("level set evolution: fast marching\n");

    	// computation variables
        boolean[] processed = new boolean[nx*ny*nz]; // note: using a int instead of boolean for the second pass
    	
        // init
    	
    	// first estimate the binary tree size
    	int size = 0;
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) if (mgdmfunctions[0][xyz]<=0.5f) size++;
		// create the tree with initial estimates of size
    	heap = new BinaryHeapPair(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeapPair.MINTREE);
    	
		for (int t=0;t<iter;t++) {
			if (debug) System.out.print("iteration "+t+"\n");

			
			if (debug) System.out.print("init ("+size+")\n");
        
			heap.reset();
			double curr, next;
			float speed = 0;
			int nbound = 0;
			for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
				// the criterion for being in the narrow band is to have a int distance to closest boundaries
				if (mgdmfunctions[0][xyz]<1.0f) {
					// add to the heap
					curr = balloonweight*balloonforces[mgdmlabels[0][xyz]][xyz];
					next = balloonweight*balloonforces[mgdmlabels[1][xyz]][xyz];
					if (next>curr) {
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
			System.out.println("initial boundary points: "+nbound);
				
			// grow the labels and functions			
			float dist = 0.0f;
			while (heap.isNotEmpty()) {
				// extract point with minimum distance
				dist = heap.getFirst();
				int xyz = heap.getFirstId1();
				int lb = heap.getFirstId2();
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
							if (mgdmlabels[0][xyzn]!=lb) {
								curr = balloonweight*balloonforces[mgdmlabels[0][xyzn]][xyzn];
								next = balloonweight*balloonforces[lb][xyzn];
								if (next>curr) {
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
			for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (!processed[xyz]) {
				mgdmfunctions[0][xyz] = dist;
			}
			// re-init
			//resetIsosurfaceBoundary();
			//fastMarchingReinitialization(false);
		}
        return;
    }
    
    private final float fastMarchingForces(int xyz, int lb) {
    	if (balloonforces[lb][xyz]>0) return 1.0f/(float)(smoothweight+balloonweight*balloonforces[lb][xyz]);
    	else return 1.0f/(float)smoothweight;
    }
}

