package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm uses the MGDM framework to evolve a labeling
 *	according to internal (curvature) and external (balloon, vector field) forces
 *	note: here we assume 2D images
 *
 *	@version    Jul 2010
 *	@author     Pierre-Louis Bazin 
 *	@author     John Bogovic
 * 	@author 	Hanlin Wan
 *		
 *
 */
 
public class Mgdm2d {
	
	// object types

	private	static	final	int	EMPTY = -1;
	
	// fast marching flags
	private final static int X = 0;
    private final static int Y = 1;
    
    // numerical quantities
	private static final	float   INF=1e15f;
	private static final	float   ZERO=1e-15f;
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private final static float SQR2 = (float) Math.sqrt(2.0f);
    private final static float diagdist = 1/(2*SQR2);
	private static final	float	UNKNOWN = -1.0f;
	
	private static int[] xoff;
    private static int[] yoff;


	// data and membership buffers
	private 	float[][] 		mgdmfunctions;  	// MGDM's pseudo level set mgdmfunctions
	private 	int[][] 		mgdmlabels;   		// MGDM's label maps
	private 	int[] 			segmentation;   	// MGDM's segmentation
	private 	float[][] 		fieldforce;  		// original image forces, indep. object (e.g. boundaries)
	private 	float[][] 		balloonforces;  	// original image forces, along object normals (e.g. from memberships)
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny;   		    // images dimensions
	private static	float 		rx,ry;   		    // images resolutions
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
		
		public final void addPoint(int xy, int[][] mgdmlb, float[][] mgdmfn) {
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
			id[currentsize] = xy;
			for (int l=0;l<nmgdm;l++) {
				labels[l][currentsize] = mgdmlb[l][xy];
				functions[l][currentsize] = mgdmfn[l][xy];
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
	public Mgdm2d(int[] init_, int nx_, int ny_, int nobj_, int nmgdm_,
						float rx_, float ry_,
						float[][] field_, float[][] balloon_, 
						float fw_, float bw_, float sw_, float pw_,
						String connectivityType_, String lutdir_) {
	
		fieldforce = field_;
		balloonforces = balloon_;
		
		fieldweight = fw_;
		balloonweight = bw_;
		smoothweight = sw_;
		pressureweight = pw_;
		
		if (debug) System.out.print("MGDM forces: "+fw_+" (field), "+bw_+" (balloon), "+sw_+" (smoothing), "+pw_+" (pressure)\n");
		
		nx = nx_;
		ny = ny_;
		
		nobj = nobj_;
		nmgdm = nmgdm_;
		
		rx = rx_;
		ry = ry_;
		
		lutdir = lutdir_;
		
		objLabel = ObjectLabeling.listOrderedLabels(init_, nx, ny);
		// note: we do expect that there are nb objects (not checked)
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0};
		yoff = new int[]{0, 0, nx, -nx};
		
		// init all the arrays
		try {
			mgdmfunctions = new float[nmgdm][nx*ny];
			mgdmlabels = new int[nmgdm][nx*ny];	
			segmentation = new int[nx*ny];	
			mask = new boolean[nx*ny];
			otherlabels = new int[nx*ny];	
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeapPair(nx*ny, BinaryHeapPair.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectivityType_.equals("26/6") || connectivityType_.equals("8/4")) lut = new CriticalPointLUT(lutdir, "critical266LUT.raw.gz",200);
			else if (connectivityType_.equals("6/26") || connectivityType_.equals("4/8")) lut = new CriticalPointLUT(lutdir, "critical626LUT.raw.gz",200);
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
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) {
			if (x>1 && x<nx-2 && y>1 && y<ny-2) mask[x+nx*y] = true;
			else mask[x+nx*y] = false;
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
	
	public final int[][] getLabels() { return mgdmlabels; }
    
	public final void fastMarchingInitializationFromSegmentation(int[] init) {
         // initialize the quantities
         for (int xy = 0; xy<nx*ny; xy++) {
         	 // mgdm functions
			for (int n = 0; n<nmgdm; n++) {
            	mgdmfunctions[n][xy] = UNKNOWN;                            
            	mgdmlabels[n][xy] = EMPTY;
            }
            // segmentation
            int nlb = EMPTY;
			for (int n=0; n<nobj; n++) {
				if (objLabel[n]==init[xy]) {
					nlb = n;
					continue;
				}
			}
			segmentation[xy] = nlb;
        }
		
        // computation variables
        int[] processed = new int[nx*ny]; // note: using a int instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
        	processed[xy] = 0;
        	// search for boundaries
        	for (int k = 0; k<4; k++) {
				int xyn = xy + xoff[k] + yoff[k];
				if (segmentation[xyn]!=segmentation[xy]) if (mask[xyn]) {
					
					// add to the heap
					heap.addValue(0.5f,xyn,segmentation[xy]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xy = heap.getFirstId1();
        	int lb = heap.getFirstId2();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xy]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			boolean done = false;
			for (int n=0; n<processed[xy]; n++)
				if (mgdmlabels[n][xy]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xy]][xy] = dist;
			mgdmlabels[processed[xy]][xy] = lb;
			processed[xy]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<4; k++) {
				int xyn = xy + xoff[k] + yoff[k];
				
				if (mask[xyn]) {
					// must be in outside the object or its processed neighborhood
					boolean isprocessed = false;
					if (segmentation[xyn]==lb) isprocessed = true;
					else {
						for (int n=0; n<processed[xyn]; n++)
							if (mgdmlabels[n][xyn]==lb) isprocessed = true;
					}
					
					if (!isprocessed) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<4; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xynb = xyn + xoff[l] + yoff[l];
							// note that there is at most one value used here
							for (int n=0; n<processed[xynb]; n++) if (mask[xynb]) if (mgdmlabels[n][xynb]==lb) {
								nbdist[l] = mgdmfunctions[n][xynb];
								nbflag[l] = true;
							}			
						}
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if (newdist<=narrowBandDist) {
							// add to the heap
							heap.addValue(newdist,xyn,lb);
						}
					}
				}
			}
		}
		// to create the MGDM functions, we need to copy the segmentation, forget the last labels
		// and compute differences between distance functions
		if (debug) BasicInfo.displayMessage("transform into MGDM functions\n");		
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
			// label permutation
			otherlabels[xy] = mgdmlabels[nmgdm-1][xy];
			for (int n=nmgdm-1;n>0;n--) {
				mgdmlabels[n][xy] = mgdmlabels[n-1][xy];
			}
			mgdmlabels[0][xy] = segmentation[xy];
			
			// distance function difference
        	for (int n = nmgdm-1; n>0; n--) {
        		mgdmfunctions[n][xy] = Numerics.max(UNKNOWN, mgdmfunctions[n][xy]
        														-mgdmfunctions[n-1][xy]);
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
        int[] processed = new int[nx*ny]; // note: using a int instead of boolean for the second pass
		float[] nbdist = new float[4];
		boolean[] nbflag = new boolean[4];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
        	// mgdm functions : reinit everywhere
			for (int n = 0; n<nmgdm; n++) {
            	if (n>0) mgdmfunctions[n][xy] = UNKNOWN;                            
            	mgdmlabels[n][xy] = EMPTY;
            }
            processed[xy] = 0;
        	// not needed, should be kept the same
        	//segmentation[xy] = mgdmlabels[0][xy];
        	// search for boundaries
        	for (int k = 0; k<4; k++) {
				int xyn = xy + xoff[k] + yoff[k];
				if (segmentation[xyn]!=segmentation[xy]) if (mask[xyn]) {
					
					// add to the heap with previous value
					heap.addValue(mgdmfunctions[0][xyn],xyn,segmentation[xy]);
                }
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xy = heap.getFirstId1();
        	int lb = heap.getFirstId2();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xy]>=nmgdm)  continue;
			
			// if there is already a label for this object, this is done
			boolean done = false;
			for (int n=0; n<processed[xy]; n++)
				if (mgdmlabels[n][xy]==lb) done = true;
			if (done) continue;
			
			// update the distance functions at the current level
			mgdmfunctions[processed[xy]][xy] = dist;
			mgdmlabels[processed[xy]][xy] = lb;
			processed[xy]++; // update the current level
 			
			// find new neighbors
			for (int k = 0; k<4; k++) {
				int xyn = xy + xoff[k] + yoff[k];
				
				if (mask[xyn]) {
					// must be in outside the object or its processed neighborhood
					boolean isprocessed = false;
					if (segmentation[xyn]==lb) isprocessed = true;
					else {
						for (int n=0; n<processed[xyn]; n++)
							if (mgdmlabels[n][xyn]==lb) isprocessed = true;
					}
					
					if (!isprocessed) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<4; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xynb = xyn + xoff[l] + yoff[l];
							// note that there is at most one value used here
							for (int n=0; n<processed[xynb]; n++) if (mask[xynb]) if (mgdmlabels[n][xynb]==lb) {
								nbdist[l] = mgdmfunctions[n][xynb];
								nbflag[l] = true;
							}			
						}
						float newdist = minimumMarchingDistance(nbdist, nbflag);
						
						if (newdist<=narrowBandDist) {
							// add to the heap
							heap.addValue(newdist,xyn,lb);
						}
					}
				}
			}
		}
		// to create the MGDM functions, we need to copy the segmentation, forget the last labels
		// and compute differences between distance functions
		if (debug) BasicInfo.displayMessage("transform into MGDM functions\n");		
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
			// label permutation
			otherlabels[xy] = mgdmlabels[nmgdm-1][xy];
			for (int n=nmgdm-1;n>0;n--) {
				mgdmlabels[n][xy] = mgdmlabels[n-1][xy];
			}
			mgdmlabels[0][xy] = segmentation[xy];
			
			// distance function difference
        	for (int n = nmgdm-1; n>0; n--) {
        		mgdmfunctions[n][xy] = Numerics.max(UNKNOWN, mgdmfunctions[n][xy]
        														-mgdmfunctions[n-1][xy]);
        	}
        }
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }

	/**
     * the Fast marching distance computation 
     * (!assumes a 4D array with opposite coordinates stacked one after the other)
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
        
        tmp = (s+(float) Math.sqrt((double) (s*s-count*(s2-1.0f))))/count;

        // The larger root
        return tmp;
    }
	/**
     * the isosurface distance computation 
     * (!assumes a 4D array with opposite coordinates stacked one after the other)
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
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedLevelSets() {
    	float[][] levelsets = new float[nobj][nx*ny];
    	for (int n=0;n<nobj;n++) {
    		for (int xy=0; xy<nx*ny; xy++) {
    			if (mgdmlabels[0][xy]==n) levelsets[n][xy] = -mgdmfunctions[0][xy];
    			else  levelsets[n][xy] = 0.0f;
    			
    			for (int l=0;l<nmgdm && mgdmlabels[l][xy]!=n;l++) {
    				levelsets[n][xy] += mgdmfunctions[l][xy];
    			}
    		}
    	}
    	return levelsets;
    }
    /**
     *	reconstruct the levelset only where it is guaranteed to be exact 
     */
    public final float[][] reconstructedExactLevelSets() {
    	float[][] levelsets = new float[nobj][nx*ny];
    	for (int n=0;n<nobj;n++) {
    		for (int xy=0; xy<nx*ny; xy++) {
    			if (mgdmlabels[0][xy]==n) levelsets[n][xy] = -mgdmfunctions[0][xy];
    			else  levelsets[n][xy] = 0.0f;
    			
    			int max=0;
    			for (int l=0;l<nmgdm && mgdmlabels[l][xy]!=n;l++) {
    				levelsets[n][xy] += mgdmfunctions[l][xy];
    				max++;
    			}
    			if (max==nmgdm) levelsets[n][xy] = UNKNOWN;
    		}
    	}
    	return levelsets;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedLevelsetForces() {
    	float[][] forces = new float[nobj][nx*ny];
    	for (int n=0;n<nobj;n++) {
    		for (int xy=0; xy<nx*ny; xy++) if (mask[xy]) {
    			forces[n][xy] = (float)levelsetForces(xy,n);
    		}
    	}
    	return forces;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[] reconstructedLevelsetForceDifference(int lv) {
    	float[] forces = new float[nx*ny];
    	for (int xy=0; xy<nx*ny; xy++) if (mask[xy]) {
    		if (mgdmlabels[lv][xy]!=EMPTY) forces[xy] = (float)levelsetForces(xy,mgdmlabels[lv][xy]);
    		else forces[xy] = 0.0f;
    		if (mgdmlabels[lv+1][xy]!=EMPTY) forces[xy]-= (float)levelsetForces(xy,mgdmlabels[lv+1][xy]);
    		
    	}
    	return forces;
    }
	/** 
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final float[][] reconstructedFastMarchingForces() {
    	float[][] forces = new float[nobj][nx*ny];
    	for (int n=0;n<nobj;n++) {
    		for (int xy=0; xy<nx*ny; xy++) if (mask[xy]) {
    			forces[n][xy] = (float)fastMarchingForces(xy,n);
    		}
    	}
    	return forces;
    }
	/**
	 *	reconstruct the level set functions with possible approximation far from contour 
	 */
    public final int[][] reconstructedLabels() {
    	int[][] labels = new int[nmgdm][nx*ny];
    	for (int n=0;n<nmgdm;n++) {
    		for (int xy=0; xy<nx*ny; xy++) {
    			if (mgdmlabels[n][xy]>-1) {
					labels[n][xy] = objLabel[mgdmlabels[n][xy]];
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
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) if (mgdmfunctions[0][xy]<narrowBandDist) size++;
		// create the narrow band with initial estimates of size
    	NarrowBand narrowband = new NarrowBand(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size));
    	BitSet landmines = new BitSet(Numerics.ceil(0.2f*size));
    	
    	if (debug) System.out.print("init ("+size+")\n");
        
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
			// the criterion for being in the narrow band is to have a short distance to closest boundaries
			if (mgdmfunctions[0][xy]<narrowBandDist) {
				narrowband.addPoint(xy, mgdmlabels, mgdmfunctions);
				// in addition, if close to the narrow band boundary, set a landmine
				if (mgdmfunctions[0][xy]>=landmineDist) {
					landmines.set(xy,true);
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
				int xy = narrowband.id[n];
				
				/* probably better; disabled for testing
				// select best label somehow (use the balloon forces with highest value (e.g. coming from memberships)
				int olb = lastNeighborCandidate(xy);
				*/
				
				// select best label from previous initialization
				int olb = otherlabels[xy];
				
				// evolve the MGDM functions
				
				// compute the forces from current levelset values, update the narrow band from it
				/*
				double curr, next;
				if (olb!=EMPTY) next = levelsetForces(xy, olb);
				else next = levelsetForces(xy, mgdmlabels[nmgdm-1][xy]);
				*/
				if (olb!=EMPTY) forces[nmgdm] = levelsetForces(xy, olb);
				else forces[nmgdm] = 0.0f;
				lbmax = nmgdm;
				lbsec = -1;
				for (int lb=nmgdm-1;lb>=0;lb--) {
					if (mgdmlabels[lb][xy]!=EMPTY) {
						forces[lb] = levelsetForces(xy, mgdmlabels[lb][xy]);
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
				for (int lb=nmgdm-1;lb>=0;lb--) if (mgdmlabels[lb][xy]!=EMPTY) {
					/*curr = levelsetForces(xy, mgdmlabels[lb][xy]);*/
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
								//if (homeomorphicLabeling(xy, olb)) {
								narrowband.labels[lb][n] = olb;
								narrowband.functions[lb][n] = -narrowband.functions[lb][n];
							} else {
								// reset to low value
								narrowband.functions[lb][n] = lowlevel;
							}
						} else if (mgdmlabels[lb+1][xy]!=EMPTY) {
							if (lb==0) {
								// check for topology here (optional)
								if (homeomorphicLabeling(xy, mgdmlabels[lb+1][xy])) {
									narrowband.labels[lb+1][n] = mgdmlabels[lb][xy];
									narrowband.labels[lb][n] = mgdmlabels[lb+1][xy];
									narrowband.functions[lb][n] = -narrowband.functions[lb][n];
									segmentation[xy] = mgdmlabels[lb+1][xy];
									// check for boundary changes in the landmines : force reinitialization
									if (landmines.get(xy)) reinitLM = true;
									// check for far labels getting mixed in: time to re-initialize
									if (narrowband.labels[0][n]==otherlabels[xy]) reinitOL = true;
								} else {
									// reset to low value
									narrowband.functions[lb][n] = lowlevel;
								}
							} else {
								//if (homeomorphicLabeling(xy, mgdmlabels[lb+1][xy])) {
									narrowband.labels[lb+1][n] = mgdmlabels[lb][xy];
									narrowband.labels[lb][n] = mgdmlabels[lb+1][xy];
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
				int xy = narrowband.id[n];
				// measure the changes at the base level
				if (mgdmlabels[0][xy]==narrowband.labels[0][n]) avgdiff += Numerics.abs(mgdmfunctions[0][xy]-narrowband.functions[0][n]);
				else avgdiff += Numerics.abs(mgdmfunctions[0][xy]+narrowband.functions[0][n]);
				
				for (int lb=0;lb<nmgdm;lb++) {
					mgdmlabels[lb][xy] = narrowband.labels[lb][n];
					mgdmfunctions[lb][xy] = narrowband.functions[lb][n];
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
				for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
					// the criterion for being in the narrow band is to have a shortdistance to closest boundaries
					if (mgdmfunctions[0][xy]<narrowBandDist) {
						narrowband.addPoint(xy, mgdmlabels, mgdmfunctions);
						if (mgdmfunctions[0][xy]>=landmineDist) {
							landmines.set(xy,true);
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
    
    /** method to select the nmgdm-th closest neighbor based on external information 
     *  (highly application-dependent, may not be usable in all cases) 
     */
    private final int lastNeighborCandidate(int xy) {
		// select best label somehow (use the balloon forces with highest value (e.g. coming from memberships)
		int olb = EMPTY;
		if (balloonweight>0) {
			float best=-INF;
			for (int n=1;n<nobj;n++) {
				boolean outside = true;
				for (int l=0;l<nmgdm;l++) if (mgdmlabels[l][xy]==n) outside = false;
				if (outside && balloonforces[n][xy]>best) olb = n;
			}
			if (olb==EMPTY) System.out.print("?");
		}
		return olb;
	}
    
	/** specific forces applied to the level sets (application dependent) */
    private final double levelsetForces(int xy, int lb) {
    	
		// simple option: rebuild level set locally
		// note: we go back to the convention of usual level sets with negative value inside, positive value outside
		float[][] phi = new float[3][3];
		
		// do the center point first
		if (mgdmlabels[0][xy]==lb) phi[1][1] = -mgdmfunctions[0][xy];
		else  phi[1][1] = 0.0f;
		for (int l=0;l<nmgdm && mgdmlabels[l][xy]!=lb;l++) {
			phi[1][1] += mgdmfunctions[l][xy];
		}
		// neighbors
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) if (i*i>0 || j*j>0) {
    		int xyn = xy + i + j*nx;

			if (mask[xyn]) {
				if (mgdmlabels[0][xyn]==lb) phi[i+1][j+1] = -mgdmfunctions[0][xyn];
				else  phi[i+1][j+1] = 0.0f;
					
				for (int l=0;l<nmgdm && mgdmlabels[l][xyn]!=lb;l++) {
					phi[i+1][j+1] += mgdmfunctions[l][xyn];
				}
			} else {
				// filling in values outside the mask?? center value
				phi[i+1][j+1] = phi[1][1];
			}
    	}
    	// first derivatives
    	double Dmx = phi[1][1] - phi[0][1];
    	double Dmy = phi[1][1] - phi[1][0];
    	
    	double Dpx = phi[2][1] - phi[1][1];
    	double Dpy = phi[1][2] - phi[1][1];
    	
    	double D0x = (phi[2][1] - phi[0][1])/2.0;
    	double D0y = (phi[1][2] - phi[1][0])/2.0;
    	
    	// if using smoothing forces:
    	double smooth = 0.0;
    	if (Numerics.abs(smoothweight)>0) {
    	
			// second derivatives
			double Dxx = phi[0][1] + phi[2][1] - 2.0*phi[1][1];
			double Dyy = phi[1][0] + phi[1][2] - 2.0*phi[1][1];
			double Dxy = (phi[0][0] + phi[2][2] - phi[0][2] - phi[2][0])/4.0;
			
			// gradient norm
			double SD0x = D0x * D0x;
			double SD0y = D0y * D0y;
			double GPhi = Math.sqrt(SD0x + SD0y);
			
			// curvature
			double K =  Dyy*SD0x + Dxx*SD0y 
						- 2.0*D0x*D0y*Dxy;
				
			if(GPhi > 0.0000001){
				double tmp = GPhi*GPhi;
				K = K/(GPhi*tmp);
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
									 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0)));
			
			double DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
									 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0)));
			
			balloon = balloonweight*stepsize*(Numerics.max(balloonforces[lb][xy], 0.0)*DeltaP + Numerics.min(balloonforces[lb][xy], 0.0)*DeltaM);
		}
		
		// external force field
		double field = 0.0;
		if (Numerics.abs(fieldweight)>0) {
			field =  fieldweight*stepsize*(Numerics.max(fieldforce[X][xy],0.0)*Dmx + Numerics.min(fieldforce[X][xy],0.0)*Dpx
										  +Numerics.max(fieldforce[Y][xy],0.0)*Dmy + Numerics.min(fieldforce[Y][xy],0.0)*Dpy);
		}
		
		// pressure forces of structures 1 voxel thick
		double pressure = 0.0;
		if (phi[1][1]<0 && Numerics.abs(pressureweight)>0) {
			for (int n=0;n<4;n++) {
				int xyn = xy+xoff[n]+yoff[n];
				// neighbor-based pressure force
				pressure = Numerics.max(pressure, Numerics.min(Numerics.max(1.0-mgdmfunctions[1][xyn],0.0),
																 Numerics.max(1.0-mgdmfunctions[0][xyn],0.0)));
			}
			pressure = -pressureweight*stepsize*pressure;
		}
		
		return - (smooth - balloon - field - pressure);
    }

    /**
	 *  critical relation detection: groups objects with relations
	 */
    private final boolean homeomorphicLabeling(int xy, int lb) {
    	// if we don't check, just exit
    	if (!checkTopology) return true;
    	
		// is the new segmentation homeomorphic ? 
		
		// inside the original object ?
		if (segmentation[xy]==lb) return true;
		
		boolean [][][] obj = new boolean[3][3][3];
		
		// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (segmentation[xy+i+j*nx]==lb) {
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
			if (segmentation[xy+i+j*nx]==segmentation[xy]) {
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
			if ( (i*i+j*j>0) 
				&& (segmentation[xy+i+j*nx]!=lb) 
				&& (segmentation[xy+i+j*nx]!=segmentation[xy]) ) {
				boolean found = false;
				for (int n=0;n<Nconfiguration;n++) 
					if (segmentation[xy+i+j*nx]==lbs[n]) { found = true; break; }
				
				if (!found) {
					lbs[Nconfiguration] = segmentation[xy+i+j*nx];
					Nconfiguration++;
				}
			}
		}
		// pairs

		for (int n=0;n<Nconfiguration;n++) {
			// in relation with previous object
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if ( (segmentation[xy+i+j*nx]==segmentation[xy])
					|| (segmentation[xy+i+j*nx]==lbs[n]) ) {
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
				if ( (segmentation[xy+i+j*nx]==lb)
					|| (segmentation[xy+i+j*nx]==lbs[n]) ) {
					obj[1+i][1+j][1+l] = true;
				} else {
					obj[1+i][1+j][1+l] = false;
				}
			}
			obj[1][1][1] = true;
			if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
		}

		// triplets: may not be needed
		for (int n=0;n<Nconfiguration;n++) {
			for (int m=n+1;m<Nconfiguration;m++) {
				// in relation with previous object
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (segmentation[xy+i+j*nx]==segmentation[xy])
						|| (segmentation[xy+i+j*nx]==lbs[n])
						|| (segmentation[xy+i+j*nx]==lbs[m]) ) {
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
					if ( (segmentation[xy+i+j*nx]==lb)
						|| (segmentation[xy+i+j*nx]==lbs[n]) 
						|| (segmentation[xy+i+j*nx]==lbs[m]) ) {
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

    	float[] nbdist = new float[4];
    	boolean[] nbflag = new boolean[4];
    	boolean boundary;
    	
		for (int n=0; n<narrowband.currentsize;n++) {
			int xy = narrowband.id[n];
			
			boundary = false;
			for (int l=0; l<4; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				int xynb = xy + xoff[l] + yoff[l];
				if (mgdmlabels[0][xynb]!=mgdmlabels[0][xy] && mask[xynb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(mgdmfunctions[0][xynb]);
					nbflag[l] = true;
					boundary = true;
				}
			}
			if (boundary) {
				narrowband.functions[0][n] = isoSurfaceDistance(mgdmfunctions[0][xy], nbdist, nbflag);
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int n=0; n<narrowband.currentsize;n++) {
			int xy = narrowband.id[n];
			mgdmfunctions[0][xy] = narrowband.functions[0][n];
		}
			
        return;
    }

    /** 
    *   isosurface distance re-initialization at the boundary
    */
    private final void resetIsosurfaceBoundary() {
    	if (debug) System.out.print("fast marching evolution: iso-surface reinit\n");

    	float[] nbdist = new float[4];
    	boolean[] nbflag = new boolean[4];
    	boolean boundary;
    	
    	float[] tmp = new float[nx*ny];
    	boolean[] processed = new boolean[nx*ny];
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
			
			boundary = false;
			for (int l=0; l<4; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				
				int xynb = xy + xoff[l] + yoff[l];
				if (segmentation[xynb]!=segmentation[xy] && mask[xynb]) {
					// compute new distance based on processed neighbors for the same object
					nbdist[l] = Numerics.abs(mgdmfunctions[0][xynb]);
					nbflag[l] = true;
					boundary = true;
				}
			}
			if (boundary) {
				tmp[xy] = isoSurfaceDistance(mgdmfunctions[0][xy], nbdist, nbflag);
				processed[xy] = true;
			}
		}
		// once all the new values are computed, copy into original GDM function (sign is not important here)
		for (int xy = 0; xy<nx*ny; xy++) {
			if (processed[xy]) mgdmfunctions[0][xy] = tmp[xy];
			else mgdmfunctions[0][xy] = UNKNOWN;
		}
			
        return;
    }
    
    /** 
    *  	Evolution using a joint fast marching scheme
    */
    public final void evolveFastMarching(int iter, int sign) {
    	
    	if (debug) System.out.print("level set evolution: fast marching\n");

    	// computation variables
        boolean[] processed = new boolean[nx*ny]; // note: using a int instead of boolean for the second pass
    	
        // init
    	
    	// first estimate the binary tree size
    	int size = 0;
		for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) if (mgdmfunctions[0][xy]<=0.5f) size++;
		// create the tree with initial estimates of size
    	heap = new BinaryHeapPair(Numerics.ceil(1.25f*size), Numerics.ceil(0.1f*size), BinaryHeap2D.MINTREE);
    	
		for (int t=0;t<iter;t++) {
			if (debug) System.out.print("iteration "+t+"\n");

			
			if (debug) System.out.print("init ("+size+")\n");
        
			heap.reset();
			double curr, next;
			float speed = 0;
			int nbound = 0;
			for (int xy = 0; xy<nx*ny; xy++) if (mask[xy]) {
				// the criterion for being in the narrow band is to have a short distance to closest boundaries
				if (mgdmfunctions[0][xy]<1.0f) {
					// add to the heap
					curr = balloonweight*balloonforces[mgdmlabels[0][xy]][xy];
					next = balloonweight*balloonforces[mgdmlabels[1][xy]][xy];
					if (next>curr) {
						speed =  (float)(1.0/smoothweight);
							 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
						else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
						else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
						
						heap.addValue(speed, xy, mgdmlabels[1][xy]);
					}
				
					nbound++;
				}
				processed[xy] = false;
			}
			System.out.println("initial boundary points: "+nbound);
				
			// grow the labels and functions			
			float dist = 0.0f;
			while (heap.isNotEmpty()) {
				// extract point with minimum distance
				dist = heap.getFirst();
				int xy = heap.getFirstId1();
				int lb = heap.getFirstId2();
				heap.removeFirst();
	
				// if more than nmgdm labels have been found already, this is done
				if (processed[xy])  continue;
				
				// check topology
				if (!homeomorphicLabeling(xy, lb)) continue; 
				
				// update the segmentation and distance function (?) at the current level
				mgdmfunctions[0][xy] = dist;
				mgdmlabels[0][xy] = lb;
				processed[xy] = true;
				segmentation[xy] = lb;
				//System.out.print(".");
				
				// find new neighbors
				for (int k = 0; k<4; k++) {
					int xyn = xy + xoff[k] + yoff[k];
					
					if (mask[xyn]) {
						// must be in outside the object or its processed neighborhood
						if (!processed[xyn]) {
							if (mgdmlabels[0][xyn]!=lb) {
								curr = balloonweight*balloonforces[mgdmlabels[0][xyn]][xyn];
								next = balloonweight*balloonforces[lb][xyn];
								if (next>curr) {
									// add to the heap
									speed =  (float)(1.0/smoothweight);
										 if (sign==0) speed = (float)(1.0/(smoothweight+0.5*(next-curr)));
									else if (sign>0 && next>0) speed = (float)(1.0/(smoothweight+next));
									else if (sign<0 && curr<0) speed = (float)(1.0/(smoothweight-curr));
									
									//System.out.print("+");
									heap.addValue(dist+speed, xyn, lb);
								}
							}
						}
					}
				}
			}
			for (int xy = 0; xy<nx*ny; xy++) if (!processed[xy]) {
				mgdmfunctions[0][xy] = dist;
			}
			// re-init
			//resetIsosurfaceBoundary();
			//fastMarchingReinitialization(false);
		}
        return;
    }
    
    private final float fastMarchingForces(int xy, int lb) {
    	if (balloonforces[lb][xy]>0) return 1.0f/(float)(smoothweight+balloonweight*balloonforces[lb][xy]);
    	else return 1.0f/(float)smoothweight;
    }
}

