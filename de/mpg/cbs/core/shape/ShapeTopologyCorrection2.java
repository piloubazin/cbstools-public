package de.mpg.cbs.core.shape;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm uses topology-preserving fast marching to perform
 *  topology correction on various image types (rewrite of original Toads-Cruise
 *  method for python wrapping)
 *
 *	@version    Mar 2018
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class ShapeTopologyCorrection2 {
	
	// object types

	private	static	final	byte	EMPTY = -1;
	private	static	final	byte	OBJ = 1;
	private	static	final	byte	BG = 0;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    	
    // numerical quantities
	private final static float	UNKNOWN=1e15f;
	
	//private static int[] xoff;
    //private static int[] yoff;
    //private static int[] zoff;

	// data and membership buffers
	private 	float[] 		levelset;  			// level set function
	private 	float[] 		intensity;  		// intensity image
	private 	float[] 		object;   	        // binary object
	private 	int[] 			initialization;   	// initialization mask
	private 	byte[] 		    segmentation;  		// corrected image
	private 	float[] 		corrected;  		// corrected image
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz, nxyz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private     String	            lutdir = null;
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	private     byte            ngb = 6;
	
	// parameters
	private float[] inputImage;
	private int[]   startImage = null;

	private float[] correctImage;
	private int[]   correctobjImage;

	private String inputType = "binary_object";
	public static final String[] inputTypes = {"binary_object", "probability_map", "signed_distance_function"}; 
	
	private String connectType = "wcs";
	public static final String[] connectTypes = {"wcs","wco", "6/18","6/26","18/6","26/6"};
	
	private String propagType = "object->background";
	public static final String[] propagTypes = {"object->background","background->object"};
	
	private float		minDistance = 0.00001f;
	
	public static final int    OBJ2BG = 1;
	public static final int    BG2OBJ = -1;
	
	private int        propagationDir = OBJ2BG;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	

	public final void setShapeImage(float[] val) { inputImage = val; }
	public final void setShapeImageType(String val) { inputType = val; }
	
	public final void setStartingObjectImage(int[] val) { startImage = val; }
	
	public final void setPropagationDirection(String val) { propagType = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setMinimumDistance(float val) { minDistance = val; }
	
	public final void setTopology(String val) { connectType = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }
	
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Shape.devel"; }
	public final String getLabel() { return "Topology Correction"; }
	public final String getName() { return "TopologyCorrection"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Corrects the topology for a binary object, a levelset surface or a probability map (from Bazin et al., 2007)"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };

	// create outputs
	public final float[] getCorrectedImage() { return correctImage; }
    public final int[] getCorrectedObjectImage() { return correctobjImage; }
	
	public final void execute() {
		
		if (propagType.equals("object->background"))
		    propagationDir = OBJ2BG;
		else
		    propagationDir = BG2OBJ;
				
		// 6-neighborhood: pre-compute the index offsets
		//xoff = new int[]{1, -1, 0, 0, 0, 0};
		//yoff = new int[]{0, 0, nx, -nx, 0, 0};
		//zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		// init the arrays : only one should be valid here
		if (inputType.equals("binary_object")) {
		    object = inputImage;
		    intensity = null;
		    levelset = null;
		} else if (inputType.equals("probability_map")) {
		    object = null;
		    intensity = inputImage;
		    levelset = null;
		} else if (inputType.equals("signed_distance_function")) {
		    object = null;
		    intensity = null;
		    levelset = inputImage;
		}
		
		initialization = startImage;
		try {
			mask = new boolean[nx*ny*nz];
            // basic mask: remove two layers off the images (for avoiding limits)
            for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
                if (x>1 && x<nx-2 && y>1 && y<ny-2 && z>1 && z<nz-2) mask[x+nx*y+nx*ny*z] = true;
                else mask[x+nx*y+nx*ny*z] = false;
            }
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectType.equals("26/6")) { lut = new CriticalPointLUT(lutdir, "critical266LUT.raw.gz",200); ngb = 26; }
			else if (connectType.equals("6/26")) { lut = new CriticalPointLUT(lutdir, "critical626LUT.raw.gz",200); ngb = 26; }
			else if (connectType.equals("18/6")) { lut = new CriticalPointLUT(lutdir, "critical186LUT.raw.gz",200); ngb = 18; }
			else if (connectType.equals("6/18")) { lut = new CriticalPointLUT(lutdir, "critical618LUT.raw.gz",200); ngb = 18; }
			else if (connectType.equals("6/6")) { lut = new CriticalPointLUT(lutdir, "critical66LUT.raw.gz",200); ngb = 26; }
			else if (connectType.equals("wcs")) {
				lut = new CriticalPointLUT(lutdir,"criticalWCLUT.raw.gz",200);
				ngb = 26;
				//checkComposed=true;
				checkComposed=false;
			}
			else if (connectType.equals("wco")) {
				lut = null;
				ngb = 26;
				checkTopology=false;
				checkComposed=true;
			}
			else if (connectType.equals("no")) {
				lut = null;
				ngb = 26;
				checkTopology=false;
			}
			else {
			    BasicInfo.displayMessage("Specified LUT type not found: "+connectType+"\n");
				lut = null;
				ngb = 26;
				checkTopology=false;
			}
			if (checkTopology) {
				if (!lut.loadCompressedPattern()) {
					System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
					BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
				} else {
					//System.out.println("LUT loaded from: "+lut.getFilename());
				}
			}
			
			// initialize from different inputs, configurations
            segmentation = new byte[nx*ny*nz];
		    corrected = new float[nx*ny*nz];
		    if (object!=null) {
		        levelset = new float[nx*ny*nz];
		        fastMarchingObjectLevelsetIntialization();
		    }
			if (initialization==null) {
			    initialization = new int[nx*ny*nz];
			    if (intensity!=null && propagationDir==OBJ2BG) 
			        initializationFromMaximumIntensity();
			    else if (levelset!=null && propagationDir==OBJ2BG)
			        initializationFromMinimumLevelset();
			    else if (propagationDir==BG2OBJ)
			        initializationFromImageBoundaries();
			}				
		} catch (OutOfMemoryError e){
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initialization\n");
		
		if (intensity!=null) fastMarchingTopologyPropagationIntensity();
		else if (levelset!=null) fastMarchingTopologyPropagationLevelset();
		
		// output
		correctImage = corrected;
		//correctImage = levelset;
		//correctImage = intensity;
		
		correctobjImage = new int[nxyz];
		if (inputType.equals("probability_map")) {
		    for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
		        if (correctImage[xyz]>0.5) correctobjImage[xyz] = 1;
		        else correctobjImage[xyz] = 0;
		    }
		} else {
		    for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
		        if (correctImage[xyz]<0) correctobjImage[xyz] = 1;
		        else correctobjImage[xyz] = 0;
		    }
		}
    }

	private final void fastMarchingObjectLevelsetIntialization() {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			levelset[xyz] = UNKNOWN;                            
           
			// segmentation
            if (object[xyz]>0) {
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
				for (byte k = 0; k<6; k++) {
					//int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
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
        while ( heap.isNotEmpty() ) {
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
 			maxdist = dist;	// keep track of distance if stopping at the narrow band
			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				//int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (segmentation[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (byte l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							//int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
							int xyznb = Ngb.neighborIndex(l, xyzn, nx,ny,nz);
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
		// clean up the masked values
		float maxlvl = 0.0f;
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (mask[xyz] && levelset[xyz]!=UNKNOWN) {
        	    maxlvl = Numerics.max(maxlvl, levelset[xyz]);
        	}
        }
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (!mask[xyz] || levelset[xyz]==UNKNOWN) {
        	    levelset[xyz] = maxlvl;
        	}
        }
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }

     private final void initializationFromMaximumIntensity() {
		if (debug) BasicInfo.displayMessage("intialize from max intensity\n");		
         // initialize to the highest value
         float best = -1e-9f;
         int idx = -1;
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
         	if (intensity[xyz]>best) {
         	    best = intensity[xyz];
         	    idx = xyz;
         	}
         }
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
             initialization[xyz] = BG;
         }
         initialization[idx] = OBJ;
     }

     private final void initializationFromMinimumLevelset() {
         if (debug) BasicInfo.displayMessage("intialize from min levelset\n");		
         // initialize to the lowest value
         float best = 1e-9f;
         int idx = -1;
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
         	if (levelset[xyz]<best) {
         	    best = levelset[xyz];
         	    idx = xyz;
         	}
         }
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
             initialization[xyz] = BG;
         }
         initialization[idx] = OBJ;
     }
     
     private final void initializationFromImageBoundaries() {
 		if (debug) BasicInfo.displayMessage("intialize from bounding box\n");		
 		// initialize just inside the mask
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
            if (x==2 || x==nx-3 || y==2 || y==ny-3 || z==2 || z==nz-3) {
                initialization[x+nx*y+nx*ny*z] = OBJ;
            } else {
                initialization[x+nx*y+nx*ny*z] = BG;
            }
        }
     }

    private final void fastMarchingTopologyPropagationLevelset() {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			corrected[xyz] = UNKNOWN;                            
           
			// segmentation: for the topology
            if (initialization[xyz]==OBJ) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
         }
		
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz];
		//boolean[] critical = new boolean[nx*ny*nz];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching correction (levelsets)\n");	
        
        // propagation direction
        if (propagationDir==OBJ2BG) heap.setMinTree();
        else heap.setMaxTree();
        heap.reset();
		// initialize the heap from boundaries
        float mainval = 0.0f;
        if (propagationDir==OBJ2BG) mainval = -1e15f;
        else if (propagationDir==BG2OBJ) mainval = 1e15f;
        float maskval = -1e15f;
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (mask[xyz]) {
        	    if (levelset[xyz]>maskval) maskval =  levelset[xyz];
        		processed[xyz] = false;
				if (initialization[xyz]==OBJ) {
				    processed[xyz] = true;
				    corrected[xyz] = levelset[xyz];
				    segmentation[xyz] = OBJ;
				    if (propagationDir==OBJ2BG && corrected[xyz]>mainval) mainval = corrected[xyz];
				    else if (propagationDir==BG2OBJ && corrected[xyz]<mainval) mainval = corrected[xyz];
				    
				    // add neighbors to the heap
				    for (byte k = 0; k<ngb; k++) {
                        //int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
                        int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
                        
                        if (mask[xyzn] && !processed[xyzn]) {
                            heap.addValue(levelset[xyzn],xyzn,OBJ);
                        }
                    }
				}
			}
        }
          
		if (debug) BasicInfo.displayMessage("init\n");		
		
		
        // grow the labels and functions
		while ( heap.isNotEmpty() ) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float val = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			/*
			// check for neighbor value constraint
			float higher = higherNeighborConstraint(val, xyz, corrected, processed);
			float lower = lowerNeighborConstraint(val, xyz, corrected, processed);
			
			if (propagationDir==OBJ2BG && val<higher) {
			    // lower (=inside) than neighbor value -> restack
			    heap.addValue(higher, xyz, OBJ);
			} else if (propagationDir==BG2OBJ && val>lower) {
			    // higher (=outside) than neighbor value -> restack
			    heap.addValue(lower, xyz, OBJ);
			} else
			*/
			if (homeomorphicLabeling(xyz, OBJ)) {
			    // all correct: update and find new neighbors
			    
			        // update the distance functions at the current level
                corrected[xyz] = val;
                processed[xyz] = true; // update the current level
                //critical[xyz] = false;
                segmentation[xyz] = OBJ;
                mainval = val;	// keep track of distance if stopping at the narrow band
            
                // find new neighbors
                for (byte k = 0; k<ngb; k++) {
                    //int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
                    int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
                    
                    if (mask[xyzn]) {
                        // must be in outside the object or its processed neighborhood
                        if (!processed[xyzn]) {
                            // add to the heap
                            if (propagationDir==OBJ2BG) {
                                heap.addValue(Numerics.max(levelset[xyzn],val+minDistance),xyzn,OBJ);
                            } else if (propagationDir==BG2OBJ) {
                                heap.addValue(Numerics.min(levelset[xyzn],val-minDistance),xyzn,OBJ);
                            }
                        }
                    }
                }
                /*
                // for adding back the critical points only
                for (byte k = 0; k<6; k++) {
                    //int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
                    int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
                    
                    if (mask[xyzn]) {
                        // must be in outside the object or its processed neighborhood
                        if (critical[xyzn]) {
                            // add to the heap
                            if (propagationDir==OBJ2BG) {
                                heap.addValue(val,xyzn,OBJ);
                            } else if (propagationDir==BG2OBJ) {
                                heap.addValue(val,xyzn,OBJ);
                            }
                        }
                    }
                }*/
            } /*else {
                critical[xyz] = true;
            }*/
		}
		// keep it out for debug
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
		    if (!mask[xyz]) corrected[xyz] = maskval;
            else if (corrected[xyz] == UNKNOWN) corrected[xyz] = mainval;
        }
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
    }

    private final void fastMarchingTopologyPropagationIntensity() {
         // initialize the quantities
         for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm functions
			corrected[xyz] = UNKNOWN;                            
          
			// segmentation: for the topology
            if (initialization[xyz]==OBJ) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
        }
		
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching correction (intensity)\n");	
        
        // propagation direction
        if (propagationDir==OBJ2BG) heap.setMaxTree();
        else heap.setMinTree();
        heap.reset();
		// initialize the heap from boundaries
        float mainval = 0.0f;
        if (propagationDir==OBJ2BG) mainval = 1e15f;
        else if (propagationDir==BG2OBJ) mainval = -1e15f;
        float maskval = 1e15f;
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
        	if (mask[xyz]) {
        		if (intensity[xyz]<maskval) maskval =  intensity[xyz];
        		processed[xyz] = false;
				if (initialization[xyz]==OBJ) {
				    processed[xyz] = true;
				    corrected[xyz] = intensity[xyz];
				    segmentation[xyz] = OBJ;
				    if (propagationDir==OBJ2BG && corrected[xyz]<mainval) mainval = corrected[xyz];
				    else if (propagationDir==BG2OBJ && corrected[xyz]>mainval) mainval = corrected[xyz];
				     
				    // add neighbors to the heap
				    for (byte k = 0; k<ngb; k++) {
                        //int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
                        int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
                        
                        if (mask[xyzn] && !processed[xyzn]) {
                            heap.addValue(intensity[xyzn],xyzn,OBJ);
                        }
                    }
				}
			}
        }
 		if (debug) BasicInfo.displayMessage("init\n");		
		
        // grow the labels and functions
		while ( heap.isNotEmpty() ) {
        	//System.out.print(".");
        	// extract point with minimum distance
        	float val = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			/* not necessary?
			if (propagationDir==OBJ2BG && val>mainval) {
			    // lower (=inside) than current value -> restack
			    //heap.addValue(mainval-minDistance, xyz, OBJ);
			} else if (propagationDir==BG2OBJ && val<mainval) {
			    // higher (=outside) than current value -> restack
			    //heap.addValue(mainval-minDistance, xyz, OBJ);
			} else 
			*/
			if (homeomorphicLabeling(xyz, OBJ)) {
			    // all correct: update and find new neighbors
			    
			    // update the distance functions at the current level
                corrected[xyz] = val;
                processed[xyz]=true; // update the current level
                segmentation[xyz] = OBJ;
                mainval = val;	// keep track of distance if stopping at the narrow band
            
                // find new neighbors
                for (byte k = 0; k<ngb; k++) {
                    //int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
                    int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
                    
                    if (mask[xyzn]) {
                        // must be in outside the object or its processed neighborhood
                        if (!processed[xyzn]) {
                            // add to the heap
                            if (propagationDir==OBJ2BG) {
                                heap.addValue(Numerics.min(intensity[xyzn],val-minDistance),xyzn,OBJ);
                            } else if (propagationDir==BG2OBJ) {
                                heap.addValue(Numerics.max(intensity[xyzn],val+minDistance),xyzn,OBJ);
                            }
                        }
                    }
                }
            }
		}
		// keep it out for debug
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
            if (!mask[xyz]) corrected[xyz] = maskval;
            else if (corrected[xyz] == UNKNOWN) corrected[xyz] = mainval;
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
		/*	
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
		*/
		// else, it works
		return true;
    }

    private final float lowerNeighborConstraint(float val, int xyz, float[] corrected, boolean[] processed) {
        float lower = val;
        for (byte k=0;k<26;k++) {
            int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
            if (processed[xyzn]) {
                if (corrected[xyzn]<lower) lower = corrected[xyzn];
            }
        }
        return lower;
    }
    private final float higherNeighborConstraint(float val, int xyz, float[] corrected, boolean[] processed) {
        float higher = val;
        for (byte k=0;k<26;k++) {
            int xyzn = Ngb.neighborIndex(k, xyz, nx,ny,nz);
            if (processed[xyzn]) {
                if (corrected[xyzn]>higher) higher = corrected[xyzn];
            }
        }
        return higher;
    }
}

