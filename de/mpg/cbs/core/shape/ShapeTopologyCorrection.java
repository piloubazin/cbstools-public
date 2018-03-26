package de.mpg.cbs.core.shape;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.io.*;
import java.util.*;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class ShapeTopologyCorrection {

	// jist containers
	private float[] inputImage;
	private int[]   startImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	// for cropping
	private int mx, my, mz;
	private int         x0,xN,y0,yN,z0,zN;
	private float		Imin, Imax;
	
	private String inputType = "binary_object";
	public static final String[] inputTypes = {"binary_object", "probability_map", "signed_distance_function"}; 
	
	private String connectType = "object->background";
	public static final String[] connectTypes = {"wsc","6/18","6/26","18/6","26/6"};
	
	private String propagType = "object->background";
	public static final String[] propagTypes = {"object->background","background->object"};
	
	private float   	highestLevel = 1.0f;
    private float   	lowestLevel  = 0.0f;
	private float       minDistance = 0.0001f;
	private boolean     useMinMax = true;
    private	boolean		thresholdStop = false;
	
	// information on the parts of the image to process
	private boolean 	cropBackground = false;
    private float       signalThreshold=0.0f;
    private int         nxmin,nymin,nzmin,nxmax,nymax,nzmax;

	private float[]     correctImage;
	private int[]       correctobjImage;
	
	// numerical quantities
	private static final	float	INVSQRT2 = (float)(1.0/FastMath.sqrt(2.0));
	private static final	float	INVSQRT3 = (float)(1.0/FastMath.sqrt(3.0));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0);

	// direction labeling		
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;

	// feature choice labeling		
	public	static	final	byte	DIST = 100;
	public	static	final	byte	PROBA = 101;
	public	static	final	byte	BIN = 102;
	
	// computation variables
	private CriticalPointLUT lut;
	private BinaryHeap3D	heap;
	private String	lutdir = null;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;

	// create inputs
	public final void setShapeImage(float[] val) { inputImage = val; }
	public final void setShapeImageType(String val) { inputType = val; }
	
	public final void setStartingObjectImage(int[] val) { startImage = val; }
	
	public final void setPropagationDirection(String val) { propagType = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setLowerThreshold(float val) { lowestLevel = val; }
	public final void setHigherThreshold(float val) { highestLevel = val; }
	public final void setNormalizeIntensityRange(boolean val) { useMinMax = val; }
	
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

				
    /**
    *	computes the topology correction for 3D images
    */
    public final void execute() {    
		float[][][]   	result=null; 
        float[][][]     image;
		boolean[][][]     objectMask;
		boolean[][][]     obj = null;
		int[][][]     label = null;
		int 		x,y,z,k,m,i,j,l;
		float   	num,den;
		float   	ngb, neighbors;
		float   	dist;
        float[][][]   	ext = null; 
        boolean[][][]     paint;
		TopologyPropagation algorithm = null;
		String info;
		float[][][] lb = null;
		float level,Dmax;
		int maxLabels,Nc;
		float[] labels;
		boolean isNew;
		float best;
        int xM,yM,zM;
		boolean stop;
                
			            
		int mod;
		int progress;

		long start_time, inner_loop_time;

		if (verbose) BasicInfo.displayMessage("start correction\n");
				
		if (verbose) BasicInfo.displayMessage("extract data..\n");
		
		/* pre-processing : expand boundaries, so that we don't have to worry for them */
        try {
            expandSize();
            image = expandBoundaries(inputImage);
            if (startImage!=null) {
               paint = expandBoundaries(startImage);
               if (verbose) BasicInfo.displayMessage("use input mask");
			}
            else paint = null;
            // dispose of imgs
            inputImage = null;
            startImage = null;
        } catch (OutOfMemoryError e) {
            BasicInfo.displayMessage("Topology Correction: Out of memory creating process buffer");
            return;
        }
		if (verbose) BasicInfo.displayMessage("set parameters..\n");

		// compute exact min and max (calcMinMax() returns integers)
		Imin = 1e30f;
		Imax = -1e30f;
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			if (image[x][y][z] < Imin) Imin = image[x][y][z];
			if (image[x][y][z] > Imax) Imax = image[x][y][z];
		}
		if (verbose) BasicInfo.displayMessage("min: "+Imin+", max:"+Imax+"\n");
		
		// use relative values for the levels
		if (useMinMax) {
			highestLevel = Imax - (1.0f-highestLevel)*(Imax-Imin);
			lowestLevel  = Imin + lowestLevel*(Imax-Imin);
		}
		if (verbose) BasicInfo.displayMessage("high: "+highestLevel+", low:"+lowestLevel+"\n");
		
		if (verbose) BasicInfo.displayMessage("create mask..\n");
		
		// process only parts where intensity > lowest 
		// (shrink the image, but doesn't use the mask for topology purposes)
		objectMask = createObjectMask(image, lowestLevel);
        
		// create smaller images for speed and memory
		if (verbose) BasicInfo.displayMessage("compute boundaries..\n");
		computeProcessingBoundaries(objectMask);
		
		if (verbose) BasicInfo.displayMessage("shrink images..\n");
		image = reduceImageSize(image); 
		objectMask = reduceImageSize(objectMask);
		
		objectMask = resetCroppedMask();
        
		if (paint!=null) paint = reduceImageSize(paint);
        
		if (verbose) BasicInfo.displayMessage("new dimensions: "+mx+"x"+my+"x"+mz+"\n");
		
		/////////////////////////////////////////
		//          MAIN ALGORITHM             //
		/////////////////////////////////////////
		if (verbose) BasicInfo.displayMessage("start algorithm\n");
 
		// record time
		 start_time = System.currentTimeMillis();
		        
		System.out.println("Start ALGORITHM");
		// correction of the topology for a scalar image			
        if (inputType.equals("probability_map") || inputType.equals("signed_distance_function")) {
			// use relative values for the levels
			if (useMinMax) {
				minDistance = minDistance*(Imax-Imin);
			}
			// distance functions: flip the sign
			if (inputType.equals("signed_distance_function")) {
			    for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
			        image[x][y][z] *= -1.0f;
			    }
			    float tmp = highestLevel;
			    highestLevel = lowestLevel;
			    lowestLevel = tmp;
			}
			// look for the starting point
			if (paint==null) {
				if (propagType.equals("object->background")) {
					// create a label image for the object (value >= middleLevel)
					obj = ObjectExtraction.objectFromImage(image, mx, my, mz,
															0.5f*(lowestLevel+highestLevel), 
															highestLevel, ObjectLabeling.SUPEQUAL, ObjectLabeling.INFEQUAL);
					
					label = ObjectLabeling.connected6Object3D(obj, mx, my, mz);
					
					// find the largest label
					obj = ObjectLabeling.largestObjectFromLabel(label, mx,my,mz);
					// set the starting point for propagation
					// if multiple points, find the closest to the center of mass
					best = lowestLevel;
					xM=0;yM=0;zM=0;
					int nbest=1;
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						if ( (obj[x][y][z]) && (image[x][y][z]>best) ) {
							best = image[x][y][z];
							xM=x; yM=y; zM=z;
							nbest=1;
						} else if ( (obj[x][y][z]) && (image[x][y][z]==best) ) {
							xM+=x; yM+=y; zM+=z;
							nbest++;
						}
					}
					System.out.println("best val: "+best+", number: "+nbest);
					if (nbest>1) {
						xM /= nbest;
						yM /= nbest;
						zM /= nbest;
						System.out.println("best center: ("+xM+","+yM+","+zM+")");
					
						if ( (!obj[xM][yM][zM]) || (image[xM][yM][zM]!=best) ) {
							// find closest candidate to center of mass
							int xB=0; int yB=0; int zB=0;
							float dbest = mx*mx+my*my+mz*mz;
							for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
								if ( (obj[x][y][z]) && (image[x][y][z]==best) ) {
									if ((x-xM)*(x-xM)+(y-yM)*(y-yM)+(z-zM)*(z-zM)<dbest) { 
										xB=x; yB=y; zB=z;
										dbest = (x-xM)*(x-xM)+(y-yM)*(y-yM)+(z-zM)*(z-zM);
									}
								}
							}
							xM = xB; yM = yB; zM = zB;
							System.out.println("best closest to center: ("+xM+","+yM+","+zM+")");
						}
					}
					
					// create a paint mask with just this point
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					paint[xM][yM][zM] = true;
					if (verbose) BasicInfo.displayMessage("starting point: "+image[xM][yM][zM]+"("+xM+","+yM+","+zM+"\n");
				} else if (propagType.equals("background->object")) {
					// create a paint mask with the outer boundary
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					for (x=0;x<mx;x++) for (y=0;y<my;y++) {
						paint[x][y][0] = true;
						paint[x][y][1] = true;
						paint[x][y][mz-1] = true;
						paint[x][y][mz-2] = true;
					}
					for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[0][y][z] = true;
						paint[1][y][z] = true;
						paint[mx-1][y][z] = true;
						paint[mx-2][y][z] = true;
					}
					for (x=0;x<mx;x++) for (z=0;z<mz;z++) {
						paint[x][0][z] = true;
						paint[x][1][z] = true;
						paint[x][my-1][z] = true;
						paint[x][my-2][z] = true;
					}
					if (verbose) BasicInfo.displayMessage("starting from a bounding box\n");
				}
            } else {
				// set all points in the mask to highestLevel
				for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
					if (paint[x][y][z]) image[x][y][z] = highestLevel;
				}
			}
			// flatten the values above and below the levels
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z] > highestLevel) image[x][y][z] = lowestLevel;
				if (image[x][y][z] < lowestLevel)  image[x][y][z] = lowestLevel;
			}
			//System.out.println("TOPOLOGY PROPOGATION");
			// start the algorithm
			if (propagType.equals("background->object")) 
				algorithm = new TopologyPropagation(image, objectMask, mx, my, mz, lowestLevel, highestLevel, 
													0.0f, minDistance, connectType, thresholdStop, 
													"paint mask", null, paint);
			else if (propagType.equals("object->background")) 
				algorithm = new TopologyPropagation(image, objectMask, mx, my, mz, lowestLevel, highestLevel, 
													0.0f, minDistance, connectType, thresholdStop, 
													"paint mask", null, paint);
			
			//System.out.println("FINISHED PROPOGATION");
			if (verbose) BasicInfo.displayMessage("initialisation ("+algorithm.isWorking()+")\n");
												
			// initialize from the paint mask
			// redundant ? algorithm.initDownwardPaintLabels();
			
			if (propagType.equals("background->object")) 
				algorithm.propagateUpwardExactSmoothing();
			else if (propagType.equals("object->background")) 
				algorithm.propagateDownwardExactSmoothing();
			
			result = extendImageSize(algorithm.exportDistance());
		
			algorithm.finalize();
			
		} else if (inputType.equals("binary_object")) {
			// extract the distance function or the object
			obj = ObjectExtraction.objectFromImage(image,mx,my,mz,lowestLevel,highestLevel,ObjectLabeling.SUPERIOR,ObjectLabeling.INFEQUAL);
			ext = new float[mx][my][mz];
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (obj[x][y][z]) ext[x][y][z] = 1.0f;
				else ext[x][y][z] = 0.0f;
			}
			// start the algorithm	
			if (propagType.equals("background->object")) {
				algorithm = new TopologyPropagation(ext, objectMask, mx, my, mz, 
											0.0f, 0.5f, 0.0f, 0.0f,
											connectType, false, 
											"intensity", null, null);
	 
				// compute distance function
				algorithm.propagateDownwardUnconstrainedGeometricDistance();
				image = algorithm.exportDistance();
				algorithm.finalize();			
			} else if (propagType.equals("object->background")) {
				algorithm = new TopologyPropagation(ext, objectMask, mx, my, mz, 
											0.5f, 1.0f, 0.0f, 0.0f,
											connectType, false, 
											"intensity", null, null);
	 
				// compute distance function
				algorithm.propagateUpwardUnconstrainedGeometricDistance();
				image = algorithm.exportDistance();
				algorithm.finalize();
			}
			// reset highest and lowest to fit the distance function
			highestLevel = 0.0f;
			lowestLevel = 0.0f;
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z]>highestLevel) highestLevel=image[x][y][z];
				if (image[x][y][z]<lowestLevel)  lowestLevel=image[x][y][z];
			}		
			// use relative values for the levels
			if (useMinMax) {
				minDistance = minDistance*(Imax-Imin);
			}
					
			// find the highest point to start from
            if (paint==null) {
				if (propagType.equals("object->background")) {
					// set the starting point for propagation
					best = 0.0f;
					xM=0;yM=0;zM=0;
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						if (image[x][y][z]>best) {
							best = image[x][y][z];
							xM=x; yM=y; zM=z;
						}
					}
					// create a paint mask with just this point
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					paint[xM][yM][zM] = true;
					if (verbose) BasicInfo.displayMessage("starting point: "+image[xM][yM][zM]+"("+xM+","+yM+","+zM+"\n");
				} else if (propagType.equals("background->object")) {
					// create a paint mask with the outer boundary
					paint = new boolean[mx][my][mz];
					for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[x][y][z] = false;
					}
					for (x=0;x<mx;x++) for (y=0;y<my;y++) {
						paint[x][y][0] = true;
						paint[x][y][1] = true;
						paint[x][y][mz-1] = true;
						paint[x][y][mz-2] = true;
					}
					for (y=0;y<my;y++) for (z=0;z<mz;z++) {
						paint[0][y][z] = true;
						paint[1][y][z] = true;
						paint[mx-1][y][z] = true;
						paint[mx-2][y][z] = true;
					}
					for (z=0;z<mz;z++) for (x=0;x<mx;x++) {
						paint[x][0][z] = true;
						paint[x][1][z] = true;
						paint[x][my-1][z] = true;
						paint[x][my-2][z] = true;
					}
					if (verbose) BasicInfo.displayMessage("starting from a bounding box\n");
										
				}
            } else {
				// set all points in the mask to highestLevel
				for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
					if (paint[x][y][z]) image[x][y][z] = highestLevel;
				}
			}
			// flatten the values above and below the levels
			for (x=0;x<mx;x++) for (y=0;y<my;y++) for (z=0;z<mz;z++) {
				if (image[x][y][z] > highestLevel) image[x][y][z] = lowestLevel;
				if (image[x][y][z] < lowestLevel)  image[x][y][z] = lowestLevel;
			}
			
			// start the algorithm
			algorithm = new TopologyPropagation(image, objectMask, mx, my, mz, lowestLevel, highestLevel, 
												0.0f, minDistance, connectType, true, 
												"paint mask", null, paint);

			if (verbose) BasicInfo.displayMessage("initialisation\n");
												
			// initialize from the paint mask
			// redundant ? algorithm.initDownwardPaintLabels();
			
			if (propagType.equals("background->object")) 
				algorithm.propagateUpwardExactSmoothing();
			else if (propagType.equals("object->background")) 
				algorithm.propagateDownwardExactSmoothing();
			
			result = extendImageSize(algorithm.exportDistance());
		
			algorithm.finalize();	
		}					

		// debug
		if (verbose) BasicInfo.displayMessage("total time: (milliseconds): " + (System.currentTimeMillis()-start_time)); 

		// extract results (segmentation and/or classes) and store them in destImage[]
		nx = nx-2;
		ny = ny-2;
		nz = nz-2;
		try {
			correctImage = new float[nx*ny*nz];
			correctobjImage = new int[nx*ny*nz];
			for (x=0;x<nx;x++) {
				for (y=0;y<ny;y++) {
					for (z=0;z<nz;z++) {
						correctImage[x + nx*y + nx*ny*z] = result[x+1][y+1][z+1];
						if (inputType.equals("binary_object")) {
						    correctImage[x + nx*y + nx*ny*z] = result[x+1][y+1][z+1];
						    if (result[x+1][y+1][z+1]>0) correctobjImage[x + nx*y + nx*ny*z] = 1;
						    else correctobjImage[x + nx*y + nx*ny*z] = 0;
						}
						else if (inputType.equals("probability_map")) {
						    correctImage[x + nx*y + nx*ny*z] = result[x+1][y+1][z+1];
						    if (result[x+1][y+1][z+1]>0.5) correctobjImage[x + nx*y + nx*ny*z] = 1;
						    else correctobjImage[x + nx*y + nx*ny*z] = 0;						    
						}
						else if (inputType.equals("signed_distance_function")) {
						    correctImage[x + nx*y + nx*ny*z] = -result[x+1][y+1][z+1];
						    if (result[x+1][y+1][z+1]>0) correctobjImage[x + nx*y + nx*ny*z] = 1;
						    else correctobjImage[x + nx*y + nx*ny*z] = 0;						    
						}
					}
				}
			}
		} catch (OutOfMemoryError e) {
			BasicInfo.displayMessage("Algorithm Topology: Out of memory creating outputs");
			return;
		}

        image = null;
		objectMask = null;
		result = null;
		lb = null;
		obj = null;
		label = null;
		System.out.println("COMPLETED");
    }
	
	/** expand boundaries for spatial comutations */
	private void expandSize() {
		nx = nx+2;
		ny = ny+2;
		nz = nz+2;
	}			
	private float[][][] expandBoundaries(float[] image) {
		int 		x,y,z;
		float[][][] 	tmp;
		
		tmp = new float[nx][ny][nz];
		for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++)
					tmp[x][y][z] = image[ (x-1) + (nx-2)*(y-1) + (nx-2)*(ny-2)*(z-1) ];
		
		return tmp;
	}
	private boolean[][][] expandBoundaries(int[] image) {
		int 		x,y,z;
		boolean[][][] 	tmp;
		
		tmp = new boolean[nx][ny][nz];
		for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++)
					tmp[x][y][z] = (image[ (x-1) + (nx-2)*(y-1) + (nx-2)*(ny-2)*(z-1) ]>0);
		
		return tmp;
	}
	
	/** creates a mask for unused data */
	private boolean[][][] createObjectMask(float[][][] image, float val) {
		int 		x,y,z;
		boolean[][][]  	objMask;

		// uses only values over the threshold, if mask used
		objMask = new boolean[nx][ny][nz];
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++)
				for (z=0;z<nz;z++) {
                    if ( (cropBackground) && (image[x][y][z] <= val) )
                        objMask[x][y][z] = false;
                    else
                        objMask[x][y][z] = true;                    
                }
		// remove the boundary from the computations
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++) {
				objMask[x][y][0] = false;
				objMask[x][y][nz-1] = false;
			}
		for (y=0;y<ny;y++)
			for (z=0;z<nz;z++) {
				objMask[0][y][z] = false;
				objMask[nx-1][y][z] = false;
			}
		for (z=0;z<nz;z++)
			for (x=0;x<nx;x++) {
				objMask[x][0][z] = false;
				objMask[x][ny-1][z] = false;
			}

		return objMask;
	} // createObjectMask
 	/** create smaller image (for saving memory) */
	private boolean[][][] resetCroppedMask() {
		boolean[][][] objMask = new boolean[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			objMask[x][y][z] = true;
		}
		for (int x=0;x<mx;x++)
			for (int y=0;y<my;y++) {
				objMask[x][y][0] = false;
				objMask[x][y][mz-1] = false;
			}
		for (int y=0;y<my;y++)
			for (int z=0;z<mz;z++) {
				objMask[0][y][z] = false;
				objMask[mx-1][y][z] = false;
			}
		for (int z=0;z<mz;z++)
			for (int x=0;x<mx;x++) {
				objMask[x][0][z] = false;
				objMask[x][my-1][z] = false;
			}
		return objMask;
	}
       
    /** sets the processing lower and upper boundaries */
    private void computeProcessingBoundaries(boolean[][][] objMask) {
		int 		x,y,z;
        
        x0 = nx;
        xN = 0;
        y0 = ny;
        yN = 0;
        z0 = nz;
        zN = 0;
        for (x=1;x<nx-1;x++)
			for (y=1;y<ny-1;y++)
				for (z=1;z<nz-1;z++) {
                    if (objMask[x][y][z]) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
		// update the smaller size parameters (include am extended boundary)
		mx = xN - x0 + 1 + 2;
		my = yN - y0 + 1 + 2;
		mz = zN - z0 + 1 + 2;
		
        // debug
        System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
	/** create smaller image (for saving memory) */
	private float[][][] reduceImageSize(float[][][] image) {
		float[][][] smaller = new float[mx][my][mz];

		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					smaller[x-x0+1][y-y0+1][z-z0+1] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
	private boolean[][][] reduceImageSize(boolean[][][] image) {
		boolean[][][] smaller = new boolean[mx][my][mz];

		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					smaller[x-x0+1][y-y0+1][z-z0+1] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
	private void reducePointPosition(float[][] pt) {
		for (int n=0;n<pt.length;n++) {
            // from image to image+boundary to reduced image
            pt[n][X] = pt[n][X] - x0 + 1 + 2;
            pt[n][Y] = pt[n][Y] - y0 + 1 + 2;
            pt[n][Z] = pt[n][Z] - z0 + 1 + 2;
		}
		return;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(float[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(int[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = (float)image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	private float[][][] extendImageSize(byte[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
					larger[x][y][z] = (float)image[x-x0+1][y-y0+1][z-z0+1];
				}
			}
		}
		return larger;
	}
    /** retrieve original size from smaller image */
	private float[][][] extendImageSize(boolean[][][] image) {
		float[][][] larger = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					larger[x][y][z] = 0.0f;
				}
			}
		}		
 		for (int x=x0-1;x<=xN+1;x++) {
            for (int y=y0-1;y<=yN+1;y++) {
                for (int z=z0-1;z<=zN+1;z++) {
                    if (image[x-x0+1][y-y0+1][z-z0+1])
                        larger[x][y][z] = 1.0f;
                    else
                        larger[x][y][z] = 0.0f;
				}
			}
		}
		return larger;
	}
}
