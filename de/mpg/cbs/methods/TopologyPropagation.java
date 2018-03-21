package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import de.mpg.cbs.structures.BinaryTree;
import de.mpg.cbs.structures.CriticalPointLUT;
import de.mpg.cbs.utilities.BasicInfo;
import de.mpg.cbs.libraries.ObjectStatistics;
 
/**
 *
 *  This algorithm propagates topology in images.
 *  <p>
 *  Different options are possible:
 *  approximate or exact, upward or downward smoothing.
 *  The propagation is based on a fast marching technique, and applies 
 *  to distance functions or grayscale images.
 *	
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class TopologyPropagation {
	
	// propagation flags, also used for checking topology
	// ! (SKELETON,OBJECT,TEST) must be > others
	public  static final byte    BACKFRONTIER = 110;
	private static final byte    BACKMASK = 	100;
	public  static final byte    OBJECT =       80;
	private static final byte    TEST =         70;
	private static final byte    BOUNDARY =     40;
	private static final byte    BACKGROUND =   30;
	private static final byte    CRITICAL =     20;
	private static final byte    MBACK =        18;
	private static final byte    MBOUND =       16;
	private static final byte    MTEST =        14;
	private static final byte    MOBJ =         12;
	private static final byte    MASK =         10;
	private static final byte    FRONTIER =     00;
	
	private static final byte	 OBJSET		=	8;
	private static final byte	 OBJTEST	=	7;
	private static final byte	 OBJBOUND	=	6;
	private static final byte	 OBJCRIT	=	5;
	private static final byte	 OBJMASK	=	3;
	
	private static final byte	 UNSET		=	4;
	
	private static final byte	 BGMASK		=	5;
	private static final byte	 BGCRIT		=	3;
	private static final byte	 BGBOUND	=	2;
	private static final byte	 BGTEST		=	1;
	private static final byte	 BGSET		=	0;
	
	// direction labeling		
	public	static	final	byte	X = 0;
	public	static	final	byte	Y = 1;
	public	static	final	byte	Z = 2;

	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-6f;
	
	// data buffers
	private 	float[][][] 	image;  			// original image
	private 	float[][][]     dist;   			// propagated function
	private 	byte[][][]      label;  			// labels for the marching progression
	private 	boolean[][][]   mask;   			// image mask: true for data points
	private 	float[][][]     geom;   			// a geometric distance term for regularisation
	private 	byte[][][]		critical;   		// a map of the critical points
	private 	BinaryTree  	boundary;   		// the binary used for the fast marching
	private static	int 		nx,ny,nz;   		// image dimensions
	
	// parameters
	private 	float   		upLevel;            // the levelset to start the upward algorithms
	private 	float   		downLevel;          // the levelset to start the downward algorithms
	private		float			maxDistance;		// the maximum change allowed from an original image
	private		boolean			thresholdStop = false;	// if we use the threshold to stop the propagation
	private     String          inputSkel;          // the type of input for the skeleton: intensity, seed point or seed region
    private     float[][]      inputSeed;          // the seed points for seed skeleton
    private     boolean[][][]   inputPaint;         // the mask for paint skeleton
	private		float			minDistance;		// the minimum distance allowed between neighboring points
    
	// computation variables
	private		int				obj;
	private		int				back;
	private		int				dz;
	
	// computation variables
	private     boolean[][][]       pattern = new boolean[3][3][3];
	private     CriticalPointLUT    lut;
	private     String	            lutdir = null;
	private     String              connectType;
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	 boolean		debug=true;
	 boolean		verbose=true;
	/**
	 *  constructor
	 */
	public TopologyPropagation(float[][][] image_, boolean [][][] mask_, 
									int nx_, int ny_, int nz_,
									float upLevel_, float downLevel_, float maxDist_, 
									float minDist_,
									String connect_, boolean tS_, 
                                    String sk_, float[][] sd_, boolean [][][] pm_) {
		image = image_;
		mask = mask_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		upLevel = upLevel_;
        downLevel = downLevel_;
		maxDistance = maxDist_;
        connectType = connect_;
		minDistance = minDist_;
		
		thresholdStop = tS_;
        inputSkel = sk_;
        inputSeed = sd_;
        inputPaint = pm_;
        
		// init all the arrays
		try {
			dist = new float[nx][ny][nz];
			geom = new float[nx][ny][nz];
			label = new byte[nx][ny][nz];
			critical = new byte[nx][ny][nz];
			boundary = new BinaryTree(nx*ny*nz, 3, BinaryTree.MINTREE);			

			// topology luts
			checkTopology=true;
			checkComposed=false;
				 if (connectType.equals("26/6")) lut = new CriticalPointLUT(lutdir, "critical266LUT.raw.gz",200);
			else if (connectType.equals("6/26")) lut = new CriticalPointLUT(lutdir, "critical626LUT.raw.gz",200);
			else if (connectType.equals("18/6")) lut = new CriticalPointLUT(lutdir, "critical186LUT.raw.gz",200);
			else if (connectType.equals("6/18")) lut = new CriticalPointLUT(lutdir, "critical618LUT.raw.gz",200);
			else if (connectType.equals("6/6")) lut = new CriticalPointLUT(lutdir, "critical66LUT.raw.gz",200);
			else if (connectType.equals("wcs")) {
				lut = new CriticalPointLUT(lutdir, "criticalWCLUT.raw.gz",200);
				checkComposed=false;
			}
			else if (connectType.equals("wco")) {
				lut = new CriticalPointLUT(lutdir, "critical66LUT.raw.gz",200);
				checkComposed=true;
			}
			else if (connectType.equals("no")) {
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
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;

		// init values
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					geom[x][y][z] = 0.0f;
					label[x][y][z] = BACKGROUND;						
					critical[x][y][z] = 0;						
					dist[x][y][z] = 0.0f;
				}
			}
		}
		if (debug) BasicInfo.displayMessage("TP:initialisation\n");
	}

	public void finalize() {
		image = null;
		dist = null;
		label = null;
		critical = null;
		geom = null;
		boundary = null;
		System.gc();
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
		label = null;
		//geom = null;
		boundary = null;
		System.gc();
	}
		
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
    public final void setUpLevel(float lv) { upLevel = lv; }
	public final void setDownLevel(float lv) { downLevel = lv; }
	public final void setSkeletonType(String sk) { inputSkel = sk; }
	public final void setImage(float [][][] img) { image = img; }
	
	public final void copyDistanceToImage() { 
		int x,y,z;
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++)
			image[x][y][z] = dist[x][y][z]; 
	}
	
	/**
	 *  higher neighborhood constraint
	 */
	private final float higherNeighborConstraint(int x, int y, int z) {
		float higher;
		int     i,j,l,xi,yj,zl;

		higher = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] == OBJECT) && (dist[xi][yj][zl] > higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  higher neighborhood constraint
	 */
	private final float higherNeighborConstraint(float val, int x, int y, int z) {
		float higher;
		int     i,j,l,xi,yj,zl;

		higher = val;
		//for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < 4) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] == OBJECT) && (dist[xi][yj][zl] > higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborConstraint(int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerNeighborBidirConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl]==OBJSET) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBackgroundConstraint(float val, int x, int y, int z) {
		float lower;
		int     i,j,l,xi,yj,zl;

		lower = val;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) 
			if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( (label[xi][yj][zl] >= OBJECT) && (dist[xi][yj][zl] < lower) )
				lower = dist[xi][yj][zl];
		}
	   
		return lower;
	}

	/**
	 *  higher neighborhood constraint
	 */
	private final float higherBoundaryConstraint(int x, int y, int z) {
		float higher,val;
		int     i,j,l,xi,yj,zl;

		val = dist[x][y][z];
		higher = INF;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < back) {
			xi=x+i; yj=y+j; zl=z+l; 
			if ( ( (label[xi][yj][zl] == BOUNDARY) || (label[xi][yj][zl] == CRITICAL) )
				&& (dist[xi][yj][zl] > val) && (dist[xi][yj][zl] < higher) )
				higher = dist[xi][yj][zl];
		}
	   
		return higher;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBoundaryConstraint(int x, int y, int z) {
		float lower,val;
		int     i,j,l,xi,yj,zl;

		lower = dist[x][y][z];
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
			if (i*i+j*j+l*l < obj) {
				xi=x+i; yj=y+j; zl=z+l; 
				//if ( (label[xi][yj][zl] == BOUNDARY) && (dist[xi][yj][zl] > lower) )
				//	lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == CRITICAL) && (dist[xi][yj][zl] < lower) )
					lower = dist[xi][yj][zl];
			}
		}
		return lower;
	}

	/**
	 *  lower neighborhood constraint
	 */
	private final float lowerBoundaryNeighbor(int x, int y, int z) {
		float lower,val;
		int     i,j,l,xi,yj,zl;

		lower = -INF;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
			if ( (i*i+j*j+l*l < 4) && (i*i+j*j+l*l > 0) ){
				xi=x+i; yj=y+j; zl=z+l; 
				if ( (label[xi][yj][zl] == BOUNDARY) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == CRITICAL) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
				if ( (label[xi][yj][zl] == BACKGROUND) && (dist[xi][yj][zl] > lower) && (dist[xi][yj][zl] <= dist[x][y][z]) )
					lower = dist[xi][yj][zl];
			}
		}
		if (lower == -INF) lower = dist[x][y][z]-0.01f;
		
		return lower;
	}

	/**
	 *  bigger neighborhood constraint
	 */
	private final int neighborCount(int x, int y, int z) {
		int   val;
		int     i,j,l,xi,yj,zl;

		val = 0;
		for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) if (i*i+j*j+l*l < obj) {
			xi=x+i; yj=y+j; zl=z+l; 
			if (label[xi][yj][zl] >= OBJECT)
				val ++;
		}
	   
		return val;
	}

	/**
	 *  regularity check
	 */
	private final boolean isRegularPoint(float[][][] img, int x, int y, int z) {
        if (!checkTopology) return true;
 		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
            if (img[x+i][y+j][z+l] > img[x][y][z]) pattern[i+1][j+1][l+1] = true;
            else if (img[x+i][y+j][z+l] < img[x][y][z]) pattern[i+1][j+1][l+1] = false;
            else pattern[i+1][j+1][l+1] = false;
        }
        if (checkComposed) if (!ObjectStatistics.isWellComposed(pattern,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(pattern,1,1,1))) return false;
		return true;
    }

	/**
	 *  regularity check
	 */
	private final boolean isRegularPoint(byte[][][] img, int x, int y, int z) {
        if (!checkTopology) return true;
 		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
            if (img[x+i][y+j][z+l]==OBJECT) pattern[i+1][j+1][l+1] = true;
            else pattern[i+1][j+1][l+1] = false;
        }
        if (checkComposed) if (!ObjectStatistics.isWellComposed(pattern,1,1,1)) return false;		
		if (!lut.get(lut.keyFromPattern(pattern,1,1,1))) return false;
		return true;
    }
 	
	/**
	 *  initialize the labels
	 */
	public final void initUpwardGeometricLabels() {
        int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( image[x][y][z] <= upLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (image[x+i][y+j][z+l] <= upLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardGeometricLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+downLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( image[x][y][z] >= downLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (image[x+i][y+j][z+l] >= downLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = -1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
    
	public final void initUpwardSmoothingLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( image[x][y][z] <= upLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
						geom[x][y][z] = 1.0f;
						dist[x][y][z] = image[x][y][z];
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (image[x+i][y+j][z+l] <= upLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initUpwardSeedLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from seed points and up
        for (int n=0;n<inputSeed.length;n++) {
            // initial object 
            int x = (int)inputSeed[n][X];
            int y = (int)inputSeed[n][Y];
            int z = (int)inputSeed[n][Z];
            label[x][y][z] = OBJECT;
            geom[x][y][z] = 1.0f;
            dist[x][y][z] = image[x][y][z];
        }
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if (label[x][y][z]==OBJECT) {
                        // do nothing
                    } else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (label[x+i][y+j][z+l] == OBJECT) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initUpwardPaintLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+upLevel+"\n");

		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = BACKFRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = BACKMASK;
					} else if ( inputPaint[x][y][z] ) {
						// initial object 
						label[x][y][z] = OBJECT;
						geom[x][y][z] = 1.0f;
						dist[x][y][z] = image[x][y][z];
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < back) {
								if ( (inputPaint[x+i][y+j][z+l]) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
						}
					}
				}
			}
		}
	}//initLabels
    
	public final void initDownwardSmoothingLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( image[x][y][z] >= downLevel ) {
						// initial object 
						label[x][y][z] = OBJECT;
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 1.0f;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (image[x+i][y+j][z+l] >= downLevel) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            //boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardSeedLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from seed points and up
        for (int n=0;n<inputSeed.length;n++) {
            // initial object 
            int x = (int)inputSeed[n][X];
            int y = (int)inputSeed[n][Y];
            int z = (int)inputSeed[n][Z];
            label[x][y][z] = OBJECT;
            geom[x][y][z] = 1.0f;
            dist[x][y][z] = image[x][y][z];
        }
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( label[x][y][z] == OBJECT ) {
						// do nothing
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (label[x+i][y+j][z+l] == OBJECT) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            //boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
	public final void initDownwardPaintLabels() {
		int[] vec = new int[3];
		if (debug) BasicInfo.displayMessage("TP:threshold "+downLevel+"\n");
		
		// init the labels from the threshold and up
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					if ( (x==0) || (y==0) || (z==0) || (x==nx-1) || (y==ny-1) || (z==nz-1) ) {
						label[x][y][z] = FRONTIER;
					} else if ( mask[x][y][z]==false ) {
						label[x][y][z] = MASK;
					} else if ( inputPaint[x][y][z] ) {
						// initial object 
						label[x][y][z] = OBJECT;
						dist[x][y][z] = image[x][y][z];
						geom[x][y][z] = 1.0f;
					} else {			
						// boundary or background ?
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
							if (i*i+j*j+l*l < obj) {
								if ( (inputPaint[x+i][y+j][z+l]) && ( mask[x+i][y+j][z+l]==true ) ) {
									label[x][y][z] = BOUNDARY;
								}					
							}
						}
						if (label[x][y][z] == BOUNDARY) {
							// boundary: stack to the tree
							dist[x][y][z] = image[x][y][z];
							geom[x][y][z] = 1.0f;
							vec[0] = x; vec[1] = y; vec[2] = z;
                            //boundary.addValue(geom[x][y][z], vec);
                            boundary.addValue(dist[x][y][z], vec);
						} else {
							label[x][y][z] = BACKGROUND;
							// for boosted version
							//dist[x][y][z] = image[x][y][z];
							dist[x][y][z] = 0.0f;
						}
					}
				}
			}
		}
	}//initLabels
			
	/**
	 *  propagate distances from a contour to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start upward geometric distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		initUpwardGeometricLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			higher = higherNeighborConstraint(x, y, z);
            
			if (higher>val) {
				val = higher;
                dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);
				critical[x][y][z] += 1;
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
               // test for update
                label[x][y][z] = TEST;

                if (isRegularPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val+1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
               }
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val+1.0f;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateUpwardGeometricDistance

	/**
	 *  propagate distances from a contour to the outside (inferior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start downward geometric distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		initDownwardGeometricLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(x, y, z);
            
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
                 vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);               
				critical[x][y][z] -= 1;
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if (isRegularPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val-1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
               }
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val-1.0f;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardGeometricDistance

	/**
	 *  propagate distances from a contour to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardUnconstrainedGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start upward unconstrained distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		initUpwardGeometricLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			higher = higherNeighborConstraint(x, y, z);
            
			if (higher>val) {
				val = higher;
                dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);
				critical[x][y][z] += 1;
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
               // test for update
                label[x][y][z] = TEST;

                //if (isRegularPoint(label,x,y,z) ) {
                if ( true ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val+1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
               }
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val+1.0f;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateUpwardGeometricDistance

	/**
	 *  propagate distances from a contour to the outside (inferior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardUnconstrainedGeometricDistance() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start downward unconstrained distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		initDownwardGeometricLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			val = dist[x][y][z];
			lower = lowerNeighborConstraint(x, y, z);
            
			if (lower<val) {
				val = lower;
                dist[x][y][z] = lower;
                 vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(dist[x][y][z], vec);               
				critical[x][y][z] -= 1;
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                //if (isRegularPoint(label,x,y,z) ) {
                if ( true ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val-1.0f;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
                    }
               } else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
               }
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val-1.0f;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardGeometricDistance

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardExactSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start upward exact distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			
			boundary.removeFirst();

			//val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if (isRegularPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                        //if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.max(image[xi][yj][zl], val+minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						//if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                        if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val+minDistance, vec);
                            //boundary.addValue(val, vec);
							//followCriticalPath(xi,yj,zl,val+minDistance);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = val;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardExactSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardExactSmoothing() {
	
		int 		t, stop;
		float   	val, dival, prec, rank;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		boolean		isCritical;
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start downward exact distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();
		
		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		dival = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			val = boundary.getFirst();
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
			if (lower<val) {
				val = lower;
				//dist[x][y][z] = lower;
				vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(lower, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary

				// test for update
				label[x][y][z] = TEST;
				
                if (isRegularPoint(label,x,y,z) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = val;

				    // search for new or critical neighbors
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l;
						if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            //boundary.addValue(val-minDistance, vec);
                            boundary.addValue(val, vec);
                        }
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = Math.min(image[xi][yj][zl], val-minDistance);
							vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(Math.min(image[xi][yj][zl], val-minDistance), vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] = 1;
				}				
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		int xb=0,yb=0,zb=0;
		float best = val;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == BACKGROUND) || (label[x][y][z] == CRITICAL) ) {
				if (thresholdStop) dist[x][y][z] = upLevel;
				else dist[x][y][z] = val;
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardExactSmoothing

    public final void updateCriticalTree(int x, int y, int z, float val) {
        // search for all critical neighbors: if none, put the point into the boundary
        boolean hasNeighbors = false;
        for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
            if ( (i*i+j*j+l*l < back) && (label[x+i][y+j][z+l] == CRITICAL) ) {
                hasNeighbors = true;
                label[x][y][z] = MASK;
                updateCriticalTree(x+i,y+j,z+l,val);
             }
        }
        hasNeighbors = false;
        if (!hasNeighbors) {    
            label[x][y][z] = BOUNDARY;
            //dist[xi][yj][zl] = val;
            int[] vec = new int[3]; 
            vec[0] = x; vec[1] = y; vec[2] = z;
            boundary.addValue(val, vec);
            vec = null;
        } else {
            // change the label back to critical
            label[x][y][z] = CRITICAL;
        }
    }
    
	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateUpwardApproxSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	higher;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
		int[] vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start upward approx distance loop\n");
		
		boundary.reset();
		boundary.setMinTree();
		if (inputSkel.equals("seed point")) initUpwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initUpwardPaintLabels();
        else initUpwardSmoothingLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;
		//time = System.currentTimeMillis();

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			boundary.removeFirst();

			//val = dist[x][y][z];
			higher = higherNeighborConstraint(val, x, y, z);
            
			if (higher>val) {
				val = higher;
                //dist[x][y][z] = higher;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val < downLevel) ) {
				// check if we're under the upward threshold if necessary
			
                // test for update
                label[x][y][z] = TEST;

                if ( (val >= image[x][y][z] + maxDistance) || (isRegularPoint(label,x,y,z) ) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = Math.min(val, image[x][y][z] + maxDistance);
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < back) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.max(image[xi][yj][zl],val+minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						// variant
						if ( (i*i+j*j+l*l < 4) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val+minDistance, vec);
                        }
                    }
				} else {
                    // critical point
                    label[x][y][z] = CRITICAL;
					critical[x][y][z] += 1;
				}
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if ( (label[x][y][z] == CRITICAL) || (label[x][y][z] == BACKGROUND) ) {
				dist[x][y][z] = Math.min(image[x][y][z]+maxDistance,val);
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

 		return;
	}//propagateUpwardApproxSmoothing

	/**
	 *  propagate a scalar field to the inside (superior values)
	 *  with a fast marching technique
	 *  and maintain the topology
	 */
	public final void propagateDownwardApproxSmoothing() {
	
		int 		t, stop;
		float   	val, prec;
		float   	lower;
		int 		xi,yj,zl;
		int 		x,y,z,k;
		long        time;
        int[]       vec = new int[3];
		
		// init: reset the boundary tree, the labels
		if (debug) BasicInfo.displayMessage("TP:start downward approx distance loop\n");
		
		boundary.reset();
		boundary.setMaxTree();
		if (inputSkel.equals("seed point")) initDownwardSeedLabels();
        else if (inputSkel.equals("paint mask")) initDownwardPaintLabels();
        else initDownwardSmoothingLabels();

		if (debug) BasicInfo.displayMessage("TP:labels initialized\n");
		
		int tmax = (nx-2)*(ny-2)*(nz-2);
		int mod = tmax/100;

		// loop on the boundary
		stop = 0;
		t = 0;
		val = 0.0f;
		while (boundary.isNotEmpty()) {
			t++;

			if ( (debug) && (t%10000==0)) BasicInfo.displayMessage(".");
			// get the next value
			x = boundary.getFirstIndex(0);
			y = boundary.getFirstIndex(1);
			z = boundary.getFirstIndex(2);
			val = boundary.getFirst();
			boundary.removeFirst();

			//val = dist[x][y][z];
			lower = lowerNeighborConstraint(val, x, y, z);
            if (lower<val) {
				// bound if necessary
				val = lower;
                //dist[x][y][z] = val;
                vec[0] = x; vec[1] = y; vec[2] = z;
                boundary.addValue(val, vec);
			} else if ( (!thresholdStop) || (val > upLevel) ) {
				// check if we're under the upward threshold if necessary
			
                label[x][y][z] = TEST;
				
                if ( (val <= image[x][y][z] - maxDistance) || (isRegularPoint(label,x,y,z) ) ) {

                    // regular point
                    label[x][y][z] = OBJECT;
                    dist[x][y][z] = Math.max(val, image[x][y][z] - maxDistance);
                    // search for new or critical neighbors
                    for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-dz;l<=dz;l++) {
                        xi=x+i; yj=y+j; zl=z+l; 
                        if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == BACKGROUND) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            dist[xi][yj][zl] = Math.min(image[xi][yj][zl],val-minDistance);
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(dist[xi][yj][zl], vec);
                        }
						// variant: keep the critical inside with lower value
						if ( (i*i+j*j+l*l < obj) && (label[xi][yj][zl] == CRITICAL) ) {
                            label[xi][yj][zl] = BOUNDARY;
                            //dist[xi][yj][zl] = val;
                            vec[0] = xi; vec[1] = yj; vec[2] = zl;
                            boundary.addValue(val-minDistance, vec);
                        }
                    }
				} else {
					// critical point
					label[x][y][z] = CRITICAL;
					critical[x][y][z] -= 1;
				}
            }
		}
		if (debug) BasicInfo.displayMessage("TP:critical points\n");

		// set the critical points
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if (label[x][y][z] == CRITICAL) {
				dist[x][y][z] = Math.max(image[x][y][z]-maxDistance, val);
			} else if (label[x][y][z] == BACKGROUND) {
				dist[x][y][z] = image[x][y][z];
			}
		}
		if (debug) BasicInfo.displayMessage("TP:end loop\n");

		return;
	}//propagateDownwardApproxSmoothing

    /** 
	 *	returns some internal functions 
	 */
	public final float[][][] exportDistance() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = dist[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances
	public final float[][][] exportGeometry() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = geom[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances
	public final float[][][] exportLabel() {
		int 	x,y,z,k;
		float[][][]	lbs = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					lbs[x][y][z] = (float)label[x][y][z];
				}
			}
		}
		return lbs;
	} // exportDistances
	public final float[][][] exportCriticalPoints() {
		int 	x,y,z,k;
		float[][][]	cpt = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					cpt[x][y][z] = (float)critical[x][y][z];
				}
			}
		}
		return cpt;
	} // exportDistances
	public final float[][][] exportImage() {
		int 	x,y,z,k;
		float[][][]	exported = new float[nx][ny][nz];
		
		for (x=0;x<nx;x++) {
			for (y=0;y<ny;y++) {
				for (z=0;z<nz;z++) {
					exported[x][y][z] = image[x][y][z];
				}
			}
		}
		return exported;
	} // exportDistances
}
