package de.mpg.cbs.core.shape;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class ShapeSimpleSkeleton {

	// jist containers
	private float[] inputImage=null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private float boundParam;
	private float distParam;
	private String featureParam = "probability_map";
	private static final String[] featureTypes = {"signed_distance","probability_map"}; 
	
	private byte[] medialImage;
	private byte[] skelImage;
	
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
	
	// computation variables
	private boolean[][][] obj = new boolean[3][3][3];
	private CriticalPointLUT lut;
	private BinaryHeap3D	heap3;
	private BinaryHeapPair	heap2;
	private String	lutdir = null;
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;

	// create inputs
	public final void setShapeImage(float[] val) { inputImage = val; }
	public final void setShapeImageType(String val) { featureParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setBoundaryThreshold(float val) {boundParam = val; }
	public final void setSkeletonThreshold(float val) {distParam = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }
	public final void setFeatureType(String val) {featureParam = val; }
			
	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Shape.devel"; }
	public final String getLabel() { return "Simple Skeleton"; }
	public final String getName() { return "SimpleSkeleton"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Create a skeleton for a levelset surface or a probability map (loosely adapted from Bouix et al., 2006)"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };

	// create outputs
	public final byte[] getMedialSurfaceImage() { return medialImage; }
	public final byte[] getMedialCurveImage() { return skelImage; }

	public void execute3D(){
		
		// import the image data into 3D arrays
		float[][][] input = new float[nx][ny][nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			input[x][y][z] = inputImage[xyz];
		}
		inputImage = null;
		
		byte mode;
		if (featureParam.equals("signed_distance")) mode = DIST;
		else if (featureParam.equals("probability_map")) mode = PROBA;
		else mode = 0;
		
		float boundary = boundParam;
		float distscale = distParam;
		
		int objsize=0;
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			// just flip sign on distance functions
			if (mode==DIST) input[x][y][z] *= -1.0f;
			// could be optimized (e.g. 1.5x boundary size)
			if (input[x][y][z]>boundary) objsize++;
		}
				
		// 1. get medial skeleton
		lut = new CriticalPointLUT(lutdir,"critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
			return;
		} else {
			//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
		}
		heap3 = new BinaryHeap3D(objsize, BinaryHeap3D.MINTREE);
		
		byte[][][] medial = medialThinning3D(input, nx, ny, nz, boundary, distscale);
		
		// 2. get curve skeleton
		byte[][][] skel = medialCurveThinning3D(input, nx, ny, nz, boundary, distscale);
		
		medialImage = new byte[nxyz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			medialImage[xyz] = medial[x][y][z];
		}
		medial=null;
		
		skelImage = new byte[nxyz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			skelImage[xyz] = skel[x][y][z];
		}
		skel=null;
	}
	
	public void execute2D(){
		
		// import the image data into 2D arrays
		float[][] input = new float[nx][ny];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) {
			int xyz = x+nx*y;
			input[x][y] = inputImage[xyz];
		}
		inputImage = null;
		
		byte mode;
		if (featureParam.equals("signed_distance")) mode = DIST;
		else if (featureParam.equals("probability_map")) mode = PROBA;
		else mode = 0;
		
		float boundary = boundParam;
		float distscale = distParam;
		
		int objsize=0;
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) {
			// just flip sign on distance functions
			if (mode==DIST) input[x][y] *= -1.0f;
			// could be optimized (e.g. 1.5x boundary size)
			if (input[x][y]>boundary) objsize++;
		}
				
		// 1. get medial skeleton
		lut = new CriticalPointLUT(lutdir,"critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			BasicInfo.displayMessage("Problem loading the algorithm's LUT from: "+lut.getFilename()+"\n");
			return;
		} else {
			//if (debug) System.out.println("LUT loaded from: "+lut.getFilename());
		}
		heap2 = new BinaryHeapPair(objsize, BinaryHeap2D.MINTREE);
		
		byte[][] medial = medialThinning2D(input, nx, ny, boundary, distscale);
		
		// 2. get curve skeleton
		byte[][] skel = medialCurveThinning2D(input, nx, ny, boundary, distscale);
		
		medialImage = new byte[nxyz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) {
			int xyz = x+nx*y;
			medialImage[xyz] = medial[x][y];
		}
		medial=null;
		
		skelImage = new byte[nxyz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) {
			int xyz = x+nx*y;
			skelImage[xyz] = skel[x][y];
		}
		skel=null;
	}
	
	public final byte[][][] medialThinning3D(float[][][] funct, int nx, int ny, int nz, float startthreshold, float distthreshold) {
		
		// from the boundary, find the middle (approx)
		if (verbose) BasicInfo.displayMessage("fast marching init\n");		
		heap3.reset();
		boolean[][][] processed = new boolean[nx][ny][nz];
		boolean[][][] medial = new boolean[nx][ny][nz];
		byte[][][] skeleton = new byte[nx][ny][nz];
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) {
			if (funct[x][y][z]>startthreshold) {
				skeleton[x][y][z] = 1;
				processed[x][y][z] = false;
				medial[x][y][z] = false;
					
				// search for boundaries
				boolean isboundary=false;
				for (int k = 0; k<6 && !isboundary; k++) {
					if (funct[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]<=startthreshold) {
						isboundary=true;
					}
				}
				// add to the heap
				if (isboundary) {
					heap3.addValue(funct[x][y][z], x, y, z);
					if (debug) BasicInfo.displayMessage(".");
				}
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
        	// thin the object
		while (heap3.isNotEmpty()) {
			// extract point with minimum distance
			float dist = heap3.getFirst();
			int x = heap3.getFirstX();
			int y = heap3.getFirstY();
			int z = heap3.getFirstZ();
			heap3.removeFirst();

			// if already processed, pass
			if (processed[x][y][z])  continue;
			
			// check for simple point
			if (isSimplePoint3D(skeleton, x, y, z)) {
				if (isMedialEndpoint3D(skeleton, x, y, z) && funct[x][y][z]>distthreshold) {
					// label as medial surface endpoint
					medial[x][y][z] = true;
					if (debug) BasicInfo.displayMessage("-");
				} else {
					if (debug) BasicInfo.displayMessage("x");
					skeleton[x][y][z] = 0;
					
					// find new neghbors: 6C
					for (int k = 0; k<6; k++) {
						if (skeleton[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]==1 && !processed[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]) {
							// add to the heap
							heap3.addValue(funct[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]], x+Ngb.x[k], y+Ngb.y[k], z+Ngb.z[k]);
						}
					}
				}
						
				// discard point
				processed[x][y][z]=true; 
			} else {
				if (debug) BasicInfo.displayMessage("!");	
			}
		}
		// label the medial endpoints
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) if (medial[x][y][z]) {
			skeleton[x][y][z] = 2;
		}
		
		if (verbose) BasicInfo.displayMessage("done\n");		
		
		return skeleton;
     }
     
	public final byte[][] medialThinning2D(float[][] funct, int nx, int ny, float startthreshold, float distthreshold) {
		
		// from the boundary, find the middle (approx)
		if (verbose) BasicInfo.displayMessage("fast marching init\n");		
		heap2.reset();
		boolean[][] processed = new boolean[nx][ny];
		boolean[][] medial = new boolean[nx][ny];
		byte[][] skeleton = new byte[nx][ny];
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) {
			if (funct[x][y]>startthreshold) {
				skeleton[x][y] = 1;
				processed[x][y] = false;
				medial[x][y] = false;
					
				// search for boundaries
				boolean isboundary=false;
				for (int k = 0; k<4 && !isboundary; k++) {
					if (funct[x+Ngb2.x[k]][y+Ngb2.y[k]]<=startthreshold) {
						isboundary=true;
					}
				}
				// add to the heap
				if (isboundary) {
					heap2.addValue(funct[x][y], x, y);
					if (debug) BasicInfo.displayMessage(".");
				}
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
        	// thin the object
		while (heap2.isNotEmpty()) {
			// extract point with minimum distance
			float dist = heap2.getFirst();
			int x = heap2.getFirstId1();
			int y = heap2.getFirstId2();
			heap2.removeFirst();

			// if already processed, pass
			if (processed[x][y])  continue;
			
			// check for simple point
			if (isSimplePoint2D(skeleton, x, y)) {
				if (isMedialEndpoint2D(skeleton, x, y) && funct[x][y]>distthreshold) {
					// label as medial surface endpoint
					medial[x][y] = true;
					if (debug) BasicInfo.displayMessage("-");
				} else {
					if (debug) BasicInfo.displayMessage("x");
					skeleton[x][y] = 0;
					
					// find new neghbors: 6C
					for (int k = 0; k<4; k++) {
						if (skeleton[x+Ngb2.x[k]][y+Ngb2.y[k]]==1 && !processed[x+Ngb2.x[k]][y+Ngb2.y[k]]) {
							// add to the heap
							heap2.addValue(funct[x+Ngb2.x[k]][y+Ngb2.y[k]], x+Ngb2.x[k], y+Ngb2.y[k]);
						}
					}
				}
						
				// discard point
				processed[x][y]=true; 
			} else {
				if (debug) BasicInfo.displayMessage("!");	
			}
		}
		// label the medial endpoints
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) if (medial[x][y]) {
			skeleton[x][y] = 2;
		}
		
		if (verbose) BasicInfo.displayMessage("done\n");		
		
		return skeleton;
     }
     
	public final byte[][][] medialCurveThinning3D(float[][][] funct, int nx, int ny, int nz, float startthreshold, float distthreshold) {
		
		// from the boundary, find the middle (approx)
		if (verbose) BasicInfo.displayMessage("fast marching init\n");		
		heap3.reset();
		boolean[][][] processed = new boolean[nx][ny][nz];
		boolean[][][] medial = new boolean[nx][ny][nz];
		byte[][][] skeleton = new byte[nx][ny][nz];
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) if (funct[x][y][z]>startthreshold) {
        		skeleton[x][y][z] = 1;
			processed[x][y][z] = false;
			medial[x][y][z] = false;
        		
			// search for boundaries
			boolean isboundary=false;
			for (int k = 0; k<6 && !isboundary; k++) {
				if (funct[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]<=startthreshold) {
					isboundary=true;
				}
			}
			// add to the heap: use the distance to boundary as driving force
			if (isboundary) {
				heap3.addValue(funct[x][y][z], x, y, z);
				if (debug) BasicInfo.displayMessage(".");
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
		// thin the object
		while (heap3.isNotEmpty()) {
			// extract point with minimum distance
			float dist = heap3.getFirst();
			int x = heap3.getFirstX();
			int y = heap3.getFirstY();
			int z = heap3.getFirstZ();
			heap3.removeFirst();

			// if already processed, pass
			if (processed[x][y][z])  continue;
			
			// check for simple point
			if (isSimplePoint3D(skeleton, x, y, z)) {
				if (isCurveEndpoint3D(skeleton, x, y, z) && funct[x][y][z]>distthreshold) {
					// label as medial surface endpoint
					medial[x][y][z] = true;
					if (debug) BasicInfo.displayMessage("-");
				} else {
					if (debug) BasicInfo.displayMessage("x");
					skeleton[x][y][z] = 0;
					
					// find new neghbors: 6C
					for (int k = 0; k<6; k++) {
						if (skeleton[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]==1 && !processed[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]]) {
							// add to the heap
							heap3.addValue(funct[x+Ngb.x[k]][y+Ngb.y[k]][z+Ngb.z[k]], x+Ngb.x[k], y+Ngb.y[k], z+Ngb.z[k]);
						}
					}
				}
						
				// discard point
				processed[x][y][z]=true; 
			} else {
				if (debug) BasicInfo.displayMessage("!");	
			}
		}
		// label the medial endpoints
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) for (short z=1;z<nz-1;z++) if (medial[x][y][z]) {
			skeleton[x][y][z] = 2;
		}
		
		if (verbose) BasicInfo.displayMessage("done\n");		
		
       return skeleton;
     }
     
	public final byte[][] medialCurveThinning2D(float[][] funct, int nx, int ny, float startthreshold, float distthreshold) {
		
		// from the boundary, find the middle (approx)
		if (verbose) BasicInfo.displayMessage("fast marching init\n");		
		heap2.reset();
		boolean[][] processed = new boolean[nx][ny];
		boolean[][] medial = new boolean[nx][ny];
		byte[][] skeleton = new byte[nx][ny];
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) if (funct[x][y]>startthreshold) {
        		skeleton[x][y] = 1;
			processed[x][y] = false;
			medial[x][y] = false;
        		
			// search for boundaries
			boolean isboundary=false;
			for (int k = 0; k<4 && !isboundary; k++) {
				if (funct[x+Ngb2.x[k]][y+Ngb2.y[k]]<=startthreshold) {
					isboundary=true;
				}
			}
			// add to the heap: use the distance to boundary as driving force
			if (isboundary) {
				heap2.addValue(funct[x][y], x, y);
				if (debug) BasicInfo.displayMessage(".");
			}
		}
		if (debug) BasicInfo.displayMessage("init\n");		
		
		// thin the object
		while (heap2.isNotEmpty()) {
			// extract point with minimum distance
			float dist = heap2.getFirst();
			int x = heap2.getFirstId1();
			int y = heap2.getFirstId2();
			heap2.removeFirst();

			// if already processed, pass
			if (processed[x][y])  continue;
			
			// check for simple point
			if (isSimplePoint2D(skeleton, x, y)) {
				if (isCurveEndpoint2D(skeleton, x, y) && funct[x][y]>distthreshold) {
					// label as medial surface endpoint
					medial[x][y] = true;
					if (debug) BasicInfo.displayMessage("-");
				} else {
					if (debug) BasicInfo.displayMessage("x");
					skeleton[x][y] = 0;
					
					// find new neghbors: 6C
					for (int k = 0; k<4; k++) {
						if (skeleton[x+Ngb.x[k]][y+Ngb.y[k]]==1 && !processed[x+Ngb.x[k]][y+Ngb.y[k]]) {
							// add to the heap
							heap2.addValue(funct[x+Ngb.x[k]][y+Ngb.y[k]], x+Ngb.x[k], y+Ngb.y[k]);
						}
					}
				}
						
				// discard point
				processed[x][y]=true; 
			} else {
				if (debug) BasicInfo.displayMessage("!");	
			}
		}
		// label the medial endpoints
		for (short x=1;x<nx-1;x++) for (short y=1;y<ny-1;y++) if (medial[x][y]) {
			skeleton[x][y] = 2;
		}
		
		if (verbose) BasicInfo.displayMessage("done\n");		
		
       return skeleton;
     }
     
	private final boolean isMedialEndpoint3D(byte[][][] image, int x, int y, int z) {
		// we assume a binary image with values 0 outside and 1 inside
		for (byte d=0;d<13;d++) {
			float val = 0.0f;
			if (d==0) {		
				val=(image[x][y-1][z] 		+image[x][y+1][z]
					+image[x][y][z-1] 		+image[x][y][z+1]
					+image[x][y-1][z-1] 	+image[x][y-1][z+1]
					+image[x][y+1][z-1]		+image[x][y+1][z+1]);
			} else if (d==1) {
				val=(image[x-1][y][z]		+image[x+1][y][z]
					+image[x][y][z-1]		+image[x][y][z+1]
					+image[x-1][y][z-1]		+image[x-1][y][z+1]
					+image[x+1][y][z-1]		+image[x+1][y][z+1]);
			} else if (d==2) { 			
				val=(image[x-1][y][z]		+image[x+1][y][z]
					 +image[x][y-1][z] 		+image[x][y+1][z]
					 +image[x-1][y-1][z]	+image[x-1][y+1][z]
					 +image[x+1][y-1][z]	+image[x+1][y+1][z]);
			} else if (d==3) { 			
				val=(image[x-1][y+1][z]		+image[x+1][y-1][z]
					 +image[x][y][z-1]		+image[x-1][y+1][z-1]
					 +image[x+1][y-1][z-1] 	+image[x][y][z+1]
					 +image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]);
			} else if (d==4) { 			
				val=(image[x][y+1][z-1]		+image[x][y-1][z+1]
					 +image[x-1][y][z]		+image[x-1][y+1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y][z]	
					 +image[x+1][y+1][z-1]	+image[x+1][y-1][z+1]);
			} else if (d==5) { 			
				val=(image[x+1][y][z-1]		+image[x-1][y][z+1]
					 +image[x][y-1][z]		+image[x+1][y-1][z-1]
					 +image[x-1][y-1][z+1]	+image[x][y+1][z]		
					 +image[x+1][y+1][z-1]	+image[x-1][y+1][z+1]);
			} else if (d==6) { 			
				val=(image[x-1][y-1][z]		+image[x+1][y+1][z]
					 +image[x][y][z-1]		+image[x-1][y-1][z-1]
					 +image[x+1][y+1][z-1]	+image[x][y][z+1]		
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==7) { 			
				val=(image[x][y-1][z-1]		+image[x][y+1][z+1]
					 +image[x-1][y][z]		+image[x-1][y-1][z-1]
					 +image[x-1][y+1][z+1]	+image[x+1][y][z]	
					 +image[x+1][y-1][z-1]	+image[x+1][y+1][z+1]);
			} else if (d==8) { 			
				val=(image[x-1][y][z-1]		+image[x+1][y][z+1]
					 +image[x][y-1][z]		+image[x-1][y-1][z-1]
					 +image[x+1][y-1][z+1]	+image[x][y+1][z]	
					 +image[x-1][y+1][z-1]	+image[x+1][y+1][z+1]);
			} else if (d==9) { 			
				val=(image[x-1][y-1][z+1]	+image[x-1][y+1][z-1]
					 +image[x+1][y-1][z-1]	+image[x+1][y+1][z-1]
					 +image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]);
			} else if (d==10) { 			
				val=(image[x+1][y-1][z+1]	+image[x+1][y+1][z-1]
					 +image[x-1][y-1][z-1]	+image[x-1][y+1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==11) { 			
				val=(image[x-1][y+1][z+1]	+image[x+1][y+1][z-1]
					 +image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]
					 +image[x-1][y-1][z+1]	+image[x+1][y+1][z+1]);
			} else if (d==12) { 			
				val=(image[x-1][y+1][z+1]	+image[x+1][y-1][z+1]
					 +image[x-1][y-1][z-1]	+image[x+1][y-1][z-1]
					 +image[x-1][y+1][z-1]	+image[x+1][y+1][z+1]);
			}
			if (val==1) return true;
		}
		return false;
	}
	
	private final boolean isMedialEndpoint2D(byte[][] image, int x, int y) {
		// we assume a binary image with values 0 outside and 1 inside
		for (byte d=0;d<4;d++) {
			float val = 0.0f;
			if (d==0) {		
				val=(image[x][y-1] 		+image[x][y+1]);
			} else if (d==1) {
				val=(image[x-1][y]		+image[x+1][y]);
			} else if (d==2) { 			
				val=(image[x-1][y+1]		+image[x+1][y-1]);
			} else if (d==6) { 			
				val=(image[x-1][y-1]		+image[x+1][y+1]);
			}
			if (val==1) return true;
		}
		return false;
	}
	
	private final boolean isCurveEndpoint3D(byte[][][] image, int x, int y, int z) {
		// we assume a binary image with values 0 outside and 1 inside
		float val = 0.0f;
		for (int n=0;n<26;n++) {
			byte[] ngb = Ngb.directionNeighbor(n);
			val += image[x+ngb[X]][y+ngb[Y]][z+ngb[Z]];
		}
		if (val==1) return true;
		else return false;
	}
	
	private final boolean isCurveEndpoint2D(byte[][] image, int x, int y) {
		// we assume a binary image with values 0 outside and 1 inside
		float val = 0.0f;
		for (int n=0;n<8;n++) {
			byte[] ngb = Ngb2.directionNeighbor(n);
			val += image[x+ngb[X]][y+ngb[Y]];
		}
		if (val==1) return true;
		else return false;
	}
	
	private final boolean isSimplePoint3D(byte[][][] image, int x, int y, int z) {
    	// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (image[x+i][y+j][z+l]==1) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		return true;
    }

	private final boolean isSimplePoint2D(byte[][] image, int x, int y) {
    	// does it change the topology of the new object ?
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
			if (image[x+i][y+j]==1) {
				obj[1+i][1+j][1+l] = true;
			} else {
				obj[1+i][1+j][1+l] = false;
			}
		}
		obj[1][1][1] = true;
		
		if (!lut.get(lut.keyFromPattern(obj,1,1,1))) return false;
			
		return true;
    }

}
