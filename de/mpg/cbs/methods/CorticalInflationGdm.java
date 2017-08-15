package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm uses the classical GDM framework to evolve a single label
 *	according to internal (curvature) and external (vector field) forces
 *	note: for now, we are assuming isotropic voxels, and boundary conditions are not enforced
 *	(not a problem unless the narrowband reaches a boundary of the image, 
 *	mirroring conditions would solve this problem if needed)
 *
 *  update: the algorithm also performs a tracking of the deformation and its inverse (Jul 2013)
 *
 *	@version    Jul 2010
 *	@author     Pierre-Louis Bazin 
 *	@author     John Bogovic
 * 	@author 	Hanlin Wan
 *		
 *
 */
 
public class CorticalInflationGdm {
	
	// object types

	private	static	final	byte	EMPTY = -1;
	private	static	final	byte	OBJ = 1;
	private	static	final	byte	BG = 0;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    private final static byte T = 3;
    	
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
	private 	float[] 		levelset;  			// level set functions
	private 	byte[] 			segmentation;   	// segmentation
	private 	float[] 		initlevelset;  		// starting levelset
	private		byte[]			labeling;
	private 	float[] 		trglevelset;  		// target levelset
	private		boolean[]		mask;				// masking regions not used in computations
	private		float[][]		mappingFD;			// mapping between original and deformed levelset (forward, direct)
	//private		float[][]		mappingFR;			// mapping between original and deformed levelset (forward, reverse)
	private		float[][]		mappingBD;			// mapping between original and deformed levelset (backward, direct)
	//private		float[][]		mappingBR;			// mapping between original and deformed levelset (backward, reverse)
	private		float[][]		mapped;
	private		boolean[]		mappedMask;
	private		int				nimg;
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		BinaryHeap2D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	
	// parameters
	private	double		smoothweight, balloonweight;
	private	double		stepsize = 0.4;
	private	float		lowlevel = 0.1f;
	private float		landmineDist = 5.0f;
	private	float		narrowBandDist = landmineDist+1.8f;
	
	//private int			maxProject = 1;
	
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
	public CorticalInflationGdm(float[] lvlsetin_,
								int nx_, int ny_, int nz_,
								float rx_, float ry_, float rz_,
								boolean[] mask_,
								float bw_, float sw_,
								String connectivityType_,
								float[][] intens_, boolean[] intensMask_, int nimg_) { 
		
		initlevelset = lvlsetin_;
		
		balloonweight = bw_;
		smoothweight = sw_;
		K0 = 0.0;
		
		mask = mask_;
		
		if (debug) System.out.print("GDM forces: "+bw_+" (balloon), "+sw_+" (smoothing)\n");
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		UNKNOWN = Numerics.max(nx,ny,nz);
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		//maxProject = maxproj_;
		mapped = intens_;
		mappedMask = intensMask_;
		nimg = nimg_;		
		
		// 6,18,26-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1};
		yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, nx, -nx, -nx, nx, -nx, nx, -nx, 0, 0, 0, 0, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, nx*ny, -nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, nx*ny, nx*ny, nx*ny, -nx*ny, -nx*ny, -nx*ny, -nx*ny};
		
		// init all the arrays
		try {
			levelset = new float[nx*ny*nz];
			segmentation = new byte[nx*ny*nz];	
			mappingFD = new float[3][nx*ny*nz];
			//mappingFR = new float[3][nx*ny*nz];
			mappingBD = new float[3][nx*ny*nz];
			//mappingBR = new float[3][nx*ny*nz];
			//mapped = new float[nx*ny*nz];
			trglevelset = new float[nx*ny*nz];
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
		// init decomposition
		fastMarchingInitializationFromInside(false);
				
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		levelset = null;
		segmentation = null;
		trglevelset = null;
		mappingFD = null;
		mappingBD = null;
		mapped = null;
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
	
	public final float[][] getMappedIntensity() { return mapped; }
	
	public final byte[] getSegmentation() { return segmentation; }
    
	public final float[] exportSegmentation() {
		float[] res = new float[nx*ny*nz];
		for (int xyz = 0; xyz<nx*ny*nz; xyz++) res[xyz] = segmentation[xyz];
		return res;
	}

	public final float[][][] exportLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = levelset[xyz];
		}
		return res;
	}

	public final float[][][] exportInitLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = initlevelset[xyz];
		}
		return res;
	}

	public final float[][][] exportTargetLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = trglevelset[xyz];
		}
		return res;
	}

	public final float[][][] generateMappedLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = ImageInterpolation.linearClosestInterpolation(initlevelset, mappingFD[X][xyz], 
																		mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz);
		}
		return res;
	}

	public final float[][][] generateInverseMappedLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = ImageInterpolation.linearClosestInterpolation(levelset, mappingBD[X][xyz], 
																		mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
		}
		return res;
	}

	public final float[][][] generateMappedIntensity(int i) {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z] = mapped[i][xyz];
		}
		return res;
	}
	
	public final byte[][][] generateInitLevelsetLabeling() {
		byte[][][] res = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz]) {
				if (initlevelset[xyz]<0) res[x][y][z] = labeling[xyz];
				else res[x][y][z] = 0;
			}
		}
		return res;
	}

	public final byte[][][] generateMappedLevelsetLabeling() {
		byte[][][] res = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz]) {
				if (ImageInterpolation.linearClosestInterpolation(initlevelset, mappingFD[X][xyz], 
															mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz)<0) 
					res[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(labeling, (byte)0, mappingFD[X][xyz], 
															mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz);
				else res[x][y][z] = 0;
			} 
		}
		return res;
	}

	public final float[][][][] exportMappingFD() {
		float[][][][] res = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z][X] = mappingFD[X][xyz];
			res[x][y][z][Y] = mappingFD[Y][xyz];
			res[x][y][z][Z] = mappingFD[Z][xyz];
		}
		return res;
	}

	public final float[][][][] exportMappingBD() {
		float[][][][] res = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z][X] = mappingBD[X][xyz];
			res[x][y][z][Y] = mappingBD[Y][xyz];
			res[x][y][z][Z] = mappingBD[Z][xyz];
		}
		return res;
	}

	public final float[][][][] generateDifferenceMapping() {
		float[][][][] res = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z][X] = mappingFD[X][xyz]-x;
			res[x][y][z][Y] = mappingFD[Y][xyz]-y;
			res[x][y][z][Z] = mappingFD[Z][xyz]-z;
		}
		return res;
	}

	public final void cleanupForwardMapping() {
		float[][] res = new float[3][nx*ny*nz];
		BitSet undef = new BitSet(nx*ny*nz);
		
		// 1. transform the original levelset with the final transform
		float[] lvl = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lvl[xyz] = ImageInterpolation.linearClosestInterpolation(initlevelset, mappingFD[X][xyz], 
																		mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz);
		}
		// 2. find regions of strong differences : only different signs??
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			//if (Numerics.abs(lvl[xyz]-levelset[xyz]) > 1.0f) undef.set(xyz);
			if (lvl[xyz]*levelset[xyz]<0) undef.set(xyz);
		}
		// 3. adjust the mapping there
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (!undef.get(xyz)) {
				res[X][xyz] = mappingFD[X][xyz];
				res[Y][xyz] = mappingFD[Y][xyz];
				res[Z][xyz] = mappingFD[Z][xyz];
			} else {
				res[X][xyz] = 0.0f;
				res[Y][xyz] = 0.0f;
				res[Z][xyz] = 0.0f;
				float nb=0.0f;
				for (int k=0;k<26;k++) if (!undef.get(xyz+xoff[k]+yoff[k]+zoff[k])) {
					// weight with distance between values?
					//float w = 1.0f-Numerics.abs(lvl[xyz+xoff[k]+yoff[k]+zoff[k]]-levelset[xyz+xoff[k]+yoff[k]+zoff[k]]);
					// weight with distance to boundary?
					float w = 1.0f/(1.0f + Numerics.abs(levelset[xyz+xoff[k]+yoff[k]+zoff[k]]));
					res[X][xyz] += w*mappingFD[X][xyz+xoff[k]+yoff[k]+zoff[k]];
					res[Y][xyz] += w*mappingFD[Y][xyz+xoff[k]+yoff[k]+zoff[k]];
					res[Z][xyz] += w*mappingFD[Z][xyz+xoff[k]+yoff[k]+zoff[k]];
					nb+= w;
				}
				if (nb==0) {
					// unchanged
					res[X][xyz] = mappingFD[X][xyz];
					res[Y][xyz] = mappingFD[Y][xyz];
					res[Z][xyz] = mappingFD[Z][xyz];
				} else {
					res[X][xyz] /= nb;
					res[Y][xyz] /= nb;
					res[Z][xyz] /= nb;
				}
				// reproject again?
				float err = levelset[xyz] - ImageInterpolation.linearInterpolation(initlevelset, initlevelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
				float ratio = 1.0f;
				float sgn = levelset[xyz]*ImageInterpolation.linearInterpolation(initlevelset, initlevelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
				//for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
				for (int e=0;e<50 && sgn<0;e++) {
					float err0 = err;
					// get levelset gradient at init location
					float dX = ImageInterpolation.linearInterpolationXderivative(initlevelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float dY = ImageInterpolation.linearInterpolationYderivative(initlevelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float dZ = ImageInterpolation.linearInterpolationZderivative(initlevelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
					
					res[X][xyz] += ratio*err/norm2*dX;
					res[Y][xyz] += ratio*err/norm2*dY;
					res[Z][xyz] += ratio*err/norm2*dZ;
					
					sgn = levelset[xyz]*ImageInterpolation.linearInterpolation(initlevelset, initlevelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					if (sgn<0) {
						err = levelset[xyz] - ImageInterpolation.linearInterpolation(initlevelset, initlevelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
						if (Numerics.abs(err)>Numerics.abs(err0)) {
							res[X][xyz] -= ratio*err/norm2*dX;
							res[Y][xyz] -= ratio*err/norm2*dY;
							res[Z][xyz] -= ratio*err/norm2*dZ;
							ratio *= 0.67f;
							err = err0;
						} else {
							ratio = 1.0f;
						}
					}
				}
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			mappingFD[X][xyz] = res[X][xyz];
			mappingFD[Y][xyz] = res[Y][xyz];
			mappingFD[Z][xyz] = res[Z][xyz];
		}
		return;
	}

	public final void cleanupBackwardMapping() {
		float[][] res = new float[3][nx*ny*nz];
		BitSet undef = new BitSet(nx*ny*nz);
		
		// 1. transform the original levelset with the final transform
		float[] lvl = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lvl[xyz] = ImageInterpolation.linearClosestInterpolation(levelset, mappingBD[X][xyz], 
																		mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
		}
		// 2. find regions of strong differences : only different signs??
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			//if (Numerics.abs(lvl[xyz]-levelset[xyz]) > 1.0f) undef.set(xyz);
			if (lvl[xyz]*initlevelset[xyz]<0) undef.set(xyz);
		}
		// 3. adjust the mapping there
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (!undef.get(xyz)) {
				res[X][xyz] = mappingBD[X][xyz];
				res[Y][xyz] = mappingBD[Y][xyz];
				res[Z][xyz] = mappingBD[Z][xyz];
			} else {
				res[X][xyz] = 0.0f;
				res[Y][xyz] = 0.0f;
				res[Z][xyz] = 0.0f;
				float nb=0.0f;
				for (int k=0;k<26;k++) if (!undef.get(xyz+xoff[k]+yoff[k]+zoff[k])) {
					// weight with distance between values?
					//float w = 1.0f-Numerics.abs(lvl[xyz+xoff[k]+yoff[k]+zoff[k]]-levelset[xyz+xoff[k]+yoff[k]+zoff[k]]);
					// weight with distance to boundary?
					float w = 1.0f/(1.0f + Numerics.abs(initlevelset[xyz+xoff[k]+yoff[k]+zoff[k]]));
					res[X][xyz] += w*mappingBD[X][xyz+xoff[k]+yoff[k]+zoff[k]];
					res[Y][xyz] += w*mappingBD[Y][xyz+xoff[k]+yoff[k]+zoff[k]];
					res[Z][xyz] += w*mappingBD[Z][xyz+xoff[k]+yoff[k]+zoff[k]];
					nb+= w;
				}
				if (nb==0) {
					// unchanged
					res[X][xyz] = mappingBD[X][xyz];
					res[Y][xyz] = mappingBD[Y][xyz];
					res[Z][xyz] = mappingBD[Z][xyz];
				} else {
					res[X][xyz] /= nb;
					res[Y][xyz] /= nb;
					res[Z][xyz] /= nb;
				}
				// reproject again?
				float err = initlevelset[xyz] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
				float ratio = 1.0f;
				float sgn = initlevelset[xyz]*ImageInterpolation.linearInterpolation(levelset, levelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
				//for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
				for (int e=0;e<50 && sgn<0;e++) {
					float err0 = err;
					// get levelset gradient at init location
					float dX = ImageInterpolation.linearInterpolationXderivative(levelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float dY = ImageInterpolation.linearInterpolationYderivative(levelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float dZ = ImageInterpolation.linearInterpolationZderivative(levelset, res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
					
					res[X][xyz] += ratio*err/norm2*dX;
					res[Y][xyz] += ratio*err/norm2*dY;
					res[Z][xyz] += ratio*err/norm2*dZ;
					
					sgn = initlevelset[xyz]*ImageInterpolation.linearInterpolation(levelset, levelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
					if (sgn<0) {
						err = initlevelset[xyz] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], res[X][xyz], res[Y][xyz], res[Z][xyz], nx, ny, nz);
						if (Numerics.abs(err)>Numerics.abs(err0)) {
							res[X][xyz] -= ratio*err/norm2*dX;
							res[Y][xyz] -= ratio*err/norm2*dY;
							res[Z][xyz] -= ratio*err/norm2*dZ;
							ratio *= 0.67f;
							err = err0;
						} else {
							ratio = 1.0f;
						}
					}
				}
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			mappingBD[X][xyz] = res[X][xyz];
			mappingBD[Y][xyz] = res[Y][xyz];
			mappingBD[Z][xyz] = res[Z][xyz];
		}
		return;
	}

	public final float[][][][] generateInverseDifferenceMapping() {
		float[][][][] res = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z][X] = mappingBD[X][xyz]-x;
			res[x][y][z][Y] = mappingBD[Y][xyz]-y;
			res[x][y][z][Z] = mappingBD[Z][xyz]-z;
		}
		return res;
	}

	public final float[][][][] generateMappingDistances() {
		float[][][][] res = new float[nx][ny][nz][4];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[x][y][z][X] = Numerics.square(mappingFD[X][xyz]-x) + Numerics.square(mappingFD[Y][xyz]-y) + Numerics.square(mappingFD[Z][xyz]-z);
			res[x][y][z][Y] = Numerics.square(mappingBD[X][xyz]-x) + Numerics.square(mappingBD[Y][xyz]-y) + Numerics.square(mappingBD[Z][xyz]-z);
			res[x][y][z][Z] = Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingBD[X], mappingFD[X][xyz], mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz)-x)
							  +Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingBD[Y], mappingFD[X][xyz], mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz)-y)
							  +Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingBD[Z], mappingFD[X][xyz], mappingFD[Y][xyz], mappingFD[Z][xyz], nx, ny, nz)-z);
			res[x][y][z][T] = Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingFD[X], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz)-x)
							  +Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingFD[Y], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz)-y)
							  +Numerics.square(ImageInterpolation.linearClosestInterpolation(mappingFD[Z], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz)-z);
			
			res[x][y][z][X] = (float)FastMath.sqrt(res[x][y][z][X]);
			res[x][y][z][Y] = (float)FastMath.sqrt(res[x][y][z][Y]);
			res[x][y][z][Z] = (float)FastMath.sqrt(res[x][y][z][Z]);
			res[x][y][z][T] = (float)FastMath.sqrt(res[x][y][z][T]);
		}
		return res;
	}

	public final float[][][][] generateGridMapping() {
		float[][][][] res = new float[nx][ny][nz][3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (Numerics.round(mappingFD[X][xyz])%10==0) res[x][y][z][X] = 1.0f;
			else res[x][y][z][X] = 0.0f;
			if (Numerics.round(mappingFD[Y][xyz])%10==0) res[x][y][z][Y] = 1.0f;
			else res[x][y][z][Y] = 0.0f;
			if (Numerics.round(mappingFD[Z][xyz])%10==0) res[x][y][z][Z] = 1.0f;
			else res[x][y][z][Z] = 0.0f;
		}
		return res;
	}

	public final void fastMarchingInitializationFromInside(boolean narrowBandOnly) {
		
        // initialize the quantities
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // gdm function
			levelset[xyz] = initlevelset[xyz];                            
           
			// segmentation
            if (initlevelset[xyz]<0) {
				segmentation[xyz] = OBJ;
			} else {
				segmentation[xyz] = BG;
			}
        }
        
        // initial mapping: identity
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	mappingFD[X][xyz] = x;
        	mappingFD[Y][xyz] = y;
        	mappingFD[Z][xyz] = z;
        	
        	mappingBD[X][xyz] = x;
        	mappingBD[Y][xyz] = y;
        	mappingBD[Z][xyz] = z;
		}
		
        // re-compute the level set (everywhere that's not masked)
        fastMarchingReinitialization(false);
        
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) {
         	 // copy the gdm function
			initlevelset[xyz] = levelset[xyz];
			trglevelset[xyz] = levelset[xyz];        	
        	//mapped[xyz] = levelset[xyz];
        }        
        
        return;
     }
     
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingReinitialization(boolean narrowBandOnly) {
        // computation variables
        boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
		// initialize the heap from boundaries
        for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
        	processed[xyz] = false;
        	// search for boundaries
        	for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
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
			for (int k = 0; k<6; k++) {
				int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
				
				if (mask[xyzn]) {
					// must be in outside the object or its processed neighborhood
					if (segmentation[xyzn]==lb && !processed[xyzn]) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xyznb = xyzn + xoff[l] + yoff[l] + zoff[l];
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
    public final float[] generateForces() {
    	float[] forces = new float[nx*ny*nz];
    	for (int xyz=0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		forces[xyz] = (float)levelsetForces(xyz);
    	}
    	return forces;
    }
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    
    public final void evolveNarrowBand(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band\n");

    	double[] vect = new double[4];
    	
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
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				prev = narrowband.levels[n];
				//narrowband.levels[n] -= Numerics.bounded(levelsetForces(xyz),-0.9,0.9);
				levelsetForceVector(xyz, vect);
				narrowband.levels[n] -= vect[3];
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						// update the segmentation
						segmentation[xyz] = OBJ;
					} else {
						narrowband.levels[n] = lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						// do nothing
						segmentation[xyz] = BG;
					} else {
						narrowband.levels[n] = -lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
			}
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
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
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBandMappingNearest(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band + mapping\n");

    	double[] vect = new double[4];
    	double[] pt = new double[3];
    	double delta;
    	//double phi, dx, dy, dz, grad;

		//float[][] newmapping = new float[3][nx*ny*nz];
		float[][] tmpmap = new float[nimg][nx*ny*nz];
		boolean[] newMappedMask = new boolean[nx*ny*nz];
		
		float[][] direct = new float[3][nx*ny*nz];
		float[][] inverse = new float[4][nx*ny*nz];
		float[][] tmpmapping = new float[3][nx*ny*nz];
		
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
			/*
			newmapping[X][xyz] = mappingFD[X][xyz];
			newmapping[Y][xyz] = mappingFD[Y][xyz];
			newmapping[Z][xyz] = mappingFD[Z][xyz];
			*/
			
			// direct mapping
			// coordinate mapping
			pt[Z] = Numerics.floor(xyz/(nx*ny));
			pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
			pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
			
			direct[X][xyz] = (float)pt[X];
			direct[Y][xyz] = (float)pt[Y];
			direct[Z][xyz] = (float)pt[Z];
			
			inverse[X][xyz] = 0.0f;
			inverse[Y][xyz] = 0.0f;
			inverse[Z][xyz] = 0.0f;
			inverse[T][xyz] = 0.0f;
		}
			
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && diff>mindiff;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinit = false;
			int nswap=0;
			double prev;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				prev = narrowband.levels[n];
				levelsetForceVector(xyz, vect);
				narrowband.levels[n] -= vect[3];
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						// update the segmentation
						segmentation[xyz] = OBJ;
					} else {
						delta = narrowband.levels[n]-lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						// do nothing
						segmentation[xyz] = BG;
					} else {
						delta = narrowband.levels[n]+lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = -lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
				// coordinate mapping
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				// move to the new location
				pt[X] += -vect[X];
				pt[Y] += -vect[Y];
				pt[Z] += -vect[Z];
				
				/*
				// project onto the corresponding levelset (i.e. find closest levelset with new phi value and project)
				// this should be a minor adjustment
				double np = 1.0;
				phi = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				while (Numerics.abs(phi - narrowband.levels[n])>0.001f && np<=maxProject) {
					dx = (ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nx, ny, nz)
						 -ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nx, ny, nz) )/2.0;
					dy = (ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nx, ny, nz)
						 -ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nx, ny, nz) )/2.0;
					dz = (ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nx, ny, nz)
						 -ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nx, ny, nz) )/2.0;
					
					grad = Math.sqrt(dx*dx + dy*dy + dz*dz);
					
					pt[X] -= (phi - narrowband.levels[n])*dx/grad;
					pt[Y] -= (phi - narrowband.levels[n])*dy/grad;
					pt[Z] -= (phi - narrowband.levels[n])*dz/grad;
				
					np++;
					phi = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				}
				*/
				
				direct[X][xyz] = (float)pt[X];
				direct[Y][xyz] = (float)pt[Y];
				direct[Z][xyz] = (float)pt[Z];
				
				/*
				// build inverse mapping? linear interpolation?
				int xyz0 = Numerics.floor(pt[X])+Numerics.floor(pt[Y])*nx+Numerics.floor(pt[Z])*nx*ny;
				double wx = pt[X]-Numerics.floor(pt[X]);
				double wy = pt[Y]-Numerics.floor(pt[Y]);
				double wz = pt[Z]-Numerics.floor(pt[Z]);
				inverse[X][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[X]);
				inverse[Y][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[Y]);
				inverse[Z][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[Z]);
				inverse[T][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz));
				inverse[X][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz)*vect[X]);
				inverse[Y][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz)*vect[Y]);
				inverse[Z][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz)*vect[Z]);
				inverse[T][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz));
				inverse[X][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz)*vect[X]);
				inverse[Y][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz)*vect[Y]);
				inverse[Z][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz)*vect[Z]);
				inverse[T][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz));
				inverse[X][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz*vect[X]);
				inverse[Y][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz*vect[Y]);
				inverse[Z][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz*vect[Z]);
				inverse[T][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz);
				inverse[X][xyz0+1+nx] += (float)( wx*wy*(1.0-wz)*vect[X]);
				inverse[Y][xyz0+1+nx] += (float)( wx*wy*(1.0-wz)*vect[Y]);
				inverse[Z][xyz0+1+nx] += (float)( wx*wy*(1.0-wz)*vect[Z]);
				inverse[T][xyz0+1+nx] += (float)( wx*wy*(1.0-wz));
				inverse[X][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz*vect[X]);
				inverse[Y][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz*vect[Y]);
				inverse[Z][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz*vect[Z]);
				inverse[T][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz);
				inverse[X][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz*vect[X]);
				inverse[Y][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz*vect[Y]);
				inverse[Z][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz*vect[Z]);
				inverse[T][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz);
				inverse[X][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz*vect[X]);
				inverse[Y][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz*vect[Y]);
				inverse[Z][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz*vect[Z]);
				inverse[T][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz);
				*/
				// build inverse mapping? maximum neighbor
				int xyz0 = Numerics.floor(pt[X])+Numerics.floor(pt[Y])*nx+Numerics.floor(pt[Z])*nx*ny;
				double wx = pt[X]-Numerics.floor(pt[X]);
				double wy = pt[Y]-Numerics.floor(pt[Y]);
				double wz = pt[Z]-Numerics.floor(pt[Z]);
				if ( (1.0-wx)*(1.0-wy)*(1.0-wz) > inverse[T][xyz0]) {
					inverse[X][xyz0] = (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[X]);
					inverse[Y][xyz0] = (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[Y]);
					inverse[Z][xyz0] = (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*vect[Z]);
					inverse[T][xyz0] = (float)( (1.0-wx)*(1.0-wy)*(1.0-wz));
				}
				if ( wx*(1.0-wy)*(1.0-wz) > inverse[T][xyz0+1]) {
					inverse[X][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz)*vect[X]);
					inverse[Y][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz)*vect[Y]);
					inverse[Z][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz)*vect[Z]);
					inverse[T][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz));
				}
				if ( (1.0-wx)*wy*(1.0-wz) > inverse[T][xyz0+nx]) {
					inverse[X][xyz0+nx] = (float)( (1.0-wx)*wy*(1.0-wz)*vect[X]);
					inverse[Y][xyz0+nx] = (float)( (1.0-wx)*wy*(1.0-wz)*vect[Y]);
					inverse[Z][xyz0+nx] = (float)( (1.0-wx)*wy*(1.0-wz)*vect[Z]);
					inverse[T][xyz0+nx] = (float)( (1.0-wx)*wy*(1.0-wz));
				}
				if ( (1.0-wx)*(1.0-wy)*wz > inverse[T][xyz0+nx*ny]) {
					inverse[X][xyz0+nx*ny] = (float)( (1.0-wx)*(1.0-wy)*wz*vect[X]);
					inverse[Y][xyz0+nx*ny] = (float)( (1.0-wx)*(1.0-wy)*wz*vect[Y]);
					inverse[Z][xyz0+nx*ny] = (float)( (1.0-wx)*(1.0-wy)*wz*vect[Z]);
					inverse[T][xyz0+nx*ny] = (float)( (1.0-wx)*(1.0-wy)*wz);
				}
				if ( wx*wy*(1.0-wz) > inverse[T][xyz0+1+nx]) {
					inverse[X][xyz0+1+nx] = (float)( wx*wy*(1.0-wz)*vect[X]);
					inverse[Y][xyz0+1+nx] = (float)( wx*wy*(1.0-wz)*vect[Y]);
					inverse[Z][xyz0+1+nx] = (float)( wx*wy*(1.0-wz)*vect[Z]);
					inverse[T][xyz0+1+nx] = (float)( wx*wy*(1.0-wz));
				}
				if ( (1.0-wx)*wy*wz > inverse[T][xyz0+nx+nx*ny]) {
					inverse[X][xyz0+nx+nx*ny] =  (float)( (1.0-wx)*wy*wz*vect[X]);
					inverse[Y][xyz0+nx+nx*ny] =  (float)( (1.0-wx)*wy*wz*vect[Y]);
					inverse[Z][xyz0+nx+nx*ny] =  (float)( (1.0-wx)*wy*wz*vect[Z]);
					inverse[T][xyz0+nx+nx*ny] =  (float)( (1.0-wx)*wy*wz);
				}
				if ( wx*(1.0-wy)*wz > inverse[T][xyz0+nx*ny+1]) {
					inverse[X][xyz0+nx*ny+1] = (float)( wx*(1.0-wy)*wz*vect[X]);
					inverse[Y][xyz0+nx*ny+1] = (float)( wx*(1.0-wy)*wz*vect[Y]);
					inverse[Z][xyz0+nx*ny+1] = (float)( wx*(1.0-wy)*wz*vect[Z]);
					inverse[T][xyz0+nx*ny+1] = (float)( wx*(1.0-wy)*wz);
				}
				if ( wx*wy*wz > inverse[T][xyz0+1+nx+nx*ny]) {
					inverse[X][xyz0+1+nx+nx*ny] = (float)( wx*wy*wz*vect[X]);
					inverse[Y][xyz0+1+nx+nx*ny] = (float)( wx*wy*wz*vect[Y]);
					inverse[Z][xyz0+1+nx+nx*ny] = (float)( wx*wy*wz*vect[Z]);
					inverse[T][xyz0+1+nx+nx*ny] = (float)( wx*wy*wz);
				}
				// build inverse mapping? most likely neighbor??
				
				/*
				// compose with current mapping
				newmapping[X][xyz] = ImageInterpolation.linearInterpolation(mappingFD[X], 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				newmapping[Y][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Y], 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				newmapping[Z][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Z], 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				*/
				
				// compose intensity
				if (mapped != null) {
					for (int i=0;i<nimg;i++) {
						tmpmap[i][xyz] = ImageInterpolation.linearInterpolation(mapped[i], mappedMask, 0.0f, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
						newMappedMask[xyz] = (tmpmap[i][xyz] != 0.0f); 
					}
				}
			}
			
			mappedMask = Arrays.copyOf(newMappedMask, nx*ny*nz);
			
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
			}
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				if (inverse[T][xyz]>0.001) {
					inverse[X][xyz] /= inverse[T][xyz];
					inverse[Y][xyz] /= inverse[T][xyz];
					inverse[Z][xyz] /= inverse[T][xyz];
				} else {
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
				}					
			}
			
			// forward mapping
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				/*
				pt[Z] = Numerics.floor(xyz/(nx*ny)) + inverse[Z][xyz];
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx) + inverse[Y][xyz];
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx) + inverse[X][xyz];
				
				tmpmapping[X][xyz] = ImageInterpolation.linearInterpolation(mappingFD[X], (float)pt[X], pt[X], pt[Y], pt[Z], nx, ny, nz);
				tmpmapping[Y][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Y], (float)pt[Y], pt[X], pt[Y], pt[Z], nx, ny, nz);
				tmpmapping[Z][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Z], (float)pt[Z], pt[X], pt[Y], pt[Z], nx, ny, nz);
				*/
				tmpmapping[X][xyz] = ImageInterpolation.linearInterpolation(mappingFD[X], direct[X][xyz], direct[X][xyz], direct[Y][xyz], direct[Z][xyz], nx, ny, nz);
				tmpmapping[Y][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Y], direct[Y][xyz], direct[X][xyz], direct[Y][xyz], direct[Z][xyz], nx, ny, nz);
				tmpmapping[Z][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Z], direct[Z][xyz], direct[X][xyz], direct[Y][xyz], direct[Z][xyz], nx, ny, nz);
				
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingFD[X][xyz] = tmpmapping[X][xyz];
				mappingFD[Y][xyz] = tmpmapping[Y][xyz];
				mappingFD[Z][xyz] = tmpmapping[Z][xyz];
			}
			// backward mapping
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				/*
				tmpmapping[X][xyz] = ImageInterpolation.linearInterpolation(direct[X], mappingBD[X][xyz], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Y][xyz] = ImageInterpolation.linearInterpolation(direct[Y], mappingBD[Y][xyz], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Z][xyz] = ImageInterpolation.linearInterpolation(direct[Z], mappingBD[Z][xyz], mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				*/
				tmpmapping[X][xyz] = mappingBD[X][xyz] + ImageInterpolation.linearInterpolation(inverse[X], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Y][xyz] = mappingBD[Y][xyz] + ImageInterpolation.linearInterpolation(inverse[Y], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Z][xyz] = mappingBD[Z][xyz] + ImageInterpolation.linearInterpolation(inverse[Z], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingBD[X][xyz] = tmpmapping[X][xyz];
				mappingBD[Y][xyz] = tmpmapping[Y][xyz];
				mappingBD[Z][xyz] = tmpmapping[Z][xyz];
			}
			// reinit the inverse values
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				inverse[X][xyz] = 0.0f;
				inverse[Y][xyz] = 0.0f;
				inverse[Z][xyz] = 0.0f;
				inverse[T][xyz] = 0.0f;				
			}
			/*
			// mapping update
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingFD[X][xyz] = newmapping[X][xyz];
				mappingFD[Y][xyz] = newmapping[Y][xyz];
				mappingFD[Z][xyz] = newmapping[Z][xyz];
			}
			*/
			// intensity mapping
			if (mapped != null) {
				for (int i=0;i<nimg;i++) {
					for (int n=0; n<narrowband.currentsize;n++) {
						xyz = narrowband.id[n];
					
						mapped[i][xyz] = tmpmap[i][xyz];
					}
				}
			}
			
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// clean up the transform in the old narrow band
				for (int n=0; n<narrowband.currentsize;n++) {
					xyz = narrowband.id[n];
					
					pt[Z] = Numerics.floor(xyz/(nx*ny));
					pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
					pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
					
					direct[X][xyz] = (float)pt[X];
					direct[Y][xyz] = (float)pt[Y];
					direct[Z][xyz] = (float)pt[Z];
					
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
					inverse[T][xyz] = 0.0f;
				}
				
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
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBandMappingClosestMean(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band + mapping\n");

    	double[] vect = new double[4];
    	double[] pt = new double[3];
    	double delta;
    	//double phi, dx, dy, dz, grad;

		//float[][] newmapping = new float[3][nx*ny*nz];
		float[][] tmpmap = new float[nimg][nx*ny*nz];
		boolean[] newMappedMask = new boolean[nx*ny*nz];
		
		float[][] direct = new float[3][nx*ny*nz];
		float[][] inverse = new float[4][nx*ny*nz];
		float[][] tmpmapping = new float[3][nx*ny*nz];
		
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
			direct[X][xyz] = 0.0f;
			direct[Y][xyz] = 0.0f;
			direct[Z][xyz] = 0.0f;
		}
			
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && diff>mindiff;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinit = false;
			int nswap=0;
			double prev;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				prev = narrowband.levels[n];
				levelsetForceVector(xyz, vect);
				narrowband.levels[n] -= vect[T];
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						// update the segmentation
						segmentation[xyz] = OBJ;
					} else {
						delta = narrowband.levels[n]-lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						// do nothing
						segmentation[xyz] = BG;
					} else {
						delta = narrowband.levels[n]+lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = -lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
				direct[X][xyz] = -(float)vect[X];
				direct[Y][xyz] = -(float)vect[Y];
				direct[Z][xyz] = -(float)vect[Z];
			}
			
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
			}
			
			// forward mapping
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				pt[X] += direct[X][xyz];
				pt[Y] += direct[Y][xyz];
				pt[Z] += direct[Z][xyz];
				
				tmpmapping[X][xyz] = ImageInterpolation.linearInterpolation(mappingFD[X], (float)pt[X], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				tmpmapping[Y][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Y], (float)pt[Y], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				tmpmapping[Z][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Z], (float)pt[Z], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				
				// compose intensity
				if (mapped != null) {
					for (int i=0;i<nimg;i++) {
						tmpmap[i][xyz] = ImageInterpolation.linearInterpolation(mapped[i], mappedMask, 0.0f, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
						newMappedMask[xyz] = (tmpmap[i][xyz] != 0.0f); 
					}
				}
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingFD[X][xyz] = tmpmapping[X][xyz];
				mappingFD[Y][xyz] = tmpmapping[Y][xyz];
				mappingFD[Z][xyz] = tmpmapping[Z][xyz];
			}
			mappedMask = Arrays.copyOf(newMappedMask, nx*ny*nz);
			
			// build the inverse transform from the direct one: get the "most representative" vector?
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				pt[X] += direct[X][xyz];
				pt[Y] += direct[Y][xyz];
				pt[Z] += direct[Z][xyz];
				
				int xyz0 = Numerics.floor(pt[X])+Numerics.floor(pt[Y])*nx+Numerics.floor(pt[Z])*nx*ny;
				double wx = pt[X]-Numerics.floor(pt[X]);
				double wy = pt[Y]-Numerics.floor(pt[Y]);
				double wz = pt[Z]-Numerics.floor(pt[Z]);
				inverse[X][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[X][xyz]);
				inverse[Y][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[Y][xyz]);
				inverse[Z][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[Z][xyz]);
				inverse[T][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz));
				inverse[X][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[X][xyz]);
				inverse[Y][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[Y][xyz]);
				inverse[Z][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[Z][xyz]);
				inverse[T][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz));
				inverse[X][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[X][xyz]);
				inverse[Y][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[Y][xyz]);
				inverse[Z][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[Z][xyz]);
				inverse[T][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz));
				inverse[X][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[X][xyz]);
				inverse[Y][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[Y][xyz]);
				inverse[Z][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[Z][xyz]);
				inverse[T][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz);
				inverse[X][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[X][xyz]);
				inverse[Y][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[Y][xyz]);
				inverse[Z][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[Z][xyz]);
				inverse[T][xyz0+1+nx] += (float)( wx*wy*(1.0-wz));
				inverse[X][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[X][xyz]);
				inverse[Y][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[Y][xyz]);
				inverse[Z][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[Z][xyz]);
				inverse[T][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz);
				inverse[X][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[X][xyz]);
				inverse[Y][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[Y][xyz]);
				inverse[Z][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[Z][xyz]);
				inverse[T][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz);
				inverse[X][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[X][xyz]);
				inverse[Y][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[Y][xyz]);
				inverse[Z][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[Z][xyz]);
				inverse[T][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz);
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				if (inverse[T][xyz]>0.001) {
					inverse[X][xyz] = inverse[X][xyz]/inverse[T][xyz];
					inverse[Y][xyz] = inverse[Y][xyz]/inverse[T][xyz];
					inverse[Z][xyz] = inverse[Z][xyz]/inverse[T][xyz];
				} else {
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
				}
				inverse[T][xyz] = 0.0f;
			}
			// second pass: finds the vector with closest value to the average
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				double corr = -(direct[X][xyz]*inverse[X][xyz]+direct[Y][xyz]*inverse[Y][xyz]+direct[Z][xyz]*inverse[Z][xyz])
							 /FastMath.sqrt(direct[X][xyz]*direct[X][xyz]+direct[Y][xyz]*direct[Y][xyz]+direct[Z][xyz]*direct[Z][xyz])
							 /FastMath.sqrt(inverse[X][xyz]*inverse[X][xyz]+inverse[Y][xyz]*inverse[Y][xyz]+inverse[Z][xyz]*inverse[Z][xyz]);
							 	
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				pt[X] += direct[X][xyz];
				pt[Y] += direct[Y][xyz];
				pt[Z] += direct[Z][xyz];
				
				int xyz0 = Numerics.floor(pt[X])+Numerics.floor(pt[Y])*nx+Numerics.floor(pt[Z])*nx*ny;
				double wx = pt[X]-Numerics.floor(pt[X]);
				double wy = pt[Y]-Numerics.floor(pt[Y]);
				double wz = pt[Z]-Numerics.floor(pt[Z]);
				if ( (1.0-wx)*(1.0-wy)*(1.0-wz)*corr > inverse[T][xyz0]) {
					inverse[X][xyz0] = -direct[X][xyz];
					inverse[Y][xyz0] = -direct[Y][xyz];
					inverse[Z][xyz0] = -direct[Z][xyz];
					inverse[T][xyz0] = (float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*corr);
				}
				if ( wx*(1.0-wy)*(1.0-wz)*corr > inverse[T][xyz0+1]) {
					inverse[X][xyz0+1] = -direct[X][xyz];
					inverse[Y][xyz0+1] = -direct[Y][xyz];
					inverse[Z][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz)*vect[Z]);
					inverse[T][xyz0+1] = (float)( wx*(1.0-wy)*(1.0-wz)*corr);
				}
				if ( (1.0-wx)*wy*(1.0-wz)*corr > inverse[T][xyz0+nx]) {
					inverse[X][xyz0+nx] = -direct[X][xyz];
					inverse[Y][xyz0+nx] = -direct[Y][xyz];
					inverse[Z][xyz0+nx] = -direct[Z][xyz];
					inverse[T][xyz0+nx] = (float)( (1.0-wx)*wy*(1.0-wz)*corr);
				}
				if ( (1.0-wx)*(1.0-wy)*wz*corr > inverse[T][xyz0+nx*ny]) {
					inverse[X][xyz0+nx*ny] = -direct[X][xyz];
					inverse[Y][xyz0+nx*ny] = -direct[Y][xyz];
					inverse[Z][xyz0+nx*ny] = -direct[Z][xyz];
					inverse[T][xyz0+nx*ny] = (float)( (1.0-wx)*(1.0-wy)*wz*corr);
				}
				if ( wx*wy*(1.0-wz)*corr > inverse[T][xyz0+1+nx]) {
					inverse[X][xyz0+1+nx] = -direct[X][xyz];
					inverse[Y][xyz0+1+nx] = -direct[Y][xyz];
					inverse[Z][xyz0+1+nx] = -direct[Z][xyz];
					inverse[T][xyz0+1+nx] = (float)( wx*wy*(1.0-wz)*corr);
				}
				if ( (1.0-wx)*wy*wz*corr > inverse[T][xyz0+nx+nx*ny]) {
					inverse[X][xyz0+nx+nx*ny] =  -direct[X][xyz];
					inverse[Y][xyz0+nx+nx*ny] =  -direct[Y][xyz];
					inverse[Z][xyz0+nx+nx*ny] =  -direct[Z][xyz];
					inverse[T][xyz0+nx+nx*ny] =  (float)( (1.0-wx)*wy*wz*corr);
				}
				if ( wx*(1.0-wy)*wz*corr > inverse[T][xyz0+nx*ny+1]) {
					inverse[X][xyz0+nx*ny+1] = -direct[X][xyz];
					inverse[Y][xyz0+nx*ny+1] = -direct[Y][xyz];
					inverse[Z][xyz0+nx*ny+1] = -direct[Z][xyz];
					inverse[T][xyz0+nx*ny+1] = (float)( wx*(1.0-wy)*wz*corr);
				}
				if ( wx*wy*wz*corr > inverse[T][xyz0+1+nx+nx*ny]) {
					inverse[X][xyz0+1+nx+nx*ny] = -direct[X][xyz];
					inverse[Y][xyz0+1+nx+nx*ny] = -direct[Y][xyz];
					inverse[Z][xyz0+1+nx+nx*ny] = -direct[Z][xyz];
					inverse[T][xyz0+1+nx+nx*ny] = (float)( wx*wy*wz*corr);
				}
			}
			// backward mapping
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				tmpmapping[X][xyz] = mappingBD[X][xyz] + ImageInterpolation.linearInterpolation(inverse[X], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Y][xyz] = mappingBD[Y][xyz] + ImageInterpolation.linearInterpolation(inverse[Y], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Z][xyz] = mappingBD[Z][xyz] + ImageInterpolation.linearInterpolation(inverse[Z], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingBD[X][xyz] = tmpmapping[X][xyz];
				mappingBD[Y][xyz] = tmpmapping[Y][xyz];
				mappingBD[Z][xyz] = tmpmapping[Z][xyz];
			}
			// reinit the inverse values
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				inverse[X][xyz] = 0.0f;
				inverse[Y][xyz] = 0.0f;
				inverse[Z][xyz] = 0.0f;
				inverse[T][xyz] = 0.0f;				
			}
			// intensity mapping
			if (mapped != null) {
				for (int i=0;i<nimg;i++) {
					for (int n=0; n<narrowband.currentsize;n++) {
						xyz = narrowband.id[n];
					
						mapped[i][xyz] = tmpmap[i][xyz];
					}
				}
			}
			
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// clean up the transform in the old narrow band
				for (int n=0; n<narrowband.currentsize;n++) {
					xyz = narrowband.id[n];
					
					pt[Z] = Numerics.floor(xyz/(nx*ny));
					pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
					pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
					
					direct[X][xyz] = (float)pt[X];
					direct[Y][xyz] = (float)pt[Y];
					direct[Z][xyz] = (float)pt[Z];
					
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
					inverse[T][xyz] = 0.0f;
				}
				
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
    
    /** 
    *  	Evolution using the narrow band scheme 
    *	(the reinitialization is incorporated)
    */
    public final void evolveNarrowBandMappingMean(int iter, float mindiff) {
    	if (debug) System.out.print("level set evolution: narrow band + mapping\n");

    	double[] vect = new double[4];
    	double[] pt = new double[3];
    	double delta;
    	//double phi, dx, dy, dz, grad;

		//float[][] newmapping = new float[3][nx*ny*nz];
		float[][] tmpmap = new float[nimg][nx*ny*nz];
		boolean[] newMappedMask = new boolean[nx*ny*nz];
		
		float[][] direct = new float[3][nx*ny*nz];
		float[][] inverse = new float[4][nx*ny*nz];
		float[][] tmpmapping = new float[3][nx*ny*nz];
		BitSet undef = new BitSet(nx*ny*nz);
		
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
			direct[X][xyz] = 0.0f;
			direct[Y][xyz] = 0.0f;
			direct[Z][xyz] = 0.0f;
		}
			
		// evolve until a landmine is closer than minDist of the boundaries
		float diff = 1.0f;
		for (int t=0;t<iter && diff>mindiff;t++) {
			if (debug) System.out.print("iteration "+t+"\n");
        			
			boolean reinit = false;
			int nswap=0;
			double prev;
			
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				// compute the forces from current levelset values, update the narrow band from it
				prev = narrowband.levels[n];
				levelsetForceVector(xyz, vect);
				narrowband.levels[n] -= vect[T];
				
				// change of sign ? bg -> obj
				if (narrowband.levels[n]<0 && prev>0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, OBJ)) {
						// update the segmentation
						segmentation[xyz] = OBJ;
					} else {
						delta = narrowband.levels[n]-lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				} else if (narrowband.levels[n]>0 && prev<0) {
					//if (debug) System.out.print(""+lb);
					nswap++;
					
					// switch labels	
					if (homeomorphicLabeling(xyz, BG)) {
						// do nothing
						segmentation[xyz] = BG;
					} else {
						delta = narrowband.levels[n]+lowlevel;
						vect[X] *= delta/vect[T];
						vect[Y] *= delta/vect[T];
						vect[Z] *= delta/vect[T];
						narrowband.levels[n] = -lowlevel;
					}
					// check for boundary changes in the landmines : force reinitialization
					if (landmines.get(xyz)) reinit = true;
				}
				direct[X][xyz] = -(float)vect[X];
				direct[Y][xyz] = -(float)vect[Y];
				direct[Z][xyz] = -(float)vect[Z];
			}
			
			diff = (nswap/(float)boundarysize);
			if (debug) System.out.print("changed labels: "+nswap+" ("+(diff*100.0f)+" % of boundary)\n");
			
			// for the case there are issues with the last label (same forces for everyone, for instance)
			//if (nswap[nmgdm-1]==0 && nswap[0]>0) reinit=true;
			
			// forward mapping
			float maxerr = 0.0f;
			float meanerr = 0.0f;
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				pt[X] += direct[X][xyz];
				pt[Y] += direct[Y][xyz];
				pt[Z] += direct[Z][xyz];
				
				/* not very effective
				// compose levelset for adjustment : locally? recursively?
				float err = narrowband.levels[n] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				// get levelset gradient at init location
				float dX = ImageInterpolation.linearInterpolationXderivative(levelset, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				float dY = ImageInterpolation.linearInterpolationYderivative(levelset, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				float dZ = ImageInterpolation.linearInterpolationZderivative(levelset, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
				for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
					pt[X] += 0.67f*err/norm2*dX;
					pt[Y] += 0.67f*err/norm2*dY;
					pt[Z] += 0.67f*err/norm2*dZ;
					
					err = narrowband.levels[n] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				}
				if (Numerics.abs(err)>Numerics.abs(maxerr)) maxerr = err;
				meanerr += Numerics.abs(err);
				*/
				tmpmapping[X][xyz] = ImageInterpolation.linearInterpolation(mappingFD[X], (float)pt[X], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				tmpmapping[Y][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Y], (float)pt[Y], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				tmpmapping[Z][xyz] = ImageInterpolation.linearInterpolation(mappingFD[Z], (float)pt[Z], (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
				
				// compose intensity
				if (mapped != null) {
					for (int i=0;i<nimg;i++) {
						tmpmap[i][xyz] = ImageInterpolation.linearInterpolation(mapped[i], mappedMask, 0.0f, (float)pt[X], (float)pt[Y], (float)pt[Z], nx, ny, nz);
						newMappedMask[xyz] = (tmpmap[i][xyz] != 0.0f); 
					}
				}
			}
			
			// compose levelsets for adjustment: globally is best
			undef.clear();
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				float err = narrowband.levels[n] - ImageInterpolation.linearInterpolation(initlevelset, initlevelset[xyz], tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
				float ratio = 1.0f;
				for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
					float err0 = err;
					// get levelset gradient at init location
					float dX = ImageInterpolation.linearInterpolationXderivative(initlevelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float dY = ImageInterpolation.linearInterpolationYderivative(initlevelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float dZ = ImageInterpolation.linearInterpolationZderivative(initlevelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
					
					tmpmapping[X][xyz] += ratio*err/norm2*dX;
					tmpmapping[Y][xyz] += ratio*err/norm2*dY;
					tmpmapping[Z][xyz] += ratio*err/norm2*dZ;
					
					err = narrowband.levels[n] - ImageInterpolation.linearInterpolation(initlevelset,initlevelset[xyz], tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					if (Numerics.abs(err)>Numerics.abs(err0)) {
						tmpmapping[X][xyz] -= ratio*err/norm2*dX;
						tmpmapping[Y][xyz] -= ratio*err/norm2*dY;
						tmpmapping[Z][xyz] -= ratio*err/norm2*dZ;
						ratio *= 0.67f;
						err = err0;
					} else {
						ratio = 1.0f;
					}
				}
				//if (Numerics.abs(err)>0.5f) undef.set(xyz);
				
				if (Numerics.abs(err)>Numerics.abs(maxerr)) maxerr = err;
				meanerr += Numerics.abs(err);
			}
			if (debug) System.out.print("max residual: "+maxerr+"\n");
			if (debug) System.out.print("mean residual: "+(meanerr/narrowband.currentsize)+"\n");
			
			/*
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				float mindiff2 = 1e9f;			 
				for (int k=0;k<6;k++) if (!undef.get(xyz+xoff[k]+yoff[k]+zoff[k])) {
					float diff2 = Numerics.square(tmpmapping[X][xyz] - tmpmapping[X][xyz+xoff[k]+yoff[k]+zoff[k]])
								 +Numerics.square(tmpmapping[Y][xyz] - tmpmapping[Y][xyz+xoff[k]+yoff[k]+zoff[k]])
								 +Numerics.square(tmpmapping[Z][xyz] - tmpmapping[Z][xyz+xoff[k]+yoff[k]+zoff[k]]);
					if (diff2<mindiff2) mindiff2 = diff2;	
				}
				if (mindiff2 > 100.0f) undef.set(xyz);
			}
			*/	
			/*
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				if (!undef.get(xyz)) {
					mappingFD[X][xyz] = tmpmapping[X][xyz];
					mappingFD[Y][xyz] = tmpmapping[Y][xyz];
					mappingFD[Z][xyz] = tmpmapping[Z][xyz];
				} else {
					pt[X] = mappingFD[X][xyz];
					pt[Y] = mappingFD[Y][xyz];
					pt[Z] = mappingFD[Z][xyz];
					mappingFD[X][xyz] = 0.0f;
					mappingFD[Y][xyz] = 0.0f;
					mappingFD[Z][xyz] = 0.0f;
					float nb=0.0f;
					for (int k=0;k<6;k++) if (!undef.get(xyz+xoff[k]+yoff[k]+zoff[k])) {
						mappingFD[X][xyz] += tmpmapping[X][xyz+xoff[k]+yoff[k]+zoff[k]];
						mappingFD[Y][xyz] += tmpmapping[Y][xyz+xoff[k]+yoff[k]+zoff[k]];
						mappingFD[Z][xyz] += tmpmapping[Z][xyz+xoff[k]+yoff[k]+zoff[k]];
						nb++;
					}
					if (nb==0) {
						// unchanged
						mappingFD[X][xyz] = (float)pt[X];
						mappingFD[Y][xyz] = (float)pt[Y];
						mappingFD[Z][xyz] = (float)pt[Z];
					} else {
						mappingFD[X][xyz] /= nb;
						mappingFD[Y][xyz] /= nb;
						mappingFD[Z][xyz] /= nb;
					}
				}
			}
			*/
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingFD[X][xyz] = tmpmapping[X][xyz];
				mappingFD[Y][xyz] = tmpmapping[Y][xyz];
				mappingFD[Z][xyz] = tmpmapping[Z][xyz];
			}
			// once all the new values are computed, copy into original MGDM functions
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				levelset[xyz] = narrowband.levels[n];
			}
			
			// intensity mapping
			if (mapped != null) {
				for (int i=0;i<nimg;i++) {
					for (int n=0; n<narrowband.currentsize;n++) {
						xyz = narrowband.id[n];
					
						mapped[i][xyz] = tmpmap[i][xyz];
					}
				}
			}
			mappedMask = Arrays.copyOf(newMappedMask, nx*ny*nz);
			
			// build the inverse transform from the direct one: get the "most representative" vector?
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				
				pt[Z] = Numerics.floor(xyz/(nx*ny));
				pt[Y] = Numerics.floor((xyz-pt[Z]*nx*ny)/nx);
				pt[X] = Numerics.floor(xyz-pt[Z]*nx*ny-pt[Y]*nx);
				
				pt[X] += direct[X][xyz];
				pt[Y] += direct[Y][xyz];
				pt[Z] += direct[Z][xyz];
				
				if (Numerics.floor(pt[X])>=0 && Numerics.floor(pt[X])<nx-1 && Numerics.floor(pt[Y])>=0 && Numerics.floor(pt[Y])<ny-1 && Numerics.floor(pt[Z])>=0 && Numerics.floor(pt[Z])<nz-1) {
					int xyz0 = Numerics.floor(pt[X])+Numerics.floor(pt[Y])*nx+Numerics.floor(pt[Z])*nx*ny;
					double wx = pt[X]-Numerics.floor(pt[X]);
					double wy = pt[Y]-Numerics.floor(pt[Y]);
					double wz = pt[Z]-Numerics.floor(pt[Z]);
					inverse[X][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[X][xyz]);
					inverse[Y][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[Y][xyz]);
					inverse[Z][xyz0] += -(float)( (1.0-wx)*(1.0-wy)*(1.0-wz)*direct[Z][xyz]);
					inverse[T][xyz0] += (float)( (1.0-wx)*(1.0-wy)*(1.0-wz));
					inverse[X][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[X][xyz]);
					inverse[Y][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[Y][xyz]);
					inverse[Z][xyz0+1] += -(float)( wx*(1.0-wy)*(1.0-wz)*direct[Z][xyz]);
					inverse[T][xyz0+1] += (float)( wx*(1.0-wy)*(1.0-wz));
					inverse[X][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[X][xyz]);
					inverse[Y][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[Y][xyz]);
					inverse[Z][xyz0+nx] += -(float)( (1.0-wx)*wy*(1.0-wz)*direct[Z][xyz]);
					inverse[T][xyz0+nx] += (float)( (1.0-wx)*wy*(1.0-wz));
					inverse[X][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[X][xyz]);
					inverse[Y][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[Y][xyz]);
					inverse[Z][xyz0+nx*ny] += -(float)( (1.0-wx)*(1.0-wy)*wz*direct[Z][xyz]);
					inverse[T][xyz0+nx*ny] += (float)( (1.0-wx)*(1.0-wy)*wz);
					inverse[X][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[X][xyz]);
					inverse[Y][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[Y][xyz]);
					inverse[Z][xyz0+1+nx] += -(float)( wx*wy*(1.0-wz)*direct[Z][xyz]);
					inverse[T][xyz0+1+nx] += (float)( wx*wy*(1.0-wz));
					inverse[X][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[X][xyz]);
					inverse[Y][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[Y][xyz]);
					inverse[Z][xyz0+nx+nx*ny] += -(float)( (1.0-wx)*wy*wz*direct[Z][xyz]);
					inverse[T][xyz0+nx+nx*ny] +=  (float)( (1.0-wx)*wy*wz);
					inverse[X][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[X][xyz]);
					inverse[Y][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[Y][xyz]);
					inverse[Z][xyz0+nx*ny+1] += -(float)( wx*(1.0-wy)*wz*direct[Z][xyz]);
					inverse[T][xyz0+nx*ny+1] += (float)( wx*(1.0-wy)*wz);
					inverse[X][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[X][xyz]);
					inverse[Y][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[Y][xyz]);
					inverse[Z][xyz0+1+nx+nx*ny] += -(float)( wx*wy*wz*direct[Z][xyz]);
					inverse[T][xyz0+1+nx+nx*ny] += (float)( wx*wy*wz);
				}
			}
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				if (inverse[T][xyz]>0.001) {
					inverse[X][xyz] = inverse[X][xyz]/inverse[T][xyz];
					inverse[Y][xyz] = inverse[Y][xyz]/inverse[T][xyz];
					inverse[Z][xyz] = inverse[Z][xyz]/inverse[T][xyz];
				} else {
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
				}
				inverse[T][xyz] = 0.0f;
			}
			// backward mapping
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				tmpmapping[X][xyz] = mappingBD[X][xyz] + ImageInterpolation.linearInterpolation(inverse[X], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Y][xyz] = mappingBD[Y][xyz] + ImageInterpolation.linearInterpolation(inverse[Y], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
				tmpmapping[Z][xyz] = mappingBD[Z][xyz] + ImageInterpolation.linearInterpolation(inverse[Z], 0.0f, mappingBD[X][xyz], mappingBD[Y][xyz], mappingBD[Z][xyz], nx, ny, nz);
			}
			// compose levelsets for adjustment: globally?
			float maxinverr = 0.0f;
			float meaninverr = 0.0f;
			undef.clear();
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
				float err = initlevelset[xyz] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
				float ratio = 1.0f;
				for (int e=0;e<50 && Numerics.abs(err)>0.01f;e++) {
					float err0 = err;
					// get levelset gradient at init location
					float dX = ImageInterpolation.linearInterpolationXderivative(levelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float dY = ImageInterpolation.linearInterpolationYderivative(levelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float dZ = ImageInterpolation.linearInterpolationZderivative(levelset, tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					float norm2 = Numerics.max(1e-3f,dX*dX+dY*dY+dZ*dZ);
					
					tmpmapping[X][xyz] += 0.67f*err/norm2*dX;
					tmpmapping[Y][xyz] += 0.67f*err/norm2*dY;
					tmpmapping[Z][xyz] += 0.67f*err/norm2*dZ;
					
					err = initlevelset[xyz] - ImageInterpolation.linearInterpolation(levelset, levelset[xyz], tmpmapping[X][xyz], tmpmapping[Y][xyz], tmpmapping[Z][xyz], nx, ny, nz);
					if (Numerics.abs(err)>Numerics.abs(err0)) {
						tmpmapping[X][xyz] -= ratio*err/norm2*dX;
						tmpmapping[Y][xyz] -= ratio*err/norm2*dY;
						tmpmapping[Z][xyz] -= ratio*err/norm2*dZ;
						ratio *= 0.67f;
						err = err0;
					} else {
						ratio = 1.0f;
					}
				}
				//if (Numerics.abs(err)>0.5f) undef.set(xyz);

				if (Numerics.abs(err)>Numerics.abs(maxinverr)) maxinverr = err;
				meaninverr += Numerics.abs(err);
			}
			if (debug) System.out.print("max inv residual: "+maxinverr+"\n");
			if (debug) System.out.print("mean inv residual: "+(meaninverr/narrowband.currentsize)+"\n");
			/*
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				float diff2 = Numerics.square(mappingFD[X][xyz] - tmpmapping[X][xyz])
							 +Numerics.square(mappingFD[Y][xyz] - tmpmapping[Y][xyz])
							 +Numerics.square(mappingFD[Z][xyz] - tmpmapping[Z][xyz]);
				if (!undef.get(xyz) && diff2 < 100.0f) {
					mappingBD[X][xyz] = tmpmapping[X][xyz];
					mappingBD[Y][xyz] = tmpmapping[Y][xyz];
					mappingBD[Z][xyz] = tmpmapping[Z][xyz];
				} else {
					pt[X] = mappingBD[X][xyz];
					pt[Y] = mappingBD[Y][xyz];
					pt[Z] = mappingBD[Z][xyz];
					mappingBD[X][xyz] = 0.0f;
					mappingBD[Y][xyz] = 0.0f;
					mappingBD[Z][xyz] = 0.0f;
					float nb=0.0f;
					for (int k=0;k<6;k++) if (!undef.get(xyz+xoff[k]+yoff[k]+zoff[k])) {
						mappingBD[X][xyz] += tmpmapping[X][xyz+xoff[k]+yoff[k]+zoff[k]];
						mappingBD[Y][xyz] += tmpmapping[Y][xyz+xoff[k]+yoff[k]+zoff[k]];
						mappingBD[Z][xyz] += tmpmapping[Z][xyz+xoff[k]+yoff[k]+zoff[k]];
						nb++;
					}
					if (nb==0) {
						mappingBD[X][xyz] = (float)pt[X];
						mappingBD[Y][xyz] = (float)pt[Y];
						mappingBD[Z][xyz] = (float)pt[Z];
					} else {
						mappingBD[X][xyz] /= nb;
						mappingBD[Y][xyz] /= nb;
						mappingBD[Z][xyz] /= nb;
					}
				}
			}
			*/
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				mappingBD[X][xyz] = tmpmapping[X][xyz];
				mappingBD[Y][xyz] = tmpmapping[Y][xyz];
				mappingBD[Z][xyz] = tmpmapping[Z][xyz];
			}
			// reinit the inverse values
			for (int n=0; n<narrowband.currentsize;n++) {
				xyz = narrowband.id[n];
			
				inverse[X][xyz] = 0.0f;
				inverse[Y][xyz] = 0.0f;
				inverse[Z][xyz] = 0.0f;
				inverse[T][xyz] = 0.0f;				
			}
			
			if (reinit) {
				if (debug) System.out.print("re-initialization\n");
        		
				resetIsosurfaceNarrowBand(narrowband);
				fastMarchingReinitialization(true);
				
				// clean up the transform in the old narrow band
				for (int n=0; n<narrowband.currentsize;n++) {
					xyz = narrowband.id[n];
					
					direct[X][xyz] = 0.0f;
					direct[Y][xyz] = 0.0f;
					direct[Z][xyz] = 0.0f;
					
					inverse[X][xyz] = 0.0f;
					inverse[Y][xyz] = 0.0f;
					inverse[Z][xyz] = 0.0f;
					inverse[T][xyz] = 0.0f;
				}
				
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
    private final void levelsetForceVector(int xyz, double[] vect) {
    	
    	xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
    	xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
    	xyzpz = xyz+nx*ny;	if (!mask[xyzpz]) xyzpz = xyz;
    	xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
    	xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
    	xyzmz = xyz-nx*ny;	if (!mask[xyzmz]) xyzmz = xyz;
    	
    	xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
    	xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
    	xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
    	xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
    	
    	xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
    	xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
    	xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
    	xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
    	
    	xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
    	xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
    	xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
    	xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;

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
		smooth = smoothweight*stepsize*(K-K0)*GPhi;
		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
								 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
								 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
		
		DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
								 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
								 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
		
		balloonforce = Numerics.bounded( trglevelset[xyz] - levelset[xyz], -1.0, 1.0);
		
		balloon = balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*DeltaP + Numerics.min(balloonforce, 0.0)*DeltaM);
		
		// levelset value
		vect[3] = Numerics.bounded(- balloon - smooth,-0.9,0.9);
		
		// coordinate transform: the levelset "force"
		vect[X] = vect[3]*D0x/(GPhi*GPhi);
		vect[Y] = vect[3]*D0y/(GPhi*GPhi);
		vect[Z] = vect[3]*D0z/(GPhi*GPhi);
		
		return;
    }
    
    
	/** specific forces applied to the level sets (application dependent) */
    private final double levelsetForces(int xyz) {
    	
    	xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
    	xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
    	xyzpz = xyz+nx*ny;	if (!mask[xyzpz]) xyzpz = xyz;
    	xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
    	xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
    	xyzmz = xyz-nx*ny;	if (!mask[xyzmz]) xyzmz = xyz;
    	
    	xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
    	xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
    	xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
    	xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
    	
    	xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
    	xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
    	xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
    	xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
    	
    	xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
    	xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
    	xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
    	xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;

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
		smooth = smoothweight*stepsize*(K-K0)*GPhi;
		
		// external object-specific balloon forces
		// we assume b > 0 inside objects, b<0 outside, as in membership functions
		DeltaP = Math.sqrt(Numerics.square(Numerics.max(Dmx, 0.0)) + Numerics.square(Numerics.min(Dpx, 0.0))
								 +Numerics.square(Numerics.max(Dmy, 0.0)) + Numerics.square(Numerics.min(Dpy, 0.0))
								 +Numerics.square(Numerics.max(Dmz, 0.0)) + Numerics.square(Numerics.min(Dpz, 0.0)));
		
		DeltaM = Math.sqrt(Numerics.square(Numerics.max(Dpx, 0.0)) + Numerics.square(Numerics.min(Dmx, 0.0))
								 +Numerics.square(Numerics.max(Dpy, 0.0)) + Numerics.square(Numerics.min(Dmy, 0.0))
								 +Numerics.square(Numerics.max(Dpz, 0.0)) + Numerics.square(Numerics.min(Dmz, 0.0)));
		
		balloonforce = Numerics.bounded( trglevelset[xyz] - levelset[xyz], -1.0, 1.0);
		
		balloon = balloonweight*stepsize*(Numerics.max(balloonforce, 0.0)*DeltaP + Numerics.min(balloonforce, 0.0)*DeltaM);
		
		return - balloon - smooth;
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
    
    public final void computeAverageCurvature() {

    	double avgK = 0.0;
    	double length = 0.0;
    	for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		
    		if (Numerics.abs(levelset[xyz])<1.0) {
    			
				xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
				xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
				xyzpz = xyz+nx*ny;	if (!mask[xyzpz]) xyzpz = xyz;
				xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
				xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
				xyzmz = xyz-nx*ny;	if (!mask[xyzmz]) xyzmz = xyz;
				
				xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
				xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
				xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
				xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
				
				xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
				xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
				xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
				xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
				
				xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
				xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
				xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
				xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
		
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
				
				// combine with simple estimate of pv:
				
				avgK += (1.0-Numerics.abs(levelset[xyz]))*K;
				length += (1.0-Numerics.abs(levelset[xyz]));
			}
		}
		K0 = avgK/length;
		if (debug) System.out.print("Avg Curvature: "+K0+"\n");
    }
    
   public final void computeCurvatureLabeling(float scale) {
   	   
   	   labeling = new byte[nx*ny*nz];
   	   
   	   float[] smoothed = ImageFilters.separableMaskedConvolution(levelset, mask, nx, ny, nz, 
    															ImageFilters.separableGaussianKernel(scale,scale,scale));
   	   
   	   for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
    		
    		if (Numerics.abs(smoothed[xyz])<1.0) {
    			
				xyzpx = xyz+1;		if (!mask[xyzpx]) xyzpx = xyz;
				xyzpy = xyz+nx;		if (!mask[xyzpy]) xyzpy = xyz;
				xyzpz = xyz+nx*ny;	if (!mask[xyzpz]) xyzpz = xyz;
				xyzmx = xyz-1;		if (!mask[xyzmx]) xyzmx = xyz;
				xyzmy = xyz-nx;		if (!mask[xyzmy]) xyzmy = xyz;
				xyzmz = xyz-nx*ny;	if (!mask[xyzmz]) xyzmz = xyz;
				
				xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
				xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
				xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
				xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
				
				xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
				xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
				xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
				xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
				
				xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
				xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
				xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
				xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
		
				// first derivatives
				Dmx = smoothed[xyz] - smoothed[xyzmx];
				Dmy = smoothed[xyz] - smoothed[xyzmy];
				Dmz = smoothed[xyz] - smoothed[xyzmz];
				
				Dpx = smoothed[xyzpx] - smoothed[xyz];
				Dpy = smoothed[xyzpy] - smoothed[xyz];
				Dpz = smoothed[xyzpz] - smoothed[xyz];
				
				D0x = (smoothed[xyzpx] - smoothed[xyzmx])/2.0;
				D0y = (smoothed[xyzpy] - smoothed[xyzmy])/2.0;
				D0z = (smoothed[xyzpz] - smoothed[xyzmz])/2.0;
				
				// second derivatives
				Dxx = smoothed[xyzmx] + smoothed[xyzpx] - 2.0*smoothed[xyz];
				Dyy = smoothed[xyzmy] + smoothed[xyzpy] - 2.0*smoothed[xyz];
				Dzz = smoothed[xyzmz] + smoothed[xyzpz] - 2.0*smoothed[xyz];
				
				Dxy = (smoothed[xyzmxmy] + smoothed[xyzpxpy] - smoothed[xyzmxpy] - smoothed[xyzpxmy])/4.0;
				Dyz = (smoothed[xyzmymz] + smoothed[xyzpypz] - smoothed[xyzmypz] - smoothed[xyzpymz])/4.0;
				Dzx = (smoothed[xyzmzmx] + smoothed[xyzpzpx] - smoothed[xyzmzpx] - smoothed[xyzpzmx])/4.0;
					
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
				
				// combine with simple estimate of pv:
				if (K>0) labeling[xyz] = 3;
				else if (K<0) labeling[xyz] = 2;
			}
		}
		return;    
   }
    
    public final void computeSmoothedLevelset(float scale, boolean rescale) {
    	// use recursive smoothing on the result levelset (not the previous target)
    	trglevelset = ImageFilters.separableMaskedConvolution(levelset, mask, nx, ny, nz, 
    															ImageFilters.separableGaussianKernel(scale,scale,scale));
    	
    	if (rescale) {
			// label the intersections
			labeling = new byte[nx*ny*nz];
			for (int xyz = 0; xyz<nx*ny*nz; xyz++) if (mask[xyz]) {
				if (levelset[xyz]>0 && trglevelset[xyz]>0) labeling[xyz] = 0;
				else if (levelset[xyz]>0 && trglevelset[xyz]<=0) labeling[xyz] = 40;
				else if (levelset[xyz]<=0 && trglevelset[xyz]>0) labeling[xyz] = 80;
				else if (levelset[xyz]<=0 && trglevelset[xyz]<=0) labeling[xyz] = 120;
				else labeling[xyz] = 0;
			}
			
			// find best levelset of new mask to maintain approx. constant surface
			float surf0 = 0.0f;
			float x0 = 0.0f, y0 = 0.0f, z0 = 0.0f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (mask[xyz]) {
					// search for neighbors
					boolean boundary=false;
					for (int l=0;l<6 && !boundary;l++) {
						int xyznb = xyz + xoff[l] + yoff[l] + zoff[l];
						if ( (mask[xyznb]) && (initlevelset[xyz]*initlevelset[xyznb]<=0.0f) ) {
							boundary=true;
						}
					}
					if (boundary) {
						surf0++;
						x0 += x;
						y0 += y;
						z0 += z;
					}
				}
			}
			x0 /= surf0;
			y0 /= surf0;
			z0 /= surf0;
			
			int nsurf = 10;
			float[] surf = new float[nsurf];
			for (int n=0;n<nsurf;n++) surf[n] = 0.0f;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (mask[xyz]) {
					for (int n=0;n<nsurf;n++) {
						float factor = 1.0f + 0.1f*(n+1);
						// scale the intensity
						int xyzT = Numerics.round(x0+(x-x0)/factor) + nx*Numerics.round(y0+(y-y0)/factor) + nx*ny*Numerics.round(z0+(z-z0)/factor);
						if (mask[xyzT]) {
							// search for neighbors
							boolean boundary=false;
							for (int l=0;l<6 && !boundary;l++) {
								int xyznb = xyzT + xoff[l] + yoff[l] + zoff[l];
								if ( (mask[xyznb]) && (trglevelset[xyzT]*trglevelset[xyznb]<=0.0f) ) {
									boundary=true;
								}
							}
							if (boundary) surf[n]++;
						}
					}
				}
			}
			
			int lbest = -1;
			for (int l=0;l<nsurf && lbest==-1;l++) if (surf[l]>surf0) lbest=l;
			if (lbest>=nsurf) lbest = nsurf-1;
			
			System.out.println("inflated surface area: "+surf0);
			for (int l=0;l<nsurf;l++) System.out.print(l+":"+surf[l]+",	");
			System.out.print("\n");
			
			float factor;
			if (lbest>0) factor = 1.0f + 0.1f*(lbest + (surf0-surf[lbest-1])/(surf[lbest]-surf[lbest-1]) );
			else factor = 1.0f + 0.1f*(lbest+1);
			
			System.out.println("closest levelset (level, input / smoothed): "+lbest+", "+surf0+", "+surf[lbest]);
			float[] tmp = new float[nx*ny*nz];
			boolean[] tmpmask = new boolean[nx*ny*nz];
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				int xyz = x+nx*y+nx*ny*z;
				// update the smoothed data and the mask
				tmp[xyz] = ImageInterpolation.linearInterpolation(trglevelset,  narrowBandDist,
								x0+(x-x0)/factor, y0+(y-y0)/factor, z0+(z-z0)/factor, nx, ny, nz);
				
				// most conservative: keep everything (must avoid the boundary, though!)
				if (mask[xyz] || ImageInterpolation.nearestNeighborInterpolation(mask,  false,
								x0+(x-x0)/factor, y0+(y-y0)/factor, z0+(z-z0)/factor, nx, ny, nz) ) 
					tmpmask[xyz] = true;
			}
			trglevelset = tmp;
			mask = tmpmask;
    	}
		return;
    }

}

