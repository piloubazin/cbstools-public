package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

import javax.vecmath.*;
//import WildMagic.LibFoundation.Mathematics.Vector3f;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm handles the creation of levelset avglevelsets
 *
 *	@version    Mar 2015
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class AverageGdm {
	
	// object types
	private	static	final	byte	EMPTY = -1;
	private	static	final	byte	OBJ = 1;
	private	static	final	byte	BG = 0;
	
	// fast marching flags
	private final static byte X = 0;
    private final static byte Y = 1;
    private final static byte Z = 2;
    	
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
	private 	float[][][] 		avglevelset;  			// avglevelset level set function
	private 	float[][][] 		stdlevelset;  			// avglevelset level set function
	private 	float[][][] 		minlevelset;  			// minimum level set function
	private 	float[][][][] 		levelsets;  		// starting group of levelsets
	private		boolean[][][]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	int			nsubj;
	private static	float 		rx,ry,rz;   		// images resolutions
	private		float		distance;			// max distance at which to process the levelsets
	private 	int 		ncx,ncy,ncz;   		// images dimensions
	private 	int 		x0,y0,z0;   		// images dimensions
	private 	int 		xN,yN,zN;   		// images dimensions
	private		BinaryHeap4D	heap;				// the heap used in fast marching
	private		CriticalPointLUT	lut;				// the LUT for critical points
	private		boolean				checkComposed;		// check if the objects are well-composed too (different LUTs)
	private		boolean				checkTopology;		// check if the objects are well-composed too (different LUTs)
	private		float			meanvoloffset = 0.0f;
	
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public AverageGdm(float[][][][] lvlsetin_, int nsubj_, float dist_, int nx_, int ny_, int nz_, float rx_, float ry_, float rz_, String connectivityType_) {
			
		levelsets = lvlsetin_;
		nsubj = nsubj_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
				
		distance = dist_;
		UNKNOWN = distance+1;
		/*
		x0 = 0, xN = nx-1; ncx = nx;
		y0 = 0, yN = ny-1; ncy = ny;
		z0 = 0, zN = nz-1; ncz = nz;
		*/
		findBoundingBox(distance+1);
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, 1, -1, 0, 0};
		zoff = new int[]{0, 0, 0, 0, 1, -1};
		
		// init all the arrays
		try {
			avglevelset = new float[nx][ny][nz];
			stdlevelset = new float[nx][ny][nz];
			minlevelset = new float[nx][ny][nz];
			mask = new boolean[nx][ny][nz];
			// initalize the heap too so we don't have to do it multiple times
			heap = new BinaryHeap4D(ncx*ncy+ncy*ncz+ncz*ncx, BinaryHeap4D.MINTREE);
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
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[x][y][z] = false;
			else mask[x][y][z] = true;
		}
		if (debug) BasicInfo.displayMessage("initialization\n");
	}
		
	public void finalize() {
		levelsets = null;
		avglevelset = null;
		stdlevelset = null;
		minlevelset = null;
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
	
	public final float[][][] getAvgLevelSet() { return avglevelset; }
	
	public final float[][][] getStdLevelSet() { return stdlevelset; }
	
	public final float[] export1DAvgLevelset() {
		float[] res = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[xyz] = avglevelset[x][y][z];
		}
		return res;
	}
	
 	public final float[][][] exportAvgLevelset() {
		float[][][] res = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			res[x][y][z] = avglevelset[x][y][z];
		}
		return res;
	}
	
	public final void findBoundingBox(float maxdist) {
		// crop the data to a reasonable size (must be at least as big as the union of objects)
		maxdist = Numerics.max(0.0f, maxdist);
		
		System.out.println("compute cropping parameters");
		x0 = nx; xN = 0;
		y0 = ny; yN = 0;
		z0 = nz; zN = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boolean keep=false;
			for (int n=0;n<nsubj;n++) if (levelsets[n][x][y][z] < maxdist) keep = true;
			if (keep) {
				x0 = Numerics.min(x0, x);
				xN = Numerics.max(xN, x);
				y0 = Numerics.min(y0, y);
				yN = Numerics.max(yN, y);
				z0 = Numerics.min(z0, z);
				zN = Numerics.max(zN, z);
			}
		}
		ncx = xN-x0+1;
		ncy = yN-y0+1;
		ncz = zN-z0+1;
		
		return;
	}
	
	public final void simpleAverageLevelset() {
		System.out.println("compute avg levelset");
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			avglevelset[x][y][z] = 0.0f;
			for (int n=0;n<nsubj;n++) avglevelset[x][y][z] += levelsets[n][x][y][z]/(float)nsubj;
		}
	}
	
	public final void minimumLevelset() {
		System.out.println("compute min levelset");
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			minlevelset[x][y][z] = levelsets[0][x][y][z];
			for (int n=1;n<nsubj;n++) {
				if (levelsets[n][x][y][z]<minlevelset[x][y][z]) {
					minlevelset[x][y][z] = levelsets[n][x][y][z];
				}
			}
		}
	}
	
	public final void stdevLevelset() {
		System.out.println("compute levelset stdev");
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			stdlevelset[x][y][z] = 0.0f;
			for (int n=0;n<nsubj;n++) stdlevelset[x][y][z] += Numerics.square(levelsets[n][x][y][z]-avglevelset[x][y][z])/(float)nsubj;
			stdlevelset[x][y][z] = (float)FastMath.sqrt(stdlevelset[x][y][z]);
		}
	}
	
	public final float[][][] generateProbaLevelset() {
		System.out.println("compute cumulative probability");
		float[][][] proba = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			NormalDistribution gaussian = new NormalDistribution(avglevelset[x][y][z], stdlevelset[x][y][z]);
			proba[x][y][z] = (float)gaussian.cumulativeProbability(0.0);
		}
		return proba;
	}
	
	public final void normalizeToMeanVolume() {
		// compute mean volume to normalize the avglevelset to the same size
		double meanvol = 0.0;
		float mindist = 1e6f;
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) if (mask[x][y][z]) {
			for (int n=0;n<nsubj;n++) if (levelsets[n][x][y][z]<0) meanvol += 1.0/nsubj;
			if (avglevelset[x][y][z]<mindist) mindist = avglevelset[x][y][z];
		}
		System.out.println("mean volume (voxels): "+meanvol);
		System.out.println("minimum avg. distance: "+mindist);
		
		// find appropriate threshold to have correct volume; should use a fast marching approach!
		heap.reset();
        heap.setMinTree();
		
		double vol = 0.0;
		boolean[][][] label = new boolean[nx][ny][nz];
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) if (mask[x][y][z]) {
			if (avglevelset[x][y][z]<0 || avglevelset[x][y][z]==mindist) {
				vol++;
				label[x][y][z] = true;
				
				// check for neighbors outside
				for (int k = 0; k<6 ; k++) {
					if (avglevelset[x+xoff[k]][y+yoff[k]][z+zoff[k]]>=0) {
						// add to the heap (multiple instances are OK)
						heap.addValue(avglevelset[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
					}
				}
			} else {
				label[x][y][z] = false;
			}
		}
		// run until the volume exceeds the mean volume
		float threshold = 0.0f;
		while (heap.isNotEmpty() && vol<meanvol) {
			threshold = heap.getFirst();
			int x = heap.getFirstX();
			int y = heap.getFirstY();
			int z = heap.getFirstZ();
			heap.removeFirst();	
			if (label[x][y][z]==false) {
				vol++;
				label[x][y][z] = true;
				// add neighbors
				for (int k = 0; k<6; k++) {
					if (label[x+xoff[k]][y+yoff[k]][z+zoff[k]]==false) {
						heap.addValue(avglevelset[x+xoff[k]][y+yoff[k]][z+zoff[k]],x+xoff[k],y+yoff[k],z+zoff[k],1);
					}
				}
			}
		}
		System.out.println("Distance offset: "+threshold);
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			avglevelset[x][y][z] -= threshold;
		}
		meanvoloffset = threshold;
	}

	public final void removeMeanVolumeOffset() {
     	System.out.println("Distance offset: "+meanvoloffset);
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			avglevelset[x][y][z] += meanvoloffset;
		}
	}
	
     /**
      *		perform joint reinitialization for all labels 
      */
     public final void fastMarchingRebuild(float dist0) {
        // computation variables
        boolean[][][] processed = new boolean[nx][ny][nz]; // note: using a byte instead of boolean for the second pass
		float[][][] recomputed = new float[nx][ny][nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
					        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
        heap.setMaxTree();
		// initialize the heap from dist0
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			processed[x][y][z] = false;
			if (minlevelset[x][y][z]>dist0) {
				recomputed[x][y][z] = avglevelset[x][y][z];
				processed[x][y][z] = true;
			}	
		}
		for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) if (mask[x][y][z]) {
        	// search for boundaries around dist0
        	if (minlevelset[x][y][z]>dist0) {
				for (int k = 0; k<6; k++) {
					int xk = x+xoff[k]; int yk = y+yoff[k]; int zk = z+zoff[k];
					if (minlevelset[xk][yk][zk]<=dist0 && mask[xk][yk][zk]) {
						// add to the heap with previous value
						heap.addValue(avglevelset[xk][yk][zk],xk, yk, zk, OBJ);
					}
				}
            }
        }
		if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float maxdist = 0.0f;
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
			int x = heap.getFirstX();
			int y = heap.getFirstY();
			int z = heap.getFirstZ();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[x][y][z])  continue;
			
			// compute new distance based on processed neighbors for the same object
			for (int l=0; l<6; l++) {
				nbdist[l] = UNKNOWN;
				nbflag[l] = false;
				int xl = x+xoff[l]; int yl = y+yoff[l]; int zl = z+zoff[l];
				// note that there is at most one value used here
				if (processed[xl][yl][zl]) if (mask[xl][yl][zl]) {
					nbdist[l] = recomputed[xl][yl][zl];
					nbflag[l] = true;
				}			
			}
			float newdist = minimumNegativeMarchingDistance(nbdist, nbflag);
			
			// update the distance functions at the current level
			recomputed[x][y][z] = newdist;
			processed[x][y][z] = true;
 			
			// find new neighbors
			for (int k = 0; k<6; k++) {
				int xk = x+xoff[k]; int yk = y+yoff[k]; int zk = z+zoff[k];
					
				if (mask[xk][yk][zk]) {
					// must be in outside the object or its processed neighborhood
					if (!processed[xk][yk][zk]) {
						// add to the heap
						heap.addValue(avglevelset[xk][yk][zk],xk,yk,zk,OBJ);
					}
				}
			}
		}
		// replace average
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			avglevelset[x][y][z] = recomputed[x][y][z];
		}
		
		if (debug) BasicInfo.displayMessage("done\n");		
		
       return;
     }

	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * 
     */
    public static final float minimumNegativeMarchingDistance(float[] val, boolean[] flag) {

        double s, s2; // s = a + b +c; s2 = a*a + b*b +c*c
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
        
        double check = (s*s-count*(s2-1.0f));
        if(check <= 0) {
        	check = 0.0001;
        }
        
        tmp = (float)(s-Math.sqrt((double) check))/count;

        // The larger root
        return tmp;
    }
   
    
}

