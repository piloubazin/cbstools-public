package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.util.FastMath;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm builds a hierarchical clustering of profile data
 *
 *	The key component is an array structure that contains the association weights for all pairs of connected points.
 *	By default, we use 6-connectivity (can be easily changed if desired). The array stores the self-weights for each
 *	cluster as well as the degree in the first element of the list. Once constructed, the association array is fully
 *	abstract (no contextual information remains).
 *
 *	@version    March 2012
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class FastAggregativeClustering {
		
	// data buffers
	private 	float[][]	image;  			// original data
	private 	int			nx,ny,nz,nc;   		// image dimensions
	private 	float		rx,ry,rz;   		// image resolutions
	private		boolean[]	mask;
	private		float		imgscale;
	private		float		basis;
	private		int			connect;
	//private		String[]	modes = {"profiles", "binary", "scalar", "noise"};
	private		float 		pvalue;
	
	private		byte		mode;
	private	static final	byte	BINARY = 2;
	private	static final	byte	SCALAR = 3;
	private	static final	byte	GDIST = 4;
	private	static final	byte	NONPARAM = 5;
	private	static final	byte	AUTOGDIST = 6;
	private	static final	byte	ROBUSTGDIST = 7;
	private	static final	byte	HNDIST = 8;
	private	static final	byte	PROFILES_DIST = 10;
	private	static final	byte	PROFILES_GDIST = 11;
	private	static final	byte	PROFILES_CORR = 12;
	private	static final	byte	PROFILES_ANGL = 13;
	private	static final	byte	PROFILES_CDIST = 14;
	private	static final	byte	PROFILES_CGDIST = 15;
	private	static final	byte	GBFH = 20;
	private	static final	byte	MAHALANOBIS = 30;
	private	static final	byte	JENSENSHANNON = 40;
	private	static final	byte	HELLINGER = 50;
	private	static final	byte	KULLBACKLEIBLER = 60;
	
	private static final double PI2 = Math.PI/2.0;
	private static final double SQPI = Math.sqrt(Math.PI);
	private static final float INF = 1e12f;
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// algorithm quantities
	private		int[]		labeling;
	private		int[]		clustering;
	private		int			nlb;
	private		int			key;
	
	//private		ArrayList<Float>	degree;
	//private		ArrayList<ArrayList<Triple>>	assoc;
	private		Triple[][]						assoc;
	private		BitSet							active;
	private		BinaryHeapPair					bintree;
	private		float[]							cost;
	private		float[]							self;
	private		float[]							other;
	private		int[]							latest;
	private		int[][]							clusterPos;
	private		int[]							invertlabeling;
	private		float[]							boundaryScore;
	
	private		Histogram		distribution;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// simple id, value class
	private static class Triple {
		public int id;
		public float weight;
		public float delta;
		public float size;
		
		public Triple(int id) {
			this.id = id;
			this.weight = 0.0f;
			this.delta = 0.0f;
			this.size = 1;
		}
		public Triple(int id, float weight) {
			this.id = id;
			this.weight = weight;
			this.delta = 0.0f;
			this.size = 1;
		}
		public Triple(int id, float weight, float size) {
			this.id = id;
			this.weight = weight;
			this.delta = 0.0f;
			this.size = size;
		}
		public Triple(int id, float weight, float delta, float size) {
			this.id = id;
			this.weight = weight;
			this.delta = delta;
			this.size = size;
		}
		
	}
	
	/* simple variable size array */
	private static class IntArray {
		public int[] val;
		public int	size;
		public int	last;
		
		public IntArray(int size) {
			this.val = new int[size];
			this.size = size;
			this.last = 0;
		}
		public void add(int id) {
			if (last>=val.length) {
				int[] tmp = new int[val.length+size];
				for (int n=0;n<last;n++) tmp[n] = val[n];
				val = tmp;
				if (debug) System.out.print("#");
			}
			val[last] = id;
			last++;
		}
		public void reset() {
			last = 0;
		}
	}
	

	/**
	 *  constructor
	 */
	public FastAggregativeClustering(float[][] img_, boolean[] msk_,
										int nx_, int ny_, int nz_, int nc_,
										float rx_, float ry_, float rz_,
										float ims_, float bas_, int conn_, String mod_) {
		image = img_;
		mask = msk_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		rx = rx_;
		ry = ry_;
		rz = rz_;		
		
		imgscale = ims_;
		basis = bas_;
		connect = conn_;
		pvalue = ims_;
		
		if (mod_.equals("profiles-distance")) 						mode = PROFILES_DIST;
		else if (mod_.equals("profiles-gauss-distance")) 			mode = PROFILES_GDIST;
		else if (mod_.equals("profiles-correlation")) 				mode = PROFILES_CORR;
		else if (mod_.equals("profiles-angle")) 					mode = PROFILES_ANGL;
		else if (mod_.equals("profiles-correlation-dist")) 		mode = PROFILES_CDIST;
		else if (mod_.equals("profiles-correlation-gauss-dist")) 	mode = PROFILES_CGDIST;
		else if (mod_.equals("binary")) 	mode = BINARY;
		else if (mod_.equals("scalar")) 	mode = SCALAR;
		else if (mod_.equals("gaussian")) 	mode = GDIST;
		else if (mod_.equals("auto-gauss")) 	mode = AUTOGDIST;
		else if (mod_.equals("robust-gauss")) 	mode = ROBUSTGDIST;
		else if (mod_.equals("halfnormal")) 	mode = HNDIST;
		else if (mod_.equals("nonparam")) 	mode = NONPARAM;
		else if (mod_.equals("gbfh")) 	mode = GBFH;
		else if (mod_.equals("mahalanobis")) 	mode = MAHALANOBIS;
		else if (mod_.equals("jensen-shannon")) 	mode = JENSENSHANNON;
		else if (mod_.equals("hellinger")) 	mode = HELLINGER;
		else if (mod_.equals("kullback-leibler")) 	mode = KULLBACKLEIBLER;

		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};

		// init: find the labels for used data points
		// pb: destroys spatial information (i.e. other data must be transferred in that space too)
		labeling = new int[nx*ny*nz];
		nlb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				nlb++;
				labeling[xyz] = nlb;
			} else {
				labeling[xyz] = 0;
			}
		}
		invertlabeling = new int[nlb+1];
		invertlabeling[0] = -1;
		int lb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				lb++;
				invertlabeling[lb] = xyz;
			}
		}
		
		if (debug) System.out.println("initial labeling ("+nlb+" voxels)");
		
	}

	final public void finalize() {
		image = null;
		labeling = null;
		assoc = null;
		System.gc();
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsAverage() {
		
		if (debug) System.out.println("-- weight initialization --");
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		Triple[] trp = new Triple[connect];
		
		if (debug) System.out.println("first pass");

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int nb=0;
				// 6-C
				if (x>0 && mask[xyz-1])	{
					trp[nb] = new Triple(labeling[xyz-1], associationWeight(xyz,xyz-1));
					nb++;
				}
				if (x<nx-1 && mask[xyz+1])	{
					trp[nb] = new Triple(labeling[xyz+1], associationWeight(xyz,xyz+1));
					nb++;
				}
				if (y>0 && mask[xyz-nx])	{
					trp[nb] = new Triple(labeling[xyz-nx], associationWeight(xyz,xyz-nx));
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					trp[nb] = new Triple(labeling[xyz+nx], associationWeight(xyz,xyz+nx));
					nb++;
				}
				if (z>0 && mask[xyz-nx*ny])	{
					trp[nb] = new Triple(labeling[xyz-nx*ny], associationWeight(xyz,xyz-nx*ny));
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					trp[nb] = new Triple(labeling[xyz+nx*ny], associationWeight(xyz,xyz+nx*ny));
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y>0 && mask[xyz-1-nx])	{
						trp[nb] = new Triple(labeling[xyz-1-nx], associationWeight(xyz,xyz-1-nx));
						nb++;
					}
					if (x<nx-1 && y>0 && mask[xyz+1-nx])	{
						trp[nb] = new Triple(labeling[xyz+1-nx], associationWeight(xyz,xyz+1-nx));
						nb++;
					}
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						trp[nb] = new Triple(labeling[xyz-1+nx], associationWeight(xyz,xyz-1+nx));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						trp[nb] = new Triple(labeling[xyz+1+nx], associationWeight(xyz,xyz+1+nx));
						nb++;
					}
					if (y>0 && z>0 && mask[xyz-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx-nx*ny], associationWeight(xyz,xyz-nx-nx*ny));
						nb++;
					}
					if (y<ny-1 && z>0 && mask[xyz+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx-nx*ny], associationWeight(xyz,xyz+nx-nx*ny));
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx+nx*ny], associationWeight(xyz,xyz-nx+nx*ny));
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx+nx*ny], associationWeight(xyz,xyz+nx+nx*ny));
						nb++;
					}
					if (z>0 && x>0 && mask[xyz-nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny-1], associationWeight(xyz,xyz-nx*ny-1));
						nb++;
					}
					if (z<nz-1 && x>0 && mask[xyz+nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny-1], associationWeight(xyz,xyz+nx*ny-1));
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny+1], associationWeight(xyz,xyz-nx*ny+1));
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny+1], associationWeight(xyz,xyz+nx*ny+1));
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z>0 && mask[xyz-1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx-nx*ny], associationWeight(xyz,xyz-1-nx-nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z>0 && mask[xyz+1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx-nx*ny], associationWeight(xyz,xyz+1-nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z>0 && mask[xyz-1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx-nx*ny], associationWeight(xyz,xyz-1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx+nx*ny], associationWeight(xyz,xyz-1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z>0 && mask[xyz+1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx-nx*ny], associationWeight(xyz,xyz+1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx+nx*ny], associationWeight(xyz,xyz-1+nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx+nx*ny], associationWeight(xyz,xyz+1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx+nx*ny], associationWeight(xyz,xyz+1+nx+nx*ny));
						nb++;
					}
				}
				// build the weight array

				// link in 6 directions
				Triple[] node = new Triple[nb+1];
				node[0] = new Triple(labeling[xyz], 1.0f);
				for (int n=0;n<nb;n++) {
					node[n+1] = trp[n];
				}
				// build the degree array
				float deg = node[0].weight/node[0].size;
				for (int n=0;n<nb;n++) {
					deg += trp[n].weight/trp[n].size/nb;
				}
				//degree.set(labeling[xyz], new Float(deg));
				node[0].delta = deg;

				// store the values
				assoc[labeling[xyz]] = node;
				
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		for (int l=1;l<=nlb;l++) {
			Triple[] node = assoc[l];
			float di = node[0].delta;
			for (int n=1;n<node.length;n++) {
				//float dj = degree.get(node.get(n).id);
				float wij = node[n].weight;
				float wii = node[0].weight;
				float wjj = assoc[node[n].id][0].weight;
				float dj = assoc[node[n].id][0].delta;
				float sij = node[n].size;
				float si = node[0].size;
				float sj = assoc[node[n].id][0].size;
				//float dval = 2.0f*node.get(n).weight/(di + dj);
				// use the formula with self-weights?
				// D = (wii+wjj+2node)/(di+dj) -wii/di -wjj/dj
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj) - wii/di - wjj/dj;
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj);
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				//node.get(n).delta = 2.0f*wij/sij*(di + dj)/(si + sj);
				
				// best so far..
				node[n].delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				
				// use the same score as the stopping criterion? slows down the process
				//node.get(n).delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij)
				//					- (1.0f- wij/sij)*wii/si*wjj/sj;
				
				//other[0] += wij/sij;
			}
			cost[0] += node[0].weight/node[0].size/node[0].delta;
			//self[0] += node.get(0).weight/node.get(0).size;
		}
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsGBFH() {
		
		if (debug) System.out.println("-- weight initialization (graph-based) --");
		
		//assoc = new ArrayList<ArrayList<Triple>>(nlb+1);
		//degree = new ArrayList(nlb+1);
		//assoc.add(0, new ArrayList<Triple>(1));
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		Triple[] trp = new Triple[connect];
		
		if (debug) System.out.println("first pass");

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int nb=0;
				// 6-C
				if (x>0 && mask[xyz-1])	{
					trp[nb] = new Triple(labeling[xyz-1], associationWeight(xyz,xyz-1));
					nb++;
				}
				if (x<nx-1 && mask[xyz+1])	{
					trp[nb] = new Triple(labeling[xyz+1], associationWeight(xyz,xyz+1));
					nb++;
				}
				if (y>0 && mask[xyz-nx])	{
					trp[nb] = new Triple(labeling[xyz-nx], associationWeight(xyz,xyz-nx));
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					trp[nb] = new Triple(labeling[xyz+nx], associationWeight(xyz,xyz+nx));
					nb++;
				}
				if (z>0 && mask[xyz-nx*ny])	{
					trp[nb] = new Triple(labeling[xyz-nx*ny], associationWeight(xyz,xyz-nx*ny));
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					trp[nb] = new Triple(labeling[xyz+nx*ny], associationWeight(xyz,xyz+nx*ny));
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y>0 && mask[xyz-1-nx])	{
						trp[nb] = new Triple(labeling[xyz-1-nx], associationWeight(xyz,xyz-1-nx));
						nb++;
					}
					if (x<nx-1 && y>0 && mask[xyz+1-nx])	{
						trp[nb] = new Triple(labeling[xyz+1-nx], associationWeight(xyz,xyz+1-nx));
						nb++;
					}
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						trp[nb] = new Triple(labeling[xyz-1+nx], associationWeight(xyz,xyz-1+nx));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						trp[nb] = new Triple(labeling[xyz+1+nx], associationWeight(xyz,xyz+1+nx));
						nb++;
					}
					if (y>0 && z>0 && mask[xyz-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx-nx*ny], associationWeight(xyz,xyz-nx-nx*ny));
						nb++;
					}
					if (y<ny-1 && z>0 && mask[xyz+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx-nx*ny], associationWeight(xyz,xyz+nx-nx*ny));
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx+nx*ny], associationWeight(xyz,xyz-nx+nx*ny));
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx+nx*ny], associationWeight(xyz,xyz+nx+nx*ny));
						nb++;
					}
					if (z>0 && x>0 && mask[xyz-nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny-1], associationWeight(xyz,xyz-nx*ny-1));
						nb++;
					}
					if (z<nz-1 && x>0 && mask[xyz+nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny-1], associationWeight(xyz,xyz+nx*ny-1));
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny+1], associationWeight(xyz,xyz-nx*ny+1));
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny+1], associationWeight(xyz,xyz+nx*ny+1));
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z>0 && mask[xyz-1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx-nx*ny], associationWeight(xyz,xyz-1-nx-nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z>0 && mask[xyz+1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx-nx*ny], associationWeight(xyz,xyz+1-nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z>0 && mask[xyz-1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx-nx*ny], associationWeight(xyz,xyz-1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx+nx*ny], associationWeight(xyz,xyz-1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z>0 && mask[xyz+1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx-nx*ny], associationWeight(xyz,xyz+1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx+nx*ny], associationWeight(xyz,xyz-1+nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx+nx*ny], associationWeight(xyz,xyz+1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx+nx*ny], associationWeight(xyz,xyz+1+nx+nx*ny));
						nb++;
					}
				}
				// build the weight array

				// link in 6 directions
				Triple[] node = new Triple[nb+1];
				node[0] = new Triple(labeling[xyz], 0.0f);
				for (int n=0;n<nb;n++) {
					node[n+1] = trp[n];
					//if (trp[n].weight==0) System.out.print("o");
				}
				// build the degree array
				float deg = node[0].weight/node[0].size;
				for (int n=0;n<nb;n++) {
					deg += trp[n].weight/trp[n].size/nb;
				}
				//degree.set(labeling[xyz], new Float(deg));
				node[0].delta = deg;

				// store the values
				assoc[labeling[xyz]] = node;
				
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		for (int l=1;l<=nlb;l++) {
			Triple[] node = assoc[l];
			float di = node[0].delta;
			for (int n=1;n<node.length;n++) {
				//float dj = degree.get(node.get(n).id);
				float wij = node[n].weight;
				float wii = node[0].weight;
				float wjj = assoc[node[n].id][0].weight;
				float dj = assoc[node[n].id][0].delta;
				float sij = node[n].size;
				float si = node[0].size;
				float sj = assoc[node[n].id][0].size;
				//float dval = 2.0f*node.get(n).weight/(di + dj);
				// use the formula with self-weights?
				// D = (wii+wjj+2node)/(di+dj) -wii/di -wjj/dj
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj) - wii/di - wjj/dj;
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj);
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				//node.get(n).delta = 2.0f*wij/sij*(di + dj)/(si + sj);
				
				node[n].delta = basis - wij;
				
				// use the same score as the stopping criterion? slows down the process
				//node.get(n).delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij)
				//					- (1.0f- wij/sij)*wii/si*wjj/sj;
				
				//other[0] += wij/sij;
			}
			cost[0] += node[0].weight/node[0].size/node[0].delta;
			//self[0] += node.get(0).weight/node.get(0).size;
		}
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsMinimax() {
		
		if (debug) System.out.println("-- weight initialization --");
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		Triple[] trp = new Triple[connect];
		
		if (debug) System.out.println("first pass");

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int nb=0;
				// 6-C
				if (x>0 && mask[xyz-1])	{
					trp[nb] = new Triple(labeling[xyz-1], associationWeight(xyz,xyz-1));
					nb++;
				}
				if (x<nx-1 && mask[xyz+1])	{
					trp[nb] = new Triple(labeling[xyz+1], associationWeight(xyz,xyz+1));
					nb++;
				}
				if (y>0 && mask[xyz-nx])	{
					trp[nb] = new Triple(labeling[xyz-nx], associationWeight(xyz,xyz-nx));
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					trp[nb] = new Triple(labeling[xyz+nx], associationWeight(xyz,xyz+nx));
					nb++;
				}
				if (z>0 && mask[xyz-nx*ny])	{
					trp[nb] = new Triple(labeling[xyz-nx*ny], associationWeight(xyz,xyz-nx*ny));
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					trp[nb] = new Triple(labeling[xyz+nx*ny], associationWeight(xyz,xyz+nx*ny));
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y>0 && mask[xyz-1-nx])	{
						trp[nb] = new Triple(labeling[xyz-1-nx], associationWeight(xyz,xyz-1-nx));
						nb++;
					}
					if (x<nx-1 && y>0 && mask[xyz+1-nx])	{
						trp[nb] = new Triple(labeling[xyz+1-nx], associationWeight(xyz,xyz+1-nx));
						nb++;
					}
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						trp[nb] = new Triple(labeling[xyz-1+nx], associationWeight(xyz,xyz-1+nx));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						trp[nb] = new Triple(labeling[xyz+1+nx], associationWeight(xyz,xyz+1+nx));
						nb++;
					}
					if (y>0 && z>0 && mask[xyz-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx-nx*ny], associationWeight(xyz,xyz-nx-nx*ny));
						nb++;
					}
					if (y<ny-1 && z>0 && mask[xyz+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx-nx*ny], associationWeight(xyz,xyz+nx-nx*ny));
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx+nx*ny], associationWeight(xyz,xyz-nx+nx*ny));
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx+nx*ny], associationWeight(xyz,xyz+nx+nx*ny));
						nb++;
					}
					if (z>0 && x>0 && mask[xyz-nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny-1], associationWeight(xyz,xyz-nx*ny-1));
						nb++;
					}
					if (z<nz-1 && x>0 && mask[xyz+nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny-1], associationWeight(xyz,xyz+nx*ny-1));
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny+1], associationWeight(xyz,xyz-nx*ny+1));
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny+1], associationWeight(xyz,xyz+nx*ny+1));
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z>0 && mask[xyz-1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx-nx*ny], associationWeight(xyz,xyz-1-nx-nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z>0 && mask[xyz+1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx-nx*ny], associationWeight(xyz,xyz+1-nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z>0 && mask[xyz-1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx-nx*ny], associationWeight(xyz,xyz-1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx+nx*ny], associationWeight(xyz,xyz-1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z>0 && mask[xyz+1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx-nx*ny], associationWeight(xyz,xyz+1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx+nx*ny], associationWeight(xyz,xyz-1+nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx+nx*ny], associationWeight(xyz,xyz+1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx+nx*ny], associationWeight(xyz,xyz+1+nx+nx*ny));
						nb++;
					}
				}
				// build the weight array
				
				// link in 6 directions
				Triple[] node = new Triple[nb+1];
				node[0] = new Triple(labeling[xyz], 1.0f);
				for (int n=0;n<nb;n++) {
					node[n+1] = trp[n];
				}
				// build the degree array
				float deg = node[0].weight/node[0].size;
				for (int n=0;n<nb;n++) {
					deg += trp[n].weight/trp[n].size/nb;
				}
				//degree.set(labeling[xyz], new Float(deg));
				node[0].delta = deg;

				// store the values
				assoc[labeling[xyz]] = node;
				
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		for (int l=1;l<=nlb;l++) {
			Triple[] node = assoc[l];
			float di = node[0].delta;
			for (int n=1;n<node.length;n++) {
				//float dj = degree.get(node.get(n).id);
				float wij = node[n].weight;
				float wii = node[0].weight;
				float wjj = assoc[node[n].id][0].weight;
				float dj = assoc[node[n].id][0].delta;
				float sij = node[n].size;
				float si = node[0].size;
				float sj = assoc[node[n].id][0].size;
				
				// best so far..
				//node[n].delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				//node[n].delta = wij*Numerics.min(wii,wjj,wij);
				//node[n].delta = wij*Numerics.min(wii,wjj,wij) - (1.0f-wij)*(float)Math.sqrt(wii*wjj);
				node[n].delta = Numerics.square(wij*Numerics.min(wii,wjj,wij)) - Numerics.square(1.0f-wij)*wii*wjj;
			}
			cost[0] += node[0].weight/node[0].size/node[0].delta;
		}
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsNormalClusters() {
		
		if (debug) System.out.println("-- weight initialization (normal clusters)--");
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		Triple[] trp = new Triple[connect];
		
		if (debug) System.out.println("first pass");

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int nb=0;
				// 6-C
				if (x>0 && mask[xyz-1])	{
					trp[nb] = new Triple(labeling[xyz-1], 0.0f, associationWeight(xyz,xyz-1), 1.0f);
					nb++;
				}
				if (x<nx-1 && mask[xyz+1])	{
					trp[nb] = new Triple(labeling[xyz+1], 0.0f, associationWeight(xyz,xyz+1), 1.0f);
					nb++;
				}
				if (y>0 && mask[xyz-nx])	{
					trp[nb] = new Triple(labeling[xyz-nx], 0.0f, associationWeight(xyz,xyz-nx), 1.0f);
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					trp[nb] = new Triple(labeling[xyz+nx], 0.0f, associationWeight(xyz,xyz+nx), 1.0f);
					nb++;
				}
				if (z>0 && mask[xyz-nx*ny])	{
					trp[nb] = new Triple(labeling[xyz-nx*ny], 0.0f, associationWeight(xyz,xyz-nx*ny), 1.0f);
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					trp[nb] = new Triple(labeling[xyz+nx*ny], 0.0f, associationWeight(xyz,xyz+nx*ny), 1.0f);
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y>0 && mask[xyz-1-nx])	{
						trp[nb] = new Triple(labeling[xyz-1-nx], 0.0f, associationWeight(xyz,xyz-1-nx), 1.0f);
						nb++;
					}
					if (x<nx-1 && y>0 && mask[xyz+1-nx])	{
						trp[nb] = new Triple(labeling[xyz+1-nx], 0.0f, associationWeight(xyz,xyz+1-nx), 1.0f);
						nb++;
					}
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						trp[nb] = new Triple(labeling[xyz-1+nx], 0.0f, associationWeight(xyz,xyz-1+nx), 1.0f);
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						trp[nb] = new Triple(labeling[xyz+1+nx], 0.0f, associationWeight(xyz,xyz+1+nx), 1.0f);
						nb++;
					}
					if (y>0 && z>0 && mask[xyz-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx-nx*ny], 0.0f, associationWeight(xyz,xyz-nx-nx*ny), 1.0f);
						nb++;
					}
					if (y<ny-1 && z>0 && mask[xyz+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx-nx*ny], 0.0f, associationWeight(xyz,xyz+nx-nx*ny), 1.0f);
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-nx+nx*ny], 0.0f, associationWeight(xyz,xyz-nx+nx*ny), 1.0f);
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+nx+nx*ny], 0.0f, associationWeight(xyz,xyz+nx+nx*ny), 1.0f);
						nb++;
					}
					if (z>0 && x>0 && mask[xyz-nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny-1], 0.0f, associationWeight(xyz,xyz-nx*ny-1), 1.0f);
						nb++;
					}
					if (z<nz-1 && x>0 && mask[xyz+nx*ny-1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny-1], 0.0f, associationWeight(xyz,xyz+nx*ny-1), 1.0f);
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz-nx*ny+1], 0.0f, associationWeight(xyz,xyz-nx*ny+1), 1.0f);
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						trp[nb] = new Triple(labeling[xyz+nx*ny+1], 0.0f, associationWeight(xyz,xyz+nx*ny+1), 1.0f);
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z>0 && mask[xyz-1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx-nx*ny], 0.0f, associationWeight(xyz,xyz-1-nx-nx*ny), 1.0f);
						nb++;
					}
					if (x<nx-1 && y>0 && z>0 && mask[xyz+1-nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx-nx*ny], 0.0f, associationWeight(xyz,xyz+1-nx-nx*ny), 1.0f);
						nb++;
					}
					if (x>0 && y<ny-1 && z>0 && mask[xyz-1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx-nx*ny], 0.0f, associationWeight(xyz,xyz-1+nx-nx*ny), 1.0f);
						nb++;
					}
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1-nx+nx*ny], 0.0f, associationWeight(xyz,xyz-1-nx+nx*ny), 1.0f);
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z>0 && mask[xyz+1+nx-nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx-nx*ny], 0.0f, associationWeight(xyz,xyz+1+nx-nx*ny), 1.0f);
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz-1+nx+nx*ny], 0.0f, associationWeight(xyz,xyz-1+nx+nx*ny), 1.0f);
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1-nx+nx*ny], 0.0f, associationWeight(xyz,xyz+1-nx+nx*ny), 1.0f);
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						trp[nb] = new Triple(labeling[xyz+1+nx+nx*ny], 0.0f, associationWeight(xyz,xyz+1+nx+nx*ny), 1.0f);
						nb++;
					}
				}
				// build the weight array

				// link in 6 directions
				Triple[] node = new Triple[nb+1];
				node[0] = new Triple(labeling[xyz], image[0][xyz], 0.0f, 1.0f);
				for (int n=0;n<nb;n++) {
					node[n+1] = trp[n];
				}
				
				// store the values
				assoc[labeling[xyz]] = node;
				
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost: not needed
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		
		return;
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final float[] computeWeightMetrics() {
		
		if (debug) System.out.println("-- metric pre-processing --");
		
		float[] distance = new float[connect/2*nx*ny*nz];
		int nb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				// 6-C
				if (x<nx-1 && mask[xyz+1])	{
					distance[nb] = dataDistance(xyz,xyz+1);
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					distance[nb] = dataDistance(xyz,xyz+nx);
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					distance[nb] = dataDistance(xyz,xyz+nx*ny);
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						distance[nb] = dataDistance(xyz,xyz-1+nx);
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						distance[nb] = dataDistance(xyz,xyz+1+nx);
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz-nx+nx*ny);
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz+nx+nx*ny);
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						distance[nb] = dataDistance(xyz,xyz-nx*ny+1);
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						distance[nb] = dataDistance(xyz,xyz+nx*ny+1);
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz-1-nx+nx*ny);
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz-1+nx+nx*ny);
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz+1-nx+nx*ny);
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						distance[nb] = dataDistance(xyz,xyz+1+nx+nx*ny);
						nb++;
					}
				}
			}
		}
		
		if (debug) System.out.println("statistics");

		Histogram hist = new Histogram(distance, 100, nb);
		System.out.println("Mean : "+hist.mean());
		System.out.println("Stdev : "+hist.stdev());
		System.out.println("Min : "+hist.min());
		System.out.println("Max : "+hist.max());
		System.out.println("Rayleigh : "+hist.rayleighParameter());
		System.out.println("10% : "+hist.percentage(0.10f));
		System.out.println("25% : "+hist.percentage(0.25f));
		System.out.println("50% : "+hist.percentage(0.50f));
		System.out.println("75% : "+hist.percentage(0.75f));
		System.out.println("90% : "+hist.percentage(0.90f));
		
		if (mode==AUTOGDIST) {
			imgscale = hist.mean();
		} else if (mode==ROBUSTGDIST) {
			imgscale = hist.percentage(imgscale);
		} else if (mode==HNDIST) {
			imgscale = (float)(0.5*hist.mean()*SQPI);
		} else if (mode==MAHALANOBIS) {
			imgscale = (float)(0.5*hist.mean()*SQPI);
		} else if (mode==JENSENSHANNON) {
			double sum = 0.0;
			for (int n=0;n<nb;n++) sum += Numerics.square(distance[n]);
			sum /= nb;
			imgscale = (float)FastMath.sqrt(0.5*sum);
		} else if (mode==HELLINGER) {
			imgscale = (float)(0.5*hist.mean()*SQPI);
		} else if (mode==KULLBACKLEIBLER) {
			imgscale = (float)(0.5*hist.mean()*SQPI);
		}
		distribution = hist;
		distribution.normalizeToMaximum();
		//distribution.normalizeToNonZeroMaximum();
		//distribution.inverseCumulativeHistogram();
		//System.out.println(distribution.printHistogram());
		System.out.println("adjusted scale: "+imgscale);
		
		return distance;
	}

	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringAverageAlt(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : average (alt)");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		
		// store cluster location dynamically? (only for new clusters?)
		clusterPos = new int[2*nlb][];

		// look for joint boundaries
		BitSet ngbcluster = new BitSet(nx*ny*nz);
		boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// self-weight & degree
					// d(uv) = d(u) + d(v)
					// w(uv,uv) = w(u,u) + w(v,v) + 2*w(uv.uv)
					aNode[0] = new Triple(id, bNode[0].weight + pNode[0].weight + 2.0f*bNode[nbest].weight, 
												bNode[0].delta + pNode[0].delta,  
												bNode[0].size + pNode[0].size + 2.0f*bNode[nbest].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += bNode[n].weight;
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += pNode[n].weight;
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lpair] = null;
					} else {
						clusterPos[id] = new int[2];
						clusterPos[id][0] = lbest;
						clusterPos[id][1] = lpair;
						// no need to check
						boundaryScore[invertlabeling[lbest]] = bcost;
						boundaryScore[invertlabeling[lpair]] = bcost;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						// use the same criterion than for stopping? makes sense, but slows down the process
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/((connect+1.0f)*nlb)) );
						//float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/((connect+1.0f))) );
						//float wsize = 1.0f;
						aNode[n].delta = wsize*( Numerics.square(aNode[n].weight/aNode[n].size)
													 *Numerics.square( ( aNode[0].weight + wngb[0].weight + 2.0f*aNode[n].weight )
																		/( aNode[0].size + wngb[0].size + 2.0f*aNode[n].size ) )
													 -Numerics.square(1.0f-aNode[n].weight/aNode[n].size)
													  *aNode[0].weight/aNode[0].size
													  *wngb[0].weight/wngb[0].size  );
						
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = 	Numerics.square(bNode[nbest].weight/bNode[nbest].size * aNode[0].weight/aNode[0].size);
					other[iter] = Numerics.square(1.0f-bNode[nbest].weight/bNode[nbest].size) * bNode[0].weight/bNode[0].size
																							* pNode[0].weight/pNode[0].size;
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
					// not a correct stopping criterion								
					if (verbose) if (bcost<0 && first) {
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
					}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);

						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}

					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					assoc[lbest][0] = new Triple(id);
					assoc[lpair][0] = new Triple(id);
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
					
					// reduce the number of clusters
					nclusters--;

					/*
					// update the neighbors too, but only up to a point 
					// (-> use both strategies to gain speed in the two extreme worst case scenarios)
					// (both strategies = this and the track updating)
					for (int b=1;b<aNode.length && b<maxlength;b++) {
						if (active.get(aNode[b].id)) {
							widthsearchtime += assoc[aNode[b].id].length;
							for (int c=1;c<assoc[aNode[b].id].length && c<maxlength;c++) {
								if (latest[assoc[aNode[b].id][c].id]==lbest || latest[assoc[aNode[b].id][c].id]==lpair)
									latest[assoc[aNode[b].id][c].id] = id;
								//if (assoc[aNode[b].id][c].id==lbest || assoc[aNode[b].id][c].id==lpair)
								//	assoc[aNode[b].id][c].id = id;
							}
						}
					}
					*/
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringMinimax(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : minimax");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		// init for max() weights
		//for (int n=0;n<2*nlb;n++) ngbWeight[n] = -INF;
		// init for min() weights
		for (int n=0;n<2*nlb;n++) ngbWeight[n] = INF;
		
		// store cluster location dynamically? (only for new clusters?)
		clusterPos = new int[2*nlb][];

		// look for joint boundaries
		BitSet ngbcluster = new BitSet(nx*ny*nz);
		boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// self-weight & degree
					// d(uv) = d(u) + d(v)
					// w(uv,uv) = w(u,u) + w(v,v) + 2*w(uv.uv)
					aNode[0] = new Triple(id, Numerics.min(bNode[0].weight, pNode[0].weight, bNode[nbest].weight), 
												Numerics.min(bNode[0].delta, pNode[0].delta),  
												bNode[0].size + pNode[0].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//bNode[n].id = lbn;
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								// use the min or max of weights along the path?
								// max corresp. to graph-based, not noise prone
								//ngbWeight[lbn] = Numerics.max(bNode[n].weight, ngbWeight[lbn]);
								// min corresp. to most conservative
								ngbWeight[lbn] = Numerics.min(bNode[n].weight, ngbWeight[lbn]);
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//pNode[n].id = lbn;
								
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								//ngbWeight[lbn] = Numerics.max(pNode[n].weight, ngbWeight[lbn]);
								ngbWeight[lbn] = Numerics.min(pNode[n].weight, ngbWeight[lbn]);
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lpair] = null;
					} else {
						clusterPos[id] = new int[2];
						clusterPos[id][0] = lbest;
						clusterPos[id][1] = lpair;
						// no need to check
						boundaryScore[invertlabeling[lbest]] = bcost;
						boundaryScore[invertlabeling[lpair]] = bcost;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						//aNode[n].delta = aNode[n].weight*Numerics.min(aNode[0].weight, wngb[0].weight, aNode[n].weight);
												
						// use the same criterion than for stopping? makes sense, but slows down the process
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/nlb) );
						//float wsize = 1.0f;
						aNode[n].delta = wsize*( Numerics.square(aNode[n].weight *Numerics.min(aNode[0].weight, wngb[0].weight, aNode[n].weight ))
											     -Numerics.square(1.0f-aNode[n].weight)*aNode[0].weight*wngb[0].weight );
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = Numerics.square(bNode[nbest].weight * aNode[0].weight);
					other[iter] = Numerics.square(1.0f-bNode[nbest].weight) * bNode[0].weight * pNode[0].weight;
					
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
																	
				   if (verbose) if (bcost<0 && first) {
				   	   // that's not the originally computed score, is it??
				   	   System.out.println(iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
				   	}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);

						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}

					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					assoc[lbest][0] = new Triple(id);
					assoc[lpair][0] = new Triple(id);
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
										
					// reduce the number of clusters
					nclusters--;
					
					/*
					// update the neighbors too, but only up to a point 
					// (-> use both strategies to gain speed in the two extreme worst case scenarios)
					// (both strategies = this and the track updating)
					for (int b=1;b<aNode.length && b<maxlength;b++) {
						if (active.get(aNode[b].id)) {
							widthsearchtime += assoc[aNode[b].id].length;
							for (int c=1;c<assoc[aNode[b].id].length && c<maxlength;c++) {
								if (latest[assoc[aNode[b].id][c].id]==lbest || latest[assoc[aNode[b].id][c].id]==lpair)
									latest[assoc[aNode[b].id][c].id] = id;
								//if (assoc[aNode[b].id][c].id==lbest || assoc[aNode[b].id][c].id==lpair)
								//	assoc[aNode[b].id][c].id = id;
							}
						}
					}
					*/
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	// perform the hierarchical clustering with the Felzenszwab-Huttenlocher method
	public final void hierarchicalClusteringGraphBasedFH(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : graph-based FH");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		for (int n=0;n<2*nlb;n++) ngbWeight[n] = INF;
		
		// store cluster location dynamically? (only for new clusters?)
		clusterPos = new int[2*nlb][];

		// look for joint boundaries
		BitSet ngbcluster = new BitSet(nx*ny*nz);
		boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// self-weight & degree
					// d(uv) = d(u) + d(v)
					// w(uv,uv) = w(u,u) + w(v,v) + 2*w(uv.uv)
					aNode[0] = new Triple(id, Numerics.max(bNode[0].weight, pNode[0].weight, bNode[nbest].weight), 
												Numerics.min(bNode[0].delta, pNode[0].delta),  
												bNode[0].size + pNode[0].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//bNode[n].id = lbn;
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] = Numerics.min(bNode[n].weight,ngbWeight[lbn]);
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//pNode[n].id = lbn;
								
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] = Numerics.min(pNode[n].weight,ngbWeight[lbn]);
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = INF;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();
					
					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lpair] = null;
					} else {
						clusterPos[id] = new int[2];
						clusterPos[id][0] = lbest;
						clusterPos[id][1] = lpair;
						// no need to check
						boundaryScore[invertlabeling[lbest]] = bcost;
						boundaryScore[invertlabeling[lpair]] = bcost;
					}

					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// score : min joint weight - max( self-weight + k / size, self-weight + k / size)
						aNode[n].delta = Numerics.min(aNode[0].weight + basis/aNode[0].size, 
														wngb[0].weight + basis/wngb[0].size) - aNode[n].weight;
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					other[iter] = bNode[nbest].weight;
					self[iter] = Numerics.min(bNode[0].weight + basis/bNode[0].size,
												pNode[0].weight + basis/pNode[0].size);
					
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "+bNode[nbest].weight
																	+", s= "+bNode[0].weight
																	+", o= "+pNode[0].weight
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
																	
				   if (verbose) if (bcost<0 && first) {
				   	   // that's not the originally computed score, is it??
				   	   System.out.println(iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
				   	}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);
						/*
						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}
						*/
					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					assoc[lbest][0] = new Triple(id);
					assoc[lpair][0] = new Triple(id);
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
										
					// reduce the number of clusters
					nclusters--;
					
					// update the neighbors too, but only up to a point 
					// (-> use both strategies to gain speed in the two extreme worst case scenarios)
					// (both strategies = this and the track updating)
					for (int b=1;b<aNode.length && b<maxlength;b++) {
						if (active.get(aNode[b].id)) {
							widthsearchtime += assoc[aNode[b].id].length;
							for (int c=1;c<assoc[aNode[b].id].length && c<maxlength;c++) {
								if (latest[assoc[aNode[b].id][c].id]==lbest || latest[assoc[aNode[b].id][c].id]==lpair)
									latest[assoc[aNode[b].id][c].id] = id;
								//if (assoc[aNode[b].id][c].id==lbest || assoc[aNode[b].id][c].id==lpair)
								//	assoc[aNode[b].id][c].id = id;
							}
						}
					}
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringMahalanobis(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : mahalanobis distances");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		
		// store cluster location dynamically? (only for new clusters?)
		clusterPos = new int[2*nlb][];

		// look for joint boundaries
		BitSet ngbcluster = new BitSet(nx*ny*nz);
		boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		double var1=0, var2=0;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// mean, variance, size: compute from data??
					
					// sum(uv) = sum(u) + sum(v)
					// sq2(uv) = sq2(u) + sq2(v) + 1/n(uv) [n(u)/n(v)*m(v) - n(v)/n(u)*m(u)]^2
					// n(uv) = n(u) + n(v)
					aNode[0] = new Triple(id, bNode[0].weight + pNode[0].weight, 
												bNode[0].delta + pNode[0].delta 
												+ Numerics.square(bNode[0].weight/bNode[0].size-pNode[0].weight/pNode[0].size)
													*bNode[0].size*pNode[0].size/(bNode[0].size+pNode[0].size),  
												bNode[0].size + pNode[0].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += bNode[n].weight;
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += pNode[n].weight;
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lpair] = null;
					} else {
						clusterPos[id] = new int[2];
						clusterPos[id][0] = lbest;
						clusterPos[id][1] = lpair;
						// no need to check
						boundaryScore[invertlabeling[lbest]] = bcost;
						boundaryScore[invertlabeling[lpair]] = bcost;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// estimate the variance with prior
						var1 = (aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis-1.0);
						var2 = (wngb[0].delta + basis*imgscale*imgscale)/(wngb[0].size+basis-1.0);
						
						// build the t-test for each possible pairing
						float tval = tTest(aNode[0].weight/aNode[0].size, wngb[0].weight/wngb[0].size,
										var1, var2, Numerics.max(aNode[0].size,2), Numerics.max(wngb[0].size,2));
						
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/((connect+1.0f)*nlb)) );
						
						// Bonferroni correction ?
						//float threshold = pvalue/(aNode[0].size+wngb[0].size-1.0f);
						
						// Sidak correction ?
						float threshold = 1.0f - (float)FastMath.pow(1.0f-pvalue, 1.0f/(aNode[0].size+wngb[0].size-1.0f));
						
						aNode[n].delta = wsize*(tval-threshold);
						
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = (float)var1;
					other[iter] = (float)var2;
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "
																	+", m= "+aNode[0].weight/aNode[0].size
																	+", s= "+Math.sqrt(aNode[0].delta/aNode[0].size)
																	+", s0= "+Math.sqrt((aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis))
																	+", n= "+aNode[0].size
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
					// not a correct stopping criterion								
					if (verbose) if (bcost<0 && first) {
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+", m= "+aNode[0].weight/aNode[0].size
																	+", s= "+Math.sqrt(aNode[0].delta/aNode[0].size)
																	+", s0= "+Math.sqrt((aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis))
																	+", n= "+aNode[0].size
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
					}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);
						/*
						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}
						*/
					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					//assoc[lbest][0] = new Triple(id);
					//assoc[lpair][0] = new Triple(id);
					assoc[lbest][0] = bNode[0];
					assoc[lpair][0] = pNode[0];
					assoc[lbest][0].id = id;
					assoc[lpair][0].id = id;
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
					
					// reduce the number of clusters
					nclusters--;
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	public final void hierarchicalClusteringNormalDistributions(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : normal distribution distances");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		
		// store cluster location dynamically? (only for new clusters?)
		clusterPos = new int[2*nlb][];

		// look for joint boundaries
		BitSet ngbcluster = new BitSet(nx*ny*nz);
		boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		double var1=0, var2=0, var12=0;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// mean, variance, size: compute from data??
					
					// sum(uv) = sum(u) + sum(v)
					// sq2(uv) = sq2(u) + sq2(v) + 1/n(uv) [n(u)/n(v)*m(v) - n(v)/n(u)*m(u)]^2
					// n(uv) = n(u) + n(v)
					aNode[0] = new Triple(id, bNode[0].weight + pNode[0].weight, 
												bNode[0].delta + pNode[0].delta 
												+ Numerics.square(bNode[0].weight/bNode[0].size-pNode[0].weight/pNode[0].size)
													*bNode[0].size*pNode[0].size/(bNode[0].size+pNode[0].size),  
												bNode[0].size + pNode[0].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += bNode[n].weight;
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += pNode[n].weight;
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyz] = bcost;
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						boundaryScore[xyz] = bcost;
						for (int k = 0; k<6; k++) {
							int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
							if (ngbcluster.get(xyzn)) {
								boundaryScore[xyzn] = bcost;
							}
						}
						clusterPos[lpair] = null;
					} else {
						clusterPos[id] = new int[2];
						clusterPos[id][0] = lbest;
						clusterPos[id][1] = lpair;
						// no need to check
						boundaryScore[invertlabeling[lbest]] = bcost;
						boundaryScore[invertlabeling[lpair]] = bcost;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// estimate the variance with prior
						var1 = (aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis-1.0);
						var2 = (wngb[0].delta + basis*imgscale*imgscale)/(wngb[0].size+basis-1.0);
						
						if (mode==KULLBACKLEIBLER) {
							var12 = (wngb[0].delta + aNode[0].delta
									+ (wngb[0].size/aNode[0].size*aNode[0].weight*aNode[0].weight
									+ aNode[0].size/wngb[0].size*wngb[0].weight*wngb[0].weight
									- 2.0*aNode[0].weight*wngb[0].weight)/(aNode[0].size+wngb[0].size)
									+ basis*imgscale*imgscale)/(aNode[0].size+wngb[0].size+basis-1.0);
						}
						
						// build the divergence metric
						float metric = (float)distributionMetric(aNode[0].weight/aNode[0].size, wngb[0].weight/wngb[0].size, 
																	var1, var2, var12, Numerics.max(aNode[0].size,2), Numerics.max(wngb[0].size,2));
						
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/((connect+1.0f)*nlb)) );
						
						
						float threshold = distributionThreshold(aNode[0].size, wngb[0].size);
						// Bonferroni correction ?
						//float threshold = pvalue/(aNode[0].size+wngb[0].size-1.0f);
						
						// Sidak correction ?
						//float threshold = 1.0f - (float)FastMath.pow(1.0f-pvalue, 1.0f/(aNode[0].size+wngb[0].size-1.0f));
						
						aNode[n].delta = wsize*(metric-threshold);
						
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = (float)var1;
					other[iter] = (float)var2;
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "
																	+", m= "+aNode[0].weight/aNode[0].size
																	+", s= "+Math.sqrt(aNode[0].delta/aNode[0].size)
																	+", s0= "+Math.sqrt((aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis))
																	+", n= "+aNode[0].size
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
					// not a correct stopping criterion								
					if (verbose) if (bcost<0 && first) {
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+", m= "+aNode[0].weight/aNode[0].size
																	+", s= "+Math.sqrt(aNode[0].delta/aNode[0].size)
																	+", s0= "+Math.sqrt((aNode[0].delta + basis*imgscale*imgscale)/(aNode[0].size+basis))
																	+", n= "+aNode[0].size
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
					}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);
						/*
						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}
						*/
					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					//assoc[lbest][0] = new Triple(id);
					//assoc[lpair][0] = new Triple(id);
					assoc[lbest][0] = bNode[0];
					assoc[lpair][0] = pNode[0];
					assoc[lbest][0].id = id;
					assoc[lpair][0].id = id;
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
					
					// reduce the number of clusters
					nclusters--;
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	public final int firstAssociationLoss(int nt) {
		if (debug) System.out.println("find largest number of clusters with high association");

		int zerolb = -1;
		for (int t=1;t<nlb;t++) {
			//if (self[t]<other[t]) {
			if (cost[t]<0) {
				zerolb = t;
				break;
			}
		}
		if (zerolb==-1) zerolb = nlb-nt-1;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	public final int firstDifferenceLoss(int nt) {
		if (debug) System.out.println("find largest number of clusters with low difference");

		int zerolb = -1;
		for (int t=1;t<nlb;t++) {
			//if (self[t]<other[t]) {
			if (cost[t]>0) {
				zerolb = t;
				break;
			}
		}
		if (zerolb==-1) zerolb = nlb-nt-1;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	public final int firstAbovePvalue(int nt) {
		if (debug) System.out.println("find largest number of clusters below P-value");

		int zerolb = -1;
		for (int t=1;t<nlb;t++) {
			//if (self[t]<other[t]) {
			if (cost[t]<0) {
				zerolb = t;
				break;
			}
		}
		if (zerolb==-1) zerolb = nlb-nt-1;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	public final int lastAssociationGain(int nt) {
		if (debug) System.out.println("find smallest number of clusters with high association");

		int zerolb = -1;
		for (int t=nlb-1;t>0;t--) {
			//if (self[t]>other[t]) {
			if (cost[t]>0) {
				zerolb = t+1;
				break;
			}
		}
		if (zerolb==-1) zerolb = 0;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	private final void 	updateLatestLabels(int maxlb) {	
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringSequential3da(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				tmp[x][y][z] = latest[labeling[xyz]];
			}
		}
		return tmp;
	}
	
	/*
	// given a clustering result, update the labeling
	public final float[][][] relabelClusteringSequential3da(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}

		int[] relabel = new int[2*nlb];
		int counter = 2;
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			if (latest[n]==n) {
				relabel[n] = counter;
				counter++;
			} else {
				relabel[n] = relabel[latest[n]];
			}
		}
		System.out.println("number of clusters found: "+(counter-2));
		
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				tmp[x][y][z] = relabel[labeling[xyz]];
			}
		}
		return tmp;
	}
	*/
	/*
	public final void relabelClustering3da(float[][][] clusters) {
		float[] lblist = ObjectLabeling.listLabels(clusters, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				boolean found = false;
				for (int n=0;n<lblist.length && !found;n++) {
					if (lblist[n]==clusters[x][y][z]) {
						clusters[x][y][z] = n;
						found = true;
					}
				}
			}
		}
		return;
	}
	*/
	// given a clustering result, update the labeling
	public final float[] exportClusteringSequential(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = latest[labeling[xyz]];
				tmp[xyz] = lb;
			}
		}
		/*
		// second pass to get smaller label values?
		float[] lblist = ObjectLabeling.listLabels(tmp, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				boolean found = false;
				for (int n=0;n<lblist.length && !found;n++) {
					if (lblist[n]==tmp[xyz]) {
						tmp[xyz] = n;
						found = true;
					}
				}
			}
		}
		*/
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final String clusteringStatisticsSequential3da(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// create volume measurements
		int[] volume = new int[2*nlb];
		int totalvol = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				volume[latest[labeling[xyz]]]++;
				totalvol++;
			}
		}
		int maxvol = Numerics.max(volume);
		int nclusters = 0;
		for (int n=0;n<2*nlb;n++) if (volume[n]>0) nclusters++; 
		String stats = "cluster number : "+nclusters+", ratio (%) = "+(100.0f*nclusters/(float)totalvol)+
		", \n volume : "+totalvol+", max cluster size (%) = "+(100.0f*maxvol/(float)totalvol);

		return stats;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] filteredClusters3da(int maxlb, float minratio, float maxratio) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// create volume measurements
		int[] volume = new int[2*nlb];
		int totalvol = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				volume[latest[labeling[xyz]]]++;
				totalvol++;
			}
		}
		
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				if (volume[latest[labeling[xyz]]]>minratio*totalvol && volume[latest[labeling[xyz]]]<maxratio*totalvol) {
					tmp[x][y][z] = latest[labeling[xyz]]+1;
				} else {
					tmp[x][y][z] = 1;
				}
			}
		}
		
		// second pass to get smaller label values?
		float[] lblist = ObjectLabeling.listLabels(tmp, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				boolean found = false;
				for (int n=0;n<lblist.length && !found;n++) {
					if (lblist[n]==tmp[x][y][z]) {
						tmp[x][y][z] = n;
						found = true;
					}
				}
			}
		}

		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] filteredClusters(int maxlb, float minratio, float maxratio) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// create volume measurements
		int[] volume = new int[2*nlb];
		int totalvol = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				volume[latest[labeling[xyz]]]++;
				totalvol++;
			}
		}
		
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				if (volume[latest[labeling[xyz]]]>minratio*totalvol && volume[latest[labeling[xyz]]]<maxratio*totalvol) {
					tmp[xyz] = latest[labeling[xyz]]+1;
				} else {
					tmp[xyz] = 1;
				}
			}
		}
		
		// second pass to get smaller label values?
		float[] lblist = ObjectLabeling.listLabels(tmp, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				boolean found = false;
				for (int n=0;n<lblist.length && !found;n++) {
					if (lblist[n]==tmp[xyz]) {
						tmp[xyz] = n;
						found = true;
					}
				}
			}
		}

		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterDegree(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				tmp[xyz] = assoc[lb][0].delta;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterSelfWeight(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				tmp[xyz] = assoc[lb][0].weight;
			}
		}
		return tmp;
	}
	// given a clustering result, update the labeling
	public final float[] exportClusterSelfWeightSequential(int maxlb) {
		
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = latest[labeling[xyz]];
				tmp[xyz] = assoc[lb][0].weight/assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusterSelfWeight3da(int maxlb) {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				tmp[x][y][z] = assoc[lb][0].weight/assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusterSelfWeight3da(float[][][] clusters) {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = (int)clusters[x][y][z];
				tmp[x][y][z] = assoc[lb][0].weight/assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterSelfWeight(float[] clusters) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = (int)clusters[xyz];
				tmp[xyz] = assoc[lb][0].weight/assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterSelfDelta(float[] clusters) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = (int)clusters[xyz];
				tmp[xyz] = assoc[lb][0].delta/assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterSize(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				tmp[xyz] = assoc[lb][0].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportMinClusterWeight(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				float w = 1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					w = Numerics.min(w, assoc[lb][l].weight/assoc[lb][l].size);
				}
				if (w==1e9f) w = 0.0f;
				tmp[xyz] = w;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportMinClusterWeight3da(int maxlb) {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				
				//lb = latest[labeling[xyz]];
				float w = 1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					w = Numerics.min(w, assoc[lb][l].weight/assoc[lb][l].size);
				}
				if (w==1e9f) w = 0.0f;
				tmp[x][y][z] = w;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportOrigMinClusterWeight3da() {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				/*
				while (assoc.get(lb).get(0).id>0 && assoc.get(lb).get(0).id!=lb 
						&& !active.get(lb) && assoc.get(lb).get(0).id<2*nlb-maxlb) {
					lb = assoc.get(lb).get(0).id;
				}
				*/
				float w = 1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					w = Numerics.min(w, assoc[lb][l].weight/assoc[lb][l].size);
				}
				if (w==1e9f) w = 0.0f;
				tmp[x][y][z] = w;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][][] exportOrigAllClusterWeight3da() {
		float[][][][] tmp = new float[nx][ny][nz][connect];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				for (int l=1;l<assoc[lb].length;l++) {
					tmp[x][y][z][l-1] = assoc[lb][l].weight/assoc[lb][l].size;
				}
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][] exportOrigAllClusterWeight() {
		float[][] tmp = new float[connect][nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				for (int l=1;l<assoc[lb].length;l++) {
					tmp[l-1][xyz] = assoc[lb][l].weight/assoc[lb][l].size;
				}
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportOrigMaxClusterDelta() {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				float w = -1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					w = Numerics.max(w, assoc[lb][l].delta);
				}
				if (w==-1e9f) w = 0.0f;
				tmp[xyz] = w;
			}
		}
		return tmp;
	}
	
	public final float[] exportMaxClusterWeight(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				float w = -1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					w = Numerics.max(w, assoc[lb][l].weight/assoc[lb][l].size);
				}
				tmp[xyz] = w;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportMinClusterDelta(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				float d = 1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					d = Numerics.min(d, assoc[lb][l].delta);
				}
				tmp[xyz] = d;
			}
		}
		return tmp;
	}
	
	public final float[] exportMaxClusterDelta(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				float d = -1e9f;
				for (int l=1;l<assoc[lb].length;l++) {
					d = Numerics.max(d, assoc[lb][l].delta);
				}
				tmp[xyz] = d;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportLabeling() {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			tmp[xyz] = labeling[xyz];
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[] exportClusterLocalization(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int lb=1;lb<=2*nlb-maxlb;lb++) {
			if (clusterPos[lb]!=null) {
				for (int n = 0;n<clusterPos[lb].length; n++) {
					int xyz = invertlabeling[clusterPos[lb][n]];
					tmp[xyz] = lb;
				}
			}
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[] getClusterBoundaryScore() {
		return boundaryScore;
	}
		
	// computes the association parameter for the image
	private final float associationWeight(int id1, int id2) {
		if (mode==PROFILES_DIST) return profileDistanceWeight(id1, id2);
		if (mode==PROFILES_GDIST) return profileGaussianDistanceWeight(id1, id2);
		if (mode==PROFILES_CORR) return correlationWeight(id1, id2);
		if (mode==PROFILES_ANGL) return correlationAngleWeight(id1, id2);
		if (mode==PROFILES_CDIST) return correlationDistanceWeight(id1, id2);
		if (mode==PROFILES_CGDIST) return correlationGaussianDistanceWeight(id1, id2);
		if (mode==NONPARAM) return nonParamDistanceWeight(id1, id2);
		if (mode==GDIST) return scalarGaussianDistanceWeight(id1, id2);
		if (mode==HNDIST) return scalarGaussianDistanceWeight(id1, id2);
		if (mode==AUTOGDIST) return scalarGaussianDistanceWeight(id1, id2);
		if (mode==ROBUSTGDIST) return scalarGaussianDistanceWeight(id1, id2);
		if (mode==SCALAR) return scalarDistanceWeight(id1, id2);
		if (mode==BINARY) return binaryWeight(id1, id2);
		if (mode==GBFH) return scalarDistance(id1, id2);
		if (mode==MAHALANOBIS) return ttestWeight(id1, id2);
		if (mode==JENSENSHANNON) return jsdivWeight(id1, id2);
		if (mode==HELLINGER) return hellingerWeight(id1, id2);
		if (mode==KULLBACKLEIBLER) return kldivWeight(id1, id2);
		return 0.5f;
	}
	private final float scalarDistanceWeight(int id1, int id2) {
		return 1.0f/( basis + Numerics.square( (image[0][id1]-image[0][id2])/imgscale ) );
	}
	private final float nonParamDistanceWeight(int id1, int id2) {
		return distribution.getHistogramCount(Numerics.abs(image[0][id1]-image[0][id2]));
	}
	private final float scalarGaussianDistanceWeight(int id1, int id2) {
		return (float)Math.exp( -0.5*Numerics.square( (image[0][id1]-image[0][id2])/imgscale ) );
	}
	private final float ttestWeight(int id1, int id2) {
		return (float)(tTest(image[0][id1],image[0][id2],imgscale*imgscale,imgscale*imgscale,2,2)-pvalue);
	}
	private final float jsdivWeight(int id1, int id2) {
		return (float)(jsDiv(image[0][id1],image[0][id2],imgscale*imgscale,imgscale*imgscale)-pvalue);
	}
	private final float hellingerWeight(int id1, int id2) {
		return (float)(hellingerCoeff(image[0][id1],image[0][id2],imgscale*imgscale,imgscale*imgscale)-pvalue);
	}
	private final float kldivWeight(int id1, int id2) {
		return (float)(klDiv(image[0][id1],image[0][id2],
								imgscale*imgscale,imgscale*imgscale,imgscale*imgscale,2,2)-pvalue);
	}
	private final float binaryWeight(int id1, int id2) {
		if (image[0][id1]==image[0][id2]) return basis;
		else return imgscale;
	}
	private final float profileDistanceWeight(int id1, int id2) {
		float dist=0.0f;
		for (int t=0;t<image.length;t++) {
			dist += Numerics.square(image[t][id1]-image[t][id2]);
		}
		return 1.0f/( basis + dist/(imgscale*imgscale) );
	}
	private final float profileGaussianDistanceWeight(int id1, int id2) {
		float dist=0.0f;
		for (int t=0;t<image.length;t++) {
			dist += Numerics.square(image[t][id1]-image[t][id2]);
		}
		return (float)Math.exp( -0.5*dist/(imgscale*imgscale) );
	}
	private final float correlationWeight(int id1, int id2) {
		float corr = 0.0f;
		float m1 = 0.0f, m2 = 0.0f;
		float v1 = 0.0f, v2 = 0.0f;
		for (int t=0;t<image.length;t++) {
			m1 += image[t][id1];
			m2 += image[t][id2];
		}
		m1 /= image.length;
		m2 /= image.length;
		for (int t=0;t<image.length;t++) {
			v1 += (image[t][id1]-m1)*(image[t][id1]-m1);
			v2 += (image[t][id2]-m2)*(image[t][id2]-m2);
		}
		v1 = (float)Math.sqrt(v1/(image.length-1.0f));
		v2 = (float)Math.sqrt(v2/(image.length-1.0f));
		for (int t=0;t<image.length;t++) {
			corr += (image[t][id1]-m1)*(image[t][id2]-m2)/(v1*v2);
		}
		return Numerics.bounded( (corr/image.length-basis)/(1.0f-basis), 0.0f, 1.0f);
	}
	
	private final float correlationAngleWeight(int id1, int id2) {
		float corr = 0.0f;
		float m1 = 0.0f, m2 = 0.0f;
		float v1 = 0.0f, v2 = 0.0f;
		for (int t=0;t<image.length;t++) {
			m1 += image[t][id1];
			m2 += image[t][id2];
		}
		m1 /= image.length;
		m2 /= image.length;
		for (int t=0;t<image.length;t++) {
			v1 += (image[t][id1]-m1)*(image[t][id1]-m1);
			v2 += (image[t][id2]-m2)*(image[t][id2]-m2);
		}
		v1 = (float)Math.sqrt(v1/(image.length-1.0f));
		v2 = (float)Math.sqrt(v2/(image.length-1.0f));
		for (int t=0;t<image.length;t++) {
			corr += (image[t][id1]-m1)*(image[t][id2]-m2)/(v1*v2);
		}
		// score = (sin^-1(corr) - theta_0)/(90-theta_0) bounded in [0,1]
		return (float)Numerics.max( (Math.asin(Numerics.bounded(corr/image.length,0.0,1.0)) - PI2*basis)/(PI2 - PI2*basis), 0.0);
	}
	
	private final float correlationDistanceWeight(int id1, int id2) {
		float corr = 0.0f;
		float m1 = 0.0f, m2 = 0.0f;
		float v1 = 0.0f, v2 = 0.0f;
		for (int t=0;t<image.length;t++) {
			m1 += image[t][id1];
			m2 += image[t][id2];
		}
		m1 /= image.length;
		m2 /= image.length;
		for (int t=0;t<image.length;t++) {
			v1 += (image[t][id1]-m1)*(image[t][id1]-m1);
			v2 += (image[t][id2]-m2)*(image[t][id2]-m2);
		}
		v1 = (float)Math.sqrt(v1/(image.length-1.0f));
		v2 = (float)Math.sqrt(v2/(image.length-1.0f));
		for (int t=0;t<image.length;t++) {
			corr += (image[t][id1]-m1)*(image[t][id2]-m2)/(v1*v2);
		}
		// score = (sin^-1(corr) - theta_0)/(90-theta_0) bounded in [0,1]
		return (float)(1.0/(1.0 + Numerics.square(Math.acos(Numerics.bounded(corr/image.length,0.0,1.0))/imgscale) ) );
	}
	
	private final float correlationGaussianDistanceWeight(int id1, int id2) {
		float corr = 0.0f;
		float m1 = 0.0f, m2 = 0.0f;
		float v1 = 0.0f, v2 = 0.0f;
		for (int t=0;t<image.length;t++) {
			m1 += image[t][id1];
			m2 += image[t][id2];
		}
		m1 /= image.length;
		m2 /= image.length;
		for (int t=0;t<image.length;t++) {
			v1 += (image[t][id1]-m1)*(image[t][id1]-m1);
			v2 += (image[t][id2]-m2)*(image[t][id2]-m2);
		}
		v1 = (float)Math.sqrt(v1/(image.length-1.0f));
		v2 = (float)Math.sqrt(v2/(image.length-1.0f));
		for (int t=0;t<image.length;t++) {
			corr += (image[t][id1]-m1)*(image[t][id2]-m2)/(v1*v2);
		}
		// score = (sin^-1(corr) - theta_0)/(90-theta_0) bounded in [0,1]
		return (float)Math.exp( -0.5*Numerics.square(Math.acos(Numerics.bounded(corr/image.length,0.0,1.0))/imgscale) );
	}

	private final float dataDistance(int id1, int id2) {
		if (mode==BINARY) return binaryDistance(id1, id2);
		if (mode==SCALAR) return scalarDistance(id1, id2);
		if (mode==GDIST) return scalarDistance(id1, id2);
		if (mode==HNDIST) return scalarDistance(id1, id2);
		if (mode==AUTOGDIST) return scalarDistance(id1, id2);
		if (mode==ROBUSTGDIST) return scalarDistance(id1, id2);
		if (mode==NONPARAM) return scalarDistance(id1, id2);
		if (mode==PROFILES_DIST) return profileEuclideanDistance(id1, id2);
		if (mode==PROFILES_GDIST) return profileEuclideanDistance(id1, id2);
		if (mode==PROFILES_CORR) return profileCorrelation(id1, id2);
		if (mode==PROFILES_ANGL) return profileProduct(id1, id2);
		if (mode==PROFILES_CDIST) return profileCorrelation(id1, id2);
		if (mode==PROFILES_CGDIST) return profileCorrelation(id1, id2);
		if (mode==GBFH) return scalarDistance(id1, id2);
		if (mode==MAHALANOBIS) return scalarDistance(id1, id2);
		if (mode==JENSENSHANNON) return scalarDistance(id1, id2);
		if (mode==HELLINGER) return scalarDistance(id1, id2);
		if (mode==KULLBACKLEIBLER) return scalarDistance(id1, id2);
		return 0.5f;
	}
	
	private final float binaryDistance(int id1, int id2) {
		if (image[0][id1]==image[0][id2]) return 0.0f;
		else return 1.0f;
	}
	
	private final float scalarDistance(int id1, int id2) {
		return Numerics.abs(image[0][id1]-image[0][id2]);
	}
	
	private final float profileEuclideanDistance(int id1, int id2) {
		float dist=0.0f;
		for (int t=0;t<image.length;t++) {
			dist += Numerics.square(image[t][id1]-image[t][id2]);
		}
		return (float)Math.sqrt(dist);
	}
	
	private final float profileProduct(int id1, int id2) {
		float dist=0.0f;
		for (int t=0;t<image.length;t++) {
			dist += image[t][id1]*image[t][id2];
		}
		return dist;
	}
	
	private final float profileCorrelation(int id1, int id2) {
		float corr = 0.0f;
		float m1 = 0.0f, m2 = 0.0f;
		float v1 = 0.0f, v2 = 0.0f;
		for (int t=0;t<image.length;t++) {
			m1 += image[t][id1];
			m2 += image[t][id2];
		}
		m1 /= image.length;
		m2 /= image.length;
		for (int t=0;t<image.length;t++) {
			v1 += (image[t][id1]-m1)*(image[t][id1]-m1);
			v2 += (image[t][id2]-m2)*(image[t][id2]-m2);
		}
		v1 = (float)Math.sqrt(v1/(image.length-1.0f));
		v2 = (float)Math.sqrt(v2/(image.length-1.0f));
		for (int t=0;t<image.length;t++) {
			corr += (image[t][id1]-m1)*(image[t][id2]-m2)/(v1*v2);
		}
		return corr;
	}
	
	private final float distributionThreshold(double n1, double n2) {
		// Sidak correction
		if (mode==MAHALANOBIS) return 1.0f - (float)FastMath.pow(1.0f-pvalue, 1.0f/(n1+n2-1.0f));
		// JS div threshold??				
		else if (mode==JENSENSHANNON) return pvalue;
		// Hellinger distance threshold??
		else if (mode==HELLINGER) return pvalue;
		// JS div threshold??				
		else if (mode==KULLBACKLEIBLER) return pvalue;
		return 0.5f;
	}

	private final double distributionMetric(final double m1, final double m2,
											  final double v1, final double v2, final double v12,
											  final double n1, final double n2) {
		// T-test
		if (mode==MAHALANOBIS) return tTest(m1, m2, v1, v2, n1, n2);
		// JS div 		
		else if (mode==JENSENSHANNON) return jsDiv(m1, m2, v1, v2);
		// Hellinger distance 
		else if (mode==HELLINGER) return hellingerCoeff(m1, m2, v1, v2);
		// Kullback-Leibler div
		else if (mode==KULLBACKLEIBLER) return klDiv(m1, m2, v1, v2, v12, n1, n2);
		return 0.5;
	}

	/**
     * Computes p-value for 2-sided, 2-sample t-test.
     * <p>
     * Does not assume subpopulation variances are equal. Degrees of freedom
     * are estimated from the data. Adapted from commons-math at apache.org</p>
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param v1 first sample variance
     * @param v2 second sample variance
     * @param n1 first sample n
     * @param n2 second sample n
     * @return p-value
     * @throws MaxCountExceededException if an error occurs computing the p-value
     */
    private float tTest(final double m1, final double m2,
                          final double v1, final double v2,
                          final double n1, final double n2) {
       
    	final double t = FastMath.abs((m1 - m2) / FastMath.sqrt((v1 / n1) + (v2 / n2)));
    	
    	
    	final double degreesOfFreedom;
    	
    	if (n1==1 || n2==1) degreesOfFreedom = 0;
    	else degreesOfFreedom = (((v1 / n1) + (v2 / n2)) * ((v1 / n1) + (v2 / n2))) /
       								((v1 * v1) / (n1 * n1 * (n1 - 1d)) + (v2 * v2) /
       									(n2 * n2 * (n2 - 1d)));
        
       	TDistribution distribution = new TDistribution(degreesOfFreedom);
        return (float)(2.0 * distribution.cumulativeProbability(-t));
    }

    private double jsDiv(final double m1, final double m2,
                          		final double v1, final double v2) {
       
    	return 1.0 - FastMath.sqrt( (m1-m2)*(m1-m2)/(v1+v2) - FastMath.log(2.0*FastMath.sqrt(v1*v2)/(v1+v2)) );
    }

    private double hellingerCoeff(final double m1, final double m2,
                          			final double v1, final double v2) {
       
    	return 1.0 - FastMath.sqrt( 1.0 - FastMath.sqrt( 2.0*FastMath.sqrt(v1*v2)/(v1+v2) )*FastMath.exp(-0.25*(m1-m2)*(m1-m2)/(v1+v2) ) );
    }

    private double klDiv(final double m1, final double m2,
							final double v1, final double v2, final double v12,
							  final double n1, final double n2) {
    	
		final double m12 = (n1*m1+n2*m2)/(n1+n2);
		
		return 2.0 - 0.5*( (m1-m12)*(m1-m12)/v12 + (m2-m12)*(m2-m12)/v12  
    						+ (v1+v2)/v12 - FastMath.log((v1/v12)*(v2/v12)) );
    }

}
