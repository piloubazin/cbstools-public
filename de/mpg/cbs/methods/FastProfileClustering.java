package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

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
 
public class FastProfileClustering {
		
	// data buffers
	private 	float[][]	image;  			// original data
	private 	int			nx,ny,nz,nc;   		// image dimensions
	private 	float		rx,ry,rz;   		// image resolutions
	private		boolean[]	mask;
	private		float		imgscale;
	private		float		basis;
	private		int			connect;
	//private		String[]	modes = {"profiles", "binary", "scalar", "noise"};
	
	private		byte		mode;
	private	static final	byte	BINARY = 2;
	private	static final	byte	SCALAR = 3;
	private	static final	byte	NOISE = 4;
	private	static final	byte	PROFILES_DIST = 10;
	private	static final	byte	PROFILES_GDIST = 11;
	private	static final	byte	PROFILES_CORR = 12;
	private	static final	byte	PROFILES_ANGL = 13;
	private	static final	byte	PROFILES_CDIST = 14;
	private	static final	byte	PROFILES_CGDIST = 15;
	
	private static final double PI2 = Math.PI/2.0;
	
	// algorithm quantities
	private		int[]		labeling;
	private		int[]		clustering;
	private		int			nlb;
	private		int			key;
	
	//private		ArrayList<Float>	degree;
	//private		ArrayList<ArrayList<Triple>>	assoc;
	private		Triple[][]						assoc;
	private		BitSet							active;
	private		BinaryHeapPair					maxtree;
	private		float[]							cost;
	private		float[]							self;
	private		float[]							other;
	private		int[]							latest;
	
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
			this.weight = 1.0f;
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
	public FastProfileClustering(float[][] img_, boolean[] msk_,
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
		
		if (mod_.equals("profiles-distance")) 						mode = PROFILES_DIST;
		else if (mod_.equals("profiles-gauss-distance")) 			mode = PROFILES_GDIST;
		else if (mod_.equals("profiles-correlation")) 				mode = PROFILES_CORR;
		else if (mod_.equals("profiles-angle")) 					mode = PROFILES_ANGL;
		else if (mod_.equals("profiles-correlation-dist")) 		mode = PROFILES_CDIST;
		else if (mod_.equals("profiles-correlation-gauss-dist")) 	mode = PROFILES_CGDIST;
		else if (mod_.equals("binary")) 	mode = BINARY;
		else if (mod_.equals("scalar")) 	mode = SCALAR;
		else if (mod_.equals("noise")) 	mode = NOISE;

		// init: find the labels for used data points
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
		if (debug) System.out.println("initial labeling ("+nlb+" voxels)");
		
	}

	final public void finalize() {
		image = null;
		labeling = null;
		assoc = null;
		System.gc();
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeights() {
		
		if (debug) System.out.println("-- weight initialization --");
		
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
				// build the weight array: only with labels of inferior value
				// and the self-weight as well
				/* duplicate all the links... ? needed for finding the new neighbors
				ArrayList<Triple> node = new ArrayList<Triple>(nl);
				node.add(0, new Triple(labeling[xyz], 1.0f, 0.0f));
				int l=1;
				for (int n=0;n<nb;n++) if (trp[n].id<labeling[xyz]) {
					node.add(l, trp[n]);
					l++;
				}
				*/
				/*
				// link in 6 directions
				ArrayList<Triple> node = new ArrayList<Triple>(nb);
				node.add(0, new Triple(labeling[xyz], 1.0f));
				for (int n=0;n<nb;n++) {
					node.add(n+1, trp[n]);
				}
				// build the degree array
				float deg = node.get(0).weight/node.get(0).size;
				for (int n=0;n<nb;n++) {
					deg += trp[n].weight/trp[n].size/nb;
				}
				//degree.set(labeling[xyz], new Float(deg));
				node.get(0).delta = deg;

				// store the values
				assoc.add(labeling[xyz], node);
				*/
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
	public final void computeWeightMetrics() {
		
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
				if (z>0 && mask[xyz-nx*ny])	{
					distance[nb] = dataDistance(xyz,xyz-nx*ny);
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
		System.out.println("50% : "+hist.percentage(0.5f));
		System.out.println("25% : "+hist.percentage(0.25f));
		System.out.println("75% : "+hist.percentage(0.75f));
	}
	
	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClustering(int k0, boolean firststop) {
		
		if (debug) System.out.println("hierarchical clustering");
		
		// 1. Build a binary tree for the delta values
		maxtree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
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
					/*
					// recompute the delta value every time
					float score = node.get(b).weight/node.get(b).delta/(node.get(0).delta+assoc.get(node.get(b).id).get(0).delta);
					if (score>bestscore) {
						//if (debug) System.out.print(">");
						best = b;
						bestscore = score;
					}
					*/
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				maxtree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		while (maxtree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			int  	lbest = maxtree.getFirstId1();
			int 	nbest = maxtree.getFirstId2();
			maxtree.removeFirst();

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
					int l=1;
					BitSet counted = new BitSet();
					counted.set(0);
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = assoc[lbn][0].id;
									/*
									int lbn2 = assoc[lbn][0].id;
									if (lbn2!=lbn) {
										// update the id if not the same? no, it destroys the hierarchy :(
										//assoc[lbn][0].id = assoc[lbn2][0].id;
										latest[lbn] = latest[lbn2];
										lbn = lbn2;
									}
									*/
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
								latest[tracks.val[lb]] = lbn;
								//latest[bNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							// retrieve the up-to-date weight linking depending on creation date
							float wnb = bNode[n].weight;
							float snb = bNode[n].size;
							// maybe not: just combine all the original weights onto a single label
							/*
							if (lbn>lbest) {
								// neighbor is more recent; must update the weights
								boolean found=false;
								for (int m=1;m<assoc.get(lbn).size() && !found;m++) {
									if (assoc.get(lbn).get(m).id==lbest) {
										wnb = assoc.get(lbn).get(m).weight;
										snb = assoc.get(lbn).get(m).size;
										found=true;
									}
								}
								if (!found) System.out.print("!");
							}
							*/
							// check on all already created values : shouldn't be needed (?)
							boolean found=false;
							/* not needed?
							for (int m=1;m<aNode.size() && !found;m++) {
								if (lbn==aNode.get(m).id && lbn!=lpair) {
									// found: add the values
									found=true;
									aNode.get(m).weight += wnb;
									aNode.get(m).size += snb;
									//if (!counted.get(m)) counted.set(m);
									//else System.out.print(":");
								}
							}
							*/
							if (!found && lbn!=lpair) {
								aNode[l] = new Triple(lbn, wnb, snb);
								l++;
							}	
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
									lbn = assoc[lbn][0].id;
									/*
									int lbn2 = assoc[lbn][0].id;
									if (lbn2!=lbn) {
										// update the id if not the same? no, it destroys the hierarchy :(
										//assoc[lbn][0].id = assoc[lbn2][0].id;
										latest[lbn] = latest[lbn2];
										lbn = lbn2;
									}
									*/
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
								latest[tracks.val[lb]] = lbn;
								//latest[pNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							// retrieve the up-to-date weight linking to lpair ?
							// or is this counting twice??
							float wnb = pNode[n].weight;
							float snb = pNode[n].size;
							/*
							if (lbn>lpair) {
								boolean found=false;
								for (int m=1;m<assoc.get(lbn).size() && !found;m++) {
									if (assoc.get(lbn).get(m).id==lpair) {
										wnb = assoc.get(lbn).get(m).weight;
										snb = assoc.get(lbn).get(m).size;
										found = true;
									}
								}
								if (!found) System.out.print("?");
							}
							*/
							// check on all already created values
							boolean found=false;
							/* not needed if using the max, not the average */
							// or just the ones from the first group
							//for (int m=1;m<aNode.size() && !found;m++) {
							int bLength = l;
							for (int m=1;m<bLength && !found;m++) {
								if (lbn==aNode[m].id && lbn!=lbest) {
									// found: add the values
									found=true;
									aNode[m].weight += wnb;
									aNode[m].size += snb;
									//if (!counted.get(m)) counted.set(m);
									//else System.out.print(";");
								}
							}
							if (!found && lbn!=lbest) {
								// create a new one
								aNode[l] = new Triple(lbn, wnb, snb);
								l++;
							}
						}
					}
					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						aNode[n].delta = aNode[n].weight/aNode[n].size
												*( aNode[0].weight + wngb[0].weight + 2.0f*aNode[n].weight )
												/( aNode[0].size + wngb[0].size + 2.0f*aNode[n].size );
												
						// use the same criterion than for stopping? makes sense, but slows down the process
						/*
						aNode.get(n).delta = aNode.get(n).weight/aNode.get(n).size
												*( aNode.get(0).weight + wngb.get(0).weight + 2.0f*aNode.get(n).weight )
												/( aNode.get(0).size + wngb.get(0).size + 2.0f*aNode.get(n).size )
											 -(1.0f-aNode.get(n).weight/aNode.get(n).size)
											 	*aNode.get(0).weight/aNode.get(0).size
											 	*wngb.get(0).weight/wngb.get(0).size;
						*/
						/*					
						aNode.get(n).delta = ( aNode.get(0).weight + wngb.get(0).weight + 2.0f*aNode.get(n).weight )
											/( aNode.get(0).delta + wngb.get(0).delta )
											- aNode.get(0).weight/aNode.get(0).delta - wngb.get(0).weight/wngb.get(0).delta;
						*/
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					aNode[0].delta = aNode[0].weight/aNode[0].size;
					for (int n=1; n<aNode.length; n++) {
						aNode[0].delta += aNode[n].weight/aNode[n].size/(aNode.length-1.0f);
					}
					/*
					cost[iter] = 	- 0.5f*bNode.get(0).weight/bNode.get(0).size/bNode.get(0).delta
									- 0.5f*pNode.get(0).weight/pNode.get(0).size/pNode.get(0).delta
									+ aNode.get(0).weight/aNode.get(0).size/aNode.get(0).delta;
					*/
					self[iter] = 	bNode[nbest].weight/bNode[nbest].size * aNode[0].weight/aNode[0].size;
					other[iter] = (1.0f-bNode[nbest].weight/bNode[nbest].size) * bNode[0].weight/bNode[0].size
																							* pNode[0].weight/pNode[0].size;
					cost[iter] = self[iter]-other[iter];
					/*
					// update the association cost
					cost[iter] = (nclusters*cost[iter-1] - bNode.get(0).weight/bNode.get(0).size/bNode.get(0).delta
											   			  - pNode.get(0).weight/pNode.get(0).size/pNode.get(0).delta
											   			  + aNode.get(0).weight/aNode.get(0).size/aNode.get(0).delta)/(nclusters-1);
					*/
					/*
					cost[iter] = cost[iter-1] - bNode.get(0).weight/bNode.get(0).size/bNode.get(0).delta
											   - pNode.get(0).weight/pNode.get(0).size/pNode.get(0).delta
											   + aNode.get(0).weight/aNode.get(0).size/aNode.get(0).delta;
		
				   if (debug) if (cost[iter]/nclusters>1) {
					   System.out.println(nclusters+": c= "
															+cost[iter]/nclusters
															+" < "+aNode.get(0).weight
															+", "+aNode.get(0).size
															+", "+aNode.get(0).delta+"> "
															+" ("+bNode.get(0).delta+", "
															+bNode.get(0).weight+" | "
															+pNode.get(0).delta+", "
															+pNode.get(0).weight+" | "
															+bNode.get(nbest).delta+", "
															+bNode.get(nbest).size+", "
															+bNode.get(nbest).weight+")");					   
				
						// brute force: recompute the association cost globally (debug only)
						cost[iter] = 0.0f;
						for (int n=1;n<assoc.size();n++) if (active.get(n)) {
							cost[iter] += assoc.get(n).get(0).weight/assoc.get(n).get(0).size/assoc.get(n).get(0).delta;
						}
						System.out.println("new cost: "+cost[iter]/nclusters);
				   }
				   */
				   if (verbose) if (iter%(nlb/100)==0) System.out.println(iter+" / "+nclusters+": c= "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");		
																	
				   if (verbose) if (self[iter]-other[iter]<0 && first) {
				   	   System.out.println(iter+" / "+nclusters+": c= "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
				   }								
					/*						   
					if (debug) if (bNode.get(nbest).weight/bNode.get(nbest).size<0.5) System.out.println(nclusters+": c= "+cost[iter]
																						+", s= "+self[iter]
																						+", o= "+other[iter]
																						+" ("+bNode.get(0).delta+", "
																						+bNode.get(0).weight+" | "
																						+pNode.get(0).delta+", "
																						+pNode.get(0).weight+" | "
																						+bNode.get(nbest).delta+", "
																						+bNode.get(nbest).size+", "
																						+bNode.get(nbest).weight+")");					   
											   
					if (debug) if (Float.isNaN(aNode.get(0).delta)) {
						System.out.println(nclusters+": c= "+cost[iter]+" ("+bNode.get(0).delta+", "
																			+bNode.get(0).weight+", "
																			+bNode.get(0).size+" | "
																			+pNode.get(0).delta+", "
																			+pNode.get(0).weight+", "
																			+pNode.get(0).size+" | "
																			+aNode.get(0).delta+", "
																			+aNode.get(0).weight+", "
																			+aNode.get(0).size+")");					   
						
						for (int n=1; n<aNode.size(); n++) {
							System.out.print("<"+aNode.get(n).delta+", "+aNode.get(n).weight+", "+aNode.get(n).size+">");
						}
						System.out.print("\n");
					}
					*/
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					/*
					// update the latest label list
					for (int n=1;n<id;n++) {
						if (latest[n]==lbest || latest[n]==lpair) latest[n] = id;
					}
					*/
					
					if (aNode.length>1) {
						int best=1;
						/*
						float bestscore = 2.0f*aNode.get(1).weight/aNode.get(1).delta
											/( aNode.get(0).delta + assoc.get(aNode.get(1).id).get(0).delta );
											
						for (int b=2;b<aNode.size();b++) {
							float score = 2.0f*aNode.get(b).weight/aNode.get(b).delta
											/( aNode.get(0).delta + assoc.get(aNode.get(b).id).get(0).delta );
											
							if (score>bestscore) {
								best = b;
								bestscore = score;
							}
						}
						*/
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						maxtree.addValue(aNode[best].delta, id, best);

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
					/*
					ArrayList<Triple> tag = new ArrayList<Triple>(1);
					tag.add(0, new Triple(id));
					assoc.set(lbest, tag);
					assoc.set(lpair, tag);
					*/
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
					
					// update the neighbors too, but only up to a point 
					// (-> use both strategies to gain speed in the two extreme worst case scenarios)
					int maxlength = 100;
					for (int b=1;b<aNode.length && b<maxlength;b++) {
						for (int c=1;c<assoc[aNode[b].id].length && c<maxlength;c++) {
							if (active.get(aNode[b].id))
								if (assoc[aNode[b].id][c].id==lbest || assoc[aNode[b].id][c].id==lpair)
									assoc[aNode[b].id][c].id = id;
						}
					}	
					
					// reduce the number of clusters
					nclusters--;
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	// select the best cost by maximizing "curvature"
	public final int bestAssociationCurvature() {
		if (debug) System.out.println("find best number of clusters");

		float maxcurv = 0.0f;
		int maxlb = -1;
		for (int t=1;t<nlb;t++) if (cost[t+1]!=0) {
			float curv = cost[t] - cost[t+1];
			if (curv>maxcurv) {
				maxcurv = curv;
				maxlb = t;
			}
		}
		if (debug) System.out.println("-> "+(nlb-maxlb));
		return (nlb-maxlb);
	}
	
	// select the best cost by maximizing "curvature"
	public final int bestAssociationScore() {
		if (debug) System.out.println("find best number of clusters");

		float maxcost = 0.0f;
		int maxlb = -1;
		for (int t=1;t<nlb;t++) {
			if (cost[t]/(nlb-t)>maxcost) {
				maxcost = cost[t]/(nlb-t);
				maxlb = t;
			}
		}
		if (debug) System.out.println("-> "+(nlb-maxlb));
		return (nlb-maxlb);
	}
	
	public final int firstAssociationLoss() {
		if (debug) System.out.println("find largest number of clusters with high association");

		int zerolb = -1;
		for (int t=1;t<nlb;t++) {
			if (self[t]<other[t]) {
				zerolb = t;
				break;
			}
		}
		if (zerolb==-1) zerolb = nlb-1;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	public final int lastAssociationGain() {
		if (debug) System.out.println("find smallest number of clusters with high association");

		int zerolb = -1;
		for (int t=nlb-1;t>0;t--) {
			if (self[t]>other[t]) {
				zerolb = t+1;
				break;
			}
		}
		if (zerolb==-1) zerolb = 0;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	/*
	// given a clustering result, update the labeling
	public final int[] generateClustering(int maxlb) {
		clustering = new int[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id!=lb) {
					lb = assoc[lb][0].id;
				}
				clustering[xyz] = lb;
			}
		}
		return clustering;
	}
	*/
	/*
	// given a clustering result, update the labeling
	public final float[] exportClustering(int maxlb) {
		float[] tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& !active.get(lb) && assoc[lb][0].id<2*nlb-maxlb) {
					lb = assoc[lb][0].id;
				}
				tmp[xyz] = lb;
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
	*/
	/*
	// given a clustering result, update the labeling
	public final float[][][] exportClustering3da(int maxlb) {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				int lb = labeling[xyz];
				
				while (assoc[lb][0].id>0 && assoc[lb][0].id!=lb 
						&& assoc[lb][0].id<2*nlb-maxlb) {
					int lb2 = assoc[lb][0].id;
					assoc[lb][0].id = assoc[lb2][0].id;
					lb = lb2;
				}
				
				//lb = latest[labeling[xyz]];
				tmp[x][y][z] = lb;
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
					if (lblist[n]==tmp[x][y][z]) {
						tmp[x][y][z] = n;
						found = true;
					}
				}
			}
		}
		*/
		/*
		return tmp;
	}
	*/
	
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
	
	// computes the association parameter for the image
	private final float associationWeight(int id1, int id2) {
		if (mode==PROFILES_DIST) return profileDistanceWeight(id1, id2);
		if (mode==PROFILES_GDIST) return profileGaussianDistanceWeight(id1, id2);
		if (mode==PROFILES_CORR) return correlationWeight(id1, id2);
		if (mode==PROFILES_ANGL) return correlationAngleWeight(id1, id2);
		if (mode==PROFILES_CDIST) return correlationDistanceWeight(id1, id2);
		if (mode==PROFILES_CGDIST) return correlationGaussianDistanceWeight(id1, id2);
		if (mode==SCALAR) return scalarDistanceWeight(id1, id2);
		if (mode==BINARY) return binaryWeight(id1, id2);
		return 0.5f;
	}
	private final float scalarDistanceWeight(int id1, int id2) {
		return 1.0f/( basis + Numerics.square( (image[0][id1]-image[0][id2])/imgscale ) );
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
		if (mode==PROFILES_DIST) return profileEuclideanDistance(id1, id2);
		if (mode==PROFILES_GDIST) return profileEuclideanDistance(id1, id2);
		if (mode==PROFILES_CORR) return profileCorrelation(id1, id2);
		if (mode==PROFILES_ANGL) return profileProduct(id1, id2);
		if (mode==PROFILES_CDIST) return profileCorrelation(id1, id2);
		if (mode==PROFILES_CGDIST) return profileCorrelation(id1, id2);
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
	
}
