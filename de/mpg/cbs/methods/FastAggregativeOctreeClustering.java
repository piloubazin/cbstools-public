package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.util.FastMath;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm builds a hierarchical clustering of octree data
 *
 *	The key component is an array structure that contains the association weights for all pairs of connected points.
 *	By default, we use 6-connectivity (can be easily changed if desired). The array stores the self-weights for each
 *	cluster as well as the degree in the first element of the list. Once constructed, the association array is fully
 *	abstract (no contextual information remains).
 *
 *	@version    March 2015
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class FastAggregativeOctreeClustering {
		
	// data buffers
	private 	OctreeCluster	octree;  			// original data
	private 	int			nx,ny,nz;   		// image dimensions
	private		boolean[]	mask;
	private		int			connect;
	private		float 		pvalue;
	private		float 		globalvar;
	
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
	private		Cluster[]						list;
	private		Ngb[][]							assoc;
	private		BitSet							active;
	private		BinaryHeapPair					maxtree;
	private		float[]							cost;
	private		int[]							latest;
	private		int[][]							clusterPos;
	private		int[]							invertlabeling;
	private		float[]							boundaryScore;
	
	static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// simple id, value class
	private static class Cluster {
		public int id;
		public float size;
		public float sum;
		public float sqd;
		// public double det; ??
		
		public Cluster(int id) {
			this.id = id;
			this.size = 1;
			this.sum = 0.0f;
			this.sqd = 0.0f;
		}
		public Cluster(int id, float sum, float size) {
			this.id = id;
			this.sum = sum;
			this.size = size;
			this.sqd = 0.0f;
		}
		public Cluster(int id, float sum, float sqd, float size) {
			this.id = id;
			this.sum = sum;
			this.sqd = sqd;
			this.size = size;
		}
		
	}
	
	// simple id, value class
	private static class Ngb {
		public int id;
		public float size;
		public float delta;
		
		public Ngb(int id, float delta) {
			this.id = id;
			this.delta = delta;
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
	public FastAggregativeOctreeClustering(OctreeCluster oct_, boolean[] msk_,
											int nx_, int ny_, int nz_,
											float sqd_, float pval_) {
		octree = oct_;
		mask = msk_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		globalvar = sqd_;
		
		pvalue = pval_;
		connect = 6;
		
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};

		// init: find the labels for used data points
		// pb: destroys spatial information (i.e. other data must be transferred in that space too)
		if (debug) { System.out.println("get octree-based labels"); System.out.flush(); }

		//labeling = octree.exportMaskedLabelingToImage(mask,nx,ny,nz);
		labeling = new int[nx*ny*nz];
		nlb = octree.createMaskedLabelingFromScale(labeling, mask, nx,ny,nz);
		
		if (debug) { System.out.println("label count: "+nlb); System.out.flush(); }

		//nlb = ObjectLabeling.countLabels(labeling, nx, ny, nz)-1;
		invertlabeling = new int[nlb+1];
		invertlabeling[0] = -1;
		int lb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (invertlabeling[labeling[xyz]]==0) {
				invertlabeling[labeling[xyz]] = xyz;
			}
		}
		
		if (debug) { System.out.println("initial labeling ("+nlb+" clusters)"); System.out.flush(); }
		
	}

	final public void finalize() {
		octree = null;
		labeling = null;
		assoc = null;
		System.gc();
	}
	
	public final int getInitClusterNumber() { return nlb; }
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initClusters() {
		
		list = new Cluster[2*nlb];
		list[0] = new Cluster(0);
		assoc = new Ngb[2*nlb][];
		assoc[0] = null;
		
		latest = new int[2*nlb];
		
		if (debug) { System.out.println("weight computation"); System.out.flush(); }

		float zerosqd = 0.0f;
		
		OctreeCluster.Cell cellX,cellY,cellZ;
		cellX = octree.traverseToLeaf(0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = octree.traverseToLeaf(cellZ,x,y,z);
					int cs = 1 << (cellZ.level);
					int xyz = x+nx*y+nx*ny*z;
					if (labeling[xyz]>0 && list[labeling[xyz]]==null) {
						list[labeling[xyz]] = new Cluster(labeling[xyz], cellZ.sum, cellZ.sqd, octree.cellsizes[cellZ.level]);
					}
					if (z+1>=cellZ.z+cs) cellZ = octree.getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = octree.getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = octree.getXplusNeighbor(cellX);
		}

		if (debug) { System.out.println("init neighbors for clustering"); System.out.flush(); }
		
		
		BitSet counted = new BitSet(nx*ny*nz);
		int maxnb = 0;
		cellX = octree.traverseToLeaf(0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = octree.traverseToLeaf(cellZ,x,y,z);
					int cs = 1 << (cellZ.level);
					int xyz = x+nx*y+nx*ny*z;
					if (labeling[xyz]>0) {
						Ngb[] ngb = new Ngb[(int)(connect*octree.celllengths[cellZ.level]*octree.celllengths[cellZ.level])];
		
						// 6-C
						int nb=0;
						OctreeCluster.Cell cellngb = octree.getXminusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell before in X direction
										int xyzngbij = (int)(xyzngb+nx*i+nx*ny*j+octree.celllengths[cellZ.level]-1);
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x+(int)(octree.celllengths[cellZ.level]-1), cellngb.y+i, cellngb.z+j);
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						cellngb = octree.getXplusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell after in X direction
										int xyzngbij = xyzngb+nx*i+nx*ny*j;
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x, cellngb.y+i, cellngb.z+j);
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						cellngb = octree.getYminusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell before in Y direction
										int xyzngbij = (int)(xyzngb+nx*ny*i+j+nx*(octree.celllengths[cellZ.level]-1));
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x+j, cellngb.y+(int)(octree.celllengths[cellZ.level]-1), cellngb.z+i);
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						cellngb = octree.getYplusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell after in Y direction
										int xyzngbij = xyzngb+nx*ny*i+j;
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x+j, cellngb.y, cellngb.z+i);
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						cellngb = octree.getZminusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell before in Z direction
										int xyzngbij = (int)(xyzngb+i+nx*j+nx*ny*(octree.celllengths[cellZ.level]-1));
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x+i, cellngb.y+j, cellngb.z+(int)(octree.celllengths[cellZ.level]-1));
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						cellngb = octree.getZplusNeighbor(cellZ);
						if (cellngb!=null) {
							int xyzngb = cellngb.x + nx*cellngb.y + nx*ny*cellngb.z;
							if (labeling[xyzngb]>0) {
								// search for all leaf neighbors
								counted.clear();
								counted.set(labeling[xyz]);
								for (int i=0;i<octree.celllengths[cellZ.level];i++) {
									for (int j=0;j<octree.celllengths[cellZ.level];j++){
										// cell before in Z direction
										int xyzngbij = xyzngb+i+nx*j;
										// check if label is already counted before going to the leaf
										if (labeling[xyzngbij]>0 && !counted.get(labeling[xyzngbij])) {
											OctreeCluster.Cell childngb = octree.traverseToLeaf(cellngb, cellngb.x+i, cellngb.y+j, cellngb.z);
											ngb[nb] = new Ngb(labeling[xyzngbij], jsDiv(cellZ.sum, childngb.sum, cellZ.sqd, childngb.sqd, octree.cellsizes[cellZ.level], octree.cellsizes[childngb.level]));
											nb++;
											counted.set(labeling[xyzngbij]);
										}
									}
								}								
							}
						}
						
						// store the values: new instance!
						assoc[labeling[xyz]] = new Ngb[nb];
						for (int n=0;n<nb;n++) {
							assoc[labeling[xyz]][n] = ngb[n];
						}
						
						// debug: check on the max size of neighborhood
						if (nb>maxnb) maxnb = nb;
				
						// store the latest active index for everything
						latest[labeling[xyz]] = labeling[xyz];
					}
					if (z+1>=cellZ.z+cs) cellZ = octree.getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = octree.getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = octree.getXplusNeighbor(cellX);
		}
		if (verbose) { System.out.println("max. neighborhood size: "+maxnb); System.out.flush(); }
		
		if (debug) { System.out.println("done"); System.out.flush(); }
		
		
		// second pass for delta, cost: not needed
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		
		return;
	}
	
	public final void hierarchicalClustering(int k0, boolean firststop) {
		
		if (debug) System.out.println("hierarchical clustering : normal distribution distances");
		
		// 1. Build a binary tree for the delta values
		maxtree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		int ntree=0;
		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Ngb[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>0) {
				int best=0;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=1;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (node[best].delta>0) {
					//if (debug) System.out.print(""+node[best].delta+",");
					maxtree.addValue(node[best].delta, list[lb].id, best);
					ntree++;
				//}
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building (start: "+ntree+")");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		
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
		
		float bcost = 0;
		while (maxtree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			bcost = maxtree.getFirst(); 
			int  	lbest = maxtree.getFirstId1();
			int 	nbest = maxtree.getFirstId2();
			maxtree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Cluster bNodeC = list[lbest];
				Ngb[] bNodeN = assoc[lbest];
				int	lpair = bNodeN[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Cluster pNodeC = list[lpair];
					Ngb[] pNodeN = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values : depends if same octree init scale
					//Ngb[] aNodeN = new Ngb[bNodeN.length+pNodeN.length-2];
					// better to define it later on
					
					// sum(uv) = sum(u) + sum(v)
					// sq2(uv) = sq2(u) + sq2(v) + 1/n(uv) [n(u)/n(v)*m(v) - n(v)/n(u)*m(u)]^2
					// n(uv) = n(u) + n(v)
					
					float sum = bNodeC.sum+pNodeC.sum;
					float size = bNodeC.size+pNodeC.size;
					float factor = bNodeC.size*pNodeC.size/size;
					float sqd = bNodeC.sqd + pNodeC.sqd 
									+ factor*(bNodeC.sum/bNodeC.size-pNodeC.sum/pNodeC.size)
											*(bNodeC.sum/bNodeC.size-pNodeC.sum/pNodeC.size);
					Cluster aNodeC = new Cluster(id, sum,  sqd, size);
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNodeN.length+pNodeN.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					ngbList.clear();
					for (int n=0;n<bNodeN.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNodeN[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[list[lbn].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNodeN[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=0;n<pNodeN.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNodeN[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[list[lbn].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNodeN[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					Ngb[] aNodeN = new Ngb[ngbList.cardinality()];
					int l=0;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNodeN[l] = new Ngb(lbn, 0.0f);
						l++;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNodeN.length) {
						Ngb[] tmp = new Ngb[l];
						for (int n=0;n<l;n++) tmp[n] = aNodeN[n];
						aNodeN = tmp;
					}
					/* optional
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
					*/
					
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					//System.out.print(">");
					
					for (int n=0; n<aNodeN.length; n++) {
						Cluster wngbC = list[aNodeN[n].id];
						
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNodeC.size+wngbC.size)/((connect+1.0f)*nlb)) );
						
						// build the divergence metric (with built-in threshold)
						float metric = (float)jsDiv(aNodeC.sum, wngbC.sum, 
													aNodeC.sqd, wngbC.sqd, 
													aNodeC.size, wngbC.size);
						
						aNodeN[n].delta = wsize*metric;
						
						//System.out.print("m="+metric+"; ");
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						float avgsqd = (float)FastMath.sqrt(aNodeC.sqd/aNodeC.size);
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "
																	+"sqd= "+avgsqd
																	+", n= "+aNodeC.size
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNodeN.length+", "+pNodeN.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						System.out.flush();
					
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
					// not a correct stopping criterion								
					if (verbose) if (bcost<0 && first) {
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+", n= "+aNodeC.size
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNodeN.length+", "+pNodeN.length+")");	
					   first=false;
					   if (firststop) stop = true;
					}								
					// add the new values to list, binary tree
					list[id] = aNodeC;
					assoc[id] = aNodeN;
					active.set(id, true);
					
					if (aNodeN.length>0) {
						int best=0;
						for (int b=1;b<aNodeN.length;b++) {
							if (aNodeN[b].delta>aNodeN[best].delta) best = b;
						}
						maxtree.addValue(aNodeN[best].delta, id, best);
					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					list[lbest].id = id;
					list[lpair].id = id;
					
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
		if (debug) System.out.println("completed ("+nclusters+"| "+bcost+"; "+maxtree.getCurrentSize()+")");
		
	}
	
	// computes the association parameter for the image
    private float jsDiv(float sum1, float sum2,
						float sqd1, float sqd2,
						float size1, float size2) {
       
		double v1 = Numerics.max(sqd1/size1, globalvar);
		double v2 = Numerics.max(sqd2/size2, globalvar);
		double v12 = 0.25*( v1 + v2 );
		
		double js = 0.25*Numerics.square(sum1/size1-sum2/size2)/v12 - 0.5*FastMath.log(v1/v12*v2/v12) + FastMath.log(2.0);
		// Gaussian model on the distances
    	return (float)FastMath.exp(-0.5*js) - pvalue;
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
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringSequential(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = list[n].id;
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
	
	// given a clustering result, update the labeling
	public final String clusteringStatisticsSequential(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = list[n].id;
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
	public final float[][][] filteredClusters(int maxlb, float minratio, float maxratio) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = list[n].id;
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
	public final float[][][] exportClusterBoundaryScore() {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			tmp[x][y][z] = boundaryScore[xyz];
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringMean(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = list[n].id;
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
				tmp[x][y][z] = list[latest[labeling[xyz]]].sum/list[latest[labeling[xyz]]].size;
			}
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringStdev(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = list[n].id;
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
				tmp[x][y][z] = (float)FastMath.sqrt(list[latest[labeling[xyz]]].sqd/Numerics.max(list[latest[labeling[xyz]]].size-1.0f,1.0f));
			}
		}
		return tmp;
	}
		
}
