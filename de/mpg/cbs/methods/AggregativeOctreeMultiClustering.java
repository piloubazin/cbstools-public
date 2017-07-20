package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.analysis.function.Gaussian;

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
 
public class AggregativeOctreeMultiClustering {
		
	// data buffers
	private 	OctreeMultiClusterSimplification		octree;  			// original data
	private		int			connect;
	private		float 		pvalue;
	private		float[] 	globalvar;
	
	private static final double PI2 = Math.PI/2.0;
	private static final double SQPI = FastMath.sqrt(Math.PI);
	private static final float INF = 1e12f;
	private static final double LN2 = FastMath.log(2.0);
	private static final double SQLN2 = FastMath.sqrt(FastMath.log(2.0));
	
	private int	nx, ny, nz, nt;
	private int[] xoff;
   	private int[] yoff;
	private int[] zoff;

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
	
	private	static	int metric;
	private	static 	final int 		JSDIV_SUM = 101;
	private	static 	final int 		JSDIV_MAX = 102;
	private	static 	final int 		JSDIV_EXACT = 103;
	private	static 	final int 		JSDIV_TABLE = 104;
	private		JensenShannonDivTable			jsdtable;
	
	private	static	int varupdate;
	private	static 	final int 		VARMAX = 201;
	private	static 	final int 		VARLIN = 202;
	private	static 	final int 		VARZERO = 203;
	private	static 	final int 		VARCONST = 204;
	private int	clustersize;

	
	
	static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// simple id, value class
	private static class Cluster {
		public int id;
		public int size;
		public float[] sum;
		public float[] sqd;
		// public double det; ??
		
		public Cluster(int id, int n) {
			this.id = id;
			this.size = 1;
			this.sum = new float[n];
			this.sqd = new float[n];
		}
		public Cluster(int id, float[] sum, int size) {
			this.id = id;
			this.size = size;
			this.sum = sum;
			this.sqd = new float[sum.length];
		}
		public Cluster(int id, float[] sum, float[] sqd, int size) {
			this.id = id;
			this.size = size;
			this.sum = sum;
			this.sqd = sqd;
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
	public AggregativeOctreeMultiClustering(OctreeMultiClusterSimplification oct_, float[] sqd_, float pval_, String metrictype, String varup, int clustsize, int conn_) {
		octree = oct_;
		
		globalvar = sqd_;
		
		pvalue = pval_;
		connect = conn_;
		
		if (metrictype.equals("joint Jensen-Shannon")) metric = JSDIV_SUM;
		else if (metrictype.equals("max Jensen-Shannon")) metric = JSDIV_MAX;
		else if (metrictype.equals("exact Jensen-Shannon")) metric = JSDIV_EXACT;
		else if (metrictype.equals("table Jensen-Shannon")) metric = JSDIV_TABLE;
		
		//if (metric==JSDIV_TABLE) jsdtable = new JensenShannonDivTable(3.0, 0.01, 0.01);
		if (metric==JSDIV_TABLE) jsdtable = new JensenShannonDivTable(6.0, 0.01, 0.01);
		
		if (varup.equals("Maximum")) varupdate = VARMAX;
		else if (varup.equals("Linear")) varupdate = VARLIN;
		else if (varup.equals("Zero")) varupdate = VARZERO;
		else varupdate = VARCONST;
		
		clustersize = clustsize;
	}

	final public void finalize() {
		octree = null;
		labeling = null;
		assoc = null;
		System.gc();
	}
	
	public final int getInitClusterNumber() { return nlb; }
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initFirstClusters(int level) {
		
		nx = octree.getDimX(level);
		ny = octree.getDimY(level);
		nz = octree.getDimZ(level);
		nt = octree.getDimT();
		
		// neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1,  0,   0, 0,          0,  1, -1,   1,  -1,     0,     0,      0,      0,     1,      1,    -1,     -1,     1,    -1,     1,    -1,      1,     -1,      1,     -1};
		yoff = new int[]{0,  0, nx, -nx, 0,          0, nx, nx, -nx, -nx,    nx,   -nx,     nx,    -nx,     0,      0,     0,      0,    nx,    nx,   -nx,   -nx,     nx,     nx,    -nx,    -nx};
		zoff = new int[]{0,  0,  0,   0, nx*ny, -nx*ny,  0,  0,   0,   0, nx*ny, nx*ny, -nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, nx*ny, nx*ny, nx*ny, -nx*ny, -nx*ny, -nx*ny, -nx*ny};
		
		float[] ngbfactor = new float[connect];
		for (int c=0;c<connect;c++) {
			if (c<6) ngbfactor[c] = 1.0f;
			else if (c<18) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(2.0f);
			else if (c<26) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(3.0f);
		}
		
		// init: find the labels for used data points
		// pb: destroys spatial information (i.e. other data must be transferred in that space too)
		if (debug) { System.out.println("get octree-based labels"); System.out.flush(); }

		//labeling = octree.exportMaskedLabelingToImage(mask,nx,ny,nz);
		nlb = octree.generateLabelsAt(level);
		labeling = octree.getLbl(level);
		
		if (debug) { System.out.println("initial labeling ("+nlb+" clusters)"); System.out.flush(); }

		list = new Cluster[2*nlb];
		list[0] = new Cluster(0,nt);
		assoc = new Ngb[2*nlb][];
		assoc[0] = null;
		
		latest = new int[2*nlb];
		
		if (debug) { System.out.println("weight computation"); System.out.flush(); }

		float zerosqd = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labeling[xyz]>0 && list[labeling[xyz]]==null) {
				list[labeling[xyz]] = new Cluster(labeling[xyz], octree.getSum(level)[xyz], octree.getSqd(level)[xyz], octree.getNpt(level)[xyz]);
			}
		}
		if (debug) { System.out.println("init neighbors for clustering"); System.out.flush(); }
		
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labeling[xyz]>0) {
				Ngb[] ngb = new Ngb[connect];
		
				// 6-C or more: change connect value
				int nb=0;
				for (int n=0;n<connect;n++) {
					int xyzngb = xyz+xoff[n]+yoff[n]+zoff[n];
					if (x+xoff[n]>=0 && x+xoff[n]<nx 
						&& nx*y+yoff[n]>=0 && nx*y+yoff[n]<nx*ny 
						&& nx*ny*z+zoff[n]>=0 && nx*ny*z+zoff[n]<nx*ny*nz) {
						if (labeling[xyzngb]>0) {
							ngb[nb] = new Ngb(labeling[xyzngb], ngbfactor[n]*jsDiv(octree.getSum(level)[xyz], octree.getSum(level)[xyzngb], 
																				  octree.getSqd(level)[xyz], octree.getSqd(level)[xyzngb], 
																				  octree.getNpt(level)[xyz], octree.getNpt(level)[xyzngb] ) );
							nb++;						
						}
					}
				}
								
				// store the values: new instance!
				assoc[labeling[xyz]] = new Ngb[nb];
				for (int n=0;n<nb;n++) {
					assoc[labeling[xyz]][n] = ngb[n];
				}
						
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			}
		}
		if (debug) { System.out.println("done"); System.out.flush(); }
		
		
		// second pass for delta, cost: not needed
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		
		return;
	}
	
	// new lists: create all the links from a previous labeling and the current octree representation
	public final void initNextClusters(int level, int prevlb, boolean propagate, float approx) {
		
		nx = octree.getDimX(level);
		ny = octree.getDimY(level);
		nz = octree.getDimZ(level);
		nt = octree.getDimT();
		
		// neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1,  0,   0, 0,          0,  1, -1,   1,  -1,     0,     0,      0,      0,     1,      1,    -1,     -1,     1,    -1,     1,    -1,      1,     -1,      1,     -1};
		yoff = new int[]{0,  0, nx, -nx, 0,          0, nx, nx, -nx, -nx,    nx,   -nx,     nx,    -nx,     0,      0,     0,      0,    nx,    nx,   -nx,   -nx,     nx,     nx,    -nx,    -nx};
		zoff = new int[]{0,  0,  0,   0, nx*ny, -nx*ny,  0,  0,   0,   0, nx*ny, nx*ny, -nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, nx*ny, nx*ny, nx*ny, -nx*ny, -nx*ny, -nx*ny, -nx*ny};
		
		float[] ngbfactor = new float[connect];
		for (int c=0;c<connect;c++) {
			if (c<6) ngbfactor[c] = 1.0f;
			else if (c<18) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(2.0f);
			else if (c<26) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(3.0f);
		}
				
		// init: find the labels for used data points
		// pb: destroys spatial information (i.e. other data must be transferred in that space too)
		if (debug) { System.out.println("get octree-based labels"); System.out.flush(); }

		// first: propagate maps from previous level
		octree.generateNextLevelMaps(level);

		// labeling: import from previous scale and add new scale elements
		if (propagate) nlb = octree.competitionWithJensenShannonDistance(level, prevlb, approx);
		else nlb = octree.generateNextLabelsAt(level, prevlb);
		labeling = octree.getLbl(level);
		
		
		if (debug) { System.out.println("initial labeling ("+nlb+" clusters)"); System.out.flush(); }

		list = new Cluster[2*nlb];
		list[0] = new Cluster(0,nt);
		assoc = new Ngb[2*nlb][];
		assoc[0] = null;
		
		latest = new int[2*nlb];
		
		if (debug) { System.out.println("weight computation"); System.out.flush(); }

		float zerosqd = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (labeling[xyz]>0 && list[labeling[xyz]]==null) {
				if (labeling[xyz]>prevlb) {
					// new labels from current scale
					list[labeling[xyz]] = new Cluster(labeling[xyz], octree.getSum(level)[xyz], octree.getSqd(level)[xyz], octree.getNpt(level)[xyz]);
				} else {
					// labels from the previous scale: use npt instead of size	
					list[labeling[xyz]] = new Cluster(labeling[xyz], octree.getSum(level)[xyz], octree.getSqd(level)[xyz], octree.getNpt(level)[xyz]);
				}
			}
		}
		if (debug) { System.out.println("init neighbors for clustering"); System.out.flush(); }
		
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			//if (labeling[xyz]>0) {
			// build relationships only for the new labels
			if (labeling[xyz]>prevlb) {
				Ngb[] ngb = new Ngb[connect];
		
				// 6-C
				int nb=0;
				for (int n=0;n<connect;n++) {
					int xyzngb = xyz+xoff[n]+yoff[n]+zoff[n];
					if (x+xoff[n]>=0 && x+xoff[n]<nx 
						&& nx*y+yoff[n]>=0 && nx*y+yoff[n]<nx*ny 
						&& nx*ny*z+zoff[n]>=0 && nx*ny*z+zoff[n]<nx*ny*nz) {
						if (labeling[xyzngb]>prevlb) {
							ngb[nb] = new Ngb(labeling[xyzngb], ngbfactor[n]*jsDiv(octree.getSum(level)[xyz], octree.getSum(level)[xyzngb], 
																				   octree.getSqd(level)[xyz], octree.getSqd(level)[xyzngb], 
																				   octree.getNpt(level)[xyz], octree.getNpt(level)[xyzngb] ) );
							nb++;						
						} else 
						// allow relationships to previous labels of course
						if (labeling[xyzngb]>0) {
							ngb[nb] = new Ngb(labeling[xyzngb], ngbfactor[n]*jsDiv(octree.getSum(level)[xyz], octree.getSum(level)[xyzngb], 
																				   octree.getSqd(level)[xyz], octree.getSqd(level)[xyzngb], 
																				   octree.getNpt(level)[xyz], octree.getNpt(level)[xyzngb] ) );
							nb++;						
						}
					}
				}
								
				// store the values: new instance!
				assoc[labeling[xyz]] = new Ngb[nb];
				for (int n=0;n<nb;n++) {
					assoc[labeling[xyz]][n] = ngb[n];
				}
						
				// store the latest active index for everything
				latest[labeling[xyz]] = labeling[xyz];
			} else {
				assoc[labeling[xyz]] = new Ngb[0];
			}
		}
		if (debug) { System.out.println("done"); System.out.flush(); }
		
		
		// second pass for delta, cost: not needed
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		
		return;
	}
	
	public final int hierarchicalClustering(int k0, boolean recompute, float recomputeratio, float sizeratio) {
		
		float mincost = - 0.1f*pvalue;
		
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
			//if (node!=null && node.length>0) {
			if (node.length>0) {
				int best=0;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=1;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				if (node[best].delta>mincost) {
					//if (debug) System.out.print(""+node[best].delta+",");
					maxtree.addValue(node[best].delta, list[lb].id, best);
					ntree++;
				} else {
					if (debug) System.out.print("o");	
				}
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
		//clusterPos = new int[2*nlb][];

		// look for joint boundaries
		//BitSet ngbcluster = new BitSet(nx*ny*nz);
		//boundaryScore = new float[nx*ny*nz];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		//boolean first = true;
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
			
				// update the link label? yes if recomputing, no otherwise
				if (!active.get(lpair)) {
					// make sure it's the most up-to-date version of the label
					tracks.reset();
					while (!active.get(lpair)) {
						tracks.add(lpair);
						lpair = latest[list[lpair].id];
						//lbn = assoc[lbn][0].id;
						depthsearchtime++;
					}
					// update the id if not up to date
					for (int lb=0;lb<tracks.last;lb++)
						latest[tracks.val[lb]] = lpair;
						//assoc[tracks.val[lb]][0].id = lbn;
					latest[bNodeN[nbest].id] = lpair;
					//lbn = latest[lbn];
				}
				if (active.get(lpair)) {
					Cluster pNodeC = list[lpair];
					Ngb[] pNodeN = assoc[lpair];

					// recompute score (?) if lower re-insert into the tree and skip
					boolean merge = true;
					if (recompute && bcost>=0) {
						// build the divergence metric (with built-in threshold and geometric factor)
						float metric = jsDiv(bNodeC.sum, pNodeC.sum, 
											 bNodeC.sqd, pNodeC.sqd, 
											 bNodeC.size, pNodeC.size);
						
						if (metric<0.0f) {
							// negative score: throw away
							merge = false;
						} else if (metric<recomputeratio*bcost) {
							// much lower score: send back into the tree
							maxtree.addValue(metric, lbest, nbest);
							merge = false;
						}
					} else if (recompute) {
						// build the divergence metric (with built-in threshold and geometric factor)
						float metric = jsDiv(bNodeC.sum, pNodeC.sum, 
											 bNodeC.sqd, pNodeC.sqd, 
											 bNodeC.size, pNodeC.size);
						
						if (metric<0.0f) {
							// stopping criterion								
							System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNodeN.length+", "+pNodeN.length+")");	
							stop = true;
							merge = false;
						} else {
							// now positive: continue
							merge = true;
						}
					} else {
						// stopping criterion								
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNodeN.length+", "+pNodeN.length+")");	
					   stop = true;
					   merge = false;
					}

	
					if (merge) {
						
						// do the regular merging
			
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
						
						int size = bNodeC.size+pNodeC.size;
						float factor = bNodeC.size*pNodeC.size/size;
						
						float[] sum = new float[nt];
						float[] sqd = new float[nt];
						for (int t=0;t<nt;t++) {
							sum[t] = bNodeC.sum[t] + pNodeC.sum[t];
							sqd[t] = bNodeC.sqd[t] + pNodeC.sqd[t] 
										+ factor*(bNodeC.sum[t]/bNodeC.size-pNodeC.sum[t]/pNodeC.size)
												*(bNodeC.sum[t]/bNodeC.size-pNodeC.sum[t]/pNodeC.size);
						}
						Cluster aNodeC = new Cluster(id, sum,  sqd, size);
						
						// new mixing weights
						// w(uv,x) = w(u,x) + w(v,x)
						// must check if the links still exist, not duplicate
						//if (bNodeN!=null && pNodeN!=null)
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
									// check for size : if the link is smaller than a ratio of current size, drop the connection
									Cluster ngbtest = list[lbn];
									// keep the connection only if neighbor is big enough
									if (ngbtest.size >= sizeratio*aNodeC.size) {
										// add to the label
										ngbList.set(lbn, true);
									}
									// if neighbor had no connection because of size, and this changes: add neighbor connection
									if (bNodeC.size < sizeratio*ngbtest.size && aNodeC.size >= sizeratio*ngbtest.size) {
										// compute score
										float metric = (float)jsDiv(aNodeC.sum, ngbtest.sum, aNodeC.sqd, ngbtest.sqd, aNodeC.size, ngbtest.size);
										//if (metric>0) {
										// must keep also links with negative scores as they can be changed
											Ngb[] tmp = new Ngb[assoc[lbn].length+1];
											for (int na=0;na<assoc[lbn].length;na++) tmp[na] = assoc[lbn][na];
											tmp[assoc[lbn].length] = new Ngb(id,metric);
											assoc[lbn] = tmp;
										//}
									}
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
									Cluster ngbtest = list[lbn];
									// keep the connection only if neighbor is big enough
									if (ngbtest.size >= sizeratio*aNodeC.size) {
										// add to the label
										ngbList.set(lbn, true);
									}
									// if neighbor had no connection because of size, and this changes: add neighbor connection
									if (pNodeC.size < sizeratio*ngbtest.size && aNodeC.size >= sizeratio*ngbtest.size) {
										// compute score
										float metric = jsDiv(aNodeC.sum, ngbtest.sum, aNodeC.sqd, ngbtest.sqd, aNodeC.size, ngbtest.size);
										//if (metric>0) {
										// must keep also links with negative scores as they can be changed
											Ngb[] tmp = new Ngb[assoc[lbn].length+1];
											for (int na=0;na<assoc[lbn].length;na++) tmp[na] = assoc[lbn][na];
											tmp[assoc[lbn].length] = new Ngb(id,metric);
											assoc[lbn] = tmp;
										//}
									}
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
	
						// make sure we don't have extra empty values (??)
						//aNode.trimToSize();
						if (l<aNodeN.length) {
							Ngb[] tmp = new Ngb[l];
							for (int n=0;n<l;n++) tmp[n] = aNodeN[n];
							aNodeN = tmp;
						}
						
						// new deltas : no need of the u,v values anymore
						// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
						//System.out.print(">");
						
						for (int n=0; n<aNodeN.length; n++) {
							Cluster wngbC = list[aNodeN[n].id];
							
							// build the divergence metric (with built-in threshold)
							float metric = jsDiv(aNodeC.sum, wngbC.sum, 
												 aNodeC.sqd, wngbC.sqd, 
												 aNodeC.size, wngbC.size);
							
							aNodeN[n].delta = metric;
							
							//System.out.print("m="+metric+"; ");
						}
						
						/// probably not needed
						// recompute the degree? (for averaged links)
						cost[iter] = bcost;
	
						if (verbose) if (iter%(Numerics.max(1,nlb/100))==0) {
							long newtime = System.currentTimeMillis();
							nstep++;
							float avgsqd = 0.0f;
							for (int t=0;t<nt;t++) avgsqd += (float)FastMath.sqrt(aNodeC.sqd[t]/aNodeC.size);
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
						/*
						// stopping criterion								
						if (bcost<0) {
						   System.out.println(iter+" / "+nclusters+": c= "+bcost
																		+", n= "+aNodeC.size
																		+" ("+maxtree.getCurrentSize()+"|"
																		+bNodeN.length+", "+pNodeN.length+")");	
						   //first=false;
						   //if (firststop) stop = true;
						   stop = true;
						}
						*/
						// add the new values to list, binary tree
						list[id] = aNodeC;
						assoc[id] = aNodeN;
						active.set(id, true);
						
						if (aNodeN.length>0) {
							int best=0;
							for (int b=1;b<aNodeN.length;b++) {
								if (aNodeN[b].delta>aNodeN[best].delta) best = b;
							}
							// always add, or stop before the first negative?
							if (aNodeN[best].delta>mincost)
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
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+"| "+bcost+"; "+maxtree.getCurrentSize()+")");
		
		return nclusters;
	}
	
	// computes the association parameter for the image
    private float jsDiv(float[] sum1, float[] sum2,
						float[] sqd1, float[] sqd2,
						float npt1, float npt2) {
       
		double js = 0.0;
		for (int t=0;t<nt;t++) {	
			double v1=0,v2=0;
			if (varupdate==VARMAX) {
				// strong prior influence
				v1 = Numerics.max(sqd1[t]/npt1, globalvar[t]);
				v2 = Numerics.max(sqd2[t]/npt2, globalvar[t]);
			} else if (varupdate==VARLIN) {
				// weaker priors: go down with number of samples
				v1 = (sqd1[t] + clustersize*globalvar[t])/(npt1 + clustersize);
				v2 = (sqd2[t] + clustersize*globalvar[t])/(npt2 + clustersize);
			} else if (varupdate==VARZERO) {
				// prior only for t=0
				if (sqd1[t]==0) v1 = globalvar[t]; else v1 = sqd1[t]/npt1;
				if (sqd2[t]==0) v2 = globalvar[t]; else v2 = sqd2[t]/npt2;
			} else if (varupdate==VARCONST) {
				// prior only for t=0
				v1 = globalvar[t];
				v2 = globalvar[t];
			}
			
			if (metric==JSDIV_EXACT) {
				double mu1 = sum1[t]/npt1;
				double mu2 = sum2[t]/npt2;
				double sq1 = FastMath.sqrt(v1);
				double sq2 = FastMath.sqrt(v2);
				
				// here we estimate the JS divergence explicitly
				double xmin = Numerics.min(mu1-3.0*sq1, mu2-3.0*sq2);
				double xmax = Numerics.max(mu1+3.0*sq1, mu2+3.0*sq2);
				Gaussian p1 = new Gaussian(mu1, sq1);
				Gaussian p2 = new Gaussian(mu2, sq2);
				
				double jsdiv = 0.0;
				double step = 0.01*(xmax-xmin);
				for (double x=xmin;x<=xmax;x+=step) {
					double p1x = p1.value(x);
					double p2x = p2.value(x);
					jsdiv += p1x*FastMath.log(2.0*p1x/(p1x+p2x)) + p2x*FastMath.log(2.0*p2x/(p1x+p2x));
				}
				js += jsdiv*step/2.0/(double)nt;
			} else if (metric==JSDIV_TABLE) {
				double mu1 = sum1[t]/npt1;
				double mu2 = sum2[t]/npt2;
				double sq1 = FastMath.sqrt(v1);
				double sq2 = FastMath.sqrt(v2);
				
				js += jsdtable.lookup(mu1,sq1,mu2,sq2)/(double)nt;
			} else {
				// here we compute explicitly the joint distribution
				double v12 = npt1*v1 + npt2*v2 + npt1*npt2/(npt1+npt2)*Numerics.square(sum1[t]/npt1-sum2[t]/npt2);
				double mu12 = sum1[t] + sum2[t];
				
				// normalize
				mu12 /= (npt1+npt2);
				v12 /= (npt1+npt2-1);
	
				double jst = 0;
				if (metric==JSDIV_SUM) {
					jst += 0.5*(Numerics.square(sum1[t]/npt1-mu12)/v12 + FastMath.log( v12/v1 ) );
					jst += 0.5*(Numerics.square(sum2[t]/npt2-mu12)/v12 + FastMath.log( v12/v2 ) );
				
					js += jst/(double)nt/2.0;
				} else {
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum1[t]/npt1-mu12)/v12 + FastMath.log( v12/v1 ) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum2[t]/npt2-mu12)/v12 + FastMath.log( v12/v2 ) ) );
				
					js = Numerics.max(js, jst);
				}
			}
		}
		// geometric progression?
		//float wsize = 1.0f/(1.0f + Numerics.square( (npt1+npt2)/(6.0f*nlb)) );
							
		// Gaussian model on the distances? may not be optimal
		//if (metric==JSDIV_GAUSS_GEOM) return (float)(FastMath.exp(-0.5*js/(2.0*nt)) - pvalue);
		//else if (metric==JSDIV_LIN_GEOM) return (float)(1.0-js/(2.0*nt) - pvalue);
		//else if (metric==JSDIV_GAUSS) return (float)(FastMath.exp(-0.5*js/nt) - pvalue);
		//else return (float)(1.0-js/nt - pvalue);
		
		//return (float)(FastMath.exp(-0.5*js) - pvalue);
		//return (float)(FastMath.exp(-js) - pvalue);
		return (float)(1.0 - FastMath.sqrt(js)/SQLN2 - pvalue);
    }

	// given a clustering result, update the labeling
	public final int exportClusteringResults(int[] lbl, float[][] sum, float[][] sqd, int[] size, int maxlb) {
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
		// get the corresponding values
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (labeling[xyz]>0) {
				for (int t=0;t<nt;t++) {
					sum[xyz][t] = list[latest[labeling[xyz]]].sum[t];
					sqd[xyz][t] = list[latest[labeling[xyz]]].sqd[t];
				}
				size[xyz] = list[latest[labeling[xyz]]].size;
			}
		}
		// remap the labels to [1, Nlb]
		int[] remap = new int[2*nlb];
		int lb=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (labeling[xyz]>0) {
				if (remap[latest[labeling[xyz]]]==0) {
					lb++;
					remap[latest[labeling[xyz]]] = lb;
				}
				lbl[xyz] = remap[latest[labeling[xyz]]];
			} else {
				lbl[xyz] = 0;
			}
		}
		
		return lb;
	}

	// given a clustering result, update the labeling
	public final int smallClusterThreshold(int ncluster) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-ncluster-1; n>0; n--) {
			int lb = list[n].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-ncluster) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// make a visitation list to make sure there's no duplication
		BitSet visited = new BitSet(2*nlb-ncluster);
		double[] scale = new double[ncluster+2];
		int nc=0;
		double maxscale = 0;
		for (int n=1;n<2*nlb-ncluster; n++) {
			if (!visited.get(latest[n])) {
				scale[nc] = list[latest[n]].size;
				if (scale[nc]>maxscale) maxscale = scale[nc];
				nc++;
				visited.set(latest[n]);
			}
		}
		float beta = ImageStatistics.robustExponentialFit(scale, true, ncluster);

		float threshold = beta*(float)FastMath.log(maxscale/beta);
		BasicInfo.displayMessage("cluster threshold "+threshold+"\n");
		
		return Numerics.ceil(threshold);
	}

}
