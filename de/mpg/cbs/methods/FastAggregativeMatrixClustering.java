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
 *  This algorithm builds a hierarchical clustering of connectivity matrix data
 *	The connectivity is assumed to be always >0; a value of -1 makes no connection
 *
 *	The key component is an array structure that contains the association weights for all pairs of connected points.
 *	cluster as well as the degree in the first element of the list. Once constructed, the association array is fully
 *	abstract (no contextual information remains).
 *
 *	@version    September 2012
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class FastAggregativeMatrixClustering {
		
	// data buffers
	private 	float[][]	image;  			// original data
	private 	int			nd;   		// image dimensions
	private		float		scale;
	private		float		basis;
	private		float 		pvalue;
	
	private		byte		mode;
	private	static final	byte	GAUSS_DIST = 10;
	private	static final	byte	NPARAM_DIST = 12;
	
	private static final double PI2 = Math.PI/2.0;
	private static final double SQPI = Math.sqrt(Math.PI);
	private static final float INF = 1e12f;
	
	// algorithm quantities
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
	public FastAggregativeMatrixClustering(float[][] img_, int nlb_, float sca_, float bas_, String mod_) {
		image = img_;
		
		nlb = nlb_;
		
		scale = sca_;
		basis = bas_;
		pvalue = sca_;
		
		if (mod_.equals("gauss-avg-distance")) 			mode = GAUSS_DIST;
		else if (mod_.equals("gauss-mmx-distance")) 		mode = GAUSS_DIST;
		
		else if (mod_.equals("np-avg-distance")) 			mode = NPARAM_DIST;
		else if (mod_.equals("np-mmx-distance")) 			mode = NPARAM_DIST;
		
	}

	final public void finalize() {
		image = null;
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
		
		Triple[] trp = new Triple[nlb-1];
		
		if (debug) System.out.println("first pass");

		for (int n1=0;n1<nlb;n1++) {
			// build the weight array
			int nb=0;
			for (int n2=0;n2<nlb;n2++) if (n2!=n1) if (image[n1][n2]>-1) {
				// note: here we can remove unnecessary nodes if desired
				trp[nb] = new Triple(n2+1, associationWeight(n1+1,n2+1));
				nb++;
			}
			
			// link in all conserved directions
			Triple[] node = new Triple[nb+1];
			node[0] = new Triple(n1+1, 1.0f);
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
			assoc[n1+1] = node;
			
			// store the latest active index for everything
			latest[n1+1] = n1+1;
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
	public final void initAllEdgeWeightsMinimax() {
		
		if (debug) System.out.println("-- weight initialization --");
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		Triple[] trp = new Triple[nlb-1];
		
		if (debug) System.out.println("first pass");

		for (int n1=0;n1<nlb;n1++) {
			// build the weight array
			int nb=0;
			for (int n2=0;n2<nlb;n2++) if (n2!=n1) {
				// note: here we can remove unnecessary nodes if desired
				trp[nb] = new Triple(n2+1, associationWeight(n1+1,n2+1));
				nb++;
			}
				
			// link in all conserved directions
			Triple[] node = new Triple[nb+1];
			node[0] = new Triple(n1+1, 1.0f);
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
			assoc[n1+1] = node;
			
			// store the latest active index for everything
			latest[n1+1] = n1+1;
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
	public final float[] computeWeightMetrics() {
		
		if (debug) System.out.println("-- metric pre-processing --");
		
		float[] distance = new float[nlb/2*nlb];
		int nb = 0;
		for (int n1=0;n1<nlb;n1++) {
			for (int n2=n1+1;n2<nlb;n2++) {
				if (image[n1][n2]>-1) {
					distance[nb] = dataDistance(n1+1,n2+1);
					nb++;
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
		
		if (mode==GAUSS_DIST) {
			scale = (float)(0.5*hist.mean()*SQPI);
		}
		distribution = hist;
		distribution.normalizeToMaximum();
		//distribution.normalizeToNonZeroMaximum();
		//distribution.inverseCumulativeHistogram();
		//System.out.println(distribution.printHistogram());
		System.out.println("adjusted scale: "+scale);
		
		return distance;
	}

	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringAverage(int k0, boolean firststop, int maxlength) {
		
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
					/* irrelevant here: skipped
					// set the new location parameters
					ngbcluster.clear();
					if (lbest>nlb && lpair>nlb && clusterPos[lbest].length>clusterPos[lpair].length) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							if (mask[invertlabeling[clusterPos[lbest][n]]]) 
								ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][clusterPos[lbest].length+n] = clusterPos[lpair][n];
							int xyz = invertlabeling[clusterPos[lpair][n]];
							if (mask[xyz]) {
								for (int k = 0; k<6; k++) {
									int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
									if (ngbcluster.get(xyzn)) {
										boundaryScore[xyz] = bcost;
										boundaryScore[xyzn] = bcost;
									}
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb && lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][n] = clusterPos[lpair][n];
							if (mask[invertlabeling[clusterPos[lpair][n]]])
								ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][clusterPos[lpair].length+n] = clusterPos[lbest][n];
							int xyz = invertlabeling[clusterPos[lbest][n]];
							if (mask[xyz]) {
								for (int k = 0; k<6; k++) {
									int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
									if (ngbcluster.get(xyzn)) {
										boundaryScore[xyz] = bcost;
										boundaryScore[xyzn] = bcost;
									}
								}
							}
						}
						clusterPos[lbest] = null;
						clusterPos[lpair] = null;
					} else if (lbest>nlb) {
						clusterPos[id] = new int[1+clusterPos[lbest].length];
						for (int n=0;n<clusterPos[lbest].length;n++) {
							clusterPos[id][n] = clusterPos[lbest][n];
							if (mask[invertlabeling[clusterPos[lbest][n]]])
								ngbcluster.set(invertlabeling[clusterPos[lbest][n]]);
						}
						clusterPos[id][clusterPos[lbest].length] = lpair;
						int xyz = invertlabeling[lpair];
						if (mask[xyz]) {
							boundaryScore[xyz] = bcost;
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyzn] = bcost;
								}
							}
						}
						clusterPos[lbest] = null;
					} else if (lpair>nlb) {
						clusterPos[id] = new int[clusterPos[lpair].length+1];
						clusterPos[id][0] = lbest;
						for (int n=0;n<clusterPos[lpair].length;n++) {
							clusterPos[id][1+n] = clusterPos[lpair][n];
							if (mask[invertlabeling[clusterPos[lpair][n]]])
								ngbcluster.set(invertlabeling[clusterPos[lpair][n]]);
						}
						int xyz = invertlabeling[lbest];
						if (mask[xyz]) {
							boundaryScore[xyz] = bcost;
							for (int k = 0; k<6; k++) {
								int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
								if (ngbcluster.get(xyzn)) {
									boundaryScore[xyzn] = bcost;
								}
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
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						// use the same criterion than for stopping? makes sense, but slows down the process
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/(nlb*nlb)) );
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
					/* skipped
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
	public final int[] exportClusteringSequential(int maxlb) {
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
		int[] tmp = new int[nlb];
		for (int n=0;n<nlb;n++) {
			tmp[n] = latest[n+1];
		}
		return tmp;
	}
	
	// computes the association parameter for the image
	private final float associationWeight(int id1, int id2) {
		if (mode==GAUSS_DIST) return gaussDistanceWeight(id1, id2);
		if (mode==NPARAM_DIST) return nonParamDistanceWeight(id1, id2);
		return 0.5f;
	}
	private final float gaussDistanceWeight(int id1, int id2) {
		return (float)FastMath.exp( -0.5*image[id1-1][id2-1]/(scale*scale) );
	}
	private final float nonParamDistanceWeight(int id1, int id2) {
		return distribution.getHistogramCount(image[id1-1][id2-1]);
	}

	private final float dataDistance(int id1, int id2) {
		if (mode==GAUSS_DIST) return gaussDistance(id1, id2);
		if (mode==NPARAM_DIST) return gaussDistance(id1, id2);
		return 0.5f;
	}
	
	private final float gaussDistance(int id1, int id2) {
		return image[id1-1][id2-1];
	}
	
}
