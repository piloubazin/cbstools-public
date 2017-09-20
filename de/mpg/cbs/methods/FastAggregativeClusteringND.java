package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.util.FastMath;

import Jama.Matrix;
import Jama.LUDecomposition;
import Jama.SingularValueDecomposition;

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
 
public class FastAggregativeClusteringND {
		
	// data buffers
	private 	float[][]	image;  			// original data
	private 	int			nx,ny,nz,nc;   		// image dimensions
	private 	float		rx,ry,rz;   		// image resolutions
	private		boolean[]	mask;
	private		float		basis;
	private		int			connect;
	//private		static 	boolean		simplified;
	//private		String[]	modes = {"profiles", "binary", "scalar", "noise"};
	private		float 		pvalue;
	private		float[][] 	globalcov;
	private		double 		det0;
	private		float		eigratio = 0.99f;
	
	private		static byte		covariance;
	private	static final	byte	SINGLE = 10;
	private	static final	byte	DIAGONAL = 11;
	private	static final	byte	FULL = 12;
	private	static final	byte	PCA = 13;
	
	private		byte		mode;
	private	static final	byte	HOTELLING = 23;
	private	static final	byte	JENSENSHANNON = 24;
	
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
		public float[] mean;
		public float[][] cov;
		public float det;
		// public double det; ??
		
		public Cluster(int id, int n) {
			this.id = id;
			this.size = 1;
			this.mean = new float[n];
			if (covariance==SINGLE) this.cov = new float[1][1];
			else if (covariance==DIAGONAL) this.cov = new float[1][n];
			else this.cov = new float[n][n];
			this.det = 0.0f;
		}
		public Cluster(int id, float[] mean, float size) {
			this.id = id;
			this.mean = mean;
			if (covariance==SINGLE) this.cov = new float[1][1];
			else if (covariance==DIAGONAL) this.cov = new float[1][mean.length];
			else this.cov = new float[mean.length][mean.length];
			this.size = size;
			this.det = 0.0f;
		}
		public Cluster(int id, float[] mean, float[][] cov, float size, float det) {
			this.id = id;
			this.mean = mean;
			this.cov = cov;
			this.size = size;
			this.det = det;
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
	public FastAggregativeClusteringND(float[][] img_, boolean[] msk_,
										int nx_, int ny_, int nz_, int nc_,
										float rx_, float ry_, float rz_,
										float ims_, float bas_, int conn_, String mod_,
										String cov_) {
		image = img_;
		mask = msk_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nc = nc_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;		
		
		pvalue = ims_;
		basis = bas_;
		connect = conn_;
		
		/*
		simplified = simplify_;
		if (simplified) covariance = DIAGONAL;
		else covariance = FULL;
		*/
		if (cov_.equals("single")) 	covariance = SINGLE;
		else if (cov_.equals("diagonal")) 	covariance = DIAGONAL;
		else if (cov_.equals("pca")) 	covariance = PCA;
		else covariance = FULL;
		
		if (mod_.equals("hotelling")) 	mode = HOTELLING;
		else if (mod_.equals("jensen-shannon")) 	mode = JENSENSHANNON;
		
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
	public final void initClusters() {
		
		if (debug) System.out.println("mean covariance estimation");
		
		globalcov = computeGlobalCovariance();
		
		list = new Cluster[2*nlb];
		list[0] = new Cluster(0,nc);
		assoc = new Ngb[2*nlb][];
		assoc[0] = null;
		
		latest = new int[2*nlb];
		
		Ngb[] ngb = new Ngb[connect];
		
		if (debug) System.out.println("weight computation");

		float[][] zerocov;
		if (covariance==SINGLE) {
			zerocov = new float[1][1];
			zerocov[0][0] = 0.0f;
		} else if (covariance==DIAGONAL) {
			zerocov = new float[1][nc];
			for (int j=0;j<nc;j++) zerocov[0][j] = 0.0f;
		} else {
			zerocov = new float[nc][nc];
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) zerocov[i][j] = 0.0f;
		}
		
		if (covariance==SINGLE) {
			det0 = globalcov[0][0]/basis;
		} else if (covariance==DIAGONAL) {
			det0 = 1.0;
			for (int j=0;j<nc;j++) det0 *= globalcov[0][j]/basis;
		} else {
			Matrix mat = new Matrix(nc,nc);
			for (int n=0;n<nc;n++) for (int m=0;m<nc;m++)
				mat.set(n,m, globalcov[n][m]/basis);
			det0 = mat.det();
		}
		
		if (debug) System.out.println("global determinant: "+det0);

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			
			float[] mean = new float[nc];
			for (int c=0;c<nc;c++) mean[c] = image[c][xyz];
			
			list[labeling[xyz]] = new Cluster(labeling[xyz], mean, zerocov, 1, (float)det0);
		}
		image = null;
		
		if (debug) System.out.println("init neighbors for clustering");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				//if (debug) System.out.print(".");
				int nb=0;
				// 6-C
				if (x>0 && mask[xyz-1])	{
					ngb[nb] = new Ngb(labeling[xyz-1], associationWeight(xyz,xyz-1));
					nb++;
				}
				if (x<nx-1 && mask[xyz+1])	{
					ngb[nb] = new Ngb(labeling[xyz+1], associationWeight(xyz,xyz+1));
					nb++;
				}
				if (y>0 && mask[xyz-nx])	{
					ngb[nb] = new Ngb(labeling[xyz-nx], associationWeight(xyz,xyz-nx));
					nb++;
				}
				if (y<ny-1 && mask[xyz+nx])	{
					ngb[nb] = new Ngb(labeling[xyz+nx], associationWeight(xyz,xyz+nx));
					nb++;
				}
				if (z>0 && mask[xyz-nx*ny])	{
					ngb[nb] = new Ngb(labeling[xyz-nx*ny], associationWeight(xyz,xyz-nx*ny));
					nb++;
				}
				if (z<nz-1 && mask[xyz+nx*ny])	{
					ngb[nb] = new Ngb(labeling[xyz+nx*ny], associationWeight(xyz,xyz+nx*ny));
					nb++;
				}
				// 18-C
				if (connect>6) {
					if (x>0 && y>0 && mask[xyz-1-nx])	{
						ngb[nb] = new Ngb(labeling[xyz-1-nx], associationWeight(xyz,xyz-1-nx));
						nb++;
					}
					if (x<nx-1 && y>0 && mask[xyz+1-nx])	{
						ngb[nb] = new Ngb(labeling[xyz+1-nx], associationWeight(xyz,xyz+1-nx));
						nb++;
					}
					if (x>0 && y<ny-1 && mask[xyz-1+nx])	{
						ngb[nb] = new Ngb(labeling[xyz-1+nx], associationWeight(xyz,xyz-1+nx));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx])	{
						ngb[nb] = new Ngb(labeling[xyz+1+nx], associationWeight(xyz,xyz+1+nx));
						nb++;
					}
					if (y>0 && z>0 && mask[xyz-nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-nx-nx*ny], associationWeight(xyz,xyz-nx-nx*ny));
						nb++;
					}
					if (y<ny-1 && z>0 && mask[xyz+nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+nx-nx*ny], associationWeight(xyz,xyz+nx-nx*ny));
						nb++;
					}
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-nx+nx*ny], associationWeight(xyz,xyz-nx+nx*ny));
						nb++;
					}
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+nx+nx*ny], associationWeight(xyz,xyz+nx+nx*ny));
						nb++;
					}
					if (z>0 && x>0 && mask[xyz-nx*ny-1])	{
						ngb[nb] = new Ngb(labeling[xyz-nx*ny-1], associationWeight(xyz,xyz-nx*ny-1));
						nb++;
					}
					if (z<nz-1 && x>0 && mask[xyz+nx*ny-1])	{
						ngb[nb] = new Ngb(labeling[xyz+nx*ny-1], associationWeight(xyz,xyz+nx*ny-1));
						nb++;
					}
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1])	{
						ngb[nb] = new Ngb(labeling[xyz-nx*ny+1], associationWeight(xyz,xyz-nx*ny+1));
						nb++;
					}
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1])	{
						ngb[nb] = new Ngb(labeling[xyz+nx*ny+1], associationWeight(xyz,xyz+nx*ny+1));
						nb++;
					}
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z>0 && mask[xyz-1-nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-1-nx-nx*ny], associationWeight(xyz,xyz-1-nx-nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z>0 && mask[xyz+1-nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+1-nx-nx*ny], associationWeight(xyz,xyz+1-nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z>0 && mask[xyz-1+nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-1+nx-nx*ny], associationWeight(xyz,xyz-1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-1-nx+nx*ny], associationWeight(xyz,xyz-1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z>0 && mask[xyz+1+nx-nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+1+nx-nx*ny], associationWeight(xyz,xyz+1+nx-nx*ny));
						nb++;
					}
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz-1+nx+nx*ny], associationWeight(xyz,xyz-1+nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+1-nx+nx*ny], associationWeight(xyz,xyz+1-nx+nx*ny));
						nb++;
					}
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny])	{
						ngb[nb] = new Ngb(labeling[xyz+1+nx+nx*ny], associationWeight(xyz,xyz+1+nx+nx*ny));
						nb++;
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
		if (debug) System.out.println("done");
		
		// second pass for delta, cost: not needed
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		
		return;
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final float[][] computeGlobalCovariance() {
		
		if (debug) System.out.println("-- global covariance estimate --");
		
		double[][] cov = new double[nc][nc];
		int nb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				// 6-C
				if (x<nx-1 && mask[xyz+1])	{ dataCovariance(cov,xyz,xyz+1); nb++; }
				if (y<ny-1 && mask[xyz+nx]) { dataCovariance(cov,xyz,xyz+nx); nb++; }
				if (z<nz-1 && mask[xyz+nx*ny])	{dataCovariance(cov,xyz,xyz+nx*ny); nb++; }
				// 18-C
				if (connect>6) {
					if (x>0 && y<ny-1 && mask[xyz-1+nx]) { dataCovariance(cov,xyz,xyz-1+nx); nb++; }
					if (x<nx-1 && y<ny-1 && mask[xyz+1+nx]) { dataCovariance(cov,xyz,xyz+1+nx); nb++; }
					if (y>0 && z<nz-1 && mask[xyz-nx+nx*ny]) { dataCovariance(cov,xyz,xyz-nx+nx*ny); nb++; }
					if (y<ny-1 && z<nz-1 && mask[xyz+nx+nx*ny]) { dataCovariance(cov,xyz,xyz+nx+nx*ny); nb++; }
					if (z>0 && x<nx-1 && mask[xyz-nx*ny+1]) { dataCovariance(cov,xyz,xyz-nx*ny+1); nb++; }
					if (z<nz-1 && x<nx-1 && mask[xyz+nx*ny+1]) { dataCovariance(cov,xyz,xyz+nx*ny+1); nb++; }
				}
				// 26-C
				if (connect>18) {
					if (x>0 && y>0 && z<nz-1 && mask[xyz-1-nx+nx*ny]) { dataCovariance(cov,xyz,xyz-1-nx+nx*ny); nb++; }
					if (x>0 && y<ny-1 && z<nz-1 && mask[xyz-1+nx+nx*ny]) { dataCovariance(cov,xyz,xyz-1+nx+nx*ny); nb++; }
					if (x<nx-1 && y>0 && z<nz-1 && mask[xyz+1-nx+nx*ny]) { dataCovariance(cov,xyz,xyz+1-nx+nx*ny); nb++; }
					if (x<nx-1 && y<ny-1 && z<nz-1 && mask[xyz+1+nx+nx*ny]) { dataCovariance(cov,xyz,xyz+1+nx+nx*ny); nb++; }
				}
			}
		}
		
		// PCA analysis?
		Matrix covMat = new Matrix(cov);
		SingularValueDecomposition svd = covMat.svd();
		if (debug) System.out.println("matrix svd: \n"+displayMatrix(svd.getSingularValues()));
		double[] eig = svd.getSingularValues();
		double sum = 0.0;
		for (int n=0;n<nc;n++) sum += eig[n];
		int nbest = -1;
		double first = 0.0;
		for (int n=0;n<nc && nbest==-1;n++) {
			first += eig[n];
			if (first >= eigratio*sum) nbest = n;
		}
		if (debug) System.out.println("pca dimension ("+eigratio+" of energy): "+nbest);
		
		if (covariance==PCA) {
			// project the data into PCA space
			double[] datapca = new double[nbest];
			double[][] basis = svd.getU().getArray();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					for (int n=0;n<nbest;n++) {
						datapca[n] = 0.0;
						for (int m=0;m<nc;m++) datapca[n] += basis[m][n]*image[m][xyz];
					}
					// replace directly the data once all coeffs are computed
					for (int n=0;n<nbest;n++) {
						image[n][xyz] = (float)datapca[n];
					}
				}
			}
			
			// trim the unused elements
			float[][] tmp = new float[nbest][];
			for (int n=0;n<nbest;n++) tmp[n] = image[n];
			image = tmp;
			
			// change the covariance matrix
			double[][] tmpcov = new double[nbest][nbest];
			for (int n=0;n<nbest;n++) {
				for (int m=0;m<nbest;m++) tmpcov[n][m] = 0.0;
				tmpcov[n][n] = eig[n];
			}
			cov = tmpcov;
			// change the data dimensions
			nc = nbest;
		}	

		float[][] meancov;
		if (covariance==SINGLE) {
			meancov = new float[1][1];
			meancov[0][0] = 0.0f;
			for (int j=0;j<nc;j++) {
				meancov[0][0] += (float)(0.5*cov[0][j]/nb/nc);	
			}
		} else if (covariance==DIAGONAL) {
			meancov = new float[1][nc];
			for (int j=0;j<nc;j++) {
				meancov[0][j] = (float)(0.5*cov[0][j]/nb);	
			}
		} else {
			meancov = new float[nc][nc];
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				meancov[i][j] = (float)(0.5*cov[i][j]/nb);	
				// for svd below
				cov[i][j] = 0.5*cov[i][j]/nb;	
			}
		}
		
		if (debug) System.out.println("Covariance matrix: \n"+displayMatrix(meancov));
		
		
		// multiply by basis to act like a sum of differences
		if (covariance==SINGLE) {
			meancov[0][0] *= basis;	
		} else if (covariance==DIAGONAL) {
			for (int j=0;j<nc;j++) {
				meancov[0][j] *= basis;	
			}
		} else {
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				meancov[i][j] *= basis;	
			}
		}
		return meancov;
	}

	private final void dataCovariance(double[][] cov, int id1, int id2) {
		if (covariance==SINGLE || covariance==DIAGONAL) {
			for (int t=0;t<nc;t++) {
				cov[0][t] += (image[t][id1]-image[t][id2])*(image[t][id1]-image[t][id2]);
			}
		} else {
			for (int t1=0;t1<nc;t1++) for (int t2=0;t2<nc;t2++) {
				cov[t1][t2] += (image[t1][id1]-image[t1][id2])*(image[t2][id1]-image[t2][id2]);
			}	
		}
		return;
	}

	public final void hierarchicalClustering(int k0, boolean firststop, int maxlength) {
		
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
					
					// new values
					Ngb[] aNodeN = new Ngb[bNodeN.length+pNodeN.length-2];
					
					// sum(uv) = sum(u) + sum(v)
					// sq2(uv) = sq2(u) + sq2(v) + 1/n(uv) [n(u)/n(v)*m(v) - n(v)/n(u)*m(u)]^2
					// n(uv) = n(u) + n(v)
					
					float[] mean = new float[nc];
					for (int n=0;n<nc;n++) mean[n] = bNodeC.mean[n]+pNodeC.mean[n];
					float size = bNodeC.size+pNodeC.size;
					float[][] cov;
					double det;
					if (covariance==SINGLE) {
						cov = new float[1][1];
						float factor = bNodeC.size*pNodeC.size/size;
						cov[0][0] = bNodeC.cov[0][0] + pNodeC.cov[0][0] 
									+ factor*(bNodeC.mean[0]/bNodeC.size-pNodeC.mean[0]/pNodeC.size)
											*(bNodeC.mean[0]/bNodeC.size-pNodeC.mean[0]/pNodeC.size);
						// determinant : need to mix with global cov					
						det = (cov[0][0]+globalcov[0][0])/(size+basis-1.0);					
					} else if (covariance==DIAGONAL) {
						cov = new float[1][nc];
						float factor = bNodeC.size*pNodeC.size/size;
						for (int m=0;m<nc;m++) {
							cov[0][m] = bNodeC.cov[0][m] + pNodeC.cov[0][m] 
									+ factor*(bNodeC.mean[m]/bNodeC.size-pNodeC.mean[m]/pNodeC.size)
											*(bNodeC.mean[m]/bNodeC.size-pNodeC.mean[m]/pNodeC.size);
						}
						// determinant : need to mix with global cov	
						det = 1.0;
						for (int m=0;m<nc;m++)	det *= (cov[0][m]+globalcov[0][m])/(size+basis-1.0);					
					} else {
						cov = new float[nc][nc];
						float factor = bNodeC.size*pNodeC.size/size;
						for (int n=0;n<nc;n++) for (int m=0;m<nc;m++) {
							cov[n][m] = bNodeC.cov[n][m] + pNodeC.cov[n][m] 
									+ factor*(bNodeC.mean[n]/bNodeC.size-pNodeC.mean[n]/pNodeC.size)
											*(bNodeC.mean[m]/bNodeC.size-pNodeC.mean[m]/pNodeC.size);
						}
						// determinant : need to mix with global cov					
						Matrix mat = new Matrix(nc, nc);
						for (int n=0;n<nc;n++) for (int m=0;m<nc;m++)
							mat.set(n,m, (cov[n][m]+globalcov[n][m])/(size+basis-1.0) );
						det = mat.det();
					}
					Cluster aNodeC = new Cluster(id, mean,  cov, size, (float)det);
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNodeN.length+pNodeN.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
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
					//System.out.print(">");
					
					/*
					// first step to get rid of bad candidates (to avoid computing too many inversions)
					BitSet discardedNgb = new BitSet(aNodeN.length);
					if ( (mode==JENSENSHANNON) && (!simplified)) {
						float best = 0.0f;
						for (int n=0; n<aNodeN.length; n++) {
							Cluster wngbC = list[aNodeN[n].id];
							float metric = (float)approxJsDiv(aNodeC.mean, wngbC.mean, 
																	aNodeC.cov, wngbC.cov, 
																	aNodeC.size, wngbC.size);
							aNodeN[n].delta = metric;
							if (metric>best) best = metric;
						}
						for (int n=0; n<aNodeN.length; n++) {
							if (aNodeN[n].delta<0.5f*best) {
								discardedNgb.set(n);
								aNodeN[n].delta = -1e6f;
							}
						}
					}
						
					for (int n=0; n<aNodeN.length; n++) if (!discardedNgb.get(n)) {
					*/
					for (int n=0; n<aNodeN.length; n++) {
						Cluster wngbC = list[aNodeN[n].id];
						
						// build the divergence metric
						float metric = (float)distributionMetric(aNodeC.mean, wngbC.mean, 
																	aNodeC.cov, wngbC.cov, 
																	aNodeC.size, wngbC.size,
																	aNodeC.det, wngbC.det);
						
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNodeC.size+wngbC.size)/((connect+1.0f)*nlb)) );
						
						
						float threshold = distributionThreshold(aNodeC.size, wngbC.size);
						
						
						aNodeN[n].delta = wsize*(metric-threshold);
						
						//System.out.print("m="+metric+"; ");
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						float avgcov = 0.0f;
						if (covariance==SINGLE) avgcov += nc*aNodeC.cov[0][0];
						else if (covariance==DIAGONAL) for (int c=0;c<nc;c++) avgcov += aNodeC.cov[0][c];
						else for (int c=0;c<nc;c++) avgcov += aNodeC.cov[c][c];
						avgcov = (float)FastMath.sqrt(avgcov/(nc*aNodeC.size));
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "
																	+"cov= "+avgcov
																	+", n= "+aNodeC.size
																	+" ("+maxtree.getCurrentSize()+"|"
																	+bNodeN.length+", "+pNodeN.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
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
	private final float associationWeight(int id1, int id2) {
		return (float)(distributionMetric(list[labeling[id1]].mean,
							 				list[labeling[id2]].mean,
							 					globalcov,globalcov,1,1,(float)det0,(float)det0)-pvalue);
	}
	
	private final float distributionThreshold(double n1, double n2) {
		// Sidak correction
		if (mode==HOTELLING) return 1.0f - (float)FastMath.pow(1.0f-pvalue, 1.0f/(n1+n2+basis-1.0f));
		// JS div threshold??				
		else if (mode==JENSENSHANNON) return pvalue;
		return 0.5f;
	}

	private final float distributionMetric(final float[] m1, final float[] m2,
											  final float[][] s1, final float[][] s2,
											  final float n1, final float n2,
											  final float det1, final float det2) {
		// T-test
		if (mode==HOTELLING) return hotellingTest(m1, m2, s1, s2, n1, n2);
		// JS div 		
		else if (mode==JENSENSHANNON) return jsDiv(m1, m2, s1, s2, n1, n2, det1, det2);
		return 0.5f;
	}

	/**
     * Computes p-value for 2-sided, 2-sample t-test.
     * <p>
     * Does not assume subpopulation variances are equal. Degrees of freedom
     * are estimated from the data. Adapted from commons-math at apache.org</p>
     *
     * @param m1 first sample mean
     * @param m2 second sample mean
     * @param s1 first sample variance x size
     * @param s2 second sample variance x size
     * @param n1 first sample n
     * @param n2 second sample n
     * @return p-value
     * @throws MaxCountExceededException if an error occurs computing the p-value
     */
    private float hotellingTest(final float[] m1, final float[] m2,
                          final float[][] s1, final float[][] s2,
                          final float n1, final float n2) {

		double t2 = 0.0;
		if (covariance==SINGLE) {
			double v12 = (s1[0][0] + s2[0][0] + globalcov[0][0])/(n1+n2+basis-1.0);
			
			for (int j=0;j<nc;j++) {
				t2 += (m1[j]/n1-m2[j]/n2)*(m1[j]/n1-m2[j]/n2)/v12;
			}
			
			t2 *= n1*n2/(n1+n2)*(n1+n2+basis-nc-1.0)/( (n1+n2+basis-2.0)*nc );
    	} else if (covariance==DIAGONAL) {
			double[] v12 = new double[nc];
			for (int j=0;j<nc;j++) {
				v12[j] = (s1[0][j] + s2[0][j] + globalcov[0][j])/(n1+n2+basis-1.0);
			}
			
			for (int j=0;j<nc;j++) {
				t2 += (m1[j]/n1-m2[j]/n2)*(m1[j]/n1-m2[j]/n2)/v12[j];
			}
			
			t2 *= n1*n2/(n1+n2)*(n1+n2+basis-nc-1.0)/( (n1+n2+basis-2.0)*nc );
    	} else {
			double[][] v12 = new double[nc][nc];
			double[][] v12inv = new double[nc][nc];
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				v12[i][j] = (s1[i][j] + s2[i][j] + globalcov[i][j])/(n1+n2+basis-1.0);
			}
			LUDecomposition lu = new LUDecomposition(new Matrix(v12));
			v12inv = lu.solve(Matrix.identity(nc,nc)).getArray();
			
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				t2 += (m1[i]/n1-m2[i]/n2)*v12inv[i][j]*(m1[j]/n1-m2[j]/n2);
			}
			
			t2 *= n1*n2/(n1+n2)*(n1+n2+basis-nc-1.0)/( (n1+n2+basis-2.0)*nc );
    	}
    	FDistribution distribution = new FDistribution(nc, n1+n2+basis-1.0-nc);
        return (float)(1.0-distribution.cumulativeProbability(t2));
    }

    private float jsDiv(final float[] m1, final float[] m2,
                          		final float[][] s1, final float[][] s2,
								  final float n1, final float n2,
								  final float det1, final float det2) {
       
		double js = 0.0;
		if (covariance==SINGLE) {
			double v12 = 0.25*( (s1[0][0] + globalcov[0][0])/(n1+basis-1.0)
								+(s2[0][0] + globalcov[0][0])/(n2+basis-1.0) );
			double det = v12;
			
			for (int j=0;j<nc;j++) {
				js += 0.25*(m1[j]/n1-m2[j]/n2)*(m1[j]/n1-m2[j]/n2)/v12;
			}
			js -= 0.5*FastMath.log(det1/det*det2/det);
			// adjust to zero
			js += nc*FastMath.log(2.0);
		} else if (covariance==DIAGONAL) {
			double[] v12 = new double[nc];
			double det = 1.0;
			for (int j=0;j<nc;j++) {
				v12[j] = 0.25*( (s1[0][j] + globalcov[0][j])/(n1+basis-1.0)
								+(s2[0][j] + globalcov[0][j])/(n2+basis-1.0) );
				det *= v12[j];
			}
			for (int j=0;j<nc;j++) {
				js += 0.25*(m1[j]/n1-m2[j]/n2)*(m1[j]/n1-m2[j]/n2)/v12[j];
			}
			js -= 0.5*FastMath.log(det1/det*det2/det);
			// adjust to zero
			js += nc*FastMath.log(2.0);
		} else {
			double[][] v12 = new double[nc][nc];
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				v12[i][j] = 0.25*( (s1[i][j] + globalcov[i][j])/(n1+basis-1.0)
									+(s2[i][j] + globalcov[i][j])/(n2+basis-1.0) );
			}
			LUDecomposition lu = new LUDecomposition(new Matrix(v12));
			double det = lu.det();
			v12 = lu.solve(Matrix.identity(nc,nc)).getArray();
			
			for (int i=0;i<nc;i++) for (int j=0;j<nc;j++) {
				js += 0.25*(m1[i]/n1-m2[i]/n2)*v12[i][j]*(m1[j]/n1-m2[j]/n2);
			}
			js -= 0.5*FastMath.log(det1/det*det2/det);
			// adjust to zero
			js += nc*FastMath.log(2.0);
		}
    	// Gaussian model on the distances
    	return (float)FastMath.exp(-0.5*js/nc);
    }

    private float approxJsDiv(final float[] m1, final float[] m2,
                          		final float[][] s1, final float[][] s2,
								  final float n1, final float n2) {
       
		double js = 0.0;
		double[] v1 = new double[nc];
		double[] v2 = new double[nc];
		double[] v12 = new double[nc];
		double det = 1.0;
		double det1 = 1.0, det2 = 1.0;
		for (int j=0;j<nc;j++) {
			v1[j] = (s1[0][j] + globalcov[0][j])/(n1+basis-1.0);
			v2[j] = (s2[0][j] + globalcov[0][j])/(n2+basis-1.0);
			v12[j] = 0.25*(v1[j] + v2[j]);
			det1 *= v1[j];
			det2 *= v2[j];
			det *= v12[j];
		}
		for (int j=0;j<nc;j++) {
			js += 0.25*(m1[j]/n1-m2[j]/n2)*(m1[j]/n1-m2[j]/n2)/v12[j];
		}
		js -= 0.5*FastMath.log(det1/det*det2/det);
		// adjust to zero
		js += nc*FastMath.log(2.0);
    	// Gaussian model on the distances
    	return (float)FastMath.exp(-0.5*js/nc);
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
	public final float[][][][] exportClusteringMean(int maxlb) {
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
		float[][][][] tmp = new float[nx][ny][nz][nc];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (mask[xyz]) {
				for (int c=0;c<nc;c++)
					tmp[x][y][z][c] = list[latest[labeling[xyz]]].mean[c]/list[latest[labeling[xyz]]].size;
			}
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][][] exportClusteringVariance(int maxlb) {
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
		float[][][][] tmp;
		if (covariance==SINGLE) {
			tmp = new float[nx][ny][nz][1];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					tmp[x][y][z][0] = (float)FastMath.sqrt(list[latest[labeling[xyz]]].cov[0][0]/Numerics.max(1.0f,(list[latest[labeling[xyz]]].size-1.0f)));
				}
			}
		} else if (covariance==DIAGONAL) {
			tmp = new float[nx][ny][nz][nc];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					for (int c=0;c<nc;c++)
						tmp[x][y][z][c] = (float)FastMath.sqrt(list[latest[labeling[xyz]]].cov[0][c]/Numerics.max(1.0f,(list[latest[labeling[xyz]]].size-1.0f)));
				}
			}
		} else {
			tmp = new float[nx][ny][nz][nc];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (mask[xyz]) {
					for (int c=0;c<nc;c++)
						tmp[x][y][z][c] = (float)FastMath.sqrt(list[latest[labeling[xyz]]].cov[c][c]/Numerics.max(1.0f,(list[latest[labeling[xyz]]].size-1.0f)));
				}
			}
		}
		return tmp;
	}
	

	private final String displayMatrix (double[][] mat) {
		String output = "";
		for (int i=0;i<mat.length;i++) {
			output += "(	";
			for (int j=0;j<mat[0].length;j++) {
				output += mat[i][j]+"	";
			}
			output += ")\n";
		}
		return output;
	}

	private final String displayMatrix (float[][] mat) {
		String output = "";
		for (int i=0;i<mat.length;i++) {
			output += "(	";
			for (int j=0;j<mat[0].length;j++) {
				output += mat[i][j]+"	";
			}
			output += ")\n";
		}
		return output;
	}
	
	private final String displayMatrix (double[] mat) {
		String output = "";
		output += "(	";
		for (int j=0;j<mat.length;j++) {
			output += mat[j]+"	";
		}
		output += ")\n";

		return output;
	}


}
