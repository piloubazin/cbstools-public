package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.analysis.function.Gaussian;

/**
 *
 *  This class handles octree-style grouping for JSDiv multi-contrast clustering
 *
 *	@version    April 2015
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class OctreeMultiClusterSimplification {

	// store the multiscale data
	private float[][][] sum;
	private float[][][] sqd;
	private BitSet[] valid;
	private int[][] npt;
	private int[][] lbl;
	
	private int[] nx;
	private int[] ny;
	private int[] nz;
	private int   nt;
	private int   nlb;
	
	private	static 	int 		LEVELS;
	public static int[]		CS;
	public static int[]		CL;
	
	private float[] sigma2;
	private float threshold;
	
	private static final double LN2 = FastMath.log(2.0);
	private static final double JS2 = 0.5*FastMath.pow(2.0,0.25);
	private static final double JS8 = 0.5*FastMath.pow(8.0,0.25);
	private static final double SQLN2 = FastMath.sqrt(FastMath.log(2.0));
	private static final double SQLN8 = FastMath.sqrt(FastMath.log(8.0));
	
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
	
	private int connect = 6;
	private int[][] ngb;
	private float[] ngbfactor;
	
	private static final boolean	debug = false;
	private static final boolean	verbose = true;	
	
	public OctreeMultiClusterSimplification(float[][] image, boolean[] mask, float[][] imgdev, int nix, int niy, int niz, int nit, float[] stdev, float devfactor, float thres, String metrictype, String varup, int clustsize, int conn_) {
		
		nt = nit;
		
		sigma2 = new float[nt];
		for (int t=0;t<nt;t++) sigma2[t] = devfactor*devfactor*stdev[t]*stdev[t];
		
		threshold  = thres;
		
		if (metrictype.equals("joint Jensen-Shannon")) metric = JSDIV_SUM;
		else if (metrictype.equals("max Jensen-Shannon")) metric = JSDIV_MAX;
		else if (metrictype.equals("exact Jensen-Shannon")) metric = JSDIV_EXACT;
		else if (metrictype.equals("table Jensen-Shannon")) metric = JSDIV_TABLE;
		
		if (metric==JSDIV_TABLE) jsdtable = new JensenShannonDivTable(6.0, 0.01, 0.01);
		
		if (varup.equals("Maximum")) varupdate = VARMAX;
		else if (varup.equals("Linear")) varupdate = VARLIN;
		else if (varup.equals("Zero")) varupdate = VARLIN;
		else varupdate = VARCONST;
		
		clustersize = clustsize;
		
		connect = conn_;
		
		LEVELS = findMinimumLevel(nix, niy, niz);
		
		CS = new int[LEVELS];
		int size = 1;
		CL = new int[LEVELS];
		int length = 1;
		for (int l=0;l<LEVELS;l++) {
			CS[l] = size;
			size *= 8;
			CL[l] = length;
			length *= 2;
		}
		
		// start the arrays
		nx = new int[LEVELS];
		ny = new int[LEVELS];
		nz = new int[LEVELS];
		
		sum = new float[LEVELS][][];
		sqd = new float[LEVELS][][];
		valid = new BitSet[LEVELS];

		npt = new int[LEVELS][];
		lbl = new int[LEVELS][];
		
		// fill them in
		nx[0] = nix;
		ny[0] = niy;
		nz[0] = niz;
		for (int l=1;l<LEVELS;l++) {
			nx[l] = Numerics.floor(0.5f*nx[l-1]);
			ny[l] = Numerics.floor(0.5f*ny[l-1]);
			nz[l] = Numerics.floor(0.5f*nz[l-1]);
		}
		
		ngb = new int[connect][LEVELS];
		for (int l=0;l<LEVELS;l++) {
			ngb[0][l] = 1;
			ngb[1][l] = -1;
			ngb[2][l] = nx[l];
			ngb[3][l] = -nx[l];
			ngb[4][l] = nx[l]*ny[l];
			ngb[5][l] = -nx[l]*ny[l];
			if (connect>6) {
				ngb[6][l] =  1+nx[l];
				ngb[7][l] = -1+nx[l];
				ngb[8][l] =  1-nx[l];
				ngb[9][l] = -1-nx[l];
				ngb[10][l] =  nx[l]+nx[l]*ny[l];
				ngb[11][l] = -nx[l]+nx[l]*ny[l];
				ngb[12][l] =  nx[l]-nx[l]*ny[l];
				ngb[13][l] = -nx[l]-nx[l]*ny[l];
				ngb[14][l] =  nx[l]*ny[l]+1;
				ngb[15][l] = -nx[l]*ny[l]+1;
				ngb[16][l] =  nx[l]*ny[l]-1;
				ngb[17][l] = -nx[l]*ny[l]-1;
			}
			if (connect>18) {
				ngb[18][l] =  1+nx[l]+nx[l]*ny[l];
				ngb[19][l] = -1+nx[l]+nx[l]*ny[l];
				ngb[20][l] =  1-nx[l]+nx[l]*ny[l];
				ngb[21][l] = -1-nx[l]+nx[l]*ny[l];
				ngb[22][l] =  1+nx[l]-nx[l]*ny[l];
				ngb[23][l] = -1+nx[l]-nx[l]*ny[l];
				ngb[24][l] =  1-nx[l]-nx[l]*ny[l];
				ngb[25][l] = -1-nx[l]-nx[l]*ny[l];
			}
		}
		ngbfactor = new float[connect];
		for (int c=0;c<connect;c++) {
			if (c<6) ngbfactor[c] = 1.0f;
			else if (c<18) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(2.0f);
			else if (c<26) ngbfactor[c] = 1.0f/(float)FastMath.sqrt(3.0f);
		}
		// work on copy
		sum[0] = new float[nx[0]*ny[0]*nz[0]][nt];;
		sqd[0] = new float[nx[0]*ny[0]*nz[0]][nt];
		npt[0] = new int[nx[0]*ny[0]*nz[0]];
		valid[0] = new BitSet(nx[0]*ny[0]*nz[0]);
		for (int x=0;x<nx[0];x++) for (int y=0;y<ny[0];y++) for (int z=0;z<nz[0];z++) {
			int xyz = x+nx[0]*y+nx[0]*ny[0]*z;
			if (mask[xyz]) {
				for (int t=0;t<nt;t++) sum[0][xyz][t] = image[xyz][t];
				valid[0].set(xyz, true);
			} else {
				for (int t=0;t<nt;t++) sum[0][xyz][t] = 0.0f;
				valid[0].set(xyz, false);
			}
			if (imgdev==null) for (int t=0;t<nt;t++) sqd[0][xyz][t] = 0.0f;
			else for (int t=0;t<nt;t++) sqd[0][xyz][t] = Numerics.square(devfactor*imgdev[xyz][t]);
			npt[0][xyz] = 1;
		}
	}
	
	/** compute the necessary octree level for storing an image of nx,ny,nz dimensions */
	public static int findMinimumLevel(int nix, int niy, int niz) {
		int level = 1;
		int size = 1;
		while ( (size<nix) || (size<niy) || (size<niz) ) {
			size = 2*size;
			level++;
		}
		return level;
	}
	
	/** group similar values in octree fashion
	 *	based on Jensen-Shannon divergence
	 */
	public int mergeToJensenShannonDistance(int toplevel) {
		boolean stop=false;
		int maxlevel = 0;
		for (int l=1;l<LEVELS && l<toplevel && !stop;l++) {
			if (verbose) System.out.print(" level "+l);
			
			// create the arrays here (so it doesn't try more than necessary)
			sum[l] = new float[nx[l]*ny[l]*nz[l]][nt];
			sqd[l] = new float[nx[l]*ny[l]*nz[l]][nt];
			npt[l] = new int[nx[l]*ny[l]*nz[l]];
			valid[l] = new BitSet(nx[l]*ny[l]*nz[l]);
			
			// must set nlb to run the metric
			nlb = valid[l-1].cardinality();
		
			stop = true;
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;

				// previous scale
				int xyzp = 2*x + nx[l-1]*2*y + nx[l-1]*ny[l-1]*2*z;
				int xyzpx = xyzp+1;
				int xyzpy = xyzp+nx[l-1];
				int xyzpz = xyzp+nx[l-1]*ny[l-1];
				int xyzpxy = xyzp+1+nx[l-1];
				int xyzpyz = xyzp+nx[l-1]+nx[l-1]*ny[l-1];
				int xyzpzx = xyzp+nx[l-1]*ny[l-1]+1;
				int xyzpxyz = xyzp+1+nx[l-1]+nx[l-1]*ny[l-1];
				
				// check if all valid at lower level
				boolean allvalid = true;
				if (!valid[l-1].get(xyzp) 
				 || !valid[l-1].get(xyzpx)
				 || !valid[l-1].get(xyzpy)
				 || !valid[l-1].get(xyzpz)
				 || !valid[l-1].get(xyzpxy)
				 || !valid[l-1].get(xyzpyz)
				 || !valid[l-1].get(xyzpzx)
				 || !valid[l-1].get(xyzpxyz) ) allvalid = false;
				
				if (allvalid) {
					//if (verbose) System.out.print(">");
			
					// compute the joint JS divergence : for each connected leaf pair!
					// assume 6C to start: 12 pairs
					boolean merge = true;
					//if (verbose) System.out.print("+");
					//System.out.flush();
					
					// groupwise JS div
					float js = octreeJSdiv(xyzp,xyzpx,xyzpy,xyzpz,xyzpxy,xyzpyz,xyzpzx,xyzpxyz, l-1);
					
					if (js<0) merge = false;
					//else if (verbose) System.out.print("js: "+js+" ");
						
					// merging: add the sum, sqd incrementally for each leaf with Chan's formula
					if (merge) {
						// add progressively all children
						float cs = CS[l-1];
						
						for (int t=0;t<nt;t++) {
							sum[l][xyz][t] = sum[l-1][xyzp][t];
							sqd[l][xyz][t] = sqd[l-1][xyzp][t];
							
							// update first the square distance to keep the valid sum and the added sum separate
							sqd[l][xyz][t] += sqd[l-1][xyzpx][t] + 1*cs*cs/(1*cs+cs)*Numerics.square(sum[l][xyz][t]/(1*cs)-sum[l-1][xyzpx][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpx][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpy][t] + 2*cs*cs/(2*cs+cs)*Numerics.square(sum[l][xyz][t]/(2*cs)-sum[l-1][xyzpy][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpy][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpz][t] + 3*cs*cs/(3*cs+cs)*Numerics.square(sum[l][xyz][t]/(3*cs)-sum[l-1][xyzpz][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpz][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpxy][t] + 4*cs*cs/(4*cs+cs)*Numerics.square(sum[l][xyz][t]/(4*cs)-sum[l-1][xyzpxy][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpxy][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpyz][t] + 5*cs*cs/(5*cs+cs)*Numerics.square(sum[l][xyz][t]/(5*cs)-sum[l-1][xyzpyz][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpyz][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpzx][t] + 6*cs*cs/(6*cs+cs)*Numerics.square(sum[l][xyz][t]/(6*cs)-sum[l-1][xyzpzx][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpzx][t];
							
							sqd[l][xyz][t] += sqd[l-1][xyzpxyz][t] + 7*cs*cs/(7*cs+cs)*Numerics.square(sum[l][xyz][t]/(7*cs)-sum[l-1][xyzpxyz][t]/cs);
							sum[l][xyz][t] += sum[l-1][xyzpxyz][t];
						}						
						valid[l].set(xyz);
						npt[l][xyz] = CS[l];
						stop = false;
						maxlevel = l;
					}
				}
			}
		}
		if (maxlevel==0) {
			int l=0;
			sum[l] = new float[nx[l]*ny[l]*nz[l]][nt];
			sqd[l] = new float[nx[l]*ny[l]*nz[l]][nt];
			npt[l] = new int[nx[l]*ny[l]*nz[l]];
			valid[l] = new BitSet(nx[l]*ny[l]*nz[l]);
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				valid[l].set(xyz);
				npt[l][xyz] = CS[l];
			}
		}
		return maxlevel;
	}
	
	
	/** group similar values in octree fashion
	 *	based on Jensen-Shannon divergence
	 */
	public int competitionWithJensenShannonDistance(int l, int nclusters, float approx) {
		// generate the merged labeling at level l
		nlb = generateNextLabelsAt(l, nclusters);
		if (verbose) System.out.println("region competition: init clusters"+nlb);
			
		// keep the cluster level data in separate (duplicate) arrays
		float[][] listsum = new float[nclusters+1][nt];
		float[][] listsqd = new float[nclusters+1][nt];
		int[] listnpt = new int[nclusters+1];
		
		// put into the tree any regions at level l with a neighbor at lvl+1
		BinaryHeapPair tree = new BinaryHeapPair(Numerics.max((nlb-nclusters)/16,20), Numerics.max((nlb-nclusters)/16, 20), BinaryHeapPair.MAXTREE);
		for (int x=1;x<nx[l]-1;x++) for (int y=1;y<ny[l]-1;y++) for (int z=1;z<nz[l]-1;z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>0 && lbl[l][xyz]<=nclusters) {
				for (int c=0;c<connect;c++) {
					int xyzn = xyz+ngb[c][l];
					
					if (lbl[l][xyzn]>nclusters) {
						float js = ngbfactor[c]*regionCompetitionJSdiv(xyz, xyzn, l);
						if (js>0) tree.addValue(js, xyzn, xyz);
					}
				}
				/*
				// look for neighbors from new scale
				int xyzpx = xyz+1;
				int xyzpy = xyz+nx[l];
				int xyzpz = xyz+nx[l]*ny[l];
				int xyzmx = xyz-1;
				int xyzmy = xyz-nx[l];
				int xyzmz = xyz-nx[l]*ny[l];
				
				if (lbl[l][xyzpx]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzpx, l);
					if (js>0) tree.addValue(js, xyzpx, xyz);
				}
				if (lbl[l][xyzmx]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzmx, l);
					if (js>0) tree.addValue(js, xyzmx, xyz);
				}
				if (lbl[l][xyzpy]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzpy, l);
					if (js>0) tree.addValue(js, xyzpy, xyz);
				}
				if (lbl[l][xyzmy]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzmy, l);
					if (js>0) tree.addValue(js, xyzmy, xyz);
				}
				if (lbl[l][xyzpz]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzpz, l);
					if (js>0) tree.addValue(js, xyzpz, xyz);
				}
				if (lbl[l][xyzmz]>nclusters) {
					float js = regionCompetitionJSdiv(xyz, xyzmz, l);
					if (js>0) tree.addValue(js, xyzmz, xyz);
				}
				*/
				// build the global label list
				for (int t=0;t<nt;t++) {
					listsum[lbl[l][xyz]][t] = sum[l][xyz][t];
					listsqd[lbl[l][xyz]][t] = sqd[l][xyz][t];
				}
				listnpt[lbl[l][xyz]] = npt[l][xyz];
			}
		}
		
		while (tree.isNotEmpty()) {
			float score = tree.getFirst();
			int xyz = tree.getFirstId1();
			int xyzngb = tree.getFirstId2();
			tree.removeFirst();
			
			if (lbl[l][xyz]>nclusters) {
				
				// update neighbor for the next js computation
				int lbngb = lbl[l][xyzngb];
				for (int t=0;t<nt;t++) {
					sqd[l][xyzngb][t] = listsqd[lbngb][t];
					sum[l][xyzngb][t] = listsum[lbngb][t];
				}
				npt[l][xyzngb] = listnpt[lbngb];
					
				// recompute the score: it may have changed
				float jsupdate = regionCompetitionJSdiv(xyz, xyzngb, l);
				if (jsupdate>=approx*score) {
				
					// merge with label
					int lb = lbl[l][xyzngb];
					lbl[l][xyz] = lb;
					
					// update the global lists
					for (int t=0;t<nt;t++) {
						listsqd[lb][t] += sqd[l][xyz][t] + listnpt[lb]*npt[l][xyz]/(listnpt[lb]+npt[l][xyz])*Numerics.square(sum[l][xyz][t]/npt[l][xyz]-listsum[lb][t]/listnpt[lb]);
						listsum[lb][t] += sum[l][xyz][t];
					}
					listnpt[lb] += npt[l][xyz];
					
					// update also current point for the next js computation
					for (int t=0;t<nt;t++) {
						sqd[l][xyz][t] = listsqd[lb][t];
						sum[l][xyz][t] = listsum[lb][t];
					}
					npt[l][xyz] = listnpt[lb];
					
					// find neighboring labels not yet relabeled
					
					// check for neighbor boundaries
					int z = Numerics.floor(xyz/(nx[l]*ny[l]));
					int y = Numerics.floor((xyz-nx[l]*ny[l]*z)/nx[l]);
					int x = xyz-nx[l]*ny[l]*z-nx[l]*y;
	
					int xyzpx = xyz+1;
					if (x+1<nx[l] && lbl[l][xyzpx]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzpx, l);
						if (js>0) tree.addValue(js, xyzpx, xyz);
					}
					int xyzmx = xyz-1;
					if (x-1>=0 && lbl[l][xyzmx]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzmx, l);
						if (js>0) tree.addValue(js, xyzmx, xyz);
					}
					int xyzpy = xyz+nx[l];
					if (y+1<ny[l] && lbl[l][xyzpy]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzpy, l);
						if (js>0) tree.addValue(js, xyzpy, xyz);
					}
					int xyzmy = xyz-nx[l];
					if (y-1>=0 && lbl[l][xyzmy]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzmy, l);
						if (js>0) tree.addValue(js, xyzmy, xyz);
					}
					int xyzpz = xyz+nx[l]*ny[l];
					if (z+1<nz[l] && lbl[l][xyzpz]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzpz, l);
						if (js>0) tree.addValue(js, xyzpz, xyz);
					}
					int xyzmz = xyz-nx[l]*ny[l];
					if (z-1>=0 && lbl[l][xyzmz]>nclusters) {
						float js = regionCompetitionJSdiv(xyz, xyzmz, l);
						if (js>0) tree.addValue(js, xyzmz, xyz);
					}
				} else if (jsupdate>0) {
					tree.addValue(jsupdate, xyz, xyzngb);
				}
			}
		}
		tree.finalize();
		tree = null;
		
		if (verbose) System.out.println("remap labels");
		
		// update data
		// and remap the leftover labels
		int maxlb = nclusters;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>nclusters) {
				maxlb++;
				// remap non-visited clusters
				lbl[l][xyz] = maxlb;
			} else if (lbl[l][xyz]>0) {
				// update the cluster values
				int lb = lbl[l][xyz];
				for (int t=0;t<nt;t++) {
					sum[l][xyz][t] = listsum[lb][t];
					sqd[l][xyz][t] = listsqd[lb][t];
				}
				npt[l][xyz] = listnpt[lb];
			}
		}
		if (verbose) System.out.println("done.\n final clusters"+maxlb);
		
		return maxlb;
	}

	float octreeJSdiv(int xyz1, int xyz2, int xyz3, int xyz4, int xyz5, int xyz6, int xyz7, int xyz8, int level) {
		// joint variance: sum or max over intensities?
		double js = 0.0;
		for (int t=0;t<nt;t++) {
			double v1=0,v2=0,v3=0,v4=0,v5=0,v6=0,v7=0,v8=0;
			if (varupdate==VARMAX) {
				// strong prior influence
				v1 = Numerics.max(sqd[level][xyz1][t]/npt[level][xyz1], sigma2[t]);
				v2 = Numerics.max(sqd[level][xyz2][t]/npt[level][xyz2], sigma2[t]);
				v3 = Numerics.max(sqd[level][xyz3][t]/npt[level][xyz3], sigma2[t]);
				v4 = Numerics.max(sqd[level][xyz4][t]/npt[level][xyz4], sigma2[t]);
				v5 = Numerics.max(sqd[level][xyz5][t]/npt[level][xyz5], sigma2[t]);
				v6 = Numerics.max(sqd[level][xyz6][t]/npt[level][xyz6], sigma2[t]);
				v7 = Numerics.max(sqd[level][xyz7][t]/npt[level][xyz7], sigma2[t]);
				v8 = Numerics.max(sqd[level][xyz8][t]/npt[level][xyz8], sigma2[t]);
			} else if (varupdate==VARLIN) {
				// weaker priors: go down with number of samples
				v1 = (sqd[level][xyz1][t] + clustersize*sigma2[t])/(npt[level][xyz1] + clustersize);
				v2 = (sqd[level][xyz2][t] + clustersize*sigma2[t])/(npt[level][xyz2] + clustersize);
				v3 = (sqd[level][xyz3][t] + clustersize*sigma2[t])/(npt[level][xyz3] + clustersize);
				v4 = (sqd[level][xyz4][t] + clustersize*sigma2[t])/(npt[level][xyz4] + clustersize);
				v5 = (sqd[level][xyz5][t] + clustersize*sigma2[t])/(npt[level][xyz5] + clustersize);
				v6 = (sqd[level][xyz6][t] + clustersize*sigma2[t])/(npt[level][xyz6] + clustersize);
				v7 = (sqd[level][xyz7][t] + clustersize*sigma2[t])/(npt[level][xyz7] + clustersize);
				v8 = (sqd[level][xyz8][t] + clustersize*sigma2[t])/(npt[level][xyz8] + clustersize);
			} else if (varupdate==VARZERO) {
				// weaker priors: go down with number of samples
				if (sqd[level][xyz1][t]==0) v1 = sigma2[t]; else v1 = sqd[level][xyz1][t]/npt[level][xyz1];
				if (sqd[level][xyz2][t]==0) v2 = sigma2[t]; else v2 = sqd[level][xyz2][t]/npt[level][xyz2];
				if (sqd[level][xyz3][t]==0) v3 = sigma2[t]; else v3 = sqd[level][xyz3][t]/npt[level][xyz3];
				if (sqd[level][xyz4][t]==0) v4 = sigma2[t]; else v4 = sqd[level][xyz4][t]/npt[level][xyz4];
				if (sqd[level][xyz5][t]==0) v5 = sigma2[t]; else v5 = sqd[level][xyz5][t]/npt[level][xyz5];
				if (sqd[level][xyz6][t]==0) v6 = sigma2[t]; else v6 = sqd[level][xyz6][t]/npt[level][xyz6];
				if (sqd[level][xyz7][t]==0) v7 = sigma2[t]; else v7 = sqd[level][xyz7][t]/npt[level][xyz7];
				if (sqd[level][xyz8][t]==0) v8 = sigma2[t]; else v8 = sqd[level][xyz8][t]/npt[level][xyz8];
			} else if (varupdate==VARCONST) {
				// weaker priors: go down with number of samples
				v1 = sigma2[t];
				v2 = sigma2[t];
				v3 = sigma2[t];
				v4 = sigma2[t];
				v5 = sigma2[t];
				v6 = sigma2[t];
				v7 = sigma2[t];
				v8 = sigma2[t];
			}
			/*
			// shouldn't we use the joint variance calculation instead?? i.e. + \sum w_ij (\mu_i-\mu_j)^2
			double v0 = 0.125*( v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 );
			
			double mu0 = 0.125*( sum[level][xyz1][t]/npt[level][xyz1] + sum[level][xyz2][t]/npt[level][xyz2] 
							   + sum[level][xyz3][t]/npt[level][xyz3] + sum[level][xyz4][t]/npt[level][xyz4] 
							   + sum[level][xyz5][t]/npt[level][xyz5] + sum[level][xyz6][t]/npt[level][xyz6] 
							   + sum[level][xyz7][t]/npt[level][xyz7] + sum[level][xyz8][t]/npt[level][xyz8] );
			*/
			if (metric==JSDIV_EXACT) {
				double mu1 = sum[level][xyz1][t]/npt[level][xyz1];
				double mu2 = sum[level][xyz2][t]/npt[level][xyz2];
				double mu3 = sum[level][xyz3][t]/npt[level][xyz3];
				double mu4 = sum[level][xyz4][t]/npt[level][xyz4];
				double mu5 = sum[level][xyz5][t]/npt[level][xyz5];
				double mu6 = sum[level][xyz6][t]/npt[level][xyz6];
				double mu7 = sum[level][xyz7][t]/npt[level][xyz7];
				double mu8 = sum[level][xyz8][t]/npt[level][xyz8];
				double sq1 = FastMath.sqrt(v1);
				double sq2 = FastMath.sqrt(v2);
				double sq3 = FastMath.sqrt(v3);
				double sq4 = FastMath.sqrt(v4);
				double sq5 = FastMath.sqrt(v5);
				double sq6 = FastMath.sqrt(v6);
				double sq7 = FastMath.sqrt(v7);
				double sq8 = FastMath.sqrt(v8);
				
				// here we estimate the JS divergence explicitly
				double xmin = Numerics.min(new double[]{mu1-3.0*sq1, mu2-3.0*sq2, mu3-3.0*sq3, mu4-3.0*sq4, mu5-3.0*sq5, mu6-3.0*sq6, mu7-3.0*sq7, mu8-3.0*sq8});
				double xmax = Numerics.max(new double[]{mu1+3.0*sq1, mu2+3.0*sq2, mu3+3.0*sq3, mu4+3.0*sq4, mu5+3.0*sq5, mu6+3.0*sq6, mu7+3.0*sq7, mu8+3.0*sq8});
				Gaussian p1 = new Gaussian(mu1, sq1);
				Gaussian p2 = new Gaussian(mu2, sq2);
				Gaussian p3 = new Gaussian(mu3, sq3);
				Gaussian p4 = new Gaussian(mu4, sq4);
				Gaussian p5 = new Gaussian(mu5, sq5);
				Gaussian p6 = new Gaussian(mu6, sq6);
				Gaussian p7 = new Gaussian(mu7, sq7);
				Gaussian p8 = new Gaussian(mu8, sq8);
				
				double jsdiv = 0.0;
				double step = 0.01*(xmax-xmin);
				for (double x=xmin;x<=xmax;x+=step) {
					double p1x = p1.value(x);
					double p2x = p2.value(x);
					double p3x = p3.value(x);
					double p4x = p4.value(x);
					double p5x = p5.value(x);
					double p6x = p6.value(x);
					double p7x = p7.value(x);
					double p8x = p8.value(x);
					double psum = (p1x+p2x+p3x+p4x+p5x+p6x+p7x+p8x)/8.0;
					jsdiv += p1x*FastMath.log(p1x/psum) + p2x*FastMath.log(p2x/psum) + p3x*FastMath.log(p3x/psum) + p4x*FastMath.log(p4x/psum)
							+ p5x*FastMath.log(p5x/psum) + p6x*FastMath.log(p6x/psum) + p7x*FastMath.log(p7x/psum) + p8x*FastMath.log(p8x/psum);
				}
				js += jsdiv*step/8.0/(double)nt;
			} else if (metric==JSDIV_TABLE) {
				double mu1 = sum[level][xyz1][t]/npt[level][xyz1];
				double mu2 = sum[level][xyz2][t]/npt[level][xyz2];
				double mu3 = sum[level][xyz3][t]/npt[level][xyz3];
				double mu4 = sum[level][xyz4][t]/npt[level][xyz4];
				double mu5 = sum[level][xyz5][t]/npt[level][xyz5];
				double mu6 = sum[level][xyz6][t]/npt[level][xyz6];
				double mu7 = sum[level][xyz7][t]/npt[level][xyz7];
				double mu8 = sum[level][xyz8][t]/npt[level][xyz8];
				double sq1 = FastMath.sqrt(v1);
				double sq2 = FastMath.sqrt(v2);
				double sq3 = FastMath.sqrt(v3);
				double sq4 = FastMath.sqrt(v4);
				double sq5 = FastMath.sqrt(v5);
				double sq6 = FastMath.sqrt(v6);
				double sq7 = FastMath.sqrt(v7);
				double sq8 = FastMath.sqrt(v8);
				
				// check for all 6-C pairs
				// ordering is: xyzp,xyzpx,xyzpy,xyzpz,xyzpxy,xyzpyz,xyzpzx,xyzpxyz
				// 				1	 2	   3	 4	   5	  6		 7		8
				double jsmax = 0.0;
				// X direction
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu2,sq2));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu3,sq3,mu5,sq5));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu4,sq4,mu7,sq7));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu6,sq6,mu8,sq8));
				// Y direction
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu3,sq3));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu5,sq5));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu4,sq4,mu6,sq6));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu7,sq7,mu8,sq8));
				// Z direction
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu4,sq4));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu7,sq7));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu3,sq3,mu6,sq6));
				jsmax = Numerics.max(jsmax, jsdtable.lookup(mu5,sq5,mu8,sq8));
				
				// more directions for 18 and 26C
				// ordering is: xyzp,xyzpx,xyzpy,xyzpz,xyzpxy,xyzpyz,xyzpzx,xyzpxyz
				// 			1	2	   3	 	4	   5	    6		 7		8
				if (connect>6) {
					// XY direction
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu5,sq5));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu4,sq4,mu8,sq8));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu3,sq3));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu7,sq7,mu6,sq6));
					// YZ direction
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu6,sq6));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu8,sq8));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu3,sq3,mu4,sq4));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu5,sq5,mu7,sq7));
					// ZX direction
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu7,sq7));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu3,sq3,mu8,sq8));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu4,sq4));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu5,sq5,mu6,sq6));
				}
				if (connect>18) {
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu1,sq1,mu8,sq8));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu2,sq2,mu6,sq6));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu3,sq3,mu7,sq7));
					jsmax = Numerics.max(jsmax, jsdtable.lookup(mu4,sq4,mu5,sq5));		
				}
				
				js += jsmax/(double)nt;
			} else {

				double npt0 = npt[level][xyz1];
				double v0 = npt0*v1;
				double mu0 = sum[level][xyz1][t];
				
				v0 += npt[level][xyz2]*v2 + npt0*npt[level][xyz2]/(npt0+npt[level][xyz2])*Numerics.square(mu0/npt0 - sum[level][xyz2][t]/npt[level][xyz2]);
				mu0 += sum[level][xyz2][t];
				npt0 += npt[level][xyz2];
				
				v0 += npt[level][xyz3]*v3 + npt0*npt[level][xyz3]/(npt0+npt[level][xyz3])*Numerics.square(mu0/npt0 - sum[level][xyz3][t]/npt[level][xyz3]);
				mu0 += sum[level][xyz3][t];
				npt0 += npt[level][xyz3];
				
				v0 += npt[level][xyz4]*v4 + npt0*npt[level][xyz4]/(npt0+npt[level][xyz4])*Numerics.square(mu0/npt0 - sum[level][xyz4][t]/npt[level][xyz4]);
				mu0 += sum[level][xyz4][t];
				npt0 += npt[level][xyz4];
				
				v0 += npt[level][xyz5]*v5 + npt0*npt[level][xyz5]/(npt0+npt[level][xyz5])*Numerics.square(mu0/npt0 - sum[level][xyz5][t]/npt[level][xyz5]);
				mu0 += sum[level][xyz5][t];
				npt0 += npt[level][xyz5];
				
				v0 += npt[level][xyz6]*v6 + npt0*npt[level][xyz6]/(npt0+npt[level][xyz6])*Numerics.square(mu0/npt0 - sum[level][xyz6][t]/npt[level][xyz6]);
				mu0 += sum[level][xyz6][t];
				npt0 += npt[level][xyz6];
				
				v0 += npt[level][xyz7]*v7 + npt0*npt[level][xyz7]/(npt0+npt[level][xyz7])*Numerics.square(mu0/npt0 - sum[level][xyz7][t]/npt[level][xyz7]);
				mu0 += sum[level][xyz7][t];
				npt0 += npt[level][xyz7];
				
				v0 += npt[level][xyz8]*v8 + npt0*npt[level][xyz8]/(npt0+npt[level][xyz8])*Numerics.square(mu0/npt0 - sum[level][xyz8][t]/npt[level][xyz8]);
				mu0 += sum[level][xyz8][t];
				npt0 += npt[level][xyz8];
				
				// normalize
				mu0 /= npt0;
				v0 /= (npt0-1);
				
				double jst = 0;
				if (metric==JSDIV_SUM) {
					jst += 0.5*(Numerics.square(sum[level][xyz1][t]/npt[level][xyz1]-mu0)/v0 + FastMath.log(v0/v1) );
					jst += 0.5*(Numerics.square(sum[level][xyz2][t]/npt[level][xyz2]-mu0)/v0 + FastMath.log(v0/v2) );
					jst += 0.5*(Numerics.square(sum[level][xyz3][t]/npt[level][xyz3]-mu0)/v0 + FastMath.log(v0/v3) );
					jst += 0.5*(Numerics.square(sum[level][xyz4][t]/npt[level][xyz4]-mu0)/v0 + FastMath.log(v0/v4) );
					jst += 0.5*(Numerics.square(sum[level][xyz5][t]/npt[level][xyz5]-mu0)/v0 + FastMath.log(v0/v5) );
					jst += 0.5*(Numerics.square(sum[level][xyz6][t]/npt[level][xyz6]-mu0)/v0 + FastMath.log(v0/v6) );
					jst += 0.5*(Numerics.square(sum[level][xyz7][t]/npt[level][xyz7]-mu0)/v0 + FastMath.log(v0/v7) );
					jst += 0.5*(Numerics.square(sum[level][xyz8][t]/npt[level][xyz8]-mu0)/v0 + FastMath.log(v0/v8) );
					
					 js += jst/(double)nt/8.0;
				} else {
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz1][t]/npt[level][xyz1]-mu0)/v0 + FastMath.log(v0/v1) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz2][t]/npt[level][xyz2]-mu0)/v0 + FastMath.log(v0/v2) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz3][t]/npt[level][xyz3]-mu0)/v0 + FastMath.log(v0/v3) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz4][t]/npt[level][xyz4]-mu0)/v0 + FastMath.log(v0/v4) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz5][t]/npt[level][xyz5]-mu0)/v0 + FastMath.log(v0/v5) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz6][t]/npt[level][xyz6]-mu0)/v0 + FastMath.log(v0/v6) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz7][t]/npt[level][xyz7]-mu0)/v0 + FastMath.log(v0/v7) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyz8][t]/npt[level][xyz8]-mu0)/v0 + FastMath.log(v0/v8) ) );
				
					js = Numerics.max(js, jst);
				}
			}
		}		
		//if (metric==JSDIV_GAUSS_GEOM) return (float)(FastMath.exp(-0.5*js/(8.0*nt)) - threshold);
		//else if (metric==JSDIV_LIN_GEOM) return (float)(1.0-js/(8.0*nt) - threshold);
		//else if (metric==JSDIV_GAUSS) return (float)(FastMath.exp(-0.5*js/nt) - threshold);
		//else return (float)(1.0-js/nt - threshold);
		
		//return (float)(FastMath.exp(-0.5*js) - threshold);
		//return (float)(JS8*FastMath.exp(-js) - threshold);
		return (float)(1.0 - FastMath.sqrt(js)/SQLN8 - threshold);
	}
	
	float regionCompetitionJSdiv(int xyzA, int xyzB, int level) {
		// joint variance
		double js = 0.0;
		for (int t=0;t<nt;t++) {
			double vA=0, vB=0;
			if (varupdate==VARMAX) {
				// strong prior influence
				vA = Numerics.max(sqd[level][xyzA][t]/npt[level][xyzA], sigma2[t]);
				vB = Numerics.max(sqd[level][xyzB][t]/npt[level][xyzB], sigma2[t]);
			} else if (varupdate==VARLIN) {
				// weaker priors: go down with number of samples
				vA = (sqd[level][xyzA][t] + clustersize*sigma2[t])/(npt[level][xyzA] + clustersize);
				vB = (sqd[level][xyzB][t] + clustersize*sigma2[t])/(npt[level][xyzB] + clustersize);
			} else if (varupdate==VARZERO) {
				// prior only for t=0
				if (sqd[level][xyzA][t]==0) vA = sigma2[t]; else vA = sqd[level][xyzA][t]/npt[level][xyzA];
				if (sqd[level][xyzB][t]==0) vB = sigma2[t]; else vB = sqd[level][xyzB][t]/npt[level][xyzB];
			}else if (varupdate==VARCONST) {
				// prior always
				vA = sigma2[t];
				vB = sigma2[t];
			}
			/*
			if (npt[level][xyzA]==1) vA = vB;
			else if (npt[level][xyzB]==1) vB = vA;
			*/
			
			if (metric==JSDIV_EXACT) {
				double muA = sum[level][xyzA][t]/npt[level][xyzA];
				double muB = sum[level][xyzB][t]/npt[level][xyzB];
				double sqA = FastMath.sqrt(vA);
				double sqB = FastMath.sqrt(vB);
				
				// here we estimate the JS divergence explicitly
				double xmin = Numerics.min(muA-3.0*sqA, muB-3.0*sqB);
				double xmax = Numerics.max(muA+3.0*sqA, muB+3.0*sqB);
				Gaussian pA = new Gaussian(muA, sqA);
				Gaussian pB = new Gaussian(muB, sqB);
				
				double jsdiv = 0.0;
				double step = 0.01*(xmax-xmin);
				for (double x=xmin;x<=xmax;x+=step) {
					double pAx = pA.value(x);
					double pBx = pB.value(x);
					jsdiv += pAx*FastMath.log(2.0*pAx/(pAx+pBx)) + pBx*FastMath.log(2.0*pBx/(pAx+pBx));
				}
				js += jsdiv*step/2.0/(double)nt;
			} else if (metric==JSDIV_TABLE) {
				double muA = sum[level][xyzA][t]/npt[level][xyzA];
				double muB = sum[level][xyzB][t]/npt[level][xyzB];
				double sqA = FastMath.sqrt(vA);
				double sqB = FastMath.sqrt(vB);
				
				js += jsdtable.lookup(muA,sqA,muB,sqB)/(double)nt;
			} else {
				// here we compute explicitly the joint distribution (!! assumes the joint distribution is normal: not true!!)
				double vAB = npt[level][xyzA]*vA + npt[level][xyzB]*vB 
							+ npt[level][xyzA]*npt[level][xyzB]/(npt[level][xyzA]+npt[level][xyzB])
							*Numerics.square(sum[level][xyzA][t]/npt[level][xyzA]-sum[level][xyzB][t]/npt[level][xyzB]);
				double muAB = sum[level][xyzA][t] + sum[level][xyzB][t];
				
				// normalize
				muAB /= (npt[level][xyzA]+npt[level][xyzB]);
				vAB /= (npt[level][xyzA]+npt[level][xyzB]-1);
	
				double jst = 0;
				if (metric==JSDIV_SUM) {
					jst += 0.5*(Numerics.square(sum[level][xyzA][t]/npt[level][xyzA]-muAB)/vAB + FastMath.log( vAB/vA ) );
					jst += 0.5*(Numerics.square(sum[level][xyzB][t]/npt[level][xyzB]-muAB)/vAB + FastMath.log( vAB/vB ) );
	
					js += jst/(double)nt/2.0;
				} else {
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyzA][t]/npt[level][xyzA]-muAB)/vAB + FastMath.log( vAB/vA ) ) );
					jst = Numerics.max(jst, 0.5*(Numerics.square(sum[level][xyzB][t]/npt[level][xyzB]-muAB)/vAB + FastMath.log( vAB/vB ) ) );
				
					js = Numerics.max(js, jst);
				}
			}
		}
		// original metric; may not be optimal... or needed at all
		//float wsize = 1.0f/(1.0f + Numerics.square( (npt[level][xyzA]+npt[level][xyzB])/(6.0f*nlb)) );
		
		//if (metric==JSDIV_GAUSS_GEOM) return (float)(FastMath.exp(-0.5*js/(2.0*nt)) - threshold);
		//else if (metric==JSDIV_LIN_GEOM) return (float)(1.0-js/(2.0*nt) - threshold);
		//else if (metric==JSDIV_GAUSS) return (float)(FastMath.exp(-0.5*js/nt) - threshold);
		//else return (float)(1.0-js/nt - threshold);
		
		//return (float)(JS2*FastMath.exp(-js) - threshold);
		//return (float)(FastMath.exp(-0.5*js) - threshold);
		return (float)(1.0 - FastMath.sqrt(js)/SQLN2 - threshold);
	}
	
	public final float[][] getSum(int lvl) { return sum[lvl]; }
	
	public final float[][] getSqd(int lvl) { return sqd[lvl]; }
	
	public final int[] getNpt(int lvl) { return npt[lvl]; }
	
	public final int[] getLbl(int lvl) { return lbl[lvl]; }
	
	public final int generateLabelsAt(int l) {
		if (lbl[l]==null) lbl[l] = new int[nx[l]*ny[l]*nz[l]];
		int label=0;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (valid[l].get(xyz)) {
				label++;
				lbl[l][xyz] = label;
			} else {
				lbl[l][xyz] = 0;
			}
		}
		return label;
	}
	
	public final int generateNextLabelsAt(int l, int prevlb) {
		if (lbl[l]==null) lbl[l] = new int[nx[l]*ny[l]*nz[l]];
		// fill in the labels from previous scale
		for (int x=0;x<nx[l+1];x++) for (int y=0;y<ny[l+1];y++) for (int z=0;z<nz[l+1];z++) {
			int xyz = x+nx[l+1]*y+nx[l+1]*ny[l+1]*z;
			if (lbl[l+1][xyz]>0) {
				int xyzn = 2*x + nx[l]*2*y + nx[l]*ny[l]*2*z;
				lbl[l][xyzn] = lbl[l+1][xyz];
				lbl[l][xyzn+1] 						= lbl[l+1][xyz];
				lbl[l][xyzn+nx[l]] 					= lbl[l+1][xyz];
				lbl[l][xyzn+nx[l]*ny[l]] 			= lbl[l+1][xyz];
				lbl[l][xyzn+1+nx[l]] 				= lbl[l+1][xyz];
				lbl[l][xyzn+nx[l]+nx[l]*ny[l]] 		= lbl[l+1][xyz];
				lbl[l][xyzn+nx[l]*ny[l]+1] 			= lbl[l+1][xyz];
				lbl[l][xyzn+1+nx[l]+nx[l]*ny[l]] 	= lbl[l+1][xyz];
			}
		}
		// use prevlb labels for the previous scale
		int label=prevlb;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			// copy down the data from previous scale
			if (lbl[l][xyz]>0) {
				// do nothing, copied above
			} else
			if (valid[l].get(xyz)) {
				// add clusters from this scale
				label++;
				lbl[l][xyz] = label;
			} else {
				lbl[l][xyz] = 0;
			}
		}
		return label;
	}
	
	public final void generateNextLevelMaps(int l) {
		
		// fill in the maps from previous scale
		for (int x=0;x<nx[l+1];x++) for (int y=0;y<ny[l+1];y++) for (int z=0;z<nz[l+1];z++) {
			int xyz = x+nx[l+1]*y+nx[l+1]*ny[l+1]*z;
			if (lbl[l+1][xyz]>0) {
				int xyzn = 2*x + nx[l]*2*y + nx[l]*ny[l]*2*z;
				for (int t=0;t<nt;t++) {
					sum[l][xyzn][t] 						= sum[l+1][xyz][t];
					sum[l][xyzn+1][t] 						= sum[l+1][xyz][t];
					sum[l][xyzn+nx[l]][t] 					= sum[l+1][xyz][t];
					sum[l][xyzn+nx[l]*ny[l]][t] 			= sum[l+1][xyz][t];
					sum[l][xyzn+1+nx[l]][t] 				= sum[l+1][xyz][t];
					sum[l][xyzn+nx[l]+nx[l]*ny[l]][t] 		= sum[l+1][xyz][t];
					sum[l][xyzn+nx[l]*ny[l]+1][t] 			= sum[l+1][xyz][t];
					sum[l][xyzn+1+nx[l]+nx[l]*ny[l]][t] 	= sum[l+1][xyz][t];
	
					sqd[l][xyzn][t] 						= sqd[l+1][xyz][t];
					sqd[l][xyzn+1][t] 						= sqd[l+1][xyz][t];
					sqd[l][xyzn+nx[l]][t] 					= sqd[l+1][xyz][t];
					sqd[l][xyzn+nx[l]*ny[l]][t] 			= sqd[l+1][xyz][t];
					sqd[l][xyzn+1+nx[l]][t] 				= sqd[l+1][xyz][t];
					sqd[l][xyzn+nx[l]+nx[l]*ny[l]][t] 		= sqd[l+1][xyz][t];
					sqd[l][xyzn+nx[l]*ny[l]+1][t] 			= sqd[l+1][xyz][t];
					sqd[l][xyzn+1+nx[l]+nx[l]*ny[l]][t] 	= sqd[l+1][xyz][t];
				}
				npt[l][xyzn] 						= npt[l+1][xyz];
				npt[l][xyzn+1] 						= npt[l+1][xyz];
				npt[l][xyzn+nx[l]] 					= npt[l+1][xyz];
				npt[l][xyzn+nx[l]*ny[l]] 			= npt[l+1][xyz];
				npt[l][xyzn+1+nx[l]] 				= npt[l+1][xyz];
				npt[l][xyzn+nx[l]+nx[l]*ny[l]] 		= npt[l+1][xyz];
				npt[l][xyzn+nx[l]*ny[l]+1] 			= npt[l+1][xyz];
				npt[l][xyzn+1+nx[l]+nx[l]*ny[l]] 	= npt[l+1][xyz];
			}
		}
		return;
	}
	
	public final int[][][] exportLabelsFrom(int lb) {
		int[][][] labeling, nextlb;
		
		labeling = new int[nx[lb]][ny[lb]][nz[lb]];
		for (int x=0;x<nx[lb];x++) for (int y=0;y<ny[lb];y++) for (int z=0;z<nz[lb];z++) {
			int xyz = x+nx[lb]*y+nx[lb]*ny[lb]*z;
			labeling[x][y][z] = lbl[lb][xyz];
		}
		for (int l=lb-1;l>=0;l--) {
			nextlb= new int[nx[l]][ny[l]][nz[l]];
			// fill in the labels from previous scale
			for (int x=0;x<nx[l+1];x++) for (int y=0;y<ny[l+1];y++) for (int z=0;z<nz[l+1];z++) {
				if (labeling[x][y][z]>0) {
					nextlb[2*x][2*y][2*z]			= labeling[x][y][z];
					nextlb[2*x+1][2*y][2*z]		= labeling[x][y][z];
					nextlb[2*x][2*y+1][2*z]		= labeling[x][y][z];
					nextlb[2*x][2*y][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y+1][2*z] 		= labeling[x][y][z];
					nextlb[2*x][2*y+1][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y+1][2*z+1] 	= labeling[x][y][z];
				}
			}
			labeling = nextlb;
		}
		return labeling;
	}
	
	public final byte[][][] exportOctreeStructure(int lb) {
		byte[][][] labeling, nextlb;
		
		labeling = new byte[nx[lb]][ny[lb]][nz[lb]];
		for (int x=0;x<nx[lb];x++) for (int y=0;y<ny[lb];y++) for (int z=0;z<nz[lb];z++) {
			int xyz = x+nx[lb]*y+nx[lb]*ny[lb]*z;
			if (npt[lb][xyz]>0) labeling[x][y][z] = (byte)(lb+1);
		}
		for (int l=lb-1;l>=0;l--) {
			nextlb= new byte[nx[l]][ny[l]][nz[l]];
			// fill in the labels from previous scale
			for (int x=0;x<nx[l+1];x++) for (int y=0;y<ny[l+1];y++) for (int z=0;z<nz[l+1];z++) {
				if (labeling[x][y][z]>0) {
					nextlb[2*x][2*y][2*z]			= labeling[x][y][z];
					nextlb[2*x+1][2*y][2*z]		= labeling[x][y][z];
					nextlb[2*x][2*y+1][2*z]		= labeling[x][y][z];
					nextlb[2*x][2*y][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y+1][2*z] 		= labeling[x][y][z];
					nextlb[2*x][2*y+1][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y][2*z+1] 		= labeling[x][y][z];
					nextlb[2*x+1][2*y+1][2*z+1] 	= labeling[x][y][z];
				}
			}
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (npt[l][xyz]>0 && nextlb[x][y][z]==0) nextlb[x][y][z] = (byte)(l+1);
			}
			labeling = nextlb;
		}
		return labeling;
	}
	
	public final int getDimX(int l) { return nx[l]; } 
	
	public final int getDimY(int l) { return ny[l]; } 
	
	public final int getDimZ(int l) { return nz[l]; } 
	
	public final int getDimT() { return nt; } 
	
	public final float[][][][] exportAvgAtLevel(int l) {
		float[][][][] avgs = new float[nx[l]][ny[l]][nz[l]][nt];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			for (int t=0;t<nt;t++)
				avgs[x][y][z][t] = sum[l][xyz][t]/Numerics.max(1.0f,npt[l][xyz]);
		}
		return avgs;
	}
	
	public final float[][][][] exportStdAtLevel(int l) {
		float[][][][] stds = new float[nx[l]][ny[l]][nz[l]][nt];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			for (int t=0;t<nt;t++)
				stds[x][y][z][t] = (float)FastMath.sqrt(sqd[l][xyz][t]/Numerics.max(1,npt[l][xyz]-1));
		}
		return stds;
	}
	
	public final int[][][] exportNptAtLevel(int l) {
		int[][][] npts = new int[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			npts[x][y][z] = npt[l][xyz];
		}
		return npts;
	}
	
	public final float[][][] exportLogNptAtLevel(int l) {
		float[][][] npts = new float[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (npt[l][xyz]>0) npts[x][y][z] = (float)(1+FastMath.log(npt[l][xyz]));
		}
		return npts;
	}
	
	public final int[][][] exportLblAtLevel(int l) {
		int[][][] lbls = new int[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			lbls[x][y][z] = lbl[l][xyz];
		}
		return lbls;
	}

	public final float[][][] exportProbaAtLevel(int l, float[][] img) {
		float[][][] probas = new float[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>0) {
				int xyz0 = x*CL[l] + nx[0]*y*CL[l] + nx[0]*ny[0]*z*CL[l];
				probas[x][y][z] = 1.0f;
				for (int t=0;t<nt;t++) {
					probas[x][y][z] *= (float)FastMath.exp( - 0.5*(img[xyz0][t]-sum[l][xyz][t]/npt[l][xyz])
																	/Numerics.max( (sqd[l][xyz][t]/Numerics.max(1,npt[l][xyz]-1)), sigma2[t]) );
				}
			}
		}
		return probas;
	}
	
	public final int removeSmallClusters(int l, int clustersize, int ncluster) {
		// finds the boundary between small and large clusters	
		int nr = 0;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>0 && npt[l][xyz]<clustersize) {
				lbl[l][xyz]=0;
				nr++;
			}				
		}
		System.out.println("removed "+nr+" clusters under "+clustersize+" voxels");
		// relabel
		int nc=0;
		int[] mapped = new int[ncluster+1];
		for (int n=0;n<ncluster+1;n++) mapped[n] = 0;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>0) {
				if (mapped[lbl[l][xyz]]==0) {
					nc++;
					mapped[lbl[l][xyz]] = nc;
				}
				lbl[l][xyz] = mapped[lbl[l][xyz]];
			}		
		}
		return nc;
	}
}
