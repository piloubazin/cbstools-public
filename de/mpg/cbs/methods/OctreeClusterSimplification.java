package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class handles octree-style grouping for JSDiv clustering
 *
 *	@version    March 2015
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class OctreeClusterSimplification {

	// store the multiscale data
	private float[][] sum;
	private float[][] sqd;
	private BitSet[] valid;
	private int[][] npt;
	private int[][] lbl;
	
	private int[] nx;
	private int[] ny;
	private int[] nz;
	
	private	static 	int 		LEVELS;
	public static int[]		CS;
	public static int[]		CL;
	
	private float sigma2;
	private float threshold;
	
	private static final boolean	debug = false;
	private static final boolean	verbose = true;	
	
	public OctreeClusterSimplification(float[] image, boolean[] mask, int nix, int niy, int niz, float stdev, float thres) {
		
		sigma2 = stdev*stdev;
		threshold  = thres;
		
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
		
		sum = new float[LEVELS][];
		sqd = new float[LEVELS][];
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
		// work on copy
		sum[0] = new float[nx[0]*ny[0]*nz[0]];
		sqd[0] = new float[nx[0]*ny[0]*nz[0]];
		npt[0] = new int[nx[0]*ny[0]*nz[0]];
		valid[0] = new BitSet(nx[0]*ny[0]*nz[0]);
		for (int x=0;x<nx[0];x++) for (int y=0;y<ny[0];y++) for (int z=0;z<nz[0];z++) {
			int xyz = x+nx[0]*y+nx[0]*ny[0]*z;
			if (mask[xyz]) {
				sum[0][xyz] = image[xyz];
				valid[0].set(xyz, true);
			} else {
				sum[0][xyz] = 0.0f;
				valid[0].set(xyz, false);
			}
			sqd[0][xyz] = 0.0f;
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
	public int mergeToJensenShannonDistance() {
		boolean stop=false;
		int maxlevel = 0;
		for (int l=1;l<LEVELS && !stop;l++) {
			if (verbose) System.out.print("level "+l);
			
			// create the arrays here (so it doesn't try more than necessary)
			sum[l] = new float[nx[l]*ny[l]*nz[l]];
			sqd[l] = new float[nx[l]*ny[l]*nz[l]];
			npt[l] = new int[nx[l]*ny[l]*nz[l]];
			valid[l] = new BitSet(nx[l]*ny[l]*nz[l]);
			
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
					
								float js = jensenShannonDivergence(xyzp, xyzpx, l-1);
					js = Numerics.min(js, jensenShannonDivergence(xyzp, xyzpy, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzp, xyzpz, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpx, xyzpxy, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpx, xyzpzx, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpy, xyzpxy, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpy, xyzpyz, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpz, xyzpyz, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpz, xyzpzx, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpxy, xyzpxyz, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpyz, xyzpxyz, l-1));
					js = Numerics.min(js, jensenShannonDivergence(xyzpzx, xyzpxyz, l-1));
					
					if (js<0) merge = false;
					//else if (verbose) System.out.print("js: "+js+" ");
						
					// merging: add the sum, sqd incrementally for each leaf with Chan's formula
					if (merge) {
						// add progressively all children
						float cs = CS[l-1];
						sum[l][xyz] = sum[l-1][xyzp];
						sqd[l][xyz] = sqd[l-1][xyzp];
						
						// update first the square distance to keep the valid sum and the added sum separate
						sqd[l][xyz] += sqd[l-1][xyzpx] + 1*cs*cs/(1*cs+cs)*Numerics.square(sum[l][xyz]/(1*cs)-sum[l-1][xyzpx]/cs);
						sum[l][xyz] += sum[l-1][xyzpx];
						
						sqd[l][xyz] += sqd[l-1][xyzpy] + 2*cs*cs/(2*cs+cs)*Numerics.square(sum[l][xyz]/(2*cs)-sum[l-1][xyzpy]/cs);
						sum[l][xyz] += sum[l-1][xyzpy];
						
						sqd[l][xyz] += sqd[l-1][xyzpz] + 3*cs*cs/(3*cs+cs)*Numerics.square(sum[l][xyz]/(3*cs)-sum[l-1][xyzpz]/cs);
						sum[l][xyz] += sum[l-1][xyzpz];
						
						sqd[l][xyz] += sqd[l-1][xyzpxy] + 4*cs*cs/(4*cs+cs)*Numerics.square(sum[l][xyz]/(4*cs)-sum[l-1][xyzpxy]/cs);
						sum[l][xyz] += sum[l-1][xyzpxy];
						
						sqd[l][xyz] += sqd[l-1][xyzpyz] + 5*cs*cs/(5*cs+cs)*Numerics.square(sum[l][xyz]/(5*cs)-sum[l-1][xyzpyz]/cs);
						sum[l][xyz] += sum[l-1][xyzpyz];
						
						sqd[l][xyz] += sqd[l-1][xyzpzx] + 6*cs*cs/(6*cs+cs)*Numerics.square(sum[l][xyz]/(6*cs)-sum[l-1][xyzpzx]/cs);
						sum[l][xyz] += sum[l-1][xyzpzx];
						
						sqd[l][xyz] += sqd[l-1][xyzpxyz] + 7*cs*cs/(7*cs+cs)*Numerics.square(sum[l][xyz]/(7*cs)-sum[l-1][xyzpxyz]/cs);
						sum[l][xyz] += sum[l-1][xyzpxyz];
						
						valid[l].set(xyz);
						npt[l][xyz] = CS[l];
						stop = false;
						maxlevel = l;
					}
				}
			}
		}
		return maxlevel;
	}
	
	
	/** group similar values in octree fashion
	 *	based on Jensen-Shannon divergence
	 */
	public int competitionWithJensenShannonDistance(int l, int nclusters) {
		// generate the merged labeling at level l
		int nlb = generateNextLabelsAt(l, nclusters);
		
		// keep the cluster level data in separate (duplicate) arrays
		float[] listsum = new float[nclusters+1];
		float[] listsqd = new float[nclusters+1];
		int[] listnpt = new int[nclusters+1];
		
		// put into the tree any regions at level l with a neighbor at lvl+1
		BinaryHeapPair tree = new BinaryHeapPair((nlb-nclusters)/16, (nlb-nclusters)/16, BinaryHeapPair.MAXTREE);
		for (int x=1;x<nx[l]-1;x++) for (int y=1;y<ny[l]-1;y++) for (int z=1;z<nz[l]-1;z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (lbl[l][xyz]>0 && lbl[l][xyz]<=nclusters) {
				// look for neighbors from new scale
				int xyzpx = xyz+1;
				int xyzpy = xyz+nx[l];
				int xyzpz = xyz+nx[l]*ny[l];
				int xyzmx = xyz-1;
				int xyzmy = xyz-nx[l];
				int xyzmz = xyz-nx[l]*ny[l];
				
				if (lbl[l][xyzpx]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzpx, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpx])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzpx, xyz);
				}
				if (lbl[l][xyzmx]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzmx, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmx])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzmx, xyz);
				}
				if (lbl[l][xyzpy]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzpy, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpy])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzpy, xyz);
				}
				if (lbl[l][xyzmy]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzmy, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmy])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzmy, xyz);
				}
				if (lbl[l][xyzpz]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzpz, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpz])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzpz, xyz);
				}
				if (lbl[l][xyzmz]>nclusters) {
					float js = jensenShannonDivergence(xyz, xyzmz, l);
					float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmz])/(6.0f*nlb)) );
					if (js>0) tree.addValue(wsize*js, xyzmz, xyz);
				}
				
				// build the global label list
				listsum[lbl[l][xyz]] = sum[l][xyz];
				listsqd[lbl[l][xyz]] = sqd[l][xyz];
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
				sqd[l][xyzngb] = listsqd[lbngb];
				sum[l][xyzngb] = listsum[lbngb];
				npt[l][xyzngb] = listnpt[lbngb];
					
				// recompute the score: it may have changed
				float wsizeup = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzngb])/(6.0f*nlb)) );
				float jsupdate = jensenShannonDivergence(xyz, xyzngb, l);
				if (wsizeup*jsupdate>=0.9f*score) {
				
					// merge with label
					int lb = lbl[l][xyzngb];
					lbl[l][xyz] = lb;
					
					// update the global lists
					listsqd[lb] += sqd[l][xyz] + listnpt[lb]*npt[l][xyz]/(listnpt[lb]+npt[l][xyz])*Numerics.square(sum[l][xyz]/npt[l][xyz]-listsum[lb]/listnpt[lb]);
					listsum[lb] += sum[l][xyz];
					listnpt[lb] += npt[l][xyz];
					
					// update also current point for the next js computation
					sqd[l][xyz] = listsqd[lb];
					sum[l][xyz] = listsum[lb];
					npt[l][xyz] = listnpt[lb];
					
					// find neighboring labels not yet relabeled
					int xyzpx = xyz+1;
					int xyzpy = xyz+nx[l];
					int xyzpz = xyz+nx[l]*ny[l];
					int xyzmx = xyz-1;
					int xyzmy = xyz-nx[l];
					int xyzmz = xyz-nx[l]*ny[l];
	
					if (lbl[l][xyzpx]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzpx, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpx])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzpx, xyz);
					}
					if (lbl[l][xyzmx]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzmx, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmx])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzmx, xyz);
					}
					if (lbl[l][xyzpy]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzpy, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpy])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzpy, xyz);
					}
					if (lbl[l][xyzmy]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzmy, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmy])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzmy, xyz);
					}
					if (lbl[l][xyzpz]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzpz, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzpz])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzpz, xyz);
					}
					if (lbl[l][xyzmz]>nclusters) {
						float js = jensenShannonDivergence(xyz, xyzmz, l);
						float wsize = 1.0f/(1.0f + Numerics.square( (npt[l][xyz]+npt[l][xyzmz])/(6.0f*nlb)) );
						if (js>0) tree.addValue(wsize*js, xyzmz, xyz);
					}
				} else if (wsizeup*jsupdate>0) {
					tree.addValue(wsizeup*jsupdate, xyz, xyzngb);
				}
			}
		}
		tree.finalize();
		tree = null;
		
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
				sum[l][xyz] = listsum[lb];
				sqd[l][xyz] = listsqd[lb];
				npt[l][xyz] = listnpt[lb];
			}
		}
		return maxlb;
	}
	
	float jensenShannonDivergence(int xyzA, int xyzB, int level) {
		// cell sizes
		//float csA = CS[level];
		//float csB = CS[level];
		float csA = npt[level][xyzA];
		float csB = npt[level][xyzB];
		//System.out.print("sizes: "+csA+","+csB+"\n");
		//System.out.flush();
		// joint variance
		double vA = Numerics.max(sqd[level][xyzA]/csA, sigma2);
		double vB = Numerics.max(sqd[level][xyzB]/csB, sigma2);
		double vAB = 0.25*( vA + vB );
		//System.out.print("variances: "+vA+","+vB+","+vAB+"\n");
		//System.out.flush();
		
		double js = 0.25*Numerics.square(sum[level][xyzA]/csA-sum[level][xyzB]/csB)/vAB 
					- 0.5*FastMath.log(vA/vAB*vB/vAB) + FastMath.log(2.0);
		//System.out.print("js dist: "+js+"\n");
		//System.out.flush();
		return (float)FastMath.exp(-0.5*js) - threshold;
		//return (float)FastMath.exp(-0.5*Numerics.square(cellA.sum-cellB.sum)/sigma2);
	}
	
	public final float[] getSum(int lvl) { return sum[lvl]; }
	
	public final float[] getSqd(int lvl) { return sqd[lvl]; }
	
	public final int[] getNpt(int lvl) { return npt[lvl]; }
	
	public final int[] getLbl(int lvl) { return lbl[lvl]; }
	
	public final int getSize(int lvl) { return CS[lvl]; }
	
	public final int[][][] exportLevelMap(int maxlevel) {
		int[][][] levels = new int[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									levels[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = l;
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return levels;
	}
	
	public final float[][][] exportSumMap(int maxlevel) {
		float[][][] sums = new float[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									sums[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = sum[l][xyz];
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return sums;
	}
	
	public final float[][][] exportSqdMap(int maxlevel) {
		float[][][] sqds = new float[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									sqds[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = sqd[l][xyz];
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return sqds;
	}
	
	public final float[][][] exportAvgMap(int maxlevel) {
		float[][][] avgs = new float[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									avgs[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = sum[l][xyz]/CS[l];
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return avgs;
	}
	
	public final float[][][] exportStdMap(int maxlevel) {
		float[][][] stds = new float[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									stds[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = (float)FastMath.sqrt(sqd[l][xyz]/CS[l]);
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return stds;
	}
		
	public final int[][][] exportLabelMap(int maxlevel) {
		int[][][] labels = new int[nx[0]][ny[0]][nz[0]];
		BitSet found = new BitSet(nx[0]*ny[0]*nz[0]);
		found.clear();
		
		int label=0;
		for (int l=maxlevel;l>=0;l--) {
			if (verbose) System.out.print("level "+l);
			
			for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
				int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
				if (valid[l].get(xyz)) {
					int xyz0 = (int)(CL[l]*x + nx[0]*CL[l]*y + nx[0]*ny[0]*CL[l]*z);
					if (!found.get(xyz0)) {
						label++;
						for (int i=0;i<CL[l];i++) {
							for (int j=0;j<CL[l];j++) { 
								for (int k=0;k<CL[l];k++) {
									labels[CL[l]*x + i][CL[l]*y + j][CL[l]*z + k] = label;
									found.set(xyz0 + i + nx[0]*j + nx[0]*ny[0]*k);
								}
							}
						}
					}
				}
			}
		}
		return labels;
	}
	/*
	public final int generateLabelsAt(int[] labels, int l) {
		int label=0;
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			if (valid[l].get(xyz)) {
				label++;
				labels[xyz] = label;
			} else {
				labels[xyz] = 0;
			}
		}
		return label;
	}
	*/
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
				sum[l][xyzn] 						= sum[l+1][xyz];
				sum[l][xyzn+1] 						= sum[l+1][xyz];
				sum[l][xyzn+nx[l]] 					= sum[l+1][xyz];
				sum[l][xyzn+nx[l]*ny[l]] 			= sum[l+1][xyz];
				sum[l][xyzn+1+nx[l]] 				= sum[l+1][xyz];
				sum[l][xyzn+nx[l]+nx[l]*ny[l]] 		= sum[l+1][xyz];
				sum[l][xyzn+nx[l]*ny[l]+1] 			= sum[l+1][xyz];
				sum[l][xyzn+1+nx[l]+nx[l]*ny[l]] 	= sum[l+1][xyz];

				sqd[l][xyzn] 						= sqd[l+1][xyz];
				sqd[l][xyzn+1] 						= sqd[l+1][xyz];
				sqd[l][xyzn+nx[l]] 					= sqd[l+1][xyz];
				sqd[l][xyzn+nx[l]*ny[l]] 			= sqd[l+1][xyz];
				sqd[l][xyzn+1+nx[l]] 				= sqd[l+1][xyz];
				sqd[l][xyzn+nx[l]+nx[l]*ny[l]] 		= sqd[l+1][xyz];
				sqd[l][xyzn+nx[l]*ny[l]+1] 			= sqd[l+1][xyz];
				sqd[l][xyzn+1+nx[l]+nx[l]*ny[l]] 	= sqd[l+1][xyz];

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
	
	
	public final int getDimX(int l) { return nx[l]; } 
	
	public final int getDimY(int l) { return ny[l]; } 
	
	public final int getDimZ(int l) { return nz[l]; } 
	
	public final float[][][] exportAvgAtLevel(int l) {
		float[][][] avgs = new float[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			avgs[x][y][z] = sum[l][xyz]/Numerics.max(1.0f,npt[l][xyz]);
		}
		return avgs;
	}
	
	public final float[][][] exportStdAtLevel(int l) {
		float[][][] stds = new float[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			stds[x][y][z] = (float)FastMath.sqrt(sqd[l][xyz]/Numerics.max(1,npt[l][xyz]-1));
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
	
	public final int[][][] exportLblAtLevel(int l) {
		int[][][] lbls = new int[nx[l]][ny[l]][nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			lbls[x][y][z] = lbl[l][xyz];
		}
		return lbls;
	}
	
	public final int[] copyLbl(int l) {
		int[] lbls = new int[nx[l]*ny[l]*nz[l]];
		
		for (int x=0;x<nx[l];x++) for (int y=0;y<ny[l];y++) for (int z=0;z<nz[l];z++) {
			int xyz = x+nx[l]*y+nx[l]*ny[l]*z;
			lbls[xyz] = lbl[l][xyz];
		}
		return lbls;
	}

}
