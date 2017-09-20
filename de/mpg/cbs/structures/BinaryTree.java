package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;

/**
 *
 *  Binary sorting trees, either min-trees or max-trees.
 *	<p>
 *  Values are sorted in a binary tree which require that each parent node is lower (resp. higher) than its children.
 *  The root of the tree is the lowest (resp. highest) value in the tree, and all operations (adding or removing a point)
 *  have <i>O( N log N )<i> complexity.
 *	<p>
 *	These trees are used principally in fast marching methods.
 *
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class BinaryTree {
	
	private float[] 	tree;
	private int[][] 	index;
	// note: alternatively, we could use ArrayLists (save some space, might be slower)
	
	private int 		N;
	private int 		Nc;
	private int			Nmax;
	private int			Nup;
	private int 		minormax;
	private static int memory;
	
	public static final int MAXTREE = 1;
	public static final int MINTREE = -1;

	public static final int FIXEDSIZE = 0;
	public static final int ADAPTATIVE = 1;

	public BinaryTree(int Nsize, int Ncoord, int type) {
		minormax = type;
		Nc = Ncoord;
		N = 0;
		Nmax = Nsize;
		Nup = Nsize;
		memory = FIXEDSIZE;
		tree = new float[Nmax];
		index = new int[Nmax][Ncoord];
	}
	
	public BinaryTree(int Nsize, int Ncoord, int type, int mem) {
		minormax = type;
		Nc = Ncoord;
		N = 0;
		Nmax = Nsize;
		Nup = Nsize;
		memory = mem;
		tree = new float[Nmax];
		index = new int[Nmax][Ncoord];
	}
	
	public void finalize() {
		tree = null;
		index = null;
	}
	
	/**
	 *  to reset the binary tree
	 */
	public final void reset() {
		N = 0;
	}
	
	/**
	 *  to set the binary tree type
	 */
	public final void setMaxTree() {
		minormax = MAXTREE;
	}
	public final void setMinTree() {
		minormax = MINTREE;
	}
	
	/**
	 *  add a new value into the binary tree
	 */
	public final void addValue(float val, int[] x) {
		int     n,n2;
		int 	stop=0;
	
		// check for size
		if ( (memory==ADAPTATIVE) && (N >= Nmax-2) ) {
			float[] tmp = new float[Nmax];
			int[][] ind = new int[Nmax][Nc];
			for (int i=1;i<=N;i++) {
				tmp[i] = tree[i];
				for (int j=0;j<Nc;j++) ind[i][j] = index[i][j];
			}
			Nmax = Nmax + Nup;
			tree = new float[Nmax];
			index = new int[Nmax][Nc];
			for (int i=1;i<=N;i++) {
				tree[i] = tmp[i];
				for (int j=0;j<Nc;j++) index[i][j] = ind[i][j];
			}
			tmp = null;
			ind = null;
		}
				
		
		// start from the new value
		n = N+1;
		n2 = n/2;
		
		// find the insert point and displace values
		if (minormax==MINTREE) {
			while (stop==0) {
				if (n2>0) {
					if ( val < tree[n2] ) {
						tree[n] = tree[n2];
						for (int m=0;m<Nc;m++)
							index[n][m] = index[n2][m];
						
						n = n2;
						n2 = n/2;
					} else {
						stop = 1;
					}
				} else {
					stop = 1;
				}
			}
		} else if (minormax==MAXTREE) {
			while (stop==0) {
				if (n2>0) {
					if ( val > tree[n2] ) {
						tree[n] = tree[n2];
						for (int m=0;m<Nc;m++)
							index[n][m] = index[n2][m];
						
						n = n2;
						n2 = n/2;
					} else {
						stop = 1;
					}
				} else {
					stop = 1;
				}
			}
		}
		
		// insert new value
		tree[n] = val;
		for (int m=0;m<Nc;m++) 
			index[n][m] = x[m];
	
		// update the size	
		N = N+1;
		
		return;
	}//addValue
	
	/**
	 *  remove the first value from the tree
	 */
	public final void removeFirst() {
		int     n,n2;
		float   val;
		int[] 	x = new int[Nc];
		int stop = 0;
	
		// replace the first value with the last attributed one
		val = tree[N];
		for (int m=0;m<Nc;m++)
			x[m] = index[N][m];
		
		// start from the new value
		n = 1;
		n2 = 2*n;
		
		// find the insert point and displace values
		if (minormax==MINTREE) {
			while (stop==0) {
				if ( n2 < N ) {
					if ( n2+1 < N ) {
						if ( (val > tree[n2]) || (val > tree[n2+1]) ){
							if ( tree[n2] > tree[n2+1] ) {
								tree[n] = tree[n2+1];
								for (int m=0;m<Nc;m++) 
									index[n][m] = index[n2+1][m];
								n = n2+1;
								n2 = 2*n;
							} else {
								tree[n] = tree[n2];
								for (int m=0;m<Nc;m++) 
									index[n][m] = index[n2][m];
								n = n2;
								n2 = 2*n;
							}
						} else {
							stop = 1;
						}
					} else {
						if ( val > tree[n2] ) {
							tree[n] = tree[n2];
							for (int m=0;m<Nc;m++) 
								index[n][m] = index[n2][m];
							n = n2;
							stop = 1;
						} else {
							stop =1;
						}
					}
				} else {
					stop = 1;
				}
			}
		} else if (minormax==MAXTREE) {
			while (stop==0) {
				if ( n2 < N ) {
					if ( n2+1 < N ) {
						if ( (val < tree[n2]) || (val < tree[n2+1]) ){
							if ( tree[n2] < tree[n2+1] ) {
								tree[n] = tree[n2+1];
								for (int m=0;m<Nc;m++) 
									index[n][m] = index[n2+1][m];
								n = n2+1;
								n2 = 2*n;
							} else {
								tree[n] = tree[n2];
								for (int m=0;m<Nc;m++) 
									index[n][m] = index[n2][m];
								n = n2;
								n2 = 2*n;
							}
						} else {
							stop = 1;
						}
					} else {
						if ( val < tree[n2] ) {
							tree[n] = tree[n2];
							for (int m=0;m<Nc;m++) 
								index[n][m] = index[n2][m];
							n = n2;
							stop = 1;
						} else {
							stop =1;
						}
					}
				} else {
					stop = 1;
				}
			}
		}
		// insert new value
		tree[n] = val;
		for (int m=0;m<Nc;m++) 
			index[n][m] = x[m];
	
		// update the number	
		N = N-1;
		
		return;
	}// removeFirstValue

	/**
	 * return the first value and its coordinates
	 */
	public final float getFirst() {
		return tree[1];
	}
	public final int[] getFirstIndices() {
		return index[1];
	}
	public final int getFirstIndex(int m) {
		return index[1][m];
	}
	
	/**
	 * various utilities
	 */
	public final void print() {
		int i;
		int n;
			
		n = 2;
		for (i=1;i<=N;i++) {
			System.out.print("  "+(int)(100*tree[i]));
			if ( ( i%n ) == n-1) {
				// odd number
				System.out.print("\n");
				n = 2*n;
			}
		}
		return;
	}//print
	
	/**
	 *  check the binary tree property
	 */
	public final boolean isNotEmpty() {
		if (N>0) return true;
		else return false;
	}
	
	/**
	 *  check the binary tree property
	 */
	public final boolean isBinaryTree() {
		int i;
		int n;
		boolean isBinary=true;
				
		if (minormax==MINTREE) {
			for (i=2;i<=N;i++) {
				n = i/2;
				if ( tree[n] > tree[i] ) {
					// odd number
					isBinary = false;
					break;
				}
			}
		} else if (minormax==MAXTREE) {
			for (i=2;i<=N;i++) {
				n = i/2;
				if ( tree[n] < tree[i] ) {
					// odd number
					isBinary = false;
					break;
				}
			}
		}
		return isBinary;
	}//isBinaryTree

}
