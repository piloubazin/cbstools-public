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
 *	These trees are used principally in fast marching methods. This specific tree has four indices for use in Toads,
 *	other versions may be needed (@see BinaryTree for a more generic implementation) 
 *
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class BinaryHeap4D {
	
	private float[] val;
	private short[] 	x;
	private short[] 	y;
	private short[] 	z;
	private short[] 	k;
	
	private int 		currentSize;
	private int			capacity;
	private int 		minormax;
	
	public static final	int MINTREE = -1;
	public static final	int MAXTREE = 1;
	
	public BinaryHeap4D(int Nsize, int type) {
		currentSize = 0;
		capacity = Nsize;
		minormax = type;
		val = new float[capacity+1];
		x = new short[capacity+1];
		y = new short[capacity+1];
		z = new short[capacity+1];
		k = new short[capacity+1];
		if (minormax==MINTREE)
			val[0] = -1e12f;
		else if (minormax==MAXTREE)
			val[0] = 1e12f;
	}
	
	public void finalize() {
		val = null;
		x = null;
		y = null;
		z = null;
		k = null;
	}
	
	/**
	 *  to reset the binary tree
	 */
	public final void reset() {
		currentSize = 0;
	}
	
	/**
	 *  to set the binary tree type
	 */
	public final void setMaxTree() {
		minormax = MAXTREE;
		val[0] = 1e12f;
	}
	public final void setMinTree() {
		minormax = MINTREE;
		val[0] = -1e12f;
	}
	
	/**
	 *  add a new value into the binary tree
	 */
	public final void addValue(float val_, int x_, int y_, int z_, int k_) {
		addValue(val_,(short)x_,(short)y_,(short)z_,(short)k_);
	}
	/**
	 *  add a new value into the binary tree
	 */
	public final void addValue(float val_, short x_, short y_, short z_, short k_) {
		// check for size
		if  (currentSize == val.length - 1) {
			float[] oldVal = val;
			short[] oldX = x;
			short[] oldY = y;
			short[] oldZ = z;
			short[] oldK = k;
			val = new float[currentSize+capacity];
			x = new short[currentSize+capacity];
			y = new short[currentSize+capacity];
			z = new short[currentSize+capacity];
			k = new short[currentSize+capacity];
			for (int i=0;i<oldVal.length;i++) {
				val[i] = oldVal[i];
				x[i] = oldX[i];
				y[i] = oldY[i];
				z[i] = oldZ[i];
				k[i] = oldK[i];
			}
		}
		// insert new  point into the proper location		
		int hole = ++currentSize;
		
		if (minormax==MINTREE) {
			for ( ; val_ < val[ hole/2 ]; hole /= 2 ) {
				val[hole] = val[hole/2];
				x[hole] = x[hole/2];
				y[hole] = y[hole/2];
				z[hole] = z[hole/2];
				k[hole] = k[hole/2];
			}
			val[hole] = val_;
			x[hole] = x_;
			y[hole] = y_;
			z[hole] = z_;
			k[hole] = k_;
		} else if (minormax==MAXTREE) {
			for ( ; val_ > val[ hole/2 ]; hole /= 2 ) {
				val[hole] = val[hole/2];
				x[hole] = x[hole/2];
				y[hole] = y[hole/2];
				z[hole] = z[hole/2];
				k[hole] = k[hole/2];
			}
			val[hole] = val_;
			x[hole] = x_;
			y[hole] = y_;
			z[hole] = z_;
			k[hole] = k_;
		}
		
		return;
	}//addValue
	
	/**
	 *  remove the first value from the tree
	 */
	public final void removeFirst() {
		int hole = 1;
		
		val[hole] = val[currentSize];
		x[hole] = x[currentSize];
		y[hole] = y[currentSize];
		z[hole] = z[currentSize];
		k[hole] = k[currentSize];
		currentSize--;
		
		int child;
		float tmp = val[hole];
		short tmpX = x[hole];
		short tmpY = y[hole];
		short tmpZ = z[hole];
		short tmpK = k[hole];
		
		if (minormax==MINTREE) {
			for ( ; hole*2 <= currentSize; hole = child ) {
				child = hole*2;
				if (child != currentSize && val[child+1]<val[child])
					child++;
				if ( val[child]<tmp ) {
					val[ hole ] = val[ child ];
					x[ hole ] = x[ child ];
					y[ hole ] = y[ child ];
					z[ hole ] = z[ child ];
					k[ hole ] = k[ child ];
				} else
					break;
			}
			val[ hole ] = tmp;
			x[ hole ] = tmpX;
			y[ hole ] = tmpY;
			z[ hole ] = tmpZ;
			k[ hole ] = tmpK;
		} else if (minormax==MAXTREE) {
			for ( ; hole*2 <= currentSize; hole = child ) {
				child = hole*2;
				if (child != currentSize && val[child+1]>val[child])
					child++;
				if ( val[child]>tmp ) {
					val[ hole ] = val[ child ];
					x[ hole ] = x[ child ];
					y[ hole ] = y[ child ];
					z[ hole ] = z[ child ];
					k[ hole ] = k[ child ];
				} else
					break;
			}
			val[ hole ] = tmp;
			x[ hole ] = tmpX;
			y[ hole ] = tmpY;
			z[ hole ] = tmpZ;
			k[ hole ] = tmpK;
		} 

		return;
	}// removeFirstValue

	/**
	 * return the first value and its coordinates
	 */
	public final float getFirst() {
		return val[1];
	}
	public final short getFirstX() {
		return x[1];
	}
	public final short getFirstY() {
		return y[1];
	}
	public final short getFirstZ() {
		return z[1];
	}
	public final short getFirstK() {
		return k[1];
	}
	
	/**
	 * various utilities
	 */
	public final void print() {
		int i;
		int n;
			
		n = 2;
		for (i=1;i<currentSize;i++) {
			System.out.print("  "+(int)(100*val[i]));
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
		return (currentSize > 0);
	}
	
	/**
	 *  check the binary tree property
	 */
	public final int getCurrentSize() {
		return currentSize;
	}
	
	/**
	 *  check the binary tree property
	 */
	public final boolean isBinaryTree() {
		int i;
		int n;
		boolean isBinary=true;
				
		if (minormax==MINTREE) {
			for (i=2;i<currentSize;i++) {
				n = i/2;
				if ( val[n] > val[i] ) {
					// odd number
					isBinary = false;
					break;
				}
			}
		} else if (minormax==MAXTREE) {
			for (i=2;i<currentSize;i++) {
				n = i/2;
				if ( val[n] < val[i] ) {
					// odd number
					isBinary = false;
					break;
				}
			}
		}
		return isBinary;
	}//isBinaryTree

}

