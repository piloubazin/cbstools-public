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

public class BinaryHeap2D {
	
	private float[] val;
	private int[] 	id;
	private byte[] 	state;
	
	private int 		currentSize;
	private int			capacity;
	private int 		minormax;
	
	public static final	int MINTREE = -1;
	public static final	int MAXTREE = 1;
	
	public BinaryHeap2D(int Nsize, int type) {
		currentSize = 0;
		capacity = Nsize;
		minormax = type;
		val = new float[capacity+1];
		id = new int[capacity+1];
		state = new byte[capacity+1];
		if (minormax==MINTREE)
			val[0] = -1e12f;
		else if (minormax==MAXTREE)
			val[0] = 1e12f;
	}
	
	public BinaryHeap2D(int Nsize, int Nincrease, int type) {
		currentSize = 0;
		capacity = Nsize;
		minormax = type;
		val = new float[capacity+1];
		id = new int[capacity+1];
		state = new byte[capacity+1];
		if (minormax==MINTREE)
			val[0] = -1e12f;
		else if (minormax==MAXTREE)
			val[0] = 1e12f;
		
		capacity = Nincrease;
	}
	
	public void finalize() {
		val = null;
		id = null;
		state = null;
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
	public final void addValue(float val_, int id_, byte state_) {
		// check for size
		if  (currentSize == val.length - 1) {
			float[] oldVal = val;
			int[] 	oldId = id;
			byte[] 	oldState = state;
			val = new float[currentSize+capacity];
			id = new int[currentSize+capacity];
			state = new byte[currentSize+capacity];
			for (int i=0;i<oldVal.length;i++) {
				val[i] = oldVal[i];
				id[i] = oldId[i];
				state[i] = oldState[i];
			}
			oldVal = null; oldId = null; oldState = null;
		}
		// insert new  point into the proper location		
		int hole = ++currentSize;
		
		if (minormax==MINTREE) {
			for ( ; val_ < val[ hole/2 ]; hole /= 2 ) {
				val[hole] = val[hole/2];
				id[hole] = id[hole/2];
				state[hole] = state[hole/2];
			}
			val[hole] = val_;
			id[hole] = id_;
			state[hole] = state_;
		} else if (minormax==MAXTREE) {
			for ( ; val_ > val[ hole/2 ]; hole /= 2 ) {
				val[hole] = val[hole/2];
				id[hole] = id[hole/2];
				state[hole] = state[hole/2];
			}
			val[hole] = val_;
			id[hole] = id_;
			state[hole] = state_;
		}
		
		return;
	}//addValue
	
	/**
	 *  remove the first value from the tree
	 */
	public final void removeFirst() {
		int hole = 1;
		
		val[hole] = val[currentSize];
		id[hole] = id[currentSize];
		state[hole] = state[currentSize];
		currentSize--;
		
		int child;
		float tmp = val[hole];
		int tmpId = id[hole];
		byte tmpState = state[hole];
		
		if (minormax==MINTREE) {
			for ( ; hole*2 <= currentSize; hole = child ) {
				child = hole*2;
				if (child != currentSize && val[child+1]<val[child])
					child++;
				if ( val[child]<tmp ) {
					val[ hole ] = val[ child ];
					id[ hole ] = id[ child ];
					state[ hole ] = state[ child ];
				} else
					break;
			}
			val[ hole ] = tmp;
			id[ hole ] = tmpId;
			state[ hole ] = tmpState;
		} else if (minormax==MAXTREE) {
			for ( ; hole*2 <= currentSize; hole = child ) {
				child = hole*2;
				if (child != currentSize && val[child+1]>val[child])
					child++;
				if ( val[child]>tmp ) {
					val[ hole ] = val[ child ];
					id[ hole ] = id[ child ];
					state[ hole ] = state[ child ];
				} else
					break;
			}
			val[ hole ] = tmp;
			id[ hole ] = tmpId;
			state[ hole ] = tmpState;
		} 

		return;
	}// removeFirstValue

	/**
	 * return the first value and its coordinates
	 */
	public final float getFirst() {
		return val[1];
	}
	public final int getFirstId() {
		return id[1];
	}
	public final byte getFirstState() {
		return state[1];
	}
	
	/**
	 *  return any value and coordinates (note that the ranking is arbitrary)
	 */
	public final float get(int n) {
		return val[n];
	}
	public final float getId(int n) {
		return id[n];
	}
	public final float getState(int n) {
		return state[n];
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

