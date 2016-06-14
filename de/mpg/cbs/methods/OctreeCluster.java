package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class handles octrees for distance-based clustering
 *	
 *	Implementation loosely adapted from:
 *	"Simple and efficient traversal methods for quadtrees and octrees"
 *	S. Frisken and R. n. Perry, MERL.
 *	
 *	(the elaborate binary operations used to speed up the system are not
 *	necessarily used here)
 *
 *	@version    March 2015
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class OctreeCluster {

	private 		Cell		root;
	private	static 	int 			LEVELS;
	private static	int 			ROOT_LEVEL;
	private static	float 			SCALE;	
	
	private static final byte		X=0;
	private static final byte		Y=1;
	private static final byte		Z=2;
	
	public static float[]		cellsizes;
	public static float[]		celllengths;
	
	private static final boolean	debug = false;
	private static final boolean	verbose = true;	
	
	public static class Cell {
		public	float		sum;		// octree sum value
		public	float		sqd;		// octree squared difference value
		public	int			x,y,z;		// coordinates
		public	int			level;		// level in the tree
		public	Cell		parent;		// parent tree
		public 	Cell[]	children;   // first children tree
		
		public Cell(int x_, int y_, int z_, int lv_) {
			sum = 0.0f;
			sqd = 0.0f;
			x = x_; y = y_; z = z_;
			level = lv_;
			parent = null;
			children = null;
		}
		
		public void finalize() {
			children = null;
		}
		
	}	
	
	public OctreeCluster(int lv_) {
		LEVELS = lv_;
		ROOT_LEVEL = lv_-1;
		SCALE = (float)FastMath.pow(2.0,ROOT_LEVEL);
		cellsizes = new float[LEVELS];
		float size = 1.0f;
		celllengths = new float[LEVELS];
		float length = 1.0f;
		for (int l=0;l<LEVELS;l++) {
			cellsizes[l] = size;
			size *= 8.0f;
			celllengths[l] = length;
			length *= 2.0f;
		}
	}
	
	/** 
	 *	traversal from a root cell to a leaf cell using x,y,z. 
	 *	<p>
	 *	The leaf cell is returned in cell
	 */
	public Cell traverseToLeaf(Cell init, int x, int y, int z) {
		Cell cell = init;
		while (cell.children!=null) {
			int cs = 1 << (cell.level-1);
			if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+","+cs+"]");
			int childIndex = 0;
			if (x>=cell.x+cs) childIndex += 1;
			if (y>=cell.y+cs) childIndex += 2;
			if (z>=cell.z+cs) childIndex += 4;
			cell = cell.children[childIndex];
		}
		if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+"]");
		return cell;
	}
	/** 
	 *	traversal from a root cell to a leaf cell using x,y,z. 
	 *	<p>
	 *	The leaf cell is returned in cell
	 */
	public Cell traverseToLeaf(int x, int y, int z) {
		return traverseToLeaf(root, x,y,z);
	}
	
	/** 
	 *	traversal from a root cell to a child cell using x,y,z down to a certain level. 
	 *	<p>
	 *	The child cell is returned in cell, and can be a leaf cell if the level is too low.
	 */
	public Cell traverseToLevel(Cell init, int x, int y, int z, int level) {
		Cell cell = init;
		int n = cell.level - level;
		while (cell.children!=null && n>0) {
			n--;
			int cs = 1 << (cell.level-1);
			if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+","+cs+"]");
			int childIndex = 0;
			if (x>=cell.x+cs) childIndex += 1;
			if (y>=cell.y+cs) childIndex += 2;
			if (z>=cell.z+cs) childIndex += 4;
			cell = cell.children[childIndex];
		}
		if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+"]");
		return cell;
	}
	
	/** 
	 *	locate the leaf cell containing point p. 
	 *	<p>
	 *	The point coordinates are in [0,1].
	 */
	public Cell locateCell(float[] p) {
		// determine the x,y,z codes for the point position
		int x = (int)(p[X]*SCALE);
		int y = (int)(p[Y]*SCALE);
		int z = (int)(p[Z]*SCALE);
	
		// follow the branching from root to leaf
		Cell cell = root;
		cell = traverseToLeaf(cell,x,y,z);
		return cell;
	}
	public Cell locateCell(float px, float py, float pz) {
		// determine the x,y,z codes for the point position
		int x = (int)(px*SCALE);
		int y = (int)(py*SCALE);
		int z = (int)(pz*SCALE);
	
		// follow the branching from root to leaf
		Cell cell = root;
		cell = traverseToLeaf(cell,x,y,z);
		return cell;
	}
	
	/** locate the left neighbor of same of larger size */
	public Cell getXminusNeighbor(Cell cell) {
		// if no left neighbor
		if (cell.x == 0) return null;
		else {
			// determine the smallest common ancestor
			Cell parentCell = cell;
			while (parentCell.x==cell.x) {
				parentCell = parentCell.parent;
				if (parentCell.parent==null) break;
			}
			
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x-1, cell.y, cell.z, cell.level);
			return parentCell;
		}
	}
	/** locate the right neighbor of same of larger size */
	public Cell getXplusNeighbor(Cell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+","+cellSize+"]");
		// if no right neighbor
		if  ( (cell.x+cellSize >= SCALE) || (cell.parent==null) ) {
			if (debug) System.out.print("no neighbor\n");
			return null;
		} else {
			// determine the smallest common ancestor
			Cell parentCell = cell.parent;
			if (debug) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			while (parentCell.x!=cell.x) {
				if (parentCell.parent==null) break;
				else parentCell = parentCell.parent;
				if (debug) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			}
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x+cellSize, cell.y, cell.z, cell.level);
			if (debug) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			if (debug) System.out.print("\n");
			return parentCell;
		}
	}
	/** locate the anterior neighbor of same of larger size */
	public Cell getYminusNeighbor(Cell cell) {
		// if no left neighbor
		if (cell.y == 0) return null;
		else {
			// determine the smallest common ancestor
			Cell parentCell = cell;
			while (parentCell.y==cell.y) {
				parentCell = parentCell.parent;
				if (parentCell.parent==null) break;
			}
			
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x, cell.y-1, cell.z, cell.level);
			return parentCell;
		}
	}
	/** locate the right neighbor of same of larger size */
	public Cell getYplusNeighbor(Cell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		// if no right neighbor
		if  ( (cell.y+cellSize >= SCALE) || (cell.parent==null) ) return null;
		else {
			// determine the smallest common ancestor
			Cell parentCell = cell.parent;
			while (parentCell.y!=cell.y) {
				if (parentCell.parent==null) break;
				else parentCell = parentCell.parent;
			}
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x, cell.y+cellSize, cell.z, cell.level);
			return parentCell;
		}
	}
	/** locate the left neighbor of same of larger size */
	public Cell getZminusNeighbor(Cell cell) {
		// if no left neighbor
		if (cell.z == 0) return null;
		else {
			// determine the smallest common ancestor
			Cell parentCell = cell;
			while (parentCell.z==cell.z) {
				parentCell = parentCell.parent;
				if (parentCell.parent==null) break;
			}
			
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x, cell.y, cell.z-1, cell.level);
			return parentCell;
		}
	}
	/** locate the right neighbor of same of larger size */
	public Cell getZplusNeighbor(Cell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		// if no right neighbor
		if  ( (cell.z+cellSize >= SCALE) || (cell.parent==null) ) return null;
		else {
			// determine the smallest common ancestor
			Cell parentCell = cell.parent;
			while (parentCell.z!=cell.z) {
				if (parentCell.parent==null) break;
				else parentCell = parentCell.parent;
			}
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x, cell.y, cell.z+cellSize, cell.level);
			return parentCell;
		}
	}
	/** 
	 *	create a new octree from the given image. 
	 *	<p>
	 *	It fills the tree with the image at finest scale 
	 *	(if the image is bigger, it gets cropped, if it's smaller it's padded with 0)
	 */
	public void createFromImage(float[] img, int nx, int ny, int nz) {
		// create the full tree recursively
		root = create(0,0,0, ROOT_LEVEL, 0, null);
		// set the image values at level 0
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,0,0,0);
					cellZ.sum = img[x+nx*y+nx*ny*z];
					int cs = 1 << (cellZ.level-1);
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level-1);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level-1);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		return;
	}
	/** 
	 *	create a new octree from the given image. 
	 *	<p>
	 *	It fills the tree with the image at the scale defined by level 
	 *	(lower levels of the tree are trimmed)
	 */
	public void createFromImage(float[] img, int nx, int ny, int nz, int level) {
		// create the full tree recursively
		root = create(0,0,0, ROOT_LEVEL, level, null);
		int incr = 1 << level;
		// set the image values at level 0
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x+=incr) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y+=incr) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z+=incr) {
					cellZ = traverseToLeaf(cellZ,0,0,0);
					cellZ.sum = img[x/incr + nx*y/incr + nx*ny*z/incr];
					int cs = 1 << (cellZ.level-1);
					if (z+incr>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level-1);
				if (y+incr>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level-1);
			if (x+incr>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		return;
	}
	/** for testing: works only with full octree */
	public void createFromImageInReverse(float[] img, int nx, int ny, int nz) {
		// create the full tree recursively
		root = create(0,0,0, ROOT_LEVEL, 0, null);
		// set the image values at level 0
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,(int)(SCALE-1),(int)(SCALE-1),(int)(SCALE-1));
		for (int x=(int)(SCALE-1);x>=0 && cellX!=null;x--) {
			cellY = cellX;
			for (int y=(int)(SCALE-1);y>=0 && cellY!=null;y--) {
				cellZ = cellY;
				for (int z=(int)(SCALE-1);z>=0 && cellZ!=null;z--) {
					if ( (x<nx) && (y<ny) && (z<nz) )
						cellZ.sum = img[x+nx*y+nx*ny*z];
					else
						cellZ.sum = 0.0f;
					cellZ = getZminusNeighbor(cellZ);
				}
				cellY = getYminusNeighbor(cellY);
			}
			cellX = getXminusNeighbor(cellX);
		}
		return;
	}
	public void createFromImageRootToLeaf(float[] img, int nx, int ny, int nz) {
		// create the full tree recursively
		root = create(0,0,0, ROOT_LEVEL, 0, null);
		// set the image values at level 0
		Cell cell;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if ( (x<SCALE) && (y<SCALE) && (z<SCALE) ) {
				if (debug) System.out.print("("+x+","+y+","+z+")->");
				cell = root;
				cell = traverseToLeaf(cell,x,y,z);
				if (debug) System.out.print("\n");
				cell.sum = img[x+nx*y+nx*ny*z];
			}
		}
		return;
	}
	
	/** 
	 *	create an empty octree fully unfolded
	 *	<p>
	 *	All values are 0
	 */
	public void createBlank() { root = create(0,0,0, ROOT_LEVEL, 0, null); }
	
	/** recursive octree building.
	 *	<p>
	 *	The octree is fully propagated from current location and level to minlevel.
	 */
	public Cell create(int x, int y, int z, int lv, int minlv, Cell parent) {
		Cell cell = new Cell(x,y,z,lv);
		cell.parent = parent;
		if (lv>minlv) {
			lv--;
			// cell size
			int cs = 1 << lv;
			cell.children = new Cell[8];
			// create the tree along proper pattern
			cell.children[0] = create(x, 	  y, 	z,    lv, minlv, cell);
			
			cell.children[1] = create(x+cs, y, 	z,    lv, minlv, cell);
			cell.children[2] = create(x,	  y+cs, z,    lv, minlv, cell);
			cell.children[4] = create(x, 	  y, 	z+cs, lv, minlv, cell);
			
			cell.children[3] = create(x+cs, y+cs, z,    lv, minlv, cell);
			cell.children[5] = create(x+cs, y, 	z+cs, lv, minlv, cell);
			cell.children[6] = create(x, 	  y+cs, z+cs, lv, minlv, cell);
			
			cell.children[7] = create(x+cs, y+cs, z+cs, lv, minlv, cell);
		}
		return cell;
	}
	
	/** compute the necessary octree level for storing an image of nx,ny,nz dimensions */
	public static int findMinimumLevel(int nx, int ny, int nz) {
		int level = 1;
		int size = 1;
		while ( (size<nx) || (size<ny) || (size<nz) ) {
			size = 2*size;
			level++;
		}
		return level;
	}
	
	/** bring the octree back into an image.
	 *	<p>
	 *	Assumes no particular structure for the tree: if might be slow.
	 *	(however it sweeps regularly rather than going back to the root)
	 */
	public float[] exportMeanToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					int cs = 1 << (cellZ.level);
					img[x+nx*y+nx*ny*z] = cellZ.sum/cellsizes[cellZ.level];
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	/** bring the octree back into an image.
	 *	<p>
	 *	Assumes no particular structure for the tree: if might be slow.
	 *	(however it sweeps regularly rather than going back to the root)
	 */
	public float[] exportSumToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					int cs = 1 << (cellZ.level);
					img[x+nx*y+nx*ny*z] = cellZ.sum;
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	/** bring the octree back into an image.
	 *	<p>
	 *	Assumes no particular structure for the tree: if might be slow.
	 *	(however it sweeps regularly rather than going back to the root)
	 */
	public float[] exportStdevToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					int cs = 1 << (cellZ.level);
					//img[x+nx*y+nx*ny*z] = cellZ.sqd;
					if (cellZ.level>0) img[x+nx*y+nx*ny*z] = (float)FastMath.sqrt( cellZ.sqd/(cellsizes[cellZ.level]-1) );
					else img[x+nx*y+nx*ny*z] = 0.0f;
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	/** prints the octree scale into an image.
	 *	<p>
	 *	Assumes no particular structure for the tree: if might be slow.
	 *	(however it sweeps regularly rather than going back to the root)
	 */
	public float[] exportScaleToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					img[x+nx*y+nx*ny*z] = cellZ.level;
					int cs = 1 << (cellZ.level);
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	/** prints the octree scale into an image.
	 *	<p>
	 *	Assumes no particular structure for the tree: if might be slow.
	 *	(however it sweeps regularly rather than going back to the root)
	 */
	public int[] exportLabelingToImage(int nx, int ny, int nz) {
		int[] img = new int[nx*ny*nz];
		int label=0;
		int EMPTY = -1;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			img[x+nx*y+nx*ny*z] = EMPTY;
		}
		
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					// only increase the label count if unlabelled cell
					if (img[x+nx*y+nx*ny*z]==EMPTY) {
						label++;
						// label all elements of the cell
						for (int i=0;i<celllengths[cellZ.level];i++) {
							for (int j=0;j<celllengths[cellZ.level];j++) { 
								for (int k=0;k<celllengths[cellZ.level];k++) {
									int id = x+i+nx*(y+j)+nx*ny*(z+k);
									if (img[id]==EMPTY) img[id] = label;
								}
							}
						}
					}
					int cs = 1 << (cellZ.level);
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	public int[] exportMaskedLabelingToImage(boolean[] mask, int nx, int ny, int nz) {
		int[] img = new int[nx*ny*nz];
		int label=0;
		int EMPTY = -1;
		int MASKED = 0;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			img[x+nx*y+nx*ny*z] = EMPTY;
		}
		
		// apparently this is very very slow :(
		
		Cell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					// only increase the label count if unlabelled cell
					if (img[x+nx*y+nx*ny*z]==EMPTY) {
						label++;
						// label all elements of the cell
						boolean ismasked = true;
						for (int i=0;i<celllengths[cellZ.level];i++) if (x+i<nx) {
							for (int j=0;j<celllengths[cellZ.level];j++) if (y+j<ny) { 
								for (int k=0;k<celllengths[cellZ.level];k++) if (z+k<nz) {
									int id = x+i+nx*(y+j)+nx*ny*(z+k);
									if (img[id]==EMPTY) img[id] = label;
									if (mask[id]) ismasked = false;
								}
							}
						}
						// if all elements are outside the mask : set to zero
						if (ismasked) {
							for (int i=0;i<celllengths[cellZ.level];i++) if (x+i<nx) {
								for (int j=0;j<celllengths[cellZ.level];j++) if (y+j<ny) { 
									for (int k=0;k<celllengths[cellZ.level];k++) if (z+k<nz) {
										int id = x+i+nx*(y+j)+nx*ny*(z+k);
										img[id] = MASKED;
									}
								}
							}
							label--;
						}
					}
					int cs = 1 << (cellZ.level);
					if (z+1>=cellZ.z+cs) cellZ = getZplusNeighbor(cellZ);
				}
				int cs = 1 << (cellY.level);
				if (y+1>=cellY.y+cs) cellY = getYplusNeighbor(cellY);
			}
			int cs = 1 << (cellX.level);
			if (x+1>=cellX.x+cs) cellX = getXplusNeighbor(cellX);
		}
		
		return img;
	}
	
	public int createMaskedLabelingFromScale(int[] labeling, boolean[] mask, int nx, int ny, int nz) {
		System.out.println("getting the octree scale"); System.out.flush();
		
		float[] scale = exportScaleToImage(nx, ny, nz);
		
		int label=0;
		int EMPTY = -1;
		int MASKED = 0;
		
		System.out.println("creating labels"); System.out.flush();
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			labeling[x+nx*y+nx*ny*z] = EMPTY;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (labeling[x+nx*y+nx*ny*z]==EMPTY) {
				label++;
				// label all elements of the cell
				boolean ismasked = true;
				for (int i=0;i<celllengths[(int)scale[x+nx*y+nx*ny*z]];i++) if (x+i<nx) {
					for (int j=0;j<celllengths[(int)scale[x+nx*y+nx*ny*z]];j++) if (y+j<ny) { 
						for (int k=0;k<celllengths[(int)scale[x+nx*y+nx*ny*z]];k++) if (z+k<nz) {
							int id = x+i+nx*(y+j)+nx*ny*(z+k);
							if (labeling[id]==EMPTY) labeling[id] = label;
							if (mask[id]) ismasked = false;
						}
					}
				}
				// if all elements are outside the mask : set to zero
				if (ismasked) {
					for (int i=0;i<celllengths[(int)scale[x+nx*y+nx*ny*z]];i++) if (x+i<nx) {
						for (int j=0;j<celllengths[(int)scale[x+nx*y+nx*ny*z]];j++) if (y+j<ny) { 
							for (int k=0;k<celllengths[(int)scale[x+nx*y+nx*ny*z]];k++) if (z+k<nz) {
								int id = x+i+nx*(y+j)+nx*ny*(z+k);
								labeling[id] = MASKED;
							}
						}
					}
					label--;
				}
			}
		}
		
		return label;
	}
	
	/** bring the octree back into an image.
	 *	<p>
	 *	Go back to the root every time (slow, for debug)
	 */
	public float[] exportToImageRootToLeaf(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if ( (x<SCALE) && (y<SCALE) && (z<SCALE) ) {
				Cell cell = traverseToLeaf(root,x,y,z);
				int cs = 1 << (cell.level);
				img[x+nx*y+nx*ny*z] = cell.sum/(cs*cs*cs);
			} else {
				img[x+nx*y+nx*ny*z] = 0.0f;
			}
		}
		
		return img;
	}

	/** display the octree structure.
	 *	<p>
	 *	Starts from cell and explores the entire tree below it.
	 */
	private String display(Cell cell) {
		if (cell.children==null) {
			return "("+cell.x+","+cell.y+","+cell.z+","+cell.level+")"+cell.sum+"\n";
		} else {
			String info = "("+cell.x+","+cell.y+","+cell.z+","+cell.level+")";
			//String info = "";
			for (int n=0;n<8;n++) {
				info += display(cell.children[n]);
			}
			return info;
		}
	}
	/** display the octree structure from the root node.
	 */
	public String display() {
		String info = "Octree: ";
		info += display(root);
		return info;
	}
	
	/** prune the octree to group similar values.
	 *	<p>
	 *	Remove child cells of value less than dist apart from their mean.
	 */
	public void pruneToDistance(float dist) {
		float mean;
		for (int l=0;l<ROOT_LEVEL;l++) {
			for (int x=0;x<SCALE;x++) for (int y=0;y<SCALE;y++) for (int z=0;z<SCALE;z++) {
				Cell cell = (traverseToLeaf(root,x,y,z)).parent;
				if (cell!=null && cell.level==l+1) {
					// check if leaf parent for all
					boolean allleaf = true;
					for (int n=0;n<8;n++) if (cell.children[n].children!=null) allleaf = false;
					
					if (allleaf) {
						// compute mean cell value
						mean = 0.0f;
						for (int n=0;n<8;n++) {
							mean += (cell.children[n]).sum/8.0f;
						}
						// compare to value
						boolean merge = true;
						for (int n=0;n<8;n++) {
							if (Numerics.abs(mean-cell.children[n].sum) > dist) {
								merge = false;
								break;
							}
						}
						// merge
						if (merge) {
							cell.sum = mean;
							if (debug) System.out.print("merge:["+cell.x+","+cell.y+","+cell.z+","+cell.level+":"+cell.sum+"]\n");
							for (int n=0;n<8;n++) cell.children[n] = null;
							cell.children = null;
						}
					}
				}
			}
		}
	}

	/** prune the octree to group similar values.
	 *	<p>
	 *	Remove child cells based on Jensen-Shannon divergence
	 */
	public void pruneToJensenShannonDistance(float sigma2, float threshold) {
		for (int l=0;l<ROOT_LEVEL;l++) {
			for (int x=0;x<SCALE;x++) for (int y=0;y<SCALE;y++) for (int z=0;z<SCALE;z++) {
				//if (verbose) System.out.print(".");
				Cell cell = (traverseToLeaf(root,x,y,z)).parent;
				if (cell!=null && cell.level==l+1) {
					//if (verbose) System.out.print("-");
				
					// check if leaf parent for all
					boolean allleaf = true;
					for (int n=0;n<8;n++) if (cell.children[n].children!=null) allleaf = false;
					
					if (allleaf) {
						//if (verbose) System.out.print(">");
				
						// compute the joint JS divergence : for each connected leaf pair!
						// assume 6C to start: 12 pairs
						boolean merge = true;
						//if (verbose) System.out.print("+");
						//System.out.flush();
						
						float js = jensenShannonDivergence(cell.children[0],cell.children[1],sigma2);
						js = Numerics.min(js, jensenShannonDivergence(cell.children[0],cell.children[2],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[0],cell.children[4],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[1],cell.children[3],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[1],cell.children[5],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[2],cell.children[3],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[2],cell.children[6],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[3],cell.children[7],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[4],cell.children[5],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[4],cell.children[6],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[5],cell.children[7],sigma2));
						js = Numerics.min(js, jensenShannonDivergence(cell.children[6],cell.children[7],sigma2));
						
						if (js<threshold) merge = false;
						//else if (verbose) System.out.print("js: "+js+" ");
							
						// merging: add the sum, sqd incrementally for each leaf with Chan's formula
						if (merge) {
							// add progressively all children
							float cs = cellsizes[cell.level-1];
							cell.sum = cell.children[0].sum;
							cell.sqd = cell.children[0].sqd;
							for (int n=1;n<8;n++) {
								// update first the square distance to keep the merged sum and the added sum separate
								cell.sqd += cell.children[n].sqd + n*cs*cs/(n*cs+cs)*Numerics.square(cell.sum/(n*cs)-cell.children[n].sum/cs);
								cell.sum += cell.children[n].sum;
							}
							
							if (debug && cell.sum>0) System.out.print("merge:["+cell.x+","+cell.y+","+cell.z+","+cell.level+":"+cell.sum+", "+cell.sqd+"]\n");
							for (int n=0;n<8;n++) cell.children[n] = null;
							cell.children = null;
						}
					}
				}
			}
		}
	}
	
	float jensenShannonDivergence(Cell cellA, Cell cellB, float sigma2) {
		// cell sizes
		float csA = cellsizes[cellA.level];
		float csB = cellsizes[cellB.level];
		//System.out.print("sizes: "+csA+","+csB+"\n");
		//System.out.flush();
		// joint variance
		double vA = Numerics.max(cellA.sqd/csA, sigma2);
		double vB = Numerics.max(cellB.sqd/csB, sigma2);
		double vAB = 0.25*( vA + vB );
		//System.out.print("variances: "+vA+","+vB+","+vAB+"\n");
		//System.out.flush();
		
		double js = 0.25*Numerics.square(cellA.sum/csA-cellB.sum/csB)/vAB - 0.5*FastMath.log(vA/vAB*vB/vAB) + FastMath.log(2.0);
		//System.out.print("js dist: "+js+"\n");
		//System.out.flush();
		return (float)FastMath.exp(-0.5*js);
		//return (float)FastMath.exp(-0.5*Numerics.square(cellA.sum-cellB.sum)/sigma2);
	}
}
