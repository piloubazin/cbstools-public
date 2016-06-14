package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;

import de.mpg.cbs.utilities.*;

/**
 *
 *  This class handles octrees and perform simple operations
 *	
 *	Implementation loosely adapted from:
 *	"Simple and efficient traversal methods for quadtrees and octrees"
 *	S. Frisken and R. n. Perry, MERL.
 *	
 *	(the elaborate binary operations used to speed up the system are not
 *	necessarily used here)
 *
 *	@version    April 2006
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Octree {

	private 		OctreeCell		root;
	private	static 	int 			LEVELS;
	private static	int 			ROOT_LEVEL;
	private static	float 			SCALE;	
	
	private static final byte		X=0;
	private static final byte		Y=1;
	private static final byte		Z=2;

	private static final boolean	debug = false;
	private static final boolean	verbose = false;	
	
	public Octree(int lv_) {
		LEVELS = lv_;
		ROOT_LEVEL = lv_-1;
		SCALE = (float)Math.pow(2.0,ROOT_LEVEL);
	}
	
	/** 
	 *	traversal from a root cell to a leaf cell using x,y,z. 
	 *	<p>
	 *	The leaf cell is returned in cell
	 */
	public OctreeCell traverseToLeaf(OctreeCell init, int x, int y, int z) {
		OctreeCell cell = init;
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
	 *	traversal from a root cell to a child cell using x,y,z down to a certain level. 
	 *	<p>
	 *	The child cell is returned in cell, and can be a leaf cell if the level is too low.
	 */
	public OctreeCell traverseToLevel(OctreeCell init, int x, int y, int z, int level) {
		OctreeCell cell = init;
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
	public OctreeCell locateCell(float[] p) {
		// determine the x,y,z codes for the point position
		int x = (int)(p[X]*SCALE);
		int y = (int)(p[Y]*SCALE);
		int z = (int)(p[Z]*SCALE);
	
		// follow the branching from root to leaf
		OctreeCell cell = root;
		cell = traverseToLeaf(cell,x,y,z);
		return cell;
	}
	public OctreeCell locateCell(float px, float py, float pz) {
		// determine the x,y,z codes for the point position
		int x = (int)(px*SCALE);
		int y = (int)(py*SCALE);
		int z = (int)(pz*SCALE);
	
		// follow the branching from root to leaf
		OctreeCell cell = root;
		cell = traverseToLeaf(cell,x,y,z);
		return cell;
	}
	
	/** locate the left neighbor of same of larger size */
	public OctreeCell getXminusNeighbor(OctreeCell cell) {
		// if no left neighbor
		if (cell.x == 0) return null;
		else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell;
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
	public OctreeCell getXplusNeighbor(OctreeCell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		if (debug) System.out.print("["+cell.x+","+cell.y+","+cell.z+","+cellSize+"]");
		// if no right neighbor
		if  ( (cell.x+cellSize >= SCALE) || (cell.parent==null) ) {
			if (debug) System.out.print("no neighbor\n");
			return null;
		} else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell.parent;
			if (debug) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			while (parentCell.x!=cell.x) {
				if (parentCell.parent==null) break;
				else parentCell = parentCell.parent;
				if (verbose) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			}
			// start from smallest ancestor and follow branching
			parentCell = traverseToLevel(parentCell, cell.x+cellSize, cell.y, cell.z, cell.level);
			if (debug) System.out.print("["+parentCell.x+","+parentCell.y+","+parentCell.z+"]");
			if (debug) System.out.print("\n");
			return parentCell;
		}
	}
	/** locate the anterior neighbor of same of larger size */
	public OctreeCell getYminusNeighbor(OctreeCell cell) {
		// if no left neighbor
		if (cell.y == 0) return null;
		else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell;
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
	public OctreeCell getYplusNeighbor(OctreeCell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		// if no right neighbor
		if  ( (cell.y+cellSize >= SCALE) || (cell.parent==null) ) return null;
		else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell.parent;
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
	public OctreeCell getZminusNeighbor(OctreeCell cell) {
		// if no left neighbor
		if (cell.z == 0) return null;
		else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell;
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
	public OctreeCell getZplusNeighbor(OctreeCell cell) {
		// get location of smallest possible right neighbor
		int cellSize = 1 << cell.level;
		// if no right neighbor
		if  ( (cell.z+cellSize >= SCALE) || (cell.parent==null) ) return null;
		else {
			// determine the smallest common ancestor
			OctreeCell parentCell = cell.parent;
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
		OctreeCell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,0,0,0);
					cellZ.value = img[x+nx*y+nx*ny*z];
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
		OctreeCell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x+=incr) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y+=incr) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z+=incr) {
					cellZ = traverseToLeaf(cellZ,0,0,0);
					cellZ.value = img[x/incr + nx*y/incr + nx*ny*z/incr];
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
		OctreeCell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,(int)(SCALE-1),(int)(SCALE-1),(int)(SCALE-1));
		for (int x=(int)(SCALE-1);x>=0 && cellX!=null;x--) {
			cellY = cellX;
			for (int y=(int)(SCALE-1);y>=0 && cellY!=null;y--) {
				cellZ = cellY;
				for (int z=(int)(SCALE-1);z>=0 && cellZ!=null;z--) {
					if ( (x<nx) && (y<ny) && (z<nz) )
						cellZ.value = img[x+nx*y+nx*ny*z];
					else
						cellZ.value = 0.0f;
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
		OctreeCell cell;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if ( (x<SCALE) && (y<SCALE) && (z<SCALE) ) {
				if (debug) System.out.print("("+x+","+y+","+z+")->");
				cell = root;
				cell = traverseToLeaf(cell,x,y,z);
				if (debug) System.out.print("\n");
				cell.value = img[x+nx*y+nx*ny*z];
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
	private OctreeCell create(int x, int y, int z, int lv, int minlv, OctreeCell parent) {
		OctreeCell cell = new OctreeCell(x,y,z,lv);
		cell.parent = parent;
		if (lv>minlv) {
			lv--;
			// cell size
			int cs = 1 << lv;
			cell.children = new OctreeCell[8];
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
	public float[] exportToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		OctreeCell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<nx && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<ny && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<nz && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,x,y,z);
					img[x+nx*y+nx*ny*z] = cellZ.value;
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
	public float[] exportScaleToImage(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		OctreeCell cellX,cellY,cellZ;
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
	
	/** bring the octree back into an image.
	 *	<p>
	 *	Go back to the root every time (slow, for debug)
	 */
	public float[] exportToImageRootToLeaf(int nx, int ny, int nz) {
		float[] img = new float[nx*ny*nz];
	
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if ( (x<SCALE) && (y<SCALE) && (z<SCALE) ) {
				OctreeCell cell = traverseToLeaf(root,x,y,z);
				img[x+nx*y+nx*ny*z] = cell.value;
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
	private String display(OctreeCell cell) {
		if (cell.children==null) {
			return "("+cell.x+","+cell.y+","+cell.z+","+cell.level+")"+cell.value+"\n";
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
		/*
		OctreeCell cellX,cellY,cellZ;
		cellX = traverseToLeaf(root,0,0,0);
		for (int x=0;(x<SCALE && cellX!=null);x++) {
			cellY = cellX;
			for (int y=0;(y<SCALE && cellY!=null);y++) {
				cellZ = cellY;
				for (int z=0;(z<SCALE && cellZ!=null);z++) {
					cellZ = traverseToLeaf(cellZ,0,0,0);
					OctreeCell cell = cellZ.parent;
					boolean merge = false;
					if (cell!=null) {
						// compute mean cell value
						mean = 0.0f;
						for (int n=0;n<8;n++) {
							mean += (cell.children[n]).value/8.0f;
						}
						// compare to value
						merge = true;
						for (int n=0;n<8;n++) {
							OctreeCell child = cell.children[n];
							if ( (child.children!=null) || (Numerics.abs(mean -child.value)>dist) ) {
								merge = false;
								break;
							}
						}
						// merge
						if (merge) {
							cell.value = mean;
							if (debug) System.out.print("merge:["+cell.x+","+cell.y+","+cell.z+","+cell.level+":"+cell.value+"]\n");
							for (int n=0;n<8;n++) cell.children[n] = null;
							cell.children = null;
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
		*/
		for (int l=0;l<ROOT_LEVEL;l++) {
			for (int x=0;x<SCALE;x++) for (int y=0;y<SCALE;y++) for (int z=0;z<SCALE;z++) {
				OctreeCell cell = (traverseToLeaf(root,x,y,z)).parent;
				if (cell!=null && cell.level==l+1) {
					// check if leaf parent for all
					boolean allleaf = true;
					for (int n=0;n<8;n++) if (cell.children[n].children!=null) allleaf = false;
					
					if (allleaf) {
						// compute mean cell value
						mean = 0.0f;
						for (int n=0;n<8;n++) {
							mean += (cell.children[n]).value/8.0f;
						}
						// compare to value
						boolean merge = true;
						for (int n=0;n<8;n++) {
							if (Numerics.abs(mean-cell.children[n].value) > dist) {
								merge = false;
								break;
							}
						}
						// merge
						if (merge) {
							cell.value = mean;
							if (debug) System.out.print("merge:["+cell.x+","+cell.y+","+cell.z+","+cell.level+":"+cell.value+"]\n");
							for (int n=0;n<8;n++) cell.children[n] = null;
							cell.children = null;
						}
					}
				}
			}
		}
	}
}
