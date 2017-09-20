package de.mpg.cbs.structures;

import java.io.*;
import java.util.*;

/**
 *
 *  This class handles octrees and perform simple operations.
 *	<p>
 *	Implementation loosely adapted from:
 *	"Simple and efficient traversal methods for quadtrees and octrees"
 *	S. Frisken and R. n. Perry, MERL.
 *
 *	@version    April 2006
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class OctreeCell {
	public	float		value;		// octree value
	public	int			x,y,z;		// coordinates
	public	int			level;		// level in the tree
	public	OctreeCell		parent;		// parent tree
	public 	OctreeCell[]	children;   // first children tree
	
	public OctreeCell(int x_, int y_, int z_, int lv_) {
		value = 0.0f;
		x = x_; y = y_; z = z_;
		level = lv_;
		parent = null;
		children = null;
	}
	
	public void finalize() {
		children = null;
	}
	
}
