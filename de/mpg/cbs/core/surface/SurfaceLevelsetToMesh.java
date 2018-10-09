package de.mpg.cbs.core.surface;


import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;

/**
 * Generates a topologically consistent iso-surface for a levelel set with a
 * specified digital connectivity rule.
 * 
 * adapted from the Connectivity-consistent Marching Cubes code
 * from Xiao Han and Blake Lucas
 * 
 */
public class SurfaceLevelsetToMesh {
    
    // image containers
    private float[] lvlImage;
    private int connectivity;
    private float level = 0.0f;
    private boolean inclusive = true;
    
    private float[] pointList;
    private int[] triangleList;
    
	private int nx=-1,ny=-1,nz=-1,nxyz=-1;
	private float rx, ry, rz;

	// create inputs
	public final void setLevelsetImage(float[] val) { lvlImage = val; }
	public final void setConnectivity(String val) { 
	         if (val.equals("6/18")) connectivity = CONNECTIVITY_6_18; 
	    else if (val.equals("18/6")) connectivity = CONNECTIVITY_18_6; 
	    else if (val.equals("6/26")) connectivity = CONNECTIVITY_6_26; 
	    else if (val.equals("26/6")) connectivity = CONNECTIVITY_26_6; 
	}
	public final void setZeroLevel(float val) { level = val; }
	public final void setInclusive(boolean val) { inclusive = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Surfaces"; }
	public final String getLabel() { return "Level Set To Mesh"; }
	public final String getName() { return "LevelSetToMesh"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Convert a levelset to a surface mesh"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	public final float[] getPointList() { return pointList; }
    public final int[] getTriangleList() { return triangleList; }
    

    // local data types
	private static final int[][] Perms = { { 0, 1, 2, 3, 4, 5, 6, 7 },
			{ 2, 0, 3, 1, 6, 4, 7, 5 }, { 3, 2, 1, 0, 7, 6, 5, 4 },
			{ 1, 3, 0, 2, 5, 7, 4, 6 }, { 4, 0, 6, 2, 5, 1, 7, 3 },
			{ 5, 4, 7, 6, 1, 0, 3, 2 }, { 1, 5, 3, 7, 0, 4, 2, 6 },
			{ 4, 5, 0, 1, 6, 7, 2, 3 }, { 6, 7, 4, 5, 2, 3, 0, 1 },
			{ 2, 3, 6, 7, 0, 1, 4, 5 }, { 1, 0, 5, 4, 3, 2, 7, 6 },
			{ 7, 6, 3, 2, 5, 4, 1, 0 }, { 4, 6, 5, 7, 0, 2, 1, 3 },
			{ 7, 5, 6, 4, 3, 1, 2, 0 }, { 2, 6, 0, 4, 3, 7, 1, 5 },
			{ 7, 3, 5, 1, 6, 2, 4, 0 }, { 0, 4, 1, 5, 2, 6, 3, 7 },
			{ 0, 2, 4, 6, 1, 3, 5, 7 }, { 3, 1, 7, 5, 2, 0, 6, 4 },
			{ 5, 1, 4, 0, 7, 3, 6, 2 }, { 3, 7, 2, 6, 1, 5, 0, 4 },
			{ 6, 4, 2, 0, 7, 5, 3, 1 }, { 5, 7, 1, 3, 4, 6, 0, 2 },
			{ 6, 2, 7, 3, 4, 0, 5, 1 } };

	private static final int HASHBIT = 5;

	private static final int HASHSIZE = (1 << (3 * HASHBIT)); /*
															 * hash table size
															 * (32768)
															 */

	private static final int MASK = ((1 << HASHBIT) - 1);

	public int HASH(int i, int j, int k) {
		return ((((((k) & MASK) << HASHBIT) | ((j) & MASK)) << HASHBIT) | ((i) & MASK));
	}

	public static class Corner {
		public int i, j, k;

		public double value;
	}

	public static class EdgeIndex {
		public int i1, j1, k1;

		public int i2, j2, k2;

		public int vid;

		public EdgeIndex next;

		public EdgeIndex(int i1, int j1, int k1, int i2, int j2, int k2,
				int vid, EdgeIndex next) {
			this.i1 = i1;
			this.j1 = j1;
			this.k1 = k1;
			this.i2 = i2;
			this.j2 = j2;
			this.k2 = k2;
			this.vid = vid;
		}
	}

	public static class TriangleIndex {
		public int vid1, vid2, vid3;
	}
	
	
	public static class KdPoint3 {
	    public float x;
	    public float y;
	    public float z;
	}

	protected static final byte COUNTER_CLOCKWISE = 9;
	protected static final byte CLOCKWISE = 3;
	protected byte direction = COUNTER_CLOCKWISE;

	protected static final int CONNECTIVITY_18_6 = 186;
	protected static final int CONNECTIVITY_6_18 = 618;
	protected static final int CONNECTIVITY_26_6 = 266;
	protected static final int CONNECTIVITY_6_26 = 626;
	
	public void setDirection(byte direction) {
		this.direction = direction;
	}

	/* local data containers */
	private float[] data;
	
	private EdgeIndex[] edgelist;

	private ArrayList<KdPoint3> verts;

	private ArrayList<TriangleIndex> trilist;

	
	/* Set vertex id for edge */
	void setedge(EdgeIndex[] table, int i1, int j1, int k1, int i2, int j2, int k2, int vid) {
		int index;
		if (i1 > i2 || (i1 == i2 && (j1 > j2 || (j1 == j2 && k1 > k2)))) {

			int t = i1;
			i1 = i2;
			i2 = t;

			t = j1;
			j1 = j2;
			j2 = t;

			t = k1;
			k1 = k2;
			k2 = t;
		}
		index = HASH(i1, j1, k1) + HASH(i2, j2, k2);
		EdgeIndex new1 = new EdgeIndex(i1, j1, k1, i2, j2, k2, vid,
				table[index]);
		new1.next = table[index];
		table[index] = new1;
	}

	/* Return vertex id for edge; return -1 if not set */
	/*
	 * Each vertex is distinguished by the edge it lies on; thus several
	 * vertices with different ID may have exactly the same coordinates if they
	 * lie on a same grid point. Hence, function value equals to iso-value is
	 * intolerable! Although we can remove triangle with two or more vertices of
	 * same coords, the comparison (of double-type values) will take time.
	 */
	private int getedge(EdgeIndex[] table, int i1, int j1, int k1, int i2, int j2, int k2)	{
		EdgeIndex q;
		if (i1 > i2 || (i1 == i2 && (j1 > j2 || (j1 == j2 && k1 > k2)))) {
			int t = i1;
			i1 = i2;
			i2 = t;
			t = j1;
			j1 = j2;
			j2 = t;
			t = k1;
			k1 = k2;
			k2 = t;
		}
		;

		q = table[HASH(i1, j1, k1) + HASH(i2, j2, k2)];
		for (; q != null; q = q.next) {
			if (q.i1 == i1 && q.j1 == j1 && q.k1 == k1 && q.i2 == i2
					&& q.j2 == j2 && q.k2 == k2) {
				return q.vid;
			}
		}
		return -1;
	}

	private int vertid(Corner c1, Corner c2)
	/*
	 * Return index for vertex on edge: c1.value and c2.value are presumed of
	 * different sign return saved index if any; else compute vertex and save
	 */
	{
		KdPoint3 v = new KdPoint3();
		KdPoint3 a = new KdPoint3(), b = new KdPoint3();
		int vid = getedge(edgelist, c1.i, c1.j, c1.k, c2.i, c2.j, c2.k);
		// System.out.printf("getid %d %d %d %d\n",vid,c1.i,c1.j,c1.k);
		if (vid != -1) {
			return vid; /* previously computed */
		}
		a.x = c1.i;
		a.y = c1.j;
		a.z = c1.k;

		b.x = c2.i;
		b.y = c2.j;
		b.z = c2.k;
		interpolate(a, b, c1.value, c2.value, v); /* position */
		// System.out.format("(%6.2f,%6.2f,%6.2f) %6.2f (%6.2f,%6.2f,%6.2f) %6.2f (%6.2f,%6.2f,%6.2f)\n",a.x,a.y,a.z,c1.value,b.x,b.y,b.z,c2.value,v.x,v.y,v.z);

		addtovertices(verts, v); /* save vertex */
		vid = verts.size() - 1;
		setedge(edgelist, c1.i, c1.j, c1.k, c2.i, c2.j, c2.k, vid);
		return vid;
	}

	int directvid(Corner c1, Corner c2)
	/*
	 * Compute zero-crossing on the line c1-c2, and directly put into vertex
	 * list
	 */
	{
		KdPoint3 v = new KdPoint3();
		KdPoint3 a = new KdPoint3(), b = new KdPoint3();
		int vid;
		a.x = c1.i;
		a.y = c1.j;
		a.z = c1.k;
		b.x = c2.i;
		b.y = c2.j;
		b.z = c2.k;
		interpolate(a, b, c1.value, c2.value, v); /* position */
		vid = directadd(verts, v); /* save vertex */
		return vid;
	}

	/* Compute the asymptote on the face formed by c1,c2,c3, and c4 */
	/* c1 and c4 are diagonal (same sign), so are c2 and c3 */
	/* Decide the polarity of the asymptote: 1 for positive; 0 otherwise */
	int asym(Corner c1, Corner c2, Corner c3, Corner c4) {
		double alpha;
		int i1, j1, k1, i2, j2, k2, ioff1, joff1, koff1, ioff2, joff2, koff2;
		double v1, v2, v3, v4;

		if (c1.value > 0)
			alpha = (c1.value) * (c4.value) - (c2.value) * (c3.value);
		else
			alpha = (c2.value) * (c3.value) - (c1.value) * (c4.value);

		if (alpha == 0) {
			ioff1 = c1.i - c2.i;
			joff1 = c1.j - c2.j;
			koff1 = c1.k - c2.k;
			ioff2 = c1.i - c3.i;
			joff2 = c1.j - c3.j;
			koff2 = c1.k - c3.k;
			i1 = c1.i + ioff1;
			j1 = c1.j + joff1;
			k1 = c1.k + koff1;
			if (i1 <= 0)
				i1 = 0;
			if (i1 >= nx)
				i1 = nx - 1;
			if (j1 <= 0)
				j1 = 0;
			if (j1 >= ny)
				j1 = ny - 1;
			if (k1 <= 0)
				k1 = 0;
			if (k1 >= nz)
				k1 = nz - 1;
			i2 = c1.i + ioff2;
			j2 = c1.j + joff2;
			k2 = c1.k + koff2;
			if (i2 <= 0)
				i2 = 0;
			if (i2 >= nx)
				i2 = nx - 1;
			if (j2 <= 0)
				j2 = 0;
			if (j2 >= ny)
				j2 = ny - 1;
			if (k2 <= 0)
				k2 = 0;
			if (k2 >= nz)
				k2 = nz - 1;
			//v1 = data[i1][j1][k1] + data[i2][j2][k2];
			v1 = data[i1+nx*j1+nx*ny*k1] + data[i2+nx*j2+nx*ny*k2];
			i1 = c4.i - ioff1;
			j1 = c4.j - joff1;
			k1 = c4.k - koff1;
			if (i1 <= 0)
				i1 = 0;
			if (i1 >= nx)
				i1 = nx - 1;
			if (j1 <= 0)
				j1 = 0;
			if (j1 >= ny)
				j1 = ny - 1;
			if (k1 <= 0)
				k1 = 0;
			if (k1 >= nz)
				k1 = nz - 1;
			i2 = c4.i - ioff2;
			j2 = c4.j - joff2;
			k2 = c4.k - koff2;
			if (i2 <= 0)
				i2 = 0;
			if (i2 >= nx)
				i2 = nx - 1;
			if (j2 <= 0)
				j2 = 0;
			if (j2 >= ny)
				j2 = ny - 1;
			if (k2 <= 0)
				k2 = 0;
			if (k2 >= nz)
				k2 = nz - 1;
			//v4 = data[i1][j1][k1] + data[i2][j2][k2];
			v4 = data[i1+nx*j1+nx*ny*k1] + data[i2+nx*j2+nx*ny*k2];

			i1 = c2.i - ioff1;
			j1 = c2.j - joff1;
			k1 = c2.k - koff1;
			if (i1 <= 0)
				i1 = 0;
			if (i1 >= nx)
				i1 = nx - 1;
			if (j1 <= 0)
				j1 = 0;
			if (j1 >= ny)
				j1 = ny - 1;
			if (k1 <= 0)
				k1 = 0;
			if (k1 >= nz)
				k1 = nz - 1;
			i2 = c2.i + ioff2;
			j2 = c2.j + joff2;
			k2 = c2.k + koff2;
			if (i2 <= 0)
				i2 = 0;
			if (i2 >= nx)
				i2 = nx - 1;
			if (j2 <= 0)
				j2 = 0;
			if (j2 >= ny)
				j2 = ny - 1;
			if (k2 <= 0)
				k2 = 0;
			if (k2 >= nz)
				k2 = nz - 1;
			//v2 = data[i1][j1][k1] + data[i2][j2][k2];
			v2 = data[i1+nx*j1+nx*ny*k1] + data[i2+nx*j2+nx*ny*k2];

			i1 = c3.i + ioff1;
			j1 = c3.j + joff1;
			k1 = c3.k + koff1;
			if (i1 <= 0)
				i1 = 0;
			if (i1 >= nx)
				i1 = nx - 1;
			if (j1 <= 0)
				j1 = 0;
			if (j1 >= ny)
				j1 = ny - 1;
			if (k1 <= 0)
				k1 = 0;
			if (k1 >= nz)
				k1 = nz - 1;
			i2 = c3.i - ioff2;
			j2 = c3.j - joff2;
			k2 = c3.k - koff2;
			if (i2 <= 0)
				i2 = 0;
			if (i2 >= nx)
				i2 = nx - 1;
			if (j2 <= 0)
				j2 = 0;
			if (j2 >= ny)
				j2 = ny - 1;
			if (k2 <= 0)
				k2 = 0;
			if (k2 >= nz)
				k2 = nz - 1;
			//v3 = data[i1][j1][k1] + data[i2][j2][k2];
			v3 = data[i1+nx*j1+nx*ny*k1] + data[i2+nx*j2+nx*ny*k2];

			if (c1.value > 0) {
				alpha = (c1.value) * v4 + (c4.value) * v1 - v2 * (c3.value)
						- v3 * (c2.value);
			} else {
				alpha = v2 * (c3.value) + v3 * (c2.value) - (c1.value) * v4
						- (c4.value) * v1;
			}
		}

		if (alpha > 0)
			return 1;
		else
			return 0;
	}

	public ArrayList<EdgeIndex> getEdges() {
		ArrayList<EdgeIndex> edges = new ArrayList<EdgeIndex>();
		for (EdgeIndex e : edgelist) {
			while (e != null) {
				edges.add(e);
				e = e.next;
			}
		}
		return edges;
	}

	/**
	 * Original implementation of iso-surface algorithm
	 * 
	 * @param vol
	 *            Levelel set
	 * @param connectivity
	 *            connectivityectivity rule
	 * @param Level
	 *            iso-levelel
	 * @param inclusive
	 *            true if algorithm should include iso-levelel in surface
	 * @return iso-surface
	 */
	public final void execute() {
	    
	    System.out.print("\nConnectivity Consistent Iso-Surface Generation");
	    
	    verts = new ArrayList<KdPoint3>();
		for (int xyz=0;xyz<nxyz;xyz++) {
	        float tmp = lvlImage[xyz] - level;
	        /* The following is for TGDM purpose */
            if (tmp > -0.1 && tmp <= 0) tmp = -0.1f;
            else if (tmp < 0.1 && tmp > 0) tmp = 0.1f;
            //Case where tmp==0 is already handled?
            if (tmp == 0) {
                if (inclusive) {
                    // Points with value exactly = isovalue is considered inside
                    lvlImage[xyz] = 0.01f;
                } else {
                    // Points with value exactly = isovalue is considered outside
                    lvlImage[xyz] = -0.01f;
                }
            } else {
                lvlImage[xyz] = tmp;
            }
        }
        data = lvlImage;
        
		edgelist = new EdgeIndex[2 * HASHSIZE];
		trilist = new ArrayList<TriangleIndex>();

	    for (int x=0;x<nx-1;x++) for (int y=0;y<ny-1;y++) for (int z=0;z<nz-1;z++) {
	         docube(x,y,z);   
	    }
	    
	    // output
	    constructMesh();
	}
	
	/*
	 * A surface is created inside a voxel (if there is one). First, an index is
	 * created. The bits in this index indicate if the Data at a corner of the
	 * voxel is higher (1) or lower (0) than <0>. To speed things up, if <idx>
	 * has only one kind of bit (i.e. it is completely inside or outside of the
	 * levelel) no further action is performed. Otherwise, the index is rotated
	 * (using a permutation of the voxel corner indices) to a standard index,
	 * which is used to create a standard surface.
	 */
	private void docube(int x, int y, int z) {

		int idx;
		int[] std = new int[] { 0 }; /* Index and standard index */
		int perm; /* Used permutation */
		int[] rev = new int[] { 0 }; /*
									 * Indicates whether <idx> was inverted
									 * before rotation
									 */

		idx = levelindex(x, y, z);
		if (idx == 0 || idx == 255) /* Inside || Outside */
			return;

		perm = rightperm(idx, std, rev);
		standardtriangles(permute(perm, std[0]), perm, rev[0], x, y, z);
	}

	/*
	 * Returns the binary index for a voxel corresponding to its polarity. A '0'
	 * means the corner of the voxel is inside the surface; a '1' means it is
	 * outside.
	 */
	int levelindex(int x, int y, int z) {
		int dz, dy, dx;
		int idx; /* The index to be returned */

		idx = 0; /* b00000000; */
		for (dz = 0; dz < 2; dz++)
			for (dy = 0; dy < 2; dy++)
				for (dx = 0; dx < 2; dx++)
					//idx |= ((data[x + dx][y + dy][z + dz] > 0) ? 1 : 0) << ((dz << 2)
					//		| (dy << 1) | (dx));
					idx |= ((data[x + dx + nx*(y + dy) + nx*ny*(z + dz)] > 0) ? 1 : 0) << ((dz << 2)
							| (dy << 1) | (dx));
		return (idx);
	}
	    
	/*
	 * Determines the right permutation for a certain index, in order to reduce
	 * building a plane in a voxel to one of 15 standard cases. If <idx> has
	 * more than 4 bits (i.e. more 1s than 0s) it is inverted first. This is
	 * returned as a flag <rev>. First all 1s are counted, then, maybe after
	 * inversion, <idx> is compared with all 256 bit-patterns and a permutation
	 * index is returned.
	 */
	/* idx: The index from the voxel */
	/* std: The corresponding standard index */
	/* rev: Flag indicating if <idx> was inverted first */
	private int rightperm(int idx, int[] std, int[] rev) {
	
		int bit, ones;

		for (ones = 0, bit = 0; bit < 8; bit++)
			ones += ((idx >> bit) & 1);
		if ((rev[0] = (ones > 4) ? 1 : 0) != 0)
			idx ^= 255; /* b11111111; */

		std[0] = idx; /* If <idx> is inverted, then <*std> != <idx> */

		switch (idx) {
		case 0: /* b00000000: */
		case 8: /* b00001000: */
		case 10: /* b00001010: */
		case 9: /* b00001001: */
		case 24: /* b00011000: */
		case 11: /* b00001011: */
		case 26: /* b00011010: */
		case 41: /* b00101001: */
		case 15: /* b00001111: */
		case 43: /* b00101011: */
		case 60: /* b00111100: */
		case 29: /* b00011101: */
		case 30: /* b00011110: */
		case 105: /* b01101001: */
		case 27: /* b00011011: */
			return (0);
		case 2: /* b00000010: */

		case 3: /* b00000011: */
		case 6: /* b00000110: */
		case 66: /* b01000010: */
		case 7: /* b00000111: */
		case 67: /* b01000011: */
		case 22: /* b00010110: */
		case 23: /* b00010111: */
		case 90: /* b01011010: */
		case 78: /* b01001110: */
		case 75: /* b01001011: */
		case 150: /* b10010110: */
		case 71: /* b01000111: */
			return (1);
		case 1: /* b00000001: */
		case 5: /* b00000101: */
		case 129: /* b10000001: */
		case 13: /* b00001101: */
		case 133: /* b10000101: */
		case 73: /* b01001001: */
		case 77: /* b01001101: */
		case 195: /* b11000011: */
		case 139: /* b10001011: */
		case 135: /* b10000111: */
		case 141: /* b10001101: */
			return (2);
		case 4: /* b00000100: */
		case 12: /* b00001100: */
		case 36: /* b00100100: */
		case 14: /* b00001110: */
		case 44: /* b00101100: */
		case 134: /* b10000110: */
		case 142: /* b10001110: */
		case 165: /* b10100101: */
		case 39: /* b00100111: */
		case 45: /* b00101101: */
		case 46: /* b00101110: */
			return (3);
		case 20: /* b00010100: */
		case 21: /* b00010101: */
		case 37: /* b00100101: */
		case 85: /* b01010101: */
		case 102: /* b01100110: */
		case 116: /* b01110100: */
		case 101: /* b01100101: */
		case 53: /* b00110101: */
			return (4);
		case 64: /* b01000000: */
		case 80: /* b01010000: */
		case 96: /* b01100000: */
		case 112: /* b01110000: */
		case 82: /* b01010010: */
		case 97: /* b01100001: */
		case 240: /* b11110000: */
		case 113: /* b01110001: */
		case 226: /* b11100010: */
		case 210: /* b11010010: */
		case 114: /* b01110010: */
			return (5);
		case 128: /* b10000000: */
		case 160: /* b10100000: */
		case 130: /* b10000010: */
		case 162: /* b10100010: */
		case 161: /* b10100001: */
		case 146: /* b10010010: */
		case 170: /* b10101010: */
		case 178: /* b10110010: */
		case 153: /* b10011001: */
		case 169: /* b10101001: */
		case 163: /* b10100011: */
			return (6);
		case 34: /* b00100010: */
		case 18: /* b00010010: */
		case 50: /* b00110010: */
		case 98: /* b01100010: */
		case 51: /* b00110011: */
		case 83: /* b01010011: */
		case 99: /* b01100011: */
			return (7);
		case 32: /* b00100000: */
		case 224: /* b11100000: */
		case 164: /* b10100100: */
		case 104: /* b01101000: */
		case 232: /* b11101000: */
		case 180: /* b10110100: */
		case 228: /* b11100100: */
			return (8);
		case 136: /* b10001000: */
		case 132: /* b10000100: */
		case 140: /* b10001100: */
		case 137: /* b10001001: */
		case 204: /* b11001100: */
		case 197: /* b11000101: */
		case 201: /* b11001001: */
			return (9);
		case 16: /* b00010000: */
		case 17: /* b00010001: */
		case 19: /* b00010011: */
		case 25: /* b00011001: */
		case 58: /* b00111010: */
		case 57: /* b00111001: */
			return (10);
		case 68: /* b01000100: */
		case 196: /* b11000100: */
		case 100: /* b01100100: */
		case 148: /* b10010100: */
		case 212: /* b11010100: */
		case 172: /* b10101100: */
		case 108: /* b01101100: */
			return (11);
		case 192: /* b11000000: */
		case 144: /* b10010000: */
		case 208: /* b11010000: */
		case 193: /* b11000001: */
		case 177: /* b10110001: */
		case 225: /* b11100001: */
		case 209: /* b11010001: */
			return (12);
		case 48: /* b00110000: */
		case 176: /* b10110000: */
		case 56: /* b00111000: */
		case 216:/* b11011000: */
		case 120: /* b01111000: */
		case 184: /* b10111000: */
			return (13);
		case 84: /* b01010100: */
		case 88: /* b01011000: */
		case 89: /* b01011001: */
		case 92: /* b01011100: */
			return (14);
		case 138: /* b10001010: */
		case 74: /* b01001010: */
		case 106: /* b01101010: */
		case 202: /* b11001010: */
			return (15);
		case 33: /* b00100001: */
		case 49: /* b00110001: */
		case 52: /* b00110100: */
		case 54: /* b00110110: */
			return (16);
		case 65: /* b01000001: */
		case 69: /* b01000101: */
		case 70: /* b01000110: */
		case 86: /* b01010110: */
			return (17);
		case 40: /* b00101000: */
		case 42: /* b00101010: */
		case 38: /* b00100110: */
		case 166: /* b10100110: */
			return (18);
		case 35: /* b00100011: */
		case 131: /* b10000011: */
		case 147: /* b10010011: */
			return (19);
		case 72: /* b01001000: */
		case 200: /* b11001000: */
		case 194: /* b11000010: */
		case 198: /* b11000110: */
			return (20);
		case 81: /* b01010001: */
		case 145: /* b10010001: */
		case 149: /* b10010101: */
			return (21);
		case 168: /* b10101000: */
		case 152: /* b10011000: */
		case 154: /* b10011010: */
			return (22);
		case 76: /* b01001100: */
		case 28: /* b00011100: */
		case 156: /* b10011100: */
			return (23);
		}

		return (0); /* Something is wrong if this line is excecuted */
	}
	    
	/* Performes a permutation on the index <idx> and returns the result. */
	/* perm: Permutation number */
	/* idx: Old index */
	private int permute(int perm, int idx) {
		int bit;
		int newidx; /* Permuted index */

		for (newidx = 0, bit = 0; bit < 8; bit++)
			newidx |= ((idx >> Perms[perm][bit]) & 1) << bit;
		return (newidx);
	}
	
	/* Add v to the list of vertices */
	void addtovertices(ArrayList<KdPoint3> vertices, KdPoint3 v) {
		// System.out.println("ADD 1 "+v);
		vertices.add(v);
	}

	/* Directly add a new vertex to the list of vertices */
	private int directadd(ArrayList<KdPoint3> vertices, KdPoint3 v) {
		// System.out.println("ADD 2 "+v);
		vertices.add(v);
		return vertices.size() - 1;
	}

	/*
	 * From two points of differing sign, interpolate to surface crossing;
	 * linear interpolation is applied
	 */
	void interpolate(KdPoint3 p1, KdPoint3 p2, double v1, double v2, KdPoint3 newp) {
		int i = 0;
		KdPoint3 pos = new KdPoint3(), neg = new KdPoint3();
		double tmp;
		if (v1 < 0) {
			pos.x = p2.x;
			pos.y = p2.y;
			pos.z = p2.z;
			neg.x = p1.x;
			neg.y = p1.y;
			neg.z = p1.z;
		} else {
			pos.x = p1.x;
			pos.y = p1.y;
			pos.z = p1.z;
			neg.x = p2.x;
			neg.y = p2.y;
			neg.z = p2.z;
			tmp = v1;
			v1 = v2;
			v2 = tmp;
		}

		tmp = v2 - v1;
		if (tmp == 0)
			return;
		newp.x = (float) ((neg.x * v2 - pos.x * v1) / tmp);
		newp.y = (float) ((neg.y * v2 - pos.y * v1) / tmp);
		newp.z = (float) ((neg.z * v2 - pos.z * v1) / tmp);
	}

	private static final float EPS = 1E-3f;

	/* Add new triangle to the list of triangles */
	/* If rev is TRUE, the order of the three vertices are reversed */
	private void addtriangle(int vid1, int vid2, int vid3, int rev)	{
		int i;
		TriangleIndex tri = new TriangleIndex(); /* *triangles; */

		if (rev == 0) {
			tri.vid1 = vid1;
			tri.vid2 = vid2;
			tri.vid3 = vid3;
		} else {
			tri.vid1 = vid3;
			tri.vid2 = vid2;
			tri.vid3 = vid1;
		}
		// System.out.format("ADD TRIANGLE %d %d %d\n",vid1,vid2,vid3);
		trilist.add(tri);
	}

	    
	/*
	 * Here the standard Triangle Strips are listed. They are just a list of
	 * edges which are intersected by the surface at <0>, indicated by its two
	 * endpoint indices. An endpoint index is an integer ranging from 0 (bin
	 * 000) to 7 (bin 111; the bits in these indices are ordered ZYX). A surface
	 * through a voxel consists of up to 4 disjoint polygons, each spanned by up
	 * to 6 vertices.
	 */
	/* idx: Binary index */
	/* perm: Permutation index */
	/* rev: Indicating whether the index is inverted */
	/* z,y,x: Origin of the voxel */
	void standardtriangles(int idx, int perm, int rev, int x, int y, int z)	{
		int e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12; /*
																	 * Vertex
																	 * ID's
																	 */
		KdPoint3 V0;
		// int n;
		int subcase;
		int alphab;
		int alpha1, alpha2, alpha3, alpha4, alpha5, alpha6; /*
															 * Asymptotes at
															 * ambiguous faces
															 */

		Corner[] corners = new Corner[8];
		for (int i = 0; i < corners.length; i++)
			corners[i] = new Corner();
		corners[0].i = x;
		corners[0].j = y;
		corners[0].k = z;

		corners[1].i = x + 1;
		corners[1].j = y;
		corners[1].k = z;

		corners[2].i = x;
		corners[2].j = y + 1;
		corners[2].k = z;

		corners[3].i = x + 1;
		corners[3].j = y + 1;
		corners[3].k = z;

		corners[4].i = x;
		corners[4].j = y;
		corners[4].k = z + 1;

		corners[5].i = x + 1;
		corners[5].j = y;
		corners[5].k = z + 1;

		corners[6].i = x;
		corners[6].j = y + 1;
		corners[6].k = z + 1;

		corners[7].i = x + 1;
		corners[7].j = y + 1;
		corners[7].k = z + 1;

		for (Corner c : corners) {
			//c.value = data[c.i][c.j][c.k];
			c.value = data[c.i+nx*c.j+nx*ny*c.k];
		}
		// System.out.format("(%d,%d,%d)\n",x,y,z);
		// System.out.format("0)%6.2f 1)%6.2f 2)%6.2f 3)%6.2f 4)%6.2f 5)%6.2f 6)%6.2f 7)%6.2f\n",
		// corners[0].value,corners[1].value,corners[2].value,corners[3].value,
		// corners[4].value,corners[5].value,corners[6].value,corners[7].value);
		// n = connectivity;

		switch (idx) {
		case 0: /* b00000000: *//* type 0 */
			break;
		case 8: /* b00001000: *//* type 1 */
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e12, e3, e2, rev);
			break;
		case 10: /* b00001010: *//* type 2 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e1, e10, e12, rev);
			addtriangle(e12, e3, e1, rev);
			break;
		case 9: /* b00001001: *//* type 3 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6)
				alpha2 = 1;
			else
				alpha2 = 0;

			if (rev != 0)
				alpha2 = 1 - alpha2;

			switch (alpha2) {
			case 0: /* Connectivityect 1 & 2 */
				addtriangle(e1, e4, e9, rev);
				addtriangle(e3, e2, e12, rev);
				break;
			case 1: /* Connectivityect 0 & 3 */
				addtriangle(e1, e2, e12, rev);
				addtriangle(e12, e9, e1, rev);
				addtriangle(e12, e3, e4, rev);
				addtriangle(e4, e9, e12, rev);
				break;
			}
			break;
		case 24: /* b00011000: *//* type 4 */
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if ((rev == 0 && connectivity == CONNECTIVITY_26_6)
					|| (rev == 1 && connectivity == CONNECTIVITY_6_26))
				alphab = 1;
			else
				alphab = 0;

			if (alphab == 1) {
				addtriangle(e12, e8, e5, rev);
				addtriangle(e12, e5, e2, rev);
				addtriangle(e2, e5, e9, rev);
				addtriangle(e2, e9, e3, rev);
				addtriangle(e3, e9, e8, rev);
				addtriangle(e3, e8, e12, rev);
			} else {
				addtriangle(e5, e9, e8, rev);
				addtriangle(e12, e3, e2, rev);
			}
			break;
		case 11: /* b00001011: *//* type 5 */
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e9, e10, e12, rev);
			addtriangle(e12, e3, e9, rev);
			addtriangle(e3, e4, e9, rev);
			break;
		case 26: /* b00011010: *//* type 6 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6)
				alpha6 = 1;
			else
				alpha6 = 0;

			if (rev != 0)
				alpha6 = 1 - alpha6;
			switch (alpha6) {
			case 0:
				addtriangle(e8, e5, e9, rev);
				addtriangle(e1, e10, e12, rev);
				addtriangle(e12, e3, e1, rev);
				break;
			case 1:
				addtriangle(e3, e1, e9, rev);
				addtriangle(e3, e9, e8, rev);
				addtriangle(e3, e8, e12, rev);
				addtriangle(e12, e8, e5, rev);
				addtriangle(e12, e5, e10, rev);
				break;
			}
			break;

		case 41: /* b00101001: *//* type 7 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e6 = vertid(corners[Perms[perm][5]], corners[Perms[perm][7]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6) {
				alpha2 = 1;
				alpha4 = 1;
				alpha6 = 1;
			} else {
				alpha2 = 0;
				alpha4 = 0;
				alpha6 = 0;
			}
			if (rev != 0) {
				alpha2 = 1 - alpha2;
				alpha4 = 1 - alpha4;
				alpha6 = 1 - alpha6;
			}
			switch (((alpha2 << 2) + (alpha4 << 1) + alpha6)) {
			case 0:
				addtriangle(e6, e10, e5, rev);
				addtriangle(e3, e2, e12, rev);
				addtriangle(e4, e9, e1, rev);
				break;
			case 1:
				addtriangle(e3, e2, e12, rev);
				addtriangle(e1, e4, e10, rev);
				addtriangle(e4, e6, e10, rev);
				addtriangle(e4, e5, e6, rev);
				addtriangle(e4, e9, e5, rev);
				break;
			case 2:
				addtriangle(e1, e4, e9, rev);
				addtriangle(e5, e6, e12, rev);
				addtriangle(e5, e12, e3, rev);
				addtriangle(e5, e3, e2, rev);
				addtriangle(e5, e2, e10, rev);
				break;
			case 3:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				addtriangle(e0, e6, e12, rev);
				addtriangle(e0, e12, e3, rev);
				addtriangle(e0, e3, e2, rev);
				addtriangle(e0, e2, e10, rev);
				addtriangle(e0, e10, e1, rev);
				addtriangle(e0, e1, e4, rev);
				addtriangle(e0, e4, e9, rev);
				addtriangle(e0, e9, e5, rev);
				addtriangle(e0, e5, e6, rev);
				break;
			case 4:
				addtriangle(e5, e6, e10, rev);
				addtriangle(e9, e1, e2, rev);
				addtriangle(e9, e2, e12, rev);
				addtriangle(e9, e12, e3, rev);
				addtriangle(e9, e3, e4, rev);
				break;
			case 5:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				addtriangle(e0, e5, e6, rev);
				addtriangle(e0, e6, e10, rev);
				addtriangle(e0, e10, e1, rev);
				addtriangle(e0, e1, e2, rev);
				addtriangle(e0, e2, e12, rev);
				addtriangle(e0, e12, e3, rev);
				addtriangle(e0, e3, e4, rev);
				addtriangle(e0, e4, e9, rev);
				addtriangle(e0, e9, e5, rev);
				break;
			case 6:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				addtriangle(e0, e1, e2, rev);
				addtriangle(e0, e2, e10, rev);
				addtriangle(e0, e10, e5, rev);
				addtriangle(e0, e5, e6, rev);
				addtriangle(e0, e6, e12, rev);
				addtriangle(e0, e12, e3, rev);
				addtriangle(e0, e3, e4, rev);
				addtriangle(e0, e4, e9, rev);
				addtriangle(e0, e9, e1, rev);

				break;
			case 7:
				addtriangle(e1, e2, e10, rev);
				addtriangle(e9, e5, e6, rev);
				addtriangle(e9, e6, e12, rev);
				addtriangle(e9, e12, e4, rev);
				addtriangle(e3, e4, e12, rev);
				break;
			}
			break;

		/* The following cases always have rev == 0 */
		case 15: /* b00001111: *//* type 8 */
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e11 = vertid(corners[Perms[perm][2]], corners[Perms[perm][6]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e9, e10, e12, rev);
			addtriangle(e12, e11, e9, rev);
			break;
		case 43: /* b00101011: *//* type 9 */
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e6 = vertid(corners[Perms[perm][5]], corners[Perms[perm][7]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e5, e6, e9, rev);
			addtriangle(e6, e12, e4, rev);
			addtriangle(e12, e3, e4, rev);
			addtriangle(e4, e9, e6, rev);
			break;
		case 60: /* b00111100: *//* type 10 */
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e6 = vertid(corners[Perms[perm][5]], corners[Perms[perm][7]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e11 = vertid(corners[Perms[perm][2]], corners[Perms[perm][6]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6) {
				addtriangle(e8, e6, e12, rev);
				addtriangle(e12, e11, e8, rev);
				addtriangle(e4, e2, e9, rev);
				addtriangle(e2, e10, e9, rev);
				// alpha3 = 1;
				// alpha4 = 1;
			} else {
				addtriangle(e9, e8, e6, rev);
				addtriangle(e6, e10, e9, rev);
				addtriangle(e12, e11, e4, rev);
				addtriangle(e4, e2, e12, rev);

				// alpha3 = 0;
				// alpha4 = 0;
			}
			break;
		case 29: /* b00011101: *//* type 11 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e11 = vertid(corners[Perms[perm][2]], corners[Perms[perm][6]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e5, e1, e2, rev);
			addtriangle(e2, e12, e11, rev);
			addtriangle(e11, e8, e5, rev);
			addtriangle(e11, e5, e2, rev);
			break;
		case 30: /* b00011110: *//* type 12 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e11 = vertid(corners[Perms[perm][2]], corners[Perms[perm][6]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6) {
				alpha3 = 1;
				alpha6 = 1;
			} else {
				alpha3 = 0;
				alpha6 = 0;
			}
			switch (((alpha3 << 1) + alpha6)) {
			case 0:
				addtriangle(e8, e5, e9, rev);
				addtriangle(e1, e10, e12, rev);
				addtriangle(e12, e11, e4, rev);
				addtriangle(e4, e1, e12, rev);
				break;
			case 1:
				addtriangle(e1, e9, e8, rev);
				addtriangle(e1, e8, e12, rev);
				addtriangle(e8, e5, e12, rev);
				addtriangle(e5, e10, e12, rev);
				addtriangle(e4, e1, e12, rev);
				addtriangle(e12, e11, e4, rev);
				break;
			case 2:
				addtriangle(e9, e4, e5, rev);
				addtriangle(e4, e12, e5, rev);
				addtriangle(e5, e12, e8, rev);
				addtriangle(e12, e11, e8, rev);
				addtriangle(e4, e1, e12, rev);
				addtriangle(e1, e10, e12, rev);
				break;
			case 3:
				addtriangle(e4, e1, e9, rev);
				addtriangle(e8, e5, e10, rev);
				addtriangle(e10, e12, e11, rev);
				addtriangle(e11, e8, e10, rev);
				break;
			}
			break;
		case 105: /* b01101001: *//* type 13 */
			e1 = vertid(corners[Perms[perm][0]], corners[Perms[perm][1]]);
			e2 = vertid(corners[Perms[perm][1]], corners[Perms[perm][3]]);
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e6 = vertid(corners[Perms[perm][5]], corners[Perms[perm][7]]);
			e7 = vertid(corners[Perms[perm][6]], corners[Perms[perm][7]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e9 = vertid(corners[Perms[perm][0]], corners[Perms[perm][4]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e11 = vertid(corners[Perms[perm][2]], corners[Perms[perm][6]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);

			if (connectivity == CONNECTIVITY_18_6
					|| connectivity == CONNECTIVITY_26_6)
				subcase = 63;
			else
				subcase = 0;

			switch (subcase) {
			case 0:
				addtriangle(e2, e12, e3, rev);
				addtriangle(e5, e6, e10, rev);
				addtriangle(e7, e8, e11, rev);
				addtriangle(e4, e9, e1, rev);
				break;
			case 1:
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						0, 'I');
				break;
			case 2:
				HandleCase13(e11, e3, e12, e7, e9, e1, e10, e5, e8, e4, e6, e2,
						0, 'I');
				break;
			case 3:
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						0, 'G');
				break;
			case 4:
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						0, 'I');
				break;
			case 5:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						e0, 'H');
				break;
			case 6:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e6, e12, e2, e10, e8, e11, e4, e9, e5, e7, e1, e3,
						e0, 'H');
				break;
			case 7:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						e0, 'L');
				break;
			case 8:
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						0, 'I');
				break;
			case 9:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e8, e9, e4, e11, e6, e10, e2, e12, e7, e5, e3, e1,
						e0, 'H');
				break;
			case 10:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						e0, 'H');
				break;
			case 11:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						e0, 'L');
				break;
			case 12:
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						0, 'G');
				break;
			case 13:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e11, e3, e12, e7, e9, e1, e10, e5, e8, e4, e6, e2,
						e0, 'E');
				break;
			case 14:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						e0, 'E');
				break;
			case 15:
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						0, 'C');
				break;
			case 16:
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						0, 'I');
				break;
			case 17:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e10, e1, e9, e5, e12, e3, e11, e7, e6, e2, e8, e4,
						e0, 'H');
				break;
			case 18:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e11, e3, e12, e7, e9, e1, e10, e5, e8, e4, e6, e2,
						e0, 'H');
				break;
			case 19:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						e0, 'E');
				break;
			case 20:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						e0, 'H');
				break;
			case 21:
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						0, 'K');
				break;
			case 22:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						e0, 'F');
				break;
			case 23:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e5, e8, e7, e6, e1, e4, e3, e2, e10, e9, e12, e11,
						e0, 'D');
				break;
			case 24:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e3, e4, e1, e2, e7, e8, e5, e6, e12, e11, e10, e9,
						e0, 'H');
				break;
			case 25:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						e0, 'F');
				break;
			case 26:
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						0, 'K');
				break;
			case 27:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						e0, 'D');
				break;
			case 28:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						e0, 'L');
				break;
			case 29:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e12, e7, e11, e3, e10, e5, e9, e1, e2, e6, e4, e8,
						e0, 'D');
				break;
			case 30:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						e0, 'D');
				break;
			case 31:
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						0, 'B');
				break;
			case 32:
				HandleCase13(e5, e8, e7, e6, e1, e4, e3, e2, e10, e9, e12, e11,
						0, 'I');
				break;
			case 33:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						e0, 'H');
				break;
			case 34:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e12, e7, e11, e3, e10, e5, e9, e1, e2, e6, e4, e8,
						e0, 'H');
				break;
			case 35:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						e0, 'E');
				break;
			case 36:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e7, e6, e5, e8, e3, e2, e1, e4, e11, e12, e9, e10,
						e0, 'H');
				break;
			case 37:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						e0, 'F');
				break;
			case 38:
				HandleCase13(e8, e9, e4, e11, e6, e10, e2, e12, e7, e5, e3, e1,
						0, 'K');
				break;
			case 39:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e3, e4, e1, e2, e7, e8, e5, e6, e12, e11, e10, e9,
						e0, 'D');
				break;
			case 40:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e5, e8, e7, e6, e1, e4, e3, e2, e10, e9, e12, e11,
						e0, 'H');
				break;
			case 41:
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						0, 'K');
				break;
			case 42:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						e0, 'F');
				break;
			case 43:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						e0, 'D');
				break;
			case 44:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						e0, 'L');
				break;
			case 45:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e11, e3, e12, e7, e9, e1, e10, e5, e8, e4, e6, e2,
						e0, 'D');
				break;
			case 46:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e10, e1, e9, e5, e12, e3, e11, e7, e6, e2, e8, e4,
						e0, 'D');
				break;
			case 47:
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						0, 'B');
				break;
			case 48:
				HandleCase13(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12,
						0, 'G');
				break;
			case 49:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e3, e12, e7, e11, e1, e10, e5, e9, e4, e2, e8, e6,
						e0, 'E');
				break;
			case 50:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						e0, 'L');
				break;
			case 51:
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						0, 'C');
				break;
			case 52:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						e0, 'E');
				break;
			case 53:
				e0 = directvid(corners[Perms[perm][2]], corners[Perms[perm][5]]);
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						e0, 'D');
				break;
			case 54:
				e0 = directvid(corners[Perms[perm][3]], corners[Perms[perm][4]]);
				HandleCase13(e8, e9, e4, e11, e6, e10, e2, e12, e7, e5, e3, e1,
						e0, 'D');
				break;
			case 55:
				HandleCase13(e4, e11, e8, e9, e2, e12, e6, e10, e1, e3, e5, e7,
						0, 'B');
				break;
			case 56:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						e0, 'E');
				break;
			case 57:
				e0 = directvid(corners[Perms[perm][0]], corners[Perms[perm][7]]);
				HandleCase13(e6, e12, e2, e10, e8, e11, e4, e9, e5, e7, e1, e3,
						e0, 'D');
				break;
			case 58:
				e0 = directvid(corners[Perms[perm][1]], corners[Perms[perm][6]]);
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						e0, 'D');
				break;
			case 59:
				HandleCase13(e2, e10, e6, e12, e4, e9, e8, e11, e3, e1, e7, e5,
						0, 'B');
				break;
			case 60:
				HandleCase13(e3, e12, e7, e11, e1, e10, e5, e9, e4, e2, e8, e6,
						0, 'C');
				break;
			case 61:
				HandleCase13(e12, e7, e11, e3, e10, e5, e9, e1, e2, e6, e4, e8,
						0, 'B');
				break;
			case 62:
				HandleCase13(e9, e5, e10, e1, e11, e7, e12, e3, e4, e8, e2, e6,
						0, 'B');
				break;
			case 63:
				addtriangle(e3, e4, e11, rev);
				addtriangle(e1, e2, e10, rev);
				addtriangle(e5, e8, e9, rev);
				addtriangle(e12, e7, e6, rev);
				break;
			}
			break;

		case 27: /* b00011011: *//* type 14 */
			e3 = vertid(corners[Perms[perm][2]], corners[Perms[perm][3]]);
			e4 = vertid(corners[Perms[perm][0]], corners[Perms[perm][2]]);
			e5 = vertid(corners[Perms[perm][4]], corners[Perms[perm][5]]);
			e8 = vertid(corners[Perms[perm][4]], corners[Perms[perm][6]]);
			e10 = vertid(corners[Perms[perm][1]], corners[Perms[perm][5]]);
			e12 = vertid(corners[Perms[perm][3]], corners[Perms[perm][7]]);
			addtriangle(e8, e5, e10, rev);
			addtriangle(e10, e12, e3, rev);
			addtriangle(e3, e4, e8, rev);
			addtriangle(e8, e10, e3, rev);
			break;
		default:
			System.out.format("Unknown index!\n");
		}
	}
	/* HandleCase13: standard triangulation of the subcases of case 13 */
	void HandleCase13(int e1, int e2, int e3, int e4, int e5, int e6, int e7,
			int e8, int e9, int e10, int e11, int e12, int e0, char subcase) {
		switch (subcase) {
		case 'B':
			addtriangle(e5, e8, e9, 0);
			addtriangle(e6, e12, e7, 0);
			addtriangle(e11, e3, e2, 0);
			addtriangle(e11, e2, e10, 0);
			addtriangle(e10, e1, e4, 0);
			addtriangle(e10, e4, e11, 0);
			break;
		case 'C':
			addtriangle(e10, e1, e4, 0);
			addtriangle(e10, e4, e11, 0);
			addtriangle(e11, e2, e10, 0);
			addtriangle(e11, e3, e2, 0);
			addtriangle(e5, e6, e12, 0);
			addtriangle(e9, e5, e12, 0);
			addtriangle(e12, e7, e9, 0);
			addtriangle(e7, e8, e9, 0);
			break;
		case 'D':
			addtriangle(e5, e8, e9, 0);
			addtriangle(e0, e6, e10, 0);
			addtriangle(e0, e10, e1, 0);
			addtriangle(e0, e7, e6, 0);
			addtriangle(e0, e12, e7, 0);
			addtriangle(e0, e2, e12, 0);
			addtriangle(e0, e3, e2, 0);
			addtriangle(e0, e11, e3, 0);
			addtriangle(e0, e4, e11, 0);
			addtriangle(e0, e1, e4, 0);
			break;
		case 'E':
			addtriangle(e0, e12, e7, 0);
			addtriangle(e0, e7, e6, 0);
			addtriangle(e0, e6, e10, 0);
			addtriangle(e0, e10, e1, 0);
			addtriangle(e0, e1, e4, 0);
			addtriangle(e0, e4, e9, 0);
			addtriangle(e0, e9, e5, 0);
			addtriangle(e0, e5, e8, 0);
			addtriangle(e0, e8, e11, 0);
			addtriangle(e0, e11, e3, 0);
			addtriangle(e0, e3, e2, 0);
			addtriangle(e0, e2, e12, 0);
			break;
		case 'F':
			addtriangle(e0, e9, e1, 0);
			addtriangle(e0, e1, e4, 0);
			addtriangle(e0, e4, e11, 0);
			addtriangle(e0, e11, e3, 0);
			addtriangle(e0, e3, e2, 0);
			addtriangle(e0, e2, e12, 0);
			addtriangle(e0, e12, e7, 0);
			addtriangle(e0, e7, e6, 0);
			addtriangle(e0, e6, e10, 0);
			addtriangle(e0, e10, e5, 0);
			addtriangle(e0, e5, e8, 0);
			addtriangle(e0, e8, e9, 0);
			break;
		case 'G':
			addtriangle(e9, e1, e2, 0);
			addtriangle(e9, e2, e12, 0);
			addtriangle(e9, e12, e4, 0);
			addtriangle(e3, e4, e12, 0);
			addtriangle(e10, e5, e8, 0);
			addtriangle(e10, e8, e11, 0);
			addtriangle(e10, e11, e6, 0);
			addtriangle(e11, e7, e6, 0);
			break;
		case 'H':
			addtriangle(e11, e7, e8, 0);
			addtriangle(e0, e9, e1, 0);
			addtriangle(e0, e1, e2, 0);
			addtriangle(e0, e2, e10, 0);
			addtriangle(e0, e10, e5, 0);
			addtriangle(e0, e5, e6, 0);
			addtriangle(e0, e6, e12, 0);
			addtriangle(e0, e12, e3, 0);
			addtriangle(e0, e3, e4, 0);
			addtriangle(e0, e4, e9, 0);
			break;
		case 'I':
			addtriangle(e11, e7, e8, 0);
			addtriangle(e10, e5, e6, 0);
			addtriangle(e3, e4, e9, 0);
			addtriangle(e9, e12, e3, 0);
			addtriangle(e9, e1, e12, 0);
			addtriangle(e1, e2, e12, 0);
			break;

		case 'K': /* Complement of F */
			/* Top, front, right seperated */
			addtriangle(e3, e2, e12, 0);
			addtriangle(e9, e5, e8, 0);
			addtriangle(e6, e10, e1, 0);
			addtriangle(e1, e7, e6, 0);
			addtriangle(e1, e4, e7, 0);
			addtriangle(e4, e11, e7, 0);
			break;

		case 'L': /* top, front, bottom seperated */
			addtriangle(e0, e2, e10, 0);
			addtriangle(e0, e10, e5, 0);
			addtriangle(e0, e5, e8, 0);
			addtriangle(e0, e8, e9, 0);
			addtriangle(e0, e9, e1, 0);
			addtriangle(e0, e1, e4, 0);
			addtriangle(e0, e4, e11, 0);
			addtriangle(e0, e11, e7, 0);
			addtriangle(e0, e7, e6, 0);
			addtriangle(e0, e6, e12, 0);
			addtriangle(e0, e12, e3, 0);
			addtriangle(e0, e3, e2, 0);
			break;
		}
	}


	private void constructMesh() {
		float[] points = new float[verts.size() * 3];
		int indices[] = new int[trilist.size() * 3];
		data = null;
		int i = 0;
		for (TriangleIndex tri : trilist) {
			// System.out.format("%9d %9d %9d\n", tri.vid1, tri.vid2, tri.vid3);
			if (direction == COUNTER_CLOCKWISE) {
				indices[i++] = tri.vid3;
				indices[i++] = tri.vid2;
				indices[i++] = tri.vid1;
			} else {
				indices[i++] = tri.vid1;
				indices[i++] = tri.vid2;
				indices[i++] = tri.vid3;
			}
		}
		i = 0;
		for (KdPoint3 pt : verts) {
		    points[3*i+0] = pt.x;
		    points[3*i+1] = pt.y;
		    points[3*i+2] = pt.z;
		    i++;
		}
		pointList = points;
		triangleList = indices;
		
		return;
	}


}
