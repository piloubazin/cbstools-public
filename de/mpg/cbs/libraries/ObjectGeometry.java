package de.mpg.cbs.libraries;


/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectGeometry {
	
	// no data: used as a library of functions
	private final static int MaxObject = 1000000;
	
	public final static int SUPERIOR =  10;
	public final static int SUPEQUAL =  11;
	public final static int INFERIOR =  12;
	public final static int INFEQUAL =  13;
	public final static int EQUAL =  	14;
	public final static int UNEQUAL =  	15;
	public final static int NONE =      25;
	public final static int AND =       26;
	public final static int OR =        27;
	public final static int XOR =       28;
	
	public final static float SQR2 = (float)Math.sqrt(2.0f);
	public final static float SQR3 = (float)Math.sqrt(3.0f);
   // object manipulation
    
	/**
     *  Extract a surface that bounds the input volume
	 *	and returns the Euler characteristic
	 *  (not working properly)
     */
    public static final int cubicSurfaceFromObject(boolean img[][][], int nx, int ny, int nz) {
		float [][] vertex;
		int [][] edge;
		int [][] face;
		int Nv = 0;
		int Ne = 0;
		int Nf = 0;
		boolean cleanUp = true;
		
		int Mv=nx*ny*nz,Me=nx*ny*nz,Mf=nx*ny*nz;
		vertex = new float[Mv][3];
		edge = new int[Me][2];
		face = new int[Mf][4];
		boolean newFace;
		int Nv0,Ne0,Nf0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				Nv0 = Nv; Ne0 = Ne; Nf0 = Nf;
				newFace = false;
				// object: check for boundary
				// check the six directions
				// +X
				if (!img[x+1][y][z]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x+0.5f; vertex[Nv+0][1] = y-0.5f; vertex[Nv+0][2] = z-0.5f;
					vertex[Nv+1][0] = x+0.5f; vertex[Nv+1][1] = y+0.5f; vertex[Nv+1][2] = z-0.5f;
					vertex[Nv+2][0] = x+0.5f; vertex[Nv+2][1] = y+0.5f; vertex[Nv+2][2] = z+0.5f;
					vertex[Nv+3][0] = x+0.5f; vertex[Nv+3][1] = y-0.5f; vertex[Nv+3][2] = z+0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				// -X
				if (!img[x-1][y][z]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x-0.5f; vertex[Nv+0][1] = y-0.5f; vertex[Nv+0][2] = z-0.5f;
					vertex[Nv+1][0] = x-0.5f; vertex[Nv+1][1] = y+0.5f; vertex[Nv+1][2] = z-0.5f;
					vertex[Nv+2][0] = x-0.5f; vertex[Nv+2][1] = y+0.5f; vertex[Nv+2][2] = z+0.5f;
					vertex[Nv+3][0] = x-0.5f; vertex[Nv+3][1] = y-0.5f; vertex[Nv+3][2] = z+0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				// +Y
				if (!img[x][y+1][z]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x-0.5f; vertex[Nv+0][1] = y+0.5f; vertex[Nv+0][2] = z-0.5f;
					vertex[Nv+1][0] = x+0.5f; vertex[Nv+1][1] = y+0.5f; vertex[Nv+1][2] = z-0.5f;
					vertex[Nv+2][0] = x+0.5f; vertex[Nv+2][1] = y+0.5f; vertex[Nv+2][2] = z+0.5f;
					vertex[Nv+3][0] = x-0.5f; vertex[Nv+3][1] = y+0.5f; vertex[Nv+3][2] = z+0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				// -Y
				if (!img[x][y-1][z]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x-0.5f; vertex[Nv+0][1] = y-0.5f; vertex[Nv+0][2] = z-0.5f;
					vertex[Nv+1][0] = x+0.5f; vertex[Nv+1][1] = y-0.5f; vertex[Nv+1][2] = z-0.5f;
					vertex[Nv+2][0] = x+0.5f; vertex[Nv+2][1] = y-0.5f; vertex[Nv+2][2] = z+0.5f;
					vertex[Nv+3][0] = x-0.5f; vertex[Nv+3][1] = y-0.5f; vertex[Nv+3][2] = z+0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				// +Z
				if (!img[x][y][z+1]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x-0.5f; vertex[Nv+0][1] = y-0.5f; vertex[Nv+0][2] = z+0.5f;
					vertex[Nv+1][0] = x+0.5f; vertex[Nv+1][1] = y-0.5f; vertex[Nv+1][2] = z+0.5f;
					vertex[Nv+2][0] = x+0.5f; vertex[Nv+2][1] = y+0.5f; vertex[Nv+2][2] = z+0.5f;
					vertex[Nv+3][0] = x-0.5f; vertex[Nv+3][1] = y+0.5f; vertex[Nv+3][2] = z+0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				// -Z
				if (!img[x][y][z-1]) {
					// boundary: create face, vertices, edges
					newFace = true;
					// vertices
					vertex[Nv+0][0] = x-0.5f; vertex[Nv+0][1] = y-0.5f; vertex[Nv+0][2] = z-0.5f;
					vertex[Nv+1][0] = x+0.5f; vertex[Nv+1][1] = y-0.5f; vertex[Nv+1][2] = z-0.5f;
					vertex[Nv+2][0] = x+0.5f; vertex[Nv+2][1] = y+0.5f; vertex[Nv+2][2] = z-0.5f;
					vertex[Nv+3][0] = x-0.5f; vertex[Nv+3][1] = y+0.5f; vertex[Nv+3][2] = z-0.5f;
					// edges
					edge[Ne+0][0] = Nv+0; edge[Ne+0][1] = Nv+1;
					edge[Ne+1][0] = Nv+1; edge[Ne+1][1] = Nv+2;
					edge[Ne+2][0] = Nv+2; edge[Ne+2][1] = Nv+3;
					edge[Ne+3][0] = Nv+3; edge[Ne+3][1] = Nv+0;
					// face
					face[Nf+0][0] = Ne+0; face[Nf+0][1] = Ne+1; face[Nf+0][2] = Ne+2; face[Nf+0][3] = Ne+3;
					// update sizes
					Nv+=4; Ne+=4; Nf+=1;
				}
				if ( (newFace) && (!cleanUp) ) {
					// cleaning up: remove duplicates
					int[] idv = new int[Nv];
					for (int n=0;n<Nv;n++) idv[n] = -1;
					int lb = 0;
					for (int n=0;n<Nv;n++) {
						if (idv[n]==-1) {
							// new one
							idv[n] = lb;
							for (int m=n+1;m<Nv;m++) {
								if ( (vertex[m][0]==vertex[n][0]) && (vertex[m][1]==vertex[n][1]) && (vertex[m][2]==vertex[n][2]) ) {
									idv[m] = lb;
								}
							}
							lb++;
						}
					}
					int nV = lb;
					float[][] Vertex = new float[nV][3];
					// order the vertices
					for (int n=0;n<Nv;n++) {
						Vertex[idv[n]][0] = vertex[n][0]; Vertex[idv[n]][1] = vertex[n][1]; Vertex[idv[n]][2] = vertex[n][2];
					}
					// update the edges
					for (int n=0;n<Ne;n++) {
						edge[n][0] = idv[edge[n][0]]; edge[n][1] = idv[edge[n][1]];
					}
					// clean up edges
					int[] ide = new int[Ne];
					for (int n=0;n<Ne;n++) ide[n] = -1;
					lb = 0;
					for (int n=0;n<Ne;n++) {
						if (ide[n]==-1) {
							// new one
							ide[n] = lb;
							for (int m=n+1;m<Ne;m++) {
								if ( ( (edge[m][0]==edge[n][0]) && (edge[m][1]==edge[n][1]) )
									|| ( (edge[m][0]==edge[n][1]) && (edge[m][1]==edge[n][0]) ) ) {
									ide[m] = lb;
								}
							}
							lb++;
						}
					}
					int nE = lb;
					int[][] Edge = new int[nE][2];
					// order the vertices
					for (int n=0;n<Ne;n++) {
						Edge[ide[n]][0] = edge[n][0]; Edge[ide[n]][1] = edge[n][1];
					}
					// update the faces
					for (int n=0;n<Nf;n++) {
						face[n][0] = ide[face[n][0]]; face[n][1] = ide[face[n][1]]; face[n][2] = ide[face[n][2]]; face[n][3] = ide[face[n][3]];
					}
					idv = null;ide = null;
					// pass it to the new values
					Nv = nV; Ne = nE;
					for (int n=0;n<Nv;n++) {
						vertex[n][0] = Vertex[n][0]; vertex[n][1] = Vertex[n][1]; vertex[n][2] = Vertex[n][2];
					}
					Vertex = null;
					for (int n=0;n<Ne;n++) {
						edge[n][0] = Edge[n][0]; edge[n][1] = Edge[n][1];
					}
					Edge = null;
					System.gc();
				}
			}
			if (Nv+24>=Mv) {
				Mv += nx*ny*nz;
				float[][] tmv = new float[Mv][3];
				for (int n=0;n<Nv;n++) { tmv[n][0] = vertex[n][0]; tmv[n][1] = vertex[n][1]; tmv[n][2] = vertex[n][2]; }
				vertex = new float[Mv][3];
				for (int n=0;n<Nv;n++) { vertex[n][0] = tmv[n][0]; vertex[n][1] = tmv[n][1]; vertex[n][2] = tmv[n][2]; }
				tmv = null;
				System.gc();
			}
			if (Ne+24>=Me) {
				Me += nx*ny*nz;
				int[][] tme = new int[Me][2];
				for (int n=0;n<Ne;n++) { tme[n][0] = edge[n][0]; tme[n][1] = edge[n][1]; }
				edge = new int[Me][2];
				for (int n=0;n<Ne;n++) { edge[n][0] = tme[n][0]; edge[n][1] = tme[n][1]; }
				tme = null;
				System.gc();
			}
			if (Nf+12>=Mf) {
				Mf += nx*ny*nz;
				int[][] tmf = new int[Mf][4];
				for (int n=0;n<Nf;n++) { tmf[n][0] = face[n][0]; tmf[n][1] = face[n][1]; tmf[n][2] = face[n][2];  tmf[n][3] = face[n][3]; }
				face = new int[Mf][4];
				for (int n=0;n<Nf;n++) { face[n][0] = tmf[n][0]; face[n][1] = tmf[n][1]; face[n][2] = tmf[n][2];  face[n][3] = tmf[n][3]; }
				tmf = null;
				System.gc();
			}
		}
		if (cleanUp) {
			// cleaning up: remove duplicates
			int[] idv = new int[Nv];
			for (int n=0;n<Nv;n++) idv[n] = -1;
			int lb = 0;
			for (int n=0;n<Nv;n++) {
				if (idv[n]==-1) {
					// new one
					idv[n] = lb;
					for (int m=n+1;m<Nv;m++) {
						if ( (vertex[m][0]==vertex[n][0]) && (vertex[m][1]==vertex[n][1]) && (vertex[m][2]==vertex[n][2]) ) {
							idv[m] = lb;
						}
					}
					lb++;
				}
			}
			int nV = lb;
			float[][] Vertex = new float[nV][3];
			// order the vertices
			for (int n=0;n<Nv;n++) {
				Vertex[idv[n]][0] = vertex[n][0]; Vertex[idv[n]][1] = vertex[n][1]; Vertex[idv[n]][2] = vertex[n][2];
			}
			// update the edges
			for (int n=0;n<Ne;n++) {
				edge[n][0] = idv[edge[n][0]]; edge[n][1] = idv[edge[n][1]];
			}
			// clean up edges
			int[] ide = new int[Ne];
			for (int n=0;n<Ne;n++) ide[n] = -1;
			lb = 0;
			for (int n=0;n<Ne;n++) {
				if (ide[n]==-1) {
					// new one
					ide[n] = lb;
					for (int m=n+1;m<Ne;m++) {
						if ( ( (edge[m][0]==edge[n][0]) && (edge[m][1]==edge[n][1]) )
							|| ( (edge[m][0]==edge[n][1]) && (edge[m][1]==edge[n][0]) ) ) {
							ide[m] = lb;
						}
					}
					lb++;
				}
			}
			int nE = lb;
			int[][] Edge = new int[nE][2];
			// order the vertices
			for (int n=0;n<Ne;n++) {
				Edge[ide[n]][0] = edge[n][0]; Edge[ide[n]][1] = edge[n][1];
			}
			// update the faces
			for (int n=0;n<Nf;n++) {
				face[n][0] = ide[face[n][0]]; face[n][1] = ide[face[n][1]]; face[n][2] = ide[face[n][2]]; face[n][3] = ide[face[n][3]];
			}
			Nv = nV; Ne = nE;
		}		
		System.out.println("Surface: "+Nv+" vertices, "+Ne+" edges, "+Nf+" faces");
		// save the surface in MIPAV's format ?
		return Nv-Ne+Nf;
	}

	/**
     *  compute the thickness of the object at each point
	 *	(not properly functional yet)
     */
    public static final float[][][] thicknessMap(boolean img[][][], int nx, int ny, int nz) {
		float [][][] map = new float[nx][ny][nz];
		float thickness;
		int n,m;
		float min;
		
		// hyp: the object is away from the image boundary
		min = nx+ny+nz;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]) {
				// check for pointwise and edgewise connections : maybe not necessary
				// (at least, it picks-up too many cases)
				/*
				// corners: 8
				if (   (img[x][y][z]==img[x+1][y+1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x+1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x+1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y+1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x-1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x-1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y-1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x+1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x+1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y+1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x+1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x+1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else  if (   (img[x][y][z]==img[x-1][y-1][z+1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z+1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x-1][y][z+1]) )
						|| ( (img[x][y][z]!=img[x][y][z+1])
						  && (img[x][y][z]!=img[x-1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x+1][y-1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x+1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x+1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x+1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y+1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y+1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y+1][z])
						  && (img[x][y][z]!=img[x-1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x-1][y+1][z]) ) ) ) map[x][y][z] = 0.0f;
						  
				else if (   (img[x][y][z]==img[x-1][y-1][z-1]) 
					&& (   ( (img[x][y][z]!=img[x-1][y][z])
						  && (img[x][y][z]!=img[x][y-1][z-1]) )
						|| ( (img[x][y][z]!=img[x][y-1][z])
						  && (img[x][y][z]!=img[x-1][y][z-1]) )
						|| ( (img[x][y][z]!=img[x][y][z-1])
						  && (img[x][y][z]!=img[x-1][y-1][z]) ) ) ) map[x][y][z] = 0.0f;
				else
				// edges: 12
				// X,Y
				if (   (img[x][y][z]==img[x+1][y+1][z]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x+1][y-1][z]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y+1][z]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y-1][z]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else
	
				// Y,Z
				if (   (img[x][y][z]==img[x][y+1][z+1]) 
					&& (img[x][y][z]!=img[x][y][z+1])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y-1][z+1]) 
					&& (img[x][y][z]!=img[x][y][z+1])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y+1][z-1]) 
					&& (img[x][y][z]!=img[x][y][z-1])
					&& (img[x][y][z]!=img[x][y+1][z]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x][y-1][z-1]) 
					&& (img[x][y][z]!=img[x][y][z-1])
					&& (img[x][y][z]!=img[x][y-1][z]) ) map[x][y][z] = 0.0f;
				else
					
				// Z,X	
				if (   (img[x][y][z]==img[x+1][y][z+1]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y][z+1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x+1][y][z-1]) 
					&& (img[x][y][z]!=img[x+1][y][z])
					&& (img[x][y][z]!=img[x][y][z-1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y][z+1]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y][z+1]) ) map[x][y][z] = 0.0f;
				else if (   (img[x][y][z]==img[x-1][y][z-1]) 
					&& (img[x][y][z]!=img[x-1][y][z])
					&& (img[x][y][z]!=img[x][y][z-1]) ) map[x][y][z] = 0.0f;
				else {
					*/
					// no 0 thick junction: find depth along all directions
					thickness = 0.0f; map[x][y][z] = nx+ny+nz;
					// compute thickness along the x,y,z directions
					n=1;m=1;
					while ( (x+n<nx) && (img[x][y][z]==img[x+n][y][z]) ) n++;
					while ( (x-m>=0) && (img[x][y][z]==img[x-m][y][z]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (y+n<ny) && (img[x][y][z]==img[x][y+n][z]) ) n++;
					while ( (y-m>=0) && (img[x][y][z]==img[x][y-m][z]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (img[x][y][z]==img[x][y][z+n]) ) n++;
					while ( (z-m>=0) && (img[x][y][z]==img[x][y][z-m]) ) m++;
					thickness = n+m-1;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					
					// compute along the edges
					/* not needed ? 
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (img[x][y][z]==img[x+n][y+n][z]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (img[x][y][z]==img[x-m][y-m][z]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y-n>=0) && (img[x][y][z]==img[x+n][y-n][z]) ) n++;
					while ( (x-m>=0) && (y+m<ny) && (img[x][y][z]==img[x-m][y+m][z]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (z+n<nz) && (img[x][y][z]==img[x+n][y][z+n]) ) n++;
					while ( (x-m>=0) && (z-m>=0) && (img[x][y][z]==img[x-m][y][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (z-n>=0) && (img[x][y][z]==img[x+n][y][z-n]) ) n++;
					while ( (x-m>=0) && (z+m<nz) && (img[x][y][z]==img[x-m][y][z+m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (y+n<ny) && (img[x][y][z]==img[x][y+n][z+n]) ) n++;
					while ( (z-m>=0) && (y-m>=0) && (img[x][y][z]==img[x][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (z+n<nz) && (y-n>=0) && (img[x][y][z]==img[x][y-n][z+n]) ) n++;
					while ( (z-m>=0) && (y+m<ny) && (img[x][y][z]==img[x][y+m][z-m]) ) m++;
					thickness = (n+m-1)*SQR2;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					
					// along corners
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (z+n<nz) && (img[x][y][z]==img[x+n][y+n][z+n]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (z-m>=0) && (img[x][y][z]==img[x-m][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y+n<ny) && (z-n>=0) && (img[x][y][z]==img[x+n][y+n][z-n]) ) n++;
					while ( (x-m>=0) && (y-m>=0) && (z+m<nz) && (img[x][y][z]==img[x-m][y-m][z+m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x+n<nx) && (y-n>=0) && (z+n<nz) && (img[x][y][z]==img[x+n][y-n][z+n]) ) n++;
					while ( (x-m>=0) && (y+m<ny) && (z-m>=0) && (img[x][y][z]==img[x-m][y+m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					n=1;m=1;
					while ( (x-n>=0) && (y+n<ny) && (z+n<nz) && (img[x][y][z]==img[x-n][y+n][z+n]) ) n++;
					while ( (x+m<nx) && (y-m>=0) && (z-m>=0) && (img[x][y][z]==img[x+m][y-m][z-m]) ) m++;
					thickness = (n+m-1)*SQR3;
					if (thickness < map[x][y][z]) map[x][y][z] = thickness;
					*/
					
				//}
				if (map[x][y][z] < min) min = map[x][y][z];
			} else {
				map[x][y][z] = 0.0f;
			}
		}
		System.out.println("minimum thickness: "+min);
		return map;
	}
	
	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectBoundary(boolean[][][] obj, int nx, int ny, int nz) {
		return objectInsideBoundary(obj,nx,ny,nz,6);
	}
	
	/**
	 *  compute the object boundary
	 */
	public static final boolean[] objectBoundary(boolean[] obj, int nx, int ny, int nz) {
		return objectInsideBoundary(obj,nx,ny,nz,6);
	}
	
	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectInsideBoundary(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		boolean[][][]		boundary = new boolean[nx][ny][nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boundary[x][y][z] = false;
			if (obj[x][y][z]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (!obj[x+i][y+j][z+l]) boundary[x][y][z] = true;
						}
					}
				}
			}
		}
		return boundary;
	}

	/**
	 *  compute the object boundary
	 */
	public static final boolean[] objectInsideBoundary(boolean[] obj, int nx, int ny, int nz, int connectivity) {
		boolean[]		boundary = new boolean[nx*ny*nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			boundary[xyz] = false;
			if (obj[xyz]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (!obj[xyz+i+nx*j+nx*ny*l]) boundary[xyz] = true;
						}
					}
				}
			}
		}
		return boundary;
	}

	/**
	 *  compute the object boundary
	 */
	public static final boolean[][][] objectOutsideBoundary(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		boolean[][][]		boundary = new boolean[nx][ny][nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boundary[x][y][z] = false;
			if (!obj[x][y][z]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (obj[x+i][y+j][z+l]) boundary[x][y][z] = true;
						}
					}
				}
			}
		}
		return boundary;
	}
	public static final boolean[] objectOutsideBoundary(boolean[] obj, int nx, int ny, int nz, int connectivity) {
		boolean[]		boundary = new boolean[nx*ny*nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			boundary[xyz] = false;
			if (!obj[xyz]) {
				// check for boundary
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<=dist) {
						if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
							if (obj[x+i+nx*(y+j)+nx*ny*(z+l)]) boundary[xyz] = true;
						}
					}
				}
			}
		}
		return boundary;
	}
	public static final boolean[] multiObjectBoundary(byte[] label, int nx, int ny, int nz, int connectivity) {
		boolean[]		boundary = new boolean[nx*ny*nz];
		
		int dist;
		if (connectivity==6) dist=1;
		else if (connectivity==18) dist=2;
		else dist=3;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			boundary[xyz] = false;
			// check for boundary
			for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
				if (i*i+j*j+l*l<=dist) {
					if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
						if (label[x+i+nx*(y+j)+nx*ny*(z+l)]!=label[xyz]) boundary[xyz] = true;
					}
				}
			}
		}
		return boundary;
	}

	/**
	 *  compute a simplistic mean curvature function on the object boundary
	 */
	public static final float[][][] objectCurvature(boolean[][][] obj, int nx, int ny, int nz) {
		float   	num, den;
		float		u,xp,xm,yp,ym,zp,zm;
		float		xpyp,xpym,xpzp,xpzm,xmyp,xmym,xmzp,xmzm,ypzp,ypzm,ymzp,ymzm;
		float		xpypzp,xpypzm,xpymzp,xpymzm,xmypzp,xmypzm,xmymzp,xmymzm;
		float   	ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		float		val;
		float[][][]		curv = new float[nx][ny][nz];
		float[][][]		level = new float[nx][ny][nz];
		int 		boundary;
		
		// create a 'levelset' map
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (obj[x][y][z]) {
				// check for boundary
				boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
					else if ( (i*i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
					else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
				}
				// attribute a levelset for the different boundaries
				if (boundary==0) level[x][y][z] = -1.0f;
				else level[x][y][z] = 0.0f;
				/*
				else if (boundary==6) level[x][y][z] = 0.0f;
				else if (boundary==18) level[x][y][z] = -SQR2+1.0f;
				else if (boundary==26) level[x][y][z] = -SQR3+1.0f;
				*/
			} else {
				level[x][y][z] = 1.0f;
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (level[x][y][z]==0) {
				// set the level set values
				u = level[x][y][z];
				
				// 6 - neighbors		
				xp = level[x+1][y][z];
				xm = level[x-1][y][z];
				yp = level[x][y+1][z];
				ym = level[x][y+1][z];
				zp = level[x][y][z+1];
				zm = level[x][y][z-1];
				
				// 18-neighbors		
				xpyp = level[x+1][y+1][z];
				xmyp = level[x-1][y+1][z];
				xpym = level[x+1][y-1][z];
				xmym = level[x-1][y-1][z];
				xpzp = level[x+1][y][z+1];
				xmzp = level[x-1][y][z+1];
				xpzm = level[x+1][y][z-1];
				xmzm = level[x-1][y][z-1];
				ypzp = level[x][y+1][z+1];
				ymzp = level[x][y-1][z+1];
				ypzm = level[x][y+1][z-1];
				ymzm = level[x][y-1][z-1];
				
				// 26-neighbors
				xpypzm = level[x+1][y+1][z-1];
				xpymzm = level[x+1][y-1][z-1];
				xmypzm = level[x-1][y+1][z-1];
				xmymzm = level[x-1][y-1][z-1];
				xpypzp = level[x+1][y+1][z+1];
				xpymzp = level[x+1][y-1][z+1];
				xmypzp = level[x-1][y+1][z+1];
				xmymzp = level[x-1][y-1][z+1];
			
				// central differences ?
				ux = 0.25f*( xp + xpyp + xpzp + xpypzp 
							-xm - xmyp - xmzp - xmypzp );
				
				uy = 0.25f*( yp + xpyp + ypzp + xpypzp 
							-ym - xpym - ymzp - xpymzp );
				
				uz = 0.25f*( zp + xpzp + ypzp + xpypzp 
							-zm - xpzm - ypzm - xpypzm );
		
				uxx = 0.0625f*( 4*xp + 2*xpyp + 2*xpzp + 2*xpym + 2*xpzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*xm + 2*xmyp + 2*xmzp + 2*xmym + 2*xmzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*yp + 2*zp + 2*ym + 2*zm + ypzp + ymzp + ypzm + ymzm);
							 
				uyy = 0.0625f*( 4*yp + 2*xpyp + 2*ypzp + 2*xmyp + 2*ypzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*ym + 2*xpym + 2*ymzp + 2*xmym + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*zp + 2*xm + 2*zm + xpzp + xmzp + xpzm + xmzm);
							 
				uzz = 0.0625f*( 4*zp + 2*xpzp + 2*ypzp + 2*xmzp + 2*ymzp + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*zm + 2*xpzm + 2*ypzm + 2*xmzm + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*yp + 2*xm + 2*ym + xpyp + xmyp + xpym + xmym);
							 
		
				uxy = 0.0625f*( 2*xpyp + 2*xmym + xpypzp + xpypzm + xmymzp + xmymzm
							  - 2*xpym - 2*xmyp - xpymzp - xpymzm - xmypzp - xmypzm );
							 
				uyz = 0.0625f*( 2*ypzp + 2*ymzm + xpypzp + xmypzp + xpymzm + xmymzm
							  - 2*ymzp - 2*ypzm - xpymzp - xmymzp - xpypzm - xmypzm );
							 
				uzx = 0.0625f*( 2*xpzp + 2*xmzm + xpypzp + xpymzp + xmypzm + xmymzm
							  - 2*xmzp - 2*xpzm - xpypzm - xpymzm - xmypzp - xmymzp );
							 
				// 3D mean curvature
				num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
				den = (float)Math.sqrt(ux*ux + uy*uy + uz*uz);
				den = 2.0f*den*den*den;
				
				if (den>0) {
					curv[x][y][z] = num/den;
				} else {
					curv[x][y][z] = 0.0f;
				}
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}

	/**
     *  compute a binary operation (and, or, xor) on two boolean images.
     */
    public static final boolean[][][] binaryOperation(boolean img1[][][], boolean img2[][][], int operator, int nx, int ny, int nz) {
		boolean[][][] res = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			res[x][y][z] = false;
			if ( (operator==AND) && (img1[x][y][z] && img2[x][y][z]) ) res[x][y][z]=true;
			else if ( (operator==OR) && (img1[x][y][z] || img2[x][y][z]) ) res[x][y][z]=true;
			else if ( (operator==XOR) && ( (img1[x][y][z] && !img2[x][y][z]) || (!img1[x][y][z] && img2[x][y][z]) ) ) res[x][y][z] = true;
		}

		return res;
	}
	
    public static final boolean[] binaryOperation(boolean img1[], boolean img2[], int operator, int nx, int ny, int nz) {
		boolean[] res = new boolean[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[xyz] = false;
			if ( (operator==AND) && (img1[xyz] && img2[xyz]) ) res[xyz]=true;
			else if ( (operator==OR) && (img1[xyz] || img2[xyz]) ) res[xyz]=true;
			else if ( (operator==XOR) && ( (img1[xyz] && !img2[xyz]) || (!img1[xyz] && img2[xyz]) ) ) res[xyz] = true;
		}

		return res;
	}
	
    public static final boolean[] binaryOperation(boolean img1[], boolean img2[], boolean img3[], int operator, int nx, int ny, int nz) {
		boolean[] res = new boolean[nx*ny*nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			res[xyz] = false;
			if ( (operator==AND) && (img1[xyz] && img2[xyz] && img3[xyz]) ) res[xyz]=true;
			else if ( (operator==OR) && (img1[xyz] || img2[xyz] || img3[xyz]) ) res[xyz]=true;
		}

		return res;
	}
	
	
	/**
     *  compute a binary operation (and, or, xor) on two boolean images.
     */
    public static final boolean[][][] inverse(boolean img[][][], int nx, int ny, int nz) {
		boolean[][][] res = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			res[x][y][z] = !img[x][y][z];
		}

		return res;
	}
	
}//ObjectProcessing class

