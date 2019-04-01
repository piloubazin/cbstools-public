package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version   Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectTransforms {
	
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
	*  compute a skeleton using a fast marching method
	*/
	public static final boolean[][][] simpleSkeleton(boolean[][][] obj, int nx, int ny, int nz, int connectivity) {
		float			val=0.0f, ngb=0.0f;
		int				x,y,z;
		int[]			vec = new int[3];
		boolean[][][]	skeleton, critical;
		float[][][]		dist;
		int max;
		BinaryTree		tree;
		CriticalPointLUT lut;
		float epsilon = 1e-3f;
		float mindist = 2.0f;
		
		/*
		dist = fastMarchingDistance(obj,nx,ny,nz,connectivity);
		float maxdist = ImageFunctions.maximum(dist,nx,ny,nz);
		*/
		
		tree = new BinaryTree(nx*ny*nz, 3, BinaryTree.MINTREE, BinaryTree.ADAPTATIVE);
		
		if (connectivity==6) lut = new CriticalPointLUT("critical626LUT.raw.gz",200);
		else if (connectivity==18) lut = new CriticalPointLUT("critical186LUT.raw.gz",200);
		else lut = new CriticalPointLUT("critical266LUT.raw.gz",200);
		if (!lut.loadCompressedPattern()) {
			System.out.println("Problem loading the algorithm's LUT from: "+lut.getFilename());
			return null;
		}
		
		/*
		if (connectivity==6) max=2;
		else if (connectivity==18) max=3; 
		else max=4;
		*/
		max=2;
		
		// no processing outside the object
		skeleton = new boolean[nx][ny][nz];
		critical = new boolean[nx][ny][nz];
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			skeleton[x][y][z] = obj[x][y][z];
			critical[x][y][z] = false;
		}
		
		// pre-compute the distance function: must use same algorithm
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			if (obj[x][y][z]) {
				// check for inside boundary
				int boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (i*i+j*j+l*l<max) {
							 if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
						else if ( (i* i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
						else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
					}
				}
				// set the distances
					 if (boundary==6)  val = 0.5f;
				else if (boundary==18) val = 0.5f*SQR2;
				else if (boundary==26) val = 0.5f*SQR3;
				// record in the sorting tree
				if (boundary>0) {
					vec[0] = x; vec[1] = y; vec[2] = z;
					tree.addValue(val,vec);
				}
			}
		}
		
		// propagate
		while (tree.isNotEmpty()) {
			// get the next value
			x = tree.getFirstIndex(0);
			y = tree.getFirstIndex(1);
			z = tree.getFirstIndex(2);
			val = tree.getFirst();
			tree.removeFirst();
			
			if (skeleton[x][y][z]) {
				
				// check the topology
				if (lut.get(lut.keyFromPattern(skeleton,x,y,z))) {
				
					// check for endpoint
					int nb = -1;
					for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
						if (skeleton[x+i][y+j][z+l]) nb++;
					}
					if (nb>1 || val<mindist) {
						
						// set the distance, update label
						//dist[x][y][z] = val;
						skeleton[x][y][z] = false;
						critical[x][y][z] = false;
		
						// find the neighbors
						for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
							if (i*i+j*j+l*l<max) {
								if (critical[x+i][y+j][z+l]) {
									ngb = val + epsilon;
									//ngb = dist[x+i][y+j][z+l] + epsilon;
									vec[0] = x+i; vec[1] = y+j; vec[2] = z+l;
									tree.addValue(ngb,vec);
								} else if (skeleton[x+i][y+j][z+l]) {
										 if (i*i+j*j+l*l==1) ngb = val + 1.0f;
									else if (i*i+j*j+l*l==2) ngb = val + SQR2;
									else if (i*i+j*j+l*l==3) ngb = val + SQR3;
									vec[0] = x+i; vec[1] = y+j; vec[2] = z+l;
									tree.addValue(ngb,vec);
								}
							}
						}
					}
				} else {
					//dist[x][y][z] = val;
					critical[x][y][z] = true;
				}
			}
		}
		// clean up
		tree.finalize(); tree = null;
		
		return skeleton;
	}//simpleSkeleton
	
	/**
	*  compute a feature transform based on Voronoi diagrams
	*	from: Maurer, Qi and Rghavan, PAMI 25:2, 2003
	*/
	public static final short[][][][] voronoiFeatureTransform(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = new short[nx][ny][nz][3];
		short[][] 		gx = new short[nx][3];
		short[][] 		gy = new short[ny][3];
		short[][] 		gz = new short[nz][3];
		
		// initialize the ft
		for (short z=0;z<nz;z++) {
			for (short y=0;y<ny;y++) {
				for (short x=0;x<nx;x++) {
					if (obj[x][y][z]) {
						ft[x][y][z][0] = x;
						ft[x][y][z][1] = y;
						ft[x][y][z][2] = z;
					} else {
						ft[x][y][z][0] = -1;
						ft[x][y][z][1] = -1;
						ft[x][y][z][2] = -1;
					}
				}
				// compute Voronoi 1D
				computePartialVoronoiDiagramX(ft, gx, y, z, nx, ny, nz);
			}
			for (short x=0;x<nx;x++) {
				// compute Voronoi 2D
				computePartialVoronoiDiagramY(ft, gy, x, z, nx, ny, nz);		
			}
		}
		for (short x=0;x<nx;x++) {
			for (short y=0;y<ny;y++) {
				// compute Voronoi 3D
				computePartialVoronoiDiagramZ(ft, gz, x, y, nx, ny, nz);
			}
		}
		gx = null;
		gy = null;
		gz = null;
		return ft;
	}
	
	private static final void computePartialVoronoiDiagramX(short[][][][] ft, short[][] g, short y, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			
			if (ft[i][y][z][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[i][y][z][0];
					g[l][1] = ft[i][y][z][1];
					g[l][2] = ft[i][y][z][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[i][y][z],xi,0)) l--;
					l++;
					g[l][0] = ft[i][y][z][0];
					g[l][1] = ft[i][y][z][1];
					g[l][2] = ft[i][y][z][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[i][y][z][0] = g[l][0];
			ft[i][y][z][1] = g[l][1];
			ft[i][y][z][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramY(short[][][][] ft, short[][] g, short x, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			
			if (ft[x][i][z][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[x][i][z][0];
					g[l][1] = ft[x][i][z][1];
					g[l][2] = ft[x][i][z][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[x][i][z],xi,1)) l--;
					l++;
					g[l][0] = ft[x][i][z][0];
					g[l][1] = ft[x][i][z][1];
					g[l][2] = ft[x][i][z][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[x][i][z][0] = g[l][0];
			ft[x][i][z][1] = g[l][1];
			ft[x][i][z][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramZ(short[][][][] ft, short[][] g, short x, short y, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			
			if (ft[x][y][i][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[x][y][i][0];
					g[l][1] = ft[x][y][i][1];
					g[l][2] = ft[x][y][i][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[x][y][i],xi,2)) l--;
					l++;
					g[l][0] = ft[x][y][i][0];
					g[l][1] = ft[x][y][i][1];
					g[l][2] = ft[x][y][i][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[x][y][i][0] = g[l][0];
			ft[x][y][i][1] = g[l][1];
			ft[x][y][i][2] = g[l][2];
		}
		return;
	}
	
	/**
	*  compute a feature transform based on Voronoi diagrams
	*	from: Maurer, Qi and Rghavan, PAMI 25:2, 2003
	*/
	public static final short[][] voronoiFeatureTransform(boolean[] obj, int nx, int ny, int nz) {
		short[][]		ft = new short[nx*ny*nz][3];
		short[][] 		gx = new short[nx][3];
		short[][] 		gy = new short[ny][3];
		short[][] 		gz = new short[nz][3];
		
		// initialize the ft
		for (short z=0;z<nz;z++) {
			for (short y=0;y<ny;y++) {
				for (short x=0;x<nx;x++) {
					int xyz = x+nx*y+nx*ny*z;
					if (obj[xyz]) {
						ft[xyz][0] = x;
						ft[xyz][1] = y;
						ft[xyz][2] = z;
					} else {
						ft[xyz][0] = -1;
						ft[xyz][1] = -1;
						ft[xyz][2] = -1;
					}
				}
				// compute Voronoi 1D
				computePartialVoronoiDiagramX(ft, gx, y, z, nx, ny, nz);
			}
			for (short x=0;x<nx;x++) {
				// compute Voronoi 2D
				computePartialVoronoiDiagramY(ft, gy, x, z, nx, ny, nz);		
			}
		}
		for (short x=0;x<nx;x++) {
			for (short y=0;y<ny;y++) {
				// compute Voronoi 3D
				computePartialVoronoiDiagramZ(ft, gz, x, y, nx, ny, nz);
			}
		}
		gx = null;
		gy = null;
		gz = null;
		return ft;
	}
	
	private static final void computePartialVoronoiDiagramX(short[][] ft, short[][] g, short y, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			int iyz = i+nx*y+nx*ny*z;
			
			if (ft[iyz][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[iyz][0];
					g[l][1] = ft[iyz][1];
					g[l][2] = ft[iyz][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[iyz],xi,0)) l--;
					l++;
					g[l][0] = ft[iyz][0];
					g[l][1] = ft[iyz][1];
					g[l][2] = ft[iyz][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nx;i++) {
			xi[0] = i;
			xi[1] = y;
			xi[2] = z;
			int iyz = i+nx*y+nx*ny*z;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[iyz][0] = g[l][0];
			ft[iyz][1] = g[l][1];
			ft[iyz][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramY(short[][] ft, short[][] g, short x, short z, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			int xiz = x+nx*i+nx*ny*z;
			
			if (ft[xiz][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[xiz][0];
					g[l][1] = ft[xiz][1];
					g[l][2] = ft[xiz][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[xiz],xi,1)) l--;
					l++;
					g[l][0] = ft[xiz][0];
					g[l][1] = ft[xiz][1];
					g[l][2] = ft[xiz][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<ny;i++) {
			xi[0] = x;
			xi[1] = i;
			xi[2] = z;
			int xiz = x+nx*i+nx*ny*z;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[xiz][0] = g[l][0];
			ft[xiz][1] = g[l][1];
			ft[xiz][2] = g[l][2];
		}
		return;
	}
	
	private static final void computePartialVoronoiDiagramZ(short[][] ft, short[][] g, short x, short y, int nx, int ny, int nz) {
		short[] xi = new short[3];
		int l=-1;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			int xyi = x+nx*y+nx*ny*i;
			
			if (ft[xyi][0]>-1) {
				if (l<1) {
					l++;
					g[l][0] = ft[xyi][0];
					g[l][1] = ft[xyi][1];
					g[l][2] = ft[xyi][2];
				} else {
					while ( (l>=1) && removeVoronoiFeature(g[l-1],g[l],ft[xyi],xi,2)) l--;
					l++;
					g[l][0] = ft[xyi][0];
					g[l][1] = ft[xyi][1];
					g[l][2] = ft[xyi][2];
				}
			}
		}
		int ns = l;
		if (ns==-1) { return; }
		l=0;
		for (short i=0;i<nz;i++) {
			xi[0] = x;
			xi[1] = y;
			xi[2] = i;
			int xyi = x+nx*y+nx*ny*i;
			
			while ( (l<ns) && ( ftDistance(xi,g[l]) > ftDistance(xi,g[l+1]) ) ) l++;
			
			ft[xyi][0] = g[l][0];
			ft[xyi][1] = g[l][1];
			ft[xyi][2] = g[l][2];
		}
		return;
	}
	
	public static final float[][][] voronoiFeatureSquaredDistance(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[][][] dist = new float[nx][ny][nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = ftDistance(pt, ft[x][y][z]);
			} else {
				dist[x][y][z] = -1;
			}
		}
		ft = null;
		
		return dist;
	}
	
	public static final float[] voronoiFeatureSquaredDistance(boolean[] obj, int nx, int ny, int nz) {
		short[][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[] dist = new float[nx*ny*nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
			if (ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = ftDistance(pt, ft[xyz]);
			} else {
				dist[xyz] = -1;
			}
		}
		ft = null;
		
		return dist;
	}
	
	public static final float[] voronoiFeatureDistance(boolean[] obj, int nx, int ny, int nz) {
		short[][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[] dist = new float[nx*ny*nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
			if (ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = (float)Math.sqrt(ftDistance(pt, ft[xyz]));
			} else {
				dist[xyz] = -1;
			}
		}
		ft = null;
		
		return dist;
	}
	
	private static final boolean removeVoronoiFeature(short[] u, short[] v, short[] w, short[] Rd, int d) {
		
		float duR = 0.0f; 
		float dvR = 0.0f; 
		float dwR = 0.0f; 
		for (int i=0;i<3;i++) if (i!=d) {
			duR += (u[i]-Rd[i])*(u[i]-Rd[i]);
			dvR += (v[i]-Rd[i])*(v[i]-Rd[i]);
			dwR += (w[i]-Rd[i])*(w[i]-Rd[i]);
		}
		return ( (w[d]-u[d])*dvR - (w[d]-v[d])*duR - (v[d]-u[d])*dwR - (w[d]-u[d])*(w[d]-v[d])*(v[d]-u[d]) > 0 );
	}
	
	private static final float ftDistance(short[] u, short[] v) {
		
		return (u[0]-v[0])*(u[0]-v[0])
			  +(u[1]-v[1])*(u[1]-v[1])
			  +(u[2]-v[2])*(u[2]-v[2]);
	}
	
	public static final float[][][] signedDistanceFunction(boolean[][][] obj, int nx, int ny, int nz) {
		short[][][][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[][][] dist = new float[nx][ny][nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (!obj[x][y][z] && ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = (float)Math.sqrt(ftDistance(pt, ft[x][y][z]))-0.5f;
			} else {
				dist[x][y][z] = 0;
			}
		}
		// same on the other part
		boolean[][][] bg = new boolean[nx][ny][nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			bg[x][y][z] = !obj[x][y][z];
		}
		ft = voronoiFeatureTransform(bg,nx,ny,nz);
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			if (!bg[x][y][z] && ft[x][y][z][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[x][y][z] = -(float)Math.sqrt(ftDistance(pt, ft[x][y][z]))+0.5f;
			} else if (!bg[x][y][z]) {
				dist[x][y][z] = 0;
			}
		}
		
		ft = null;
		
		return dist;
	}
	
	public static final float[] signedDistanceFunction(boolean[] obj, int nx, int ny, int nz) {
		short[][]	ft = voronoiFeatureTransform(obj,nx,ny,nz);
		float[] dist = new float[nx*ny*nz];
		
		// compute the distances
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (obj[xyz] && ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = (float)Math.sqrt(ftDistance(pt, ft[xyz]));
			} else {
				dist[xyz] = 0;
			}
		}
		// same on the other part
		boolean[] bg = new boolean[nx*ny*nz];
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			bg[xyz] = !obj[xyz];
		}
		ft = voronoiFeatureTransform(bg,nx,ny,nz);
		for (short x=0;x<nx;x++) for (short y=0;y<ny;y++) for (short z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (bg[xyz] && ft[xyz][0]>-1) {
				short[] pt = new short[]{x,y,z};
				dist[xyz] = -(float)Math.sqrt(ftDistance(pt, ft[xyz]));
			} else {
				dist[xyz] = 0;
			}
		}
		
		ft = null;
		
		return dist;
	}

	public static final float[] fastMarchingDistanceFunction(boolean[] object, int nx, int ny, int nz)  {
        // computation variables
        float[] levelset = new float[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] mask = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		BinaryHeap2D heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
				        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        //if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
        // initialize mask and processing domain
		float maxlvl = Numerics.max(nx/2.0f,ny/2.0f,nz/2.0f);
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
        	if (object[xyz]) levelset[xyz]=-0.5f;
        	else levelset[xyz] = 0.5f;
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (!mask[xyz]) { // inside the masked region: either fully inside or fully outside
				if (object[xyz]) levelset[xyz] = -maxlvl;
				else levelset[xyz] = maxlvl;
			}
		}
		// initialize the heap from boundaries
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	// search for boundaries
        	for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				if (object[xyzn]!=object[xyz]) {
					// we assume the levelset value is correct at the boundary
					
					// add to the heap with previous value
					byte lb = 0;
					if (object[xyzn]) lb = 1;
					heap.addValue(Numerics.abs(levelset[xyzn]),xyzn,lb);
                }
            }
        }
		//if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float maxdist = 0.0f;
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==1) levelset[xyz] = -dist;
			else levelset[xyz] = dist;
			processed[xyz]=true; // update the current level
 			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn] && !processed[xyzn]) if (object[xyzn]==object[xyz]) {
					// compute new distance based on processed neighbors for the same object
					for (byte l=0; l<6; l++) {
						nbdist[l] = -1.0f;
						nbflag[l] = false;
						int xyznb = Ngb.neighborIndex(l, xyzn, nx, ny, nz);
						// note that there is at most one value used here
						if (mask[xyznb] && processed[xyznb]) if (object[xyznb]==object[xyz]) {
							nbdist[l] = Numerics.abs(levelset[xyznb]);
							nbflag[l] = true;
						}			
					}
					float newdist = minimumMarchingDistance(nbdist, nbflag);
					
					// add to the heap
					heap.addValue(newdist,xyzn,lb);
				}
			}			
		}
		//if (debug) BasicInfo.displayMessage("done\n");		
		
       return levelset;
     }

	
	public static final float[] fastMarchingDistanceFunction(float[] levelset, int nx, int ny, int nz)  {
        // computation variables
        boolean[] object = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] mask = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		BinaryHeap2D heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
				        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        //if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
        // initialize mask and processing domain
        float maxlvl = Numerics.max(nx/2.0f,ny/2.0f,nz/2.0f);
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
        	object[xyz] = (levelset[xyz]<=0);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (!mask[xyz]) { // inside the masked region: either fully inside or fully outside
				if (object[xyz]) levelset[xyz] = -maxlvl;
				else levelset[xyz] = maxlvl;
			}
		}
		// initialize the heap from boundaries
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	processed[xyz] = false;
        	// search for boundaries
        	for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				if (object[xyzn]!=object[xyz]) {
					// we assume the levelset value is correct at the boundary
					
					// add to the heap with previous value
					byte lb = 0;
					if (object[xyzn]) lb = 1;
					heap.addValue(Numerics.abs(levelset[xyzn]),xyzn,lb);
                }
            }
        }
		//if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float maxdist = 0.0f;
        while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==1) levelset[xyz] = -dist;
			else levelset[xyz] = dist;
			processed[xyz]=true; // update the current level
 			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn] && !processed[xyzn]) if (object[xyzn]==object[xyz]) {
					// compute new distance based on processed neighbors for the same object
					for (byte l=0; l<6; l++) {
						nbdist[l] = -1.0f;
						nbflag[l] = false;
						int xyznb = Ngb.neighborIndex(l, xyzn, nx, ny, nz);
						// note that there is at most one value used here
						if (mask[xyznb] && processed[xyznb]) if (object[xyznb]==object[xyz]) {
							nbdist[l] = Numerics.abs(levelset[xyznb]);
							nbflag[l] = true;
						}			
					}
					float newdist = minimumMarchingDistance(nbdist, nbflag);
					
					// add to the heap
					heap.addValue(newdist,xyzn,lb);
				}
			}			
		}
		//if (debug) BasicInfo.displayMessage("done\n");		
		
       return levelset;
     }

	public static final float[] fastMarchingDistanceFunction(float[] levelset, float maxdist, int nx, int ny, int nz)  {
        // computation variables
        boolean[] object = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] processed = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		boolean[] mask = new boolean[nx*ny*nz]; // note: using a byte instead of boolean for the second pass
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		BinaryHeap2D heap = new BinaryHeap2D(nx*ny+ny*nz+nz*nx, BinaryHeap2D.MINTREE);
				        		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
        //if (debug) BasicInfo.displayMessage("fast marching\n");		
        heap.reset();
        // initialize mask and processing domain
        //float maxlvl = Numerics.max(nx/2.0f,ny/2.0f,nz/2.0f);
		float maxlvl = maxdist + 1.0f;
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
        	object[xyz] = (levelset[xyz]<=0);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) mask[x+nx*y+nx*ny*z] = true;
			else mask[x+nx*y+nx*ny*z] = false;
			if (!mask[xyz]) { // inside the masked region: either fully inside or fully outside
				if (object[xyz]) levelset[xyz] = -maxlvl;
				else levelset[xyz] = maxlvl;
			}
		}
		// initialize the heap from boundaries
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
        	int xyz = x+nx*y+nx*ny*z;
        	processed[xyz] = false;
        	// search for boundaries
        	for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				if (object[xyzn]!=object[xyz]) {
					// we assume the levelset value is correct at the boundary
					
					// add to the heap with previous value
					byte lb = 0;
					if (object[xyzn]) lb = 1;
					heap.addValue(Numerics.abs(levelset[xyzn]),xyzn,lb);
                }
            }
        }
		//if (debug) BasicInfo.displayMessage("init\n");		

        // grow the labels and functions
        float dist = 0.0f;
        while (heap.isNotEmpty() && dist<maxdist) {
        	// extract point with minimum distance
        	dist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	byte lb = heap.getFirstState();
			heap.removeFirst();

			// if more than nmgdm labels have been found already, this is done
			if (processed[xyz])  continue;
			
			// update the distance functions at the current level
			if (lb==1) levelset[xyz] = -dist;
			else levelset[xyz] = dist;
			processed[xyz]=true; // update the current level
 			
			// find new neighbors
			for (byte k = 0; k<6; k++) {
				int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
				
				// must be in outside the object or its processed neighborhood
				if (mask[xyzn] && !processed[xyzn]) if (object[xyzn]==object[xyz]) {
					// compute new distance based on processed neighbors for the same object
					for (byte l=0; l<6; l++) {
						nbdist[l] = -1.0f;
						nbflag[l] = false;
						int xyznb = Ngb.neighborIndex(l, xyzn, nx, ny, nz);
						// note that there is at most one value used here
						if (mask[xyznb] && processed[xyznb]) if (object[xyznb]==object[xyz]) {
							nbdist[l] = Numerics.abs(levelset[xyznb]);
							nbflag[l] = true;
						}			
					}
					float newdist = minimumMarchingDistance(nbdist, nbflag);
					
					// add to the heap
					heap.addValue(newdist,xyzn,lb);
				}
			}			
		}
		// set unprocessed values to +/-maxdist
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (!processed[xyz]) {
				if (object[xyz]) levelset[xyz] = -maxlvl;
				else levelset[xyz] = maxlvl;
			}
		}
		//if (debug) BasicInfo.displayMessage("done\n");		
		
       return levelset;
     }

	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * 
     */
    private static final float minimumMarchingDistance(float[] val, boolean[] flag) {

        float s, s2; // s = a + b +c; s2 = a*a + b*b +c*c
        float tmp;
        int count;
        s = 0;
        s2 = 0;
        count = 0;

        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.min(val[n], val[n+1]); // Take the smaller one if both are processed
				s += tmp;
				s2 += tmp*tmp;
				count++;
			} else if (flag[n]) {
				s += val[n]; // Else, take the processed one
				s2 += val[n]*val[n];
				count++;
			} else if (flag[n+1]) {
				s += val[n+1];
				s2 += val[n+1]*val[n+1];
				count++;
			}
		}
         // count must be greater than zero since there must be at least one processed pt in the neighbors
        if (count==0) System.err.print("!");
        if (s*s-count*(s2-1.0)<0) {
        	System.err.print(":");
        	tmp = 0;
        	for (int n=0;n<6;n++) if (flag[n]) tmp = Numerics.max(tmp,val[n]);
        	for (int n=0;n<6;n++) if (flag[n]) tmp = Numerics.min(tmp,val[n]);
        } else {
			tmp = (s + (float)FastMath.sqrt((double) (s*s-count*(s2-1.0))))/count;
		}
        // The larger root
        return tmp;
    }

}//ObjectTransforms class

