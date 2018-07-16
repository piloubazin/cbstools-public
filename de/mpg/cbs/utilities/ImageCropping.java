package de.mpg.cbs.utilities;

import java.io.*;
import java.util.*;
import java.lang.*;


/**
 *
 *  Utility for cropping images before processing.
 *  <p>
 *	Include dilating, cropping with rigid transform, etc.
 *	It is meant to replace the many cropping routines in 
 * 	Algorithm classes.
 *  <p>
 *  When adding a border and cropping a background value,
 *  you have to do it in that particular order (1. add border, 2. crop).
 *
 *	@version    March 2006
 *	@author     Pilou Bazin
 *
*/
public class ImageCropping {

    private     int  		nx,ny,nz;			// original dimensions
	private		int			Nx,Ny,Nz;			// with added border
    private     int         x0,xN,y0,yN,z0,zN;	// cropping coordinates
    private     int         mx,my,mz;			// cropped dimensions
	private		int			borderSize;			// added border
	private		float		maskingValue;		// value for masked areas
	private		float		maskingThreshold;	// threshold for masking
	private		String		maskingMode;		// type of masking
	
	private     boolean		verbose = true;
	private     boolean		debug = true;
	
    /**	Constructor  */
	public ImageCropping(int nx_, int ny_, int nz_) {
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		x0 = 0; y0 = 0; z0 = 0;
		xN = nx-1; yN = ny-1; zN = nz-1;
		borderSize = 0;
		Nx = nx; Ny = ny; Nz = nz;

		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        
		maskingValue = 0.0f;
		maskingThreshold = 0.0f;
	}
	
	/** access the dimensions */
	public int nx() { return nx; }
	public int ny() { return ny; }
	public int nz() { return nz; }

	public int Nx() { return Nx; }
	public int Ny() { return Ny; }
	public int Nz() { return Nz; }

	public int mx() { return mx; }
	public int my() { return my; }
	public int mz() { return mz; }
	
	public int border() { return borderSize; }
	
	public int[] getOriginalExtents() {
		int[] extent= new int[3];
		extent[0] = nx;
		extent[1] = ny;
		extent[2] = nz;
		return extent;
	}
	
	public int[] getOriginalExtents(int nt) {
		int[] extent= new int[4];
		extent[0] = nx;
		extent[1] = ny;
		extent[2] = nz;
		extent[3] = nt;
		return extent;
	}
	
	public int[] boundingBox() {
		int[] box = new int[6];
		box[0] = x0; box[1] = xN;
		box[2] = y0; box[3] = yN;
		box[4] = z0; box[5] = zN;
		return box;
	}
	
	/** set the boundaries at a defined value */
	public void setBoundingBox(int[] box) {
		x0 = box[0]; xN = box[1];
		y0 = box[2]; yN = box[3];
		z0 = box[4]; zN = box[5];
		// update the cropped dimensions
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
	}

	/** set up parameters */
	public void setMaskingValue(float v) { maskingValue = v; }
	
	/** set up parameters */
	public void setBorderSize(int b) { 
		borderSize = b; 
		
		// update dependent parameters
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;        
	}
	
	/** conversion between 1D and 3D arrays */
	public float[][][] convertArray(float[] image) {
		float[][][] 	tmp;
		
		tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[x][y][z] = image[ x + nx*y + nx*ny*z ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static float[][][] convertArray(float[] image, int nx, int ny, int nz) {
		float[][][] 	tmp;
		
		tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[x][y][z] = image[ x + nx*y + nx*ny*z ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static float[][][] convertArray(float[] image, int[] n) {
		float[][][] 	tmp;
		
		tmp = new float[n[0]][n[1]][n[2]];
		for (int x=0;x<n[0];x++) for (int y=0;y<n[1];y++) for (int z=0;z<n[2];z++)
			tmp[x][y][z] = image[ x + n[0]*y + n[0]*n[1]*z ];
		
		return tmp;
	}
	/** conversion between 1D and 3D vector arrays */
	public float[][][][] convertArray(float[] image, int nv) {
		float[][][][] 	tmp;
		
		tmp = new float[nv][nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int v=0;v<nv;v++)
			tmp[v][x][y][z] = image[ x + nx*y + nx*ny*z +nx*ny*nz*v ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays (import directly from Mipav format)*/
	/*
	public float[][][] convertFloatImage(ModelImage image) {
		float[] buffer;
		try {
			buffer = new float[nx*ny*nz];
			image.exportData(0,nx*ny*nz, buffer); // locks and releases lock
			float[][][] tmp = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
				tmp[x][y][z] = buffer[ x + nx*y + nx*ny*z ];
			
			return tmp;
		} catch (IOException error) {
			buffer = null;
			return null;
		} catch (OutOfMemoryError e){
			buffer = null;
			return null;
		}
	}
	/** conversion between 1D and 3D arrays (import directly from Mipav format)*/
	/*
	public float[][][][] convertFloatImage(ModelImage image, int nt) {
		float[] buffer;
		try {
			buffer = new float[nx*ny*nz*nt];
			image.exportData(0,nx*ny*nz*nt, buffer); // locks and releases lock
			float[][][][] tmp = new float[nt][nx][ny][nz];
			for (int t=0;t<nt;t++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
				tmp[t][x][y][z] = buffer[ x + nx*y + nx*ny*z + nx*ny*nz*t];
			
			return tmp;
		} catch (IOException error) {
			buffer = null;
			return null;
		} catch (OutOfMemoryError e){
			buffer = null;
			return null;
		}
	}
	/** conversion between 1D and 3D arrays (import directly from Mipav format)*/
	/*
	public boolean[][][][] convertBooleanImage(ModelImage image, int nt) {
		BitSet buffer;
		try {
			buffer = new BitSet(nx*ny*nz*nt);
			image.exportData(0,nx*ny*nz*nt, buffer); // locks and releases lock
			boolean[][][][] tmp = new boolean[nt][nx][ny][nz];
			for (int t=0;t<nt;t++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
				tmp[t][x][y][z] = buffer.get(x + nx*y + nx*ny*z + nx*ny*nz*t);
			
			return tmp;
		} catch (IOException error) {
			buffer = null;
			return null;
		} catch (OutOfMemoryError e){
			buffer = null;
			return null;
		}
	}
	/** conversion between 1D and 3D arrays (import directly from Mipav format)*/
	/*
	public byte[][][] convertByteImage(ModelImage image) {
		byte[] buffer;
		try {
			buffer = new byte[nx*ny*nz];
			image.exportData(0,nx*ny*nz, buffer); // locks and releases lock
			byte[][][] tmp = new byte[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
				tmp[x][y][z] = buffer[ x + nx*y + nx*ny*z ];
			
			return tmp;
		} catch (IOException error) {
			buffer = null;
			return null;
		} catch (OutOfMemoryError e){
			buffer = null;
			return null;
		}
	}
	/** conversion between 1D and 3D arrays (import directly from Mipav format)*/
	/*
	public int[][][] convertIntImage(ModelImage image) {
		int[] buffer;
		try {
			buffer = new int[nx*ny*nz];
			image.exportData(0,nx*ny*nz, buffer); // locks and releases lock
			int[][][] tmp = new int[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
				tmp[x][y][z] = buffer[ x + nx*y + nx*ny*z ];
			
			return tmp;
		} catch (IOException error) {
			buffer = null;
			return null;
		} catch (OutOfMemoryError e){
			buffer = null;
			return null;
		}
	}
	
	/** conversion between 1D and 3D arrays */
	public float[] convertArray(float[][][] image) {
		float[] 	tmp;
		
		tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static float[] convertArray(float[][][] image, int nx, int ny, int nz) {
		float[] 	tmp;
		
		tmp = new float[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static float[] convertArray(float[][][] image, int[] n) {
		float[] 	tmp;
		
		tmp = new float[n[0]*n[1]*n[2]];
		for (int x=0;x<n[0];x++) for (int y=0;y<n[1];y++) for (int z=0;z<n[2];z++)
			tmp[ x + n[0]*y + n[0]*n[1]*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public float[] convertArray(float[][][][] image, int nv) {
		float[] 	tmp;
		
		tmp = new float[nx*ny*nz*nv];
		for (int v=0;v<nv;v++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z + nx*ny*nz*v ] = image[v][x][y][z];
		
		return tmp;
	}
	
	/** conversion between 1D and 3D arrays */
	public static float[] convertArrayVector(float[][][][] image, int[] n) {
		float[] 	tmp;
		
		tmp = new float[n[0]*n[1]*n[2]*n[3]];
		for (int x=0;x<n[0];x++) for (int y=0;y<n[1];y++) for (int z=0;z<n[2];z++) for (int v=0;v<n[3];v++) 
			tmp[ x + n[0]*y + n[0]*n[1]*z + n[0]*n[1]*n[2]*v ] = image[x][y][z][v];
		
		return tmp;
	}
	
	public float[][] convertVectorArray(float[][][][] image, int nv) {
		float[][] 	tmp;
		
		tmp = new float[nv][nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int v=0;v<nv;v++) 
			tmp[ v ][ x + nx*y + nx*ny*z ] = image[x][y][z][v];
		
		return tmp;
	}
	
	public float[][] convertArrayVector(float[][][][] image, int nv) {
		float[][] 	tmp;
		
		tmp = new float[nx*ny*nz][nv];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int v=0;v<nv;v++) 
			tmp[ x + nx*y + nx*ny*z ][ v ] = image[x][y][z][v];
		
		return tmp;
	}
	
	public float[][] convertArrayVector(float[][][][][] image, int nv, int nw) {
		float[][] 	tmp;
		
		tmp = new float[nx*ny*nz][nv*nw];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int v=0;v<nv;v++)  for (int w=0;w<nw;w++) 
			tmp[ x + nx*y + nx*ny*z ][ v + nv*w ] = image[x][y][z][v][w];
		
		return tmp;
	}
	
	/** conversion between 1D and 3D arrays */
	public byte[][][] convertArray(byte[] image) {
		byte[][][] 	tmp;
		
		tmp = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[x][y][z] = image[ x + nx*y + nx*ny*z ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static byte[][][] convertArray(byte[] image, int[] n) {
		byte[][][] 	tmp;
		
		tmp = new byte[n[0]][n[1]][n[2]];
		for (int x=0;x<n[0];x++) for (int y=0;y<n[1];y++) for (int z=0;z<n[2];z++)
			tmp[x][y][z] = image[ x + n[0]*y + n[0]*n[1]*z ];
		
		return tmp;
	}
	
	/** conversion between 1D and 3D arrays */
	public short[] convertArray(short[][][] image) {
		short[] 	tmp;
		
		tmp = new short[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public byte[] convertArray(byte[][][] image) {
		byte[] 	tmp;
		
		tmp = new byte[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public boolean[] convertArray(boolean[][][] image) {
		boolean[] 	tmp;
		
		tmp = new boolean[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public boolean[] convertArray(boolean[][][][] image, int nv) {
		boolean[] 	tmp;
		
		tmp = new boolean[nx*ny*nz*nv];
		for (int v=0;v<nv;v++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z + nx*ny*nz*v ] = image[v][x][y][z];
		
		return tmp;
	}
	
	/** conversion between 1D and 3D arrays */
	public static boolean[] convertArray(boolean[][][] image, int nx, int ny, int nz) {
		boolean[] 	tmp;
		
		tmp = new boolean[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static boolean[][][] convertArray(boolean[] image, int nx, int ny, int nz) {
		boolean[][][] 	tmp;
		
		tmp = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[x][y][z] = image[ x + nx*y + nx*ny*z ];
		
		return tmp;
	}
	public static byte[] convertArray(byte[][][] image, int nx, int ny, int nz) {
		byte[] 	tmp;
		
		tmp = new byte[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static byte[] convertArray(byte[][][][] image, int nx, int ny, int nz, int nt) {
		byte[] 	tmp;
		
		tmp = new byte[nx*ny*nz*nt];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++)
			tmp[ x + nx*y + nx*ny*z + nx*ny*nz*t ] = image[x][y][z][t];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static short[] convertArray(short[][][] image, int nx, int ny, int nz) {
		short[] 	tmp;
		
		tmp = new short[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static byte[] convertArray(byte[][][] image, int[] n) {
		byte[] 	tmp;
		
		tmp = new byte[n[0]*n[1]*n[2]];
		for (int x=0;x<n[0];x++) for (int y=0;y<n[1];y++) for (int z=0;z<n[2];z++)
			tmp[ x + n[0]*y + n[0]*n[1]*z ] = image[x][y][z];
		
		return tmp;
	}
	
	/** conversion between 1D and 3D arrays */
	public int[][][] convertArray(int[] image) {
		int[][][] 	tmp;
		
		tmp = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[x][y][z] = image[ x + nx*y + nx*ny*z ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public int[][][][] convertArray(int[] image, int nv) {
		int[][][][] 	tmp;
		
		tmp = new int[nv][nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int v=0;v<nv;v++)
			tmp[v][x][y][z] = image[ x + nx*y + nx*ny*z + nx*ny*nz*v ];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public int[] convertArray(int[][][] image) {
		int[] 	tmp;
		
		tmp = new int[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}
	/** conversion between 1D and 3D arrays */
	public static int[] convertArray(int[][][] image, int nx, int ny, int nz) {
		int[] 	tmp;
		
		tmp = new int[nx*ny*nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++)
			tmp[ x + nx*y + nx*ny*z ] = image[x][y][z];
		
		return tmp;
	}

	/** add border */
	public float[][][] addBorder(float[][][] image, int border) {
		float[][][] 	tmp;
		
		borderSize = border;
		Nx = nx+2*borderSize;
		Ny = ny+2*borderSize;
		Nz = nz+2*borderSize;
		tmp = new float[Nx][Ny][Nz];
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			if ( (x<borderSize) || (x>=nx+borderSize)
				|| (y<borderSize) || (y>=ny+borderSize)
				|| (z<borderSize) || (z>=nz+borderSize) )
				tmp[x][y][z] = maskingValue;
			else
				tmp[x][y][z] = image[x-borderSize][y-borderSize][z-borderSize];
		}
		return tmp;
	}
	/** add border */
	public float[][][][] addBorder(float[][][][] image, int nv, int border) {
		float[][][][] 	tmp;
		
		borderSize = border;
		Nx = nx+2*borderSize;
		Ny = ny+2*borderSize;
		Nz = nz+2*borderSize;
		tmp = new float[nv][Nx][Ny][Nz];
		for (int v=0;v<nv;v++) for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			if ( (x<borderSize) || (x>=nx+borderSize)
				|| (y<borderSize) || (y>=ny+borderSize)
				|| (z<borderSize) || (z>=nz+borderSize) )
				tmp[v][x][y][z] = maskingValue;
			else
				tmp[v][x][y][z] = image[v][x-borderSize][y-borderSize][z-borderSize];
		}
		return tmp;
	}
	/** add border */
	public byte[][][] addBorder(byte[][][] image, int bgLabel, int border) {
		byte[][][] 	tmp;
		
		borderSize = border;
		Nx = nx+2*borderSize;
		Ny = ny+2*borderSize;
		Nz = nz+2*borderSize;
		tmp = new byte[Nx][Ny][Nz];
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			if ( (x<borderSize) || (x>=nx+borderSize)
				|| (y<borderSize) || (y>=ny+borderSize)
				|| (z<borderSize) || (z>=nz+borderSize) )
				tmp[x][y][z] = (byte)bgLabel;
			else
				tmp[x][y][z] = image[x-borderSize][y-borderSize][z-borderSize];
		}
		return tmp;
	}
	
	/** add border */
	public int[][][] addBorder(int[][][] image, int bgLabel, int border) {
		int[][][] 	tmp;
		
		borderSize = border;
		Nx = nx+2*borderSize;
		Ny = ny+2*borderSize;
		Nz = nz+2*borderSize;
		tmp = new int[Nx][Ny][Nz];
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			if ( (x<borderSize) || (x>=nx+borderSize)
				|| (y<borderSize) || (y>=ny+borderSize)
				|| (z<borderSize) || (z>=nz+borderSize) )
				tmp[x][y][z] = bgLabel;
			else
				tmp[x][y][z] = image[x-borderSize][y-borderSize][z-borderSize];
		}
		return tmp;
	}
		
		
	/** remove border */
	public float[][][] removeBorder(float[][][] image) {
		float[][][] 	tmp;
		
		tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = image[x+borderSize][y+borderSize][z+borderSize];
		}
		return tmp;
	}
		
	/** remove border */
	public float[][][][] removeBorder(float[][][][] image, int nv) {
		float[][][][] 	tmp;
		
		tmp = new float[nv][nx][ny][nz];
		for (int v=0;v<nv;v++) for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[v][x][y][z] = image[v][x+borderSize][y+borderSize][z+borderSize];
		}
		return tmp;
	}
		
	/** remove border */
	public int[][][] removeBorder(int[][][] image) {
		int[][][] 	tmp;
		
		tmp = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = image[x+borderSize][y+borderSize][z+borderSize];
		}
		return tmp;
	}

	/** remove border */
	public byte[][][] removeBorder(byte[][][] image) {
		byte[][][] 	tmp;
		
		tmp = new byte[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = image[x+borderSize][y+borderSize][z+borderSize];
		}
		return tmp;
	}

	/**
	 *	Update the border size after the image has been cropped
	 */
	public float[][][] updateBorder(float[][][] image, int border, int prevBorder) {
		float[][][] 	tmp;
		int borderDiff = border-prevBorder;
		
		borderSize = border;
		Nx = nx+2*borderSize;
		Ny = ny+2*borderSize;
		Nz = nz+2*borderSize;
		x0 = x0 + borderDiff;
		y0 = y0 + borderDiff;
		z0 = z0 + borderDiff;
		xN = xN + borderDiff;
		yN = yN + borderDiff;
		zN = zN + borderDiff;
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        tmp = new float[mx][my][mz];
		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			if ( (borderDiff>0) && 
				 ( (x<borderDiff) || (x>=mx-borderDiff)
				|| (y<borderDiff) || (y>=my-borderDiff)
				|| (z<borderDiff) || (z>=mz-borderDiff) ) )
				tmp[x][y][z] = maskingValue;
			else
				tmp[x][y][z] = image[x-borderDiff][y-borderDiff][z-borderDiff];
		}
		return tmp;
	}
    /** sets the bounding box for cropping */
    public void findCroppingBoundaries(boolean[][][] mask) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
                    if (mask[x][y][z]) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
        mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
    /** sets the bounding box for cropping, including two images */
    public void findCroppingBoundaries(boolean[][][] mask1, boolean[][][] mask2) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
                    if ( (mask1[x][y][z]) || (mask2[x][y][z]) ) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
        mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
    /** sets the bounding box for cropping, including two images */
    public void findCroppingBoundaries(boolean[] mask1, boolean[] mask2) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
				    int xyz  = x+nx*y+nx*ny*z;
                    if ( (mask1[xyz]) || (mask2[xyz]) ) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
        mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }

    /** sets the bounding box for cropping */
    public void findCroppingBoundaries(float[][][] image, float val) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
		maskingThreshold = val;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
                    if (image[x][y][z]>maskingThreshold) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }

    /** sets the bounding box for cropping */
    public void findCroppingBoundaries(float[] image, float val) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
		maskingThreshold = val;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
					int xyz = x+nx*y+nx*ny*z;
                    if (image[xyz]>maskingThreshold) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }

    /** sets the bounding box for cropping from a vector image */
    public void findCroppingBoundaries(float[][][][] image, float val, int nt) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
		maskingThreshold = val;
		for (int t=0;t<nt;t++) 
			for (int x=borderSize;x<Nx-borderSize;x++)
				for (int y=borderSize;y<Ny-borderSize;y++)
					for (int z=borderSize;z<Nz-borderSize;z++) {
						if (image[t][x][y][z]>maskingThreshold) {
							if (x < x0) x0 = x;
							if (x > xN) xN = x;
							if (y < y0) y0 = y;
							if (y > yN) yN = y;
							if (z < z0) z0 = z;
							if (z > zN) zN = z;
						}
					}
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
    /** sets the bounding box for cropping */
    public void findCroppingBoundaries(byte[][][] label, int bgLabel) {
		x0 = Nx;
        xN = 0;
        y0 = Ny;
        yN = 0;
        z0 = Nz;
        zN = 0;
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
                    if (label[x][y][z]>bgLabel) {
                        if (x < x0) x0 = x;
                        if (x > xN) xN = x;
                        if (y < y0) y0 = y;
                        if (y > yN) yN = y;
                        if (z < z0) z0 = z;
                        if (z > zN) zN = z;
                    }
                }
        mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        // debug
        if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
        return;
    }
	
    /** changes the original image boundaries so that it includes the complete cropped image */
	public int[] updateOriginalBoundaries() {
		if (debug) System.out.print("boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n");
        
		if (x0<borderSize) {
			int offset = borderSize-x0;
			if (debug) System.out.println("x0 +"+offset+"\n");
			nx += offset;
			Nx += offset;
			x0 += offset;
			xN += offset;
		}
		if (xN>=Nx-borderSize) {
			int offset = xN-Nx+borderSize;
			if (debug) System.out.println("xN +"+offset+"\n");
			nx += offset;
			Nx += offset;
		}
		if (y0<borderSize) {
			int offset = borderSize-y0;
			if (debug) System.out.println("y0 +"+offset+"\n");
			ny += offset;
			Ny += offset;
			y0 += offset;
			yN += offset;
		}
		if (yN>=Ny-borderSize) {
			int offset = yN-Ny+borderSize;
			if (debug) System.out.println("yN +"+offset+"\n");
			ny += offset;
			Ny += offset;
		}
		if (z0<borderSize) {
			int offset = borderSize-z0;
			if (debug) System.out.println("z0 +"+offset+"\n");
			nz += offset;
			Nz += offset;
			z0 += offset;
			zN += offset;
        }
		if (zN>=Nz-borderSize) {
			int offset = zN-Nz+borderSize;
			if (debug) System.out.println("zN +"+offset+"\n");
			nz += offset;
			Nz += offset;
		}
		// debug
        if (debug) System.out.print("new dimensions: ["+nx+","+ny+","+nz+"]\n");
        
		int[] dim = new int[3];
		dim[0] = nx; dim[1] = ny; dim[2] = nz;
        return dim;
    }

    /** display the changes to the original image boundaries so that it includes the complete cropped image */
	public String displayOriginalBoundaryUpdate() {
		String info = "boundary update: \n";
		if (x0<borderSize) {
			int offset = borderSize-x0;
			info += "x0 +"+offset+"\n";
		}
		if (xN>=Nx-borderSize) {
			int offset = xN-Nx+borderSize;
			info += "xN +"+offset+"\n";
		}
		if (y0<borderSize) {
			int offset = borderSize-y0;
			info += "y0 +"+offset+"\n";
		}
		if (yN>=Ny-borderSize) {
			int offset = yN-Ny+borderSize;
			info += "yN +"+offset+"\n";
		}
		if (z0<borderSize) {
			int offset = borderSize-z0;
			info += "z0 +"+offset+"\n";
        }
		if (zN>=Nz-borderSize) {
			int offset = zN-Nz+borderSize;
			info += "zN +"+offset+"\n";
		}
        return info;
    }

	/** creates a mask for cropping */
	public boolean[][][] createMask(float[][][] image, float val) {
		boolean[][][]  	mask;
        
		// uses only values over the threshold, if mask used
		maskingThreshold = val;
		mask = new boolean[nx][ny][nz];
        for (int x=borderSize;x<Nx-borderSize;x++)
			for (int y=borderSize;y<Ny-borderSize;y++)
				for (int z=borderSize;z<Nz-borderSize;z++) {
					if (image[x][y][z] <= maskingThreshold) mask[x][y][z] = false;
                    else mask[x][y][z] = true;
                }
		// remove the boundary from the computations
		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int i=0;i<borderSize;i++) {
			mask[x][y][i] = false;
			mask[x][y][Nz-i] = false;
		}
		for (int i=0;i<borderSize;i++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			mask[i][y][z] = false;
			mask[Nx-i][y][z] = false;
		}
		for (int x=0;x<Nx;x++) for (int i=0;i<borderSize;i++) for (int z=0;z<Nz;z++) {
			mask[x][i][z] = false;
			mask[x][Ny-i][z] = false;
		}

		return mask;
	} // createMask
        
    /** crop image into smaller one */
	public float[][][] cropImage(float[][][] image) {
		float[][][] smaller = new float[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public float[] cropImage(float[] image) {
		float[] smaller = new float[mx*my*mz];

		for (int xyz=0;xyz<mx*my*mz;xyz++) {
			smaller[xyz] = maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)] = image[x+nx*y+nx*ny*z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public float[] cropImage(float[] image, int nv) {
		float[] smaller = new float[mx*my*mz*nv];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) for (int n=0;n<nv;n++) {
			smaller[x+mx*y+mx*my*z+mx*my*mz*n] = maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					for (int n=0;n<nv;n++) {
						smaller[x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize) + mx*my*mz*n] = image[x+nx*y+nx*ny*z+nx*ny*nz*n];
					}
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public float[][] cropImage(float[][] image, int nv) {
		float[][] smaller = new float[nv][mx*my*mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) for (int n=0;n<nv;n++) {
			smaller[n][x+mx*y+mx*my*z] = maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					for (int n=0;n<nv;n++) {
						smaller[n][x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)] = image[n][x+nx*y+nx*ny*z];
					}
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public int[][][] cropImage(int[][][] image, int bgLabel) {
		int[][][] smaller = new int[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = bgLabel;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public short[][][] cropImage(short[][][] image) {
		short[][][] smaller = new short[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = (short)maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public byte[][][] cropImage(byte[][][] image) {
		byte[][][] smaller = new byte[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = (byte)maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public byte[][][] cropImage(byte[][][] image, int bgLabel) {
		byte[][][] smaller = new byte[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = (byte)bgLabel;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public boolean[][][] cropImage(boolean[][][] image) {
		boolean[][][] smaller = new boolean[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			smaller[x][y][z] = false;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize] = image[x][y][z];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public boolean[] cropImage(boolean[] image) {
		boolean[] smaller = new boolean[mx*my*mz];

		for (int xyz=0;xyz<mx*my*mz;xyz++) {
			smaller[xyz] = false;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
                    int xyz = x+nx*y+nx*ny*z;
                    int xyzs = x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize);
					smaller[xyzs] = image[xyz];
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public float[][][][] cropImage(float[][][][] image, int nv) {
		float[][][][] smaller = new float[mx][my][mz][nv];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) for (int n=0;n<nv;n++) {
			smaller[x][y][z][n] = maskingValue;
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					for (int n=0;n<nv;n++) {
						smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize][n] = image[x][y][z][n];
					}
				}
			}
		}
		return smaller;
	}
    /** crop image into smaller one */
	public float[][][][][] cropImage(float[][][][][] image, int nv, int nw) {
		float[][][][][] smaller = new float[mx][my][mz][nv][nw];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			for (int n=0;n<nv;n++) for (int w=0;w<nw;w++) {
				smaller[x][y][z][n][w] = maskingValue;
			}
		}
		for (int x=x0;x<=xN;x++) {
            for (int y=y0;y<=yN;y++) {
                for (int z=z0;z<=zN;z++) {
					for (int n=0;n<nv;n++) {
						for (int w=0;w<nw;w++) {
							smaller[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize][n][w] = image[x][y][z][n][w];
						}
					}
				}
			}
		}
		return smaller;
	}

    /** recrop image using the original bounding box and the bounding differences.
	 *	The differences are >0 if the new boundix box is outside the original one. */
	public float[][][] updateCroppedImage(float[][][] image, int[] box, int[] diff) {
		x0 = box[0] - diff[0];
		xN = box[1] + diff[1];
		y0 = box[2] - diff[2];
		yN = box[3] + diff[3];
		z0 = box[4] - diff[4];
		zN = box[5] + diff[5];
		mx = xN-x0+1+2*borderSize;
        my = yN-y0+1+2*borderSize;
        mz = zN-z0+1+2*borderSize;
        float[][][] recropped = new float[mx][my][mz];

		for (int x=0;x<mx;x++) for (int y=0;y<my;y++) for (int z=0;z<mz;z++) {
			if (   (x<diff[0]+borderSize) || (x>=mx-diff[1]-borderSize) 
				|| (y<diff[2]+borderSize) || (y>=my-diff[3]-borderSize) 
				|| (z<diff[4]+borderSize) || (z>=mz-diff[5]-borderSize) )
				recropped[x][y][z] = maskingValue;
			else
				recropped[x][y][z] = image[x-diff[0]][y-diff[2]][z-diff[4]];
		}
		return recropped;
	}
	
	/** retrieve original size from smaller image */
	public float[][][] uncropImage(float[][][] image) {
		float[][][] larger = new float[Nx][Ny][Nz];

		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			larger[x][y][z] = maskingValue;
		}		
 		for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
            for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
                for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
					larger[x][y][z] = image[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	public float[] uncropImage(float[] image) {
		float[] larger = new float[Nx*Ny*Nz];

		for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) {
			larger[xyz] = maskingValue;
		}		
 		for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
            for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
                for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
					larger[x+Nx*y+Nx*Ny*z] = image[x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	public byte[] uncropImage(byte[] image) {
		byte[] larger = new byte[Nx*Ny*Nz];

		for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) {
			larger[xyz] = (byte)maskingValue;
		}		
 		for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
            for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
                for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
					larger[x+Nx*y+Nx*Ny*z] = image[x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)];
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	public float[][][][] uncropImage(float[][][][] image, int nv) {
		float[][][][] larger = new float[nv][Nx][Ny][Nz];

		for (int v=0;v<nv;v++) for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			larger[v][x][y][z] = maskingValue;
		}		
 		for (int v=0;v<nv;v++) {
			for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
				for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
					for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
						larger[v][x][y][z] = image[v][x-x0+borderSize][y-y0+borderSize][z-z0+borderSize];
					}
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	public float[][] uncropImage(float[][] image, int nv) {
		float[][] larger = new float[nv][Nx*Ny*Nz];

		for (int v=0;v<nv;v++) for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) {
			larger[v][xyz] = maskingValue;
		}		
 		for (int v=0;v<nv;v++) {
			for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
				for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
					for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
						larger[v][x+Nx*y+Nx*Ny*z] = image[v][x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)];
					}
				}
			}
		}
		return larger;
	}
	/** retrieve original size from smaller image */
	public float[] uncropImage(float[] image, int nv) {
		float[] larger = new float[nv*Nx*Ny*Nz];

		for (int v=0;v<nv;v++) for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) {
			larger[xyz+Nx*Ny*Nz*v] = maskingValue;
		}		
 		for (int v=0;v<nv;v++) {
			for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
				for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
					for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
						larger[x+Nx*y+Nx*Ny*z+Nx*Ny*Nz*v] = image[x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize) + mx*my*mz*v];
					}
				}
			}
		}
		return larger;
	}
	public byte[][][] uncropImage(byte[][][] image) {
		byte[][][] larger = new byte[Nx][Ny][Nz];

		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			larger[x][y][z] = (byte)maskingValue;
		}		
 		for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
            for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
                for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
					larger[x][y][z] = image[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize];
				}
			}
		}
		return larger;
	}
	public int[][][] uncropImage(int[][][] image, int bgValue) {
		int[][][] larger = new int[Nx][Ny][Nz];

		for (int x=0;x<Nx;x++) for (int y=0;y<Ny;y++) for (int z=0;z<Nz;z++) {
			larger[x][y][z] = bgValue;
		}		
 		for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
            for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
                for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
					larger[x][y][z] = image[x-x0+borderSize][y-y0+borderSize][z-z0+borderSize];
				}
			}
		}
		return larger;
	}
	public int[][][] uncropImage(int[][][] image) {
		return uncropImage(image, 0);
	}
	/** retrieve original size from smaller image */
	public float[] uncropAndBuffer(float[][] image, int nv) {
		float[] larger = new float[nv*Nx*Ny*Nz];

		for (int xyz=0;xyz<nv*Nx*Ny*Nz;xyz++) {
			larger[xyz] = maskingValue;
		}		
 		for (int v=0;v<nv;v++) {
			for (int x=Numerics.max(x0-borderSize,0);x<=Numerics.min(xN+borderSize,Nx-1);x++) {
				for (int y=Numerics.max(y0-borderSize,0);y<=Numerics.min(yN+borderSize,Ny-1);y++) {
					for (int z=Numerics.max(z0-borderSize,0);z<=Numerics.min(zN+borderSize,Nz-1);z++) {
						larger[x+Nx*y+Nx*Ny*z+Nx*Ny*Nz*v] = image[v][x-x0+borderSize + mx*(y-y0+borderSize) + mx*my*(z-z0+borderSize)];
					}
				}
			}
		}
		return larger;
	}

	/** 
	* passes a rigid transform matrix from the cropped image space
	* back into the original image space.
	* The matrices are in Mipav's convention.
	*/
	public double[][] uncropTransformMatrix(double[][] mat) {
		double[][] larger = new double[4][4];

		for (int i=0;i<4;i++) for (int j=0;j<4;j++) 
            larger[i][j] = mat[i][j];
        // add the offset: (Id-R^T)P0
		larger[0][3] += ( (x0-borderSize) - mat[0][0]*(x0-borderSize) - mat[0][1]*(y0-borderSize) - mat[0][2]*(z0-borderSize));
		larger[1][3] += ( (y0-borderSize) - mat[1][0]*(x0-borderSize) - mat[1][1]*(y0-borderSize) - mat[1][2]*(z0-borderSize));
		larger[2][3] += ( (z0-borderSize) - mat[2][0]*(x0-borderSize) - mat[2][1]*(y0-borderSize) - mat[2][2]*(z0-borderSize));
               
		return larger;
	}
	public double[][] uncropTransformMatrix(float[][] mat) {
		double[][] larger = new double[4][4];

		if (mat.length==3) {
			for (int i=0;i<3;i++) for (int j=0;j<4;j++) 
				larger[i][j] = mat[i][j];
			for (int j=0;j<3;j++) larger[3][j] = 0.0f;
				larger[3][3] = 1.0f;
        } else {
			for (int i=0;i<4;i++) for (int j=0;j<4;j++) 
				larger[i][j] = mat[i][j];
		}
        // add the offset: (Id-R^T)P0
		larger[0][3] += ( (x0-borderSize) - mat[0][0]*(x0-borderSize) - mat[0][1]*(y0-borderSize) - mat[0][2]*(z0-borderSize));
		larger[1][3] += ( (y0-borderSize) - mat[1][0]*(x0-borderSize) - mat[1][1]*(y0-borderSize) - mat[1][2]*(z0-borderSize));
		larger[2][3] += ( (z0-borderSize) - mat[2][0]*(x0-borderSize) - mat[2][1]*(y0-borderSize) - mat[2][2]*(z0-borderSize));
               
		return larger;
	}
	public double[][] uncropTransformMatrix(double[][] mat, float rx, float ry, float rz) {
		double[][] larger = new double[4][4];

		if (mat.length==3) {
			for (int i=0;i<3;i++) for (int j=0;j<4;j++) 
				larger[i][j] = mat[i][j];
			for (int j=0;j<3;j++) larger[3][j] = 0.0f;
				larger[3][3] = 1.0f;
        } else {
			for (int i=0;i<4;i++) for (int j=0;j<4;j++) 
				larger[i][j] = mat[i][j];
		}
        // add the offset: (Id-R^T)P0
		larger[0][3] += ( (x0-borderSize)*rx - mat[0][0]*(x0-borderSize)*rx - mat[0][1]*(y0-borderSize)*ry - mat[0][2]*(z0-borderSize)*rz);
		larger[1][3] += ( (y0-borderSize)*ry - mat[1][0]*(x0-borderSize)*rx - mat[1][1]*(y0-borderSize)*ry - mat[1][2]*(z0-borderSize)*rz);
		larger[2][3] += ( (z0-borderSize)*rz - mat[2][0]*(x0-borderSize)*rx - mat[2][1]*(y0-borderSize)*ry - mat[2][2]*(z0-borderSize)*rz);
               
		return larger;
	}
	public String displayParameters() {
		String line = "cropping parameters: \n";
		line +=  "original size: ["+nx+","+ny+","+nz+"]\n";
		line +=  "border: "+borderSize+"\n";
		line +=  "extended size: ["+Nx+","+Ny+","+Nz+"]\n";
		line +=  "boundaries: ["+x0+","+xN+"] ["+y0+","+yN+"] ["+z0+","+zN+"]\n";
		line +=  "cropped size: ["+mx+","+my+","+mz+"]\n";
		return line;
	}
	
	public final int[][] mapColorsToLabels(int[] buffer) {
		ArrayList<int[]> map = new ArrayList();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int[] color = new int[4];
			color[0] = buffer[0+4*x+4*nx*y+4*nx*ny*z];
			color[1] = buffer[1+4*x+4*nx*y+4*nx*ny*z];
			color[2] = buffer[2+4*x+4*nx*y+4*nx*ny*z];
			color[3] = buffer[3+4*x+4*nx*y+4*nx*ny*z];
		
			boolean isFound = false;
			for (int n=0;n<map.size();n++) {
				if ( color[0]==map.get(n)[0] && color[1]==map.get(n)[1] && color[2]==map.get(n)[2] && color[3]==map.get(n)[3] )
					isFound = true;
			}
			if (!isFound) map.add(color);
		}
		int[][] finalmap = new int[map.size()][4];
		for (int n=0;n<map.size();n++) for (int i=0;i<4;i++) {
			finalmap[n][i] = map.get(n)[i];
		}
		return finalmap;
	}
	
	public final int[][][] convertColorToLabelImage(int[] buffer, int[][] map) {
		int[][][] img = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			img[x][y][z] = 0;
			
			int[] color = new int[4];
			color[0] = buffer[0+4*x+4*nx*y+4*nx*ny*z];
			color[1] = buffer[1+4*x+4*nx*y+4*nx*ny*z];
			color[2] = buffer[2+4*x+4*nx*y+4*nx*ny*z];
			color[3] = buffer[3+4*x+4*nx*y+4*nx*ny*z];
		
			for (int n=0;n<map.length;n++) {
				if ( color[0]==map[n][0] && color[1]==map[n][1] && color[2]==map[n][2] && color[3]==map[n][3] )
					img[x][y][z] = n+1;
			}
		}
		return img;
	}
}