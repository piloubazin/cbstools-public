package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This class computes various image interpolation routines.
 *	
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ImageInterpolation {
	
	// no data: used as a library of functions
	
    // simple image functions
    
    public static final byte X = 0;
    public static final byte Y = 1;
    public static final byte Z = 2;
    public static final byte T = 3;

	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float nearestNeighborInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return 0.0f;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float nearestNeighborInterpolation(float[][][] image, float zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float nearestNeighborInterpolation(float[][][][] image, float zero, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0][c];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float nearestNeighborInterpolation(float[] image, float zero, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0 + nx*y0 + nx*ny*z0 + nx*ny*nz*c];
	}
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return 0;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][] image, byte zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[][][][] image, byte zero, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0][y0][z0][c];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static boolean nearestNeighborInterpolation(boolean[] image, boolean zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static byte nearestNeighborInterpolation(byte[] image, byte zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static short nearestNeighborInterpolation(short[] image, short zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static int nearestNeighborInterpolation(int[] image, int zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float nearestNeighborInterpolation(float[] image, float zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) {
			return zero;
        }
		
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		return image[x0+nx*y0+nx*ny*z0];
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float nearestNeighborInterpolation(float[] image, boolean[] mask, float zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) {
			return zero;
        }
		
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		if (mask[x0+nx*y0+nx*ny*z0]) return image[x0+nx*y0+nx*ny*z0];
		else return zero;
	}
	public static float nearestNeighborInterpolation(float[] image, int offset, boolean[] mask, float zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) {
			return zero;
        }
		
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		if (mask[x0+nx*y0+nx*ny*z0]) return image[x0+nx*y0+nx*ny*z0+offset];
		else return zero;
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static float[] nearestNeighborInterpolation(float[][] image, float x, float y, float z, int nx, int ny, int nz, int nd) {
		int x0,y0,z0;
		float[] val = new float[3];

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) {
			for (int d=0;d<nd;d++) val[d] = 0.0f;
            return val;
        }
		
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

		for (int d=0;d<nd;d++) val[d] = image[d][x0+nx*y0+nx*ny*z0];
		return val;
	}
	/**
	 *	linear interpolation, with given value outside the image
	 */
	public static int nearestNeighborInterpolation(int[][][] image, int zero, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-1) || (y<0) || (y>ny-1) || (z<0) || (z>nz-1) ) 
            return zero;
        
		x0 = Numerics.round(x);
		y0 = Numerics.round(y);
		z0 = Numerics.round(z);

        return image[x0][y0][z0];
	}
	
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float nearestNeighborClosestInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.round(x),0,nx-1);
		y0 = Numerics.bounded(Numerics.round(y),0,ny-1);
		z0 = Numerics.bounded(Numerics.round(z),0,nz-1);

		return image[x0][y0][z0];
	}
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float nearestNeighborClosestInterpolation(float[][][][] image, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.round(x),0,nx-1);
		y0 = Numerics.bounded(Numerics.round(y),0,ny-1);
		z0 = Numerics.bounded(Numerics.round(z),0,nz-1);

		return image[x0][y0][z0][c];
	}
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float nearestNeighborClosestInterpolation(float[] image, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.round(x),0,nx-1);
		y0 = Numerics.bounded(Numerics.round(y),0,ny-1);
		z0 = Numerics.bounded(Numerics.round(z),0,nz-1);

		return image[x0 + nx*y0 + nx*ny*z0 + nx*ny*nz*c];
	}
	/**
	 *	linear interpolation, with 0 outside the image
	 */
	public static float linearInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[][][] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);
		
		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[] image, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);
		
		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0 + nx*y0 + nx*ny*z0 + nx*ny*nz*c] 
			+ alpha*nbeta*ngamma*image[(x0+1) + nx*y0 + nx*ny*z0 + nx*ny*nz*c] 
			+ nalpha*beta*ngamma*image[x0 + nx*(y0+1) + nx*ny*z0 + nx*ny*nz*c] 
			+ nalpha*nbeta*gamma*image[x0 + nx*y0 + nx*ny*(z0+1) + nx*ny*nz*c] 
			+ alpha*beta*ngamma*image[(x0+1) + nx*(y0+1) + nx*ny*z0 + nx*ny*nz*c] 
			+ nalpha*beta*gamma*image[x0 + nx*(y0+1) + nx*ny*(z0+1) + nx*ny*nz*c] 
			+ alpha*nbeta*gamma*image[(x0+1) + nx*y0 + nx*ny*(z0+1) + nx*ny*nz*c] 
			+ alpha*beta*gamma*image[(x0+1) + nx*(y0+1) + nx*ny*(z0+1) + nx*ny*nz*c];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[][][][] image, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, use closest point
		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);
		
		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0][c] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0][c] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0][c] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1][c] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0][c] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1][c] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1][c] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1][c];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[][][][] image, float value, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0][c] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0][c] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0][c] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1][c] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0][c] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1][c] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1][c] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1][c];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        // if out of boundary, replace all with zero
        if ( (x0<0) || (x0>nx-2) || (y0<0) || (y0>ny-2) || (z0<0) || (z0>nz-2) ) 
            return value;
        
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[] image, float value, float x, float y, float z, int c, int nx, int ny, int nz, int nc) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        // if out of boundary, replace all with zero
        if ( (x0<0) || (x0>nx-2) || (y0<0) || (y0>ny-2) || (z0<0) || (z0>nz-2) ) 
            return value;
        
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny + c*nx*ny*nz;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[] image, int offset, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        // if out of boundary, replace all with zero
        if ( (x0<0) || (x0>nx-2) || (y0<0) || (y0>ny-2) || (z0<0) || (z0>nz-2) ) 
            return value;
        
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[offset+xyz0] 
			+ alpha*nbeta*ngamma*image[offset+xyz0+1] 
			+ nalpha*beta*ngamma*image[offset+xyz0+nx] 
			+ nalpha*nbeta*gamma*image[offset+xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[offset+xyz0+1+nx] 
			+ nalpha*beta*gamma*image[offset+xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[offset+xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[offset+xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearInterpolation(float[] image, float value, double x, double y, double z, int nx, int ny, int nz) {
		double alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        // if out of boundary, replace all with zero
        if ( (x0<0) || (x0>nx-2) || (y0<0) || (y0>ny-2) || (z0<0) || (z0>nz-2) ) 
            return value;
        
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return (float)val;
	}
	/**
	 *	linear interpolation, with value outside the image
	 */
	public static float linearClosestInterpolation(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

		x0 = Numerics.bounded(Numerics.floor(x),0,nx-2);
		y0 = Numerics.bounded(Numerics.floor(y),0,ny-2);
		z0 = Numerics.bounded(Numerics.floor(z),0,nz-2);

		alpha = Numerics.bounded(x - x0, 0.0f, 1.0f);
		nalpha = 1.0f - alpha;

		beta = Numerics.bounded(y - y0, 0.0f, 1.0f);
		nbeta = 1.0f - beta;

		gamma = Numerics.bounded(z - z0, 0.0f, 1.0f);
		ngamma = 1.0f - gamma;

		int xyz0 = x0 + y0*nx + z0*nx*ny;
		
		val = nalpha*nbeta*ngamma*image[xyz0] 
			+ alpha*nbeta*ngamma*image[xyz0+1] 
			+ nalpha*beta*ngamma*image[xyz0+nx] 
			+ nalpha*nbeta*gamma*image[xyz0+nx*ny] 
			+ alpha*beta*ngamma*image[xyz0+1+nx] 
			+ nalpha*beta*gamma*image[xyz0+nx+nx*ny] 
			+ alpha*nbeta*gamma*image[xyz0+1+nx*ny] 
			+ alpha*beta*gamma*image[xyz0+1+nx+nx*ny];

		return val;
	}
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float linearInterpolation(float[][][] image, boolean[][][] mask, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val,den;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		den = 0.0f;
		
		if (mask[x0][y0][z0]) {
			val += nalpha*nbeta*ngamma*image[x0][y0][z0];
			den += nalpha*nbeta*ngamma;
		}
		if (mask[x0+1][y0][z0]) {
			val += alpha*nbeta*ngamma*image[x0+1][y0][z0];
			den += alpha*nbeta*ngamma;
		}
		if (mask[x0][y0+1][z0]) {
			val += nalpha*beta*ngamma*image[x0][y0+1][z0];
			den += nalpha*beta*ngamma;
		}
		if (mask[x0][y0][z0+1]) {
			val += nalpha*nbeta*gamma*image[x0][y0][z0+1];
			den += nalpha*nbeta*gamma;
		}
		if (mask[x0+1][y0+1][z0]) {
			val += alpha*beta*ngamma*image[x0+1][y0+1][z0];
			den += alpha*beta*ngamma;
		}
		if (mask[x0][y0+1][z0+1]) {
			val += nalpha*beta*gamma*image[x0][y0+1][z0+1];
			den += nalpha*beta*gamma;
		}
		if (mask[x0+1][y0][z0+1]) {
			val += alpha*nbeta*gamma*image[x0+1][y0][z0+1];
			den += alpha*nbeta*gamma;
		}
		if (mask[x0+1][y0+1][z0+1]) {
			val += alpha*beta*gamma*image[x0+1][y0+1][z0+1];
			den += alpha*beta*gamma;
		}
		if (den>0) {
			return val/den;
		} else {
			return value;
		}
	}
	
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float linearInterpolation(float[] image, boolean[] mask, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val,den;
		int x0,y0,z0;
		int xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;
		
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		den = 0.0f;
		
		if (mask[xyz0]) {
			val += nalpha*nbeta*ngamma*image[xyz0];
			den += nalpha*nbeta*ngamma;
		}
		if (mask[xyz0+1]) {
			val += alpha*nbeta*ngamma*image[xyz0+1];
			den += alpha*nbeta*ngamma;
		}
		if (mask[xyz0+nx]) {
			val += nalpha*beta*ngamma*image[xyz0+nx];
			den += nalpha*beta*ngamma;
		}
		if (mask[xyz0+nx*ny]) {
			val += nalpha*nbeta*gamma*image[xyz0+nx*ny];
			den += nalpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx]) {
			val += alpha*beta*ngamma*image[xyz0+1+nx];
			den += alpha*beta*ngamma;
		}
		if (mask[xyz0+nx+nx*ny]) {
			val += nalpha*beta*gamma*image[xyz0+nx+nx*ny];
			den += nalpha*beta*gamma;
		}
		if (mask[xyz0+1+nx*ny]) {
			val += alpha*nbeta*gamma*image[xyz0+1+nx*ny];
			den += alpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx+nx*ny]) {
			val += alpha*beta*gamma*image[xyz0+1+nx+nx*ny];
			den += alpha*beta*gamma;
		}
		if (den>0) {
			return val/den;
		} else {
			return value;
		}
	}
	
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float linearInterpolation(float[] image, boolean[] mask, float value, double x, double y, double z, int nx, int ny, int nz) {
		double alpha,beta,gamma,nalpha,nbeta,ngamma,val,den;
		int x0,y0,z0;
		int xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;
		
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		den = 0.0f;
		
		if (mask[xyz0]) {
			val += nalpha*nbeta*ngamma*image[xyz0];
			den += nalpha*nbeta*ngamma;
		}
		if (mask[xyz0+1]) {
			val += alpha*nbeta*ngamma*image[xyz0+1];
			den += alpha*nbeta*ngamma;
		}
		if (mask[xyz0+nx]) {
			val += nalpha*beta*ngamma*image[xyz0+nx];
			den += nalpha*beta*ngamma;
		}
		if (mask[xyz0+nx*ny]) {
			val += nalpha*nbeta*gamma*image[xyz0+nx*ny];
			den += nalpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx]) {
			val += alpha*beta*ngamma*image[xyz0+1+nx];
			den += alpha*beta*ngamma;
		}
		if (mask[xyz0+nx+nx*ny]) {
			val += nalpha*beta*gamma*image[xyz0+nx+nx*ny];
			den += nalpha*beta*gamma;
		}
		if (mask[xyz0+1+nx*ny]) {
			val += alpha*nbeta*gamma*image[xyz0+1+nx*ny];
			den += alpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx+nx*ny]) {
			val += alpha*beta*gamma*image[xyz0+1+nx+nx*ny];
			den += alpha*beta*gamma;
		}
		if (den>0) {
			return (float) (val/den);
		} else {
			return (float) value;
		}
	}
	
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float linearInterpolation(float[] image, int offset, boolean[] mask, float value, double x, double y, double z, int nx, int ny, int nz) {
		double alpha,beta,gamma,nalpha,nbeta,ngamma,val,den;
		int x0,y0,z0;
		int xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;
		
		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		den = 0.0f;
		
		if (mask[xyz0]) {
			val += nalpha*nbeta*ngamma*image[xyz0+offset];
			den += nalpha*nbeta*ngamma;
		}
		if (mask[xyz0+1]) {
			val += alpha*nbeta*ngamma*image[xyz0+offset+1];
			den += alpha*nbeta*ngamma;
		}
		if (mask[xyz0+nx]) {
			val += nalpha*beta*ngamma*image[xyz0+offset+nx];
			den += nalpha*beta*ngamma;
		}
		if (mask[xyz0+nx*ny]) {
			val += nalpha*nbeta*gamma*image[xyz0+offset+nx*ny];
			den += nalpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx]) {
			val += alpha*beta*ngamma*image[xyz0+offset+1+nx];
			den += alpha*beta*ngamma;
		}
		if (mask[xyz0+nx+nx*ny]) {
			val += nalpha*beta*gamma*image[xyz0+offset+nx+nx*ny];
			den += nalpha*beta*gamma;
		}
		if (mask[xyz0+1+nx*ny]) {
			val += alpha*nbeta*gamma*image[xyz0+offset+1+nx*ny];
			den += alpha*nbeta*gamma;
		}
		if (mask[xyz0+1+nx+nx*ny]) {
			val += alpha*beta*gamma*image[xyz0+offset+1+nx+nx*ny];
			den += alpha*beta*gamma;
		}
		if (den>0) {
			return (float) (val/den);
		} else {
			return (float) value;
		}
	}
	
	/**
	 *	linear interpolation, with value outside the image and mask
	 */
	public static float minimumInterpolation(float[] image, boolean[] mask, float value, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		int xyz0;
		float UNSET = 1e12f;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return value;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;
		
		float val = UNSET;
		if (mask[xyz0] && image[xyz0]<val) val = image[xyz0];
		if (mask[xyz0+1] && image[xyz0+1]<val) val = image[xyz0+1];
		
		if (mask[xyz0+nx] && image[xyz0+nx]<val) val = image[xyz0+nx];
		if (mask[xyz0+nx*ny] && image[xyz0+nx*ny]<val) val = image[xyz0+nx*ny];
		
		if (mask[xyz0+1+nx] && image[xyz0+1+nx]<val) val = image[xyz0+1+nx];
		if (mask[xyz0+nx+nx*ny] && image[xyz0+nx+nx*ny]<val) val = image[xyz0+nx+nx*ny];
	
		if (mask[xyz0+1+nx*ny] && image[xyz0+1+nx*ny]<val) val = image[xyz0+1+nx*ny];
		if (mask[xyz0+1+nx+nx*ny] && image[xyz0+1+nx+nx*ny]<val) val = image[xyz0+1+nx+nx*ny];
		
		if (val!=UNSET) return val;
		else return value;
	}
	
	
	/**
	 *	linear interpolation, with +/- value outside the image
	 *	useful for interpolating gradients (suppose the gradient is +value away from the image)
	 */
	public static float linearGradientInterpolation(float[][][] image, float value, int dim, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with consistent value
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) {
			// dim corresponds to X,Y or Z
				 if ( (dim==0) && (x<0) )    return -value;
			else if ( (dim==0) && (x>nx-2) ) return value;
			else if ( (dim==1) && (y<0) )    return -value;
			else if ( (dim==1) && (y>ny-2) ) return value;
			else if ( (dim==2) && (z<0) )    return -value;
			else if ( (dim==2) && (z>nz-2) ) return value;
			else return 0.0f;
        }
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 *	linear interpolation, with distance to boundary outside image
	 */
	public static float linearDistanceInterpolation(float[][][] image, float value, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		int x0,y0,z0;
		float dist = 0.0f;

        // if out of boundary, replace all with consistent value
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) {
			// dim corresponds to X,Y or Z
			if (x<0) dist += (0-x)*value;
			else if (x>nx-2) dist+= (x-nx+2)*value; 
			if (y<0) dist += (0-y)*value;
			else if (y>ny-2) dist+= (y-ny+2)*value; 
			if (z<0) dist += (0-z)*value;
			else if (z>nz-2) dist+= (z-nz+2)*value; 
			return dist;
        }
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*nbeta*ngamma*image[x0][y0][z0] 
			+ alpha*nbeta*ngamma*image[x0+1][y0][z0] 
			+ nalpha*beta*ngamma*image[x0][y0+1][z0] 
			+ nalpha*nbeta*gamma*image[x0][y0][z0+1] 
			+ alpha*beta*ngamma*image[x0+1][y0+1][z0] 
			+ nalpha*beta*gamma*image[x0][y0+1][z0+1] 
			+ alpha*nbeta*gamma*image[x0+1][y0][z0+1] 
			+ alpha*beta*gamma*image[x0+1][y0+1][z0+1];

		return val;
	}
	/**
	 * interpolate a computation mask:
	 * where there is some data, the mask is set to true
	 */
	public static boolean maskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		return inclusiveMaskInterpolation(image, x, y, z, nx, ny, nz);
	}
		
	/**
	 * interpolate a computation mask:
	 * the mask is true if the linear interpolant is above 0.5
	 */
	public static boolean linearMaskInterpolation(boolean[][][] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		if (mask[x0][y0][z0])       val +=nalpha*nbeta*ngamma;
		if (mask[x0+1][y0][z0])     val +=alpha*nbeta*ngamma;
		if (mask[x0][y0+1][z0])     val +=nalpha*beta*ngamma;
		if (mask[x0][y0][z0+1])     val +=nalpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0])   val +=alpha*beta*ngamma;
		if (mask[x0][y0+1][z0+1])   val +=nalpha*beta*gamma;
		if (mask[x0+1][y0][z0+1])   val +=alpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0+1]) val +=alpha*beta*gamma;

		if (val > 0.5) return true;
		else return false;
	}

	/**
	 * interpolate a computation mask:
	 * the mask is true if the linear interpolant is above 0.5
	 */
	public static boolean linearMaskInterpolation(boolean[] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		int xyz0;
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		if (mask[xyz0])       		val +=nalpha*nbeta*ngamma;
		if (mask[xyz0+1])     		val +=alpha*nbeta*ngamma;
		if (mask[xyz0+nx])     		val +=nalpha*beta*ngamma;
		if (mask[xyz0+nx*ny])     	val +=nalpha*nbeta*gamma;
		if (mask[xyz0+1+nx])   		val +=alpha*beta*ngamma;
		if (mask[xyz0+nx+nx*ny])   val +=nalpha*beta*gamma;
		if (mask[xyz0+1+nx*ny])   	val +=alpha*nbeta*gamma;
		if (mask[xyz0+1+nx+nx*ny]) val +=alpha*beta*gamma;

		if (val > 0.5) return true;
		else return false;
	}

	/**
	 * interpolate a computation mask:
	 * the mask is true if any of the interpolant values is true
	 */
	public static boolean minimumMaskInterpolation(boolean[] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		int xyz0;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		xyz0 = x0 + nx*y0 + nx*ny*z0;

		if (mask[xyz0])       		return true;
		if (mask[xyz0+1])     		return true;
		if (mask[xyz0+nx])     		return true;
		if (mask[xyz0+nx*ny])     	return true;
		if (mask[xyz0+1+nx])   		return true;
		if (mask[xyz0+nx+nx*ny])   return true;
		if (mask[xyz0+1+nx*ny])   	return true;
		if (mask[xyz0+1+nx+nx*ny]) return true;

		return false;
	}
	
	/**
	 * interpolate a computation mask:
	 * returns the linearly interpolated value in [0,1]
	 */
	public static float linearMaskInterpolationValue(boolean[][][] mask, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;
		float alpha,beta,gamma,nalpha,nbeta,ngamma,val;
		
        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = 0.0f;
		if (mask[x0][y0][z0])       val +=nalpha*nbeta*ngamma;
		if (mask[x0+1][y0][z0])     val +=alpha*nbeta*ngamma;
		if (mask[x0][y0+1][z0])     val +=nalpha*beta*ngamma;
		if (mask[x0][y0][z0+1])     val +=nalpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0])   val +=alpha*beta*ngamma;
		if (mask[x0][y0+1][z0+1])   val +=nalpha*beta*gamma;
		if (mask[x0+1][y0][z0+1])   val +=alpha*nbeta*gamma;
		if (mask[x0+1][y0+1][z0+1]) val +=alpha*beta*gamma;

		return val;
	}

	/**
	 *  interpolate a computation mask:
	 *  where there is some data, the mask is set to true 
	 *  (generates a mask that contain the transformed image)
	 */
	public static boolean inclusiveMaskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        if (image[x0][y0][z0]) return true;
		else if (image[x0+1][y0][z0]) return true; 
        else if (image[x0][y0+1][z0]) return true;
        else if (image[x0][y0][z0+1]) return true; 
        else if (image[x0+1][y0+1][z0]) return true;
        else if (image[x0][y0+1][z0+1]) return true;
        else if (image[x0+1][y0][z0+1]) return true; 
        else if (image[x0+1][y0+1][z0+1]) return true; 

		return false;
	}
	/**
	 *  interpolate a computation mask:
	 *  where there is no data, the mask is set to false
	 *  (generates a mask that is contained in the transformed image)
	 */
	public static boolean exclusiveMaskInterpolation(boolean[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return false;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

        if (!image[x0][y0][z0]) return false;
		else if (!image[x0+1][y0][z0]) return false; 
        else if (!image[x0][y0+1][z0]) return false;
        else if (!image[x0][y0][z0+1]) return false; 
        else if (!image[x0+1][y0+1][z0]) return false;
        else if (!image[x0][y0+1][z0+1]) return false;
        else if (!image[x0+1][y0][z0+1]) return false; 
        else if (!image[x0+1][y0+1][z0+1]) return false; 

		return true;
	}

	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationXderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float beta,gamma,nbeta,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nbeta*ngamma*(image[x0+1][y0][z0]-image[x0][y0][z0]) 
			+ beta*ngamma*(image[x0+1][y0+1][z0]-image[x0][y0+1][z0])
			+ nbeta*gamma*(image[x0+1][y0][z0+1]-image[x0][y0][z0+1])
			+ beta*gamma*(image[x0+1][y0+1][z0+1]-image[x0][y0+1][z0+1]);

		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationYderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,gamma,nalpha,ngamma,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*ngamma*(image[x0][y0+1][z0]-image[x0][y0][z0])
			+ alpha*ngamma*(image[x0+1][y0+1][z0]-image[x0+1][y0][z0])
			+ nalpha*gamma*(image[x0][y0+1][z0+1]-image[x0][y0][z0+1])
			+ alpha*gamma*(image[x0+1][y0+1][z0+1]-image[x0+1][y0][z0+1]);
		
		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationZderivative(float[][][] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,nalpha,nbeta,val;
		int x0,y0,z0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		val = nalpha*nbeta*(image[x0][y0][z0+1]-image[x0][y0][z0])
			+ alpha*nbeta*(image[x0+1][y0][z0+1]-image[x0+1][y0][z0])
			+ nalpha*beta*(image[x0][y0+1][z0+1]-image[x0][y0+1][z0])
			+ alpha*beta*(image[x0+1][y0+1][z0+1]-image[x0+1][y0+1][z0]);
		
		return val;
	}
    
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationXderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float beta,gamma,nbeta,ngamma,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);
		
		xyz0 = x0 + nx*y0 + nx*ny*z0;

		beta = y - y0;
		nbeta = 1.0f - beta;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nbeta*ngamma*(image[xyz0+1]-image[xyz0]) 
			+ beta*ngamma*(image[xyz0+1+nx]-image[xyz0+nx])
			+ nbeta*gamma*(image[xyz0+1+nx*ny]-image[xyz0+nx*ny])
			+ beta*gamma*(image[xyz0+1+nx+nx*ny]-image[xyz0+nx+nx*ny]);

		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationYderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,gamma,nalpha,ngamma,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		xyz0 = x0 + nx*y0 + nx*ny*z0;

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		gamma = z - z0;
		ngamma = 1.0f - gamma;

		val = nalpha*ngamma*(image[xyz0+nx]-image[xyz0])
			+ alpha*ngamma*(image[xyz0+1+nx]-image[xyz0+1])
			+ nalpha*gamma*(image[xyz0+nx+nx*ny]-image[xyz0+nx*ny])
			+ alpha*gamma*(image[xyz0+1+nx+nx*ny]-image[xyz0+1+nx*ny]);
		
		return val;
	}
	/**
	 *	exact derivatives of linear interpolation
	 */
	public static float linearInterpolationZderivative(float[] image, float x, float y, float z, int nx, int ny, int nz) {
		float alpha,beta,nalpha,nbeta,val;
		int x0,y0,z0,xyz0;

        // if out of boundary, replace all with zero
        if ( (x<0) || (x>nx-2) || (y<0) || (y>ny-2) || (z<0) || (z>nz-2) ) 
            return 0.0f;
        
		x0 = Numerics.floor(x);
		y0 = Numerics.floor(y);
		z0 = Numerics.floor(z);

		xyz0 = x0 + nx*y0 + nx*ny*z0;

		alpha = x - x0;
		nalpha = 1.0f - alpha;

		beta = y - y0;
		nbeta = 1.0f - beta;

		val = nalpha*nbeta*(image[xyz0+nx*ny]-image[xyz0])
			+ alpha*nbeta*(image[xyz0+1+nx*ny]-image[xyz0+1])
			+ nalpha*beta*(image[xyz0+nx+nx*ny]-image[xyz0+nx])
			+ alpha*beta*(image[xyz0+1+nx+nx*ny]-image[xyz0+1+nx]);
		
		return val;
	}
    
    /**
     *    Setup 3D cubic Lagrangian
     *    @return	the cubic Lagrangian interpolation polynomial kernel
     */
    public static float[][] setup3DCubicLagrangianInterpolation() {
        int i;
        double arg;
		float[][] wt;

        wt = new float[4][1000];
        for ( i = 0; i < 1000; i++ ) {
            arg = i / 999.0;
            wt[0][i] = (float) ( arg * ( 1.0 - arg ) * ( arg - 2.0 ) / 6.0 );
            wt[1][i] = (float) ( ( arg + 1.0 ) * ( arg - 1.0 ) * ( arg - 2.0 ) * 0.5 );
            wt[2][i] = (float) ( arg * ( arg + 1.0 ) * ( 2.0 - arg ) * 0.5 );
            wt[3][i] = (float) ( arg * ( arg + 1.0 ) * ( arg - 1.0 ) / 6.0 );
        }
		return wt;
    }
	
    /**
     *    3D cubic Lagrangian function
	 *	  @param wt		 the interpolation kernel
     *	  @param min	 minimum threshold for interpolated image
     *	  @param max	 maximum threshold for interpolated image
     *    @param x       float point index
     *    @param y       float point index
     *    @param z       float point index
     *    @return        the cubicLagrangian3D interpolated data point
     */
    public static final float cubicLagrangianInterpolation3D(float[][][] image, float [][] wt, float min, float max, float x, float y, float z, int nx, int ny, int nz ) {

        int xbase, ybase, zbase;
        int j0, j1, j2;
        int l0, l1, l2;
        int ix, iy, iz;
        int indexX, indexY, indexZ;
        float diffX, diffY, diffZ;
        float sum;
        float ySum, zSum;
        int offset;

		
		int xD = nx-1;
		int yD = ny-1;
		int zD = nz-1;
		
        xbase = Numerics.floor(x);
        ybase = Numerics.floor(y);
        zbase = Numerics.floor(z);
        diffX = x - xbase;
        diffY = y - ybase;
        diffZ = z - zbase;
        indexX = Numerics.floor( 999.0 * diffX );
        indexY = Numerics.floor( 999.0 * diffY );
        indexZ = Numerics.floor( 999.0 * diffZ );

        // 15% - 20% faster since Math.max and Math.min are function calls
        // I also replaced the Math.abs but saw no speed improvement.
        sum = 0.0f;
        for ( ix = 0, j0 = xbase - 1; j0 <= xbase + 2; ix++, j0++ ) {
            l0 = xD;
            if ( j0 < l0 ) {
                l0 = j0;
            }
            if ( l0 < 0 ) {
                l0 = 0;
            }
            ySum = 0.0f;
            for ( iy = 0, j1 = ybase - 1; j1 <= ybase + 2; iy++, j1++ ) {
                l1 = yD;
                if ( j1 < l1 ) {
                    l1 = j1;
                }
                if ( l1 < 0 ) {
                    l1 = 0;
                }
                zSum = 0.0f;
                for ( iz = 0, j2 = zbase - 1; j2 <= zbase + 2; iz++, j2++ ) {
                    l2 = zD;
                    if ( j2 < l2 ) {
                        l2 = j2;
                    }
                    if ( l2 < 0 ) {
                        l2 = 0;
                    }
                    zSum += wt[iz][indexZ] * image[l0][l1][l2];
                } // for (iz = 0, j2 = zbase - 1; j2 <= zbase + 2;iz++, j2++)
                ySum += wt[iy][indexY] * zSum;
            } // for (iy = 0,j1 = ybase - 1; j1 <= ybase + 2;iy++, j1++)
            sum += wt[ix][indexX] * ySum;
        } // for (ix = 0,j0 = xbase - 1; j0 <= xbase + 2;ix++, j0++)
        if (sum>max) sum = max;
		if (sum<min) sum = min;
		 
        return sum;
    }
    /**
     *    3D cubic Lagrangian function
     *	  @param wt		 the interpolation kernel
     *	  @param min	 minimum threshold for interpolated image
     *	  @param max	 maximum threshold for interpolated image
     *    @param mask    boolean mask: true if the image data is to be used in the interpolation
     *    @param x       float point index
     *    @param y       float point index
     *    @param z       float point index
     *    @return        the cubicLagrangian3D interpolated data point
     */
    public static final float cubicLagrangianInterpolation3D(float[][][] image, boolean[][][] mask, float [][] wt, float min, float max, float x, float y, float z, int nx, int ny, int nz ) {

        int xbase, ybase, zbase;
        int j0, j1, j2;
        int l0, l1, l2;
        int ix, iy, iz;
        int indexX, indexY, indexZ;
        float diffX, diffY, diffZ;
        float sum,den;
        float ySum, zSum, yDen, zDen;
        int offset;

		
		int xD = nx-1;
		int yD = ny-1;
		int zD = nz-1;
		
        xbase = Numerics.floor(x);
        ybase = Numerics.floor(y);
        zbase = Numerics.floor(z);
        diffX = x - xbase;
        diffY = y - ybase;
        diffZ = z - zbase;
        indexX = Numerics.floor( 999.0 * diffX );
        indexY = Numerics.floor( 999.0 * diffY );
        indexZ = Numerics.floor( 999.0 * diffZ );

        // 15% - 20% faster since Math.max and Math.min are function calls
        // I also replaced the Math.abs but saw no speed improvement.
        sum = 0.0f;
		den = 0.0f;
        for ( ix = 0, j0 = xbase - 1; j0 <= xbase + 2; ix++, j0++ ) {
            l0 = xD;
            if ( j0 < l0 ) {
                l0 = j0;
            }
            if ( l0 < 0 ) {
                l0 = 0;
            }
            ySum = 0.0f; yDen = 0.0f;
            for ( iy = 0, j1 = ybase - 1; j1 <= ybase + 2; iy++, j1++ ) {
                l1 = yD;
                if ( j1 < l1 ) {
                    l1 = j1;
                }
                if ( l1 < 0 ) {
                    l1 = 0;
                }
                zSum = 0.0f; zDen = 0.0f;
                for ( iz = 0, j2 = zbase - 1; j2 <= zbase + 2; iz++, j2++ ) {
                    l2 = zD;
                    if ( j2 < l2 ) {
                        l2 = j2;
                    }
                    if ( l2 < 0 ) {
                        l2 = 0;
                    }
					if (mask[l0][l1][l2]) {
						zSum += wt[iz][indexZ] * image[l0][l1][l2];
						zDen += wt[iz][indexZ];
					}
                } // for (iz = 0, j2 = zbase - 1; j2 <= zbase + 2;iz++, j2++)
                ySum += wt[iy][indexY] * zSum;
				yDen += wt[iy][indexY] * zDen;
            } // for (iy = 0,j1 = ybase - 1; j1 <= ybase + 2;iy++, j1++)
            sum += wt[ix][indexX] * ySum;
            den += wt[ix][indexX] * yDen;
        } // for (ix = 0,j0 = xbase - 1; j0 <= xbase + 2;ix++, j0++)
		if (den>0) {
			sum = sum/den;
		} else {
			sum = min;
		}
        if (sum>max) sum = max;
		if (sum<min) sum = min;
		 
        return sum;
    }//cubicLagrangian3D
	
	/**
	 *	scale down by a factor
	 */
	public static float[][][] subsample(float[][][] image, int nx, int ny, int nz, int factor) {
		int nsx,nsy,nsz;
		float[][][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx][nsy][nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x][y][z] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x][y][z] += image[x*factor+i][y*factor+j][z*factor+l];
			}
			sub[x][y][z] /= scale;
		}
		return sub;
	}
	/**
	 *	scale down by a factor
	 */
	public static float[] subsample(float[] image, int nx, int ny, int nz, int factor) {
		int nsx,nsy,nsz;
		float[] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x + nsx*y + nsx*nsy*z] += image[x*factor+i + nx*(y*factor+j) + nx*ny*(z*factor+l)];
			}
			sub[x + nsx*y + nsx*nsy*z] /= scale;
		}
		return sub;
	}
	/**
	 *	scale by a value, with linear interpolation
	 */
	public static float[] scaleLinear(float[] image, boolean[] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		float[] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = linearInterpolation(image, mask, 0.0f, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	/**
	 *	scale back by a value, with linear interpolation (inverse from previous)
	 *  nx, ny, nz are the target resolutions
	 */
	public static float[] rescaleLinear(float[] image, boolean[] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		float[] sup;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sup = new float[nx*ny*nz];
			
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			sup[x + nx*y + nx*ny*z] = linearInterpolation(image, mask, 0.0f, x/factor, y/factor, z/factor, nsx, nsy, nsz);
		}
		return sup;
	}
	/**
	 *	scale by a value, with linear interpolation
	 */
	public static boolean[] scaleLinear(boolean[] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		boolean[] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new boolean[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = linearMaskInterpolation(mask, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	
	/**
	 *	scale by a value, with linear interpolation
	 */
	public static float[][][] scaleLinear(float[][][] image, boolean[][][] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		float[][][] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx][nsy][nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x][y][z] = linearInterpolation(image, mask, 0.0f, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	/**
	 *	scale back by a value, with linear interpolation (inverse from previous)
	 *  nx, ny, nz are the target resolutions
	 */
	public static float[][][] rescaleLinear(float[][][] image, boolean[][][] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		float[][][] sup;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sup = new float[nx][ny][nz];
			
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			sup[x][y][z] = linearInterpolation(image, mask, 0.0f, x/factor, y/factor, z/factor, nsx, nsy, nsz);
		}
		return sup;
	}
	/**
	 *	scale by a value, with linear interpolation
	 */
	public static boolean[][][] scaleLinear(boolean[][][] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		boolean[][][] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new boolean[nsx][nsy][nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x][y][z] = linearMaskInterpolation(mask, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	/**
	 *	scale by a value, keeping the minimum value
	 */
	public static float[] scaleMinimum(float[] image, boolean[] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		float[] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = minimumInterpolation(image, mask, 0.0f, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	/**
	 *	scale by a value, keeping the minimum values (i.e. masking only entirely masked regions)
	 */
	public static boolean[] scaleMinimum(boolean[] mask, int nx, int ny, int nz, float factor) {
		int nsx,nsy,nsz;
		boolean[] sub;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new boolean[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = minimumMaskInterpolation(mask, x*factor, y*factor, z*factor, nx, ny, nz);
		}
		return sub;
	}
	
	/**
	 *	scale down a vector image by a factor
	 */
	public static float[][][][] subsample(float[][][][] image, int nx, int ny, int nz, int nv, int factor) {
		int nsx,nsy,nsz;
		float[][][][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx][nsy][nsz][nv];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int v=0;v<nv;v++) {
			sub[x][y][z][v] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x][y][z][v] += image[x*factor+i][y*factor+j][z*factor+l][v];
			}
			sub[x][y][z][v] /= scale;
		}
		return sub;
	}
	/**
	 *	scale down a vector image by a factor
	 */
	public static float[][] subsample(float[][] image, int nx, int ny, int nz, int nv, int factor) {
		int nsx,nsy,nsz;
		float[][] sub;
		float scale = factor*factor*factor;
		
		nsx = Numerics.floor(nx/(float)factor);
		nsy = Numerics.floor(ny/(float)factor);
		nsz = Numerics.floor(nz/(float)factor);
		sub = new float[nsx*nsy*nsz][nv];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) for (int v=0;v<nv;v++) {
			sub[x+nsx*y+nsx*nsy*z][v] = 0.0f;
			for (int i=0;i<factor;i++) for (int j=0;j<factor;j++) for (int l=0;l<factor;l++) {
				sub[x+nsx*y+nsx*nsy*z][v] += image[x*factor+i + nx*(y*factor+j) + nx*ny*(z*factor+l)][v];
			}
			sub[x+nsx*y+nsx*nsy*z][v] /= scale;
		}
		return sub;
	}
	
	/**
	 *	scale down by a factor
	 */
	public static float[] subsample(float[] image, int nx, int ny, int nz, int sx, int sy, int sz) {
		int nsx = nx/sx;
		int nsy = ny/sy;
		int nsz = nz/sz;
		System.out.println("subsample new dimensions: "+nsx+", "+nsy+", "+nsz);
		float scale = sx*sy*sz;
		float[] sub = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x + nsx*y + nsx*nsy*z] = 0.0f;
			for (int i=0;i<sx;i++) for (int j=0;j<sy;j++) for (int l=0;l<sz;l++) {
				sub[x + nsx*y + nsx*nsy*z] += image[x*sx+i + nx*(y*sy+j) + nx*ny*(z*sz+l)];
			}
			sub[x + nsx*y + nsx*nsy*z] /= scale;
		}
		return sub;
	}
	
	/**
	 *	scale up by a factor
	 */
	public static float[] supersample(float[] image, int nx, int ny, int nz, int sx, int sy, int sz) {
		int nsx = nx*sx;
		int nsy = ny*sy;
		int nsz = nz*sz;
		float[] sup = new float[nsx*nsy*nsz];
			
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (int i=0;i<sx;i++) for (int j=0;j<sy;j++) for (int l=0;l<sz;l++) {
				sup[x*sx+i + nsx*(y*sy+j) + nsx*nsy*(z*sz+l)] = image[x + nx*y + nx*ny*z];
			}
		}
		return sup;
	}
	
	/**
	 *	scale up by a factor
	 */
	public static float[] supersample(float[] image, int nx, int ny, int nz, int nv, int sx, int sy, int sz) {
		int nsx = nx*sx;
		int nsy = ny*sy;
		int nsz = nz*sz;
		float[] sup = new float[nv*nsx*nsy*nsz];
			
		for (int v=0;v<nv;v++) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				for (int i=0;i<sx;i++) for (int j=0;j<sy;j++) for (int l=0;l<sz;l++) {
					sup[x*sx+i + nsx*(y*sy+j) + nsx*nsy*(z*sz+l) + nsx*nsy*nsz*v] = image[x + nx*y + nx*ny*z + nx*ny*nz*v];
				}
			}
		}
		return sup;
	}
	
	/**
	 *	scale up by a factor
	 */
	public static byte[] supersample(byte[] image, int nx, int ny, int nz, int sx, int sy, int sz) {
		int nsx = nx*sx;
		int nsy = ny*sy;
		int nsz = nz*sz;
		System.out.println("supersample new dimensions: "+nsx+", "+nsy+", "+nsz);
		byte[] sup = new byte[nsx*nsy*nsz];
			
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (int i=0;i<sx;i++) for (int j=0;j<sy;j++) for (int l=0;l<sz;l++) {
				sup[x*sx+i + nsx*(y*sy+j) + nsx*nsy*(z*sz+l)] = image[x + nx*y + nx*ny*z];
			}
		}
		return sup;
	}

	/**
	 *	scale up by a factor
	 */
	public static boolean[] supersample(boolean[] image, int nx, int ny, int nz, int sx, int sy, int sz) {
		int nsx = nx*sx;
		int nsy = ny*sy;
		int nsz = nz*sz;
		//System.out.println("supersample new dimensions: "+nsx+", "+nsy+", "+nsz);
		boolean[] sup = new boolean[nsx*nsy*nsz];
			
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			for (int i=0;i<sx;i++) for (int j=0;j<sy;j++) for (int l=0;l<sz;l++) {
				sup[x*sx+i + nsx*(y*sy+j) + nsx*nsy*(z*sz+l)] = image[x + nx*y + nx*ny*z];
			}
		}
		return sup;
	}
	
}
