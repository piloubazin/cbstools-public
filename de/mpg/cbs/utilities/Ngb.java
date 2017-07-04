package de.mpg.cbs.utilities;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class deals with various way to work with neighbor indices
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 */

public class Ngb {
	
	public	static	final	byte	pX = 0;
	public	static	final	byte	pY = 1;
	public	static	final	byte	pZ = 2;
	public	static	final	byte	mX = 3;
	public	static	final	byte	mY = 4;
	public	static	final	byte	mZ = 5;
	public	static 	final 	byte 	pXpY = 6;
	public	static 	final 	byte 	pYpZ = 7;
	public	static 	final 	byte 	pZpX = 8;
	public	static 	final 	byte 	pXmY = 9;
	public	static 	final 	byte 	pYmZ = 10;
	public	static 	final 	byte 	pZmX = 11;
	public	static 	final 	byte 	mXpY = 12;
	public	static 	final 	byte 	mYpZ = 13;
	public	static 	final 	byte 	mZpX = 14;
	public	static 	final 	byte 	mXmY = 15;
	public	static 	final 	byte 	mYmZ = 16;
	public	static 	final 	byte 	mZmX = 17;
	public	static 	final 	byte 	pXpYpZ = 18;
	public	static 	final 	byte 	pXmYmZ = 19;
	public	static 	final 	byte 	pXmYpZ = 20;
	public	static 	final 	byte 	pXpYmZ = 21;
	public	static 	final 	byte 	mXpYpZ = 22;
	public	static 	final 	byte 	mXmYmZ = 23;
	public	static 	final 	byte 	mXmYpZ = 24;
	public	static 	final 	byte 	mXpYmZ = 25;
	public	static 	final 	byte 	NGB = 26;

	public static final byte[] x = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	public static final byte[] y = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	public static final byte[] z = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	private static final	float	INVSQRT2 = (float)(1.0/FastMath.sqrt(2.0));
	private static final	float	INVSQRT3 = (float)(1.0/FastMath.sqrt(3.0));
	private static final	float	SQRT2 = (float)FastMath.sqrt(2.0);
	private static final	float	SQRT3 = (float)FastMath.sqrt(3.0);

	public static final int neighborIndex(byte d, int id, int nx, int ny, int nz) {
		switch (d) {
			case pX		: 	return id+1; 		
			case mX		:	return id-1;
			case pY		:	return id+nx;
			case mY		:	return id-nx;
			case pZ		:	return id+nx*ny;
			case mZ		:	return id-nx*ny;
			case pXpY	:	return id+1+nx;
			case mXpY	:	return id-1+nx;
			case pYpZ	:	return id+nx+nx*ny;
			case mYpZ	:	return id-nx+nx*ny;
			case pZpX	:	return id+nx*ny+1;	
			case mZpX	:	return id-nx*ny+1;
			case pXmY	:	return id+1-nx;	
			case mXmY	:	return id-1-nx;
			case pYmZ	:	return id+nx-nx*ny;
			case mYmZ	:	return id-nx-nx*ny;
			case pZmX	:	return id+nx*ny-1;
			case mZmX	:	return id-nx*ny-1;
			case pXpYpZ	:	return id+1+nx+nx*ny;
			case mXpYpZ	:	return id-1+nx+nx*ny;
			case pXmYmZ	:	return id+1-nx-nx*ny; 
			case mXmYmZ	:	return id-1-nx-nx*ny;
			case pXmYpZ	:	return id+1-nx+nx*ny;
			case mXmYpZ	:	return id-1-nx+nx*ny;
			case pXpYmZ	:	return id+1+nx-nx*ny; 
			case mXpYmZ	:	return id-1+nx-nx*ny;
			default		:	return id;
		}
	}

	public static final float neighborDistance(byte d) {
		if (d<6) return 1.0f;
		else if (d<18) return SQRT2;
		else return SQRT3;
	}

	public static final float neighborDistanceRatio(byte d) {
		if (d<6) return 1.0f;
		else if (d<18) return 1.0f/SQRT2;
		else return 1.0f/SQRT3;
	}

	public static final int neighborhoodDifferences6C(double[] sample, float[][][] img, boolean[][][] mask, int x, int y, int z) {
		int n=0;
		for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) if (mask[x+dx][y+dy][z]) {
			if (mask[x+dx][y+dy][z-1]) sample[n] = Numerics.abs(img[x+dx][y+dy][z-1]-img[x+dx][y+dy][z]); n++;
			if (mask[x+dx][y+dy][z+1]) sample[n] = Numerics.abs(img[x+dx][y+dy][z]-img[x+dx][y+dy][z+1]); n++;
		}
		for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) if (mask[x][y+dy][z+dz]) {
			if (mask[x-1][y+dy][z+dz]) sample[n] = Numerics.abs(img[x-1][y+dy][z+dz]-img[x][y+dy][z+dz]); n++;
			if (mask[x+1][y+dy][z+dz]) sample[n] = Numerics.abs(img[x][y+dy][z+dz]-img[x+1][y+dy][z+dz]); n++;
		}
		for (int dz=-1;dz<=1;dz++) for (int dx=-1;dx<=1;dx++) if (mask[x+dx][y][z+dz]) {
			if (mask[x+dx][y-1][z+dz]) sample[n] = Numerics.abs(img[x+dx][y-1][z+dz]-img[x+dx][y][z+dz]); n++;
			if (mask[x+dx][y+1][z+dz]) sample[n] = Numerics.abs(img[x+dx][y][z+dz]-img[x+dx][y+1][z+dz]); n++;
		}
		return n;	
	}

	public static final int neighborhoodDifferences18C(double[] sample, float[][][] img, boolean[][][] mask, int x, int y, int z) {
		// use the code for 6C directly
		int n=neighborhoodDifferences6C(sample, img, mask, x,y,z);
		for (int dx=-1;dx<=1;dx++) {
			if (mask[x+dx][y-1][z-1] && mask[x+dx][y][z]) sample[n] = Numerics.abs(img[x+dx][y-1][z-1]-img[x+dx][y][z]); n++;
			if (mask[x+dx][y][z-1] && mask[x+dx][y+1][z]) sample[n] = Numerics.abs(img[x+dx][y][z-1]-img[x+dx][y+1][z]); n++;
			if (mask[x+dx][y-1][z] && mask[x+dx][y][z+1]) sample[n] = Numerics.abs(img[x+dx][y-1][z]-img[x+dx][y][z+1]); n++;
			if (mask[x+dx][y][z] && mask[x+dx][y+1][z+1]) sample[n] = Numerics.abs(img[x+dx][y][z]-img[x+dx][y+1][z+1]); n++;
			
			if (mask[x+dx][y+1][z-1] && mask[x+dx][y][z]) sample[n] = Numerics.abs(img[x+dx][y+1][z-1]-img[x+dx][y][z]); n++;
			if (mask[x+dx][y][z-1] && mask[x+dx][y-1][z]) sample[n] = Numerics.abs(img[x+dx][y][z-1]-img[x+dx][y-1][z]); n++;
			if (mask[x+dx][y+1][z] && mask[x+dx][y][z+1]) sample[n] = Numerics.abs(img[x+dx][y+1][z]-img[x+dx][y][z+1]); n++;
			if (mask[x+dx][y][z] && mask[x+dx][y-1][z+1]) sample[n] = Numerics.abs(img[x+dx][y][z]-img[x+dx][y-1][z+1]); n++;
		}
		for (int dy=-1;dy<=1;dy++) {
			if (mask[x-1][y+dy][z-1] && mask[x][y+dy][z]) sample[n] = Numerics.abs(img[x-1][y+dy][z-1]-img[x][y+dy][z]); n++;
			if (mask[x][y+dy][z-1] && mask[x+1][y+dy][z]) sample[n] = Numerics.abs(img[x][y+dy][z-1]-img[x+1][y+dy][z]); n++;
			if (mask[x-1][y+dy][z] && mask[x][y+dy][z+1]) sample[n] = Numerics.abs(img[x-1][y+dy][z]-img[x][y+dy][z+1]); n++;
			if (mask[x][y+dy][z] && mask[x+1][y+dy][z+1]) sample[n] = Numerics.abs(img[x][y+dy][z]-img[x+1][y+dy][z+1]); n++;
			
			if (mask[x+1][y+dy][z-1] && mask[x][y+dy][z]) sample[n] = Numerics.abs(img[x+1][y+dy][z-1]-img[x][y+dy][z]); n++;
			if (mask[x][y+dy][z-1] && mask[x-1][y+dy][z]) sample[n] = Numerics.abs(img[x][y+dy][z-1]-img[x-1][y+dy][z]); n++;
			if (mask[x+1][y+dy][z] && mask[x][y+dy][z+1]) sample[n] = Numerics.abs(img[x+1][y+dy][z]-img[x][y+dy][z+1]); n++;
			if (mask[x][y+dy][z] && mask[x-1][y+dy][z+1]) sample[n] = Numerics.abs(img[x][y+dy][z]-img[x-1][y+dy][z+1]); n++;
		}
		for (int dz=-1;dz<=1;dz++) {
			if (mask[x-1][y-1][z+dz] && mask[x][y][z+dz]) sample[n] = Numerics.abs(img[x-1][y-1][z+dz]-img[x][y][z+dz]); n++;
			if (mask[x][y-1][z+dz] && mask[x+1][y][z+dz]) sample[n] = Numerics.abs(img[x][y-1][z+dz]-img[x+1][y][z+dz]); n++;
			if (mask[x-1][y][z+dz] && mask[x][y+1][z+dz]) sample[n] = Numerics.abs(img[x-1][y][z+dz]-img[x][y+1][z+dz]); n++;
			if (mask[x][y][z+dz] && mask[x+1][y+1][z+dz]) sample[n] = Numerics.abs(img[x][y][z+dz]-img[x+1][y+1][z+dz]); n++;
			
			if (mask[x+1][y-1][z+dz] && mask[x][y][z+dz]) sample[n] = Numerics.abs(img[x+1][y-1][z+dz]-img[x][y][z+dz]); n++;
			if (mask[x][y-1][z+dz] && mask[x-1][y][z+dz]) sample[n] = Numerics.abs(img[x][y-1][z+dz]-img[x-1][y][z+dz]); n++;
			if (mask[x+1][y][z+dz] && mask[x][y+1][z+dz]) sample[n] = Numerics.abs(img[x+1][y][z+dz]-img[x][y+1][z+dz]); n++;
			if (mask[x][y][z+dz] && mask[x-1][y+1][z+dz]) sample[n] = Numerics.abs(img[x][y][z+dz]-img[x-1][y+1][z+dz]); n++;
		}
		return n;	
	}

	public static final int neighborhoodDifferences26C(double[] sample, float[][][] img, boolean[][][]mask, int x, int y, int z) {
		int n=neighborhoodDifferences18C(sample, img, mask, x,y,z);
		if (mask[x][y][z]) {
			if (mask[x-1][y-1][z-1]) sample[n] = Numerics.abs(img[x-1][y-1][z-1]-img[x][y][z]); n++;
			if (mask[x+1][y-1][z-1]) sample[n] = Numerics.abs(img[x+1][y-1][z-1]-img[x][y][z]); n++;
			if (mask[x-1][y+1][z-1]) sample[n] = Numerics.abs(img[x-1][y+1][z-1]-img[x][y][z]); n++;
			if (mask[x-1][y-1][z+1]) sample[n] = Numerics.abs(img[x-1][y-1][z+1]-img[x][y][z]); n++;
			if (mask[x-1][y+1][z+1]) sample[n] = Numerics.abs(img[x-1][y+1][z+1]-img[x][y][z]); n++;
			if (mask[x+1][y-1][z+1]) sample[n] = Numerics.abs(img[x+1][y-1][z+1]-img[x][y][z]); n++;
			if (mask[x+1][y+1][z-1]) sample[n] = Numerics.abs(img[x+1][y+1][z-1]-img[x][y][z]); n++;
			if (mask[x+1][y+1][z+1]) sample[n] = Numerics.abs(img[x+1][y+1][z+1]-img[x][y][z]); n++;
		}	
		return n;	
	}
	
	public static final float[] directionVector(int d) {
		if (d==pX) return new float[]{1.0f, 0.0f, 0.0f};
		else if (d==pY) return new float[]{0.0f, 1.0f, 0.0f};
		else if (d==pZ) return new float[]{0.0f, 0.0f, 1.0f};
		else if (d==pXpY) return new float[]{INVSQRT2, INVSQRT2, 0.0f};
		else if (d==pYpZ) return new float[]{0.0f, INVSQRT2, INVSQRT2};
		else if (d==pZpX) return new float[]{INVSQRT2, 0.0f, INVSQRT2};
		else if (d==pXmY) return new float[]{INVSQRT2, -INVSQRT2, 0.0f};
		else if (d==pYmZ) return new float[]{0.0f, INVSQRT2, -INVSQRT2};
		else if (d==pZmX) return new float[]{-INVSQRT2, 0.0f, INVSQRT2};
		else if (d==pXpYpZ) return new float[]{INVSQRT3, INVSQRT3, INVSQRT3};
		else if (d==pXmYmZ) return new float[]{INVSQRT3, -INVSQRT3, -INVSQRT3};
		else if (d==pXmYpZ) return new float[]{INVSQRT3, -INVSQRT3, INVSQRT3};
		else if (d==pXpYmZ) return new float[]{INVSQRT3, INVSQRT3, -INVSQRT3};
		else if (d==mX) return new float[]{-1.0f, 0.0f, 0.0f};
		else if (d==mY) return new float[]{0.0f, -1.0f, 0.0f};
		else if (d==mZ) return new float[]{0.0f, 0.0f, -1.0f};
		else if (d==mXpY) return new float[]{-INVSQRT2, INVSQRT2, 0.0f};
		else if (d==mYpZ) return new float[]{0.0f, -INVSQRT2, INVSQRT2};
		else if (d==mZpX) return new float[]{INVSQRT2, 0.0f, -INVSQRT2};
		else if (d==mXmY) return new float[]{-INVSQRT2, -INVSQRT2, 0.0f};
		else if (d==mYmZ) return new float[]{0.0f, -INVSQRT2, -INVSQRT2};
		else if (d==mZmX) return new float[]{-INVSQRT2, 0.0f, -INVSQRT2};
		else if (d==mXpYpZ) return new float[]{-INVSQRT3, INVSQRT3, INVSQRT3};
		else if (d==mXmYmZ) return new float[]{-INVSQRT3, -INVSQRT3, -INVSQRT3};
		else if (d==mXmYpZ) return new float[]{-INVSQRT3, -INVSQRT3, INVSQRT3};
		else if (d==mXpYmZ) return new float[]{-INVSQRT3, INVSQRT3, -INVSQRT3};
		else return new float[]{0.0f, 0.0f, 0.0f};
	}
	public static final byte[] directionNeighbor(int d) {
		if (d==pX) return new byte[]{1, 0, 0};
		else if (d==pY) return new byte[]{0, 1, 0};
		else if (d==pZ) return new byte[]{0, 0, 1};
		else if (d==pXpY) return new byte[]{1, 1, 0};
		else if (d==pYpZ) return new byte[]{0, 1, 1};
		else if (d==pZpX) return new byte[]{1, 0, 1};
		else if (d==pXmY) return new byte[]{1, -1, 0};
		else if (d==pYmZ) return new byte[]{0, 1, -1};
		else if (d==pZmX) return new byte[]{-1, 0, 1};
		else if (d==pXpYpZ) return new byte[]{1, 1, 1};
		else if (d==pXmYmZ) return new byte[]{1, -1, -1};
		else if (d==pXmYpZ) return new byte[]{1, -1, 1};
		else if (d==pXpYmZ) return new byte[]{1, 1, -1};
		else if (d==mX) return new byte[]{-1, 0, 0};
		else if (d==mY) return new byte[]{0, -1, 0};
		else if (d==mZ) return new byte[]{0, 0, -1};
		else if (d==mXpY) return new byte[]{-1, 1, 0};
		else if (d==mYpZ) return new byte[]{0, -1, 1};
		else if (d==mZpX) return new byte[]{1, 0, -1};
		else if (d==mXmY) return new byte[]{-1, -1, 0};
		else if (d==mYmZ) return new byte[]{0, -1, -1};
		else if (d==mZmX) return new byte[]{-1, 0, -1};
		else if (d==mXpYpZ) return new byte[]{-1, 1, 1};
		else if (d==mXmYmZ) return new byte[]{-1, -1, -1};
		else if (d==mXmYpZ) return new byte[]{-1, -1, 1};
		else if (d==mXpYmZ) return new byte[]{-1, 1, -1};
		else return new byte[]{0, 0, 0};
	}
}
