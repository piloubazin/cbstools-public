package de.mpg.cbs.utilities;

import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class deals with various way to work with neighbor indices
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 */

public class Ngb2 {
	
	public	static	final	byte	pX = 0;
	public	static	final	byte	pY = 1;
	public	static	final	byte	mX = 2;
	public	static	final	byte	mY = 3;
	public	static 	final 	byte 	pXpY = 4;
	public	static 	final 	byte 	pXmY = 5;
	public	static 	final 	byte 	mXpY = 6;
	public	static 	final 	byte 	mXmY = 7;
	public	static 	final 	byte 	NGB = 8;

	public static final byte[] x = {+1,  0, -1,  0,  +1, +1, -1, -1};
	public static final byte[] y = { 0, +1,  0, -1,  +1, -1, +1, -1};
	
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
			case pXpY	:	return id+1+nx;
			case mXpY	:	return id-1+nx;
			case pXmY	:	return id+1-nx;	
			case mXmY	:	return id-1-nx;
			default		:	return id;
		}
	}

	public static final float neighborDistance(byte d) {
		if (d<6) return 1.0f;
		else return SQRT2;
	}

	public static final float neighborDistanceRatio(byte d) {
		if (d<6) return 1.0f;
		else return 1.0f/SQRT2;
	}

	public static final float[] directionVector(int d) {
		if (d==pX) return new float[]{1.0f, 0.0f};
		else if (d==pY) return new float[]{0.0f, 1.0f};
		else if (d==pXpY) return new float[]{INVSQRT2, INVSQRT2};
		else if (d==pXmY) return new float[]{INVSQRT2, -INVSQRT2};
		else if (d==mX) return new float[]{-1.0f, 0.0f};
		else if (d==mY) return new float[]{0.0f, -1.0f};
		else if (d==mXpY) return new float[]{-INVSQRT2, INVSQRT2};
		else if (d==mXmY) return new float[]{-INVSQRT2, -INVSQRT2};
		else return new float[]{0.0f, 0.0f};
	}
	public static final byte[] directionNeighbor(int d) {
		if (d==pX) return new byte[]{1, 0};
		else if (d==pY) return new byte[]{0, 1};
		else if (d==pXpY) return new byte[]{1, 1};
		else if (d==pXmY) return new byte[]{1, -1};
		else if (d==mX) return new byte[]{-1, 0};
		else if (d==mY) return new byte[]{0, -1};
		else if (d==mXpY) return new byte[]{-1, 1};
		else if (d==mXmY) return new byte[]{-1, -1};
		else return new byte[]{0, 0, 0};
	}
}
