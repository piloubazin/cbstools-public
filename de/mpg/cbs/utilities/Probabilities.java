package de.mpg.cbs.utilities;

import java.io.*;
import java.util.*;

/**
 *
 *  This class computes various basic numerical functions.
 *	<p> 
 *	Includes min, max, gaussian random numbers, special functions, 
 *	mathematical constants. This duplicates some of Java's Math functions
 *	where these have limited type support (e.g. double, but not float) or
 *	perform multiple casting operations (slow).
 *
 *	@version    May 17, 2005
 *	@author     Pierre-Louis Bazin
 */

public class Probabilities {
	
	/** mathematical constants */
	public static final float ZERO = 1e-30f;
	/** mathematical constants */
	public static final float INF = 1e+30f;
	/** mathematical constants */
	public static final float PI = 3.1416f;
	
	/**
	 *	Gaussian random numbers (from the Box-Muller method, as in NRC)
	 *  mean 0, variance 1
	 */
    public static final double randomGaussian() {
		double v1,v2,rsq;
        do {
			v1 = 2.0*Math.random()-1.0;
			v2 = 2.0*Math.random()-1.0;
			rsq = v1*v1+v2*v2;
		} while ( (rsq<=0) || (rsq>=1) );
		return v1*Math.sqrt( -2.0*Math.log(rsq)/rsq);
    }
	
	/**
	 *	Student T 3rd order random numbers
	 *  mean 0, variance given
	 *  note: the formula should be the inverse of the probability integral (wrong here!)
	 */
    public static final double randomT3(double sigma) {
		double v1,v2,rsq;
        do {
			v1 = 2.0*Math.random()-1.0;
			v2 = 2.0*Math.random()-1.0;
			rsq = v1*v1+v2*v2;
		} while ( (rsq<=0) || (rsq>=1) );
		if (v1>0) return Math.sqrt(sigma*sigma + Math.sqrt(2.0*sigma*sigma*sigma/Math.PI/rsq));
		else return -Math.sqrt(sigma*sigma + Math.sqrt(2.0*sigma*sigma*sigma/Math.PI/rsq));
    }
	
}
