package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;
import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.FastMath;
import Jama.Matrix;

/**
 *
 *  This class deals with various simple polynomial computations
 *	
 *	@version    April 2018
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Polynomials {
	
	// no data: used as a library of functions
	
    // simple tensor functions
	
	// constants
	private static final double SQ2 = Math.sqrt(2.0);
    private static final double ZERO = 1e-30;
    
    public static final double[] normalizedChebyshevFit(double[] data, double[] location, int degree) {
        int nsample = data.length;
		double[][] chebcoeff = new double[nsample][degree+1];
		double[][] mdata = new double[nsample][1];
		double min = location[0], max = location[0];
		for (int n=1;n<nsample;n++) {
		    if (location[n]<min) min = location[n];
		    if (location[n]>max) max = location[n];
		}
		for (int n=0;n<nsample;n++) {
		    computeChebyshevCoefficients(chebcoeff[n], (location[n]-min)/(max-min), degree);
		    mdata[n][0] = data[0];
		}
		// invert the linear model
		Matrix mtx = new Matrix(chebcoeff);
		Matrix smp = new Matrix(mdata);
		Matrix val = mtx.solve(smp);
		
		double[] coeffs = new double[degree+1];
		for (int d=0;d<=degree;d++) {
			coeffs[d] = val.get(d,0);
		}
		return coeffs;
	}

    public static final double[] chebyshevFit(double[] data, double[] location, int degree) {
        int nsample = data.length;
		double[][] chebcoeff = new double[nsample][degree+1];
		double[][] mdata = new double[nsample][1];
		for (int n=0;n<nsample;n++) {
		    computeChebyshevCoefficients(chebcoeff[n], location[n], degree);
		    mdata[n][0] = data[0];
		}
		// invert the linear model
		Matrix mtx = new Matrix(chebcoeff);
		Matrix smp = new Matrix(mdata);
		Matrix val = mtx.solve(smp);
		
		double[] coeffs = new double[degree+1];
		for (int d=0;d<=degree;d++) {
			coeffs[d] = val.get(d,0);
		}
		return coeffs;
	}

	public static final void computeChebyshevCoefficients(double[] pol, double val, int deg) {
		for (int n=0;n<=deg;n++) {
			pol[n] = 0.0;
		}
		double r = (2.0*(double)val-1.0);
		
		if (deg>=0) pol[0] = 1.0;
		if (deg>=1) pol[1] = r;
		if (deg>=2) pol[2] = 2.0*r*r - 1.0;
		if (deg>=3) pol[3] = 4.0*r*r*r - 3.0*r;
		if (deg>=4) pol[4] = 8.0*r*r*r*r - 8.0*r*r + 1.0;
		if (deg>=5) pol[5] = 16.0*r*r*r*r*r - 20.0*r*r*r + 5.0*r;
		if (deg>=6) pol[6] = 32.0*r*r*r*r*r*r - 48.0*r*r*r*r + 18.0*r*r - 1.0;	
		if (deg==7) pol[7] = 64.0*r*r*r*r*r*r*r - 112.0*r*r*r*r*r + 56.0*r*r*r - 7.0*r;
		return;
	}
	 	 
}
