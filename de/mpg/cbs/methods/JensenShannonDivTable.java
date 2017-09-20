package de.mpg.cbs.methods;

import de.mpg.cbs.utilities.Numerics;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.analysis.function.Gaussian;

/**
 *
 *	@version    February 2016
 *	@author     Pierre-Louis Bazin
 */

/**
 * convenience class to compute Jensen-Shannon divergence of Gaussian mixtures quickly.
 * the values are pre-computed and then accessed in the computations.
 * non-covered cases use the regular estimation function instead.
 */
public class JensenShannonDivTable {
    private double range;			// max relative distance between means [0, 3]
    private double resolution;		// discretization step of the parameters (affects matrix size)
    private double sampling;			// discretization step of the JSD estimator
    private int maxdiff;
    private int maxratio;

    private double[][] jsd;

	private final double	ZERO = 1e-30;
    private final double	INF = 1e+30;
	
    // setting up the table: max range should be 3, with a table size of 1/resolution x 3/resolution
    public JensenShannonDivTable(double range_, double resolution_, double sampling_){
   	   range = range_;
   	   resolution = resolution_;
   	   sampling = sampling_;
   	   
   	   maxdiff = Numerics.ceil(range/resolution);
   	   maxratio = Numerics.ceil(1.0/resolution);
   	   jsd = new double[maxdiff][maxratio];
   	   
   	   for (int d=0;d<maxdiff;d++) for (int r=1;r<maxratio;r++) {
   	   	   jsd[d][r] = computeJSD(d*resolution, r*resolution);
   	   }
   }
   
   public final double lookup(double mu1, double sq1, double mu2, double sq2) {	
		// normalize
		int diff, ratio;
		if (sq1>sq2) {
			diff = Numerics.round(Numerics.abs(mu1-mu2)/sq1/resolution);
			ratio = Numerics.round(sq2/sq1/resolution);
		} else {
			diff = Numerics.round(Numerics.abs(mu1-mu2)/sq2/resolution);
			ratio = Numerics.round(sq1/sq2/resolution);
		}
     	if (diff<maxdiff && ratio<maxratio) {
     		return jsd[diff][ratio];
     	} else {
     		// simpler: use worst case scenario
     		System.out.print(";");
     		//return computeJSD(diff*resolution, ratio*resolution);
     		return jsd[maxdiff-1][maxratio-1];
     	}
   }
	
   private final double computeJSD(double diff, double ratio) {
		// here we estimate the JS divergence explicitly
		double xmin = Numerics.min(-3.0, diff-3.0*ratio);
		double xmax = Numerics.max( 3.0, diff+3.0*ratio);
		Gaussian p1 = new Gaussian(0.0, 1.0);
		Gaussian p2 = new Gaussian(diff, ratio);
		
		double jsdiv = 0.0;
		double step = sampling*(xmax-xmin);
		for (double x=xmin;x<=xmax;x+=step) {
			double p1x = p1.value(x);
			double p2x = p2.value(x);
			jsdiv += p1x*FastMath.log(2.0*p1x/(p1x+p2x)) + p2x*FastMath.log(2.0*p2x/(p1x+p2x));
		}
		return jsdiv*step/2.0;
	}

}
