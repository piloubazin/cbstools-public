package de.mpg.cbs.utilities;


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

public class Matrix3D {
	
	/** mathematical constants */
	public static final float ZERO = 1e-30f;
	/** mathematical constants */
	public static final float INF = 1e+30f;
	/** mathematical constants */
	public static final float PI = 3.1416f;
	
	/**
	 *  Close-form 3D matrix determinant
	 */
	public static final float determinant(float m[][]) {
		float det;
		det = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2]
			- m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];

		return det;
	}
	
	/**
	 *  Close-form 3D matrix determinant
	 */
	public static final double determinant(double m[][]) {
		double det;
		det = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2]
			- m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];

		return det;
	}
	
	/**
	 *  Close-form 3D matrix inversion
	 */
	public static final void invert(float m[][]) {
		float[] d = new float[9];
		float det;
		det = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2]
			- m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];
		det = 1.0f/det;

		d[0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1])*det;
		d[1] = (m[2][1]*m[0][2]-m[2][2]*m[0][1])*det;
		d[2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1])*det;
		d[3] = (m[1][2]*m[2][0]-m[1][0]*m[2][2])*det;
		d[4] = (m[0][0]*m[2][2]-m[0][2]*m[2][0])*det;
		d[5] = (m[1][0]*m[0][2]-m[1][2]*m[0][0])*det;
		d[6] = (m[1][0]*m[2][1]-m[1][1]*m[2][0])*det;
		d[7] = (m[0][1]*m[2][0]-m[0][0]*m[2][1])*det;
		d[8] = (m[0][0]*m[1][1]-m[0][1]*m[1][0])*det;

		m[0][0] = d[0];
		m[1][0] = d[1];
		m[2][0] = d[2];
		m[0][1] = d[3];
		m[1][1] = d[4];
		m[2][1] = d[5];
		m[0][2] = d[6];
		m[1][2] = d[7];
		m[2][2] = d[8];

		return;
	}
	
	/**
	 * This method compute the inverse of a 3D affine transform matrix 
	 * @param m : the given matrix which should be a square array of floats 
	 * @param inverse : is the result of inversion
	 * @return : determinent of the matrix
	 */
	
	public static final float[][] invertAffineTransform(float m[][]){
		
		float det = 0.0f;
		//float[][] inverse = new float[nColumn][nColumn];
		
		float[][] inverse = new float[3][4];
		
		det = m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2]
			- m[2][0]*m[1][1]*m[0][2] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[1][0]*m[0][1];
		if (det != 0.0f) {
			inverse[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;
			inverse[0][1]= (m[2][1]*m[0][2]-m[2][2]*m[0][1])/det;
			inverse[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;
			inverse[1][0] = (m[1][2]*m[2][0]-m[1][0]*m[2][2])/det;
			inverse[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;
			inverse[1][2] = (m[1][0]*m[0][2]-m[1][2]*m[0][0])/det;
			inverse[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;
			inverse[2][1]= (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det;
			inverse[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;
			
			inverse[0][3] = -inverse[0][0]*m[0][3]-inverse[0][1]*m[1][3]-inverse[0][2]*m[2][3];
			inverse[1][3] = -inverse[1][0]*m[0][3]-inverse[1][1]*m[1][3]-inverse[1][2]*m[2][3];
			inverse[2][3] = -inverse[2][0]*m[0][3]-inverse[2][1]*m[1][3]-inverse[2][2]*m[2][3];
		} else {
			System.out.println("Matrix is singular");
		}
		return inverse;
	}
	
	public static final float[] eigenvalues(float m[][]) {
		float[] values = new float[3];
		
		// polynomial coefficients X^3 + a2 X^2 + a1 X + a0 = 0
		
		double a2 = - (m[0][0]+m[1][1]+m[2][2]);
		
		double a1 =   (m[0][0]*m[1][1]-m[0][1]*m[1][0])
					+(m[1][1]*m[2][2]-m[1][2]*m[2][1])
					+(m[2][2]*m[0][0]-m[2][0]*m[0][2]);
		 
		double a0 = -( m[0][0]*m[1][1]*m[2][2]
					  +m[0][1]*m[1][2]*m[2][0]
					  +m[1][0]*m[0][2]*m[2][1]
					  -m[0][1]*m[1][0]*m[2][2]
					  -m[1][2]*m[2][1]*m[0][0]
					  -m[2][0]*m[0][2]*m[1][1] );
	
		double Q = ( 3.0*a1-a2*a2 )/9.0;
		double R = ( 9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2 )/54.0;
		
		// discriminant
		double D = Q*Q*Q + R*R;
		
		if ( D > 0 ) { 
			// positive discriminant: imaginary solutions
			//System.out.print("positive discriminant");
			values[0] = 0.0f;
			values[1] = 0.0f;
			values[2] = 0.0f;
		} else if ( D < 0 ) {
			// usual case
			double theta = Math.acos( R / Math.sqrt( -Q*Q*Q ) );
			double sQ = 2.0*Math.sqrt(-Q);
			
			values[0] = (float)(sQ*Math.cos(theta/3.0) - a2/3.0 );
			values[1] = (float)(sQ*Math.cos((theta+2.0*Math.PI)/3.0) - a2/3.0 );
			values[2] = (float)(sQ*Math.cos((theta+4.0*Math.PI)/3.0) - a2/3.0 );
			
			// re-ordering by decreasing size : |values[0]| >= |values[1]| >= |values[2]|
			if (Numerics.abs(values[0])<Numerics.abs(values[1])) { float tmp = values[0]; values[0] = values[1]; values[1] = tmp; }
			if (Numerics.abs(values[0])<Numerics.abs(values[2])) { float tmp = values[0]; values[0] = values[2]; values[2] = tmp; }
			if (Numerics.abs(values[1])<Numerics.abs(values[2])) { float tmp = values[1]; values[1] = values[2]; values[2] = tmp; }					
		} else {
			// multiple roots
			double S = Math.cbrt(R);
			
			values[0] = (float)(-a2/3.0 + 2.0*S);
			values[1] = (float)(-a2/3.0 - S);
			values[2] = values[1];
			
			// re-ordering by decreasing size : |values[0]| >= |values[1]| >= |values[2]|
			if (Numerics.abs(values[0])<Numerics.abs(values[1])) { float tmp = values[0]; values[0] = values[2]; values[2] = tmp; }
		}
		return values;
	}
	
	/** 
	 *	Simple geometric technique to get eigenvectors of 3D matrices.
	 *	<p>
	 *	This may not work well all the time: sometimes you have a 2 or 3D sub-space with one eigenvalue.
	 *	We assume here that the highest eigenvalues give the most reliable directions, 
	 *  and lesser eigenvectors are projected on orthogonal subspaces (seems wrong: to check).
	 *	We also assume ordered, positive eigenvalues. 
	 */
	public static final float[][] eigenvectors(float m[][], float values[]) {
		float[][] vector = new float[3][3];
		
		double[] ArIx = new double[3];
		double[] ArIy = new double[3];
		double[] ArIz = new double[3];
		
		if (values[0]==0) {
			vector[0][0] = 1.0f;
			vector[1][0] = 0.0f;			
			vector[2][0] = 0.0f;
			vector[0][1] = 0.0f;
			vector[1][1] = 1.0f;			
			vector[2][1] = 0.0f;
			vector[0][2] = 0.0f;
			vector[1][2] = 0.0f;			
			vector[2][2] = 1.0f;
			return vector;
		}			
		
		for (int i=0;i<3;i++) {
		
			// first eigenvalue
			ArIx[0] = m[0][0]-values[i];
			ArIx[1] = m[1][0];
			ArIx[2] = m[2][0];
			double normx2 = ArIx[0]*ArIx[0]+ArIx[1]*ArIx[1]+ArIx[2]*ArIx[2];
			
			ArIy[0] = m[0][1];
			ArIy[1] = m[1][1]-values[i];
			ArIy[2] = m[2][1];
			double normy2 = ArIy[0]*ArIy[0]+ArIy[1]*ArIy[1]+ArIy[2]*ArIy[2];
			
			ArIz[0] = m[0][2];
			ArIz[1] = m[1][2];
			ArIz[2] = m[2][2]-values[i];
			double normz2 = ArIz[0]*ArIz[0]+ArIz[1]*ArIz[1]+ArIz[2]*ArIz[2];
			
			if ( (normx2<normy2) && (normx2<normz2) ) {
				vector[0][i] = (float)(ArIy[1]*ArIz[2]-ArIy[2]*ArIz[1]);
				vector[1][i] = (float)(ArIy[2]*ArIz[0]-ArIy[0]*ArIz[2]);			
				vector[2][i] = (float)(ArIy[0]*ArIz[1]-ArIy[1]*ArIz[0]);
			} else if ( (normy2<normz2) && (normy2<normx2) ) {
				vector[0][i] = (float)(ArIz[1]*ArIx[2]-ArIz[2]*ArIx[1]);
				vector[1][i] = (float)(ArIz[2]*ArIx[0]-ArIz[0]*ArIx[2]);			
				vector[2][i] = (float)(ArIz[0]*ArIx[1]-ArIz[1]*ArIx[0]);
			} else {
				vector[0][i] = (float)(ArIx[1]*ArIy[2]-ArIx[2]*ArIy[1]);
				vector[1][i] = (float)(ArIx[2]*ArIy[0]-ArIx[0]*ArIy[2]);			
				vector[2][i] = (float)(ArIx[0]*ArIy[1]-ArIx[1]*ArIy[0]);
			}
			// orthogonalize? makes sens for symmetric matrices, but not in general case
			/*
			if (i==1) {
				// v0 is already normalized
				float v0v1 = vector[0][0]*vector[0][1]+vector[1][0]*vector[1][1]+vector[2][0]*vector[2][1];
				vector[0][1] = vector[0][1] - v0v1*vector[0][0];
				vector[1][1] = vector[1][1] - v0v1*vector[1][0];	
				vector[2][1] = vector[2][1] - v0v1*vector[2][0];
			} else if (i==2) {
				// v0, v1 normalized
				// (note: it would be tempting to take the cross product, but it may not be well defined)
				float v0v2 = vector[0][0]*vector[0][2]+vector[1][0]*vector[1][2]+vector[2][0]*vector[2][2];
				float v1v2 = vector[0][1]*vector[0][2]+vector[1][1]*vector[1][2]+vector[2][1]*vector[2][2];
				vector[0][2] = vector[0][2] - v0v2*vector[0][0] - v1v2*vector[0][1];
				vector[1][2] = vector[1][2] - v0v2*vector[1][0] - v1v2*vector[1][1];	
				vector[2][2] = vector[2][2] - v0v2*vector[2][0] - v1v2*vector[2][1];
			}
			*/
			// normalize
			double norm = Math.sqrt(vector[0][i]*vector[0][i]+vector[1][i]*vector[1][i]+vector[2][i]*vector[2][i]);
			vector[0][i] = (float)(vector[0][i]/norm);
			vector[1][i] = (float)(vector[1][i]/norm);	
			vector[2][i] = (float)(vector[2][i]/norm);			
		}

		/* debug: comment out for faster code 
		// check for orthogonality, same value ?
		float prodxy = vector[0][0]*vector[0][1] + vector[1][0]*vector[1][1] + vector[2][0]*vector[2][1];
		float prodyz = vector[0][1]*vector[0][2] + vector[1][1]*vector[1][2] + vector[2][1]*vector[2][2];
		float prodzx = vector[0][2]*vector[0][0] + vector[1][2]*vector[1][0] + vector[2][2]*vector[2][0];
			
		if ( (prodxy > 0.1) || (prodxy < -0.1) 
			|| (prodyz > 0.1) || (prodyz < -0.1) 
			|| (prodzx > 0.1) || (prodzx < -0.1) ) {
			System.out.print("pb:"+prodxy+"|"+prodyz+"|"+prodzx+"|");
			System.out.print("l:"+values[0]+"|"+values[1]+"|"+values[2]+"|");
		}
		*/
		return vector;
	}
	
	/** 
	 *	Close-form 3D matrix eigenvalues.
	 *	<p>
	 *	This may not work well all the time with arbitrary matrices.
	 */
	public static final double[] eigenvalues(double m[][]) {
		double[] values = new double[3];
		
		// polynomial coefficients X^3 + a2 X^2 + a1 X + a0 = 0
		
		double a2 = - (m[0][0]+m[1][1]+m[2][2]);
		
		double a1 =   (m[0][0]*m[1][1]-m[0][1]*m[1][0])
					+(m[1][1]*m[2][2]-m[1][2]*m[2][1])
					+(m[2][2]*m[0][0]-m[2][0]*m[0][2]);
		 
		double a0 = -( m[0][0]*m[1][1]*m[2][2]
					  +m[0][1]*m[1][2]*m[2][0]
					  +m[1][0]*m[0][2]*m[2][1]
					  -m[0][1]*m[1][0]*m[2][2]
					  -m[1][2]*m[2][1]*m[0][0]
					  -m[2][0]*m[0][2]*m[1][1] );
	
		double Q = ( 3.0*a1-a2*a2 )/9.0;
		double R = ( 9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2 )/54.0;
		
		// discriminant
		double D = Q*Q*Q + R*R;
		
		if ( D > 0 ) { 
			// positive discriminant: imaginary solutions
			//System.out.print("positive discriminant");
			values[0] = 0.0f;
			values[1] = 0.0f;
			values[2] = 0.0f;
		} else if ( D < 0 ) {
			// usual case
			double theta = Math.acos( R / Math.sqrt( -Q*Q*Q ) );
			double sQ = 2.0*Math.sqrt(-Q);
			
			values[0] = (sQ*Math.cos(theta/3.0) - a2/3.0 );
			values[1] = (sQ*Math.cos((theta+2.0*Math.PI)/3.0) - a2/3.0 );
			values[2] = (sQ*Math.cos((theta+4.0*Math.PI)/3.0) - a2/3.0 );
			
			// re-ordering by decreasing size : |values[0]| >= |values[1]| >= |values[2]|
			if (Numerics.abs(values[0])<Numerics.abs(values[1])) { double tmp = values[0]; values[0] = values[1]; values[1] = tmp; }
			if (Numerics.abs(values[0])<Numerics.abs(values[2])) { double tmp = values[0]; values[0] = values[2]; values[2] = tmp; }
			if (Numerics.abs(values[1])<Numerics.abs(values[2])) { double tmp = values[1]; values[1] = values[2]; values[2] = tmp; }					
		} else {
			// multiple roots
			double S = Math.cbrt(R);
			
			values[0] = (-a2/3.0 + 2.0*S);
			values[1] = (-a2/3.0 - S);
			values[2] = values[1];
			
			// re-ordering by decreasing size : |values[0]| >= |values[1]| >= |values[2]|
			if (Numerics.abs(values[0])<Numerics.abs(values[1])) { double tmp = values[0]; values[0] = values[2]; values[2] = tmp; }
		}
		return values;
	}
	
	/** 
	 *	Simple geometric technique to get eigenvectors of 3D matrices.
	 *	<p>
	 *	This may not work well all the time: sometimes you have a 2 or 3D sub-space with one eigenvalue.
	 *	We assume here that the highest eigenvalues give the most reliable directions, and lesser eigenvectors
	 *	are projected on orthogonal subspaces.
	 *	We also assume ordered, positive eigenvalues (needed for orthogonal eigenvectors). 
	 */
	public static final double[][] eigenvectors(double m[][], double values[]) {
		double[][] vector = new double[3][3];
		
		double[] ArIx = new double[3];
		double[] ArIy = new double[3];
		double[] ArIz = new double[3];
		
		if (values[0]==0) {
			vector[0][0] = 1.0f;
			vector[1][0] = 0.0f;			
			vector[2][0] = 0.0f;
			vector[0][1] = 0.0f;
			vector[1][1] = 1.0f;			
			vector[2][1] = 0.0f;
			vector[0][2] = 0.0f;
			vector[1][2] = 0.0f;			
			vector[2][2] = 1.0f;
			return vector;
		}			
		
		for (int i=0;i<3;i++) {
		
			// first eigenvalue
			ArIx[0] = m[0][0]-values[i];
			ArIx[1] = m[1][0];
			ArIx[2] = m[2][0];
			double normx2 = ArIx[0]*ArIx[0]+ArIx[1]*ArIx[1]+ArIx[2]*ArIx[2];
			
			ArIy[0] = m[0][1];
			ArIy[1] = m[1][1]-values[i];
			ArIy[2] = m[2][1];
			double normy2 = ArIy[0]*ArIy[0]+ArIy[1]*ArIy[1]+ArIy[2]*ArIy[2];
			
			ArIz[0] = m[0][2];
			ArIz[1] = m[1][2];
			ArIz[2] = m[2][2]-values[i];
			double normz2 = ArIz[0]*ArIz[0]+ArIz[1]*ArIz[1]+ArIz[2]*ArIz[2];
			
			if ( (normx2<normy2) && (normx2<normz2) ) {
				vector[0][i] = (ArIy[1]*ArIz[2]-ArIy[2]*ArIz[1]);
				vector[1][i] = (ArIy[2]*ArIz[0]-ArIy[0]*ArIz[2]);			
				vector[2][i] = (ArIy[0]*ArIz[1]-ArIy[1]*ArIz[0]);
			} else if ( (normy2<normz2) && (normy2<normx2) ) {
				vector[0][i] = (ArIz[1]*ArIx[2]-ArIz[2]*ArIx[1]);
				vector[1][i] = (ArIz[2]*ArIx[0]-ArIz[0]*ArIx[2]);			
				vector[2][i] = (ArIz[0]*ArIx[1]-ArIz[1]*ArIx[0]);
			} else {
				vector[0][i] = (ArIx[1]*ArIy[2]-ArIx[2]*ArIy[1]);
				vector[1][i] = (ArIx[2]*ArIy[0]-ArIx[0]*ArIy[2]);			
				vector[2][i] = (ArIx[0]*ArIy[1]-ArIx[1]*ArIy[0]);
			}
			// orthogonalize? not correct for the general case
			/*
			if (i==1) {
				// v0 is already normalized
				double v0v1 = vector[0][0]*vector[0][1]+vector[1][0]*vector[1][1]+vector[2][0]*vector[2][1];
				vector[0][1] = vector[0][1] - v0v1*vector[0][0];
				vector[1][1] = vector[1][1] - v0v1*vector[1][0];	
				vector[2][1] = vector[2][1] - v0v1*vector[2][0];
			} else if (i==2) {
				// v0, v1 normalized
				// (note: it would be tempting to take the cross product, but it may not be well defined)
				double v0v2 = vector[0][0]*vector[0][2]+vector[1][0]*vector[1][2]+vector[2][0]*vector[2][2];
				double v1v2 = vector[0][1]*vector[0][2]+vector[1][1]*vector[1][2]+vector[2][1]*vector[2][2];
				vector[0][2] = vector[0][2] - v0v2*vector[0][0] - v1v2*vector[0][1];
				vector[1][2] = vector[1][2] - v0v2*vector[1][0] - v1v2*vector[1][1];	
				vector[2][2] = vector[2][2] - v0v2*vector[2][0] - v1v2*vector[2][1];
			}
			*/
			// normalize
			double norm = Math.sqrt(vector[0][i]*vector[0][i]+vector[1][i]*vector[1][i]+vector[2][i]*vector[2][i]);
			vector[0][i] = (vector[0][i]/norm);
			vector[1][i] = (vector[1][i]/norm);	
			vector[2][i] = (vector[2][i]/norm);			
		}

		return vector;
	}
	
}
