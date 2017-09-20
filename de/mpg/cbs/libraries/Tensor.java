package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This class computes various basic DTI tensor functions, e.g.
 *  fractioanl anisotropy, decomposition, log-Euclidean operators, etc.
 *	
 *	@version    February 2008
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class Tensor {
	
	// no data: used as a library of functions
	
    // simple tensor functions
	
	// constants
	private static final double SQ2 = Math.sqrt(2.0);
    private static final double ZERO = 1e-30;
    
    // numbering for the eigenvalue, eigenvector representation
	public static final int L1 = 0;
	public static final int L2 = 1;
	public static final int L3 = 2;
	public static final int V1 = 3;
    public static final int V1X = 3;
    public static final int V1Y = 4;
    public static final int V1Z = 5;
    public static final int V2 = 6;
    public static final int V2X = 6;
    public static final int V2Y = 7;
    public static final int V2Z = 8;
	public static final int V3 = 9;
    public static final int V3X = 9;
    public static final int V3Y = 10;
    public static final int V3Z = 11;

	/**
	 *	Tensor mapping from vector to matrix form
	 */
	public static final double[][] matrixForm(float [] tensor) {
		
		double[][] mat = new double[3][3];
		
		mat[0][0] = tensor[0];
		mat[1][1] = tensor[1];
		mat[2][2] = tensor[2];
		
		mat[0][1] = tensor[3];
		mat[1][2] = tensor[4];
		mat[2][0] = tensor[5];
		
		mat[1][0] = tensor[3];
		mat[2][1] = tensor[4];
		mat[0][2] = tensor[5];
		
		return mat;
	}
	/**
	 *	Tensor mapping from matrix to vector form
	 */
	public static final float[] vectorForm(double[][] mat) {
		
		float[] vect = new float[6];
		
		vect[0] = (float)mat[0][0];
		vect[1] = (float)mat[1][1];
		vect[2] = (float)mat[2][2];
		
		vect[3] = 0.5f*(float)(mat[0][1]+mat[1][0]);
		vect[4] = 0.5f*(float)(mat[1][2]+mat[2][1]);
		vect[5] = 0.5f*(float)(mat[2][0]+mat[0][2]);
		
		return vect;
	}
	/**
	 *	Tensor mapping from diagonal vector (assumed convention) to triangular vector form (DTIStudio convention)
	 */
	public static final float[] triangularForm(float [] tensor) {
		
		float[] triangular = new float[6];
		
		triangular[0] = tensor[0];
		triangular[1] = tensor[3];
		triangular[2] = tensor[5];
		triangular[3] = tensor[1];
		triangular[4] = tensor[4];
		triangular[5] = tensor[2];
		
		return triangular;
	}
	/**
	 *	Tensor mapping from triangular vector form (DTIStudio convention) to diagonal vector (assumed convention)
	 */
	public static final float[] diagonalForm(float [] tensor) {
		
		float[] diagonal = new float[6];
		
		diagonal[0] = tensor[0];
		diagonal[1] = tensor[3];
		diagonal[2] = tensor[5];
		diagonal[3] = tensor[1];
		diagonal[4] = tensor[4];
		diagonal[5] = tensor[2];
		
		return diagonal;
	}

	/**
	 *	Tensor decomposition into {Eigenvalues,Eigenvectors}
	 *	assuming the tensor is given as a 3x3 matrix
	 */
	public static final double[] eigenDecomposition(double[][] tensor) {
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double[][] vect = Matrix3D.eigenvectors(tensor,lambda);

		// build new vector
		double[] eigen = new double[12];
		for (int n=0;n<3;n++) {
			eigen[n] = lambda[n];
			for (int m=0;m<3;m++) eigen[3+3*n+m] = vect[m][n];
		}
		return eigen;
	}
	
	/**
	 *	Tensor recomposition from {Eigenvalues,Eigenvectors}
	 */
	public static final double[][] eigenTensor(double [] eigen) {
		double[][] tensor = new double[3][3];
		
		for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
			tensor[n][m] = 0.0f;
			for (int p=0;p<3;p++) tensor[n][m] += eigen[3+3*n+p]*eigen[p]*eigen[3+3*m+p];
		}
		
		return tensor;
	}
	
	/**
	 *	inverse Tensor recomposition from {Eigenvalues,Eigenvectors}
	 */
	public static final double[][] inverseEigenTensor(double [] eigen) {
		double[][] tensor = new double[3][3];
		
		for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
			tensor[n][m] = 0.0f;
			for (int p=0;p<3;p++) tensor[n][m] += eigen[3+3*p+n]/eigen[p]*eigen[3+3*p+m];
		}
		
		return tensor;
	}
	
	/**
	 *	Tensor logarithm
	 *	assuming the tensor is given as a 3x3 matrix
	 */
	public static final double[][] log(double [][] tensor) {
		double[][] logt = new double[3][3];
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double[][] vect = Matrix3D.eigenvectors(tensor,lambda);

		for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
			logt[n][m] = 0.0f;
			for (int p=0;p<3;p++) logt[n][m] += (vect[p][n]*Math.log(lambda[p])*vect[p][m]);
		}
		return logt;
	}
	
	/**
	 *	Tensor exponential
	 *	assuming the tensor is given as a 3x3 matrix
	 */
	public static final double[][] exp(double [][] tensor) {
		double[][] expt = new double[3][3];
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double[][] vect = Matrix3D.eigenvectors(tensor,lambda);

		for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
			expt[n][m] = 0.0f;
			for (int p=0;p<3;p++) expt[n][m] += (vect[p][n]*Math.exp(lambda[p])*vect[p][m]);
		}
		return expt;
	}
	
	/**
	 *	Tensor transform into Log-Euclidean vector
	 *	assuming the tensor is given as a 3x3 matrix
	 */
	public static  final float[] logEuclideanVector(double[][] tensor) {
		float[] lev = new float[6];
		
		double[][] logt = log(tensor);
		
		lev[0] = (float)logt[0][0];
		lev[1] = (float)logt[1][1];
		lev[2] = (float)logt[2][2];
		lev[3] = (float)(SQ2*0.5*(logt[0][1]+logt[1][0]));
		lev[4] = (float)(SQ2*0.5*(logt[1][2]+logt[2][1]));
		lev[5] = (float)(SQ2*0.5*(logt[2][0]+logt[0][2]));

		return lev;
	}
	
	/**
	 *	transform Log-Euclidean vector back to tensor
	 */
	public static  final double[][] logEuclideanTensor(float[] lev) {
		double[][] logt = new double[3][3];
		
		logt[0][0] = lev[0];
		logt[1][1] = lev[1];
		logt[2][2] = lev[2];
		
		logt[0][1] = lev[3]/SQ2;
		logt[1][2] = lev[4]/SQ2;
		logt[2][0] = lev[5]/SQ2;
		
		logt[1][0] = lev[3]/SQ2;
		logt[2][1] = lev[4]/SQ2;
		logt[0][2] = lev[5]/SQ2;
		
		return exp(logt);
	}
	
	/**
	 *	Tensor transform into Log-Euclidean vector
	 *	assuming the tensor is given as its eigendecomposition
	 */
	public static  final float[] logEuclideanVector(double[] eigen) {
		float[] lev = new float[6];
		
		double[][] logt = new double[3][3];
		for (int n=0;n<3;n++) for (int m=0;m<=n;m++) {
			logt[n][m] = 0.0f;
			for (int p=0;p<3;p++) logt[n][m] += eigen[3+3*n+p]*Math.log(eigen[p])*eigen[3+3*m+p];
		}
		
		lev[0] = (float)logt[0][0];
		lev[1] = (float)logt[1][1];
		lev[2] = (float)logt[2][2];
		lev[3] = (float)(SQ2*0.5*(logt[1][0]+logt[0][1]));
		lev[4] = (float)(SQ2*0.5*(logt[2][1]+logt[1][2]));
		lev[5] = (float)(SQ2*0.5*(logt[2][0]+logt[0][2]));

		return lev;
	}
	
	/**
	 *	transform Log-Euclidean vector back to tensor eigendecomposition
	 */
	public static  final double[] logEuclideanEigenTensor(double[] lev) {
		double[][] logt = new double[3][3];
		
		logt[0][0] = lev[0];
		logt[1][1] = lev[1];
		logt[2][2] = lev[2];
		
		logt[0][1] = lev[3]/SQ2;
		logt[1][2] = lev[4]/SQ2;
		logt[2][0] = lev[5]/SQ2;
		
		logt[1][0] = lev[3]/SQ2;
		logt[2][1] = lev[4]/SQ2;
		logt[0][2] = lev[5]/SQ2;
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(logt);
		double[][] vect = Matrix3D.eigenvectors(logt,lambda);

		// build new vector
		double[] eigen = new double[12];
		for (int n=0;n<3;n++) {
			eigen[n] = Math.exp(lambda[n]);
			for (int m=0;m<3;m++) eigen[3+3*n+m] = vect[m][n];
		}
		return eigen;
	}
	
	/**
	 *	extract the principla direction from a tensor
	 *	assuming the tensor is given as a 3x3 matrix
	 */
	public static final float[] principalDirection(double [][] tensor) {
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double[][] vect = Matrix3D.eigenvectors(tensor,lambda);

		// build new vector
		float[] dir1 = new float[3];
		for (int n=0;n<3;n++) {
			dir1[n] = (float)vect[n][0];
		}
		return dir1;
	}
	
	/**
	 *	extract the principla direction from a tensor
	 *	assuming the tensor is given by its eigendecomposition
	 */
	public static final float[] principalDirection(double [] eigen) {
		// return main direction vector
		float[] dir1 = new float[3];
		for (int n=0;n<3;n++) {
			dir1[n] = (float)eigen[3+n];
		}
		return dir1;
	}
	
	/**
	 *	extract the principla direction from a tensor
	 *	assuming the tensor is given by its eigendecomposition
	 */
	public static final float[] secondDirection(double [] eigen) {
		// return main direction vector
		float[] dir2 = new float[3];
		for (int n=0;n<3;n++) {
			dir2[n] = (float)eigen[6+n];
		}
		return dir2;
	}
	
	/**
	 *	extract the principla direction from a tensor
	 *	assuming the tensor is given by its eigendecomposition
	 */
	public static final float[] thirdDirection(double [] eigen) {
		// return main direction vector
		float[] dir3 = new float[3];
		for (int n=0;n<3;n++) {
			dir3[n] = (float)eigen[9+n];
		}
		return dir3;
	}
	
	/**
	 *	extract direction colormap values from a tensor
	 *	assuming the tensor is given as a 3x3 matrix
	 *	The color is normalized in [-127,128]
	 */
	public static final byte[] directionColor(double [][] tensor) {
		
		// decomposition
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double[][] vect = Matrix3D.eigenvectors(tensor,lambda);

		// build new vector
		byte[] color = new byte[3];
		for (int n=0;n<3;n++) {
			color[n] = (byte)(255*Math.abs(vect[n][0])-127);
		}
		return color;
	}
	
	/**
	 *	changes directions so that the directions are always 
	 *	in positive half-space
	 */
	public static final void alignDirections(double [] eigen) {
		for (int n=0;n<3;n++) {
			int best = Numerics.argmax(eigen[3+3*n]*eigen[3+3*n],eigen[4+3*n]*eigen[4+3*n],eigen[5+3*n]*eigen[5+3*n]);
		
			if (eigen[3+best+3*n]<0) {
				for (int d=0;d<3;d++) {
					eigen[3+d+3*n] = - eigen[3+d+3*n];
				}
			}
		}
	}
	
	/**
	 *	measure fractional anisotropy from tensor
	 */
	public static  final float fractionalAnisotropy(double[][] tensor) {
		double[] lambda = Matrix3D.eigenvalues(tensor);
		double lambdaM = (lambda[0]+lambda[1]+lambda[2])/3.0f;
		double lambdaP = (lambda[0]*lambda[0]+lambda[1]*lambda[1]+lambda[2]*lambda[2]);
		float fa;
		if (lambdaP>0) {
			fa = (float)Math.sqrt(3.0f*( (lambda[0]-lambdaM)*(lambda[0]-lambdaM)
									     +(lambda[1]-lambdaM)*(lambda[1]-lambdaM)
									     +(lambda[2]-lambdaM)*(lambda[2]-lambdaM) )
										  /(2.0f*lambdaP));
		} else {
			fa = 0;
		}
		return fa;
	}		
	
	/**
	 *	measure fractional anisotropy from eigendecomposition
	 */
	public static  final float fractionalAnisotropy(double[] eigen) {
		double lambdaM = (eigen[0]+eigen[1]+eigen[2])/3.0f;
		double lambdaP = (eigen[0]*eigen[0]+eigen[1]*eigen[1]+eigen[2]*eigen[2]);
		float fa;
		if (lambdaP>0) {
			fa = (float)Math.sqrt(3.0f*( (eigen[0]-lambdaM)*(eigen[0]-lambdaM)
									     +(eigen[1]-lambdaM)*(eigen[1]-lambdaM)
									     +(eigen[2]-lambdaM)*(eigen[2]-lambdaM) )
										  /(2.0f*lambdaP) );
		} else {
			fa = 0;
		}
		return fa;
	}		
	
	/**
	 *	measure mean diffusivity from tensor
	 */
	public static final float meanDiffusivity(double[][] tensor) {
		return (float)(tensor[0][0]+tensor[1][1]+tensor[2][2])/3.0f;
	}		
	
	/**
	 *	measure mean diffusivity from eigendecomposition
	 */
	public static final float meanDiffusivity(double[] eigen) {
		return (float)(eigen[0]+eigen[1]+eigen[2])/3.0f;
	}		
	
	/**
	 *	measure mean diffusivity from eigendecomposition
	 */
	public static final float parallelDiffusivity(double[] eigen) {
		return (float)(eigen[0]);
	}		
	
	/**
	 *	measure mean diffusivity from eigendecomposition
	 */
	public static final float perpendicularDiffusivity(double[] eigen) {
		return (float)(eigen[1]+eigen[2])/2.0f;
	}		
	
	/**
	 *	measure fractional anisotropy from eigendecomposition
	 */
	public static  final float[] eigenvalues(double[] eigen) {
		float[] val = new float[3];
		for (int n=0;n<3;n++) val[n] = (float)eigen[n];
		return val;
	}
	public static  final float principalEigenvalue(double[] eigen) {
		return (float)eigen[0];
	}
	public static  final float secondEigenvalue(double[] eigen) {
		return (float)eigen[1];
	}
	public static  final float thirdEigenvalue(double[] eigen) {
		return (float)eigen[2];
	}
	public static  final float principalFraction(double[] eigen) {
		return (float)(eigen[0]/(eigen[0]+eigen[1]+eigen[2]));
	}
	public static  final float secondFraction(double[] eigen) {
		return (float)(eigen[1]/(eigen[0]+eigen[1]+eigen[2]));
	}
	public static  final float thirdFraction(double[] eigen) {
		return (float)(eigen[2]/(eigen[0]+eigen[1]+eigen[2]));
	}
	public static  final float absolutePrincipalFraction(double[] eigen) {
		return (float)(Numerics.abs(eigen[0])/(Numerics.max(ZERO,Numerics.abs(eigen[0])+Numerics.abs(eigen[1])+Numerics.abs(eigen[2]))));
	}
	public static  final float absoluteSecondFraction(double[] eigen) {
		return (float)(Numerics.abs(eigen[1])/(Numerics.max(ZERO,Numerics.abs(eigen[0])+Numerics.abs(eigen[1])+Numerics.abs(eigen[2]))));
	}
	public static  final float absoluteThirdFraction(double[] eigen) {
		return (float)(Numerics.abs(eigen[2])/(Numerics.max(ZERO,Numerics.abs(eigen[0])+Numerics.abs(eigen[1])+Numerics.abs(eigen[2]))));
	}
	
	/**
	 *	check for positivity of tensor (needed in Log-Euclidean formulas)
	 */
	 public static final boolean isPositive(double[][] tensor) {
		 double[] lambda = Matrix3D.eigenvalues(tensor);
		 if (lambda[0]<=0 || lambda[1]<=0 || lambda[2]<=0) return false;
		 else return true;
	 }
	 
	/**
	 *	check for positivity of tensor (needed in Log-Euclidean formulas)
	 */
	 public static final boolean isPositive(double[] eigen) {
		 if (eigen[0]<=0 || eigen[1]<=0 || eigen[2]<=0) return false;
		 else return true;
	 }
	 
	/**
	 *	enforce positivity of tensor (arbitrary normalization, not good at preserving orientation)
	 */
	 public static final double[][] makePositive(double[][] tensor) {
		 double[] lambda = Matrix3D.eigenvalues(tensor);
		 double[][] vect = Matrix3D.eigenvectors(tensor,lambda);
		 lambda[0] = Numerics.abs(lambda[0]);
		 lambda[1] = Numerics.abs(lambda[1]);
		 lambda[2] = Numerics.abs(lambda[2]);
		 
		 double[][] newtensor = new double[3][3];
		 for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
			newtensor[n][m] = 0.0f;
			for (int p=0;p<3;p++) newtensor[n][m] += vect[p][n]*lambda[p]*vect[p][m];
		}
		return newtensor;
	 }
	 
	/**
	 *	compute the tensor distance in Log-Euclidean framework
	 */
	 public static final float logEuclideanDistance(double[] eigen1, double[] eigen2) {
		float[] lev1 = Tensor.logEuclideanVector(eigen1);
		float[] lev2 = Tensor.logEuclideanVector(eigen2);
		
		float frobenius = 0.0f;
		for (int i=0;i<6;i++) frobenius += (lev1[i]-lev2[i])*(lev1[i]-lev2[i]);

		return frobenius;
	 }
	 
	/**
	 *	compute the tensor distance in Log-Euclidean framework
	 */
	 public static final float logEuclideanDistance(double[][] t1, double[][] t2) {
		float[] lev1 = Tensor.logEuclideanVector(t1);
		float[] lev2 = Tensor.logEuclideanVector(t2);
		
		float frobenius = 0.0f;
		for (int i=0;i<6;i++) frobenius += (lev1[i]-lev2[i])*(lev1[i]-lev2[i]);

		return frobenius;
	 }
	 
	/**
	 *	compute the KL-divergence tensor distance
	 */
	 public static final float KLdivergenceDistance(double[] eigen1, double[] eigen2) {
		float div = 0.0f;
		
		double[][] t1 = eigenTensor(eigen1);
		double[][] t2 = eigenTensor(eigen2);
		double[][] it1 = inverseEigenTensor(eigen1);
		double[][] it2 = inverseEigenTensor(eigen2);
		
		float[] mtx = new float[3];
		for (int n=0;n<3;n++) {
			mtx[n] = -2.0f;
			
			for (int p=0;p<3;p++) mtx[n] += (float)(t1[n][p]*it2[p][n] + it1[n][p]*t2[p][n]);
				
		}
		div = 0.5f*(mtx[0]+mtx[1]+mtx[2]);
		
		return div;
	 }
	 
	 /**
	 *	compute the KL-divergence tensor distance
	 */
	 public static final float KLdivergenceDistance(double[][] t1, double[][] t2) {
		float div = 0.0f;
		
		double[][] it1 = inverseEigenTensor(eigenDecomposition(t1));
		double[][] it2 = inverseEigenTensor(eigenDecomposition(t2));
		
		float[] mtx = new float[3];
		for (int n=0;n<3;n++) {
			mtx[n] = -2.0f;
			
			for (int p=0;p<3;p++) mtx[n] += (float)(t1[n][p]*it2[p][n] + it1[n][p]*t2[p][n]);
				
		}
		div = 0.5f*(mtx[0]+mtx[1]+mtx[2]);
		
		return div;
	 }
	 	 
}
