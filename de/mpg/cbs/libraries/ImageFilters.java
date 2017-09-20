package de.mpg.cbs.libraries;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This class computes various basic image filters.
 *	
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ImageFilters {
	
	// no data: used as a library of functions
	
    // simple image functions
    
    public static final byte X = 0;
    public static final byte Y = 1;
    public static final byte Z = 2;
    public static final byte T = 3;

	/**
	 *	convolution
	 */
	public static float[][][] convolution(float[][][] image, int nx, int ny, int nz, float[][][] kernel, int kx, int ky, int kz) {
		float[][][] result = new float[nx][ny][nz];
		int xi,yj,zl;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
				xi = x+i; yj = y+j; zl = z+l;
				if ( (xi>=0) && (xi<nx) && (yj>=0) && (yj<ny) && (zl>=0) && (zl<nz) ) {
					result[x][y][z] += image[xi][yj][zl]*kernel[kx+i][ky+j][kz+l];
				}
			}
		}

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	 */
	public static float[][][] separableConvolution(float[][][] image, int nx, int ny, int nz, float[][] kernel, int kx, int ky, int kz) {
		float[][][] result = new float[nx][ny][nz];
		float[][][] temp = new float[nx][ny][nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[x][y][z] += image[x+i][y][z]*kernel[0][kx+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			temp[x][y][z] = 0.0f;	
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[x][y][z] += result[x][y+i][z]*kernel[1][ky+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[x][y][z] += temp[x][y][z+i]*kernel[2][kz+i];
				}
			}
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	 */
	public static float[] separableConvolution(float[] image, int nx, int ny, int nz, float[][] kernel, int kx, int ky, int kz) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[xyz] += image[xyz+i]*kernel[0][kx+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			temp[xyz] = 0.0f;	
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[xyz] += result[xyz+i*nx]*kernel[1][ky+i];
				}
			}
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[xyz] += temp[xyz+i*nx*ny]*kernel[2][kz+i];
				}
			}
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	 */
	public static float[] separableConvolution(float[] image, int nx, int ny, int nz, float[][] kernel) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		
		int kx = (kernel[X].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[xyz] += image[xyz+i]*kernel[0][kx+i];
				}
			}
		}
		int ky = (kernel[Y].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			temp[xyz] = 0.0f;	
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[xyz] += result[xyz+i*nx]*kernel[1][ky+i];
				}
			}
		}
		int kz = (kernel[Z].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[xyz] += temp[xyz+i*nx*ny]*kernel[2][kz+i];
				}
			}
		}
		temp = null;

		return result;
	}
		
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	*	this method compensates for boundaries 
	 */
	public static float[] separableBoundedConvolution(float[] image, int nx, int ny, int nz, float[][] kernel) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		int kx = (kernel[X].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;
			float den = 0.0f;
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) ) {
					result[xyz] += image[xyz+i]*kernel[X][kx+i];
					den += kernel[X][kx+i];
				}
			}
			result[xyz] /= den;
		}
		int ky = (kernel[Y].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			temp[xyz] = 0.0f;	
			float den = 0.0f;
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) ) {
					temp[xyz] += result[xyz+i*nx]*kernel[Y][ky+i];
					den += kernel[Y][ky+i];
				}
			}
			result[xyz] /= den;
		}
		int kz = (kernel[Z].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			float den = 0.0f;
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) ) {
					result[xyz] += temp[xyz+i*nx*ny]*kernel[Z][kz+i];
					den += kernel[Z][kz+i];
				}
			}
			result[xyz] /= den;
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	*	this method compensates for masked regions 
	 */
	public static float[] separableMaskedConvolution(float[] image, boolean[] mask, int nx, int ny, int nz, float[][] kernel) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		int kx = (kernel[X].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;
			double num = 0.0;
			double den = 0.0;
			
			if (mask[xyz]){
				for (int i=-kx;i<=kx;i++) {
					if ( (x+i>=0) && (x+i<nx) && mask[xyz+i] ) {
						num += image[xyz+i]*kernel[X][kx+i];
						den += kernel[X][kx+i];
					}
				}
				if (den*den>0) num /= den;
				result[xyz] = (float)num;
			}
		}
		
		int ky = (kernel[Y].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			temp[xyz] = 0.0f;
			double num = 0.0;
			double den = 0.0;
			
			if (mask[xyz]){
				for (int i=-ky;i<=ky;i++) {
					if ( (y+i>=0) && (y+i<ny) && mask[xyz+nx*i] ) {
						num += result[xyz+i*nx]*kernel[Y][ky+i];
						den += kernel[Y][ky+i];
					}
				}
				if (den*den>0) num /= den;
				temp[xyz] = (float)num;
			}
		}
		
		int kz = (kernel[Z].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			double num = 0.0;
			double den = 0.0;
			
			if (mask[xyz]){
				for (int i=-kz;i<=kz;i++) {
					if ( (z+i>=0) && (z+i<nz) && mask[xyz+nx*ny*i] ) {
						num += temp[xyz+i*nx*ny]*kernel[Z][kz+i];
						den += kernel[Z][kz+i];
					}
				}
				if (den*den>0) num /= den;
				result[xyz] = (float)num;
			}
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	*	this method compensates for masked regions 
	 */
	public static float[][][] separableMaskedConvolution(float[][][] image, boolean[][][] mask, int nx, int ny, int nz, float[][] kernel) {
		float[][][] result = new float[nx][ny][nz];
		float[][][] temp = new float[nx][ny][nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		int kx = (kernel[X].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;
			double num = 0.0;
			double den = 0.0;
			
			if (mask[x][y][z]){
				for (int i=-kx;i<=kx;i++) {
					if ( (x+i>=0) && (x+i<nx) && mask[x+i][y][z] ) {
						num += image[x+i][y][z]*kernel[X][kx+i];
						den += kernel[X][kx+i];
					}
				}
				if (den*den>0) num /= den;
				result[x][y][z] = (float)num;
			}
		}
		
		int ky = (kernel[Y].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			temp[x][y][z] = 0.0f;
			double num = 0.0;
			double den = 0.0;
			
			if (mask[x][y][z]){
				for (int i=-ky;i<=ky;i++) {
					if ( (y+i>=0) && (y+i<ny) && mask[x][y+i][z] ) {
						num += result[x][y+i][z]*kernel[Y][ky+i];
						den += kernel[Y][ky+i];
					}
				}
				if (den*den>0) num /= den;
				temp[x][y][z] = (float)num;
			}
		}
		
		int kz = (kernel[Z].length-1)/2;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = 0.0f;	
			double num = 0.0;
			double den = 0.0;
			
			if (mask[x][y][z]){
				for (int i=-kz;i<=kz;i++) {
					if ( (z+i>=0) && (z+i<nz) && mask[x][y][z+i] ) {
						num += temp[x][y][z+i]*kernel[Z][kz+i];
						den += kernel[Z][kz+i];
					}
				}
				if (den*den>0) num /= den;
				result[x][y][z] = (float)num;
			}
		}
		temp = null;

		return result;
	}
		
	/**
	*	convolution with a separable kernel (the kernel is 3x{kx,ky,kz})
	*	this method compensates for masked regions 
	 */
	public static float[] separableMaskedConvolution(float[] image, boolean[] mask, int nx, int ny, int nz, float[][] kernel, int kx, int ky, int kz) {
		float[] result = new float[nx*ny*nz];
		float[] temp = new float[nx*ny*nz];
		int xi,yj,zl;
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			double num = 0.0;
			double den = 0.0;
			for (int i=-kx;i<=kx;i++) {
				if ( (x+i>=0) && (x+i<nx) && mask[xyz+i] ) {
					num += image[xyz+i]*kernel[X][kx+i];
					den += kernel[X][kx+i];
				}
			}
			if (den*den>0) num /= den;
			result[xyz] = (float)num;
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			double num = 0.0;
			double den = 0.0;
			for (int i=-ky;i<=ky;i++) {
				if ( (y+i>=0) && (y+i<ny) && mask[xyz+nx*i] ) {
					num += result[xyz+i*nx]*kernel[Y][ky+i];
					den += kernel[Y][ky+i];
				}
			}
			if (den*den>0) num /= den;
			temp[xyz] = (float)num;
		}
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;	
			double num = 0.0;
			double den = 0.0;
			for (int i=-kz;i<=kz;i++) {
				if ( (z+i>=0) && (z+i<nz) && mask[xyz+nx*ny*i] ) {
					num += temp[xyz+i*nx*ny]*kernel[Z][kz+i];
					den += kernel[Z][kz+i];
				}
			}
			if (den*den>0) num /= den;
			result[xyz] = (float)num;
		}
		temp = null;

		return result;
	}
		
	/**
	 *	Gaussian kernel
	 */
	public static float[][][] gaussianKernel(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Numerics.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Numerics.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Numerics.max(3.0f*sz-0.5f,0.0f));
		
		// create the kernel
		float[][][] kernel = new float[2*kx+1][2*ky+1][2*kz+1];
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
			kernel[kx+i][ky+j][kz+l] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) + (j*j)/(sy*sy) + (l*l)/(sz*sz) );
			sum += kernel[kx+i][ky+j][kz+l];
		}
		// normalize
		for (int i=-kx;i<=kx;i++) for (int j=-ky;j<=ky;j++) for (int l=-kz;l<=kz;l++) {
			kernel[kx+i][ky+j][kz+l] = kernel[kx+i][ky+j][kz+l]/sum;
		}

		return kernel;
	}
	/**
	 *	Gaussian kernel for separable convolution
	 */
	public static float[][] separableGaussianKernel(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[0] = new float[2*kx+1]; 
		kernel[1] = new float[2*ky+1]; 
		kernel[2] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[0][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[0][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[0][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[1][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[1][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[1][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[2][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[2][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[2][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivX(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		// derivative
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] *= -i/(sx*sx);
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivY(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		// derivative
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] *= -j/(sy*sy);
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivZ(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		// derivative
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] *= -l/(sz*sz);
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivX2(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		// second derivative
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] *= (i*i/(sx*sx)-1.0f)/(sx*sx);
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivY2(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		// second derivative
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] *= (j*j/(sy*sy)-1.0f)/(sy*sy);
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivZ2(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		// second derivative
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] *= (l*l/(sz*sz)-1.0f)/(sz*sz);
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivXY(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		// derivative
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] *= -i/(sx*sx);
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		// derivative
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] *= -j/(sy*sy);
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivYZ(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		// derivative
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] *= -j/(sy*sy);
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		// derivative
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] *= -l/(sz*sz);
		
		return kernel;
	}
	
	public static float[][] separableGaussianKernelDerivZX(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[X] = new float[2*kx+1]; 
		kernel[Y] = new float[2*ky+1]; 
		kernel[Z] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[X][kx+i] = (float)Math.exp( - 0.5f*(i*i)/(sx*sx) );
			//MedicUtilPublic.displayMessage("exp("+( - 0.5f*(i*i)/(sx*sx) )+") = "+kernel[0][kx+i]+"\n");
			sum += kernel[X][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[Y][ky+j] = (float)Math.exp( - 0.5f*(j*j)/(sy*sy) );
			sum += kernel[Y][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[Y][ky+j] /= sum;
		// derivative
		for (int i=-kx;i<=kx;i++) kernel[X][kx+i] *= -i/(sx*sx);
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[Z][kz+l] = (float)Math.exp( - 0.5f*(l*l)/(sz*sz) );
			sum += kernel[Z][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] /= sum;
		// derivative
		for (int l=-kz;l<=kz;l++) kernel[Z][kz+l] *= -l/(sz*sz);
		
		return kernel;
	}

	public static float[][][] surfaceEdgeFilter(float[][][] image, boolean[][][] mask, int nx, int ny, int nz) {
		float[][][] result = new float[nx][ny][nz];
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			result[x][y][z] = 0.0f;
			
			if (mask[x][y][z]) {
				float val, best = 0.0f;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z]-image[x+1][y][z]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z]-image[x+1][y-1][z]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z]-image[x+1][y+1][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y][z-1]-image[x+1][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y][z+1]-image[x+1][y][z+1]
					 +2.0f*image[x][y-1][z-1]-image[x-1][y-1][z-1]-image[x+1][y-1][z-1]
					 +2.0f*image[x][y-1][z+1]-image[x-1][y-1][z+1]-image[x+1][y-1][z+1]
					 +2.0f*image[x][y+1][z-1]-image[x-1][y+1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x][y+1][z+1]-image[x-1][y+1][z+1]-image[x+1][y+1][z+1])/9.0f;			
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z]-image[x][y+1][z]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z]-image[x-1][y+1][z]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z]-image[x+1][y+1][z]
					 +2.0f*image[x][y][z-1]-image[x][y-1][z-1]-image[x][y+1][z-1]
					 +2.0f*image[x][y][z+1]-image[x][y-1][z+1]-image[x][y+1][z+1]
					 +2.0f*image[x-1][y][z-1]-image[x-1][y-1][z-1]-image[x-1][y+1][z-1]
					 +2.0f*image[x-1][y][z+1]-image[x-1][y-1][z+1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y][z-1]-image[x+1][y-1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x+1][y][z+1]-image[x+1][y-1][z+1]-image[x+1][y+1][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y][z-1]-image[x][y][z+1]
					 +2.0f*image[x-1][y][z]-image[x-1][y][z-1]-image[x-1][y][z+1]
					 +2.0f*image[x+1][y][z]-image[x+1][y][z-1]-image[x+1][y][z+1]
					 +2.0f*image[x][y-1][z]-image[x][y-1][z-1]-image[x][y-1][z+1]
					 +2.0f*image[x][y+1][z]-image[x][y+1][z-1]-image[x][y+1][z+1]
					 +2.0f*image[x-1][y-1][z]-image[x-1][y-1][z-1]-image[x-1][y-1][z+1]
					 +2.0f*image[x-1][y+1][z]-image[x-1][y+1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y-1][z]-image[x+1][y-1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x+1][y+1][z]-image[x+1][y+1][z-1]-image[x+1][y+1][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z]-image[x+1][y+1][z]
					 +2.0f*image[x-1][y+1][z]-image[x-2][y][z]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z]-image[x][y-2][z]-image[x+2][y][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y-1][z-1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]-image[x][y+2][z-1]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y-1][z+1]-image[x+1][y+1][z+1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]-image[x][y+2][z+1]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z-1]-image[x][y+1][z+1]
					 +2.0f*image[x][y+1][z-1]-image[x][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x][y-1][z+1]-image[x][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]-image[x-1][y+2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]-image[x-1][y][z+2]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]-image[x+1][y+2][z]
					 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]-image[x+1][y][z+2])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z-1]-image[x+1][y][z+1]
					 +2.0f*image[x+1][y][z-1]-image[x][y][z-2]-image[x+2][y][z]
					 +2.0f*image[x-1][y][z+1]-image[x-2][y][z]-image[x][y][z+2]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]-image[x+2][y-1][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]-image[x][y-1][z+2]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]-image[x+2][y+1][z]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]-image[x][y+1][z+2])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y+1][z]-image[x+1][y-1][z]
					 +2.0f*image[x-1][y-1][z]-image[x-2][y][z]-image[x][y-2][z]
					 +2.0f*image[x+1][y+1][z]-image[x][y-2][z]-image[x+2][y][z]
					 +2.0f*image[x][y][z-1]-image[x-1][y+1][z-1]-image[x+1][y-1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]-image[x][y-2][z-1]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]-image[x+2][y][z-1]
					 +2.0f*image[x][y][z+1]-image[x-1][y+1][z+1]-image[x+1][y-1][z+1]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]-image[x][y-2][z+1]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]-image[x+2][y][z+1])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x][y-1][z+1]-image[x][y+1][z-1]
					 +2.0f*image[x][y-1][z-1]-image[x][y-2][z]-image[x][y][z-2]
					 +2.0f*image[x][y+1][z+1]-image[x][y][z+2]-image[x][y+2][z]
					 +2.0f*image[x-1][y][z]-image[x-1][y-1][z+1]-image[x-1][y+1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]-image[x-1][y][z-2]
					 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]-image[x-1][y+2][z]
					 +2.0f*image[x+1][y][z]-image[x+1][y-1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]-image[x+1][y][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]-image[x+1][y+2][z])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y][z+1]-image[x+1][y][z-1]
					 +2.0f*image[x-1][y][z-1]-image[x-2][y][z]-image[x][y][z-2]
					 +2.0f*image[x+1][y][z+1]-image[x][y][z+2]-image[x+2][y][z]
					 +2.0f*image[x][y-1][z]-image[x-1][y-1][z+1]-image[x+1][y-1][z-1]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]-image[x][y-1][z-2]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]-image[x+2][y-1][z]
					 +2.0f*image[x][y+1][z]-image[x-1][y+1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]-image[x][y+1][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]-image[x+2][y+1][z])/9.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z-1]-image[x+1][y+1][z+1]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]-image[x+2][y][z]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]-image[x+2][y+2][z]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]-image[x][y+2][z+2]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]-image[x+2][y][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x+1][y-1][z-1]-image[x-1][y+1][z+1]
					 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]-image[x][y][z+2]
					 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]-image[x][y+2][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]-image[x-2][y][z]
					 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]-image[x-2][y+2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]-image[x-2][y][z+2]
					 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]-image[x][y+2][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y+1][z-1]-image[x+1][y-1][z+1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]-image[x][y][z+2]
					 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]-image[x+2][y][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]-image[x][y-2][z]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]-image[x+2][y-2][z]
					 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]-image[x][y-2][z+2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]-image[x+2][y][z+2])/7.0f;
				if (val*val>best*best) best = val;
				
				val =(2.0f*image[x][y][z]-image[x-1][y-1][z+1]-image[x+1][y+1][z-1]
					 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]-image[x][y+2][z]
					 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]-image[x+2][y][z]
					 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]-image[x][y][z-2]
					 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]-image[x+2][y][z-2]
					 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]-image[x][y+2][z-2]
					 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]-image[x+2][y+2][z])/7.0f;
				if (val*val>best*best) best = val;
				
				result[x][y][z] = 0.5f*best;
			}
		}
		return result;
	}
	
	public static float[] surfaceEdgeFilter(float[] image, boolean[] mask, int nx, int ny, int nz) {
		float[] result = new float[nx*ny*nz];
		
		int X = 1;
		int Y = nx;
		int Z = nx*ny;
		int X2 = 2;
		int Y2 = 2*nx;
		int Z2 = 2*nx*ny;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x+nx*y+nx*ny*z;
			result[xyz] = 0.0f;
			
			boolean process = true;
			for (int i=-2;i<=2 && process;i++) for (int j=-2;j<=2 && process;j++) for (int l=-2;l<=2 && process;l++) {
				if (i*i+j*j+l*l<=8) {
					process = mask[xyz+i+j*nx+l*nx*ny];
				}
			}
			if (process) {
				float val, best = 0.0f, second = 0.0f;
				
				val =(2.0f*image[xyz]		-image[xyz-X]		-image[xyz+X]
					 +2.0f*image[xyz-Y]		-image[xyz-X-Y]		-image[xyz+X-Y]
					 +2.0f*image[xyz+Y]		-image[xyz-X+Y]		-image[xyz+X+Y]
					 +2.0f*image[xyz-Z]		-image[xyz-X-Z]		-image[xyz+X-Z]
					 +2.0f*image[xyz+Z]		-image[xyz-X+Z]		-image[xyz+X+Z]
					 +2.0f*image[xyz-Y-Z]	-image[xyz-X-Y-Z]	-image[xyz+X-Y-Z]
					 +2.0f*image[xyz-Y+Z]	-image[xyz-X-Y+Z]	-image[xyz+X-Y+Z]
					 +2.0f*image[xyz+Y-Z]	-image[xyz-X+Y-Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz+Y+Z]	-image[xyz-X+Y+Z]	-image[xyz+X+Y+Z])/9.0f;			
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}

				val =(2.0f*image[xyz]		-image[xyz-Y]		-image[xyz+Y]
					 +2.0f*image[xyz-X]		-image[xyz-X-Y]		-image[xyz-X+Y]
					 +2.0f*image[xyz+X]		-image[xyz+X-Y]		-image[xyz+X+Y]
					 +2.0f*image[xyz-Z]		-image[xyz-Y-Z]		-image[xyz+Y-Z]
					 +2.0f*image[xyz+Z]		-image[xyz-Y+Z]		-image[xyz+Y+Z]
					 +2.0f*image[xyz-X-Z]	-image[xyz-X-Y-Z]	-image[xyz-X+Y-Z]
					 +2.0f*image[xyz-X+Z]	-image[xyz-X-Y+Z]	-image[xyz-X+Y+Z]
					 +2.0f*image[xyz+X-Z]	-image[xyz+X-Y-Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz+X+Z]	-image[xyz+X-Y+Z]	-image[xyz+X+Y+Z])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]		-image[xyz-Z]		-image[xyz+Z]
					 +2.0f*image[xyz-X]		-image[xyz-X-Z]		-image[xyz-X+Z]
					 +2.0f*image[xyz+X]		-image[xyz+X-Z]		-image[xyz+X+Z]
					 +2.0f*image[xyz-Y]		-image[xyz-Y-Z]		-image[xyz-Y+Z]
					 +2.0f*image[xyz+Y]		-image[xyz+Y-Z]		-image[xyz+Y+Z]
					 +2.0f*image[xyz-X-Y]	-image[xyz-X-Y-Z]	-image[xyz-X-Y+Z]
					 +2.0f*image[xyz-X+Y]	-image[xyz-X+Y-Z]	-image[xyz-X+Y+Z]
					 +2.0f*image[xyz+X-Y]	-image[xyz+X-Y-Z]	-image[xyz+X-Y+Z]
					 +2.0f*image[xyz+X+Y]	-image[xyz+X+Y-Z]	-image[xyz+X+Y+Z])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X-Y]		-image[xyz+X+Y]
					 +2.0f*image[xyz-X+Y]		-image[xyz-X2]		-image[xyz+Y2]
					 +2.0f*image[xyz+X-Y]		-image[xyz-Y2]		-image[xyz+X2]
					 +2.0f*image[xyz-Z]			-image[xyz-X-Y-Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-X2-Z]	-image[xyz+Y2-Z]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz-Y2-Z]	-image[xyz+X2-Z]
					 +2.0f*image[xyz+Z]			-image[xyz-X-Y+Z]	-image[xyz+X+Y+Z]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X2+Z]	-image[xyz+Y2+Z]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz-Y2+Z]	-image[xyz+X2+Z])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-Y-Z]		-image[xyz+Y+Z]
					 +2.0f*image[xyz+Y-Z]		-image[xyz-Z2]		-image[xyz+Y2]
					 +2.0f*image[xyz-Y+Z]		-image[xyz-Y2]		-image[xyz+Z2]
					 +2.0f*image[xyz-X]			-image[xyz-X-Y-Z]	-image[xyz-X+Y+Z]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-X-Z2]	-image[xyz-X+Y2]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-X-Y2]	-image[xyz-X+Z2]
					 +2.0f*image[xyz+X]			-image[xyz+X-Y-Z]	-image[xyz+X+Y+Z]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz+X-Z2]	-image[xyz+X+Z2]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz+X-Y2]	-image[xyz+X+Z2])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X-Z]		-image[xyz+X+Z]
					 +2.0f*image[xyz+X-Z]		-image[xyz-Z2]		-image[xyz+X2]
					 +2.0f*image[xyz-X+Z]		-image[xyz-X2]		-image[xyz+Z2]
					 +2.0f*image[xyz-Y]			-image[xyz-X-Y-Z]	-image[xyz+X-Y+Z]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz-Y-Z2]	-image[xyz+X2-Y]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-X2-Y]	-image[xyz-Y+Z2]
					 +2.0f*image[xyz+Y]			-image[xyz-X+Y-Z]	-image[xyz+X+Y+Z]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz+Y-Z2]	-image[xyz+X2+Y]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X2+Y]	-image[xyz+Y+Z2])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X+Y]		-image[xyz+X-Y]
					 +2.0f*image[xyz-X-Y]		-image[xyz-X2]		-image[xyz-Y2]
					 +2.0f*image[xyz+X+Y]		-image[xyz-Y2]		-image[xyz+X2]
					 +2.0f*image[xyz-Z]			-image[xyz-X+Y-Z]	-image[xyz+X-Y-Z]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-X2-Z]	-image[xyz-Y2-Z]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz-Y2-Z]	-image[xyz+X2-Z]
					 +2.0f*image[xyz+Z]			-image[xyz-X+Y+Z]	-image[xyz+X-Y+Z]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-X2+Z]	-image[xyz-Y2+Z]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz-Y2+Z]	-image[xyz+X2+Z])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-Y+Z]		-image[xyz+Y-Z]
					 +2.0f*image[xyz-Y-Z]		-image[xyz-Y2]		-image[xyz-Z2]
					 +2.0f*image[xyz+Y+Z]		-image[xyz+Z2]		-image[xyz+Y2]
					 +2.0f*image[xyz-X]			-image[xyz-X-Y+Z]	-image[xyz-X+Y-Z]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-X-Y2]	-image[xyz-X-Z2]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X+Z2]	-image[xyz-X+Y2]
					 +2.0f*image[xyz+X]			-image[xyz+X-Y+Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz+X-Y2]	-image[xyz+X-Z2]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz+X+Z2]	-image[xyz+X+Y2])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X+Z]		-image[xyz+X-Z]
					 +2.0f*image[xyz-X-Z]		-image[xyz-X2]		-image[xyz-Z2]
					 +2.0f*image[xyz+X+Z]		-image[xyz+Z2]		-image[xyz+X2]
					 +2.0f*image[xyz-Y]			-image[xyz-X-Y+Z]	-image[xyz+X-Y-Z]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-X2-Y]	-image[xyz-Y-Z2]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz-Y+Z2]	-image[xyz+X2-Y]
					 +2.0f*image[xyz+Y]			-image[xyz-X+Y+Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-X2+Y]	-image[xyz+Y-Z2]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz+Y+Z2]	-image[xyz+X2+Y])/9.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X-Y-Z]	-image[xyz+X+Y+Z]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-X2-Y2]	-image[xyz+Z2]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-X2-Z2]	-image[xyz+Y2]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz-Y2-Z2]	-image[xyz+X2]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz-Z2]		-image[xyz+X2+Y2]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X2]		-image[xyz+Y2+Z2]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz-Y2]		-image[xyz+X2+Z2])/7.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz+X-Y-Z]	-image[xyz-X+Y+Z]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz+X2-Y2]	-image[xyz+Z2]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz+X2-Z2]	-image[xyz+Y2]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-Y2-Z2]	-image[xyz-X2]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-Z2]		-image[xyz-X2+Y2]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-Y2]		-image[xyz-X2+Z2]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz+X2]		-image[xyz+Y2+Z2])/7.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X+Y-Z]	-image[xyz+X-Y+Z]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X2+Y2]	-image[xyz+Z2]
					 +2.0f*image[xyz+X+Y-Z]		-image[xyz+Y2-Z2]	-image[xyz+X2]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-X2-Z2]	-image[xyz-Y2]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz-Z2]		-image[xyz+X2-Y2]
					 +2.0f*image[xyz-X-Y+Z]		-image[xyz-X2]		-image[xyz-Y2+Z2]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz+Y2]		-image[xyz+X2+Z2])/7.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(2.0f*image[xyz]			-image[xyz-X-Y+Z]	-image[xyz+X+Y-Z]
					 +2.0f*image[xyz-X+Y+Z]		-image[xyz-X2+Z2]	-image[xyz+Y2]
					 +2.0f*image[xyz+X-Y+Z]		-image[xyz-Y2+Z2]	-image[xyz+X2]
					 +2.0f*image[xyz-X-Y-Z]		-image[xyz-X2-Y2]	-image[xyz-Z2]
					 +2.0f*image[xyz+X-Y-Z]		-image[xyz-Y2]		-image[xyz+X2-Z2]
					 +2.0f*image[xyz-X+Y-Z]		-image[xyz-X2]		-image[xyz+Y2-Z2]
					 +2.0f*image[xyz+X+Y+Z]		-image[xyz+Z2]		-image[xyz+X2+Y2])/7.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				//result[xyz] = 0.25f*best*(1.0f - Numerics.abs(second/best));
				result[xyz] = 0.25f*best;
			}
		}
		return result;
	}
	
	public static final float planeScore(float[][][] image, int x, int y, int z, int d) {
		float val = 0.0f;
		if (d==0) {		
			val =(2.0f*image[x][y][z]		-image[x-1][y][z]		-image[x+1][y][z]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z-1]	-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x][y-1][z+1]	-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x][y+1][z-1]	-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x][y+1][z+1]	-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==1) {
			val =(2.0f*image[x][y][z]		-image[x][y-1][z]		-image[x][y+1][z]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1]
				 +2.0f*image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1]
				 +2.0f*image[x-1][y][z-1]	-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y][z+1]	-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y][z-1]	-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y][z+1]	-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==2) { 			
			val =(2.0f*image[x][y][z]		-image[x][y][z-1]		-image[x][y][z+1]
				 +2.0f*image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1]
				 +2.0f*image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1]
				 +2.0f*image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1]
				 +2.0f*image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1]
				 +2.0f*image[x-1][y-1][z]	-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z]	-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z]	-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y+1][z]	-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1])/18.0f;
		} else if (d==3) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z]		-image[x+1][y+1][z]
				 +2.0f*image[x-1][y+1][z]	-image[x-2][y][z]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z]	-image[x][y-2][z]		-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-1]		-image[x][y+2][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-1]		-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+1]		-image[x][y+2][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+1]		-image[x+2][y][z+1])/18.0f;
		} else if (d==4) { 			
			val =(2.0f*image[x][y][z]		-image[x][y-1][z-1]		-image[x][y+1][z+1]
				 +2.0f*image[x][y+1][z-1]	-image[x][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x][y-1][z+1]	-image[x][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-1][y][z-2]		-image[x-1][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-1][y-2][z]		-image[x-1][y][z+2]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x+1][y][z-2]		-image[x+1][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x+1][y-2][z]		-image[x+1][y][z+2])/18.0f;
		} else if (d==5) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y][z-1]		-image[x+1][y][z+1]
				 +2.0f*image[x+1][y][z-1]	-image[x][y][z-2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y][z+1]	-image[x-2][y][z]		-image[x][y][z+2]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-1][z-2]		-image[x+2][y-1][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-1][z]		-image[x][y-1][z+2]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+1][z-2]		-image[x+2][y+1][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+1][z]		-image[x][y+1][z+2])/18.0f;
		} else if (d==6) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y+1][z]		-image[x+1][y-1][z]
				 +2.0f*image[x-1][y-1][z]	-image[x-2][y][z]		-image[x][y-2][z]
				 +2.0f*image[x+1][y+1][z]	-image[x][y-2][z]		-image[x+2][y][z]
				 +2.0f*image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-1]		-image[x][y-2][z-1]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y-2][z-1]		-image[x+2][y][z-1]
				 +2.0f*image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z+1]		-image[x][y-2][z+1]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y-2][z+1]		-image[x+2][y][z+1])/18.0f;
		} else if (d==7) { 			
			val =(2.0f*image[x][y][z]		-image[x][y-1][z+1]		-image[x][y+1][z-1]
				 +2.0f*image[x][y-1][z-1]	-image[x][y-2][z]		-image[x][y][z-2]
				 +2.0f*image[x][y+1][z+1]	-image[x][y][z+2]		-image[x][y+2][z]
				 +2.0f*image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-1][y-2][z]		-image[x-1][y][z-2]
				 +2.0f*image[x-1][y+1][z+1]-image[x-1][y][z+2]		-image[x-1][y+2][z]
				 +2.0f*image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x+1][y-1][z-1]-image[x+1][y-2][z]		-image[x+1][y][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+1][y][z+2]		-image[x+1][y+2][z])/18.0f;
		} else if (d==8) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y][z+1]		-image[x+1][y][z-1]
				 +2.0f*image[x-1][y][z-1]	-image[x-2][y][z]		-image[x][y][z-2]
				 +2.0f*image[x+1][y][z+1]	-image[x][y][z+2]		-image[x+2][y][z]
				 +2.0f*image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-1][z]		-image[x][y-1][z-2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-1][z+2]		-image[x+2][y-1][z]
				 +2.0f*image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y+1][z]		-image[x][y+1][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+1][z+2]		-image[x+2][y+1][z])/18.0f;
		} else if (d==9) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z+1]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z-2]		-image[x+2][y][z]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y][z-2]		-image[x+2][y+2][z]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z]		-image[x][y+2][z+2]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z]		-image[x+2][y][z+2])/14.0f;
		} else if (d==10) { 			
			val =(2.0f*image[x][y][z]		-image[x+1][y-1][z-1]	-image[x-1][y+1][z+1]
				 +2.0f*image[x+1][y-1][z+1]-image[x+2][y-2][z]		-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x+2][y][z-2]		-image[x][y+2][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x][y-2][z-2]		-image[x-2][y][z]
				 +2.0f*image[x-1][y+1][z-1]-image[x][y][z-2]		-image[x-2][y+2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x][y-2][z]		-image[x-2][y][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x+2][y][z]		-image[x][y+2][z+2])/14.0f;
		} else if (d==11) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z+1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y+2][z]		-image[x][y][z+2]
				 +2.0f*image[x+1][y+1][z-1]-image[x][y+2][z-2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y][z-2]		-image[x][y-2][z]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y][z-2]		-image[x+2][y-2][z]
				 +2.0f*image[x-1][y-1][z+1]-image[x-2][y][z]		-image[x][y-2][z+2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y+2][z]		-image[x+2][y][z+2])/14.0f;
		} else if (d==12) { 			
			val =(2.0f*image[x][y][z]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z-1]
				 +2.0f*image[x-1][y+1][z+1]-image[x-2][y][z+2]		-image[x][y+2][z]
				 +2.0f*image[x+1][y-1][z+1]-image[x][y-2][z+2]		-image[x+2][y][z]
				 +2.0f*image[x-1][y-1][z-1]-image[x-2][y-2][z]		-image[x][y][z-2]
				 +2.0f*image[x+1][y-1][z-1]-image[x][y-2][z]		-image[x+2][y][z-2]
				 +2.0f*image[x-1][y+1][z-1]-image[x-2][y][z]		-image[x][y+2][z-2]
				 +2.0f*image[x+1][y+1][z+1]-image[x][y][z+2]		-image[x+2][y+2][z])/14.0f;
		}
		return 0.5f*val;
	}

	public static float[] tubularEdgeFilter(float[] image, boolean[] mask, int nx, int ny, int nz) {
		float[] result = new float[nx*ny*nz];
		
		int X = 1;
		int Y = nx;
		int Z = nx*ny;
		int X2 = 2;
		int Y2 = 2*nx;
		int Z2 = 2*nx*ny;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x+nx*y+Z*z;
			result[xyz] = 0.0f;
			
			boolean process = true;
			for (int i=-2;i<=2 && process;i++) for (int j=-2;j<=2 && process;j++) for (int l=-2;l<=2 && process;l++) {
				process = mask[xyz+i+j*nx+l*nx*ny];
			}
			if (process) {
				float val, best = 0.0f, second = 0.0f;
				
				val =(8.0f*(image[xyz]			+image[xyz-X]		+image[xyz+X])
					 		-image[xyz-Y]		-image[xyz-X-Y]		-image[xyz+X-Y]
					 		-image[xyz+Y]		-image[xyz-X+Y]		-image[xyz+X+Y]
					 		-image[xyz-Z]		-image[xyz-Z-X]		-image[xyz-Z+X]
					 		-image[xyz+Z]		-image[xyz+Z-X]		-image[xyz+Z+X]
					 		-image[xyz-Y-Z]		-image[xyz-X-Y-Z]	-image[xyz+X-Y-Z]
					 		-image[xyz-Y+Z]		-image[xyz-X-Y+Z]	-image[xyz+X-Y+Z]
					 		-image[xyz+Y-Z]		-image[xyz-X+Y-Z]	-image[xyz+X+Y-Z]
					 		-image[xyz+Y+Z]		-image[xyz-X+Y+Z]	-image[xyz+X+Y+Z])/24.0f;			
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				val =(8.0f*(image[xyz]			+image[xyz-Y]		+image[xyz+Y])
					 		-image[xyz-X]		-image[xyz-X-Y]		-image[xyz-X+Y]
					 		-image[xyz+X]		-image[xyz+X-Y]		-image[xyz+X+Y]
					 		-image[xyz-Z]		-image[xyz-Y-Z]		-image[xyz+Y-Z]
					 		-image[xyz+Z]		-image[xyz-Y+Z]		-image[xyz+Y+Z]
					 		-image[xyz-Z-X]		-image[xyz-X-Y-Z]	-image[xyz-X+Y-Z]
							-image[xyz+Z-X]		-image[xyz-X-Y+Z]	-image[xyz-X+Y+Z]
							-image[xyz-Z+X]		-image[xyz+X-Y-Z]	-image[xyz+X+Y-Z]
							-image[xyz+Z+X]		-image[xyz+X-Y+Z]	-image[xyz+X+Y+Z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-Z]		+image[xyz+Z])
					 		-image[xyz-X]		-image[xyz-Z-X]		-image[xyz+Z-X]
					 		-image[xyz+X]		-image[xyz-Z+X]		-image[xyz+Z+X]
					 		-image[xyz-Y]		-image[xyz-Y-Z]		-image[xyz-Y+Z]
					 		-image[xyz+Y]		-image[xyz+Y-Z]		-image[xyz+Y+Z]
					 		-image[xyz-X-Y]		-image[xyz-X-Y-Z]	-image[xyz-X-Y+Z]
					 		-image[xyz-X+Y]		-image[xyz-X+Y-Z]	-image[xyz-X+Y+Z]
					 		-image[xyz+X-Y]		-image[xyz+X-Y-Z]	-image[xyz+X-Y+Z]
					 		-image[xyz+X+Y]		-image[xyz+X+Y-Z]	-image[xyz+X+Y+Z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-X-Y]		+image[xyz+X+Y])
					 		-image[xyz-X+Y]		-image[xyz-X2]		-image[xyz+Y2]
					 		-image[xyz+X-Y]		-image[xyz-Y2]		-image[xyz+X2]
					 		-image[xyz-Z]		-image[xyz-X-Y-Z]	-image[xyz+X+Y-Z]
					 		-image[xyz-X+Y-Z]	-image[xyz-X2-Z]	-image[xyz+Y2-Z]
					 		-image[xyz+X-Y-Z]	-image[xyz-Y2-Z]	-image[xyz+X2-Z]
					 		-image[xyz+Z]		-image[xyz-X-Y+Z]	-image[xyz+X+Y+Z]
					 		-image[xyz-X+Y+Z]	-image[xyz-X2+Z]	-image[xyz+Y2+Z]
					 		-image[xyz+X-Y+Z]	-image[xyz-Y2+Z]	-image[xyz+X2+Z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-Y-Z]		+image[xyz+Y+Z])
					 		-image[xyz+Y-Z]		-image[xyz-Z2]		-image[xyz+Y2]
					 		-image[xyz-Y+Z]		-image[xyz-Y2]		-image[xyz+Z2]
					 		-image[xyz-X]		-image[xyz-X-Y-Z]	-image[xyz-X+Y+Z]
					 		-image[xyz-X+Y-Z]	-image[xyz-X-Z2]	-image[xyz-X+Y2]
					 		-image[xyz-X-Y+Z]	-image[xyz-X-Y2]	-image[xyz-X+Z2]
					 		-image[xyz+X]		-image[xyz+X-Y-Z]	-image[xyz+X+Y+Z]
					 		-image[xyz+X+Y-Z]	-image[xyz+X-Z2]	-image[xyz+X+Z2]
					 		-image[xyz+X-Y+Z]	-image[xyz+X-Y2]	-image[xyz+X+Z2])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-Z-X]		+image[xyz+Z+X])
					 		-image[xyz-Z+X]		-image[xyz-Z2]		-image[xyz+X2]
					 		-image[xyz+Z-X]		-image[xyz-X2]		-image[xyz+Z2]
					 		-image[xyz-Y]		-image[xyz-X-Y-Z]	-image[xyz+X-Y+Z]
					 		-image[xyz+X-Y-Z]	-image[xyz-Y-Z2]	-image[xyz+X2-Y]
					 		-image[xyz-X-Y+Z]	-image[xyz-X2-Y]	-image[xyz-Y+Z2]
					 		-image[xyz+Y]		-image[xyz-X+Y-Z]	-image[xyz+X+Y+Z]
					 		-image[xyz+X+Y-Z]	-image[xyz+Y-Z2]	-image[xyz+X2+Y]
					 		-image[xyz-X+Y+Z]	-image[xyz-X2+Y]	-image[xyz+Y+Z2])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-X+Y]		+image[xyz+X-Y])
					 		-image[xyz-X-Y]		-image[xyz-X2]		-image[xyz-Y2]
					 		-image[xyz+X+Y]		-image[xyz-Y2]		-image[xyz+X2]
					 		-image[xyz-Z]		-image[xyz-X+Y-Z]	-image[xyz+X-Y-Z]
					 		-image[xyz-X-Y-Z]	-image[xyz-X2-Z]	-image[xyz-Y2-Z]
					 		-image[xyz+X+Y-Z]	-image[xyz-Y2-Z]	-image[xyz+X2-Z]
					 		-image[xyz+Z]		-image[xyz-X+Y+Z]	-image[xyz+X-Y+Z]
					 		-image[xyz-X-Y+Z]	-image[xyz-X2+Z]	-image[xyz-Y2+Z]
					 		-image[xyz+X+Y+Z]	-image[xyz-Y2+Z]	-image[xyz+X2+Z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz-Y+Z]		+image[xyz+Y-Z])
					 		-image[xyz-Y-Z]		-image[xyz-Y2]		-image[xyz-Z2]
					 		-image[xyz+Y+Z]		-image[xyz+Z2]		-image[xyz+Y2]
					 		-image[xyz-X]		-image[xyz-X-Y+Z]	-image[xyz-X+Y-Z]
					 		-image[xyz-X-Y-Z]	-image[xyz-X-Y2]	-image[xyz-X-Z2]
					 		-image[xyz-X+Y+Z]	-image[xyz-X+Z2]	-image[xyz-X+Y2]
					 		-image[xyz+X]		-image[xyz+X-Y+Z]	-image[xyz+X+Y-Z]
					 		-image[xyz+X-Y-Z]	-image[xyz+X-Y2]	-image[xyz+X-Z2]
					 		-image[xyz+X+Y+Z]	-image[xyz+X+Z2]	-image[xyz+X+Y2])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[xyz]			+image[xyz+Z-X]		+image[xyz-Z+X])
					 		-image[xyz-Z-X]		-image[xyz-X2]		-image[xyz-Z2]
					 		-image[xyz+Z+X]		-image[xyz+Z2]		-image[xyz+X2]
					 		-image[xyz-Y]		-image[xyz-X-Y+Z]	-image[xyz+X-Y-Z]
					 		-image[xyz-X-Y-Z]	-image[xyz-X2-Y]	-image[xyz-Y-Z2]
					 		-image[xyz+X-Y+Z]	-image[xyz-Y+Z2]	-image[xyz+X2-Y]
					 		-image[xyz+Y]		-image[xyz-X+Y+Z]	-image[xyz+X+Y-Z]
					 		-image[xyz-X+Y-Z]	-image[xyz-X2+Y]	-image[xyz+Y-Z2]
					 		-image[xyz+X+Y+Z]	-image[xyz+Y+Z2]	-image[xyz+X2+Y])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[xyz]			+image[xyz-X-Y-Z]	+image[xyz+X+Y+Z])
							-image[xyz-X-Y+Z]	-image[xyz-X2-Y2]	-image[xyz+Z2]
							-image[xyz-X+Y-Z]	-image[xyz-X2-Z2]	-image[xyz+Y2]
							-image[xyz+X-Y-Z]	-image[xyz-Y2-Z2]	-image[xyz+X2]
							-image[xyz+X+Y-Z]	-image[xyz-Z2]		-image[xyz+X2+Y2]
							-image[xyz-X+Y+Z]	-image[xyz-X2]		-image[xyz+Y2+Z2]
							-image[xyz+X-Y+Z]	-image[xyz-Y2]		-image[xyz+X2+Z2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[xyz]			+image[xyz+X-Y-Z]	+image[xyz-X+Y+Z])
					 		-image[xyz+X-Y+Z]	-image[xyz+X2-Y2]	-image[xyz+Z2]
					 		-image[xyz+X+Y-Z]	-image[xyz+X2-Z2]	-image[xyz+Y2]
					 		-image[xyz-X-Y-Z]	-image[xyz-Y2-Z2]	-image[xyz-X2]
					 		-image[xyz-X+Y-Z]	-image[xyz-Z2]		-image[xyz-X2+Y2]
					 		-image[xyz-X-Y+Z]	-image[xyz-Y2]		-image[xyz-X2+Z2]
					 		-image[xyz+X+Y+Z]	-image[xyz+X2]		-image[xyz+Y2+Z2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[xyz]			+image[xyz-X+Y-Z]	+image[xyz+X-Y+Z])
							 -image[xyz-X+Y+Z]	-image[xyz-X2+Y2]	-image[xyz+Z2]
							 -image[xyz+X+Y-Z]	-image[xyz+Y2-Z2]	-image[xyz+X2]
							 -image[xyz-X-Y-Z]	-image[xyz-X2-Z2]	-image[xyz-Y2]
							 -image[xyz+X-Y-Z]	-image[xyz-Z2]		-image[xyz+X2-Y2]
							 -image[xyz-X-Y+Z]	-image[xyz-X2]		-image[xyz-Y2+Z2]
							 -image[xyz+X+Y+Z]	-image[xyz+Y2]		-image[xyz+X2+Z2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[xyz]			+image[xyz-X-Y+Z]	+image[xyz+X+Y-Z])
					 		-image[xyz-X+Y+Z]	-image[xyz-X2+Z2]	-image[xyz+Y2]
					 		-image[xyz+X-Y+Z]	-image[xyz-Y2+Z2]	-image[xyz+X2]
					 		-image[xyz-X-Y-Z]	-image[xyz-X2-Y2]	-image[xyz-Z2]
					 		-image[xyz+X-Y-Z]	-image[xyz-Y2]		-image[xyz+X2-Z2]
					 		-image[xyz-X+Y-Z]	-image[xyz-X2]		-image[xyz+Y2-Z2]
					 		-image[xyz+X+Y+Z]	-image[xyz+Z2]		-image[xyz+X2+Y2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				//result[xyz] = 0.5f*best*(1.0f - Numerics.abs(second/best));
				result[xyz] = 0.5f*best;
			}
		}
		return result;
	}
	public static float[][][] tubularEdgeFilter(float[][][] image, boolean[][][] mask, int nx, int ny, int nz) {
		float[][][] result = new float[nx][ny][nz];
		
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x+nx*y+Z*z;
			result[x][y][z] = 0.0f;
			
			boolean process = true;
			for (int i=-2;i<=2 && process;i++) for (int j=-2;j<=2 && process;j++) for (int l=-2;l<=2 && process;l++) {
				process = mask[x+i][y+j][z+l];
			}
			if (process) {
				float val, best = 0.0f, second = 0.0f;
				
				val =(8.0f*(image[x][y][z]			+image[x-1][y][z]		+image[x+1][y][z])
					 		-image[x][y-1][z]		-image[x-1][y-1][z]		-image[x+1][y-1][z]
					 		-image[x][y+1][z]		-image[x-1][y+1][z]		-image[x+1][y+1][z]
					 		-image[x][y][z-1]		-image[x-1][y][z-1]		-image[x+1][y][z-1]
					 		-image[x][y][z+1]		-image[x-1][y][z+1]		-image[x+1][y][z+1]
					 		-image[x][y-1][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z-1]
					 		-image[x][y-1][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z+1]
					 		-image[x][y+1][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z-1]
					 		-image[x][y+1][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z+1])/24.0f;			
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				val =(8.0f*(image[x][y][z]			+image[x][y-1][z]		+image[x][y+1][z])
					 		-image[x-1][y][z]		-image[x-1][y-1][z]		-image[x-1][y+1][z]
					 		-image[x+1][y][z]		-image[x+1][y-1][z]		-image[x+1][y+1][z]
					 		-image[x][y][z-1]		-image[x][y-1][z-1]		-image[x][y+1][z-1]
					 		-image[x][y][z+1]		-image[x][y-1][z+1]		-image[x][y+1][z+1]
					 		-image[x-1][y][z-1]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z-1]
							-image[x-1][y][z+1]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z+1]
							-image[x+1][y][z-1]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z-1]
							-image[x+1][y][z+1]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z+1])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x][y][z-1]		+image[x][y][z+1])
					 		-image[x-1][y][z]		-image[x-1][y][z-1]		-image[x-1][y][z+1]
					 		-image[x+1][y][z]		-image[x+1][y][z-1]		-image[x+1][y][z+1]
					 		-image[x][y-1][z]		-image[x][y-1][z-1]		-image[x][y-1][z+1]
					 		-image[x][y+1][z]		-image[x][y+1][z-1]		-image[x][y+1][z+1]
					 		-image[x-1][y-1][z]		-image[x-1][y-1][z-1]	-image[x-1][y-1][z+1]
					 		-image[x-1][y+1][z]		-image[x-1][y+1][z-1]	-image[x-1][y+1][z+1]
					 		-image[x+1][y-1][z]		-image[x+1][y-1][z-1]	-image[x+1][y-1][z+1]
					 		-image[x+1][y+1][z]		-image[x+1][y+1][z-1]	-image[x+1][y+1][z+1])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x-1][y-1][z]		+image[x+1][y+1][z])
					 		-image[x-1][y+1][z]		-image[x-2][y][z]		-image[x][y+2][z]
					 		-image[x+1][y-1][z]		-image[x][y-2][z]		-image[x+2][y][z]
					 		-image[x][y][z-1]		-image[x-1][y-1][z-1]	-image[x+1][y+1][z-1]
					 		-image[x-1][y+1][z-1]	-image[x-2][y][z-1]		-image[x][y+2][z-1]
					 		-image[x+1][y-1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
					 		-image[x][y][z+1]		-image[x-1][y-1][z+1]	-image[x+1][y+1][z+1]
					 		-image[x-1][y+1][z+1]	-image[x-2][y][z+1]		-image[x][y+2][z+1]
					 		-image[x+1][y-1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x][y-1][z-1]		+image[x][y+1][z+1])
					 		-image[x][y+1][z-1]		-image[x][y][z-2]		-image[x][y+2][z]
					 		-image[x][y-1][z+1]		-image[x][y-2][z]		-image[x][y][z+2]
					 		-image[x-1][y][z]		-image[x-1][y-1][z-1]	-image[x-1][y+1][z+1]
					 		-image[x-1][y+1][z-1]	-image[x-1][y][z-2]		-image[x-1][y-2][z]
					 		-image[x-1][y-1][z+1]	-image[x-1][y-2][z]		-image[x-1][y][z+2]
					 		-image[x+1][y][z]		-image[x+1][y-1][z-1]	-image[x+1][y+1][z+1]
					 		-image[x+1][y+1][z-1]	-image[x+1][y][z-2]		-image[x+1][y][z+2]
					 		-image[x+1][y-1][z+1]	-image[x+1][y-2][z]		-image[x+1][y][z+2])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x-1][y][z-1]		+image[x+1][y][z+1])
					 		-image[x+1][y][z-1]		-image[x][y][z-2]		-image[x+2][y][z]
					 		-image[x-1][y][z+1]		-image[x-2][y][z]		-image[x][y][z+2]
					 		-image[x][y-1][z]		-image[x-1][y-1][z-1]	-image[x+1][y-1][z+1]
					 		-image[x+1][y-1][z-1]	-image[x][y-1][z-2]		-image[x+2][y-1][z]
					 		-image[x-1][y-1][z+1]	-image[x-2][y-1][z]		-image[x][y-1][z+2]
					 		-image[x][y+1][z]		-image[x-1][y+1][z-1]	-image[x+1][y+1][z+1]
					 		-image[x+1][y+1][z-1]	-image[x][y+1][z-2]		-image[x+2][y+1][z]
					 		-image[x-1][y+1][z+1]	-image[x-2][y+1][z]		-image[x][y+1][z+2])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x-1][y+1][z]		+image[x+1][y-1][z])
					 		-image[x-1][y-1][z]		-image[x-2][y][z]		-image[x][y-2][z]
					 		-image[x+1][y+1][z]		-image[x][y-2][z]		-image[x+2][y][z]
					 		-image[x][y][z-1]		-image[x-1][y+1][z-1]	-image[x+1][y-1][z-1]
					 		-image[x-1][y-1][z-1]	-image[x-2][y][z-1]		-image[x][y-2][z-1]
					 		-image[x+1][y+1][z-1]	-image[x][y-2][z-1]		-image[x+2][y][z-1]
					 		-image[x][y][z+1]		-image[x-1][y+1][z+1]	-image[x+1][y-1][z+1]
					 		-image[x-1][y-1][z+1]	-image[x-2][y][z+1]		-image[x][y-2][z+1]
					 		-image[x+1][y+1][z+1]	-image[x][y-2][z+1]		-image[x+2][y][z+1])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x][y-1][z+1]		+image[x][y+1][z-1])
					 		-image[x][y-1][z-1]		-image[x][y-2][z]		-image[x][y][z-2]
					 		-image[x][y+1][z+1]		-image[x][y][z+2]		-image[x][y+2][z]
					 		-image[x-1][y][z]		-image[x-1][y-1][z+1]	-image[x-1][y+1][z-1]
					 		-image[x-1][y-1][z-1]	-image[x-1][y-2][z]		-image[x-1][y][z-2]
					 		-image[x-1][y+1][z+1]	-image[x-1][y][z+2]		-image[x-1][y+2][z]
					 		-image[x+1][y][z]		-image[x+1][y-1][z+1]	-image[x+1][y+1][z-1]
					 		-image[x+1][y-1][z-1]	-image[x+1][y-2][z]		-image[x+1][y][z-2]
					 		-image[x+1][y+1][z+1]	-image[x+1][y][z+2]		-image[x+1][y+2][z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(8.0f*(image[x][y][z]			+image[x-1][y][z+1]		+image[x+1][y][z-1])
					 		-image[x-1][y][z-1]		-image[x-2][y][z]		-image[x][y][z-2]
					 		-image[x+1][y][z+1]		-image[x][y][z+2]		-image[x+2][y][z]
					 		-image[x][y-1][z]		-image[x-1][y-1][z+1]	-image[x+1][y-1][z-1]
					 		-image[x-1][y-1][z-1]	-image[x-2][y-1][z]		-image[x][y-1][z-2]
					 		-image[x+1][y-1][z+1]	-image[x][y-1][z+2]		-image[x+2][y-1][z]
					 		-image[x][y+1][z]		-image[x-1][y+1][z+1]	-image[x+1][y+1][z-1]
					 		-image[x-1][y+1][z-1]	-image[x-2][y+1][z]		-image[x][y+1][z-2]
					 		-image[x+1][y+1][z+1]	-image[x][y+1][z+2]		-image[x+2][y+1][z])/24.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z-1]	+image[x+1][y+1][z+1])
							-image[x-1][y-1][z+1]	-image[x-2][y-2][z]		-image[x][y][z+2]
							-image[x-1][y+1][z-1]	-image[x-2][y][z-2]		-image[x][y+2][z]
							-image[x+1][y-1][z-1]	-image[x][y-2][z-2]		-image[x+2][y][z]
							-image[x+1][y+1][z-1]	-image[x][y][z-2]		-image[x+2][y+2][z]
							-image[x-1][y+1][z+1]	-image[x-2][y][z]		-image[x][y+2][z+2]
							-image[x+1][y-1][z+1]	-image[x][y-2][z]		-image[x+2][y][z+2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[x][y][z]			+image[x+1][y-1][z-1]	+image[x-1][y+1][z+1])
					 		-image[x+1][y-1][z+1]	-image[x+2][y-2][z]		-image[x][y][z+2]
					 		-image[x+1][y+1][z-1]	-image[x+2][y][z-2]		-image[x][y+2][z]
					 		-image[x-1][y-1][z-1]	-image[x][y-2][z-2]		-image[x-2][y][z]
					 		-image[x-1][y+1][z-1]	-image[x][y][z-2]		-image[x-2][y+2][z]
					 		-image[x-1][y-1][z+1]	-image[x][y-2][z]		-image[x-2][y][z+2]
					 		-image[x+1][y+1][z+1]	-image[x+2][y][z]		-image[x][y+2][z+2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[x][y][z]			+image[x-1][y+1][z-1]	+image[x+1][y-1][z+1])
							 -image[x-1][y+1][z+1]	-image[x-2][y+2][z]		-image[x][y][z+2]
							 -image[x+1][y+1][z-1]	-image[x][y+2][z-2]		-image[x+2][y][z]
							 -image[x-1][y-1][z-1]	-image[x-2][y][z-2]		-image[x][y-2][z]
							 -image[x+1][y-1][z-1]	-image[x][y][z-2]		-image[x+2][y-2][z]
							 -image[x-1][y-1][z+1]	-image[x-2][y][z]		-image[x][y-2][z+2]
							 -image[x+1][y+1][z+1]	-image[x][y+2][z]		-image[x+2][y][z+2])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(6.0f*(image[x][y][z]			+image[x-1][y-1][z+1]	+image[x+1][y+1][z-1])
					 		-image[x-1][y+1][z+1]	-image[x-2][y][z-2]		-image[x][y+2][z]
					 		-image[x+1][y-1][z+1]	-image[x][y-2][z+2]		-image[x+2][y][z]
					 		-image[x-1][y-1][z-1]	-image[x-2][y-2][z]		-image[x][y][z-2]
					 		-image[x+1][y-1][z-1]	-image[x][y-2][z]		-image[x+2][y][z-2]
					 		-image[x-1][y+1][z-1]	-image[x-2][y][z]		-image[x][y+2][z-2]
					 		-image[x+1][y+1][z+1]	-image[x][y][z+2]		-image[x+2][y+2][z])/18.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				//result[x][y][z] = 0.5f*best*(1.0f - Numerics.abs(second/best));
				result[x][y][z] = 0.5f*best;
			}
		}
		return result;
	}
		
	public static float[] sphereEdgeFilter(float[] image, boolean[] mask, int nx, int ny, int nz) {
		float[] result = new float[nx*ny*nz];
		
		int X = 1;
		int Y = nx;
		int Z = nx*ny;
		int X2 = 2;
		int Y2 = 2*nx;
		int Z2 = 2*nx*ny;
		for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x+nx*y+Z*z;
			result[xyz] = 0.0f;
			
			boolean process = true;
			for (int i=-2;i<=2 && process;i++) for (int j=-2;j<=2 && process;j++) for (int l=-2;l<=2 && process;l++) {
				process = mask[xyz+i+j*nx+l*nx*ny];
			}
			if (process) {
				float val, best = 0.0f, second = 0.0f;
				
				val =(26.0f*(image[xyz])
							-image[xyz-X]		-image[xyz+X]
					 		-image[xyz-Y]		-image[xyz-X-Y]		-image[xyz+X-Y]
					 		-image[xyz+Y]		-image[xyz-X+Y]		-image[xyz+X+Y]
					 		-image[xyz-Z]		-image[xyz-Z-X]		-image[xyz-Z+X]
					 		-image[xyz+Z]		-image[xyz+Z-X]		-image[xyz+Z+X]
					 		-image[xyz-Y-Z]		-image[xyz-X-Y-Z]	-image[xyz+X-Y-Z]
					 		-image[xyz-Y+Z]		-image[xyz-X-Y+Z]	-image[xyz+X-Y+Z]
					 		-image[xyz+Y-Z]		-image[xyz-X+Y-Z]	-image[xyz+X+Y-Z]
					 		-image[xyz+Y+Z]		-image[xyz-X+Y+Z]	-image[xyz+X+Y+Z])/26.0f;			
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(26.0f*(image[xyz])
							-image[xyz-X-Y]		-image[xyz+X+Y]
					 		-image[xyz-X+Y]		-image[xyz-X2]		-image[xyz+Y2]
					 		-image[xyz+X-Y]		-image[xyz-Y2]		-image[xyz+X2]
					 		-image[xyz-Z]		-image[xyz-X-Y-Z]	-image[xyz+X+Y-Z]
					 		-image[xyz-X+Y-Z]	-image[xyz-X2-Z]	-image[xyz+Y2-Z]
					 		-image[xyz+X-Y-Z]	-image[xyz-Y2-Z]	-image[xyz+X2-Z]
					 		-image[xyz+Z]		-image[xyz-X-Y+Z]	-image[xyz+X+Y+Z]
					 		-image[xyz-X+Y+Z]	-image[xyz-X2+Z]	-image[xyz+Y2+Z]
					 		-image[xyz+X-Y+Z]	-image[xyz-Y2+Z]	-image[xyz+X2+Z])/26.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				val =(20.0f*(image[xyz])
							-image[xyz-X-Y-Z]	-image[xyz+X+Y+Z]
							-image[xyz-X-Y+Z]	-image[xyz-X2-Y2]	-image[xyz+Z2]
							-image[xyz-X+Y-Z]	-image[xyz-X2-Z2]	-image[xyz+Y2]
							-image[xyz+X-Y-Z]	-image[xyz-Y2-Z2]	-image[xyz+X2]
							-image[xyz+X+Y-Z]	-image[xyz-Z2]		-image[xyz+X2+Y2]
							-image[xyz-X+Y+Z]	-image[xyz-X2]		-image[xyz+Y2+Z2]
							-image[xyz+X-Y+Z]	-image[xyz-Y2]		-image[xyz+X2+Z2])/20.0f;
				if (val*val>best*best) {
					second = best;
					best = val;
				} else if (val*val>second*second) {
					second = val;
				}
				
				//result[xyz] = 0.5f*best*(1.0f - Numerics.abs(second/best));
				result[xyz] = 0.5f*best;
			}
		}
		return result;
	}
	
	public static float[][] separableRandomKernel(float sx, float sy, float sz) {
		int kx,ky,kz;
		float sum;
		
		// kernel size
		kx = Numerics.ceil(Math.max(3.0f*sx-0.5f,0.0f));
		ky = Numerics.ceil(Math.max(3.0f*sy-0.5f,0.0f));
		kz = Numerics.ceil(Math.max(3.0f*sz-0.5f,0.0f));
		
		//MedicUtilPublic.displayMessage("kernel size: "+kx+"x"+ky+"x"+kz+"\n");
		//MedicUtilPublic.displayMessage("scale: "+sx+"x"+sy+"x"+sz+"\n");
		// create the kernel
		float[][] kernel = new float[3][];
		kernel[0] = new float[2*kx+1]; 
		kernel[1] = new float[2*ky+1]; 
		kernel[2] = new float[2*kz+1]; 
		
		sum = 0.0f;
		for (int i=-kx;i<=kx;i++) {
			kernel[0][kx+i] = (float)FastMath.random();
			sum += kernel[0][kx+i];
		}
		for (int i=-kx;i<=kx;i++) kernel[0][kx+i] /= sum;
		
		sum = 0.0f;
		for (int j=-ky;j<=ky;j++) {
			kernel[1][ky+j] = (float)FastMath.random();
			sum += kernel[1][ky+j];
		}
		for (int j=-ky;j<=ky;j++) kernel[1][ky+j] /= sum;
		
		sum = 0.0f;
		for (int l=-kz;l<=kz;l++) {
			kernel[2][kz+l] = (float)FastMath.random();
			sum += kernel[2][kz+l];
		}
		for (int l=-kz;l<=kz;l++) kernel[2][kz+l] /= sum;
		
		return kernel;
	}
	
}
