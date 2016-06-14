package de.mpg.cbs.libraries;

import Jama.Matrix;
import Jama.EigenvalueDecomposition;
import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;



/**
 *
 *  This class computes various basic image geometric features.
 *	
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ImageGeometry {
	
	// no data: used as a library of functions
	
    // simple image functions
    
    public static final byte X = 0;
    public static final byte Y = 1;
    public static final byte Z = 2;
    public static final byte T = 3;

	/**
	 *	Horn-like derivatives: use a 2x2 forward square as averaging
	 *	the gradient images are passed from above
	 */
    public static void computeGradientHorn(float[][][] img, float[][][] dx, float[][][] dy, float[][][] dz, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=0;x<nx-1;x++) for (y=0;y<ny-1;y++) for (z=0;z<nz-1;z++) {
			dx[x][y][z] = ( (img[x+1][y][z] + img[x+1][y+1][z] + img[x+1][y][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x][y+1][z] + img[x][y][z+1] + img[x][y+1][z+1]) )/4.0f;
			
			dy[x][y][z] = ( (img[x][y+1][z] + img[x+1][y+1][z] + img[x][y+1][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x+1][y][z] + img[x][y][z+1] + img[x+1][y][z+1]) )/4.0f;
			
			dz[x][y][z] = ( (img[x][y][z+1] + img[x+1][y][z+1] + img[x][y+1][z+1] + img[x+1][y+1][z+1])
						  -(img[x][y][z] + img[x+1][y][z] + img[x][y+1][z] + img[x+1][y+1][z]) )/4.0f;
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
			dx[x][y][nz-1] = 0.0f;
			dy[x][y][nz-1] = 0.0f;
			dz[x][y][nz-1] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
			dx[x][ny-1][z] = 0.0f;
			dy[x][ny-1][z] = 0.0f;
			dz[x][ny-1][z] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			dx[nx-1][y][z] = 0.0f;
			dy[nx-1][y][z] = 0.0f;
			dz[nx-1][y][z] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple gradient (central differences)
	 */
    public static void computeGradient(float[] img, float[][] J, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			int xyz = x+nx*y+nx*ny*z;
			J[0][xyz] = 0.5f*(img[xyz+1] - img[xyz-1]);
			J[1][xyz] = 0.5f*(img[xyz+nx] - img[xyz-nx]);
			J[2][xyz] = 0.5f*(img[xyz+nx*ny] - img[xyz-nx*ny]);
		}
		return;
    }
	
	/**
	 *	simple gradient (central differences)
	 */
    public static void computeGradient(float[][][] img, float[][][] dx, float[][][] dy, float[][][] dz, int nx, int ny, int nz) {
        int x,y,z;
               
        for (x=0;x<nx-1;x++) for (y=0;y<ny-1;y++) for (z=0;z<nz-1;z++) {
			dx[x][y][z] = 0.5f*(img[x+1][y][z] - img[x-1][y][z]);
			dy[x][y][z] = 0.5f*(img[x][y+1][z] - img[x][y-1][z]);
			dz[x][y][z] = 0.5f*(img[x][y][z+1] - img[x][y][z-1]);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
			dx[x][y][nz-1] = 0.0f;
			dy[x][y][nz-1] = 0.0f;
			dz[x][y][nz-1] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) {
			dx[x][ny-1][z] = 0.0f;
			dy[x][ny-1][z] = 0.0f;
			dz[x][ny-1][z] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			dx[nx-1][y][z] = 0.0f;
			dy[nx-1][y][z] = 0.0f;
			dz[nx-1][y][z] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessian(float[][][] img, float[][][][] hessian, int nx, int ny, int nz) {
        int x,y,z,n;
               
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			hessian[x][y][z][0] = img[x+1][y][z] + img[x-1][y][z] - 2.0f*img[x][y][z];
			hessian[x][y][z][1] = img[x][y+1][z] + img[x][y-1][z] - 2.0f*img[x][y][z];
			hessian[x][y][z][2] = img[x][y][z+1] + img[x][y][z-1] - 2.0f*img[x][y][z];
			
			hessian[x][y][z][3] = img[x+1][y+1][z] + img[x-1][y-1][z] - img[x+1][y-1][z] - img[x-1][y+1][z];
			hessian[x][y][z][4] = img[x][y+1][z+1] + img[x][y-1][z-1] - img[x][y+1][z-1] - img[x][y-1][z+1];
			hessian[x][y][z][5] = img[x+1][y][z+1] + img[x-1][y][z-1] - img[x-1][y][z+1] - img[x+1][y][z-1];
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<6;n++) {
			hessian[x][y][0][n] = 0.0f;
			hessian[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<6;n++) {
			hessian[x][0][z][n] = 0.0f;
			hessian[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<6;n++) {
			hessian[0][y][z][n] = 0.0f;
			hessian[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessianAt(float[] img, int xyz, double[][] hessian, int nx, int ny, int nz, float rx, float ry, float rz) {
       
    	hessian[0][0] = (img[xyz-1]			+img[xyz+1]			-2.0*img[xyz])/rx/rx;
		hessian[1][1] = (img[xyz-nx]		+img[xyz+nx]		-2.0*img[xyz])/ry/ry;
		hessian[2][2] = (img[xyz-nx*ny]	+img[xyz+nx*ny]		-2.0*img[xyz])/rz/rz;
		
		hessian[0][1] = (img[xyz-1-nx]		+img[xyz+1+nx]		-img[xyz+1-nx]		-img[xyz-1+nx])/4.0/rx/ry;
		hessian[1][2] = (img[xyz-nx-nx*ny]	+img[xyz+nx+nx*ny]	-img[xyz+nx-nx*ny]	-img[xyz-nx+nx*ny])/4.0/ry/rz;
		hessian[2][0] = (img[xyz-nx*ny-1]	+img[xyz+nx*ny+1]	-img[xyz+nx*ny-1]	-img[xyz-nx*ny+1])/4.0/rz/rx;
		
		hessian[1][0] = hessian[0][1];
		hessian[2][1] = hessian[1][2];
		hessian[0][2] = hessian[2][0];					
					
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessianAt(float[] img, int xyz, double[][] hessian, int nx, int ny, int nz) {
       
    	hessian[0][0] = (img[xyz-1]			+img[xyz+1]			-2.0*img[xyz]);
		hessian[1][1] = (img[xyz-nx]		+img[xyz+nx]		-2.0*img[xyz]);
		hessian[2][2] = (img[xyz-nx*ny]	+img[xyz+nx*ny]		-2.0*img[xyz]);
		
		hessian[0][1] = (img[xyz-1-nx]		+img[xyz+1+nx]		-img[xyz+1-nx]		-img[xyz-1+nx])/4.0;
		hessian[1][2] = (img[xyz-nx-nx*ny]	+img[xyz+nx+nx*ny]	-img[xyz+nx-nx*ny]	-img[xyz-nx+nx*ny])/4.0;
		hessian[2][0] = (img[xyz-nx*ny-1]	+img[xyz+nx*ny+1]	-img[xyz+nx*ny-1]	-img[xyz-nx*ny+1])/4.0;
		
		hessian[1][0] = hessian[0][1];
		hessian[2][1] = hessian[1][2];
		hessian[0][2] = hessian[2][0];					
					
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessianOrder2At(float[] img, int xyz, double[][] hessian, int nx, int ny, int nz) {
       
    	hessian[0][0] = (-1.0/12.0*img[xyz-2]			+4.0/3.0*img[xyz-1]			-5.0/2.0*img[xyz]	+4.0/3.0*img[xyz+1]			-1.0/12.0*img[xyz+2]);
		hessian[1][1] = (-1.0/12.0*img[xyz-2*nx]		+4.0/3.0*img[xyz-nx]		-5.0/2.0*img[xyz]	+4.0/3.0*img[xyz+nx]		-1.0/12.0*img[xyz+2*nx]);
		hessian[2][2] = (-1.0/12.0*img[xyz-2*nx*ny]	+4.0/3.0*img[xyz-nx*ny]		-5.0/2.0*img[xyz]	+4.0/3.0*img[xyz+nx*ny]		-1.0/12.0*img[xyz+2*nx*ny]);
		
		hessian[0][1] = (1.0/144.0*img[xyz-2-2*nx]		-1.0/18.0*img[xyz-2-nx] 	+1.0/18.0*img[xyz-2+nx]		-1.0/144.0*img[xyz-2+2*nx]
						  -1.0/18.0*img[xyz-1-2*nx]		 +4.0/9.0*img[xyz-1-nx]		 -4.0/9.0*img[xyz-1+nx]		 +1.0/18.0*img[xyz-1+2*nx]
						  +1.0/18.0*img[xyz+1-2*nx]		 -4.0/9.0*img[xyz+1-nx]		 +4.0/9.0*img[xyz+1+nx]		 -1.0/18.0*img[xyz+1+2*nx]
						 -1.0/144.0*img[xyz+2-2*nx]		+1.0/18.0*img[xyz+2-nx] 	-1.0/18.0*img[xyz+2+nx]		+1.0/144.0*img[xyz+2+2*nx]);
		
		hessian[1][2] = (1.0/144.0*img[xyz-2*nx-2*nx*ny]		-1.0/18.0*img[xyz-2*nx-nx*ny] 	+1.0/18.0*img[xyz-2*nx+nx*ny]		-1.0/144.0*img[xyz-2*nx+2*nx*ny]
						  -1.0/18.0*img[xyz-1*nx-2*nx*ny]		 +4.0/9.0*img[xyz-1*nx-nx*ny]	 -4.0/9.0*img[xyz-1*nx+nx*ny]		 +1.0/18.0*img[xyz-1*nx+2*nx*ny]
						  +1.0/18.0*img[xyz+1*nx-2*nx*ny]		 -4.0/9.0*img[xyz+1*nx-nx*ny]	 +4.0/9.0*img[xyz+1*nx+nx*ny]		 -1.0/18.0*img[xyz+1*nx+2*nx*ny]
						 -1.0/144.0*img[xyz+2*nx-2*nx*ny]		+1.0/18.0*img[xyz+2*nx-nx*ny] 	-1.0/18.0*img[xyz+2*nx+nx*ny]		+1.0/144.0*img[xyz+2*nx+2*nx*ny]);
		
		hessian[2][0] = (1.0/144.0*img[xyz-2-2*nx*ny]		-1.0/18.0*img[xyz-2-nx*ny] 	+1.0/18.0*img[xyz-2+nx*ny]		-1.0/144.0*img[xyz-2+2*nx*ny]
						  -1.0/18.0*img[xyz-1-2*nx*ny]		 +4.0/9.0*img[xyz-1-nx*ny]	 -4.0/9.0*img[xyz-1+nx*ny]		 +1.0/18.0*img[xyz-1+2*nx*ny]
						  +1.0/18.0*img[xyz+1-2*nx*ny]		 -4.0/9.0*img[xyz+1-nx*ny]	 +4.0/9.0*img[xyz+1+nx*ny]		 -1.0/18.0*img[xyz+1+2*nx*ny]
						 -1.0/144.0*img[xyz+2-2*nx*ny]		+1.0/18.0*img[xyz+2-nx*ny] 	-1.0/18.0*img[xyz+2+nx*ny]		+1.0/144.0*img[xyz+2+2*nx*ny]);
		
		hessian[1][0] = hessian[0][1];
		hessian[2][1] = hessian[1][2];
		hessian[0][2] = hessian[2][0];					
					
		return;
    }
	
	/**
	 *	simple Hessian
	 */
    public static void computeHessianEigenvalues(float[][][] img, float[][][][] values, int nx, int ny, int nz) {
        int x,y,z,n;
		float[][] hessian = new float[3][3];
		       
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			hessian[0][0] = img[x+1][y][z] + img[x-1][y][z] - 2.0f*img[x][y][z];
			hessian[1][1] = img[x][y+1][z] + img[x][y-1][z] - 2.0f*img[x][y][z];
			hessian[2][2] = img[x][y][z+1] + img[x][y][z-1] - 2.0f*img[x][y][z];
			
			hessian[0][1] = img[x+1][y+1][z] + img[x-1][y-1][z] - img[x+1][y-1][z] - img[x-1][y+1][z];
			hessian[1][2] = img[x][y+1][z+1] + img[x][y-1][z-1] - img[x][y+1][z-1] - img[x][y-1][z+1];
			hessian[2][0] = img[x+1][y][z+1] + img[x-1][y][z-1] - img[x-1][y][z+1] - img[x+1][y][z-1];
			
			hessian[1][0] = hessian[0][1];
			hessian[2][1] = hessian[1][2];
			hessian[0][2] = hessian[2][0];
			
			// compute eigenvalues
			values[x][y][z] = Matrix3D.eigenvalues(hessian);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<3;n++) {
			values[x][y][0][n] = 0.0f;
			values[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[x][0][z][n] = 0.0f;
			values[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[0][y][z][n] = 0.0f;
			values[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
	 *	gradient shape tensor: return the eigenvalues at each point
	 */
    public static void computeGSTeigenvalues(float[][][] img, float[][][][] values, int nx, int ny, int nz) {
        int x,y,z,n;
		float[][] gst = new float[3][3];
		       
        for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			// IxIx
			gst[0][0] = 0.25f*(img[x+1][y][z]-img[x-1][y][z])*(img[x+1][y][z]-img[x-1][y][z]);
			// IyIy
			gst[1][1] = 0.25f*(img[x][y+1][z]-img[x][y-1][z])*(img[x][y+1][z]-img[x][y-1][z]);
			// IzIz
			gst[2][2] = 0.25f*(img[x][y][z+1]-img[x][y][z-1])*(img[x][y][z+1]-img[x][y][z-1]);
			// IxIy
			gst[0][1] = 0.25f*(img[x+1][y][z]-img[x-1][y][z])*(img[x][y+1][z]-img[x][y-1][z]);
			gst[1][0] = gst[0][1];
			// IyIz
			gst[1][2] = 0.25f*(img[x][y+1][z]-img[x][y-1][z])*(img[x][y][z+1]-img[x][y][z-1]);
			gst[2][1] = gst[1][2];
			// Izx ?
			gst[2][0] = 0.25f*(img[x][y][z+1]-img[x][y][z-1])*(img[x+1][y][z]-img[x-1][y][z]);
			gst[0][2] = gst[2][0];
			
			// compute eigenvalues
			values[x][y][z] = Matrix3D.eigenvalues(gst);
		}
		// boundaries: zero
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (n=0;n<3;n++) {
			values[x][y][0][n] = 0.0f;
			values[x][y][nz-1][n] = 0.0f;
		}
		for (x=0;x<nx;x++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[x][0][z][n] = 0.0f;
			values[x][ny-1][z][n] = 0.0f;
		}
		for (y=0;y<ny;y++) for (z=0;z<nz;z++) for (n=0;n<3;n++) {
			values[0][y][z][n] = 0.0f;
			values[nx-1][y][z][n] = 0.0f;
		}
		return;
    }
	
	/**
     *    GVF-style diffusion on scalar images
     *    @param 	image	the original image, assumed to be in [0,1]
	 *    @param 	ratio	amount of diffusion
	 *    @param 	steps	number of diffusion steps
	 *    @param 	dt		time interval
	 *	  @return 			the diffused image
     */
    public static final float[][][] heatDiffusion(float[][][] image, float ratio, int steps, float dt, int nx, int ny, int nz ) {
		float[][][] diff = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			diff[x][y][z] = image[x][y][z];
		}
		for (int t=0;t<steps;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				diff[x][y][z] += ratio*dt/6.0f*(diff[x-1][y][z]+diff[x+1][y][z]
											   +diff[x][y-1][z]+diff[x][y+1][z]
											   +diff[x][y][z-1]+diff[x][y][z+1]
											   -6.0f*diff[x][y][z])
								 + image[x][y][z]*image[x][y][z]*dt*(image[x][y][z]-diff[x][y][z]);
			}
		}
		return diff;
	}
	/**
     *    GVF-style diffusion on scalar images
     *    @param 	image	the original image
	 *	  @param	scale	maximum image value
	 *    @param 	ratio	amount of diffusion
	 *    @param 	steps	number of diffusion steps
	 *    @param 	dt		time interval
	 *	  @return 			the diffused image
     */
    public static final float[][][] heatDiffusion(float[][][] image, float scale, float ratio, int steps, float dt, int nx, int ny, int nz ) {
		float[][][] diff = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			diff[x][y][z] = image[x][y][z];
		}
		for (int t=0;t<steps;t++) {
			for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
				diff[x][y][z] += ratio*dt/6.0f*(diff[x-1][y][z]+diff[x+1][y][z]
											   +diff[x][y-1][z]+diff[x][y+1][z]
											   +diff[x][y][z-1]+diff[x][y][z+1]
											   -6.0f*diff[x][y][z])
								 + (image[x][y][z]*image[x][y][z])/(scale*scale)*dt*(image[x][y][z]-diff[x][y][z]);
			}
		}
		return diff;
	}
		
	/**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][][] meanCurvature(float[][][] image, int nx, int ny, int nz ) {
		float[][][] curv = new float[nx][ny][nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// central differences, no spatial smoothing
			ux = 0.5*(image[x+1][y][z]-image[x-1][y][z]);
			uy = 0.5*(image[x][y+1][z]-image[x][y-1][z]);
			uz = 0.5*(image[x][y][z+1]-image[x][y][z-1]);
				
			uxx = image[x+1][y][z] - 2.0*image[x][y][z] + image[x-1][y][z];
			uyy = image[x][y+1][z] - 2.0*image[x][y][z] + image[x][y-1][z];
			uzz = image[x][y][z+1] - 2.0*image[x][y][z] + image[x][y][z-1];
			
			uxy = 0.25*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			uyz = 0.25*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			uzx = 0.25*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				curv[x][y][z] = (float)(num/den);
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}
    
	/**
     *    squared curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][][] squaredMeanCurvature(float[][][] image, int nx, int ny, int nz ) {
		float[][][] sqcurv = new float[nx][ny][nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// central differences, no spatial smoothing
			ux = 0.5*(image[x+1][y][z]-image[x-1][y][z]);
			uy = 0.5*(image[x][y+1][z]-image[x][y-1][z]);
			uz = 0.5*(image[x][y][z+1]-image[x][y][z-1]);
				
			uxx = image[x+1][y][z] - 2.0*image[x][y][z] + image[x-1][y][z];
			uyy = image[x][y+1][z] - 2.0*image[x][y][z] + image[x][y-1][z];
			uzz = image[x][y][z+1] - 2.0*image[x][y][z] + image[x][y][z-1];
			
			uxy = 0.25*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			uyz = 0.25*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			uzx = 0.25*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				sqcurv[x][y][z] = (float)(num/den)*(float)(num/den);
			} else {
				sqcurv[x][y][z] = 0.0f;
			}
		}
		return sqcurv;
	}
    
	
	/**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[] meanCurvature(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[] curv = new float[nx*ny*nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			// central differences, no spatial smoothing
			ux = 0.5*(image[xyzpx]-image[xyzmx]);
			uy = 0.5*(image[xyzpy]-image[xyzmy]);
			uz = 0.5*(image[xyzpz]-image[xyzmz]);
				
			uxx = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			uyy = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			uzz = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			uxy = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			uyz = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			uzx = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				curv[xyz] = (float)(num/den);
			} else {
				curv[xyz] = 0.0f;
			}
		}
		return curv;
	}
    
    /**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[] squaredMeanCurvature(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[] sqcurv = new float[nx*ny*nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			// central differences, no spatial smoothing
			ux = 0.5*(image[xyzpx]-image[xyzmx]);
			uy = 0.5*(image[xyzpy]-image[xyzmy]);
			uz = 0.5*(image[xyzpz]-image[xyzmz]);
				
			uxx = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			uyy = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			uzz = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			uxy = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			uyz = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			uzx = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				sqcurv[xyz] = (float)(num/den)*(float)(num/den);
			} else {
				sqcurv[xyz] = 0.0f;
			}
		}
		return sqcurv;
	}
    
	
    public static final double[] squaredMeanCurvature(float[] image, int nx, int ny, int nz, boolean[] mask) {
		double[] sqcurv = new double[nx*ny*nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			// central differences, no spatial smoothing
			ux = 0.5*(image[xyzpx]-image[xyzmx]);
			uy = 0.5*(image[xyzpy]-image[xyzmy]);
			uz = 0.5*(image[xyzpz]-image[xyzmz]);
				
			uxx = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			uyy = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			uzz = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			uxy = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			uyz = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			uzx = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				sqcurv[xyz] = (num/den)*(num/den);
			} else {
				sqcurv[xyz] = 0.0;
			}
		}
		return sqcurv;
	}
    
    /**
     *    squared curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[] squaredMeanCurvature(float[] image, int nx, int ny, int nz ) {
		float[] sqcurv = new float[nx*ny*nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			//if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			//if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		//if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		//if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		//if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		//if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		//if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		//if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		//if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		//if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		//if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		//if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		//if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		//if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		//if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		//if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		//if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		//if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			// central differences, no spatial smoothing
			ux = 0.5*(image[xyzpx]-image[xyzmx]);
			uy = 0.5*(image[xyzpy]-image[xyzmy]);
			uz = 0.5*(image[xyzpz]-image[xyzmz]);
				
			uxx = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			uyy = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			uzz = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			uxy = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			uyz = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			uzx = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
							 
			// 3D mean curvature
			num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
			den = FastMath.sqrt(ux*ux + uy*uy + uz*uz);
			den = 2.0*den*den*den;
				
			if (den>0) {
				sqcurv[xyz] = (float)(num/den)*(float)(num/den);
			} else {
				sqcurv[xyz] = 0.0f;
			}
		}
		return sqcurv;
	}
    
	/**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][][] gaussCurvature(float[][][] image, int nx, int ny, int nz ) {
		float[][][] curv = new float[nx][ny][nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// central differences, no spatial smoothing
			ux = 0.5*(image[x+1][y][z]-image[x-1][y][z]);
			uy = 0.5*(image[x][y+1][z]-image[x][y-1][z]);
			uz = 0.5*(image[x][y][z+1]-image[x][y][z-1]);
				
			uxx = image[x+1][y][z] - 2.0*image[x][y][z] + image[x-1][y][z];
			uyy = image[x][y+1][z] - 2.0*image[x][y][z] + image[x][y-1][z];
			uzz = image[x][y][z+1] - 2.0*image[x][y][z] + image[x][y][z-1];
			
			uxy = 0.25*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			uyz = 0.25*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			uzx = 0.25*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
							 
			// 3D gaussian curvature
			num = ux*ux*(uyy*uzz-uyz*uyz) + uy*uy*(uzz*uxx-uzx*uzx) + uz*uz*(uxx*uyy-uxy*uxy)
					+ 2.0*ux*uy*(uyz*uzx-uxy*uzz) + 2.0*uy*uz*(uxy*uzx-uxx*uyz) + 2.0*uz*ux*(uxy*uyz-uyy*uzx);
			den = ux*ux + uy*uy + uz*uz;
			den = den*den;
				
			if (den>0) {
				curv[x][y][z] = (float)(num/den);
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}

	/**
     *    curvature estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[] gaussCurvature(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[] curv = new float[nx*ny*nz];
		double ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		double num,den;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			// central differences, no spatial smoothing
			ux = 0.5*(image[xyzpx]-image[xyzmx]);
			uy = 0.5*(image[xyzpy]-image[xyzmy]);
			uz = 0.5*(image[xyzpz]-image[xyzmz]);
				
			uxx = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			uyy = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			uzz = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			uxy = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			uyz = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			uzx = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
							 
			// 3D gaussian curvature
			num = ux*ux*(uyy*uzz-uyz*uyz) + uy*uy*(uzz*uxx-uzx*uzx) + uz*uz*(uxx*uyy-uxy*uxy)
					+ 2.0*ux*uy*(uyz*uzx-uxy*uzz) + 2.0*uy*uz*(uxy*uzx-uxx*uyz) + 2.0*uz*ux*(uxy*uyz-uyy*uzx);
			den = ux*ux + uy*uy + uz*uz;
			den = den*den;
				
			if (den>0) {
				curv[xyz] = (float)(num/den);
			} else {
				curv[xyz] = 0.0f;
			}
		}
		return curv;
	}
	
	/**
     *    curvature direction estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][][][] principalCurvatureDirections(float[][][] image, int nx, int ny, int nz ) {
		float[][][][] curv = new float[8][nx][ny][nz];
		double num,den;
		double[][] hessian = new double[3][3];
		double[] eigen;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			
			hessian[0][0] = image[x+1][y][z] - 2.0*image[x][y][z] + image[x-1][y][z];
			hessian[1][1] = image[x][y+1][z] - 2.0*image[x][y][z] + image[x][y-1][z];
			hessian[2][2] = image[x][y][z+1] - 2.0*image[x][y][z] + image[x][y][z-1];
			
			hessian[0][1] = 0.25*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			hessian[1][2] = 0.25*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			hessian[2][0] = 0.25*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
			
			hessian[1][0] = 0.25*(image[x+1][y+1][z] + image[x-1][y-1][z] - image[x+1][y-1][z] - image[x-1][y+1][z]);
			hessian[2][1] = 0.25*(image[x][y+1][z+1] + image[x][y-1][z-1] - image[x][y+1][z-1] - image[x][y-1][z+1]);
			hessian[0][2] = 0.25*(image[x+1][y][z+1] + image[x-1][y][z-1] - image[x-1][y][z+1] - image[x+1][y][z-1]);
				
			eigen = Matrix3D.eigenvalues(hessian);
			hessian = Matrix3D.eigenvectors(hessian,eigen);
			
			curv[0][x][y][z] = (float)eigen[0];
			curv[1][x][y][z] = (float)hessian[0][0];
			curv[2][x][y][z] = (float)hessian[1][0];
			curv[3][x][y][z] = (float)hessian[2][0];
			curv[4][x][y][z] = (float)eigen[1];
			curv[5][x][y][z] = (float)hessian[0][1];
			curv[6][x][y][z] = (float)hessian[1][1];
			curv[7][x][y][z] = (float)hessian[2][1];
		}
		return curv;
	}

	/**
     *    curvature direction estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][] principalCurvatures(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[][] curv = new float[2][nx*ny*nz];
		double num,den;
		double[][] hessian = new double[3][3];
		double[] eigen;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			hessian[0][0] = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			hessian[1][1] = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			hessian[2][2] = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			hessian[0][1] = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			hessian[1][2] = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			hessian[2][0] = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
			
			hessian[1][0] = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			hessian[2][1] = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			hessian[0][2] = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
				
			eigen = Matrix3D.eigenvalues(hessian);
			//hessian = Matrix3D.eigenvectors(hessian,eigen);
			
			curv[0][xyz] = (float)eigen[0];
			curv[1][xyz] = (float)eigen[1];
		}
		return curv;
	}

	/**
     *    curvature direction estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][] principalCurvatureDirectionsJama(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[][] curv = new float[8][nx*ny*nz];
		double num,den;
		Matrix hessian = new Matrix(3,3);
		double[] eigen;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			hessian.set(0,0, image[xyzpx] - 2.0*image[xyz] + image[xyzmx]);
			hessian.set(1,1, image[xyzpy] - 2.0*image[xyz] + image[xyzmy]);
			hessian.set(2,2, image[xyzpz] - 2.0*image[xyz] + image[xyzmz]);
			
			hessian.set(0,1, 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]) );
			hessian.set(1,2, 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]) );
			hessian.set(2,0, 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]) );
			
			hessian.set(1,0, 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]) );
			hessian.set(2,1, 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]) );
			hessian.set(0,2, 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]) );
				
			EigenvalueDecomposition eig = hessian.eig();
			eigen = eig.getRealEigenvalues();
			hessian = eig.getV();
			int first = Numerics.argmaxmag(eigen[0],eigen[1],eigen[2]);
			int second = Numerics.argsecmag(eigen[0],eigen[1],eigen[2]);
			
			curv[0][xyz] = (float)eigen[first];
			curv[1][xyz] = (float)hessian.get(0,first);
			curv[2][xyz] = (float)hessian.get(1,first);
			curv[3][xyz] = (float)hessian.get(2,first);
			curv[4][xyz] = (float)eigen[second];
			curv[5][xyz] = (float)hessian.get(0,second);
			curv[6][xyz] = (float)hessian.get(1,second);
			curv[7][xyz] = (float)hessian.get(2,second);
		}
		return curv;
	}

	/**
     *    curvature direction estimation on a smooth distance function
	 *	  with no additionnal smoothing.
     *    @param 	image	the original distance image
	 *	  @return 			the diffused image
     */
    public static final float[][] principalCurvatureDirections(float[] image, boolean[] mask, int nx, int ny, int nz ) {
		float[][] curv = new float[8][nx*ny*nz];
		double num,den;
		double[][] hessian = new double[3][3];
		double[] eigen;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyz = x + nx*y + nx*ny*z;

			int xyzpx = xyz+1; 			if (!mask[xyzpx]) xyzpx = xyz;
			int xyzmx = xyz-1; 			if (!mask[xyzmx]) xyzmx = xyz;
			int xyzpy = xyz+nx; 		if (!mask[xyzpy]) xyzpy = xyz;
			int xyzmy = xyz-nx; 		if (!mask[xyzmy]) xyzmy = xyz;
			int xyzpz = xyz+nx*ny;		if (!mask[xyzpz]) xyzpz = xyz;
			int xyzmz = xyz-nx*ny; 		if (!mask[xyzmz]) xyzmz = xyz;
			
			int xyzpxpy = xyz+1+nx;		if (!mask[xyzpxpy]) xyzpxpy = xyz;
			int xyzmxpy = xyz-1+nx;		if (!mask[xyzmxpy]) xyzmxpy = xyz;
			int xyzpxmy = xyz+1-nx;		if (!mask[xyzpxmy]) xyzpxmy = xyz;
			int xyzmxmy = xyz-1-nx;		if (!mask[xyzmxmy]) xyzmxmy = xyz;
			
			int xyzpypz = xyz+nx+nx*ny;		if (!mask[xyzpypz]) xyzpypz = xyz;
			int xyzmypz = xyz-nx+nx*ny;		if (!mask[xyzmypz]) xyzmypz = xyz;
			int xyzpymz = xyz+nx-nx*ny;		if (!mask[xyzpymz]) xyzpymz = xyz;
			int xyzmymz = xyz-nx-nx*ny;		if (!mask[xyzmymz]) xyzmymz = xyz;
			
			int xyzpzpx = xyz+nx*ny+1;		if (!mask[xyzpzpx]) xyzpzpx = xyz;
			int xyzmzpx = xyz-nx*ny+1;		if (!mask[xyzmzpx]) xyzmzpx = xyz;
			int xyzpzmx = xyz+nx*ny-1;		if (!mask[xyzpzmx]) xyzpzmx = xyz;
			int xyzmzmx = xyz-nx*ny-1;		if (!mask[xyzmzmx]) xyzmzmx = xyz;
			
			hessian[0][0] = image[xyzpx] - 2.0*image[xyz] + image[xyzmx];
			hessian[1][1] = image[xyzpy] - 2.0*image[xyz] + image[xyzmy];
			hessian[2][2] = image[xyzpz] - 2.0*image[xyz] + image[xyzmz];
			
			hessian[0][1] = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			hessian[1][2] = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			hessian[2][0] = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
			
			hessian[1][0] = 0.25*(image[xyzpxpy] + image[xyzmxmy] - image[xyzpxmy] - image[xyzmxpy]);
			hessian[2][1] = 0.25*(image[xyzpypz] + image[xyzmymz] - image[xyzpymz] - image[xyzmypz]);
			hessian[0][2] = 0.25*(image[xyzpzpx] + image[xyzmzmx] - image[xyzpzmx] - image[xyzmzpx]);
				
			eigen = Matrix3D.eigenvalues(hessian);
			hessian = Matrix3D.eigenvectors(hessian,eigen);
			
			curv[0][xyz] = (float)eigen[0];
			curv[1][xyz] = (float)hessian[0][0];
			curv[2][xyz] = (float)hessian[1][0];
			curv[3][xyz] = (float)hessian[2][0];
			curv[4][xyz] = (float)eigen[1];
			curv[5][xyz] = (float)hessian[0][1];
			curv[6][xyz] = (float)hessian[1][1];
			curv[7][xyz] = (float)hessian[2][1];
		}
		return curv;
	}

	/**
     *    curvature estimation on a levelset distance function
	 *	  (no additionnal smoothing).
     *    @param 	image	the levelset distance function
	 *	  @return 			the curvature at that location
     */
    public static final float[] curvatureMetrics(float[] image, int xyz, int nx, int ny, int nz) {
		// first derivatives
		double D0x = (image[xyz+1] 		- image[xyz-1])/2.0f;
		double D0y = (image[xyz+nx] 	- image[xyz-nx])/2.0f;
		double D0z = (image[xyz+nx*ny] - image[xyz-nx*ny])/2.0f;
		
		// second derivatives
		double Dxx = image[xyz-1] + image[xyz+1] - 2.0*image[xyz];
		double Dyy = image[xyz-nx] + image[xyz+nx] - 2.0*image[xyz];
		double Dzz = image[xyz-nx*ny] + image[xyz+nx*ny] - 2.0*image[xyz];
	
		double Dxy = (image[xyz-1-nx] + image[xyz+1+nx] - image[xyz-1+nx] - image[xyz+1-nx])/4.0;
		double Dyz = (image[xyz-nx-nx*ny] + image[xyz+nx+nx*ny] - image[xyz-nx+nx*ny] - image[xyz+nx-nx*ny])/4.0;
		double Dzx = (image[xyz-nx*ny-1] + image[xyz+nx*ny+1] - image[xyz-nx*ny+1] - image[xyz+nx*ny-1])/4.0;
	
		// gradient norm
		double SD0x = D0x * D0x;
		double SD0y = D0y * D0y;
		double SD0z = D0z * D0z;
		double GPhi = FastMath.sqrt(SD0x + SD0y + SD0z);
	
		// mean curvature : K = k1 + k2;
		double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
					- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
		
		// gaussian curvature: G = k1*k2;
		double G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
					+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));

		double k1 = 0, k2 = 0, cn = 0, si = 0;			
		if(GPhi > 0.0000001){
			double tmp = GPhi*GPhi;
			K = K/(GPhi*tmp);
			G = G/(tmp*tmp);
			tmp = K*K-4.0*G;
			if (tmp>0) {
				tmp = FastMath.sqrt(tmp);
				double tmp2 = FastMath.sqrt(K*K-2.0*G);
				cn = tmp2/FastMath.sqrt(2.0);
				si = 2.0/FastMath.PI*FastMath.atan2(K,tmp);
				if (K>0) {
					k2 = K/2.0 + tmp;
					k1 = K/2.0 - tmp;
				} else {
					k2 = K/2.0 - tmp;
					k1 = K/2.0 + tmp;
				}
			}
		}
		return new float[]{(float)k1, (float)k2, (float)cn, (float)si};
	}
	
	public static final float gaussCurvatureValue(float[] image, int xyz, int nx, int ny, int nz) {
		// first derivatives
		double D0x = (image[xyz+1] 		- image[xyz-1])/2.0f;
		double D0y = (image[xyz+nx] 	- image[xyz-nx])/2.0f;
		double D0z = (image[xyz+nx*ny] - image[xyz-nx*ny])/2.0f;
		
		// second derivatives
		double Dxx = image[xyz-1] + image[xyz+1] - 2.0*image[xyz];
		double Dyy = image[xyz-nx] + image[xyz+nx] - 2.0*image[xyz];
		double Dzz = image[xyz-nx*ny] + image[xyz+nx*ny] - 2.0*image[xyz];
	
		double Dxy = (image[xyz-1-nx] + image[xyz+1+nx] - image[xyz-1+nx] - image[xyz+1-nx])/4.0;
		double Dyz = (image[xyz-nx-nx*ny] + image[xyz+nx+nx*ny] - image[xyz-nx+nx*ny] - image[xyz+nx-nx*ny])/4.0;
		double Dzx = (image[xyz-nx*ny-1] + image[xyz+nx*ny+1] - image[xyz-nx*ny+1] - image[xyz+nx*ny-1])/4.0;
	
		// gradient norm
		double SD0x = D0x * D0x;
		double SD0y = D0y * D0y;
		double SD0z = D0z * D0z;
		double GPhi = FastMath.sqrt(SD0x + SD0y + SD0z);
	
		// gaussian curvature
		double G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
					+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));
		
		double den = GPhi*GPhi*GPhi*GPhi;
		
		if (den>0) return (float)(G/den);
		else return 0.0f;
	}
		
	public static final float meanCurvatureValue(float[] image, int xyz, int nx, int ny, int nz) {
		// first derivatives
		double D0x = (image[xyz+1] 		- image[xyz-1])/2.0f;
		double D0y = (image[xyz+nx] 	- image[xyz-nx])/2.0f;
		double D0z = (image[xyz+nx*ny] - image[xyz-nx*ny])/2.0f;
		
		// second derivatives
		double Dxx = image[xyz-1] + image[xyz+1] - 2.0*image[xyz];
		double Dyy = image[xyz-nx] + image[xyz+nx] - 2.0*image[xyz];
		double Dzz = image[xyz-nx*ny] + image[xyz+nx*ny] - 2.0*image[xyz];
	
		double Dxy = (image[xyz-1-nx] + image[xyz+1+nx] - image[xyz-1+nx] - image[xyz+1-nx])/4.0;
		double Dyz = (image[xyz-nx-nx*ny] + image[xyz+nx+nx*ny] - image[xyz-nx+nx*ny] - image[xyz+nx-nx*ny])/4.0;
		double Dzx = (image[xyz-nx*ny-1] + image[xyz+nx*ny+1] - image[xyz-nx*ny+1] - image[xyz+nx*ny-1])/4.0;
	
		// gradient norm
		double SD0x = D0x * D0x;
		double SD0y = D0y * D0y;
		double SD0z = D0z * D0z;
		double GPhi = FastMath.sqrt(SD0x + SD0y + SD0z);
	
		// mean curvature
		double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
					- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
		
		double den = 2.0*GPhi*GPhi*GPhi;
		
		if (den>0) return (float)(K/den);
		else return 0.0f;
	}
	
	public static final float curvednessValue(float[] image, int xyz, int nx, int ny, int nz) {
		// first derivatives
		double D0x = (image[xyz+1] 		- image[xyz-1])/2.0f;
		double D0y = (image[xyz+nx] 	- image[xyz-nx])/2.0f;
		double D0z = (image[xyz+nx*ny] - image[xyz-nx*ny])/2.0f;
		
		// second derivatives
		double Dxx = image[xyz-1] + image[xyz+1] - 2.0*image[xyz];
		double Dyy = image[xyz-nx] + image[xyz+nx] - 2.0*image[xyz];
		double Dzz = image[xyz-nx*ny] + image[xyz+nx*ny] - 2.0*image[xyz];
	
		double Dxy = (image[xyz-1-nx] + image[xyz+1+nx] - image[xyz-1+nx] - image[xyz+1-nx])/4.0;
		double Dyz = (image[xyz-nx-nx*ny] + image[xyz+nx+nx*ny] - image[xyz-nx+nx*ny] - image[xyz+nx-nx*ny])/4.0;
		double Dzx = (image[xyz-nx*ny-1] + image[xyz+nx*ny+1] - image[xyz-nx*ny+1] - image[xyz+nx*ny-1])/4.0;
	
		// gradient norm
		double SD0x = D0x * D0x;
		double SD0y = D0y * D0y;
		double SD0z = D0z * D0z;
		double GPhi = FastMath.sqrt(SD0x + SD0y + SD0z);

		if (GPhi==0) return 0.0f;
		
		// mean curvature
		double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
					- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
		
		double denK = 2.0*GPhi*GPhi*GPhi;
		
		// gaussian curvature
		double G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
					+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));

		double denG = GPhi*GPhi*GPhi*GPhi;
		
		if (denK>0 && denG>0) {		
			return (float)FastMath.sqrt(Numerics.abs(2.0*K/denK*K/denK-G/denG));
		} else {
			return 0.0f;
		}
	}
	
	
	public static final float shapeIndexValue(float[] image, int xyz, int nx, int ny, int nz) {
		// first derivatives
		double D0x = (image[xyz+1] 		- image[xyz-1])/2.0f;
		double D0y = (image[xyz+nx] 	- image[xyz-nx])/2.0f;
		double D0z = (image[xyz+nx*ny] - image[xyz-nx*ny])/2.0f;
		
		// second derivatives
		double Dxx = image[xyz-1] + image[xyz+1] - 2.0*image[xyz];
		double Dyy = image[xyz-nx] + image[xyz+nx] - 2.0*image[xyz];
		double Dzz = image[xyz-nx*ny] + image[xyz+nx*ny] - 2.0*image[xyz];
	
		double Dxy = (image[xyz-1-nx] + image[xyz+1+nx] - image[xyz-1+nx] - image[xyz+1-nx])/4.0;
		double Dyz = (image[xyz-nx-nx*ny] + image[xyz+nx+nx*ny] - image[xyz-nx+nx*ny] - image[xyz+nx-nx*ny])/4.0;
		double Dzx = (image[xyz-nx*ny-1] + image[xyz+nx*ny+1] - image[xyz-nx*ny+1] - image[xyz+nx*ny-1])/4.0;
	
		// gradient norm
		double SD0x = D0x * D0x;
		double SD0y = D0y * D0y;
		double SD0z = D0z * D0z;
		double GPhi = FastMath.sqrt(SD0x + SD0y + SD0z);

		if (GPhi==0) return 0.0f;
		
		// mean curvature
		double K =  (Dyy + Dzz)*SD0x + (Dxx + Dzz)*SD0y + (Dxx + Dyy)*SD0z 
					- 2.0*(D0x*D0y*Dxy + D0y*D0z*Dyz + D0z*D0x*Dzx);
		
		double denK = 2.0*GPhi*GPhi*GPhi;
		
		// gaussian curvature
		double G = (Dyy*Dzz - Dyz*Dyz)*SD0x + (Dzz*Dxx - Dzx*Dzx)*SD0y + (Dxx*Dyy - Dxy*Dxy)*SD0z 
					+ 2.0*(D0x*D0y*(Dyz*Dzx - Dxy*Dzz) + D0z*D0x*(Dxy*Dyz - Dzx*Dyy) + D0y*D0z*(Dxy*Dzx - Dyz*Dxx));

		double denG = GPhi*GPhi*GPhi*GPhi;
		
		if (denK>0 && denG>0) {		
			return (float)(2.0/FastMath.PI*FastMath.atan2(K/denK, FastMath.sqrt(Numerics.abs(K/denK*K/denK-G/denG))));
		} else {
			return 0.0f;
		}
	}
	
	/** compute a quadric approximation to a levelset surface, assuming a centered paraboloid */
	public static final double[] quadricLevelsetApproximation(float[][][] image, boolean[][][] mask, int x0, int y0, int z0, int size, int nx, int ny, int nz) {
		// 1. build the linear estimate's vector and matrix
		Matrix X = new Matrix(9,9);
		Matrix Y = new Matrix(9,1);
		double[] vect = new double[9];
		double sig2 = (size/3.0)*(size/3.0);
		
		for (int n=0;n<9;n++) {
			for (int m=0;m<9;m++) X.set(n,m, 0);
			Y.set(n,0, 0);
		}
		
		for (int i=-size;i<=size;i++) {
			for (int j=-size;j<=size;j++) {
				for (int l=-size;l<=size;l++) {
					if (mask[x0+i][y0+j][z0+l]) {
						vect[0] = i*i; vect[1] = j*j; vect[2] = l*l;
						vect[3] = i*j; vect[4] = j*l; vect[5] = l*i;
						vect[6] = 2*i; vect[7] = 2*j; vect[8] = 2*l;
						
						double dist = (image[x0+i][y0+j][z0+l]-image[x0][y0][z0]);
						// distance weighting: lower the contribution of value further from the surface
						// and the distance to center point as well?
						double w = FastMath.exp( -0.5*dist*dist/sig2)*FastMath.exp( -0.5*(i*i+j*j+l*l)/sig2);
							
						for (int n=0;n<9;n++)  {
							for (int m=0;m<9;m++) X.set(n,m, X.get(n,m) + w*vect[n]*vect[m]);
							Y.set(n,0, Y.get(n,0) + w*vect[n]*dist);	
						}
					}
				}
			}
		}
		
		// 2. Get the parameters
		Matrix A = X.solve(Y);
		for (int n=0;n<9;n++) {
			vect[n] = A.get(n,0);
		}
		return vect;
	}
	
	/** get curvatures from a quadric approximation */
	public static final float[] quadricApproximationCurvature(double[] quadric) {
		// build the quadratic term matrix
		Matrix A = new Matrix(3,3);
		A.set(0,0, quadric[0]);
		A.set(1,1, quadric[1]);
		A.set(2,2, quadric[2]);
		
		A.set(0,1, quadric[3]);
		A.set(1,2, quadric[4]);
		A.set(2,0, quadric[5]);
		
		A.set(1,0, quadric[3]);
		A.set(2,1, quadric[4]);
		A.set(0,2, quadric[5]);
		
		// decompose
		EigenvalueDecomposition eig = A.eig();
		double[] eigen = eig.getRealEigenvalues();
		A = eig.getV();
		
		// find the normal vector: best alignment with linear term
		int norm = 0;
		if ( Numerics.abs(A.get(0,1)*quadric[6]+A.get(1,1)*quadric[7]+A.get(2,1)*quadric[8])
			>Numerics.abs(A.get(0,norm)*quadric[6]+A.get(1,norm)*quadric[7]+A.get(2,norm)*quadric[8]) ) norm=1;
		
		if ( Numerics.abs(A.get(0,2)*quadric[6]+A.get(1,2)*quadric[7]+A.get(2,2)*quadric[8])
			>Numerics.abs(A.get(0,norm)*quadric[6]+A.get(1,norm)*quadric[7]+A.get(2,norm)*quadric[8]) ) norm=2;
		
		if (eigen[norm]>0.1) System.out.print("!");
		eigen[norm] = 0.0;
		
		int first = Numerics.argmaxmag(eigen[0],eigen[1],eigen[2]);
		int second = Numerics.argsecmag(eigen[0],eigen[1],eigen[2]);
		
		float[] curv = new float[8];
		curv[0] = (float)eigen[first];
		curv[1] = (float)A.get(0,first);
		curv[2] = (float)A.get(1,first);
		curv[3] = (float)A.get(2,first);
		curv[4] = (float)eigen[second];
		curv[5] = (float)A.get(0,second);
		curv[6] = (float)A.get(1,second);
		curv[7] = (float)A.get(2,second);
		
		return curv;
	}
	
	/** compute a quadric approximation to a levelset surface, assuming a centered paraboloid */
	public static final float[][][][] quadricCurvatureEstimates(float[][][] image, boolean[][][] mask, int size, int nx, int ny, int nz) {
		Matrix X = new Matrix(9,9);
		Matrix Y = new Matrix(9,1);
		double[] vect = new double[9];
		double sig2 = (size/3.0)*(size/3.0);
		float[][][][] curv = new float[8][nx][ny][nz];
		
		for (int x=size;x<nx-size;x++) for (int y=size;y<ny-size;y++) for (int z=size;z<nz-size;z++) if (mask[x][y][z]) {
			// 1. build the linear estimate's vector and matrix
			for (int n=0;n<9;n++) {
				for (int m=0;m<9;m++) X.set(n,m, 0);
				Y.set(n,0, 0);
			}
			
			for (int i=-size;i<=size;i++) {
				for (int j=-size;j<=size;j++) {
					for (int l=-size;l<=size;l++) {
						if (mask[x+i][y+j][z+l]) {
							vect[0] = i*i; vect[1] = j*j; vect[2] = l*l;
							vect[3] = i*j; vect[4] = j*l; vect[5] = l*i;
							vect[6] = 2*i; vect[7] = 2*j; vect[8] = 2*l;
							
							double dist = (image[x+i][y+j][z+l]-image[x][y][z]);
							// distance weighting: lower the contribution of value further from the surface
							// and the distance to center point as well?
							double w = FastMath.exp( -0.5*dist*dist/sig2)*FastMath.exp( -0.5*(i*i+j*j+l*l)/sig2);
							
							for (int n=0;n<9;n++)  {
								for (int m=0;m<9;m++) X.set(n,m, X.get(n,m) + w*vect[n]*vect[m]);
								Y.set(n,0, Y.get(n,0) + w*vect[n]*dist);	
							}
						}
					}
				}
			}
		
			// 2. Get the parameters
			Matrix Q = X.solve(Y);

			// 3. build the quadratic term matrix
			Matrix A = new Matrix(3,3);
			A.set(0,0, Q.get(0,0));
			A.set(1,1, Q.get(1,0));
			A.set(2,2, Q.get(2,0));
			
			A.set(0,1, Q.get(3,0));
			A.set(1,2, Q.get(4,0));
			A.set(2,0, Q.get(5,0));
			
			A.set(1,0, Q.get(3,0));
			A.set(2,1, Q.get(4,0));
			A.set(0,2, Q.get(5,0));
			
			// decompose
			EigenvalueDecomposition eig = A.eig();
			double[] eigen = eig.getRealEigenvalues();
			A = eig.getV();
			
			// find the normal vector: best alignment with linear term
			int norm = 0;
			if ( Numerics.abs(A.get(0,1)*Q.get(6,0)+A.get(1,1)*Q.get(7,0)+A.get(2,1)*Q.get(8,0))
				>Numerics.abs(A.get(0,norm)*Q.get(6,0)+A.get(1,norm)*Q.get(7,0)+A.get(2,norm)*Q.get(8,0)) ) norm=1;
			
			if ( Numerics.abs(A.get(0,2)*Q.get(6,0)+A.get(1,2)*Q.get(7,0)+A.get(2,2)*Q.get(8,0))
				>Numerics.abs(A.get(0,norm)*Q.get(6,0)+A.get(1,norm)*Q.get(7,0)+A.get(2,norm)*Q.get(8,0)) ) norm=2;
			
			double neig = eigen[norm];
			eigen[norm] = 0.0;
			
			int first = Numerics.argmaxmag(eigen[0],eigen[1],eigen[2]);
			int second = Numerics.argsecmag(eigen[0],eigen[1],eigen[2]);
			
			if (Numerics.abs(eigen[norm]/eigen[first])>0.1) System.out.print("!");
			
			curv[0][x][y][z] = (float)eigen[first];
			curv[1][x][y][z] = (float)A.get(0,first);
			curv[2][x][y][z] = (float)A.get(1,first);
			curv[3][x][y][z] = (float)A.get(2,first);
			curv[4][x][y][z] = (float)eigen[second];
			curv[5][x][y][z] = (float)A.get(0,second);
			curv[6][x][y][z] = (float)A.get(1,second);
			curv[7][x][y][z] = (float)A.get(2,second);
		}
		return curv;
	}
	
	/** compute a quadric approximation to a levelset surface, assuming a centered paraboloid */
	public static final float[][] quadricCurvatureEstimates(float[] image, boolean[] sampling, boolean[] mask, int size, int nx, int ny, int nz) {
		Matrix X = new Matrix(9,9);
		Matrix Y = new Matrix(9,1);
		double[] vect = new double[9];
		double sig2 = (size/3.0)*(size/3.0);
		float[][] curv = new float[8][nx*ny*nz];
		
		for (int x=size;x<nx-size;x++) for (int y=size;y<ny-size;y++) for (int z=size;z<nz-size;z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (sampling[xyz] && mask[xyz]) {
				// 1. build the linear estimate's vector and matrix
				for (int n=0;n<9;n++) {
					for (int m=0;m<9;m++) X.set(n,m, 0);
					Y.set(n,0, 0);
				}
				
				for (int i=-size;i<=size;i++) {
					for (int j=-size;j<=size;j++) {
						for (int l=-size;l<=size;l++) {
							if (mask[xyz+i+nx*j+nx*ny*l]) {
								vect[0] = i*i; vect[1] = j*j; vect[2] = l*l;
								vect[3] = i*j; vect[4] = j*l; vect[5] = l*i;
								vect[6] = 2*i; vect[7] = 2*j; vect[8] = 2*l;
								
								double dist = (image[xyz+i+nx*j+nx*ny*l]-image[xyz]);
								// distance weighting: lower the contribution of value further from the surface
								// and the distance to center point as well?
								double w = FastMath.exp( -0.5*dist*dist/sig2)*FastMath.exp( -0.5*(i*i+j*j+l*l)/sig2);
								
								for (int n=0;n<9;n++)  {
									for (int m=0;m<9;m++) X.set(n,m, X.get(n,m) + w*vect[n]*vect[m]);
									Y.set(n,0, Y.get(n,0) + w*vect[n]*dist);	
								}
							}
						}
					}
				}
			
				// 2. Get the parameters
				Matrix Q = X.solve(Y);
	
				// 3. build the quadratic term matrix
				Matrix A = new Matrix(3,3);
				A.set(0,0, Q.get(0,0));
				A.set(1,1, Q.get(1,0));
				A.set(2,2, Q.get(2,0));
				
				A.set(0,1, Q.get(3,0));
				A.set(1,2, Q.get(4,0));
				A.set(2,0, Q.get(5,0));
				
				A.set(1,0, Q.get(3,0));
				A.set(2,1, Q.get(4,0));
				A.set(0,2, Q.get(5,0));
				
				// decompose
				EigenvalueDecomposition eig = A.eig();
				double[] eigen = eig.getRealEigenvalues();
				A = eig.getV();
				
				// find the normal vector: best alignment with linear term
				int norm = 0;
				if ( Numerics.abs(A.get(0,1)*Q.get(6,0)+A.get(1,1)*Q.get(7,0)+A.get(2,1)*Q.get(8,0))
					>Numerics.abs(A.get(0,norm)*Q.get(6,0)+A.get(1,norm)*Q.get(7,0)+A.get(2,norm)*Q.get(8,0)) ) norm=1;
				
				if ( Numerics.abs(A.get(0,2)*Q.get(6,0)+A.get(1,2)*Q.get(7,0)+A.get(2,2)*Q.get(8,0))
					>Numerics.abs(A.get(0,norm)*Q.get(6,0)+A.get(1,norm)*Q.get(7,0)+A.get(2,norm)*Q.get(8,0)) ) norm=2;
				
				double neig = eigen[norm];
				eigen[norm] = 0.0;
			
				int first = Numerics.argmaxmag(eigen[0],eigen[1],eigen[2]);
				int second = Numerics.argsecmag(eigen[0],eigen[1],eigen[2]);
				
				//if (Numerics.abs(eigen[norm]/eigen[first])>0.1) System.out.print("!");
			
				curv[0][xyz] = (float)eigen[first];
				curv[1][xyz] = (float)A.get(0,first);
				curv[2][xyz] = (float)A.get(1,first);
				curv[3][xyz] = (float)A.get(2,first);
				curv[4][xyz] = (float)eigen[second];
				curv[5][xyz] = (float)A.get(0,second);
				curv[6][xyz] = (float)A.get(1,second);
				curv[7][xyz] = (float)A.get(2,second);
			}
		}
		return curv;
	}
	

}
