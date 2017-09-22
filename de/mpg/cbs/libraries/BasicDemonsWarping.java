package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
 
/**
 *
 *  This algorithm handles registration algorithms
 *	of images with the Demons algorithm
 *	and a multiscale / multigrid technique
 *	(original Demons method, not robust)
 *
 *
 *	@version    November 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class BasicDemonsWarping {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// convenience tags
	private static final	int   X=0;
	private static final	int   Y=1;
	private static final	int   Z=2;
	private static final	int   T=3;
	
	// flag labels
	public static final	int   FLUID 		= 301;
	public static final	int   DIFFUSION 	= 302;
	public static final	int   MIXED	 		= 303;
	
	// data buffers
	private 	float[]		source;  		// source image
	private 	float[]		target;			// target image
	private		float[][]		c;				// added transform from target space to original space
	private		float[][]		s;				// added transform from target space to original space
	private		float[][]		u;				// added transform update from target space to original space
	private static	int		nix,niy,niz;   	// image dimensions
	private static	int		nsx,nsy,nsz;   	// scaled transform dimensions
	private static	float 	scale;
	private static	float	rix,riy,riz;   	// image resolutions
	private 		float[][]	transform;		// prior transform matrix (for instance from prior alignment)		
	
	// parameters
	private		float		smoothing; 	// smoothing kernel size
	private		int		regType;
	
	// computation variables
	private		float			sigma2; 		// scale of transformation
	

	private		float[][]		kernel;
	private		int				gx,gy,gz;
    
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=false;
	static final boolean		verbose=true;
    
	/**
	 *  constructor
     *  note: the membership function gain_ must have sum = 0 on the masked areas
	 */
	public BasicDemonsWarping(float[] src_, float[] trg_,
								int nx_, int ny_, int nz_,
								float rx_, float ry_, float rz_,
								float smoothing_,
								float spScale_,
								float scale_,
								int reg_,
								float[][] trans_) {
	
		source = src_;
		target = trg_;
		
		scale = scale_;
		
		smoothing = smoothing_;
		sigma2 = spScale_*spScale_;
		
		transform = trans_;
		
		nix = nx_;
		niy = ny_;
		niz = nz_;
		
		rix = rx_;
		riy = ry_;
		riz = rz_;
					
		nsx = Numerics.ceil(nix/scale);
		nsy = Numerics.ceil(niy/scale);
		nsz = Numerics.ceil(niz/scale);
		
		regType = reg_;

		//initialize the smoothing kernel
		kernel = ImageFilters.separableGaussianKernel(smoothing/rix,smoothing/riy,smoothing/riz);
		gx = (kernel[X].length-1)/2;
		gy = (kernel[Y].length-1)/2;
		gz = (kernel[Z].length-1)/2;
				
		isWorking = true;
		if (debug) BasicInfo.displayMessage("Basic Demons:initialisation\n");
	}

	final public void finalize() {
		s = null; u = null; c = null;
		System.gc();
	}
    
   public final float[][] getCurrentTransform() {
	   return s;
   }
   public final float[][] getCurrentUpdate() {
		return u;
	}
   public final float[] getSource() {
	   return source;
   }
   public final float[] getTarget() {
	   return target;
   }
   
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
	
    /** initialize the transform from previous estimate, if exists */
    public final void initializeTransform() {
				
		// initialization: start from prior transform or zero
		s = new float[3][nsx*nsy*nsz];
		if (transform==null) {
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				s[X][xyz] = x*scale;
				s[Y][xyz] = y*scale;
				s[Z][xyz] = z*scale;
			}
		} else {
			if (debug) BasicInfo.displayMessage("transform matrix: \n ["+transform[X][X]+", "+transform[X][Y]+", "+transform[X][Z]+", "+transform[X][T]+"]\n ["
															   +transform[Y][X]+", "+transform[Y][Y]+", "+transform[Y][Z]+", "+transform[Y][T]+"]\n ["
															   +transform[Z][X]+", "+transform[Z][Y]+", "+transform[Z][Z]+", "+transform[Z][T]+"]\n");
			
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				s[X][xyz] = transform[X][X]*x*scale + transform[X][Y]*y*scale + transform[X][Z]*z*scale + transform[X][T];
				s[Y][xyz] = transform[Y][X]*x*scale + transform[Y][Y]*y*scale + transform[Y][Z]*z*scale + transform[Y][T];
				s[Z][xyz] = transform[Z][X]*x*scale + transform[Z][Y]*y*scale + transform[Z][Z]*z*scale + transform[Z][T];
			}
		}
				
		// update always zero
		u = new float[3][nsx*nsy*nsz];
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			u[X][xyz] = 0.0f;
			u[Y][xyz] = 0.0f;
			u[Z][xyz] = 0.0f;
		}
		
		// composite update always zero
		c = new float[3][nsx*nsy*nsz];
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			c[X][xyz] = 0.0f;
			c[Y][xyz] = 0.0f;
			c[Z][xyz] = 0.0f;
		}
		
		return;
    }
    
	
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public void registerImageToTarget() {
		
		if (debug) BasicInfo.displayMessage("update: ");
		
		float meanDiff = 0.0f;
		for (int x=1;x<nsx-1;x++) for (int y=1;y<nsy-1;y++) for (int z=1;z<nsz-1;z++) {
			int xyz = x+nsx*y+nsx*nsy*z;
		
			// compute the update field
			float xs = s[X][xyz];
			float ys = s[Y][xyz];
			float zs = s[Z][xyz];
		
			float xsmx = s[X][xyz-1];
			float ysmx = s[Y][xyz-1];
			float zsmx = s[Z][xyz-1];
	
			float xspx = s[X][xyz+1];
			float yspx = s[Y][xyz+1];
			float zspx = s[Z][xyz+1];
	
			float xsmy = s[X][xyz-nsx];
			float ysmy = s[Y][xyz-nsx];
			float zsmy = s[Z][xyz-nsx];
	
			float xspy = s[X][xyz+nsx];
			float yspy = s[Y][xyz+nsx];
			float zspy = s[Z][xyz+nsx];
	
			float xsmz = s[X][xyz-nsx*nsy];
			float ysmz = s[Y][xyz-nsx*nsy];
			float zsmz = s[Z][xyz-nsx*nsy];
	
			float xspz = s[X][xyz+nsx*nsy];
			float yspz = s[Y][xyz+nsx*nsy];
			float zspz = s[Z][xyz+nsx*nsy];
			
			float xt = x*scale;
			float yt = y*scale;
			float zt = z*scale;
			
			u[X][xyz] = 0.0f;
			u[Y][xyz] = 0.0f;
			u[Z][xyz] = 0.0f;
			float den = 0.0f;
				
			float diff = 0.0f;
			float diffpx = 0.0f; float diffpy = 0.0f; float diffpz = 0.0f;
			float diffmx = 0.0f; float diffmy = 0.0f; float diffmz = 0.0f;
			
			diff = Numerics.max(diff,ImageInterpolation.linearClosestInterpolation(source,xs,ys,zs,nix,niy,niz));
			diffpx = Numerics.max(diffpx,ImageInterpolation.linearClosestInterpolation(source,xspx,yspx,zspx,nix,niy,niz));
			diffpy = Numerics.max(diffpy,ImageInterpolation.linearClosestInterpolation(source,xspy,yspy,zspy,nix,niy,niz));
			diffpz = Numerics.max(diffpz,ImageInterpolation.linearClosestInterpolation(source,xspz,yspz,zspz,nix,niy,niz));
			diffmx = Numerics.max(diffmx,ImageInterpolation.linearClosestInterpolation(source,xsmx,ysmx,zsmx,nix,niy,niz));
			diffmy = Numerics.max(diffmy,ImageInterpolation.linearClosestInterpolation(source,xsmy,ysmy,zsmy,nix,niy,niz));
			diffmz = Numerics.max(diffmz,ImageInterpolation.linearClosestInterpolation(source,xsmz,ysmz,zsmz,nix,niy,niz));
			
			diff = (ImageInterpolation.linearClosestInterpolation(target,xt,yt,zt,nix,niy,niz) - diff);
			
			float Jx, Jy, Jz;
			// symmetric ? no
			Jx = 0.5f/rix*(diffpx-diffmx);
			Jy = 0.5f/riy*(diffpy-diffmy);
			Jz = 0.5f/riz*(diffpz-diffmz);
		
			float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
			den = Numerics.max(ZERO, diff*diff/sigma2 + J2);
			
			u[X][xyz] += diff*Jx/den;
			u[Y][xyz] += diff*Jy/den;
			u[Z][xyz] += diff*Jz/den;
			
			meanDiff += diff;
		}
		meanDiff /= (nix*niy*niz);
		
		if (regType==FLUID || regType==MIXED) {
			if (debug) BasicInfo.displayMessage("GAUSS_FLUID regularization \n");
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFilters.separableConvolution(u[X],nsx,nsy,nsz,kernel,gx,gy,gz);
			u[Y] = ImageFilters.separableConvolution(u[Y],nsx,nsy,nsz,kernel,gx,gy,gz);
			u[Z] = ImageFilters.separableConvolution(u[Z],nsx,nsy,nsz,kernel,gx,gy,gz);
		}
		
		// compose the transformations
		if (debug) BasicInfo.displayMessage("compose with current transform \n");
						
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x+nsx*y+nsx*nsy*z;
			
			float xu = x+u[X][xyz];
			float yu = y+u[Y][xyz];
			float zu = z+u[Z][xyz];
			
			// note: if outside, extrapolate as X+u
			c[X][xyz] = ImageInterpolation.linearClosestInterpolation(s[X],xu,yu,zu,nsx,nsy,nsz) - x*scale*rix/rix;
			c[Y][xyz] = ImageInterpolation.linearClosestInterpolation(s[Y],xu,yu,zu,nsx,nsy,nsz) - y*scale*riy/riy;
			c[Z][xyz] = ImageInterpolation.linearClosestInterpolation(s[Z],xu,yu,zu,nsx,nsy,nsz) - z*scale*riz/riz;
		}
		
		if (regType==DIFFUSION || regType==MIXED) {
			if (debug) BasicInfo.displayMessage("GAUSS_DIFFUSION regularization \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFilters.separableConvolution(c[X],nsx,nsy,nsz,kernel,gx,gy,gz);
			c[Y] = ImageFilters.separableConvolution(c[Y],nsx,nsy,nsz,kernel,gx,gy,gz);
			c[Z] = ImageFilters.separableConvolution(c[Z],nsx,nsy,nsz,kernel,gx,gy,gz);
		}
					
		float meanC = 0.0f;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x+nsx*y+nsx*nsy*z;
			
			s[X][xyz] = x*scale*rix/rix + c[X][xyz];
			s[Y][xyz] = y*scale*riy/riy + c[Y][xyz];
			s[Z][xyz] = z*scale*riz/riz + c[Z][xyz];
			
			meanC += c[X][xyz]*c[X][xyz]+c[Y][xyz]*c[Y][xyz]+c[Z][xyz]*c[Z][xyz];
		}
		meanC /= (nsx*nsy*nsz);
		
		if (debug) BasicInfo.displayMessage("convergence "+meanC+" -> "+meanDiff+"\n");

        return;
    } // 
    
    /** 
	 *	returns the transformed coordinates
	 */
	public final void getCurrentMapping(float[] Xs, float x, float y, float z) {
		float xs = x/scale;
		float ys = y/scale;
		float zs = z/scale;
		
		Xs[X] = ImageInterpolation.linearClosestInterpolation(s[X],xs,ys,zs,nsx,nsy,nsz);
		Xs[Y] = ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys,zs,nsx,nsy,nsz);
		Xs[Z] = ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys,zs,nsx,nsy,nsz);
		return;
	}// getCurrentMapping

    /** 
	 *	rotates a transformed unit vector
	 */
	public final void mapDirection(float[] dir, int x, int y, int z) {
		float xs = x/scale;
		float ys = y/scale;
		float zs = z/scale;

		// Jacobian
		float[][] Ds = new float[3][3];
		Ds[X][X] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[X],xs+1,ys,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xs-1,ys,zs,nsx,nsy,nsz));
		Ds[Y][X] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[X],xs,ys+1,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xs,ys-1,zs,nsx,nsy,nsz));
		Ds[Z][X] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[X],xs,ys,zs+1,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xs,ys,zs-1,nsx,nsy,nsz));

		Ds[X][Y] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Y],xs+1,ys,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xs-1,ys,zs,nsx,nsy,nsz));
		Ds[Y][Y] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys+1,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys-1,zs,nsx,nsy,nsz));
		Ds[Z][Y] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys,zs+1,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys,zs-1,nsx,nsy,nsz));

		Ds[X][Z] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Z],xs+1,ys,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xs-1,ys,zs,nsx,nsy,nsz));
		Ds[Y][Z] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys+1,zs,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys-1,zs,nsx,nsy,nsz));
		Ds[Z][Z] = 0.5f*(ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys,zs+1,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys,zs-1,nsx,nsy,nsz));
		
		// transform the vector
		float[] vec = new float[3];
		vec[X] = Ds[X][X]*dir[X] + Ds[X][Y]*dir[Y] + Ds[X][Z]*dir[Z];
		vec[Y] = Ds[Y][X]*dir[X] + Ds[Y][Y]*dir[Y] + Ds[Y][Z]*dir[Z];
		vec[Z] = Ds[Z][X]*dir[X] + Ds[Z][Y]*dir[Y] + Ds[Z][Z]*dir[Z];
		float den = Numerics.max(ZERO, (float)Math.sqrt(vec[X]*vec[X] + vec[Y]*vec[Y] + vec[Z]*vec[Z]));
		
		// result
		dir[X] /= den;
		dir[Y] /= den;
		dir[Z] /= den;
		
		return;
	}// mapDirection

    /** 
	 *	returns the transformed image
	 */
	public final float[] exportTransformedSource() {
		float 	xs,ys,zs;
		float[]	src;
		
		src = new float[nix*niy*niz];
		
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
		
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			// compute interpolated values
			src[xyz] = ImageInterpolation.linearClosestInterpolation(source,xs,ys,zs,nix,niy,niz);
		}
		return src;
	} // exportTransformedImage

    /** 
	 *	returns the transform field v defined as X' = X+v (=s)
	 */
	public final float[][] exportTransformField() {
		float 	xs,ys,zs;
		float[][]	vec = new float[3][nix*niy*niz];
		
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
			
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			vec[X][xyz] = xs-x;
			vec[Y][xyz] = ys-y;
			vec[Z][xyz] = zs-z;			
 		}
		return vec;
	} // exportTransformedImage

}
