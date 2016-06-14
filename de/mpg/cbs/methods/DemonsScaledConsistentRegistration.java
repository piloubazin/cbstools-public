package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;

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
 *	@version    nay 2010
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class DemonsScaledConsistentRegistration {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// convenience tags
	private static final	int   X=0;
	private static final	int   Y=1;
	private static final	int   Z=2;
	private static final	int   T=3;
	
	// flag labels
	public static final	int   COMPOSITIVE 		= 101;
	public static final	int   DIFFEOMORPHIC 	= 102;
	
	public static final	int   FIXED 			= 201;
	public static final	int   MOVING 			= 202;
	public static final	int   SYMMETRIC 		= 203;
	
	public static final	int   GAUSS_FLUID 		= 301;
	public static final	int   GAUSS_DIFFUSION 	= 302;
	public static final	int   GAUSS_MIXED	 	= 303;
	
	public static final int	  RAW				= 401;
	public static final int	  NORMALIZED		= 402;
	public static final int	  SCALED			= 403;
	public static final int	  SEGMENTATION		= 404;
	public static final int	  FCM				= 405;
	public static final int	  MULTICHANNEL		= 406;
	
	// data buffers
	private 	float[][]		image;  		// source image
	private 	float[][]	    target;			// target image: best memberships
	private		float[][]		c;				// added transform from target space to original space (multiscale)
	private		float[][]		s;				// added transform from target space to original space (multiscale)
	private		float[][]		u;				// added transform update from target space to original space (multiscale)
	private		float[][]		ic;				// added transform from target space to original space (multiscale)
	private		float[][]		is;				// added transform from target space to original space (multiscale)
	private		float[][]		iu;				// added transform update from target space to original space (multiscale)
	private static	int		nix,niy,niz;   	// image dimensions 
	private static	int		ntx,nty,ntz;   	// target dimensions
	private static	int		nsx,nsy,nsz;   	// forward deformation dimensions
	private static	int		nisx,nisy,nisz;   	// backward deformation dimensions
	private static	float 	scale;
	private static	float	rix,riy,riz;   	// image resolutions (no pyramid)
	private static	float	rtx,rty,rtz;   	// target resolutions (no pyramid)
	private static 	int		nc;
	private		float		imin, imax;
	private		float		tmin, tmax;
	private		float[]		lb;
	private 		float[][]	transform;		// prior transform matrix (for instance from prior alignment)		
	private 		float[][]	itransform;		// prior transform inverse matrix (for instance from prior alignment)		
	
	// parameters
	private		float		smoothingKernel; 	// smoothing kernel size
	private		float		spatialScale; 		// scale of transformation
	private		boolean		useMask;
	private		float  		maskingThreshold;			// background threshold
	private		int		regType;
	private		int		forceType;
	private		int		fieldType;
	private		int		preType;
	private		int		Niter;
	
	// computation variables
	private		float			sigma2; 		// scale of transformation
	

	private		float[][]		gaussKernel;
	private		int				gx,gy,gz;
	private		float[][]		igaussKernel;
	private		int				igx,igy,igz;
    
	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=true;
	static final boolean		verbose=true;
    
	/**
	 *  constructor
     *  note: the membership function mems_ must have sum = 0 on the masked areas
	 */
	public DemonsScaledConsistentRegistration(float[] image_, float[] target_,
									int nix_, int niy_, int niz_,
									float rix_, float riy_, float riz_,
									int ntx_, int nty_, int ntz_,
									float rtx_, float rty_, float rtz_,
									float smoothing_,
									float spScale_,
									boolean useMask_,
									float maskVal_,
									float scale_, int Ni_, 
									int reg_, int force_, int field_,
									int pre_,
									float[][] trans_) {
	
		scale = scale_;
		
		smoothingKernel = smoothing_;
		spatialScale = spScale_;
		sigma2 = spatialScale*spatialScale;
		
		useMask = useMask_;
		maskingThreshold = maskVal_;
		Niter = Ni_;

		transform = trans_;
		//if (transform!=null) itransform = Numerics.invert3Dtransform(transform);
		if (transform!=null) itransform = Matrix3D.invertAffineTransform(transform);
		else itransform = null;
		
		nix = nix_;
		niy = niy_;
		niz = niz_;
		
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		ntx = ntx_;
		nty = nty_;
		ntz = ntz_;
		
		rtx = rtx_;
		rty = rty_;
		rtz = rtz_;
			
		nsx = Numerics.ceil(ntx/scale);
		nsy = Numerics.ceil(nty/scale);
		nsz = Numerics.ceil(ntz/scale);
		
		nisx = Numerics.ceil(nix/scale);
		nisy = Numerics.ceil(niy/scale);
		nisz = Numerics.ceil(niz/scale);
		
		regType = reg_;
		forceType = force_;
		fieldType = field_;
		preType = pre_;
		
		//initialize the smoothing kernel
		gaussKernel = ImageFilters.separableGaussianKernel(smoothingKernel/rtx,smoothingKernel/rty,smoothingKernel/rtz);
		gx = (gaussKernel[X].length-1)/2;
		gy = (gaussKernel[Y].length-1)/2;
		gz = (gaussKernel[Z].length-1)/2;
				
		igaussKernel = ImageFilters.separableGaussianKernel(smoothingKernel/rix,smoothingKernel/riy,smoothingKernel/riz);
		igx = (igaussKernel[X].length-1)/2;
		igy = (igaussKernel[Y].length-1)/2;
		igz = (igaussKernel[Z].length-1)/2;
				
		try {			
			preprocessInputImages(image_, target_);
		} catch (OutOfMemoryError e){
			isWorking = false;
            finalize();
			System.out.println(e.getMessage());
			return;
		}
		isWorking = true;
		if (debug) BasicInfo.displayMessage("Demons:initialisation\n");
	}

	final public void finalize() {
		image = null;
		target = null;
		s = null; u = null; c = null;
		is = null; iu = null; ic = null;
		System.gc();
	}
    
   public final float[][] getCurrentTransform() {
	   return s;
   }
   public final float[][] getCurrentUpdate() {
		return u;
	}
    
   public final float[][] getCurrentInverseTransform() {
	   return is;
   }
   public final float[][] getCurrentInverseUpdate() {
		return iu;
	}
    
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
	
	private final void preprocessInputImages(float[] img, float[] trg) {
		if (preType==RAW) {
			nc = 1;
			image = new float[nc][];
			target = new float[nc][];
			image[0] = img;
			target[0] = trg;
		} else if (preType==NORMALIZED) {
			// find min, max
			imin = ImageStatistics.robustMinimum(img, 0.01f, 2, nix, niy, niz);
			imax = ImageStatistics.robustMaximum(img, 0.01f, 2, nix, niy, niz);
			
			tmin = ImageStatistics.robustMinimum(trg, 0.01f, 2, ntx, nty, ntz);
			tmax = ImageStatistics.robustMaximum(trg, 0.01f, 2, ntx, nty, ntz);

			if (debug) BasicInfo.displayMessage("Image min / max: "+imin+" / "+imax+"\n");
			if (debug) BasicInfo.displayMessage("Target min / max: "+tmin+" / "+tmax+"\n");
			
			// allocate
			nc = 1;
			image = new float[nc+1][];
			target = new float[nc][];
			image[0] = ImageStatistics.normalize(img, nix,niy,niz, imin, imax);
			target[0] = ImageStatistics.normalize(trg, ntx,nty,ntz, tmin, tmax);
			// keep original image for results
			image[1] = img;
		} else if (preType==SCALED) {
			// find min, max
			imin = ImageStatistics.robustMinimum(img, 0.01f, 2, nix, niy, niz);
			imax = ImageStatistics.robustMaximum(img, 0.01f, 2, nix, niy, niz);
			
			tmin = ImageStatistics.robustMinimum(trg, 0.01f, 2, ntx, nty, ntz);
			tmax = ImageStatistics.robustMaximum(trg, 0.01f, 2, ntx, nty, ntz);
			
			if (debug) BasicInfo.displayMessage("Image min / max: "+imin+" / "+imax+"\n");
			if (debug) BasicInfo.displayMessage("Target min / max: "+tmin+" / "+tmax+"\n");
			
			ImageStatistics.rescale(img, nix,niy,niz, imin, imax);
			ImageStatistics.rescale(trg, ntx,nty,ntz, tmin, tmax);
			
			// allocate
			nc = 1;
			image = new float[nc][];
			target = new float[nc][];
			image[0] = img;
			target[0] = trg;
			
		} else if (preType==SEGMENTATION) {
			// find max number of objects and list them; assume the target has proper labeling
			lb = ObjectLabeling.listLabels(trg,ntx,nty,ntz);
			nc = lb.length;
			if (useMask) nc--;
			
			image = new float[nc][];
			target = new float[nc][];
			
			int c=0;
			for (int n=0;n<lb.length;n++) {
				if ( (!useMask) || (lb[n]!=maskingThreshold) ) {
					image[c] = ObjectExtraction.floatObjectFromImage(img,nix,niy,niz,lb[n],ObjectLabeling.EQUAL);
					target[c] = ObjectExtraction.floatObjectFromImage(trg,nix,niy,niz,lb[n],ObjectLabeling.EQUAL);
					c++;
				}			
			}
		} else if (preType==FCM) {
			// perform basic FCM segmentation
			nc = 3;
			//...
		} else if (preType==MULTICHANNEL) {
			nc = Numerics.round(img.length/(nix*niy*niz));
			image = new float[nc][nix*niy*niz];
			target = new float[nc][ntx*nty*ntz];
			for (int n=0;n<nc;n++) for (int xyz=0;xyz<nix*niy*niz;xyz++) {
				image[n][xyz] = img[xyz + n*nix*niy*niz];
			}
			for (int n=0;n<nc;n++) for (int xyz=0;xyz<ntx*nty*ntz;xyz++) {
				target[n][xyz] = trg[xyz + n*ntx*nty*ntz];
			}
			img = null;
			trg = null;
		} else {
			nc = 1;
			image = new float[nc][];
			target = new float[nc][];
			image[0] = img;
			target[0] = trg;
		}	
	}
	
    /** initialize the transform from previous estimate, if exists */
    public final void initializeTransform() {
				
		// initialization: start from prior transform or zero
		s = new float[3][nsx*nsy*nsz];
		if (transform==null) {
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				s[X][xyz] = x*scale*rtx/rix;
				s[Y][xyz] = y*scale*rty/riy;
				s[Z][xyz] = z*scale*rtz/riz;
			}
		} else {
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
    
    /** initialize the transform from previous estimate, if exists */
    public final void initializeInverseTransform() {
				
		// initialization: start from prior transform or zero
		is = new float[3][nisx*nisy*nisz];
		if (transform==null) {
			for (int x=0;x<nisx;x++) for (int y=0;y<nisy;y++) for (int z=0;z<nisz;z++) {
				int xyz = x+nisx*y+nisx*nisy*z;
				is[X][xyz] = x*scale*rix/rtx;
				is[Y][xyz] = y*scale*riy/rty;
				is[Z][xyz] = z*scale*riz/rtz;
			}
		} else {
			for (int x=0;x<nisx;x++) for (int y=0;y<nisy;y++) for (int z=0;z<nisz;z++) {
				int xyz = x+nisx*y+nisx*nisy*z;
				s[X][xyz] = itransform[X][X]*x*scale + itransform[X][Y]*y*scale + itransform[X][Z]*z*scale + itransform[X][T];
				s[Y][xyz] = itransform[Y][X]*x*scale + itransform[Y][Y]*y*scale + itransform[Y][Z]*z*scale + itransform[Y][T];
				s[Z][xyz] = itransform[Z][X]*x*scale + itransform[Z][Y]*y*scale + itransform[Z][Z]*z*scale + itransform[Z][T];
			}
		}
				
		// update always zero
		iu = new float[3][nisx*nisy*nisz];
		for (int xyz=0;xyz<nisx*nisy*nisz;xyz++) {
			iu[X][xyz] = 0.0f;
			iu[Y][xyz] = 0.0f;
			iu[Z][xyz] = 0.0f;
		}
		
		// composite update always zero
		ic = new float[3][nisx*nisy*nisz];
		for (int xyz=0;xyz<nisx*nisy*nisz;xyz++) {
			ic[X][xyz] = 0.0f;
			ic[Y][xyz] = 0.0f;
			ic[Z][xyz] = 0.0f;
		}
		
		return;
    }
    
	
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerImageToTargetUpdate() {
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
			
			float xsi = xs/scale, ysi = ys/scale, zsi = zs/scale;
			float xsmxi = xsmx/scale, ysmxi = ysmx/scale, zsmxi = zsmx/scale;
			float xspxi = xspx/scale, yspxi = yspx/scale, zspxi = zspx/scale;
			float xsmyi = xsmy/scale, ysmyi = ysmy/scale, zsmyi = zsmy/scale;
			float xspyi = xspy/scale, yspyi = yspy/scale, zspyi = zspy/scale;
			float xsmzi = xsmz/scale, ysmzi = ysmz/scale, zsmzi = zsmz/scale;
			float xspzi = xspz/scale, yspzi = yspz/scale, zspzi = zspz/scale;
			
			u[X][xyz] = 0.0f;
			u[Y][xyz] = 0.0f;
			u[Z][xyz] = 0.0f;
			float den = 0.0f;
			for (int n=0;n<nc;n++) {
				// image differences
				float diff = (ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt,ntx,nty,ntz) 
							- ImageInterpolation.linearClosestInterpolation(image[n],xs,ys,zs,nix,niy,niz));				
				
				// image gradient
				float Jx, Jy, Jz;
				if (forceType==FIXED) {
					Jx = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xt+scale,yt,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt-scale,yt,zt,ntx,nty,ntz));
					Jy = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt+scale,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt-scale,zt,ntx,nty,ntz));
					Jz = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt+scale,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt-scale,ntx,nty,ntz));
				} else if (forceType==MOVING) {
					Jx = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xspx,yspx,zspx,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmx,ysmx,zsmx,nix,niy,niz));
					Jy = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xspy,yspy,zspy,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmy,ysmy,zsmy,nix,niy,niz));
					Jz = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xspz,yspz,zspz,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmz,ysmz,zsmz,nix,niy,niz));
				} else {
					Jx = 0.25f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xt+scale,yt,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt-scale,yt,zt,ntx,nty,ntz));
					Jy = 0.25f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt+scale,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt-scale,zt,ntx,nty,ntz));
					Jz = 0.25f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt+scale,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt-scale,ntx,nty,ntz));

					Jx += 0.25f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xspx,yspx,zspx,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmx,ysmx,zsmx,nix,niy,niz));
					Jy += 0.25f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xspy,yspy,zspy,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmy,ysmy,zsmy,nix,niy,niz));
					Jz += 0.25f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xspz,yspz,zspz,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmz,ysmz,zsmz,nix,niy,niz));
				}
				
				// inverse consistent coupling: backprojected coordinates
				float iTx = xt - ImageInterpolation.linearClosestInterpolation(is[X],xsi,ysi,zsi,nisx,nisy,nisz);
				float iTy = yt - ImageInterpolation.linearClosestInterpolation(is[Y],xsi,ysi,zsi,nisx,nisy,nisz);
				float iTz = zt - ImageInterpolation.linearClosestInterpolation(is[Z],xsi,ysi,zsi,nisx,nisy,nisz);
				
				// erse consistent coupling: inverse transform jacobian
				float[][] diT = new float[3][3];
				diT[X][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[X],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[X][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[X],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[X][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[X],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				diT[Y][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[Y],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[Y][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[Y],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[Y][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[Y],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				diT[Z][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[Z],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[Z][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[Z],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[Z][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[Z],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				// putting it all together
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				
				float diT2 = Matrix3D.determinant(diT);
				diT2 *= diT2;				
				
				float sigmaI = diff*diff/sigma2;
				
				den += sigmaI*(1.0f+diT2) + J2;
				
				u[X][xyz] += diff*Jx;
				u[Y][xyz] += diff*Jy;
				u[Z][xyz] += diff*Jz;
				
				u[X][xyz] += sigmaI*( iTx*diT[X][X] + iTy*diT[X][Y] + iTz*diT[X][Z] );
				u[Y][xyz] += sigmaI*( iTx*diT[Y][X] + iTy*diT[Y][Y] + iTz*diT[Y][Z] );
				u[Z][xyz] += sigmaI*( iTx*diT[Z][X] + iTy*diT[Z][Y] + iTz*diT[Z][Z] );
				
				meanDiff += Numerics.abs(diff);
			}
			u[X][xyz] /= Numerics.max(ZERO,den);
			u[Y][xyz] /= Numerics.max(ZERO,den);
			u[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (ntx*nty*ntz);
		
		if (debug) BasicInfo.displayMessage("mean intensity difference "+meanDiff+"\n");

        return meanDiff;
    } // 

    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerImageToTargetUpdateDecoupled() {
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
			for (int n=0;n<nc;n++) {
				// image differences
				float diff = (ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt,ntx,nty,ntz) 
							- ImageInterpolation.linearClosestInterpolation(image[n],xs,ys,zs,nix,niy,niz));				
				
				// image gradient
				float Jx, Jy, Jz;
				if (forceType==FIXED) {
					Jx = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xt+scale,yt,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt-scale,yt,zt,ntx,nty,ntz));
					Jy = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt+scale,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt-scale,zt,ntx,nty,ntz));
					Jz = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt+scale,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt-scale,ntx,nty,ntz));
				} else if (forceType==MOVING) {
					Jx = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xspx,yspx,zspx,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmx,ysmx,zsmx,nix,niy,niz));
					Jy = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xspy,yspy,zspy,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmy,ysmy,zsmy,nix,niy,niz));
					Jz = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xspz,yspz,zspz,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmz,ysmz,zsmz,nix,niy,niz));
				} else {
					Jx = 0.25f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xt+scale,yt,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt-scale,yt,zt,ntx,nty,ntz));
					Jy = 0.25f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt+scale,zt,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt-scale,zt,ntx,nty,ntz));
					Jz = 0.25f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt+scale,ntx,nty,ntz) - ImageInterpolation.linearClosestInterpolation(target[n],xt,yt,zt-scale,ntx,nty,ntz));

					Jx += 0.25f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xspx,yspx,zspx,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmx,ysmx,zsmx,nix,niy,niz));
					Jy += 0.25f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xspy,yspy,zspy,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmy,ysmy,zsmy,nix,niy,niz));
					Jz += 0.25f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xspz,yspz,zspz,nix,niy,niz)-ImageInterpolation.linearClosestInterpolation(image[n],xsmz,ysmz,zsmz,nix,niy,niz));
				}
				
				// putting it all together
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				
				den += diff*diff/sigma2 + J2;
				
				u[X][xyz] += diff*Jx;
				u[Y][xyz] += diff*Jy;
				u[Z][xyz] += diff*Jz;
				
				meanDiff += Numerics.abs(diff);
			}
			u[X][xyz] /= Numerics.max(ZERO,den);
			u[Y][xyz] /= Numerics.max(ZERO,den);
			u[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (ntx*nty*ntz);
		
		if (debug) BasicInfo.displayMessage("mean intensity difference "+meanDiff+"\n");

        return meanDiff;
    } // 

    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerImageToTargetCompose() {		
		if (regType==GAUSS_FLUID || regType==GAUSS_MIXED) {
			//if (debug) BasicInfo.displayMessage("GAUSS_FLUID regularization \n");
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFilters.separableConvolution(u[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Y] = ImageFilters.separableConvolution(u[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Z] = ImageFilters.separableConvolution(u[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		}
		
		// compose the transformations
		if (fieldType==COMPOSITIVE) {
			//if (debug) BasicInfo.displayMessage("compose with current transform \n");
							
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				
				float xu = x+u[X][xyz];
				float yu = y+u[Y][xyz];
				float zu = z+u[Z][xyz];
				
				// note: if outside, extrapolate as X+u
				c[X][xyz] = ImageInterpolation.linearClosestInterpolation(s[X],xu,yu,zu,nsx,nsy,nsz) - x*scale;
				c[Y][xyz] = ImageInterpolation.linearClosestInterpolation(s[Y],xu,yu,zu,nsx,nsy,nsz) - y*scale;
				c[Z][xyz] = ImageInterpolation.linearClosestInterpolation(s[Z],xu,yu,zu,nsx,nsy,nsz) - z*scale;
			}
		}
		
		if (regType==GAUSS_DIFFUSION || regType==GAUSS_MIXED) {
			//if (debug) BasicInfo.displayMessage("GAUSS_DIFFUSION regularization \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFilters.separableConvolution(c[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Y] = ImageFilters.separableConvolution(c[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Z] = ImageFilters.separableConvolution(c[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		}
					
		float meanC = 0.0f;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x+nsx*y+nsx*nsy*z;
			
			s[X][xyz] = x*scale + c[X][xyz];
			s[Y][xyz] = y*scale + c[Y][xyz];
			s[Z][xyz] = z*scale + c[Z][xyz];
			
			meanC += c[X][xyz]*c[X][xyz]+c[Y][xyz]*c[Y][xyz]+c[Z][xyz]*c[Z][xyz];
		}
		meanC /= (nsx*nsy*nsz);
		
		//if (debug) BasicInfo.displayMessage("mean update size "+meanC+"\n");

        return meanC;
    } // 
    
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerTargetToImageUpdate() {
		float meanDiff = 0.0f;
		for (int x=1;x<nisx-1;x++) for (int y=1;y<nisy-1;y++) for (int z=1;z<nisz-1;z++) {
			int xyz = x+nisx*y+nisx*nisy*z;
		
			// compute the update field
			float xs = is[X][xyz];
			float ys = is[Y][xyz];
			float zs = is[Z][xyz];
		
			float xsmx = is[X][xyz-1];
			float ysmx = is[Y][xyz-1];
			float zsmx = is[Z][xyz-1];
	
			float xspx = is[X][xyz+1];
			float yspx = is[Y][xyz+1];
			float zspx = is[Z][xyz+1];
	
			float xsmy = is[X][xyz-nisx];
			float ysmy = is[Y][xyz-nisx];
			float zsmy = is[Z][xyz-nisx];
	
			float xspy = is[X][xyz+nisx];
			float yspy = is[Y][xyz+nisx];
			float zspy = is[Z][xyz+nisx];
	
			float xsmz = is[X][xyz-nisx*nisy];
			float ysmz = is[Y][xyz-nisx*nisy];
			float zsmz = is[Z][xyz-nisx*nisy];
	
			float xspz = is[X][xyz+nisx*nisy];
			float yspz = is[Y][xyz+nisx*nisy];
			float zspz = is[Z][xyz+nisx*nisy];
			
			float xsi = xs/scale, ysi = ys/scale, zsi = zs/scale;
			float xsmxi = xsmx/scale, ysmxi = ysmx/scale, zsmxi = zsmx/scale;
			float xspxi = xspx/scale, yspxi = yspx/scale, zspxi = zspx/scale;
			float xsmyi = xsmy/scale, ysmyi = ysmy/scale, zsmyi = zsmy/scale;
			float xspyi = xspy/scale, yspyi = yspy/scale, zspyi = zspy/scale;
			float xsmzi = xsmz/scale, ysmzi = ysmz/scale, zsmzi = zsmz/scale;
			float xspzi = xspz/scale, yspzi = yspz/scale, zspzi = zspz/scale;
			
			float xt = x*scale;
			float yt = y*scale;
			float zt = z*scale;
			
			iu[X][xyz] = 0.0f;
			iu[Y][xyz] = 0.0f;
			iu[Z][xyz] = 0.0f;
			float den = 0.0f;
			for (int n=0;n<nc;n++) {
				// image differences
				float diff = (ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt,nix,niy,niz) 
							- ImageInterpolation.linearClosestInterpolation(target[n],xs,ys,zs,ntx,nty,ntz));				
				
				// image gradient
				float Jx, Jy, Jz;
				if (forceType==FIXED) {
					Jx = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xt+scale,yt,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt-scale,yt,zt,nix,niy,niz));
					Jy = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt+scale,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt-scale,zt,nix,niy,niz));
					Jz = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt+scale,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt-scale,nix,niy,niz));
				} else if (forceType==MOVING) {
					Jx = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xspx,yspx,zspx,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmx,ysmx,zsmx,ntx,nty,ntz));
					Jy = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xspy,yspy,zspy,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmy,ysmy,zsmy,ntx,nty,ntz));
					Jz = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xspz,yspz,zspz,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmz,ysmz,zsmz,ntx,nty,ntz));
				} else {
					Jx = 0.25f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xt+scale,yt,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt-scale,yt,zt,nix,niy,niz));
					Jy = 0.25f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt+scale,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt-scale,zt,nix,niy,niz));
					Jz = 0.25f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt+scale,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt-scale,nix,niy,niz));
					
					Jx += 0.25f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xspx,yspx,zspx,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmx,ysmx,zsmx,ntx,nty,ntz));
					Jy += 0.25f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xspy,yspy,zspy,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmy,ysmy,zsmy,ntx,nty,ntz));
					Jz += 0.25f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xspz,yspz,zspz,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmz,ysmz,zsmz,ntx,nty,ntz));
				}
			
				// inverse consistent coupling: backprojected coordinates
				float iTx = xt - ImageInterpolation.linearClosestInterpolation(s[X],xsi,ysi,zsi,nsx,nsy,nsz);
				float iTy = yt - ImageInterpolation.linearClosestInterpolation(s[Y],xsi,ysi,zsi,nsx,nsy,nsz);
				float iTz = zt - ImageInterpolation.linearClosestInterpolation(s[Z],xsi,ysi,zsi,nsx,nsy,nsz);
				
				// erse consistent coupling: inverse transform jacobian
				float[][] diT = new float[3][3];
				diT[X][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[X],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[X][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[X],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[X][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[X],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));

				diT[Y][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[Y],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[Y][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[Y],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[Y][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[Y],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));

				diT[Z][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[Z],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[Z][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[Z],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[Z][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[Z],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));
				
				// putting all together
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				
				float diT2 = Matrix3D.determinant(diT);
				diT2 *= diT2;				
				
				float sigmaI = diff*diff/sigma2;
				
				den += sigmaI*(1.0f+diT2) + J2;
				
				iu[X][xyz] += diff*Jx;
				iu[Y][xyz] += diff*Jy;
				iu[Z][xyz] += diff*Jz;
				
				iu[X][xyz] += sigmaI*( iTx*diT[X][X] + iTy*diT[X][Y] + iTz*diT[X][Z] );
				iu[Y][xyz] += sigmaI*( iTx*diT[Y][X] + iTy*diT[Y][Y] + iTz*diT[Y][Z] );
				iu[Z][xyz] += sigmaI*( iTx*diT[Z][X] + iTy*diT[Z][Y] + iTz*diT[Z][Z] );
				
				meanDiff += Numerics.abs(diff);
			}
			iu[X][xyz] /= Numerics.max(ZERO,den);
			iu[Y][xyz] /= Numerics.max(ZERO,den);
			iu[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (nix*niy*niz);

		if (debug) BasicInfo.displayMessage("mean intensity difference "+meanDiff+"\n");

        return meanDiff;
    } // 
		
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerTargetToImageUpdateDecoupled() {
		
		float meanDiff = 0.0f;
		for (int x=1;x<nisx-1;x++) for (int y=1;y<nisy-1;y++) for (int z=1;z<nisz-1;z++) {
			int xyz = x+nisx*y+nisx*nisy*z;
		
			// compute the update field
			float xs = is[X][xyz];
			float ys = is[Y][xyz];
			float zs = is[Z][xyz];
		
			float xsmx = is[X][xyz-1];
			float ysmx = is[Y][xyz-1];
			float zsmx = is[Z][xyz-1];
	
			float xspx = is[X][xyz+1];
			float yspx = is[Y][xyz+1];
			float zspx = is[Z][xyz+1];
	
			float xsmy = is[X][xyz-nisx];
			float ysmy = is[Y][xyz-nisx];
			float zsmy = is[Z][xyz-nisx];
	
			float xspy = is[X][xyz+nisx];
			float yspy = is[Y][xyz+nisx];
			float zspy = is[Z][xyz+nisx];
	
			float xsmz = is[X][xyz-nisx*nisy];
			float ysmz = is[Y][xyz-nisx*nisy];
			float zsmz = is[Z][xyz-nisx*nisy];
	
			float xspz = is[X][xyz+nisx*nisy];
			float yspz = is[Y][xyz+nisx*nisy];
			float zspz = is[Z][xyz+nisx*nisy];
			
			float xt = x*scale;
			float yt = y*scale;
			float zt = z*scale;
			
			iu[X][xyz] = 0.0f;
			iu[Y][xyz] = 0.0f;
			iu[Z][xyz] = 0.0f;
			float den = 0.0f;
			for (int n=0;n<nc;n++) {
				// image differences
				float diff = (ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt,nix,niy,niz) 
							- ImageInterpolation.linearClosestInterpolation(target[n],xs,ys,zs,ntx,nty,ntz));				
				
				// image gradient
				float Jx, Jy, Jz;
				if (forceType==FIXED) {
					Jx = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xt+scale,yt,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt-scale,yt,zt,nix,niy,niz));
					Jy = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt+scale,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt-scale,zt,nix,niy,niz));
					Jz = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt+scale,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt-scale,nix,niy,niz));
				} else if (forceType==MOVING) {
					Jx = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xspx,yspx,zspx,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmx,ysmx,zsmx,ntx,nty,ntz));
					Jy = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xspy,yspy,zspy,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmy,ysmy,zsmy,ntx,nty,ntz));
					Jz = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xspz,yspz,zspz,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmz,ysmz,zsmz,ntx,nty,ntz));
				} else {
					Jx = 0.25f/rix*(ImageInterpolation.linearClosestInterpolation(image[n],xt+scale,yt,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt-scale,yt,zt,nix,niy,niz));
					Jy = 0.25f/riy*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt+scale,zt,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt-scale,zt,nix,niy,niz));
					Jz = 0.25f/riz*(ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt+scale,nix,niy,niz) - ImageInterpolation.linearClosestInterpolation(image[n],xt,yt,zt-scale,nix,niy,niz));
					
					Jx += 0.25f/rtx*(ImageInterpolation.linearClosestInterpolation(target[n],xspx,yspx,zspx,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmx,ysmx,zsmx,ntx,nty,ntz));
					Jy += 0.25f/rty*(ImageInterpolation.linearClosestInterpolation(target[n],xspy,yspy,zspy,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmy,ysmy,zsmy,ntx,nty,ntz));
					Jz += 0.25f/rtz*(ImageInterpolation.linearClosestInterpolation(target[n],xspz,yspz,zspz,ntx,nty,ntz)-ImageInterpolation.linearClosestInterpolation(target[n],xsmz,ysmz,zsmz,ntx,nty,ntz));
				}
			
				// putting all together
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				
				den += diff*diff/sigma2 + J2;
				
				iu[X][xyz] += diff*Jx;
				iu[Y][xyz] += diff*Jy;
				iu[Z][xyz] += diff*Jz;
				
				meanDiff += Numerics.abs(diff);
			}
			iu[X][xyz] /= Numerics.max(ZERO,den);
			iu[Y][xyz] /= Numerics.max(ZERO,den);
			iu[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (nix*niy*niz);

		if (debug) BasicInfo.displayMessage("mean intensity difference "+meanDiff+"\n");

        return meanDiff;
    } // 
		
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerTargetToImageCompose() {
		
		if (regType==GAUSS_FLUID || regType==GAUSS_MIXED) {
			//if (debug) BasicInfo.displayMessage("GAUSS_FLUID regularization \n");
		
			// smooth the result with a gaussian kernel
			iu[X] = ImageFilters.separableConvolution(iu[X],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
			iu[Y] = ImageFilters.separableConvolution(iu[Y],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
			iu[Z] = ImageFilters.separableConvolution(iu[Z],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
		}
		
		// compose the transformations
		if (fieldType==COMPOSITIVE) {
			//if (debug) BasicInfo.displayMessage("compose with current transform \n");
							
			for (int x=0;x<nisx;x++) for (int y=0;y<nisy;y++) for (int z=0;z<nisz;z++) {
				int xyz = x+nisx*y+nisx*nisy*z;
				
				float xu = x+iu[X][xyz];
				float yu = y+iu[Y][xyz];
				float zu = z+iu[Z][xyz];
				
				// note: if outside, extrapolate as X+u
				ic[X][xyz] = ImageInterpolation.linearClosestInterpolation(is[X],xu,yu,zu,nisx,nisy,nisz) - x*scale;
				ic[Y][xyz] = ImageInterpolation.linearClosestInterpolation(is[Y],xu,yu,zu,nisx,nisy,nisz) - y*scale;
				ic[Z][xyz] = ImageInterpolation.linearClosestInterpolation(is[Z],xu,yu,zu,nisx,nisy,nisz) - z*scale;
			}
		}
		
		if (regType==GAUSS_DIFFUSION || regType==GAUSS_MIXED) {
			//if (debug) BasicInfo.displayMessage("GAUSS_DIFFUSION regularization \n");
			
			// smooth the result with a gaussian kernel
			ic[X] = ImageFilters.separableConvolution(ic[X],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
			ic[Y] = ImageFilters.separableConvolution(ic[Y],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
			ic[Z] = ImageFilters.separableConvolution(ic[Z],nisx,nisy,nisz,igaussKernel,igx,igy,igz);
		}
					
		float meanC = 0.0f;
		for (int x=0;x<nisx;x++) for (int y=0;y<nisy;y++) for (int z=0;z<nisz;z++) {
			int xyz = x+nisx*y+nisx*nisy*z;
			
			is[X][xyz] = x*scale + ic[X][xyz];
			is[Y][xyz] = y*scale + ic[Y][xyz];
			is[Z][xyz] = z*scale + ic[Z][xyz];
			
			meanC += ic[X][xyz]*ic[X][xyz]+ic[Y][xyz]*ic[Y][xyz]+ic[Z][xyz]*ic[Z][xyz];
		}
		meanC /= (nisx*nisy*nisz);
		
		//if (debug) BasicInfo.displayMessage("mean update size "+meanC+"\n");

        return meanC;
    } // 
    
   /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerImageToTargetInverseConsistent() {
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
			
			float xsi = xs/scale, ysi = ys/scale, zsi = zs/scale;
			float xsmxi = xsmx/scale, ysmxi = ysmx/scale, zsmxi = zsmx/scale;
			float xspxi = xspx/scale, yspxi = yspx/scale, zspxi = zspx/scale;
			float xsmyi = xsmy/scale, ysmyi = ysmy/scale, zsmyi = zsmy/scale;
			float xspyi = xspy/scale, yspyi = yspy/scale, zspyi = zspy/scale;
			float xsmzi = xsmz/scale, ysmzi = ysmz/scale, zsmzi = zsmz/scale;
			float xspzi = xspz/scale, yspzi = yspz/scale, zspzi = zspz/scale;
			
			u[X][xyz] = 0.0f;
			u[Y][xyz] = 0.0f;
			u[Z][xyz] = 0.0f;
			float den = 0.0f;
			for (int n=0;n<nc;n++) {
				// inverse consistent coupling: backprojected coordinates
				float iTx = xt - ImageInterpolation.linearClosestInterpolation(is[X],xsi,ysi,zsi,nisx,nisy,nisz);
				float iTy = yt - ImageInterpolation.linearClosestInterpolation(is[Y],xsi,ysi,zsi,nisx,nisy,nisz);
				float iTz = zt - ImageInterpolation.linearClosestInterpolation(is[Z],xsi,ysi,zsi,nisx,nisy,nisz);
				
				// inverse consistent coupling: inverse transform jacobian
				float[][] diT = new float[3][3];
				diT[X][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[X],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[X][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[X],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[X][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[X],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[X],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				diT[Y][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[Y],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[Y][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[Y],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[Y][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[Y],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Y],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				diT[Z][X] = 0.5f/rix*(ImageInterpolation.linearClosestInterpolation(is[Z],xspxi,yspxi,zspxi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmxi,ysmxi,zsmxi,nisx,nisy,nisz));
				diT[Z][Y] = 0.5f/riy*(ImageInterpolation.linearClosestInterpolation(is[Z],xspyi,yspyi,zspyi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmyi,ysmyi,zsmyi,nisx,nisy,nisz));
				diT[Z][Z] = 0.5f/riz*(ImageInterpolation.linearClosestInterpolation(is[Z],xspzi,yspzi,zspzi,nisx,nisy,nisz)-ImageInterpolation.linearClosestInterpolation(is[Z],xsmzi,ysmzi,zsmzi,nisx,nisy,nisz));

				// putting it all together
				float diT2 = Matrix3D.determinant(diT);
				diT2 *= diT2;				
				
				float iT2 = iTx*iTx + iTy*iTy + iTz*iTz;
				
				den += iT2 + diT2;
				
				u[X][xyz] += 0.5f*( iTx*diT[X][X] + iTy*diT[X][Y] + iTz*diT[X][Z] );
				u[Y][xyz] += 0.5f*( iTx*diT[Y][X] + iTy*diT[Y][Y] + iTz*diT[Y][Z] );
				u[Z][xyz] += 0.5f*( iTx*diT[Z][X] + iTy*diT[Z][Y] + iTz*diT[Z][Z] );
				
				meanDiff += (float)Math.sqrt(iT2);
			}
			u[X][xyz] /= Numerics.max(ZERO,den);
			u[Y][xyz] /= Numerics.max(ZERO,den);
			u[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (ntx*nty*ntz);

		if (debug) BasicInfo.displayMessage("inverse consistency: "+meanDiff+"\n");

        return meanDiff;
    } // 
    
    /**
	 * compute the image position given the target
	 * for a given level l
     * performs only one iteration
	 */
    final public float registerTargetToImageInverseConsistent() {
		float meanDiff = 0.0f;
		for (int x=1;x<nisx-1;x++) for (int y=1;y<nisy-1;y++) for (int z=1;z<nisz-1;z++) {
			int xyz = x+nisx*y+nisx*nisy*z;
		
			// compute the update field
			float xs = is[X][xyz];
			float ys = is[Y][xyz];
			float zs = is[Z][xyz];
		
			float xsmx = is[X][xyz-1];
			float ysmx = is[Y][xyz-1];
			float zsmx = is[Z][xyz-1];
	
			float xspx = is[X][xyz+1];
			float yspx = is[Y][xyz+1];
			float zspx = is[Z][xyz+1];
	
			float xsmy = is[X][xyz-nisx];
			float ysmy = is[Y][xyz-nisx];
			float zsmy = is[Z][xyz-nisx];
	
			float xspy = is[X][xyz+nisx];
			float yspy = is[Y][xyz+nisx];
			float zspy = is[Z][xyz+nisx];
	
			float xsmz = is[X][xyz-nisx*nisy];
			float ysmz = is[Y][xyz-nisx*nisy];
			float zsmz = is[Z][xyz-nisx*nisy];
	
			float xspz = is[X][xyz+nisx*nisy];
			float yspz = is[Y][xyz+nisx*nisy];
			float zspz = is[Z][xyz+nisx*nisy];
			
			float xsi = xs/scale, ysi = ys/scale, zsi = zs/scale;
			float xsmxi = xsmx/scale, ysmxi = ysmx/scale, zsmxi = zsmx/scale;
			float xspxi = xspx/scale, yspxi = yspx/scale, zspxi = zspx/scale;
			float xsmyi = xsmy/scale, ysmyi = ysmy/scale, zsmyi = zsmy/scale;
			float xspyi = xspy/scale, yspyi = yspy/scale, zspyi = zspy/scale;
			float xsmzi = xsmz/scale, ysmzi = ysmz/scale, zsmzi = zsmz/scale;
			float xspzi = xspz/scale, yspzi = yspz/scale, zspzi = zspz/scale;
			
			float xt = x*scale;
			float yt = y*scale;
			float zt = z*scale;
			
			iu[X][xyz] = 0.0f;
			iu[Y][xyz] = 0.0f;
			iu[Z][xyz] = 0.0f;
			float den = 0.0f;
			for (int n=0;n<nc;n++) {
				// inverse consistent coupling: backprojected coordinates
				float iTx = xt - ImageInterpolation.linearClosestInterpolation(s[X],xsi,ysi,zsi,nsx,nsy,nsz);
				float iTy = yt - ImageInterpolation.linearClosestInterpolation(s[Y],xsi,ysi,zsi,nsx,nsy,nsz);
				float iTz = zt - ImageInterpolation.linearClosestInterpolation(s[Z],xsi,ysi,zsi,nsx,nsy,nsz);
				
				// erse consistent coupling: inverse transform jacobian
				float[][] diT = new float[3][3];
				diT[X][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[X],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[X][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[X],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[X][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[X],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[X],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));

				diT[Y][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[Y],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[Y][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[Y],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[Y][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[Y],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Y],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));

				diT[Z][X] = 0.5f/rtx*(ImageInterpolation.linearClosestInterpolation(s[Z],xspxi,yspxi,zspxi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmxi,ysmxi,zsmxi,nsx,nsy,nsz));
				diT[Z][Y] = 0.5f/rty*(ImageInterpolation.linearClosestInterpolation(s[Z],xspyi,yspyi,zspyi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmyi,ysmyi,zsmyi,nsx,nsy,nsz));
				diT[Z][Z] = 0.5f/rtz*(ImageInterpolation.linearClosestInterpolation(s[Z],xspzi,yspzi,zspzi,nsx,nsy,nsz)-ImageInterpolation.linearClosestInterpolation(s[Z],xsmzi,ysmzi,zsmzi,nsx,nsy,nsz));
				
				// putting all together
				float diT2 = Matrix3D.determinant(diT);
				diT2 *= diT2;				
				
				float iT2 = iTx*iTx + iTy*iTy + iTz*iTz;
				
				den += iT2 + diT2;
				
				iu[X][xyz] += 0.5f*( iTx*diT[X][X] + iTy*diT[X][Y] + iTz*diT[X][Z] );
				iu[Y][xyz] += 0.5f*( iTx*diT[Y][X] + iTy*diT[Y][Y] + iTz*diT[Y][Z] );
				iu[Z][xyz] += 0.5f*( iTx*diT[Z][X] + iTy*diT[Z][Y] + iTz*diT[Z][Z] );
				
				meanDiff += (float)Math.sqrt(iT2);
			}
			iu[X][xyz] /= Numerics.max(ZERO,den);
			iu[Y][xyz] /= Numerics.max(ZERO,den);
			iu[Z][xyz] /= Numerics.max(ZERO,den);			
		}
		meanDiff /= (nix*niy*niz);

		if (debug) BasicInfo.displayMessage("inverse consistency "+meanDiff+"\n");

        return meanDiff;
    } // 
    
    /** 
	 *	runs the algorithm
	 */
	public final void runRegistration() {        
        // top level
        if (debug) BasicInfo.displayMessage("initialize all parameters \n");
        initializeTransform();
        initializeInverseTransform();
		
		float dist=0, dist0=0, prec=0;
		for (int t=0;t<Niter;t++) {
			if (debug) BasicInfo.displayMessage("iteration "+(t+1)+"\n");
			prec = dist;
			
			dist = registerImageToTargetUpdateDecoupled();
			dist += registerTargetToImageUpdateDecoupled();
			
			if (t==0) {
				dist0 = dist;
			} else {
				BasicInfo.displayMessage("score: "+(dist/dist0)+", variation: "+( (dist-prec)/dist0 )+"\n");
			}
			registerImageToTargetCompose();
			registerTargetToImageCompose();
			
			registerImageToTargetInverseConsistent();
			registerTargetToImageInverseConsistent();

			registerImageToTargetCompose();
			registerTargetToImageCompose();
		}
    }   
	
	/** 
	 *	returns the transformed coordinates
	 */
	public final float[] getCurrentMapping(int x, int y, int z) {
		float xs = x/scale;
		float ys = y/scale;
		float zs = z/scale;
		
		float[] Xs = new float[3];
		Xs[X] = ImageInterpolation.linearClosestInterpolation(s[X],xs,ys,zs,nsx,nsy,nsz);
		Xs[Y] = ImageInterpolation.linearClosestInterpolation(s[Y],xs,ys,zs,nsx,nsy,nsz);
		Xs[Z] = ImageInterpolation.linearClosestInterpolation(s[Z],xs,ys,zs,nsx,nsy,nsz);
		return Xs;
	}// getCurrentMapping

	/** 
	 *	returns the transformed coordinates
	 */
	public final float[] getCurrentInverseMapping(int x, int y, int z) {
		float xs = x/scale;
		float ys = y/scale;
		float zs = z/scale;
		
		float[] Xs = new float[3];
		Xs[X] = ImageInterpolation.linearClosestInterpolation(is[X],xs,ys,zs,nisx,nisy,nisz);
		Xs[Y] = ImageInterpolation.linearClosestInterpolation(is[Y],xs,ys,zs,nisx,nisy,nisz);
		Xs[Z] = ImageInterpolation.linearClosestInterpolation(is[Z],xs,ys,zs,nisx,nisy,nisz);
		return Xs;
	}// getCurrentMapping

    /** 
	 *	returns the transformed image
	 */
	public final float[][] exportTransformedImage() {
		float 	xs,ys,zs;
		float[][]	img;
		
		if (preType==MULTICHANNEL) img = new float[nc][ntx*nty*ntz];
		else img = new float[1][ntx*nty*ntz];
		
		if (preType==SCALED) ImageStatistics.unscale(image[0],nix,niy,niz,imin,imax);
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
				// compute interpolated values
			if (preType==RAW) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(image[0],xs,ys,zs,nix,niy,niz);
			} else if (preType==NORMALIZED) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(image[1],xs,ys,zs,nix,niy,niz);
			} else if (preType==SCALED) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(image[0],xs,ys,zs,nix,niy,niz);
			} else if (preType==SEGMENTATION) {
				float nbest = maskingThreshold;
				float best = 0.0f;
				int c = 0;
				for (int n=0;n<lb.length;n++) if ( (!useMask) || (lb[n]!=maskingThreshold) ) {
					float score = ImageInterpolation.linearClosestInterpolation(image[c],xs,ys,zs,nix,niy,niz);
					if (score>best) {
						best = score;
						nbest = lb[n];
					}
					c++;
				}
				img[0][xyz] = nbest;
			} else if (preType==MULTICHANNEL) {
				for (int n=0;n<nc;n++) {
					img[n][xyz] = ImageInterpolation.linearClosestInterpolation(image[n],xs,ys,zs,nix,niy,niz);
				}
			}
		}
		return img;
	} // exportTransformedImage

     /** 
	 *	returns the transformed image
	 */
	public final float[][] exportInverseTransformedTarget() {
		float 	xs,ys,zs;
		float[][]	img;
		
		if (preType==MULTICHANNEL) img = new float[nc][nix*niy*niz];
		else img = new float[1][nix*niy*niz];
		
		if (preType==SCALED) ImageStatistics.unscale(target[0],ntx,nty,ntz,tmin,tmax);
		
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
			xs = ImageInterpolation.linearClosestInterpolation(is[X],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			ys = ImageInterpolation.linearClosestInterpolation(is[Y],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			zs = ImageInterpolation.linearClosestInterpolation(is[Z],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			
				// compute interpolated values
			if (preType==RAW) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(target[0],xs,ys,zs,ntx,nty,ntz);
			} else if (preType==NORMALIZED) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(target[1],xs,ys,zs,ntx,nty,ntz);
			} else if (preType==SCALED) {
				img[0][xyz] = ImageInterpolation.linearClosestInterpolation(target[0],xs,ys,zs,ntx,nty,ntz);
			} else if (preType==SEGMENTATION) {
				float nbest = maskingThreshold;
				float best = 0.0f;
				int c = 0;
				for (int n=0;n<lb.length;n++) if ( (!useMask) || (lb[n]!=maskingThreshold) ) {
					float score = ImageInterpolation.linearClosestInterpolation(target[c],xs,ys,zs,ntx,nty,ntz);
					if (score>best) {
						best = score;
						nbest = lb[n];
					}
					c++;
				}
				img[0][xyz] = nbest;
			} else if (preType==MULTICHANNEL) {
				for (int n=0;n<nc;n++) {
					img[n][xyz] = ImageInterpolation.linearClosestInterpolation(target[n],xs,ys,zs,ntx,nty,ntz);
				}
			}
		}
		return img;
	} // exportTransformedImage

    /** 
	 *	returns the transform field v defined as X' = X+v (=s)
	 */
	public final float[][] exportTransformMapping() {
		float 	xs,ys,zs;
		float[][]	vec = new float[3][ntx*nty*ntz];
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
			
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			vec[X][xyz] = xs;
			vec[Y][xyz] = ys;
			vec[Z][xyz] = zs;			
 		}
		return vec;
	} // exportTransformedImage

    /** 
	 *	returns the transform field v defined as X' = X+v (=s)
	 */
	public final float[][] exportInverseTransformMapping() {
		float 	xs,ys,zs;
		float[][]	vec = new float[3][nix*niy*niz];
		
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
			
			xs = ImageInterpolation.linearClosestInterpolation(is[X],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			ys = ImageInterpolation.linearClosestInterpolation(is[Y],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			zs = ImageInterpolation.linearClosestInterpolation(is[Z],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			
			vec[X][xyz] = xs;
			vec[Y][xyz] = ys;
			vec[Z][xyz] = zs;			
 		}
		return vec;
	} // exportTransformedImage

	/** 
	 *	returns the image residuals
	 */
	public final float[][] exportResiduals() {
		float 	xs,ys,zs;
		float[][]	img;
		
		if (preType==MULTICHANNEL) img = new float[nc][ntx*nty*ntz];
		else img = new float[1][ntx*nty*ntz];
		
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
    		xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
				// compute interpolated values
			if (preType==RAW || preType==NORMALIZED || preType==SCALED) {
				img[0][xyz] = Numerics.abs(target[0][xyz] - ImageInterpolation.linearClosestInterpolation(image[0],xs,ys,zs,nix,niy,niz));
			} else if (preType==SEGMENTATION) {
				img[0][xyz] = 0.0f;
				for (int n=0;n<nc;n++) {
					img[0][xyz] += Numerics.square(target[n][xyz]-ImageInterpolation.linearClosestInterpolation(image[n],xs,ys,zs,nix,niy,niz));
				}
				img[0][xyz] = (float)Math.sqrt(img[0][xyz]);
			} else if (preType==MULTICHANNEL) {
				for (int n=0;n<nc;n++) {
					img[n][xyz] = Numerics.abs(target[n][xyz]-ImageInterpolation.linearClosestInterpolation(image[n],xs,ys,zs,nix,niy,niz));
				}
			}
		}
		return img;
	} // exportResiduals

	/** 
	 *	returns the image residuals
	 */
	public final float[][] exportInverseResiduals() {
		float 	xs,ys,zs;
		float[][]	img;
		
		if (preType==MULTICHANNEL) img = new float[nc][nix*niy*niz];
		else img = new float[1][nix*niy*niz];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
    		xs = ImageInterpolation.linearClosestInterpolation(is[X],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			ys = ImageInterpolation.linearClosestInterpolation(is[Y],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			zs = ImageInterpolation.linearClosestInterpolation(is[Z],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			
				// compute interpolated values
			if (preType==RAW || preType==NORMALIZED || preType==SCALED) {
				img[0][xyz] = Numerics.abs(image[0][xyz] - ImageInterpolation.linearClosestInterpolation(target[0],xs,ys,zs,ntx,nty,ntz));
			} else if (preType==SEGMENTATION) {
				img[0][xyz] = 0.0f;
				for (int n=0;n<nc;n++) {
					img[0][xyz] += Numerics.square(image[n][xyz]-ImageInterpolation.linearClosestInterpolation(target[n],xs,ys,zs,ntx,nty,ntz));
				}
				img[0][xyz] = (float)Math.sqrt(img[0][xyz]);
			} else if (preType==MULTICHANNEL) {
				for (int n=0;n<nc;n++) {
					img[n][xyz] = Numerics.abs(image[n][xyz]-ImageInterpolation.linearClosestInterpolation(target[n],xs,ys,zs,ntx,nty,ntz));
				}
			}
		}
		return img;
	} // exportResiduals

	/** 
	 *	returns the transform residuals
	 */
	public final float[] exportForwardBackwardDistance() {
		float 	xs,ys,zs;
		float 	xsi,ysi,zsi;
		float[]	img;
		
		img = new float[ntx*nty*ntz];
		
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
    		xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			xsi = ImageInterpolation.linearClosestInterpolation(is[X],xs/scale,ys/scale,zs/scale,nisx,nisy,nisz);
			ysi = ImageInterpolation.linearClosestInterpolation(is[Y],xs/scale,ys/scale,zs/scale,nisx,nisy,nisz);
			zsi = ImageInterpolation.linearClosestInterpolation(is[Z],xs/scale,ys/scale,zs/scale,nisx,nisy,nisz);
			
			img[xyz] = (float)Math.sqrt( Numerics.square(x-xsi) + Numerics.square(y-ysi) + Numerics.square(z-zsi) );
		}
		return img;
	} // exportResiduals

	
	/** 
	 *	returns the transform residuals
	 */
	public final float[] exportBackwardForwardDistance() {
		float 	xs,ys,zs;
		float 	xsi,ysi,zsi;
		float[]	img;
		
		img = new float[nix*niy*niz];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			int xyz = x+nix*y+nix*niy*z;
    		xs = ImageInterpolation.linearClosestInterpolation(is[X],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			ys = ImageInterpolation.linearClosestInterpolation(is[Y],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			zs = ImageInterpolation.linearClosestInterpolation(is[Z],x/scale,y/scale,z/scale,nisx,nisy,nisz);
			
			xsi = ImageInterpolation.linearClosestInterpolation(s[X],xs/scale,ys/scale,zs/scale,nsx,nsy,nsz);
			ysi = ImageInterpolation.linearClosestInterpolation(s[Y],xs/scale,ys/scale,zs/scale,nsx,nsy,nsz);
			zsi = ImageInterpolation.linearClosestInterpolation(s[Z],xs/scale,ys/scale,zs/scale,nsx,nsy,nsz);
			
			img[xyz] = (float)Math.sqrt( Numerics.square(x-xsi) + Numerics.square(y-ysi) + Numerics.square(z-zsi) );
		}
		return img;
	} // exportResiduals


	
}
