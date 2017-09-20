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
 *	@version    November 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class ShapeAtlasWarping {
		
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
	public static final	int   DIV_FLUID 		= 304;
	public static final	int   DIV_DIFFUSION 	= 305;
	public static final	int   GAUSS_DIV_FLUID 		= 306;
	public static final	int   GAUSS_DIV_DIFFUSION 	= 307;
	
	// data buffers
	private 	float[][]		atlas;  		// source image
	private		int				na;				// number of atlas images
	private 	float[]		    famap;			// target image: FA map
	private 	float[][]	    gain;			// target image: best memberships
	private 	short[][]	   	label;			// target image: best labels
	private		int				nb;				// number of best targets
	private		boolean[]		registered;		// labels used for registration
	private		float[][]		c;				// added transform from target space to original space (multiscale)
	private		float[][]		s;				// added transform from target space to original space (multiscale)
	private		float[][]		u;				// added transform update from target space to original space (multiscale)
	private static	int		nax,nay,naz;   	// image dimensions (pyramid)
	private static	int		ntx,nty,ntz;   	// target dimensions (pyramid)
	private static	int		nsx,nsy,nsz;   	// target dimensions (pyramid)
	private static	float 	scale;
	private static	float	rax,ray,raz;   	// image resolutions (no pyramid)
	private static	float	rtx,rty,rtz;   	// target resolutions (no pyramid)
	private 		float[][]	transform;		// prior transform matrax (for instance from prior alignment)		
	
	// parameters
	private		float		smoothingKernel; 	// smoothing kernel size
	private		float		spatialScale; 		// scale of transformation
	private		boolean		useMask;
	private		float  		maskingThreshold;			// background threshold
	private		int		regType;
	private		int		forceType;
	private		int		fieldType;
	
	private		float		atlasScale = 5.0f;
	
	// computation variables
	private		float			sigma2; 		// scale of transformation
	

	private		float[][]		gaussKernel;
	private		int				gx,gy,gz;
    
	private		float[][]		divKernelX, divKernelY, divKernelZ;

	// computation flags
	private 	boolean 		isWorking;
	private 	boolean 		isCompleted;
	
	// for debug and display
	static final boolean		debug=true;
	static final boolean		verbose=true;
    
	/**
	 *  constructor
     *  note: the membership function gain_ must have sum = 0 on the masked areas
	 */
	public ShapeAtlasWarping(float[][] img_atlas_, int na_, 
									float[] trg_fa_, float[][] trg_gain_, short[][] trg_lbgain_, int nb_,
									boolean[] reged_,
									int nax_, int nay_, int naz_,
									float rax_, float ray_, float raz_,
									int ntx_, int nty_, int ntz_,
									float rtx_, float rty_, float rtz_,
									float smoothing_,
									float spScale_,
									boolean useMask_,
									float maskVal_,
									float scale_, int Ni_, int Nt_,
									int reg_, int force_, int field_,
									float[][] trans_) {
	
		atlas = img_atlas_;
		na = na_;
		famap = trg_fa_;
		gain = trg_gain_;
		label = trg_lbgain_;
		nb = nb_;
		registered = reged_;
		
		scale = scale_;
		
		smoothingKernel = smoothing_;
		spatialScale = spScale_;
		sigma2 = spatialScale*spatialScale;
		
		useMask = useMask_;
		maskingThreshold = maskVal_;

		transform = trans_;
		
		nax = nax_;
		nay = nay_;
		naz = naz_;
		
		rax = rax_;
		ray = ray_;
		raz = raz_;
		
		ntx = ntx_;
		nty = nty_;
		ntz = ntz_;
		
		rtx = rtx_;
		rty = rty_;
		rtz = rtz_;
			
		nsx = Numerics.ceil(ntx/scale);
		nsy = Numerics.ceil(nty/scale);
		nsz = Numerics.ceil(ntz/scale);
		
		regType = reg_;
		forceType = force_;
		fieldType = field_;

		//initialize the smoothing kernel
		gaussKernel = ImageFilters.separableGaussianKernel(smoothingKernel/rtx,smoothingKernel/rty,smoothingKernel/rtz);
		gx = (gaussKernel[X].length-1)/2;
		gy = (gaussKernel[Y].length-1)/2;
		gz = (gaussKernel[Z].length-1)/2;
		
		divKernelX = new float[3][];
		divKernelX[X] = new float[3];
		divKernelX[Y] = new float[1];
		divKernelX[Z] = new float[1];
		divKernelX[X][0] = 1.0f/4.0f;
		divKernelX[X][1] = 1.0f/2.0f;
		divKernelX[X][2] = 1.0f/4.0f;
		divKernelX[Y][0] = 1.0f;
		divKernelX[Z][0] = 1.0f;
				
		divKernelY = new float[3][];
		divKernelY[X] = new float[1];
		divKernelY[Y] = new float[3];
		divKernelY[Z] = new float[1];
		divKernelY[X][0] = 1.0f;
		divKernelY[Y][0] = 1.0f/4.0f;
		divKernelY[Y][1] = 1.0f/2.0f;
		divKernelY[Y][2] = 1.0f/4.0f;
		divKernelY[Z][0] = 1.0f;
				
		divKernelZ = new float[3][];
		divKernelZ[X] = new float[1];
		divKernelZ[Y] = new float[1];
		divKernelZ[Z] = new float[3];
		divKernelZ[X][0] = 1.0f;
		divKernelZ[Y][0] = 1.0f;
		divKernelZ[Z][0] = 1.0f/4.0f;
		divKernelZ[Z][1] = 1.0f/2.0f;
		divKernelZ[Z][2] = 1.0f/4.0f;

		
		isWorking = true;
		if (debug) BasicInfo.displayMessage("Demons:initialisation\n");
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
    
	public final boolean isWorking() { return isWorking; }
	public final boolean isCompleted() { return isCompleted; }
    
	
    /** initialize the transform from previous estimate, if exists */
    public final void initializeTransform() {
				
		// initialization: start from prior transform or zero
		s = new float[3][nsx*nsy*nsz];
		if (transform==null) {
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				s[X][xyz] = x*scale*rtx/rax;
				s[Y][xyz] = y*scale*rty/ray;
				s[Z][xyz] = z*scale*rtz/raz;
			}
		} else {
			/*
			MedicUtil.displayMessage("transform matrix: \n ["+transform[X][X]+", "+transform[X][Y]+", "+transform[X][Z]+", "+transform[X][T]+"]\n ["
															   +transform[Y][X]+", "+transform[Y][Y]+", "+transform[Y][Z]+", "+transform[Y][T]+"]\n ["
															   +transform[Z][X]+", "+transform[Z][Y]+", "+transform[Z][Z]+", "+transform[Z][T]+"]\n");
			*/
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
				
			// case by case
			if (famap!=null) {
				float diff = 0.0f;
				float diffpx = 0.0f; float diffpy = 0.0f; float diffpz = 0.0f;
				float diffmx = 0.0f; float diffmy = 0.0f; float diffmz = 0.0f;
				for (int n=0;n<na;n++) {
					diff = Numerics.max(diff,ImageInterpolation.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz));
					diffpx = Numerics.max(diffpx,ImageInterpolation.linearClosestInterpolation(atlas[n],xspx,yspx,zspx,nax,nay,naz));
					diffpy = Numerics.max(diffpy,ImageInterpolation.linearClosestInterpolation(atlas[n],xspy,yspy,zspy,nax,nay,naz));
					diffpz = Numerics.max(diffpz,ImageInterpolation.linearClosestInterpolation(atlas[n],xspz,yspz,zspz,nax,nay,naz));
					diffmx = Numerics.max(diffmx,ImageInterpolation.linearClosestInterpolation(atlas[n],xsmx,ysmx,zsmx,nax,nay,naz));
					diffmy = Numerics.max(diffmy,ImageInterpolation.linearClosestInterpolation(atlas[n],xsmy,ysmy,zsmy,nax,nay,naz));
					diffmz = Numerics.max(diffmz,ImageInterpolation.linearClosestInterpolation(atlas[n],xsmz,ysmz,zsmz,nax,nay,naz));
				}
				diff = (ImageInterpolation.linearClosestInterpolation(famap,xt,yt,zt,ntx,nty,ntz) - diff);
				
				float Jx, Jy, Jz;
				// symmetric ? no
				Jx = 0.5f/rax*(diffpx-diffmx);
				Jy = 0.5f/ray*(diffpy-diffmy);
				Jz = 0.5f/raz*(diffpz-diffmz);
			
				float J2 = Jx*Jx+Jy*Jy+Jz*Jz;
				den = Numerics.max(ZERO, diff*diff/sigma2 + J2);
				
				u[X][xyz] += diff*Jx/den;
				u[Y][xyz] += diff*Jy/den;
				u[Z][xyz] += diff*Jz/den;
			} else {
				float mem, diff;
				float Jx, Jy, Jz, J2;
				
				float sum = 0.0f;
				for (int n=0;n<nb;n++) {
					if (ImageInterpolation.linearClosestInterpolation(gain[n],xt,yt,zt,ntx,nty,ntz)>0) sum += ImageInterpolation.linearClosestInterpolation(gain[n],xt,yt,zt,ntx,nty,ntz);
				}
				sum = Numerics.max(ZERO,sum);
				for (int n=0;n<nb;n++) {
					short nbest = ImageInterpolation.nearestNeighborInterpolation(label[n],(short)0,xt,yt,zt,ntx,nty,ntz);
					
					if (nbest==-1) {
						// stop the loop
						n = nb;
					} else {
						if (nbest<100) {
							if (registered[nbest]) {
								mem = 0.5f + 0.5f*(2.0f*ImageInterpolation.linearClosestInterpolation(gain[n],xt,yt,zt,ntx,nty,ntz) - sum)/sum;
								diff = (mem - mappedAtlasInterpolation(nbest,xs,ys,zs));
						
								Jx = 0.5f/rax*(mappedAtlasInterpolation(nbest,xspx,yspx,zspx)-mappedAtlasInterpolation(nbest,xsmx,ysmx,zsmx));
								Jy = 0.5f/ray*(mappedAtlasInterpolation(nbest,xspy,yspy,zspy)-mappedAtlasInterpolation(nbest,xsmy,ysmy,zsmy));
								Jz = 0.5f/raz*(mappedAtlasInterpolation(nbest,xspz,yspz,zspz)-mappedAtlasInterpolation(nbest,xsmz,ysmz,zsmz));
								 
								J2 = Jx*Jx+Jy*Jy+Jz*Jz;
								den += diff*diff/sigma2 + J2;
								
								u[X][xyz] += diff*Jx;
								u[Y][xyz] += diff*Jy;
								u[Z][xyz] += diff*Jz;
									
								meanDiff += Numerics.abs(diff);
							}
						} else {
							int nbest2 = nbest%100;
							int nbest1 = (nbest-nbest2)/100;
							
							mem = 0.5f + 0.5f*(2.0f*ImageInterpolation.linearClosestInterpolation(gain[n],xt,yt,zt,ntx,nty,ntz) - sum)/sum;
								
							// first tract
							if (registered[nbest1]) {
								diff = (mem - mappedAtlasInterpolation(nbest1,xs,ys,zs));
						
								Jx = 0.5f/rax*(mappedAtlasInterpolation(nbest1,xspx,yspx,zspx)-mappedAtlasInterpolation(nbest1,xsmx,ysmx,zsmx));
								Jy = 0.5f/ray*(mappedAtlasInterpolation(nbest1,xspy,yspy,zspy)-mappedAtlasInterpolation(nbest1,xsmy,ysmy,zsmy));
								Jz = 0.5f/raz*(mappedAtlasInterpolation(nbest1,xspz,yspz,zspz)-mappedAtlasInterpolation(nbest1,xsmz,ysmz,zsmz));
								 
								J2 = Jx*Jx+Jy*Jy+Jz*Jz;
								den += diff*diff/sigma2 + J2;
								
								u[X][xyz] += diff*Jx;
								u[Y][xyz] += diff*Jy;
								u[Z][xyz] += diff*Jz;
									
								meanDiff += Numerics.abs(diff);
							}
							// second tract
							if (registered[nbest2]) {
								diff = (mem - mappedAtlasInterpolation(nbest2,xs,ys,zs));
						
								Jx = 0.5f/rax*(mappedAtlasInterpolation(nbest2,xspx,yspx,zspx)-mappedAtlasInterpolation(nbest2,xsmx,ysmx,zsmx));
								Jy = 0.5f/ray*(mappedAtlasInterpolation(nbest2,xspy,yspy,zspy)-mappedAtlasInterpolation(nbest2,xsmy,ysmy,zsmy));
								Jz = 0.5f/raz*(mappedAtlasInterpolation(nbest2,xspz,yspz,zspz)-mappedAtlasInterpolation(nbest2,xsmz,ysmz,zsmz));
								 
								J2 = Jx*Jx+Jy*Jy+Jz*Jz;
								den += diff*diff/sigma2 + J2;
								
								u[X][xyz] += diff*Jx;
								u[Y][xyz] += diff*Jy;
								u[Z][xyz] += diff*Jz;
									
								meanDiff += Numerics.abs(diff);
							}
						}							
					}
				}
				u[X][xyz] /= Numerics.max(ZERO,den);
				u[Y][xyz] /= Numerics.max(ZERO,den);
				u[Z][xyz] /= Numerics.max(ZERO,den);
			}
		}
		meanDiff /= (ntx*nty*ntz);
		
		if (regType==GAUSS_FLUID || regType==GAUSS_MIXED || regType==GAUSS_DIV_FLUID) {
			if (debug) BasicInfo.displayMessage("GAUSS_FLUID regularazation \n");
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFilters.separableConvolution(u[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Y] = ImageFilters.separableConvolution(u[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			u[Z] = ImageFilters.separableConvolution(u[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		} else if (regType==DIV_FLUID || regType==GAUSS_DIV_FLUID) {
		
			// smooth the result with a gaussian kernel
			u[X] = ImageFilters.separableConvolution(u[X],nsx,nsy,nsz,divKernelX,1,0,0);
			u[Y] = ImageFilters.separableConvolution(u[Y],nsx,nsy,nsz,divKernelY,0,1,0);
			u[Z] = ImageFilters.separableConvolution(u[Z],nsx,nsy,nsz,divKernelZ,0,0,1);
		}
		
		// compose the transformations
		if (fieldType==COMPOSITIVE) {
			if (debug) BasicInfo.displayMessage("compose with current transform \n");
							
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				int xyz = x+nsx*y+nsx*nsy*z;
				
				float xu = x+u[X][xyz];
				float yu = y+u[Y][xyz];
				float zu = z+u[Z][xyz];
				
				// note: if outside, extrapolate as X+u
				c[X][xyz] = ImageInterpolation.linearClosestInterpolation(s[X],xu,yu,zu,nsx,nsy,nsz) - x*scale*rtx/rax;
				c[Y][xyz] = ImageInterpolation.linearClosestInterpolation(s[Y],xu,yu,zu,nsx,nsy,nsz) - y*scale*rty/ray;
				c[Z][xyz] = ImageInterpolation.linearClosestInterpolation(s[Z],xu,yu,zu,nsx,nsy,nsz) - z*scale*rtz/raz;
			}
		}
		
		if (regType==GAUSS_DIFFUSION || regType==GAUSS_MIXED || regType==GAUSS_DIV_DIFFUSION) {
			if (debug) BasicInfo.displayMessage("GAUSS_DIFFUSION regularazation \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFilters.separableConvolution(c[X],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Y] = ImageFilters.separableConvolution(c[Y],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
			c[Z] = ImageFilters.separableConvolution(c[Z],nsx,nsy,nsz,gaussKernel,gx,gy,gz);
		} else if (regType==DIV_DIFFUSION || regType==GAUSS_DIV_DIFFUSION) {
			if (debug) BasicInfo.displayMessage("GAUSS_DIFFUSION regularazation \n");
			
			// smooth the result with a gaussian kernel
			c[X] = ImageFilters.separableConvolution(c[X],nsx,nsy,nsz,divKernelX,1,0,0);
			c[Y] = ImageFilters.separableConvolution(c[Y],nsx,nsy,nsz,divKernelY,0,1,0);
			c[Z] = ImageFilters.separableConvolution(c[Z],nsx,nsy,nsz,divKernelZ,0,0,1);
		}
					
		float meanC = 0.0f;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyz = x+nsx*y+nsx*nsy*z;
			
			s[X][xyz] = x*scale*rtx/rax + c[X][xyz];
			s[Y][xyz] = y*scale*rty/ray + c[Y][xyz];
			s[Z][xyz] = z*scale*rtz/raz + c[Z][xyz];
			
			meanC += c[X][xyz]*c[X][xyz]+c[Y][xyz]*c[Y][xyz]+c[Z][xyz]*c[Z][xyz];
		}
		meanC /= (nsx*nsy*nsz);
		
		if (debug) BasicInfo.displayMessage("convergence "+meanC+" -> "+meanDiff+"\n");

        return;
    } // 
    
    /** 
	 *	returns the transformed coordinates
	 */
	public final void getCurrentMapping(float[] Xs, int x, int y, int z) {
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
	public final float[][] exportTransformedAtlas() {
		float 	xs,ys,zs;
		float[][]	atl;
		
		atl = new float[na][ntx*nty*ntz];
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
		
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			// compute interpolated values
			for (int n=0;n<na;n++) {
				atl[n][xyz] = mappedAtlasInterpolation(n,xs,ys,zs);
			}
		}
		return atl;
	} // exportTransformedImage

    /** 
	 *	returns the transform field v defined as X' = X+v (=s)
	 */
	public final float[][] exportTransformField() {
		float 	xs,ys,zs;
		float[][]	vec = new float[3][ntx*nty*ntz];
		
        for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			int xyz = x+ntx*y+ntx*nty*z;
			
			xs = ImageInterpolation.linearClosestInterpolation(s[X],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			ys = ImageInterpolation.linearClosestInterpolation(s[Y],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			zs = ImageInterpolation.linearClosestInterpolation(s[Z],x/scale,y/scale,z/scale,nsx,nsy,nsz);
			
			vec[X][xyz] = xs-x*rtx/rax;
			vec[Y][xyz] = ys-y*rty/ray;
			vec[Z][xyz] = zs-z*rtz/raz;			
 		}
		return vec;
	} // exportTransformedImage

	private final float mappedAtlasInterpolation(int n, float xs, float ys, float zs) {
		return Numerics.bounded(atlasScale*(ImageInterpolation.linearClosestInterpolation(atlas[n],xs,ys,zs,nax,nay,naz)-0.5f)+0.5f,0.0f,1.0f);
	}	
}
