package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;

/**
 *
 *  This class handles basic rigid registration
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class BasicRigidRegistration {
	
	// image quantities
	private		float[]				source;
	private		float[]				target;
	
	private 	int 				nx,ny,nz; 			// image dimensions
	private 	float 				rx,ry,rz; 			// image resolutions
	
	// spatial transformations : only rigid
	private		float[]					transform;		// the transform parameters to get into image space
	private		float[][]				rotation;		// the associated rotation matrix
	private		ParametricTransform		transformModel;	// the type of transform (from possible ones below)
	private		RotationMatrix			rmat;
	private		float[][]				dRa, dRb, dRc;
	private		final int				Nd;
	
	
	// gradient descent parameters
	private		float		lambda	=	1.0f;		// the Levenberg-Marquardt parameter for Levenberg Marquardt estimation
	private		float		lfactor = 1.5f;			// the multiplicative factor for the adaptive term
	private		int			itSupp = 10;				// maximum of steps if the cost function is not improving
	private		int			itMax,itPlus,Nturn;		// counters for various loops
	private static final	float   INIT_LAMBDA = 1;
	private		float		minEdiff = 1e-6f;		// the minimum variation of energy to require a better alignment
	private		float		minLambda = 0.001f;		// the minimum variation of energy to require a better alignment
	private		int			subsample = 3;			// scale for the registration: just subsample the volume 
	private		int			offset = 0;				// offset used in subsampling (cyclic)
	private		float		maxdiff;
	
	
	// preset computation arrays for speed up
	private		float[] 	hessian, gradient;
	private		float[]		trial;	
	private		float[] 	Xi;
    private		float[][] 	dXi;
    private		float[]		dsource;
		
	// constants
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private static final	float	ISQRT2 =  (float)(1.0/Math.sqrt(2.0f));
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=false;
	
	
	// numerics
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;

	// convenience tags
	private static final	int   X=0;
	private static final	int   Y=1;
	private static final	int   Z=2;
	private static final	int   T=3;

	/**
	 *	constructor: load the atlas information from a file.
	 *	<p>
	 *	The atlas files follow a certain template; 
	 *  the separator between numbers is a tab, not a space.
	 */
	public BasicRigidRegistration(float[] src_, float[] trg_, 
									int nx_, int ny_, int nz_, float rx_, float ry_, float rz_,
									int iter_, float diff_, int sub_) {
		
		source = src_;
		target = trg_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
		
		itMax = iter_;
		maxdiff = diff_;
		subsample = sub_;
		
		// transform
		transformModel = new ParametricTransform("rigid", nx/2.0f,ny/2.0f,nz/2.0f, rx,ry,rz, nx, ny, nz, nx/2.0f,ny/2.0f,nz/2.0f, rx,ry,rz, nx,ny,nz);
		
		Nd = transformModel.getDimension();
		
		// init parameters
		transform = new float[Nd];
		trial = new float[Nd];
		gradient = new float[Nd];
		hessian = new float[Nd];
		Xi = new float[3];
		dXi = new float[3][Nd];
		dsource = new float[Nd];
			
		for (int n=0;n<Nd;n++) transform[n] = 0.0f;
		rotation = transformModel.computeRotation(transform);
	}
	
	public final void finalize() {
		source = null;
		target = null;
	}
	
	public final float[] getSource() { return source; }
	
	public final float[] getTarget() { return target; }
	
	public final float[] getTransform() { return transform; }
	
	public final float[][] exportTransformMatrix() {
		float[][] matrix = new float[3][4];
		rotation = transformModel.computeRotation(transform);
		transformModel.precomputeImageToTemplateMatrix(matrix, transform, rotation, 1.0f);
		return matrix;
	}
	
    /** 
	 *	runs a Levenberg-Marquardt step
	 */
	public final void register() {    
		boolean stop;
		
        // one level
        if (debug) BasicInfo.displayMessage("registration\n");
        lambda = INIT_LAMBDA;
		Nturn = 0; itPlus = 0; stop = false;
		float E0 = computeRegistrationCoefficients(transform);
		float E;
		float diff = 1.0f;
		for (int n=0;n<itMax && diff>minEdiff && lambda>minLambda;n++) {
			if (verbose) System.out.print("iteration "+n);
			E = registerGradientDescent();
			diff = (E-E0)/E0;
			if (verbose) System.out.print(" -> E = "+E+" ("+diff+")\n");
			E0 = E;
		}
		// update the rotation coefficients
		rotation = transformModel.computeRotation(transform);
    }
  
	/**
	 * compute the new transform using gradient descent
	 * performs only one iteration
	 */
    final private float registerGradientDescent() {
		float		E0,E,Eprev;
		boolean		stop = false;
		boolean		admissible = true;
		
		// compute the coefficients at current transform
		E0 = computeRegistrationCoefficients(transform);
		E = E0;
		
		// search along the line
		int iter = 0;
		while (!stop) {
			Eprev = E;
			
			// new values for the coeffs 
			for (int n=0; n<Nd; n++) 
				trial[n] = transform[n] + lambda/Numerics.max(ZERO,hessian[n])*gradient[n];
			
			// is the energy better ?
			E = computeRegistrationEnergy(trial);
			
			if (debug) System.out.print( "a: "+lambda+"("+E0+")->("+E+")\n");

			// test on energy value : maximisation
			if ( E > E0 ) {
				// better value: changing the transform and the scale
				for (int n=0; n<Nd; n++) 
					transform[n] = trial[n];
				
				lambda = lambda*lfactor;
				stop = true;
				if (verbose) BasicInfo.displayMessage("E = "+E+" ("+displayTransform(transform)+")\n");
			} 
			else {
				lambda = lambda/lfactor;
				iter++;
				E = Eprev;
				if (iter>itSupp) {
					stop = true;
					if (debug) System.out.print( "stop search\n");
				}
			}
		}
		
		return E;
	} // registerGradientDescent

	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeRegistrationCoefficients(float[] trans) {
        float xP=0,yP=0,zP=0;
        float dPx,dPy,dPz,sourceT,trg;
     	int				Npt;
		float			weight;
		float			cost,norm;
	
        // init the coefficients
		cost = 0.0f;
		norm = 0.0f;
		for (int i=0;i<Nd;i++) {
			hessian[i] = 0.0f;
			gradient[i] = 0.0f;
		}
		
        // set up rotation parameters
        rmat = transformModel.computeRotationMatrix(trans);
		rotation = rmat.getMatrix();
		dRa = rmat.derivatives(1.0f, 0.0f, 0.0f);
		dRb = rmat.derivatives(0.0f, 1.0f, 0.0f);
		dRc = rmat.derivatives(0.0f, 0.0f, 1.0f);
		
		if (debug) System.out.println(displayTransform(trans));	
			
		// main loop
		for (int x=offset;x<nx;x+=subsample) for (int y=offset;y<ny;y+=subsample) for (int z=offset;z<nz;z+=subsample) {
            // compute the local position
			transformModel.imageToTemplate(Xi,x,y,z,trans,rotation,1.0f);
				
			// compute interpolated values
			sourceT = ImageInterpolation.linearInterpolation(source,0.0f,Xi[0],Xi[1],Xi[2],nx,ny,nz);
				
			int xyz = x + y*nx + z*nx*ny;
			trg = target[xyz];
			
			if (sourceT>0 && trg>0) {
				// data terms
				weight = imageRegistrationWeight(trg,sourceT);
						
				// spatial derivatives
				dPx = ImageInterpolation.linearInterpolationXderivative(source,Xi[0],Xi[1],Xi[2],nx,ny,nz);
				dPy = ImageInterpolation.linearInterpolationYderivative(source,Xi[0],Xi[1],Xi[2],nx,ny,nz);
				dPz = ImageInterpolation.linearInterpolationZderivative(source,Xi[0],Xi[1],Xi[2],nx,ny,nz);
						
				// coordinate derivatives
				transformModel.imageToTemplateDerivatives(dXi,x,y,z,trans,rotation,dRa,dRb,dRc,1.0f);
								
				// assemble everything
				for (int i=0;i<Nd;i++) {
					dsource[i] = dPx*dXi[0][i] + dPy*dXi[1][i] + dPz*dXi[2][i];
				}
				cost += registrationCost(weight,trg,sourceT);
				norm += registrationNorm(weight,trg,sourceT);
	
				for (int i=0;i<Nd;i++) {
					gradient[i] += registrationCostGradient(weight,trg,sourceT,dsource,i);
					hessian[i] += registrationCostHessian(weight,trg,sourceT,dsource,i);
				}
			}
       }
		if (cost>ZERO) {
			for (int i=0;i<Nd;i++) {
				gradient[i] = gradient[i]/norm;
				hessian[i]  = hessian[i]/norm;
			}
		} else {
			for (int i=0;i<Nd;i++) {
				gradient[i] = 0.0f;
				hessian[i]  = 1.0f;
			}
		}
		if (debug) {
			System.out.println(displayVector(gradient));
			System.out.println(displayVector(hessian));
		}
        return cost/norm;
    } // computeRegistrationCoefficients
    
	/** transformations: how to update the transform using memberships */
	/**
	 * compute the image position given the membership functions
	 * for a given level l
     * performs only one iteration
	 */
    final private float computeRegistrationEnergy(float[] trans) {
        float weight;
        float dPx,dPy,dPz,sourceT;
        float			cost,norm;
        float trg;
		
        // init the coefficients
		cost = 0.0f;
		norm = 0.0f;
		
        // set up rotation parameters
       rmat = transformModel.computeRotationMatrix(trans);
	   rotation = rmat.getMatrix();
		
	   if (debug) System.out.println(displayTransform(trans));	
		
		// main loop
		for (int x=offset;x<nx;x+=subsample) for (int y=offset;y<ny;y+=subsample) for (int z=offset;z<nz;z+=subsample) {
            // compute the local position
			transformModel.imageToTemplate(Xi,x,y,z,trans,rotation,1.0f);
			
			// compute interpolated values
			sourceT = ImageInterpolation.linearInterpolation(source,0.0f,Xi[0],Xi[1],Xi[2],nx,ny,nz);
			
			int xyz = x + y*nx + z*nx*ny;
			trg = target[xyz];
			// check if the region is zero: no calculation needed then
			if (sourceT>0 && trg>0) {
					
				// data term : function of the memberships
				weight = imageRegistrationWeight(trg,sourceT);
				
				cost += registrationCost(weight,trg,sourceT);
				norm += registrationNorm(weight,trg,sourceT);
			}
        }
 
        return cost/norm;
    } // computeRegistrationEnergy
    
	final private float registrationCost(float w, float t, float sT) {
		return w*t*t*sT*sT;
		//return w*(t-sT)*(t-sT);
	}
	
	final private float registrationCostGradient(float w, float t, float sT, float[] dsT, int i) {
		return 2.0f*w*t*t*sT*dsT[i];
		//return -2.0f*w*(t-sT)*dsT[i];
	}
	
	final private float registrationCostHessian(float w, float t, float sT, float[] dsT, int i) {
		return 2.0f*w*t*t*dsT[i]*dsT[i];
		//return 2.0f*w*dsT[i]*dsT[i];
	}
	
	final private float registrationNorm(float w, float t, float sT) {
		return w;
	}
   	
	/**
	 * gives the weighting factor for the image data
	 */
    final private float imageRegistrationWeight(float trg, float src) {
		return Numerics.max(trg,src);
	}
	
	public final String displayTransform(float[] trans) {
		String info = "transform: (";
		for (int n=0;n<Nd-1;n++) info += trans[n]+", ";
		info += trans[Nd-1]+")\n";
		
		return info;
	}
	
	/** display the atlas data */
	public final String displayVector(float[] vect) {
		String info = "vector: (";
		for (int n=0;n<vect.length-1;n++) info += vect[n]+", ";
		info += vect[vect.length-1]+")\n";
		
		return info;
	}
	
	/** 
	 *  mapp the source to target with given transformation
	 */
    final public float[] transformedSource() {
		float[] img = new float[nx*ny*nz];
		float[] XP=new float[3];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			transformModel.imageToTemplate(XP,x,y,z,transform,rotation, 1.0f);
			
			img[x+y*nx+z*nx*ny] = ImageInterpolation.linearInterpolation(source,0.0f,XP[0],XP[1],XP[2],nx,ny,nz);	
		}
		return img;
	}

}
