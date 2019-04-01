package de.mpg.cbs.methods;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm computes basic total variation minimizers
 *  from Chambolle's algorithm
 *
 *	@version    July 2014
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class TotalVariation {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// data buffers
	private 	float[][][]			image;  			// original image
	private 	float[][][][]		proj;				// tv-projection function
	private 	boolean[][][]		mask;   			// image mask: true for data points
	private		static	int			nx,ny,nz;   		// image dimensions
	private		float 				Imin, Imax;
	private		float				Isize;				// non-masked region size
	
	// parameters
	private 	float 		lambdaScale = 0.05f;		// scaling parameter
	private 	float 		tauStep = 0.125f;		// internal step parameter (default 1/4)
	private 	float 		maxdist = 0.0001f;		// maximum error for stopping
	private 	int 		maxiter = 100;		// maximum number of iterations
	private     float       wrapping = 1.0f;    // wrapping value, as a function of [Imin, Imax]
			
	// computation variables
	private 	float[][][]			field;  // image difference field
	
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	
	// for debug
	static final boolean		debug=true;
	static final boolean		verbose=true;
	
	/**
	 *  constructor
	 *	note: all images passed to the algorithm are just linked, not copied
	 */
	 
	public TotalVariation(float[][][] image_, boolean [][][] mask_, 
					int nx_, int ny_, int nz_) {
		
		image = image_;
		mask = mask_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		// init all the new arrays
		try {
		    //image = new float[nx][ny][nz];
		    //mask = new boolean[nx][ny][nz];
			proj = new float[3][nx][ny][nz];
			field = new float[nx][ny][nz];
		} catch (OutOfMemoryError e){
			finalize();
			System.out.println(e.getMessage());
			return;
		}
		/*
		// copy data on the boundary?
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    image[x][y][z] = image_[Numerics.bounded(x-1, 0, nx-3)][Numerics.bounded(y-1, 0, ny-3)][Numerics.bounded(z-1, 0, nz-3)];
		    mask[x][y][z] = mask_[Numerics.bounded(x-1, 0, nx-3)][Numerics.bounded(y-1, 0, ny-3)][Numerics.bounded(z-1, 0, nz-3)];
		}*/
		    
		// init values
		Imin = INF;
		Imax = -INF;
		Isize = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];

			proj[X][x][y][z] = 0.0f;
			proj[Y][x][y][z] = 0.0f;
			proj[Z][x][y][z] = 0.0f;
			
			if (mask[x][y][z]) Isize++;
		}
		// rescale intensities in [0,1]
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			image[x][y][z] = (image[x][y][z]-Imin)/(Imax-Imin);
		}
		// scale parameter: invariant by intensity scaling
		//lambdaScale *= (Imax-Imin)*(Imax-Imin);

		// the [Imin, Imax] part is taken care of
		wrapping = 1.0f/lambdaScale;

		if (debug) BasicInfo.displayMessage("TV:initialisation\n");
	}

	public TotalVariation(float[][][] image_, boolean [][][] mask_, 
					int nx_, int ny_, int nz_,
					float scale_, float step_, float dist_, int iter_) {
		
		image = image_;
		mask = mask_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		lambdaScale = scale_;
		tauStep = step_;
		maxdist = dist_;
		maxiter = iter_;

		// init all the new arrays
		try {
			//image = new float[nx][ny][nz];
		    //mask = new boolean[nx][ny][nz];
			proj = new float[3][nx][ny][nz];
			field = new float[nx][ny][nz];
		} catch (OutOfMemoryError e){
			finalize();
			System.out.println(e.getMessage());
			return;
		}
		/*
		// copy data on the boundary?
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    // copy the previous slice? or mirror?
		    //image[x][y][z] = image_[Numerics.bounded(x-1, 0, nx-3)][Numerics.bounded(y-1, 0, ny-3)][Numerics.bounded(z-1, 0, nz-3)];
		    //mask[x][y][z] = mask_[Numerics.bounded(x-1, 0, nx-3)][Numerics.bounded(y-1, 0, ny-3)][Numerics.bounded(z-1, 0, nz-3)];
		    int x0 = x-1;
		    if (x0<0) x0 = -x0;
		    if (x0>nx-3) x0 = nx-3 - (x0-(nx-3));
		    
		    int y0 = y-1;
		    if (y0<0) y0 = -y0;
		    if (y0>ny-3) y0 = ny-3 - (y0-(ny-3));
		    
		    int z0 = z-1;
		    if (z0<0) z0 = -z0;
		    if (z0>nz-3) z0 = nz-3 - (z0-(nz-3));
		    
		    image[x][y][z] = image_[x0][y0][z0];
		    mask[x][y][z] = mask_[x0][y0][z0];
		}*/
		
		// init values
		Imin = INF;
		Imax = -INF;
		Isize = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];

			proj[X][x][y][z] = 0.0f;
			proj[Y][x][y][z] = 0.0f;
			proj[Z][x][y][z] = 0.0f;
			
			if (mask[x][y][z]) Isize++;
		}
		// rescale intensities in [0,1]
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			image[x][y][z] = (image[x][y][z]-Imin)/(Imax-Imin);
		}
		//lambdaScale *= (Imax-Imin)*(Imax-Imin);

		// the [Imin, Imax] part is taken care of
		wrapping = 1.0f/lambdaScale;
		
		if (debug) BasicInfo.displayMessage("TV:initialisation\n");
	}

	/** clean-up: destroy membership and centroid arrays */
	public final void finalize() {
		proj = null;
		field = null;
		System.gc();
	}
	
    /** accessor for computed data */ 
    public final float[][][][] getProjection() { return proj; }
    
    /** 
	 *  compute the new TV projection function
	 */
    final public double computeProjection() {
        double distance, dist;
        int x,y,z;
        double norm, divp;
        double[] grad = new double[3];
        double[] prev = new double[3];
        double num = 0.0;
        
		distance = 0.0;
		// we ignore the image boundary for convenience
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (mask[x][y][z]) {
			
			// first compute the intermediate field
			divp = 0.0;
			if (x==0) divp += (proj[X][x][y][z]);
			else if (x==nx-1) divp += (-proj[X][x-1][y][z]);
			else divp += (proj[X][x][y][z]-proj[X][x-1][y][z]);
			
			if (y==0) divp += (proj[Y][x][y][z]);
			else if (y==ny-1) divp += (-proj[Y][x][y-1][z]);
			else divp += (proj[Y][x][y][z]-proj[Y][x][y-1][z]);
			
			if (z==0) divp += (proj[Z][x][y][z]);
			else if (z==nz-1) divp += (-proj[Z][x][y][z-1]);
			else divp += (proj[Z][x][y][z]-proj[Z][x][y][z-1]);
			
			field[x][y][z] = (float)(divp - image[x][y][z]/lambdaScale);
		}
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (mask[x][y][z]) {
			// second compute the projection update
			prev[X] = proj[X][x][y][z];
			prev[Y] = proj[Y][x][y][z];
			prev[Z] = proj[Z][x][y][z];

			/*
			if (x==0) grad[X] = field[x+1][y][z];
			else if (x==nx-1) grad[X] = -field[x][y][z];
            else grad[X] = field[x+1][y][z]-field[x][y][z];
			
            if (y==0) grad[Y] = field[x][y+1][z];
			else if (y==ny-1) grad[Y] = -field[x][y][z];
			else grad[Y] = field[x][y+1][z]-field[x][y][z];
			
			if (z==0) grad[Z] = field[x][y][z+1];
			else if (z==nz-1) grad[Z] = -field[x][y][z];
			else grad[Z] = field[x][y][z+1]-field[x][y][z];
			*/
			if (x+1<nx) grad[X] = (field[x+1][y][z]-field[x][y][z]);
			else grad[X] = 0.0;
            if (y+1<ny) grad[Y] = (field[x][y+1][z]-field[x][y][z]);
			else grad[Y] = 0.0;
            if (z+1<nz) grad[Z] = (field[x][y][z+1]-field[x][y][z]);
			else grad[Z] = 0.0;

			norm = 1.0 + tauStep*FastMath.sqrt(grad[X]*grad[X]+grad[Y]*grad[Y]+grad[Z]*grad[Z]);
			
			proj[X][x][y][z] = (float)( (prev[X] + tauStep*grad[X])/norm);
			proj[Y][x][y][z] = (float)( (prev[Y] + tauStep*grad[Y])/norm);
			proj[Z][x][y][z] = (float)( (prev[Z] + tauStep*grad[Z])/norm);
			
			dist = 	 (proj[X][x][y][z]-prev[X])*(proj[X][x][y][z]-prev[X])
					+(proj[Y][x][y][z]-prev[Y])*(proj[Y][x][y][z]-prev[Y])
					+(proj[Z][x][y][z]-prev[Z])*(proj[Z][x][y][z]-prev[Z]);
					
			//if (dist>distance) distance = dist;
			distance += dist;
			num++;
		}

        return distance/num;
    } // computeMemberships
    
    /** 
	 *  compute the new TV projection function
	 */
    final public double computeProjectionWrapped() {
        double distance, dist;
        int x,y,z;
        double norm, divp;
        double[] grad = new double[3];
        double[] prev = new double[3];
        double num = 0.0;
        
		distance = 0.0;
		// we ignore the image boundary for convenience
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (mask[x][y][z]) {
			
			// first compute the intermediate field
			divp = 0.0;
			if (x==0) divp += Numerics.modulo(proj[X][x][y][z], wrapping);
			else if (x==nx-1) divp += Numerics.modulo(-proj[X][x-1][y][z], wrapping);
			else divp += Numerics.modulo(proj[X][x][y][z]-proj[X][x-1][y][z], wrapping);
			
			if (y==0) divp += Numerics.modulo(proj[Y][x][y][z], wrapping);
			else if (y==ny-1) divp += Numerics.modulo(-proj[Y][x][y-1][z], wrapping);
			else divp += Numerics.modulo(proj[Y][x][y][z]-proj[Y][x][y-1][z], wrapping);
			
			if (z==0) divp += Numerics.modulo(proj[Z][x][y][z], wrapping);
			else if (z==nz-1) divp += Numerics.modulo(-proj[Z][x][y][z-1], wrapping);
			else divp += Numerics.modulo(proj[Z][x][y][z]-proj[Z][x][y][z-1], wrapping);
			//if (x-1>=0) divp += Numerics.modulo(proj[X][x][y][z]-proj[X][x-1][y][z], wrapping); 
			//if (y-1>=0) divp += Numerics.modulo(proj[Y][x][y][z]-proj[Y][x][y-1][z], wrapping);
			//if (z-1>=0) divp += Numerics.modulo(proj[Z][x][y][z]-proj[Z][x][y][z-1], wrapping);
			//divp = (proj[X][x][y][z]-proj[X][x-1][y][z]) 
			//	 + (proj[Y][x][y][z]-proj[Y][x][y-1][z]) 
			//	 + (proj[Z][x][y][z]-proj[Z][x][y][z-1]);
			
			//field[x][y][z] = (float)(divp - image[x][y][z]/lambdaScale);
			field[x][y][z] = (float)Numerics.modulo(divp - image[x][y][z]/lambdaScale, wrapping);
		}
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (mask[x][y][z]) {
			// second compute the projection update
			prev[X] = proj[X][x][y][z];
			prev[Y] = proj[Y][x][y][z];
			prev[Z] = proj[Z][x][y][z];
			
			/*
			if (x==0) grad[X] = Numerics.modulo(field[x+1][y][z], wrapping);
			else if (x==nx-1) grad[X] = Numerics.modulo(-field[x][y][z], wrapping);
            else grad[X] = Numerics.modulo(field[x+1][y][z]-field[x][y][z], wrapping);
			
            if (y==0) grad[Y] = Numerics.modulo(field[x][y+1][z], wrapping);
			else if (y==ny-1) grad[Y] = Numerics.modulo(-field[x][y][z], wrapping);
			else grad[Y] = Numerics.modulo(field[x][y+1][z]-field[x][y][z], wrapping);
			
			if (z==0) grad[Z] = Numerics.modulo(field[x][y][z+1], wrapping);
			else if (z==nz-1) grad[Z] = Numerics.modulo(-field[x][y][z], wrapping);
			else grad[Z] = Numerics.modulo(field[x][y][z+1]-field[x][y][z], wrapping);
			*/
			if (x+1<nx) grad[X] = Numerics.modulo(field[x+1][y][z]-field[x][y][z], wrapping);
			else grad[X] = 0.0;
            if (y+1<ny) grad[Y] = Numerics.modulo(field[x][y+1][z]-field[x][y][z], wrapping);
			else grad[Y] = 0.0;
            if (z+1<nz) grad[Z] = Numerics.modulo(field[x][y][z+1]-field[x][y][z], wrapping);
			else grad[Z] = 0.0;
            //grad[X] = (field[x+1][y][z]-field[x][y][z]);
			//grad[Y] = (field[x][y+1][z]-field[x][y][z]);
			//grad[Z] = (field[x][y][z+1]-field[x][y][z]);
			norm = 1.0 + tauStep*FastMath.sqrt(grad[X]*grad[X]+grad[Y]*grad[Y]+grad[Z]*grad[Z]);
			
			proj[X][x][y][z] = (float)( (prev[X] + tauStep*grad[X])/norm);
			proj[Y][x][y][z] = (float)( (prev[Y] + tauStep*grad[Y])/norm);
			proj[Z][x][y][z] = (float)( (prev[Z] + tauStep*grad[Z])/norm);
			
			dist = 	 Numerics.modulo(proj[X][x][y][z]-prev[X], wrapping)*Numerics.modulo(proj[X][x][y][z]-prev[X], wrapping)
					+Numerics.modulo(proj[Y][x][y][z]-prev[Y], wrapping)*Numerics.modulo(proj[Y][x][y][z]-prev[Y], wrapping)
					+Numerics.modulo(proj[Z][x][y][z]-prev[Z], wrapping)*Numerics.modulo(proj[Z][x][y][z]-prev[Z], wrapping);
					
			//if (dist>distance) distance = dist;
			distance += dist;
			num++;
		}

        return distance/num;
    } // computeMemberships
    
    /**
	 * denoising algorithm
	 */
    final public void solve() {
    	// loop until distance is minimized
    	double distance = 1e10;
    	int t=0;
    	while ((distance>maxdist || t<0) && t<maxiter) {
    		t++;
    		if (verbose) BasicInfo.displayMessage("iter "+t);
			// get the new projection
			distance = computeProjection();
			
			if (verbose) BasicInfo.displayMessage(": d="+distance+"\n");
        }
        return;
    }
        
    /**
	 * denoising algorithm
	 */
    final public void solveWrapped() {
    	// loop until distance is minimized
    	double distance = 1e10;
    	int t=0;
    	while ((distance>maxdist || t<0) && t<maxiter) {
    		t++;
    		if (verbose) BasicInfo.displayMessage("iter "+t);
			// get the new projection
			distance = computeProjectionWrapped();
			
			if (verbose) BasicInfo.displayMessage(": d="+distance+"\n");
        }
        return;
    }
        
    /**
	 * denoising algorithm
	 */
    final public void denoiseImage(float stdev, boolean adaptive) {
    	double fn;
    	float ratio = 1.0f;
    	
    	// loop until distance is minimized
    	double distance = 1e10;
    	int t=0;
    	double sqdist = FastMath.sqrt(maxdist);
    	while ((distance>maxdist || t<0) && t<maxiter) {
    		t++;
    		if (verbose) BasicInfo.displayMessage("iter "+t);
			// get the new projection
			distance = computeProjection();
			
			// get the scaling factor
			fn = 0.0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
				double divp = 0.0;
                if (x==0) divp += (proj[X][x][y][z]);
                else if (x==nx-1) divp += (-proj[X][x-1][y][z]);
                else divp += (proj[X][x][y][z]-proj[X][x-1][y][z]);
                
                if (y==0) divp += (proj[Y][x][y][z]);
                else if (y==ny-1) divp += (-proj[Y][x][y-1][z]);
                else divp += (proj[Y][x][y][z]-proj[Y][x][y-1][z]);
                
                if (z==0) divp += (proj[Z][x][y][z]);
                else if (z==nz-1) divp += (-proj[Z][x][y][z-1]);
                else divp += (proj[Z][x][y][z]-proj[Z][x][y][z-1]);
			
				fn += divp*divp;
			}
			fn = FastMath.sqrt(fn/Isize)*lambdaScale;
			//fn *= lambdaScale/Isize;
			ratio = (float)(stdev/fn);
	
			// change scaling only after some level of convergence
			if (adaptive && distance<sqdist) {
				// update the lambda
				lambdaScale *= ratio;
			}
			if (verbose) BasicInfo.displayMessage(": d="+distance+", f="+fn+", r="+ratio+"\n");
        }
        return;
    }
        
	/** 
	 *	export membership functions 
	 */
	public final float[][][] exportResult() {
		int 	x,y,z,k;
		/*
		float[][][]	res = new float[nx-2][ny-2][nz-2];
		float divp;
       
		for (x=2;x<nx-1;x++) for (y=2;y<ny-1;y++) for (z=2;z<nz-1;z++) {
			divp = (proj[X][x-1][y-1][z-1]-proj[X][x-2][y-1][z-1]) 
				 + (proj[Y][x-1][y-1][z-1]-proj[Y][x-1][y-2][z-1]) 
				 + (proj[Z][x-1][y-1][z-1]-proj[Z][x-1][y-1][z-2]);
				 
			res[x-1][y-1][z-1] = image[x][y][z] - lambdaScale*divp;
		}*/
		
		float[][][]	res = new float[nx][ny][nz];
		double divp;
       
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
		    divp = 0.0;
			if (x==0) divp += (proj[X][x][y][z]);
			else if (x==nx-1) divp += (-proj[X][x-1][y][z]);
			else divp += (proj[X][x][y][z]-proj[X][x-1][y][z]);
			
			if (y==0) divp += (proj[Y][x][y][z]);
			else if (y==ny-1) divp += (-proj[Y][x][y-1][z]);
			else divp += (proj[Y][x][y][z]-proj[Y][x][y-1][z]);
			
			if (z==0) divp += (proj[Z][x][y][z]);
			else if (z==nz-1) divp += (-proj[Z][x][y][z-1]);
			else divp += (proj[Z][x][y][z]-proj[Z][x][y][z-1]);
							 
			res[x][y][z] = (float)(Imin + (image[x][y][z] - lambdaScale*divp)*(Imax-Imin));
		}
		return res;
	} // exportMemberships
        
	/** 
	 *	export membership functions 
	 */
	public final float[][][] exportResultWrapped() {
		int 	x,y,z,k;
		/*
		float[][][]	res = new float[nx-2][ny-2][nz-2];
		double divp;
       
		for (x=2;x<nx-1;x++) for (y=2;y<ny-1;y++) for (z=2;z<nz-1;z++) {
			divp = Numerics.modulo(proj[X][x-1][y-1][z-1]-proj[X][x-2][y-1][z-1], wrapping) 
				 + Numerics.modulo(proj[Y][x-1][y-1][z-1]-proj[Y][x-1][y-2][z-1], wrapping) 
				 + Numerics.modulo(proj[Z][x-1][y-1][z-1]-proj[Z][x-1][y-1][z-2], wrapping);
			//divp = (proj[X][x][y][z]-proj[X][x-1][y][z]) 
			//	 + (proj[Y][x][y][z]-proj[Y][x][y-1][z]) 
			//	 + (proj[Z][x][y][z]-proj[Z][x][y][z-1]);
				 
			//res[x][y][z] = (float)(image[x][y][z] - lambdaScale*divp);
			//res[x][y][z] = (float)Numerics.modulo(image[x][y][z] - lambdaScale*divp, Imax-Imin);
			res[x-1][y-1][z-1] = (float)(-lambdaScale*divp);
		}*/
		float[][][]	res = new float[nx][ny][nz];
		double divp;
       
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
			divp = 0.0;
			if (x==0) divp += Numerics.modulo(proj[X][x][y][z], wrapping);
			else if (x==nx-1) divp += Numerics.modulo(-proj[X][x-1][y][z], wrapping);
			else divp += Numerics.modulo(proj[X][x][y][z]-proj[X][x-1][y][z], wrapping);
			
			if (y==0) divp += Numerics.modulo(proj[Y][x][y][z], wrapping);
			else if (y==ny-1) divp += Numerics.modulo(-proj[Y][x][y-1][z], wrapping);
			else divp += Numerics.modulo(proj[Y][x][y][z]-proj[Y][x][y-1][z], wrapping);
			
			if (z==0) divp += Numerics.modulo(proj[Z][x][y][z], wrapping);
			else if (z==nz-1) divp += Numerics.modulo(-proj[Z][x][y][z-1], wrapping);
			else divp += Numerics.modulo(proj[Z][x][y][z]-proj[Z][x][y][z-1], wrapping);
				 
			//res[x][y][z] = (float)(Imin + (image[x][y][z] - lambdaScale*divp)*(Imax-Imin));
			res[x][y][z] = (float)( - lambdaScale*divp);
		}
		return res;
	} // exportMemberships

}
