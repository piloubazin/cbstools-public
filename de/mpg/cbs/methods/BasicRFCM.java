package de.mpg.cbs.methods;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/**
 *
 *  This algorithm handles the main Fuzzy C-Means operations
 *  for 3D data in FANTASM (membership, centroid computations)
 *	with limited options
 *
 *	@version    March 2014
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class BasicRFCM {
		
	// numerical quantities
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;
	
	// data buffers
	private 	float[][][]			image;  			// original image
	private 	float[][][][]		mems;				// membership function
	private 	float[]				centroids;			// cluster centroids
	private 	boolean[][][]		mask;   			// image mask: true for data points
	private		static	int			nx,ny,nz;   		// image dimensions
	private		float 				Imin, Imax;
	
	// parameters
	private 	int 		classes;    // number of classes
	private 	float 		smoothing;	// MRF smoothing
			
	// computation variables
	private		float[]		prev;	// previous membership values
	
	// for debug
	static final boolean		debug=true;
	static final boolean		verbose=true;
	
	/**
	 *  constructor
	 *	note: all images passed to the algorithm are just linked, not copied
	 */
	public BasicRFCM(float[][][] image_, boolean [][][] mask_, 
					int nx_, int ny_, int nz_,
					int classes_, float smoothing_) {
		
		image = image_;
		mask = mask_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		classes = classes_;
		smoothing = smoothing_;

		// init all the new arrays
		try {
			//mems = new float[classes][nx][ny][nz];
			mems = new float[nx][ny][nz][classes];
			centroids = new float[classes];
			prev = new float[classes];
		} catch (OutOfMemoryError e){
			finalize();
			System.out.println(e.getMessage());
			return;
		}
		
		// init values
		Imin = INF;
		Imax = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];

			for (int k=0;k<classes;k++) {
				mems[x][y][z][k] = 0.0f;
			}
		}
		for (int k=0;k<classes;k++) {
			centroids[k] = 0.0f;
		}

		smoothing *= (Imax-Imin);
		
		if (debug) BasicInfo.displayMessage("FCM:initialisation\n");
	}

	/** clean-up: destroy membership and centroid arrays */
	public final void finalize() {
		mems = null;
		centroids = null;
		prev = null;
		System.gc();
	}
	
    /** accessor for computed data */ 
    public final float[][][][] getMemberships() { return mems; }
    /** accessor for computed data */ 
    public final float[] getCentroids() { return centroids; }
    /** accessor for computed data */ 
	final public void setCentroids(float[] cent) { 
		centroids = cent; 
		if (debug) {
			BasicInfo.displayMessage("centroids: ("+centroids[0]);
			for (int k=1;k<classes;k++) BasicInfo.displayMessage(", "+centroids[k]);
			BasicInfo.displayMessage(")\n");
		}

	}
	
	/** estimate the starting values for the centroids based on data range (simple, but not robust) */
	final public void initCentroidsRange() {
		centroids[0] = Imin + 0.5f*(Imax-Imin)/(float)classes;
		for (int k=1;k<classes;k++) {
			centroids[k] = centroids[k-1] + (Imax-Imin)/(float)classes;
		}
	}
	
    /** 
	 *  compute the FCM membership functions given the centroids
	 */
    final public float computeMemberships() {
        float distance,dist;
        int x,y,z,k,m;
        float den,num;
        float neighbors, ngb;
        
		distance = 0.0f;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) if ( mask[x][y][z] ) {
			den = 0;
			// remember the previous values
			for (k=0;k<classes;k++) {
				prev[k] = mems[x][y][z][k];
			}
			for (k=0;k<classes;k++) {
					
				// data term
				num = (image[x][y][z]-centroids[k])*(image[x][y][z]-centroids[k]);
					
				// spatial smoothing
				if (smoothing > 0.0f) { 
					ngb = 0.0f;  
					neighbors = 0.0f;
					// 6-connected neighbors
					if (mask[x+1][y][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x+1][y][z][m]*mems[x+1][y][z][m];
						neighbors ++;
					}
					if (mask[x-1][y][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x-1][y][z][m]*mems[x-1][y][z][m];
						neighbors ++;
					}
					if (mask[x][y+1][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y+1][z][m]*mems[x][y+1][z][m];
						neighbors ++;
					}
					if (mask[x][y-1][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y-1][z][m]*mems[x][y-1][z][m];
						neighbors ++;
					}
					if (mask[x][y][z+1]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y][z+1][m]*mems[x][y][z+1][m];
						neighbors ++;
					}
					if (mask[x][y][z-1]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y][z-1][m]*mems[x][y][z-1][m];
						neighbors ++;
					}
					if (neighbors>0.0) num = num + smoothing*ngb/neighbors;
				}
				
				// invert the result
				if (num>ZERO) num = 1.0f/num;
				else num = INF;

				mems[x][y][z][k] = num;
				den += num;
			}

			// normalization
			for (k=0;k<classes;k++) {
				mems[x][y][z][k] = mems[x][y][z][k]/den;

				// compute the maximum distance
				dist = Numerics.abs(mems[x][y][z][k]-prev[k]);
				if (dist > distance) distance = dist;
			}
		} else {
			for (k=0;k<classes;k++) 
				mems[x][y][z][k] = 0.0f;
		}

        return distance;
    } // computeMemberships
    
    /**
	 * compute the centroids given the membership functions
	 */
    final public void computeCentroids() {
        int x,y,z,k;
		float num,den;
		
		for (k=0;k<classes;k++) {
			num = 0;
			den = 0;
			for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) if (mask[x][y][z]) {
				num += mems[x][y][z][k]*mems[x][y][z][k]*image[x][y][z];
				den += mems[x][y][z][k]*mems[x][y][z][k];
			}
			if (den>0.0) {
				centroids[k] = num/den;
			} else {
				centroids[k] = 0.0f;
			}
        }
        if (verbose) {
			BasicInfo.displayMessage("centroids: ("+centroids[0]);
			for (k=1;k<classes;k++) BasicInfo.displayMessage(", "+centroids[k]);
			BasicInfo.displayMessage(")\n");
		}   
        return;
    } // computeCentroids
        
   /** 
	 *  compute the FCM membership functions given the centroids everywhere (incl. outside of the mask)
	 */
    final public float computeMembershipsEverywhere() {
        float distance,dist;
        int x,y,z,k,m;
        float den,num;
        float neighbors, ngb;
        
		distance = 0.0f;
		for (x=1;x<nx-1;x++) for (y=1;y<ny-1;y++) for (z=1;z<nz-1;z++) {
			den = 0;
			// remember the previous values
			for (k=0;k<classes;k++) {
				prev[k] = mems[x][y][z][k];
			}
			for (k=0;k<classes;k++) {
					
				// data term
				num = (image[x][y][z]-centroids[k])*(image[x][y][z]-centroids[k]);
					
				// spatial smoothing
				if (smoothing > 0.0f) { 
					ngb = 0.0f;  
					neighbors = 0.0f;
					// 6-connected neighbors
					if (mask[x+1][y][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x+1][y][z][m]*mems[x+1][y][z][m];
						neighbors ++;
					}
					if (mask[x-1][y][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x-1][y][z][m]*mems[x-1][y][z][m];
						neighbors ++;
					}
					if (mask[x][y+1][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y+1][z][m]*mems[x][y+1][z][m];
						neighbors ++;
					}
					if (mask[x][y-1][z]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y-1][z][m]*mems[x][y-1][z][m];
						neighbors ++;
					}
					if (mask[x][y][z+1]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y][z+1][m]*mems[x][y][z+1][m];
						neighbors ++;
					}
					if (mask[x][y][z-1]) for (m=0;m<classes;m++) if (m!=k) {
						ngb += mems[x][y][z-1][m]*mems[x][y][z-1][m];
						neighbors ++;
					}
					if (neighbors>0.0) num = num + smoothing*ngb/neighbors;
				}
				
				// invert the result
				if (num>ZERO) num = 1.0f/num;
				else num = INF;

				mems[x][y][z][k] = num;
				den += num;
			}

			// normalization
			for (k=0;k<classes;k++) {
				mems[x][y][z][k] = mems[x][y][z][k]/den;

				// compute the maximum distance
				dist = Numerics.abs(mems[x][y][z][k]-prev[k]);
				if (dist > distance) distance = dist;
			}
		}

        return distance;
    } // computeMembershipsEverywhere
    
	/** 
	 *	returns the hard classification (max_{clusters}(Mems)) 
	 */
	public final byte[][][] exportHardClassification() {
		int 	x,y,z,t,k,best;
		byte[][][]	classification = new byte[nx][ny][nz];
        float bestmem;
        
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++) {                    
			if (mask[x][y][z]) {
				best = -1;bestmem = 0.0f;
				for (k=0;k<classes;k++) {
					if (mems[x][y][z][k] > bestmem) {
						best = k;
						bestmem = mems[x][y][z][k];
					}
				}   
				if (best==-1)
					classification[x][y][z] = 0;
				else
					classification[x][y][z] = (byte)(best+1);
			} else {
				classification[x][y][z] = 0;
			}
		}
		return classification;
	} // exportHardClassification

	/** 
	 *	export membership functions 
	 */
	public final float[][][][] exportMemberships() {
		int 	x,y,z,k;
		float[][][][]	Mems = new float[nx][ny][nz][classes];
       
		for (x=0;x<nx;x++) for (y=0;y<ny;y++) for (z=0;z<nz;z++)  for (k=0;k<classes;k++) {
			Mems[x][y][z][k] = mems[x][y][z][k];
		}
		return Mems;
	} // exportMemberships

}
