package de.mpg.cbs.libraries;

import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;


/**
 *
 *  This class computes various basic image measures.
 *	
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ImageStatistics {
	
	// no data: used as a library of functions
	
    // simple image functions
    
    public static final byte X = 0;
    public static final byte Y = 1;
    public static final byte Z = 2;
    public static final byte T = 3;

	/**
	 *	mean value of the image
	 */
    public static float mean(float[] img, int nx, int ny, int nz) {
		float mean = 0.0f;
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			mean += img[xyz];
		}
		return mean/(float)(nx*ny*nz);
	}
    
	/**
	 *	mean value of the image
	 */
    public static float mean(float[] img, boolean[] mask, int nx, int ny, int nz) {
		float mean = 0.0f;
		int count = 0;
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (mask[xyz]) {
				mean += img[xyz];
				count++;
			}
		}
		return mean/(float)(count);
	}
    
    public static double mean(double[] img, boolean[] mask, int nx, int ny, int nz) {
		double mean = 0.0;
		int count = 0;
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (mask[xyz]) {
				mean += img[xyz];
				count++;
			}
		}
		return mean/count;
	}
    
    /**
	 *	minimum value of the image
	 */
    public static float minimum(float[][][] img, int nx, int ny, int nz) {
		float min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
    /**
	 *	minimum value of the image
	 */
    public static float minimum(float[][][][] img, int nx, int ny, int nz, int nt) {
		float min = img[0][0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
			if (img[x][y][z][t]<min) min = img[x][y][z][t];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static float maximum(float[][][] img, int nx, int ny, int nz) {
		float max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static float maximum(float[][][][] img, int nx, int ny, int nz, int nt) {
		float max = img[0][0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
			if (img[x][y][z][t]>max) max = img[x][y][z][t];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static float minimum(float[] img, int nx, int ny, int nz) {
		float min = img[0];
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]<min) min = img[xyz];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static float maximum(float[] img, int nx, int ny, int nz) {
		float max = img[0];
        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]>max) max = img[xyz];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static float minimum(float[][][] img, boolean[][][] mask, int nx, int ny, int nz) {
		float min = 0;
		boolean started=false;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
        	if (!started) {
        		min = img[x][y][z];
        		started = true;
        	} else if (img[x][y][z]<min) {
        		min = img[x][y][z];
        	}
		}
		return min;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static float minimum(float[] img, boolean[] mask, int nx, int ny, int nz) {
		float min = 0;
		boolean started=false;
        for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
        	if (!started) {
        		min = img[xyz];
        		started = true;
        	} else if (img[xyz]<min) {
        		min = img[xyz];
        	}
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static float maximum(float[][][] img, boolean[][][] mask, int nx, int ny, int nz) {
		float max = 0;
		boolean started=false;
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
        	if (!started) {
        		max = img[x][y][z];
        		started = true;
        	} else if (img[x][y][z]>max) {
        		max = img[x][y][z];
        	}
		}
		return max;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static float maximum(float[] img, boolean[] mask, int nx, int ny, int nz) {
		float max = 0;
		boolean started=false;
        for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
        	if (!started) {
        		max = img[xyz];
        		started = true;
        	} else if (img[xyz]>max) {
        		max = img[xyz];
        	}
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static int minimum(int[][][] img, int nx, int ny, int nz) {
		int min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static int maximum(int[][][] img, int nx, int ny, int nz) {
		int max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	minimum value of the image
	 */
    public static byte minimum(byte[][][] img, int nx, int ny, int nz) {
		byte min = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]<min) min = img[x][y][z];
		}
		return min;
	}
	
	/**
	 *	maximum value of the image
	 */
    public static byte maximum(byte[][][] img, int nx, int ny, int nz) {
		byte max = img[0][0][0];
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>max) max = img[x][y][z];
		}
		return max;
	}
	
	/**
	 *	normalizes values of the image in [0,1]
	 */
    public static float[] normalize(float[] img, int nx, int ny, int nz, float min, float max) {
		float[] res = new float[nx*ny*nz];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]<min) res[xyz] = min;
			else if (img[xyz]>max) res[xyz] = max;
			else res[xyz] = (img[xyz]-min)/(max-min);
		}
		return res;
	}
	
	/**
	 *	rescale values of the image in [0,1], but allow for values beyond [min,max]
	 */
    public static void rescale(float[] img, int nx, int ny, int nz, float min, float max) {
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			img[xyz] = (img[xyz]-min)/(max-min);
		}
		return;
	}
	
	/**
	 *	rescale values of the image in [0,1], but allow for values beyond [min,max]
	 */
    public static void unscale(float[] img, int nx, int ny, int nz, float min, float max) {
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			img[xyz] = img[xyz]*(max-min) + min;
		}
		return;
	}
	
	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[][][] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0][0][0];
		Imax = image[0][0][0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n>0) ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[][][] image, boolean[][][] mask, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin=0; 
		Imax=0;
		boolean started=false;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			if (!started) {
				Imin = image[x][y][z];
				Imax = image[x][y][z];
				started = true;
			} else {
				if (image[x][y][z]>Imax) Imax = image[x][y][z];
				if (image[x][y][z]<Imin) Imin = image[x][y][z];
			}
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n>0) ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust minimum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[] image, boolean[] mask, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin=0; 
		Imax=0;
		boolean started=false;
		for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
			if (!started) {
				Imin = image[xyz];
				Imax = image[xyz];
				started = true;
			} else {
				if (image[xyz]>Imax) Imax = image[xyz];
				if (image[xyz]<Imin) Imin = image[xyz];
			}
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n>0) ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust minimum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n>0) ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust minimum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points below or equal to the minimum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust minimum
	 *	  @return 			the robust minimum value	
     */
    public static final float robustMinimum(float[] image, float ratio, int scales, int nx, int ny, int nz, int sub) {
		float Imin,Imax,Rmin;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		int xyz;
		
		// ratio: global value
		ratio = ratio*nx/sub*ny/sub*nz/sub;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int x=0;x<nx;x+=sub) for (int y=0;y<ny;y+=sub) for (int z=0;z<nz;z+=sub) {
			xyz = x+nx*y+nx*ny*z;
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmin = Imin;
		
		for (int t=0;t<scales;t++) {
			
			Rmin = Imin;
		
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x+=sub) for (int y=0;y<ny;y+=sub) for (int z=0;z<nz;z+=sub) {
				xyz = x+nx*y+nx*ny*z;
			
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=0;
			while ( (count < ratio) && (n<Nbins) ) {
				count +=bins[n];
				n=n+1;
			}
			Rmin = Imin + (float)(n-0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust minimum: "+Rmin+" ("+n+", "+ratio+", "+count+")\n");
		
			// new boundaries
			float I0 = Imin + (float)(n-1)/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n>0) ratio = ratio - (count-bins[n-1]);		
		}
		
		return Rmin;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[][][] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0][0][0];
		Imax = image[0][0][0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[][][] image, boolean[][][] mask, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = 0;
		Imax = 0;
		boolean started=false;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			if (!started) {
				Imin = image[x][y][z];
				Imax = image[x][y][z];
				started = true;
			} else {
				if (image[x][y][z]>Imax) Imax = image[x][y][z];
				if (image[x][y][z]<Imin) Imin = image[x][y][z];
			}
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
				// first one include both boundaries 
				if (  (image[x][y][z] >= Imin )
					&&(image[x][y][z] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[x][y][z] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[x][y][z] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug 
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[] image, boolean[] mask, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = 0;
		Imax = 0;
		boolean started=false;
		for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
			if (!started) {
				Imin = image[xyz];
				Imax = image[xyz];
				started = true;
			} else {
				if (image[xyz]>Imax) Imax = image[xyz];
				if (image[xyz]<Imin) Imin = image[xyz];
			}
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) if (mask[xyz]) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug 
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[] image, float ratio, int scales, int nx, int ny, int nz ) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		
		// ratio: global value
		ratio = ratio*nx*ny*nz;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nx*ny*nz;xyz++) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Robust maximum estimation
     *    @param 	ratio	float fraction in [0,1]: the minimum number of points above or equal to the maximum over the total volume
	 *    @param 	scales	int: the number of times the scale is refined for finding the robust maximum
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @param	sub			subsampling factor for speed
	 *	  @return 			the robust maximum value	
     */
    public static final float robustMaximum(float[] image, float ratio, int scales, int nx, int ny, int nz, int sub) {
		float Imin,Imax,Rmax;
		int Nbins = 10;
		float[] bins = new float[Nbins];
		float count;
		int n;
		int xyz;
		
		// ratio: global value
		ratio = ratio*nx/sub*ny/sub*nz/sub;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		for (int x=0;x<nx;x+=sub) for (int y=0;y<ny;y+=sub) for (int z=0;z<nz;z+=sub) {
			xyz = x+nx*y+nx*ny*z;
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
		}
		Rmax = Imax;
		
		for (int t=0;t<scales;t++) {
			
			Rmax = Imax;
			
			// compute coarse histogram
			for (n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int x=0;x<nx;x+=sub) for (int y=0;y<ny;y+=sub) for (int z=0;z<nz;z+=sub) {
				xyz = x+nx*y+nx*ny*z;
				
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0f/(float)Nbins*(Imax-Imin) ) ) bins[0]++;
				for (n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (float)n/(float)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (float)(n+1)/(float)Nbins*(Imax-Imin) ) ) bins[n]++;	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = 0;
			n=Nbins;
			while ( (count < ratio) && (n>0) ) {
				n=n-1;
				count +=bins[n];
			}
			Rmax = Imin + (float)(n+0.5f)/(float)Nbins*(Imax-Imin);
			
			//System.out.print("robust maximum: "+Rmax+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			float I0 = Imin + (float)n/(float)Nbins*(Imax-Imin);
			float I1 = Imin + (float)(n+1)/(float)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);			
		}
		
		return Rmax;
	}

	/**
     *    Percentile estimation when samples have weights
     */
    public static final double weightedPercentile(double[] image, double[] weight, double ratio, int nxyz) {
		double Imin,Imax,value;
		int Nbins = 10;
		int scales = 5;
		double[] bins = new double[Nbins];
		double count;
		
		// find first min, max
		Imin = image[0];
		Imax = image[0];
		double Nweight = 0.0;
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (image[xyz]>Imax) Imax = image[xyz];
			if (image[xyz]<Imin) Imin = image[xyz];
			Nweight += weight[xyz];
		}
		// ratio: global value
		ratio = ratio*Nweight/100.0;
		
		value = Imin;
		
		for (int t=0;t<scales;t++) {
			
			value = Imin;
			
			// compute coarse histogram
			for (int n=0;n<Nbins;n++) bins[n] = 0;
			
			for (int xyz=0;xyz<nxyz;xyz++) {
				// first one include both boundaries 
				if (  (image[xyz] >= Imin )
					&&(image[xyz] <= Imin + 1.0/(double)Nbins*(Imax-Imin) ) ) bins[0]+=weight[xyz];
				for (int n=1;n<Nbins;n++) {
					if (  (image[xyz] >  Imin + (double)n/(double)Nbins*(Imax-Imin) )
						&&(image[xyz] <= Imin + (double)(n+1)/(double)Nbins*(Imax-Imin) ) ) bins[n]+=weight[xyz];	
				}
			}
			/* debug
			System.out.print("histogram: \n");
			for (n=0;n<Nbins;n++) System.out.print(" | "+bins[n]);
			System.out.print("\n");
			*/
			
			// find the value corresponding to the ratio
			count = bins[0];
			int n = 0;
			while ( (count < ratio) && (n<Nbins-1) ) {
				n++;
				count +=bins[n];
			}
			value = Imin + (double)(n+0.5)/(double)Nbins*(Imax-Imin);
			
			//System.out.print("estimated value: "+value+" ("+n+", "+ratio+", "+count+")\n");		

			// new boundaries
			double I0 = Imin + (double)n/(double)Nbins*(Imax-Imin);
			double I1 = Imin + (double)(n+1)/(double)Nbins*(Imax-Imin);
			
			Imin = I0;
			Imax = I1;
			
			// new ratio
			if (n<Nbins) ratio = ratio - (count-bins[n]);
		}
		
		return value;
	}

	/**
     *    Robust half gaussian distribution fit
     *    @param 	outliers	boolean: iterates the estimation to account for outliers (uniformly distributed)
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the half gaussian estimated variance	
     */
	public static final float robustHalfGaussianFit(float[][][] image, boolean[][][] mask, boolean outliers, int nx, int ny, int nz) {
		int nb=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) nb++;
		double[] data = new double[nb];
		int n=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			data[n] = image[x][y][z];
			n++;
		}
		return ImageStatistics.robustHalfGaussianFit(data, outliers,nb);
	}
	/*	
    public static final float robustHalfGaussianFit(float[][][] image, boolean[][][] mask, boolean outliers, int nx, int ny, int nz) {
	
    	int nb=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) nb++;
		double[] data = new double[nb];
		float max = 0.0f;
		
		int n=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			data[n] = image[x][y][z];
			if (image[x][y][z]>max) max = image[x][y][z];
			n++;
		}
		Percentile measure = new Percentile();
		measure.setData(data);
		double median = measure.evaluate(50.0);
		
		// med = sigma x sqrt(2) x erf-1(1/2)
		double sigma2 = median*median/0.45493642;
		
		BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
		
		if (outliers) {
			double sigma2prev;
			for (int t=0;t<20;t++) {
				BasicInfo.displayMessage("iteration "+(t+1)+"\n");

				sigma2prev = sigma2;
				
				// compute the global ratio of outliers as the average of outlier probas
				double normg = FastMath.sqrt(2.0/(sigma2*FastMath.PI));
				double ratio = 0.0f;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
					double pg = normg*FastMath.exp(-image[x][y][z]*image[x][y][z]/(2.0*sigma2));
					ratio += (1.0/max/( 1.0/max+pg));	
				}
				ratio /= nb;
				
				BasicInfo.displayMessage("outlier ratio: "+ratio+"\n");
				
				median = measure.evaluate(50.0*(1.0-ratio));
				sigma2 = median*median/0.45493642;
				
				BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
				
				// stopping criterion: less than 1% change in relative value
				if (2.0*Numerics.abs(FastMath.sqrt(sigma2)-FastMath.sqrt(sigma2prev))
					/(FastMath.sqrt(sigma2)+FastMath.sqrt(sigma2prev))<0.01) t=1000;
			}
		}
		return (float)sigma2;
	}
	*/	
	/**
     *    Robust half gaussian distribution fit
     *    @param 	outliers	boolean: iterates the estimation to account for outliers (uniformly distributed)
	 *	  @param	nb	array dimensions
	 *	  @return 			the half gaussian estimated variance	
     */
    public static final float robustHalfGaussianFit(double[] data, boolean outliers, int nb) {
	
    		double max = 0.0;
		for (int n=0;n<nb;n++) if (data[n]>max) max = data[n];
		//if (max==0) return 0.0f;
		
		Percentile measure = new Percentile();
		measure.setData(data, 0, nb);
		double median = measure.evaluate(50.0);
		
		// med = sigma x sqrt(2) x erf-1(1/2)
		double sigma2 = median*median/0.45493642;
		
		//BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
		
		if (outliers) {
			double sigma2prev;
			for (int t=0;t<20 && sigma2>0;t++) {
				//BasicInfo.displayMessage("iteration "+(t+1)+"\n");

				sigma2prev = sigma2;
				
				// compute the global ratio of outliers as the average of outlier probas
				double normg = FastMath.sqrt(2.0/(sigma2*FastMath.PI));
				double ratio = 0.0f;
				for (int n=0;n<nb;n++) {
					double pg = normg*FastMath.exp(-data[n]*data[n]/(2.0*sigma2));
					ratio += (1.0/max/( 1.0/max+pg));	
				}
				ratio /= nb;
				
				if (ratio<0.9) {
					//BasicInfo.displayMessage("outlier ratio: "+ratio+"\n");
					
					median = measure.evaluate(50.0*(1.0-ratio));
					sigma2 = median*median/0.45493642;
					
					//BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
					
					// stopping criterion: less than 1% change in relative value
					if (2.0*Numerics.abs(FastMath.sqrt(sigma2)-FastMath.sqrt(sigma2prev))
						/(FastMath.sqrt(sigma2)+FastMath.sqrt(sigma2prev))<0.01) t=1000;
				} else {
					sigma2 = 0.0;
				}
			}
		}
		return (float)FastMath.sqrt(sigma2);
	}
	
		
	/**
     *    Robust exponential distribution fit
     *    @param 	outliers	boolean: iterates the estimation to account for outliers (uniformly distributed)
	 *	  @param	nx,ny,nz	image dimensions
	 *	  @return 			the exponential estimated beta parameter	
     */
	public static final float robustExponentialFit(float[][][] image, boolean[][][] mask, boolean outliers, int nx, int ny, int nz) {
		int nb=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) nb++;
		double[] data = new double[nb];
		int n=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (mask[x][y][z]) {
			data[n] = image[x][y][z];
			n++;
		}
		return ImageStatistics.robustExponentialFit(data, outliers,nb);
	}
	/**
     *    Robust exponential distribution fit
     *    @param 	outliers	boolean: iterates the estimation to account for outliers (uniformly distributed)
	 *	  @param	nb	array dimensions
	 *	  @return 			the exponential estimated beta parameter	
     */
    public static final float robustExponentialFit(double[] data, boolean outliers, int nb) {
	
    	double max = 0.0;
		for (int n=0;n<nb;n++) if (data[n]>max) max = data[n];
		
		Percentile measure = new Percentile();
		measure.setData(data, 0, nb);
		double median = measure.evaluate(50.0);
		
		// med = beta x log(2)
		double beta = median/FastMath.log(2.0);
		
		//BasicInfo.displayMessage("parameter estimates: beta "+beta+"\n");
		
		if (outliers) {
			double betaprev;
			for (int t=0;t<20 && beta>0;t++) {
				//BasicInfo.displayMessage("iteration "+(t+1)+"\n");

				betaprev = beta;
				
				// compute the global ratio of outliers as the average of outlier probas
				double ratio = 0.0f;
				for (int n=0;n<nb;n++) {
					double pe = FastMath.exp( -data[n]/beta)/beta;
					ratio += (1.0/max/( 1.0/max+pe));	
				}
				ratio /= nb;
				
				//BasicInfo.displayMessage("outlier ratio: "+ratio+"\n");
				if (ratio<0.9) {
					median = measure.evaluate(50.0*(1.0-ratio));
					beta = median/FastMath.log(2.0);
					
					//BasicInfo.displayMessage("parameter estimates: beta "+beta+"\n");
					
					// stopping criterion: less than 1% change in relative value
					if (2.0*Numerics.abs(beta-betaprev)/(beta+betaprev)<0.01) t=1000;
				} else {
					beta = 0.0;
				}
			}
		}
		return (float)beta;
	}
	

}
