package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

/**
 *
 *  This algorithm uses an image-based diffusion technique to regularize images
 *	in regions of low contrast or high noise
 *
 *	@version    June 2017
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class SpatialImageDiffusion {
	
	
	// image data
	private 	float[][]		image;  			// original images
	private 	int				nx,ny,nz, nxyzi; // image dimensions
	//private 	float			rix,riy,riz;   		// image resolutions
	private		float[][]		imvar;				// image intensity scaling
	private		int				nc;					// number of channels
	
	// labeling parameters
	private		float[]			signalproba;		// signal quality probability function
	private		boolean[]		mask;				// masking regions not used in computations
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public SpatialImageDiffusion(float[][] img_, float[][] var_, float[] proba_, byte[] mask_, int nc_,
									int nx_, int ny_, int nz_) {
	
		image = img_;
		imvar = var_;
		signalproba = proba_;
		
		nc = nc_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nxyzi = nx*ny*nz;
		
		// init all the arrays in atlas space
		try {
			//segmentation = new byte[nx*ny*nz];	
			mask = new boolean[nx*ny*nz];	
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial probability decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyzi = x+nx*y+nx*ny*z;
			mask[xyzi] = false;
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				if (mask_==null) {
					for (int c=0;c<nc;c++) if (image[c][xyzi]!=0) mask[xyzi] = true;	
				} else {
					mask[xyzi] = (mask_[xyzi]>0);
				}
			}
		}
	}
	
	public void finalize() {
		//segmentation = null;
	}
	
	public final float[][] getImage() { return image; }
	
	public final float[] getProba() { return signalproba; }
	
     public final void diffuseCertainty(int iter, float scale, float factor, int ngbsize, float mincertainty, float diffratio) {
    	
		//mix with the neighbors?
		float[] certainty = new float[nx*ny*nz];   	
		float[][] newimg = new float[nc][nx*ny*nz];
		float[] newproba = new float[nx*ny*nz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nx*ny*nz][26];   	
		
		// SOR-scheme?
		//float sorfactor = 1.95f;
		float sorfactor = 1.0f;
		
		// compute the functional factor
		//float certaintyfactor = (float)(FastMath.log(0.5)/FastMath.log(factor));
		float certaintyfactor = factor;
		BasicInfo.displayMessage("certainty exponent "+certaintyfactor+"\n");
		
		// rescale the certainty threshold so that it maps to the computed certainty values
		mincertainty = certaintyFunction(mincertainty,certaintyfactor);
		BasicInfo.displayMessage("minimum certainty threshold"+mincertainty+"\n");
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			int xyzi = x+nx*y+nx*ny*z;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
				imgweight[xyzi][j] = sorfactor*diffusionImageWeightFunction(xyzi,xyzj,scale)/ngbsize;
			}
		}
		for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
			for (int c=0;c<nc;c++) newimg[c][xyzi] = image[c][xyzi];
			newproba[xyzi] = signalproba[xyzi];
		}
		
		float maxdiff = 1.0f;
		float meandiff = 1.0f;
		float t0diff = 1.0f;
		for (int t=0;t<iter && meandiff>diffratio*t0diff;t++) {
		//for (int t=0;t<iter && maxdiff>0.1f;t++) {
			BasicInfo.displayMessage("iter "+(t+1));
			maxdiff = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			float nflip = 0.0f;
			float nproc = 0.0f;

			// main loop
			BasicInfo.displayMessage(".");

			// get the gain ; normalize
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (mask[xyzi]) {
				certainty[xyzi] = certaintyFunction(signalproba[xyzi],certaintyfactor);
			}
			// propagate the values : diffusion
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) if (mask[xyzi]) {					
				float den = certainty[xyzi];
				float[] num = new float[nc];
				for (int c=0;c<nc;c++) num[c] = den*image[c][xyzi];
				float nump = den*signalproba[xyzi];
				
				float prev = signalproba[xyzi];
				
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nx, ny, nz);
					ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
				}
				byte[] rank = Numerics.argmax(ngbweight, ngbsize);
				for (int l=0;l<ngbsize;l++) {
					int xyzl = Ngb.neighborIndex(rank[l], xyzi, nx, ny, nz);
					for (int c=0;c<nc;c++) num[c] += ngbweight[rank[l]]*image[c][xyzl];
					nump += ngbweight[rank[l]]*signalproba[xyzl];
					den += ngbweight[rank[l]];
				}
				if (den>1e-9f) {
					for (int c=0;c<nc;c++) num[c] /= den;
					nump /= den;
				}
				for (int c=0;c<nc;c++) newimg[c][xyzi] = num[c];
				newproba[xyzi] = nump;
						
				meandiff += Numerics.abs(nump-prev);
				ndiff++;
				maxdiff = Numerics.max(maxdiff, Numerics.abs(nump-prev));
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int xyzi=0;xyzi<nx*ny*nz;xyzi++) {
				for (int c=0;c<nc;c++) image[c][xyzi] = newimg[c][xyzi];
				signalproba[xyzi] = newproba[xyzi];
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiff+"\n");
		}
		newproba = null;
		newimg = null;
    }
       
    private final float diffusionImageWeightFunction(int xyz, int ngb, float scale) {
    	float maxdiff = 0.0f;
    	for (int c=0;c<nc;c++) {
			float diff = Numerics.abs((image[c][xyz] - image[c][ngb])/imvar[c][xyz]);
			// argmax
			if (diff>maxdiff) maxdiff = diff;
		}
    	return 1.0f/(1.0f+Numerics.square( maxdiff/scale));
    }
    
    private final float certaintyFunction(float delta, float scale) {
    	//return 1.0f - 1.0f/(1.0f+Numerics.square(delta/scale));	
    	//return (float)FastMath.pow(1.0+Numerics.abs(delta), scale);	
    	//return (1.0f+scale)*(1.0f-1.0f/(1.0f+Numerics.abs(delta/scale)));
    	//return (1.0f+Numerics.quad(scale))*(1.0f - 1.0f/(1.0f+Numerics.quad(delta/scale)));	
    	return Numerics.min(delta/scale,1.0f);
    }
    
}

