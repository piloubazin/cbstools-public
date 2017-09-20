package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

/**
 *
 *  This algorithm grows geodesic regions of interest for cortical sheets
 *
 *	@version    Mar 2012
 *	@author     Christine Lucas Tardif 
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class CorticalRegionSmoothing {
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private		float[]			data;
    private 	float[] 		lvlset;  			// level set functions
	private		boolean[]		mask;				// masking regions not used in computations
	private static	int 		nx,ny,nz;   		// images dimensions
	private static	float 		rx,ry,rz;   		// images resolutions
	private		float			sigma;
	private		float			sigmasqr;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
	// new constructor
	public CorticalRegionSmoothing(float[] data_, float fwhm_, 
			int nx_, int ny_, int nz_,
			float rx_, float ry_, float rz_,
			boolean[] mask_) {
			
		mask = mask_;
		data = data_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		
		rx = rx_;
		ry = ry_;
		rz = rz_;
				
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nx, -nx, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny};
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) mask[x+nx*y+nx*ny*z] = false;
		}
		
		sigma = (fwhm_ / rx_) / (2.0f*(float)Math.sqrt(2.0f*(float)Math.log(2.0f)));
		sigmasqr = sigma*sigma;
		
		if (debug) BasicInfo.displayMessage("Cortical Region Smoothing (fwhm: "+fwhm_+"mm | sigma: "+sigma+"voxels)");		
        
	}
	
		
	public void finalize() {
		data = null;
		mask = null;
		
	}
	
	/**
	 *	clean up the computation arrays
	 */
	public final void cleanUp() {
		System.gc();
	}


	public final float dilateFromPoint(int xyz0) {
		
		float dist = 0.0f;
        double sdata = 0.0;
        double sumWeight = 0.0;
        double gaussWeight = 1.0;
        
        if (!mask[xyz0]) return 0.0f;
        
        BitSet boundary = new BitSet(nx*ny*nz);
        BitSet next = new BitSet(nx*ny*nz);
        BitSet used = new BitSet(nx*ny*nz);
        
        boundary.set(xyz0);
        
        int dmax = Numerics.ceil(3.05f*sigma);
		for (int d=0;d<=dmax;d++) {
			// constant ratio for each shell
			gaussWeight = Math.exp(-0.5f*d*d/sigmasqr);
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				// add to average
				sumWeight += gaussWeight;
				sdata += gaussWeight*data[xyz];
 			
				// build the next ring of distances
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask[xyzn] && !used.get(xyzn)) next.set(xyzn);
				}
				// add to the processed list
				used.set(xyz);
			}
			// replace active boundary
			BitSet tmp = boundary;
			boundary = next;
			next = tmp;
			next.clear();
		}
		//if (gaussWeight > 0.01) System.out.print("!");
		
		sdata /= sumWeight;
		
       return (float)sdata;
    }
	
	
}

