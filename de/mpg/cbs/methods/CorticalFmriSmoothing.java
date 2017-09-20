package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

import org.apache.commons.math3.util.FastMath;

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
 
public class CorticalFmriSmoothing {
	
	private static int[] xoff;
    private static int[] yoff;
    private static int[] zoff;

	// data and membership buffers
	private		float[][][][]		data;
	private		float[][][][]		mapping;
	private		byte[][][]			labeling;
	private		BitSet				mask;				// masking regions not used in computations
	private static	int 		nax,nay,naz;   		// images dimensions
	private static	float 		rax,ray,raz;   		// images resolutions
	private static	int 		nfx,nfy,nfz;   		// images dimensions
	private static	float 		rfx,rfy,rfz;   		// images resolutions
	private static int			nfd;					// 4th dimension
	private		float			sigma;
	private		float			sigmasqr;
	
	// inner variables
	private		BitSet		boundary;
	private		BitSet		next;
	private		BitSet		used;
	private		BitSet		imgbound;
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	
	private	byte				interp;	
	private static final byte NEAREST=10;
	private static final byte LINEAR=20;
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=true;
	
	// new constructor
	public CorticalFmriSmoothing(float[][][][] data_, float[][][][] map_, boolean[][][] msk_, float fwhm_, String interp_,
									int nax_, int nay_, int naz_, float rax_, float ray_, float raz_, 
									int nfx_, int nfy_, int nfz_, int nfd_, float rfx_, float rfy_, float rfz_) {
			
		data = data_;
		mapping = map_;
		
		if (interp_.equals("linear")) {
			interp = LINEAR;
		} else {
			interp = NEAREST;
		}
		
		nax = nax_;
		nay = nay_;
		naz = naz_;
		
		rax = rax_;
		ray = ray_;
		raz = raz_;
				
		nfx = nfx_;
		nfy = nfy_;
		nfz = nfz_;
		nfd = nfd_;
		
		rfx = rfx_;
		rfy = rfy_;
		rfz = rfz_;
				
		// 6-neighborhood: pre-compute the index offsets
		xoff = new int[]{1, -1, 0, 0, 0, 0};
		yoff = new int[]{0, 0, nax, -nax, 0, 0};
		zoff = new int[]{0, 0, 0, 0, nax*nay, -nax*nay};
		
		// create 1D mask
		mask = new BitSet(nax*nay*naz);
		for (int x=0; x<nax; x++) for (int y=0; y<nay; y++) for (int z = 0; z<naz; z++) {
			if (msk_[x][y][z]) mask.set(x+nax*y+nax*nay*z, true);
		}
		
		// basic mask: remove two layers off the images (for avoiding limits)
		imgbound = new BitSet(nax*nay*naz);
		for (int x=0; x<nax; x++) for (int y=0; y<nay; y++) for (int z = 0; z<naz; z++) {
			if (x<=1 || x>=nax-2 || y<=1 || y>=nay-2 || z<=1 || z>=naz-2) {
				mask.set(x+nax*y+nax*nay*z, false);
				imgbound.set(x+nax*y+nax*nay*z, true);
			}
		}
		
		sigma = (fwhm_ / rax_) / (2.0f*(float)Math.sqrt(2.0f*(float)Math.log(2.0f)));
		sigmasqr = sigma*sigma;
		
		boundary = new BitSet(nax*nay*naz);
        next = new BitSet(nax*nay*naz);
        used = new BitSet(nax*nay*naz);
       
        if (debug) BasicInfo.displayMessage("Cortical Region Smoothing (fwhm: "+fwhm_+"mm | sigma: "+sigma+"voxels)");		
        
	}
	
	public CorticalFmriSmoothing(float[][][][] data_, float[][][][] map_, byte[][][] label_, boolean[][][] msk_, float fwhm_, String interp_,
									int nax_, int nay_, int naz_, float rax_, float ray_, float raz_, 
									int nfx_, int nfy_, int nfz_, int nfd_, float rfx_, float rfy_, float rfz_) {
		this(data_, map_, msk_, fwhm_, interp_, nax_, nay_, naz_, rax_, ray_, raz_, nfx_, nfy_, nfz_, nfd_, rfx_, rfy_, rfz_);
		labeling = label_;	
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


	public final void dilateFromPoint(float[] smoothed, int x0, int y0, int z0) {
		
		float dist = 0.0f;
        double[] sdata = new double[nfd];
        double sumWeight = 0.0;
        double gaussWeight = 1.0;
        int x,y,z;
        int xyz0 = x0+nax*y0+nax*nay*z0;
        
        if (!mask.get(xyz0)) return;
        
        boundary.clear();
        used.clear();
        next.clear();
        
        boundary.set(xyz0);
        
        int dmax = Numerics.ceil(3.05f*sigma);
		for (int d=0;d<=dmax;d++) {
			// constant ratio for each shell
			gaussWeight = FastMath.exp(-0.5f*d*d/sigmasqr);
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				// retrieve coordinates
				z = Numerics.floor(xyz/(nax*nay));
				y = Numerics.floor((xyz-z*nax*nay)/nax);
				x = Numerics.floor(xyz-z*nax*nay-y*nax);
				
				// get mapping
				
				// add to average
				sumWeight += gaussWeight;
				for (int n=0;n<nfd;n++) {
					if (interp==LINEAR) sdata[n] += gaussWeight*ImageInterpolation.linearInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
					else sdata[n] += gaussWeight*ImageInterpolation.nearestNeighborInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
				}
				// build the next ring of distances
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask.get(xyzn) && !used.get(xyzn) && !boundary.get(xyzn)) next.set(xyzn);
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
		
		for (int n=0;n<nfd;n++) {
			smoothed[n] = (float)(sdata[n]/sumWeight);
		}
       return;
    }
	
	public final void dilateFromPoint(float[] smoothed, int x0, int y0, int z0, byte lb) {
		
		float dist = 0.0f;
        double[] sdata = new double[nfd];
        double sumWeight = 0.0;
        double gaussWeight = 1.0;
        int x,y,z;
        int xyz0 = x0+nax*y0+nax*nay*z0;
        
        if (!mask.get(xyz0)) return;
        
        boundary.clear();
        used.clear();
        next.clear();
        
        boundary.set(xyz0);
        
        int dmax = Numerics.ceil(3.05f*sigma);
		for (int d=0;d<=dmax;d++) {
			// constant ratio for each shell
			gaussWeight = FastMath.exp(-0.5f*d*d/sigmasqr);
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				// retrieve coordinates
				z = Numerics.floor(xyz/(nax*nay));
				y = Numerics.floor((xyz-z*nax*nay)/nax);
				x = Numerics.floor(xyz-z*nax*nay-y*nax);
				
				// get mapping
				
				// add to average
				sumWeight += gaussWeight;
				for (int n=0;n<nfd;n++) {
					if (interp==LINEAR) sdata[n] += gaussWeight*ImageInterpolation.linearInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
					else sdata[n] += gaussWeight*ImageInterpolation.nearestNeighborInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
				}
				// build the next ring of distances
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (mask.get(xyzn) && !used.get(xyzn) && !boundary.get(xyzn) && labeling[x][y][z]==lb) next.set(xyzn);
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
		
		for (int n=0;n<nfd;n++) {
			smoothed[n] = (float)(sdata[n]/sumWeight);
		}
       return;
    }
	
	public final void dilateVolumetrically(float[] smoothed, int x0, int y0, int z0) {
		
		float dist = 0.0f;
        double[] sdata = new double[nfd];
        double sumWeight = 0.0;
        double gaussWeight = 1.0;
        int x,y,z;
        int xyz0 = x0+nax*y0+nax*nay*z0;
        
        if (!mask.get(xyz0)) return;
        
        boundary.clear();
        next.clear();
        used.clear();
       
        boundary.set(xyz0);
        
        int dmax = Numerics.ceil(3.05f*sigma);
		for (int d=0;d<=dmax;d++) {
			// constant ratio for each shell
			gaussWeight = FastMath.exp(-0.5f*d*d/sigmasqr);
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				// retrieve coordinates
				z = Numerics.floor(xyz/(nax*nay));
				y = Numerics.floor((xyz-z*nax*nay)/nax);
				x = Numerics.floor(xyz-z*nax*nay-y*nax);
				
				// get mapping
				
				// add to average only if inside the mask
				if (mask.get(xyz)) {
					sumWeight += gaussWeight;
					for (int n=0;n<nfd;n++) {
						if (interp==LINEAR) sdata[n] += gaussWeight*ImageInterpolation.linearInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
						else sdata[n] += gaussWeight*ImageInterpolation.nearestNeighborInterpolation(data, 0.0f, mapping[x][y][z][X],mapping[x][y][z][Y],mapping[x][y][z][Z],n, nfx,nfy,nfz,nfd);
					}
				}
				// build the next ring of distances
				for (int k = 0; k<6; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (!used.get(xyzn) && !boundary.get(xyzn) && !imgbound.get(xyzn)) next.set(xyzn);
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
		
		for (int n=0;n<nfd;n++) {
			smoothed[n] = (float)(sdata[n]/sumWeight);
		}
       return;
    }
	
	public final void volumetricSmoothing(float[] smoothed, int x0, int y0, int z0) {
		
		float dist = 0.0f;
        double[] sdata = new double[nfd];
        double sumWeight = 0.0;
        double gaussWeight = 1.0;
        int x,y,z;
        int xyz0 = x0+nax*y0+nax*nay*z0;
        
        if (!mask.get(xyz0)) return;
        
        int dmax = Numerics.ceil(3.05f*sigma);
        for (int i=-dmax;i<=dmax;i++) for (int j=-dmax;j<=dmax;j++) for (int l=-dmax;l<=dmax;l++) {
        	if (mask.get(xyz0+i+nax*j+nax*nay*l)) {
				gaussWeight = FastMath.exp(-0.5f*(i*i+j*j+l*l)/sigmasqr);
				sumWeight += gaussWeight;
				for (int n=0;n<nfd;n++) {
					if (interp==LINEAR) sdata[n] += gaussWeight*ImageInterpolation.linearInterpolation(data, 0.0f, mapping[x0+i][y0+j][z0+l][X],mapping[x0+i][y0+j][z0+l][Y],mapping[x0+i][y0+j][z0+l][Z],n, nfx,nfy,nfz,nfd);
					else sdata[n] += gaussWeight*ImageInterpolation.linearInterpolation(data, 0.0f, mapping[x0+i][y0+j][z0+l][X],mapping[x0+i][y0+j][z0+l][Y],mapping[x0+i][y0+j][z0+l][Z],n, nfx,nfy,nfz,nfd);
				}
			}
		}
		//if (gaussWeight > 0.01) System.out.print("!");
		
		for (int n=0;n<nfd;n++) {
			smoothed[n] = (float)(sdata[n]/sumWeight);
		}
       return;
    }
	
}

