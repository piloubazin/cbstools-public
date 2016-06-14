package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm builds tools for analyzing the intrinsic scale
 *	of (large) images based on various models
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class ScaleAnalysis {
		
	// data buffers
	private 	float[]		image;  			// original images
	private 	int			nix,niy,niz;   		// image dimensions
	private 	float		rix,riy,riz;   		// image resolutions
	private		float		range;
	
	private		int			nsx,nsy,nsz;		// scaled image dimensions
	private		byte[]		map;				// scale map
	private		float[]		min,max;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
			
	/**
	 *  constructor
	 */
	public ScaleAnalysis(float[] img_, int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_) {
		image = img_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		nsx = Numerics.floor(nix/2);
		nsy = Numerics.floor(niy/2);
		nsz = Numerics.floor(niz/2);
		
		map = new byte[nsx*nsy*nsz];
		min = new float[nsx*nsy*nsz];
		max = new float[nsx*nsy*nsz];
					
		float Imin = ImageStatistics.robustMinimum(image, 0.01f, 2, nix, niy, niz, 4);
		float Imax = ImageStatistics.robustMaximum(image, 0.01f, 2, nix, niy, niz, 4);

		range = Imax-Imin;
	}

	final public void finalize() {
		image = null;
		System.gc();
	}
	
	final public byte[] getMap() { return map; }

	/** 
	 *	compute an octree hierarchy based on a threshold of maximum difference
	 * 	(note that many shortcuts are possible here because the min and max have nice scaling properties)
	 */
	public final void octreeMaxDiffScale(float diff, int maxlvl) {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		boolean scaling = false;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			map[xyzs] = 0;
			
			min[xyzs] = image[xyzi];
			max[xyzs] = image[xyzi];
			
			// check the neighborhood
			for (int n=1;n<8;n++) {
				max[xyzs] = Numerics.max(max[xyzs], image[xyzi+nbi[n]]);
				min[xyzs] = Numerics.min(min[xyzs], image[xyzi+nbi[n]]);
			}
			if (max[xyzs]-min[xyzs]<diff*range) {
				map[xyzs] = 1;
				scaling = true;
			}
		}
			
		int step=2;
		byte lv=2;
		int[] nbs = {0,1,nsx,nsx*nsy,1+nsx,nsx+nsx*nsy,nsx*nsy+1,1+nsx+nsx*nsy};
		 
		while (scaling && lv<maxlvl) {
			System.out.print("scale: "+step);
			scaling = false;
			float Imin, Imax;
			boolean process;
			for (int x=0;x<nsx-step;x+=step) for (int y=0;y<nsy-step;y+=step) for (int z=0;z<nsz-step;z+=step) {
				xyzs = x+nsx*y+nsx*nsy*z;
				if (map[xyzs]==lv-1) {
					Imin = min[xyzs];
					Imax = max[xyzs];
					process = true;
					for (int n=1;n<8 && process;n++) {
						if (map[xyzs+nbs[n]*step/2]!=lv-1) {
							process=false;
						} else {
							Imax = Numerics.max(Imax, max[xyzs+nbs[n]*step/2]);	
							Imin = Numerics.min(Imin, min[xyzs+nbs[n]*step/2]);
						}
					}
					if (Imax-Imin<diff*range && process) {
						for (int i=0;i<step;i++) for (int j=0;j<step;j++) for (int k=0;k<step;k++) {
							map[xyzs+i+nsx*j+nsx*nsy*k] = lv;
							min[xyzs+i+nsx*j+nsx*nsy*k] = Imin;
							max[xyzs+i+nsx*j+nsx*nsy*k] = Imax;
						}
						scaling = true;
					}
				}
			}
			step*=2;
			lv++;
		}
		
		return;
	}
	
	/** 
	 *	compute the minimum threshold for each transition of scale
	 */
	public final void scaleMaxDiffOctree(int level) {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			map[xyzs] = 0;
			
			min[xyzs] = image[xyzi];
			max[xyzs] = image[xyzi];
			
			// check the neighborhood
			for (int n=1;n<8;n++) {
				max[xyzs] = Numerics.max(max[xyzs], image[xyzi+nbi[n]]);
				min[xyzs] = Numerics.min(min[xyzs], image[xyzi+nbi[n]]);
			}
			map[xyzs] = 1;
		}
			
		int step=2;
		byte lv=2;
		int[] nbs = {0,1,nsx,nsx*nsy,1+nsx,nsx+nsx*nsy,nsx*nsy+1,1+nsx+nsx*nsy};
		 
		while (lv<=level) {
			System.out.print("scale: "+step);
			float Imin, Imax;
			for (int x=0;x<nsx-step;x+=step) for (int y=0;y<nsy-step;y+=step) for (int z=0;z<nsz-step;z+=step) {
				xyzs = x+nsx*y+nsx*nsy*z;
				if (map[xyzs]==lv-1) {
					Imin = min[xyzs];
					Imax = max[xyzs];
					for (int n=1;n<8;n++) {
						Imax = Numerics.max(Imax, max[xyzs+nbs[n]*step/2]);	
						Imin = Numerics.min(Imin, min[xyzs+nbs[n]*step/2]);
					}
					for (int i=0;i<step;i++) for (int j=0;j<step;j++) for (int k=0;k<step;k++) {
						map[xyzs+i+nsx*j+nsx*nsy*k] = lv;
						min[xyzs+i+nsx*j+nsx*nsy*k] = Imin;
						max[xyzs+i+nsx*j+nsx*nsy*k] = Imax;
					}
				}
			}
			step*=2;
			lv++;
		}
		
		return;
	}
	
	public final float[] generateApproxImage() {
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			min[xyz] = 0.5f*(max[xyz]+min[xyz]);
		}
		return min;
	}
	
	public final float[] exportMap() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		// for temporary output (not optimal, obviously)
		float[] imap = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			for (int n=0;n<8;n++) {
				imap[xyzi+nbi[n]] = map[xyzs];
			}
		}
		return imap;
	}
	
	public final float[] exportMin() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		// for temporary output (not optimal, obviously)
		float[] imin = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			for (int n=0;n<8;n++) {
				imin[xyzi+nbi[n]] = min[xyzs];
			}
		}
		return imin;
	}
	
	public final float[] exportMax() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		// for temporary output (not optimal, obviously)
		float[] imax = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			for (int n=0;n<8;n++) {
				imax[xyzi+nbi[n]] = max[xyzs];
			}
		}
		return imax;
	}
	
	public final float[] exportMinMaxRatio() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};
		
		int xyzs,xyzi;
		// for temporary output (not optimal, obviously)
		float[] iratio = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			for (int n=0;n<8;n++) {
				iratio[xyzi+nbi[n]] = (max[xyzs]-min[xyzs])/range;
			}
		}
		return iratio;
	}
	
	public final float[] exportOctreeApproxImage() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};

		int xyzs,xyzi;
		float[] img = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			if (map[xyzs]==0) {
				for (int n=0;n<8;n++) {
					img[xyzi+nbi[n]] = image[xyzi+nbi[n]];
				}
			} else {
				for (int n=0;n<8;n++) {
					img[xyzi+nbi[n]] = 0.5f*(max[xyzs]+min[xyzs]);
				}
			}
		}
		return img;
	}
	
	public final float[] exportOctreeApproxError() {
		int[] nbi = {0,1,nix,nix*niy,1+nix,nix+nix*niy,nix*niy+1,1+nix+nix*niy};

		int xyzs,xyzi;
		float[] img = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = 2*x+nix*2*y+nix*niy*2*z;
			if (map[xyzs]==0) {
				for (int n=0;n<8;n++) {
					img[xyzi+nbi[n]] = 0.0f;
				}
			} else {
				for (int n=0;n<8;n++) {
					img[xyzi+nbi[n]] = Numerics.abs(0.5f*(max[xyzs]+min[xyzs])-image[xyzi+nbi[n]]);
				}
			}
		}
		return img;
	}
	
	// get average value, RMS errors, etc?
}
