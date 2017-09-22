package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

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
 
public class MinMaxFiltering {
		
	// data buffers
	private 	float[][]	images;  			// original images
	private		int			nc;
	private 	int			nix,niy,niz;   		// image dimensions
	private 	float		rix,riy,riz;   		// image resolutions
	private		float[]		Imin, Imax;
	
	private		float			scale;				// minimum scale used here
	private		int				nsx,nsy,nsz;		// scaled image dimensions
	private		float[][]		min,max;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
			

	/**
	 *  constructor
	 */
	public MinMaxFiltering(float[][] img_, int nc_, int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_) {
		images = img_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		Imin = new float[nc];
		Imax = new float[nc];
		for (int c=0;c<nc;c++) {
			Imin[c] = ImageStatistics.robustMinimum(images[c], 0.01f, 2, nix, niy, niz, 4);
			Imax[c] = ImageStatistics.robustMaximum(images[c], 0.01f, 2, nix, niy, niz, 4);
		}
		if (debug) System.out.println("image: ["+Imin+", "+Imax+"]");
	}

	/**
	 *  constructor
	 */
	public MinMaxFiltering(float[] img_, int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_) {
		images = new float[1][];
		images[0] = img_;
		nc = 1;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		
		Imin = new float[nc];
		Imax = new float[nc];
		Imin[0] = ImageStatistics.robustMinimum(images[0], 0.001f, 4, nix, niy, niz, 4);
		Imax[0] = ImageStatistics.robustMaximum(images[0], 0.001f, 4, nix, niy, niz, 4);
		if (debug) System.out.println("image: ["+Imin[0]+", "+Imax[0]+"]");
	}

	final public void finalize() {
		images = null;
		System.gc();
	}
	
	/** 
	 *	compute a scaled map of the image with (min, max)
	 */
	public final void buildScaledMaps(float sc, boolean is2d) {
		
		scale = sc;
		float scalez = scale;
		if (is2d) scalez = 1;
		
		nsx = Numerics.floor(nix/scale);
		nsy = Numerics.floor(niy/scale);
		nsz = Numerics.floor(niz/scalez);
		
		if (debug) System.out.println("scaled dimensions: "+nsx+" x "+nsy+" x "+nsz);
		
		min = new float[nc][nsx*nsy*nsz];
		max = new float[nc][nsx*nsy*nsz];
					
		int xyzs,xyzi;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);
			
			for (int c=0;c<nc;c++) {
				min[c][xyzs] = images[c][xyzi];
				max[c][xyzs] = images[c][xyzi];
				
				// check the neighborhood
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					max[c][xyzs] = Numerics.max(max[c][xyzs], images[c][xyzi+i+nix*j+nix*niy*k]);
					min[c][xyzs] = Numerics.min(min[c][xyzs], images[c][xyzi+i+nix*j+nix*niy*k]);
				}
			}
		}
		
		return;
	}
	
	private boolean[] rescale(boolean[] image, float sc, boolean is2d) {
		float scalez = scale;
		float scz = sc;
		if (is2d) { scz = 1.0f; scalez = 1.0f; }
		
		int nrx = Numerics.floor(nix/sc);
		int nry = Numerics.floor(niy/sc);
		int nrz = Numerics.floor(niz/scz);
		
		boolean[] sup = new boolean[nrx*nry*nrz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	int xyzr = Numerics.floor(x*scale/sc) + nrx*Numerics.floor(y*scale/sc) + nrx*nry*Numerics.floor(z*scalez/scz);
			for (int i=0;i<scale/sc;i++) for (int j=0;j<scale/sc;j++) for (int l=0;l<scalez/scz;l++) {
				sup[xyzr+i+nrx*j+nrx*nry*l] = image[x + nsx*y + nsx*nsy*z];
			}
		}
		return sup;
	}

	private byte[] rescale(byte[] image, float sc,boolean is2d) {
		float scalez = scale;
		float scz = sc;
		if (is2d) { scz = 1.0f; scalez = 1.0f; }
		
		int nrx = Numerics.floor(nix/sc);
		int nry = Numerics.floor(niy/sc);
		int nrz = Numerics.floor(niz/scz);
		
		byte[] sup = new byte[nrx*nry*nrz];
			
        for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
        	int xyzr = Numerics.floor(x*scale/sc) + nrx*Numerics.floor(y*scale/sc) + nrx*nry*Numerics.floor(z*scalez/scz);
			for (int i=0;i<scale/sc;i++) for (int j=0;j<scale/sc;j++) for (int l=0;l<scalez/scz;l++) {
				sup[xyzr+i+nrx*j+nrx*nry*l] = image[x + nsx*y + nsx*nsy*z];
			}
		}
		return sup;
	}

	public final float[] growRegion(float[] minval, float[] maxval, float[] rngval, float sc0, int lv, boolean is2d) {
	
		// 1. build the maps at sc0
		buildScaledMaps(sc0,is2d);
		float[][] rng = generateApproxRange();
		
		// 2. find largest region with properties
		boolean[] obj=null,mask;
		for (int c=0;c<nc;c++) {
			boolean[] minobj = ObjectExtraction.objectFromImage(min[c], nsx,nsy,nsz, minval[c]*Imax[c], ObjectExtraction.SUPEQUAL);
			System.out.println("min obj size: "+ObjectStatistics.volume(minobj,nsx,nsy,nsz));
			boolean[] maxobj = ObjectExtraction.objectFromImage(max[c], nsx,nsy,nsz, maxval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("max obj size: "+ObjectStatistics.volume(maxobj,nsx,nsy,nsz));
			boolean[] rngobj = ObjectExtraction.objectFromImage(rng[c], nsx,nsy,nsz, rngval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("rng obj size: "+ObjectStatistics.volume(rngobj,nsx,nsy,nsz));
			
			if (c==0) {
				obj = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
			} else {
				mask = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
				obj = ObjectGeometry.binaryOperation(obj,mask, ObjectGeometry.AND, nsx,nsy,nsz);
			}
		}
		System.out.println("and obj size: "+ObjectStatistics.volume(obj,nsx,nsy,nsz));
		obj = ObjectLabeling.largestObject(obj, nsx,nsy,nsz, 6); 
		obj = ObjectLabeling.removeHoles(obj, nsx,nsy,nsz, 6); 
		
		byte[] label = new byte[nsx*nsy*nsz];
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			if (obj[xyz]) label[xyz] = (byte)scale;
			else label[xyz] = 0;
		}
		
		System.out.println("main obj size: "+ObjectStatistics.volume(obj,nsx,nsy,nsz));
		boolean[] bound = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6); 
		System.out.println("boundary size: "+ObjectStatistics.volume(bound,nsx,nsy,nsz));
		
		// 3. enlarge region at finer scales
		boolean[] neighbor,grown;
		int l=1;
		while (scale>2 && l<lv) {
			obj = rescale(obj, scale/2, is2d);
			bound = rescale(bound, scale/2, is2d);
			label = rescale(label, scale/2, is2d);
			
			// new scale now, nsxyz are different
			buildScaledMaps(scale/2, is2d);
			rng = generateApproxRange();
			
			grown = new boolean[nsx*nsy*nsz];
			// candidates for boundary must be 6-c to object within the larger boundary
			neighbor = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6);
			boolean done = false;
			while (!done) {
				done = true;
				for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
					int xyz = x+nsx*y+nsx*nsy*z;
					if (neighbor[xyz] && bound[xyz]) {
						grown[xyz] = true;
						for (int c=0;c<nc;c++) {
							if (min[c][xyz]<minval[c] || max[c][xyz]>maxval[c] || rng[c][xyz]>rngval[c]) {
								grown[xyz] = false;
							}
						}
						if (grown[xyz]) done = false;
					} else {
						grown[xyz] = obj[xyz];
					}
				}
				obj = grown;
				obj = ObjectLabeling.removeHoles(obj, nsx,nsy,nsz, 6); 
				neighbor = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6);
			}
			bound = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6); 
			for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
				if (label[xyz]==0 && obj[xyz]) label[xyz] = (byte)scale;
			}
			l++;
		}	
		
		// 4. finest scale: image
		//obj = rescale(obj, scale/2);
		//bound = rescale(bound, scale/2);
		float scalez = scale;
		if (is2d) scalez = 1.0f;
		
		grown = new boolean[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyzs = x+nsx*y+nsx*nsy*z;
			int xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);
			
			if (obj[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					grown[ngbi] = true;
				}
			} else {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					grown[ngbi] = false;
				}
			}
		}
		neighbor = ObjectGeometry.objectOutsideBoundary(grown, nix,niy,niz, 6);
		
		float[] result = new float[nix*niy*niz];
		int nobj = 0;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyzs = x+nsx*y+nsx*nsy*z;
			int xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);
			
			if (scale<=2 && bound[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					if (neighbor[xyzi]) {
						result[ngbi] = 1.0f;
						for (int c=0;c<nc;c++) {
							if (images[c][ngbi]<minval[c] || images[c][ngbi]>maxval[c]) {
								result[ngbi] = 0.0f;
							}
						}
						if (result[ngbi]==1) nobj++;
					} else {
						result[ngbi] = 0.0f;
					}
				}
			} else if (obj[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
					nobj++;
				}
			} else {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
				}
			}
		}
		System.out.println("region size: "+nobj);
		return result;
	}
	
	public final float[] growIncreasingRegion(float[] minval, float[] maxval, float[] rngval, float sc0, int lv, boolean is2d) {
	
		// 1. build the maps at sc0
		buildScaledMaps(sc0, is2d);
		float[][] rng = generateApproxRange();
		
		// 2. find largest region with properties
		boolean[] obj=null,mask;
		for (int c=0;c<nc;c++) {
			boolean[] minobj = ObjectExtraction.objectFromImage(min[c], nsx,nsy,nsz, minval[c], ObjectExtraction.SUPEQUAL);
			System.out.println("min obj size: "+ObjectStatistics.volume(minobj,nsx,nsy,nsz));
			boolean[] maxobj = ObjectExtraction.objectFromImage(max[c], nsx,nsy,nsz, maxval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("max obj size: "+ObjectStatistics.volume(maxobj,nsx,nsy,nsz));
			boolean[] rngobj = ObjectExtraction.objectFromImage(rng[c], nsx,nsy,nsz, rngval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("rng obj size: "+ObjectStatistics.volume(rngobj,nsx,nsy,nsz));
			
			if (c==0) {
				obj = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
			} else {
				mask = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
				obj = ObjectGeometry.binaryOperation(obj,mask, ObjectGeometry.AND, nsx,nsy,nsz);
			}
		}
		System.out.println("and obj size: "+ObjectStatistics.volume(obj,nsx,nsy,nsz));
		obj = ObjectLabeling.largestObject(obj, nsx,nsy,nsz, 6); 
		obj = ObjectLabeling.removeHoles(obj, nsx,nsy,nsz, 6); 
		
		byte[] label = new byte[nsx*nsy*nsz];
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			if (obj[xyz]) label[xyz] = (byte)scale;
			else label[xyz] = 0;
		}
		
		System.out.println("main obj size: "+ObjectStatistics.volume(obj,nsx,nsy,nsz));
		boolean[] bound = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6); 
		System.out.println("boundary size: "+ObjectStatistics.volume(bound,nsx,nsy,nsz));
		
		// 3. enlarge region at finer scales
		boolean[] neighbor,grown;
		int l=1;
		while (scale>2 && l<lv) {
			obj = rescale(obj, scale/2, is2d);
			bound = rescale(bound, scale/2, is2d);
			label = rescale(label, scale/2, is2d);
			
			// new scale now, nsxyz are different
			buildScaledMaps(scale/2, is2d);
			rng = generateApproxRange();
			
			grown = new boolean[nsx*nsy*nsz];
			// candidates for boundary must be 6-c to object within the larger boundary
			neighbor = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6);
			boolean done = false;
			while (!done) {
				done = true;
				for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
					int xyz = x+nsx*y+nsx*nsy*z;
					if (neighbor[xyz] && bound[xyz]) {
						grown[xyz] = false;
						// find all possibly previous neighbors ?
						for (byte n=0;n<6 && !grown[xyz];n++) {
							int ngb = Ngb.neighborIndex(n, xyz, nsx, nsy, nsz);
							if (obj[ngb]) {
								grown[xyz] = true;
								for (int c=0;c<nc;c++) {
									if (min[c][xyz]<min[c][ngb] || max[c][xyz]<max[c][ngb]) {
										grown[xyz] = false;
									}
								}
							}
						}
						if (grown[xyz]) done = false;
					} else {
						grown[xyz] = obj[xyz];
					}
				}
				obj = grown;
				obj = ObjectLabeling.removeHoles(obj, nsx,nsy,nsz, 6); 
				neighbor = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6);
			}
			bound = ObjectGeometry.objectOutsideBoundary(obj, nsx,nsy,nsz, 6); 
			for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
				if (label[xyz]==0 && obj[xyz]) label[xyz] = (byte)scale;
			}
			l++;
		}	
		
		// 4. finest scale: image
		//obj = rescale(obj, scale/2);
		//bound = rescale(bound, scale/2);
		float scalez = scale;
		if (is2d) scalez = 1.0f;
		
		grown = new boolean[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyzs = x+nsx*y+nsx*nsy*z;
			int xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);
			
			if (obj[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					grown[ngbi] = true;
				}
			} else {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					grown[ngbi] = false;
				}
			}
		}
		neighbor = ObjectGeometry.objectOutsideBoundary(grown, nix,niy,niz, 6);
		
		float[] result = new float[nix*niy*niz];
		int nobj = 0;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyzs = x+nsx*y+nsx*nsy*z;
			int xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);

			if (scale<=2 && bound[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					if (neighbor[xyzi]) {
						result[ngbi] = 1.0f;
						for (int c=0;c<nc;c++) {
							if (images[c][ngbi]<images[c][xyzi]) {
								result[ngbi] = 0.0f;
							}
						}
						if (result[ngbi]==1) nobj++;
					} else {
						result[ngbi] = 0.0f;
					}
				}
			} else if (obj[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
					nobj++;
				}
			} else {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
				}
			}
		}
		System.out.println("region size: "+nobj);
		return result;
	}
	
	
	public final float[] initialRegion(float[] minval, float[] maxval, float[] rngval, float sc0, int lv, boolean is2d) {
	
		// 1. build the maps at sc0
		buildScaledMaps(sc0, is2d);
		float[][] rng = generateApproxRange();
		
		// 2. find largest region with properties
		boolean[] obj=null,mask;
		for (int c=0;c<nc;c++) {
			boolean[] minobj = ObjectExtraction.objectFromImage(min[c], nsx,nsy,nsz, minval[c]*Imax[c], ObjectExtraction.SUPEQUAL);
			System.out.println("min obj size: "+ObjectStatistics.volume(minobj,nsx,nsy,nsz));
			boolean[] maxobj = ObjectExtraction.objectFromImage(max[c], nsx,nsy,nsz, maxval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("max obj size: "+ObjectStatistics.volume(maxobj,nsx,nsy,nsz));
			boolean[] rngobj = ObjectExtraction.objectFromImage(rng[c], nsx,nsy,nsz, rngval[c]*Imax[c], ObjectExtraction.INFEQUAL);
			System.out.println("rng obj size: "+ObjectStatistics.volume(rngobj,nsx,nsy,nsz));
			
			if (c==0) {
				obj = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
			} else {
				mask = ObjectGeometry.binaryOperation(minobj,maxobj,rngobj, ObjectGeometry.AND, nsx,nsy,nsz);
				obj = ObjectGeometry.binaryOperation(obj,mask, ObjectGeometry.AND, nsx,nsy,nsz);
			}
		}
		System.out.println("and obj size: "+ObjectStatistics.volume(obj,nsx,nsy,nsz));
		obj = ObjectLabeling.largestObject(obj, nsx,nsy,nsz, 6); 
		obj = ObjectLabeling.removeHoles(obj, nsx,nsy,nsz, 6); 
		
		byte[] label = new byte[nsx*nsy*nsz];
		for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			if (obj[xyz]) label[xyz] = (byte)scale;
			else label[xyz] = 0;
		}
		float scalez = scale;
		if (is2d) scalez = scale;
				
		float[] result = new float[nix*niy*niz];
		int nobj = 0;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			int xyzs = x+nsx*y+nsx*nsy*z;
			int xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scalez*z);
			
			if (obj[xyzs]) {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
					nobj++;
				}
			} else {
				for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scalez;k++) {
					int ngbi = xyzi+i+nix*j+nix*niy*k;
					result[ngbi] = label[xyzs];
				}
			}
		}
		System.out.println("region size: "+nobj);
		return result;
	}
	
	
	public final float[][] getMin() { return min; }
		
	public final float[][] getMax() { return max; }
		
	public final float[][] generateApproxImages() {
		float[][] img = new float[nc][nsx*nsy*nsz];
		for (int c=0;c<nc;c++) for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			img[c][xyz] = 0.5f*(max[c][xyz]+min[c][xyz]);
		}
		return img;
	}
	
	public final float[][] generateApproxRange() {
		float[][] img = new float[nc][nsx*nsy*nsz];
		for (int c=0;c<nc;c++) for (int xyz=0;xyz<nsx*nsy*nsz;xyz++) {
			img[c][xyz] = (max[c][xyz]-min[c][xyz]);
		}
		return img;
	}
	
	public final float[] exportMin(int c) {
		int xyzs,xyzi;
		
		float[] imin = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scale*z);
			
			for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scale;k++) {
				imin[xyzi+i+nix*j+nix*niy*k] = min[c][xyzs];
			}
		}
		return imin;
	}
	
	public final float[] exportMax(int c) {
		int xyzs,xyzi;
		
		float[] imax = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scale*z);
			
			for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scale;k++) {
				imax[xyzi+i+nix*j+nix*niy*k] = max[c][xyzs];
			}
		}
		return imax;
	}
	
	public final float[] exportError(int c) {
		int xyzs,xyzi;
		
		float[] ierr = new float[nix*niy*niz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			xyzs = x+nsx*y+nsx*nsy*z;
			xyzi = Numerics.floor(scale*x)+nix*Numerics.floor(scale*y)+nix*niy*Numerics.floor(scale*z);
			
			for (int i=0;i<scale;i++) for (int j=0;j<scale;j++) for (int k=0;k<scale;k++) {
				ierr[xyzi+i+nix*j+nix*niy*k] = images[c][xyzi+i+nix*j+nix*niy*k] - 0.5f*(max[c][xyzs]+min[c][xyzs]);
			}
		}
		return ierr;
	}


}
