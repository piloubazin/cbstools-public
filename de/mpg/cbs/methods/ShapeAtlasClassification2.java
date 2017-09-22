package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

import Jama.Matrix;

/**
 *
 *  This algorithm performs a simple alignment of the image with a topology template
 *	<p>
 *	The algorithm handles all the little things needed for image cropping, conversion
 * 	to an image buffer, padding, etc.
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class ShapeAtlasClassification2 {
		
	// data buffers
	private 	float[][]		images;  			// original images
	private		int				nc;					// number of images
	private 	int				nix,niy,niz;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		String[]		modality;			// image modality / contrast
	private		float[]			imrange;			// image normalization factor
	private		float[][]		scaling;			// image standard deviation scaling factor
	
	// atlas parameters
	private		SimpleShapeAtlas2	atlas;
	private 	int 				nobj;    	// number of shapes
	private 	int				nax,nay,naz;   			// shape dimensions
	private 	float			rax,ray,raz;   			// shape resolutions
	
	// output results
	private		float[][]		mems;
	private		float[][]		centroids;
	private		float[][]		tissues;
	private		int				nclass;
	private		boolean[]		processed;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// constants
	private static final	float   ISQRT2 = (float)(1.0/Math.sqrt(2.0f));
	private static final	float   ZERO = 1E-20f;
	private static final	float   INF = 1E20f;
		
	public static final	byte BG=0;
	public static final	byte CSF=1;
	public static final	byte GM=2;
	public static final	byte WM=3;
	public static final	byte OUT=4;
	public static final byte nt=5;
		
	
	/**
	 *  constructor
	 */
	public ShapeAtlasClassification2(float[][] img_, String[] mod_, int nc_,
										int nix_, int niy_, int niz_,
										float rix_, float riy_, float riz_,
										SimpleShapeAtlas2 atlas_) {
		images = img_;
		modality = mod_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		imrange = new float[nc];
		
		atlas = atlas_;
		nobj = atlas.getNumber();
		nax = atlas.getShapeDim()[0];
		nay = atlas.getShapeDim()[1];
		naz = atlas.getShapeDim()[2];
		rax = atlas.getShapeRes()[0];
		ray = atlas.getShapeRes()[1];
		raz = atlas.getShapeRes()[2];		
	}

	final public void finalize() {
		images = null;
		System.gc();
	}

	public final int getNclass() { return nclass; }
	
	public final float[][] getMemberships() { return mems; }
	
	public final float[][] getCentroids() { return centroids; }
	
	public final float[] getImageRanges() { return imrange; }
	
	/**
	 *	compute the intensity priors given the atlas
	 */
    public final float[][] atlasIntensityPriors() {
        float num,den;
		float[][] intensity;
		float Imin, Imax;
		
		// get centroids from the atlas
		intensity = atlas.getIntensityPriors(modality, nc);
		
		if (debug) {
			String info = "";
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				info += "image "+c+" intensity priors: ("+intensity[c][0];
				for (int k=1;k<nobj;k++) 
					info += ", "+intensity[c][k];
						
				info += ")\n";
			}
			BasicInfo.displayMessage(info);
		}

		// adjust the centroids to image scale
		for (int c=0;c<nc;c++) {
			/* may be needed for different MP2RAGE implementations
			if (atlas.isQuantitativeContrast(atlas.contrastId(modality[c]))) {
				Imin = 0.0f;
				Imax = 4000.0f; // default setting from our scanner (so it is equivalent on others if the boundaries are different)
				
				BasicInfo.displayMessage("image "+c+" min, max: "+Imin+", "+Imax+"\n");
				
				for (int k=0;k<nobj;k++) {
					intensity[c][k] = Imin + intensity[c][k]*(Imax-Imin);
				}
				imrange[c] = Imax-Imin;
			} else 
			*/
			if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				//Imin = ImageStatistics.robustMinimum(images[c], 0.01f, 3, nix, niy, niz, 5);
				Imin = 0.0f;
				Imax = ImageStatistics.robustMaximum(images[c], 0.0001f, 5, nix, niy, niz, 5);
				
				BasicInfo.displayMessage("image "+c+" min, max: "+Imin+", "+Imax+"\n");
				
				if (atlas.isNormalizedIntensityContrast(atlas.contrastId(modality[c]))) {
					for (int k=0;k<nobj;k++) {
						intensity[c][k] = Imin + intensity[c][k]*(Imax-Imin);
					}
				}
				imrange[c] = Imax-Imin;
			} else {
				imrange[c] = 1.0f;
			}
		}
		if (debug) {
			String info = "";
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				info += "image "+c+" intensity priors: ("+intensity[c][0];
				for (int k=1;k<nobj;k++) 
					info += ", "+intensity[c][k];
						
				info += ")\n";
			}
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				info += "image "+c+" range: ("+imrange[c]+")\n";
			}
			BasicInfo.displayMessage(info);
		}
		return intensity;
    }
    

	/** 
	 *  compute the FCM membership functions for initial tissue types
	 */
    final public void initialAtlasTissueCentroids() {
		tissues = new float[nc][nt];
		scaling = new float[nc][nt];
		float[]	count = new float[nt];
		
		float img,den;
		float[] num = new float[nt];
		
		float[][] intensity = atlasIntensityPriors();
		
		// get intensity mappings : average over all tissue types (and respective sizes from atlas)
		for (int c=0;c<nc;c++) for (int t=0;t<nt;t++) {
			tissues[c][t] = 0.0f;
			scaling[c][t] = 1.0f;
		}
		for (int t=0;t<nt;t++) count[t] = 0.0f;

        	for (int xyz=0;xyz<nax*nay*naz;xyz++) {
			for (int n=0;n<nobj;n++) {
				if (atlas.getObjectType()[n].startsWith("csf")) {
					for (int c=0;c<nc;c++) {
						int contrastId = atlas.contrastId(modality[c]);
						if (atlas.isIntensityContrast(contrastId)) {
							tissues[c][CSF] += atlas.getShapes()[n][xyz]*intensity[c][n];
							scaling[c][CSF] += atlas.getShapes()[n][xyz]*atlas.getMap(contrastId)[n][1]*imrange[c];
						}
					}
					count[CSF]+=atlas.getShapes()[n][xyz];
				} else if (atlas.getObjectType()[n].startsWith("gm")) {
					for (int c=0;c<nc;c++) {
						int contrastId = atlas.contrastId(modality[c]);
						if (atlas.isIntensityContrast(contrastId)) {
							tissues[c][GM] += atlas.getShapes()[n][xyz]*intensity[c][n];
							scaling[c][GM] += atlas.getShapes()[n][xyz]*atlas.getMap(contrastId)[n][1]*imrange[c];
						}
					}
					count[GM]+=atlas.getShapes()[n][xyz];
				} else if (atlas.getObjectType()[n].startsWith("wm")) {
					for (int c=0;c<nc;c++) {
						int contrastId = atlas.contrastId(modality[c]);
						if (atlas.isIntensityContrast(contrastId)) {
							tissues[c][WM] += atlas.getShapes()[n][xyz]*intensity[c][n];
							scaling[c][WM] += atlas.getShapes()[n][xyz]*atlas.getMap(contrastId)[n][1]*imrange[c];
						}
					}
					count[WM]+=atlas.getShapes()[n][xyz];
				} else if (atlas.getObjectType()[n].startsWith("out")) {
					for (int c=0;c<nc;c++) {
						int contrastId = atlas.contrastId(modality[c]);
						if (atlas.isIntensityContrast(contrastId)) {
							tissues[c][OUT] += atlas.getShapes()[n][xyz]*intensity[c][n];
							scaling[c][OUT] += atlas.getShapes()[n][xyz]*atlas.getMap(contrastId)[n][1]*imrange[c];
						}
					}
					count[OUT]+=atlas.getShapes()[n][xyz];
				}
			}
		}
		processed = new boolean[nt];
		for (int t=0;t<nt;t++) {
			if (count[t]>0) {
				processed[t] = true;
				for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
					tissues[c][t] /= count[t];
					scaling[c][t] /= count[t];
				}
			} else {
				processed[t] = false;
			}
		}
		
		// use bg : set to zero
		for (int n=0;n<nobj;n++) {
			if (atlas.getObjectType()[n].startsWith("mask")) {
				processed[BG] = true;
				for (int c=0;c<nc;c++)  if (atlas.isIntensityContrast(atlas.contrastId(modality[c])))
					tissues[c][BG] = 0.0f;
			}
		}
		
		if (verbose) {
			String info = "";
			for (int c=0;c<nc;c++)  if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				info += "image "+c+" tissue priors: ("+tissues[c][0];
				for (int k=1;k<nt;k++) 
					info += ", "+tissues[c][k];
						
				info += ")\n";
			}
			
			BasicInfo.displayMessage(info);
		}
		
		// set as classification
		centroids = new float[nc][nt];
		for (int c=0;c<nc;c++) for (int t=0;t<nt;t++) {
			centroids[c][t] = tissues[c][t];
		}
		nclass = nt;
		mems = new float[nt][nax*nay*naz];
	}

	/** 
	 *  compute the FCM membership functions for initial tissue types
	 */
    final public void computeMemberships() {
    	float den,img;
    	float[] num = new float[nclass];
    	
    	float outdist = estimateOutlierDistance();
    	
		// perform a simple tissue classification with a FCM model on the atlas space
		float[] Xi = new float[3];
        for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
        	
        	atlas.shapeToImageCoordinates(Xi, x, y, z);
        	den = 0;
			for (int k=0;k<nclass;k++) num[k] = 0;
			
			// normalization over image range
			for (int c=0;c<nc;c++)  if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				img =  ImageInterpolation.linearClosestInterpolation(images[c], Xi[0], Xi[1], Xi[2], nix, niy, niz); 
				
				// special case for bg: low variance
				if (processed[BG]) {
					num[BG] += (img-centroids[c][BG])/(0.01f*scaling[c][BG]*imrange[c])*(img-centroids[c][BG])/(0.01f*scaling[c][BG]*imrange[c]);
				}
				for (int k=1;k<nclass-1;k++) if (processed[k]) {
					num[k] += (img-centroids[c][k])/(scaling[c][k]*imrange[c])*(img-centroids[c][k])/(scaling[c][k]*imrange[c]);
				}
			}
			// special case for outliers: constant distance
			if (processed[OUT]) {
				num[OUT] = outdist;
			}
			for (int k=0;k<nclass;k++) if (processed[k]) {
				// invert the result
				if (num[k]>ZERO) num[k] = 1.0f/num[k];
				else num[k] = INF;

				mems[k][x+nax*y+nax*nay*z] = num[k];
				den += num[k];
			}
			
			// normalization
			if (den>0.0f) {
				for (int k=0;k<nclass;k++) {
					mems[k][x+nax*y+nax*nay*z] /= den;
				}
			}
		}
        return;
    }// initialMemberships

    final public float computeCentroids() {
    	float[] den = new float[nclass];
    	float[][] cent = new float[nc][nclass];
    	float img;
		// estimates the new centroids from classification
		for (int k=1;k<nclass;k++) {
			for (int c=0;c<nc;c++)  if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) cent[c][k] = 0.0f;
			den[k] = 0.0f;
		}
		float[] Xi = new float[3];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
        	
        	atlas.shapeToImageCoordinates(Xi, x, y, z);
        	
        	for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				img =  ImageInterpolation.linearClosestInterpolation(images[c], Xi[0], Xi[1], Xi[2], nix, niy, niz); 
					
				for (int k=1;k<nclass;k++) if (processed[k]) {
					cent[c][k] += mems[k][x+nax*y+nax*nay*z]*mems[k][x+nax*y+nax*nay*z]*img;
				}
			}
			for (int k=1;k<nclass;k++) if (processed[k]) {
				den[k] += mems[k][x+nax*y+nax*nay*z]*mems[k][x+nax*y+nax*nay*z];
			}
		}
		float maxdiff = 0.0f;
		for (int k=1;k<nclass;k++) if (processed[k]) {
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				cent[c][k] /= den[k];
				if (Numerics.abs(cent[c][k]-centroids[c][k])/imrange[c] > maxdiff) {
					maxdiff = Numerics.abs(cent[c][k]-centroids[c][k])/imrange[c];
				}
				centroids[c][k] = cent[c][k];
			}
		}
		
		if (debug) {
			String info = "";
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				info += "image "+c+" centroids: ("+centroids[c][0];
				for (int k=1;k<nclass;k++) 
					info += ", "+centroids[c][k];
						
				info += ")\n";
			}
			
			BasicInfo.displayMessage(info);
		}

        return maxdiff;
    }// initialMemberships
    
    final public void estimateIntensityTransform() {
    	String info = "intensity transform: \n";
    	
    	for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
    		float[] coeffs = new float[3];
    	
			double[][] 	tissuepriors = new double[5][3];
			double[][] 	imagevalues = new double[5][1];
			
			// compute a simple quadratic fit?
			for (int k=1;k<=3;k++) if (processed[k]) {
				tissuepriors[k][0] = 1.0;
				tissuepriors[k][1] = tissues[c][k];
				tissuepriors[k][2] = 2.0*tissues[c][k]*tissues[c][k]-1.0;
				imagevalues[k][0] = centroids[c][k];
			}
			// add 0,0 and 1,1 for robustness
			tissuepriors[0][0] = 1.0;
			tissuepriors[0][1] = 0.0;
			tissuepriors[0][2] = -1.0;
			imagevalues[0][0] = 0.0;
			
			tissuepriors[4][0] = 1.0;
			tissuepriors[4][1] = imrange[c];
			tissuepriors[4][2] = 2.0*imrange[c]*imrange[c]-1.0;
			imagevalues[4][0] = imrange[c];
			
			// invert
			//Matrix3D.invert(tissuepriors);
			Matrix matA = new Matrix(tissuepriors);
			Matrix vecb = new Matrix(imagevalues);
			Matrix vecx = matA.solve(vecb);
			
			for (int i=0;i<3;i++) {
				coeffs[i] = (float)vecx.get(i,0);
			}
			
			// corrected priors (sanity check)
			info += "image "+c+" centroids: (";
			for (int k=1;k<nclass;k++) if (processed[k]) {
				float val = coeffs[0] + coeffs[1]*tissues[c][k] + coeffs[2]*(2.0f*tissues[c][k]*tissues[c][k]-1.0f);
				
				info += val+", ";
			}		
			info += ")\n";
			
			info += "from image "+c+" tissue priors: ("+tissues[c][0];
				for (int k=1;k<nt;k++) 
					info += ", "+tissues[c][k];
			info += ")\n";
			
			//BasicInfo.displayMessage(info);
			if (verbose) BasicInfo.displayMessage(info);

			// change the atlas priors afterwards?
			//if (modality[c].equals("T1map") || modality[c].equals("Iso") || modality[c].equals("FlatImg")) {
			if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				int modal = atlas.contrastId(modality[c]);
				for (int n=0;n<nobj;n++) if (!atlas.getObjectType()[n].startsWith("mask")) {
					for (int t=0;t<atlas.getMap(modal)[n].length;t+=3) {
						float val = atlas.getMap(modal)[n][t];
						val = coeffs[0] + coeffs[1]*val + coeffs[2]*(2.0f*val*val-1.0f);
						atlas.setMap(modal,n, t, val);
						float val1 = atlas.getMap(modal)[n][t]
										-0.5f*atlas.getMap(modal)[n][t+1];
						float val2 = atlas.getMap(modal)[n][t]
										+0.5f*atlas.getMap(modal)[n][t+1];
						val1 = coeffs[0] + coeffs[1]*val1 + coeffs[2]*(2.0f*val1*val1-1.0f);
						val2 = coeffs[0] + coeffs[1]*val2 + coeffs[2]*(2.0f*val2*val2-1.0f);
						atlas.setMap(modal, n, t+1, val2-val1);
					}
				}
			}
		}

		if (debug) {
			for (int c=0;c<nc;c++) if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				BasicInfo.displayMessage(atlas.displayMapIntensity(atlas.contrastId(modality[c])));
			}
		}
        return;
    }
    
 	/** 
	 *  compute the average FCM distance among tissues to generate an outlier distance estimate
	 */
    final public float estimateOutlierDistance() {
    	// if no outliers, return a large value
    	if (!processed[OUT]) return INF;
    	
    	float img;
    	float[] num = new float[nclass-1];
    	
    	float avg = 0.0f;
		float den = 0.0f;
		// perform a simple tissue classification with a FCM model on the atlas space
		float[] Xi = new float[3];
        for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
        	
        	atlas.shapeToImageCoordinates(Xi, x, y, z);
        	den = 0;
			for (int k=0;k<nclass-1;k++) num[k] = 0;
			
			// normalization over image range
			for (int c=0;c<nc;c++)  if (atlas.isIntensityContrast(atlas.contrastId(modality[c]))) {
				img =  ImageInterpolation.linearClosestInterpolation(images[c], Xi[0], Xi[1], Xi[2], nix, niy, niz); 
				
				// special case for bg: low variance
				if (processed[BG]) {
					num[BG] += (img-centroids[c][BG])/(0.01f*scaling[c][BG]*imrange[c])*(img-centroids[c][BG])/(0.01f*scaling[c][BG]*imrange[c]);
				} else {
					num[BG] += INF;
				}
				for (int k=1;k<nclass-1;k++) {
					if (processed[k]) {
						num[k] += (img-centroids[c][k])/(scaling[c][k]*imrange[c])*(img-centroids[c][k])/(scaling[c][k]*imrange[c]);
					} else {
						num[k] += INF;
					}
				}				
			}
			// if best is not BG, then accumulate
			byte id = Numerics.argmin(num);
			if (id!=BG) {
				avg += num[id];
				den++;
			}
		}
		if (den>0) avg /= den;
		// set outliers square distance beyond the average best distance
		return 4.0f*avg;
    }

 
    final public float[] getSegmentation(int id) {
    	float[] seg = new float[nax*nay*naz];
    	
		// perform a simple tissue classification with a FCM model on the atlas space
		float max;
		int best;
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
        	best = 0;
        	max = mems[0][xyz];
        	for (int k=1;k<nclass;k++) {
        		if (mems[k][xyz]>max) {
        			best = k;
        			max = mems[k][xyz];
        		}
        	}
        	if (best==id) seg[xyz] = 1.0f;
        	else seg[xyz] = 0.0f;
        }
        return seg;
    }

    final public float[] getSegmentation(int pid, int nid) {
    	float[] seg = new float[nax*nay*naz];
    	
		// perform a simple tissue classification with a FCM model on the atlas space
		float max;
		int best;
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
        	best = 0;
        	max = mems[0][xyz];
        	for (int k=1;k<nclass;k++) {
        		if (mems[k][xyz]>max) {
        			best = k;
        			max = mems[k][xyz];
        		}
        	}
        	if (best==pid) seg[xyz] = 1.0f;
        	else if (best==nid) seg[xyz] = 2.0f;
        	else seg[xyz] = 0.0f;
        }
        return seg;
    }

   final public float[] getDifferentialSegmentation(int pid, int nid) {
    	float[] seg = new float[nax*nay*naz];
    	
		// perform a simple tissue classification with a FCM model on the atlas space
		float max;
		int best;
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
        	best = 0;
        	max = mems[0][xyz];
        	for (int k=1;k<nclass;k++) {
        		if (mems[k][xyz]>max) {
        			best = k;
        			max = mems[k][xyz];
        		}
        	}
        	if (best==pid) seg[xyz] = 2.0f;
        	else if (best==nid) seg[xyz] = 0.0f;
        	else seg[xyz] = 1.0f;
        }
        return seg;
    }

    final public float[] getClassification() {
    	float[] seg = new float[nax*nay*naz];
    	
		// perform a simple tissue classification with a FCM model on the atlas space
		float max;
		int best;
        for (int xyz=0;xyz<nax*nay*naz;xyz++) {
        	best = 0;
        	max = mems[0][xyz];
        	for (int k=1;k<nclass;k++) {
        		if (mems[k][xyz]>max) {
        			best = k;
        			max = mems[k][xyz];
        		}
        	}
        	if (best==WM) seg[xyz] = 3.0f;
        	else if (best==GM) seg[xyz] = 2.0f;
        	else if (best==CSF) seg[xyz] = 1.0f;
        	else seg[xyz] = 0.0f;
        }
        return seg;
    }

    public final float[][] exportMemberships() {
    	float[][] result = new float[nclass][nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		for (int k=0;k<nclass;k++) {
				result[k][x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(mems[k], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
			}
		}
		return result;
    }
    
    public final float[] exportMembership(int k) {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		result[x+nix*y+nix*niy*z] = ImageInterpolation.linearClosestInterpolation(mems[k], Xs[0], Xs[1], Xs[2], nax, nay, naz); 
		}
		return result;
    }
    
    public final float[] exportClassification() {
    	float[] result = new float[nix*niy*niz];
    	
    	float[] Xs = new float[3];
    	float val,best;
    	int id;
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		
    		best = ImageInterpolation.linearInterpolation(mems[0], 0.0f, Xs[0], Xs[1], Xs[2], nax, nay, naz); 
    		id = 0;
    		for (int k=1;k<nclass;k++) {
    			val = ImageInterpolation.linearInterpolation(mems[k], 0.0f, Xs[0], Xs[1], Xs[2], nax, nay, naz);
    			if (val>best) { best = val; id = k; }
    		}
    		result[x+nix*y+nix*niy*z] = id;
		}
		return result;
    }
    
    public final byte[][][] exportClassificationByte() {
    	byte[][][] result = new byte[nix][niy][niz];
    	
    	float[] Xs = new float[3];
    	float val,best;
    	byte id;
        for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
    		atlas.imageToShapeCoordinates(Xs, x, y, z);
    		
    		best = ImageInterpolation.linearInterpolation(mems[0], 0.0f, Xs[0], Xs[1], Xs[2], nax, nay, naz); 
    		id = 0;
    		for (byte k=1;k<nclass;k++) {
    			val = ImageInterpolation.linearInterpolation(mems[k], 0.0f, Xs[0], Xs[1], Xs[2], nax, nay, naz);
    			if (val>best) { best = val; id = k; }
    		}
    		result[x][y][z] = id;
		}
		return result;
    }
}
