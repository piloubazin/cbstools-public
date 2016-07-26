package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm uses an image-based diffusion technique to reduce uncertainty
 *	in multi-atlas or statistical labelings
 *
 *	@version    August 2016
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class StatisticalUncertaintyReduction {
	
	
	// image data
	private 	float[][]		image;  			// original images
	private 	int				nix,niy,niz, nxyzi;   		// image dimensions
	private 	float			rix,riy,riz;   		// image resolutions
	private		boolean[]		imused;				// check if image modality / contrast is used
	private		float[]			imscale;			// image intensity scaling
	private		int				nc;					// number of channels

	// labeling parameters
	private 	int 			nobj;    			// number of shapes
	private 	int[]			objlabel;			// label values in the original image
	private		float[][]		bestproba;			// best probability function
	private		byte[][]		bestlabel;			// corresponding labels
	private static	byte 	  	nbest;				// total number of probability functions
	private 	byte[] 			segmentation;   	// MGDM's segmentation (object indices 1 to N)
	private		boolean[]		mask;				// masking regions not used in computations
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public StatisticalUncertaintyReduction(float[][] img_, boolean[] used_, float[] sca_, int nc_,
									int nix_, int niy_, int niz_, float rix_, float riy_, float riz_,
									byte nbest_) {
	
		image = img_;
		imscale = sca_;
		imused = used_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		nxyzi = nix*niy*niz;
		rix = rix_;
		riy = riy_;
		riz = riz_;
		nbest = nbest_;
		
		// init all the arrays in atlas space
		try {
			segmentation = new byte[nix*niy*niz];	
			mask = new boolean[nix*niy*niz];	
		} catch (OutOfMemoryError e){
			 finalize();
			System.out.println(e.getMessage());
			return;
		}
		if (debug) BasicInfo.displayMessage("initial probability decomposition\n");		
		
		// basic mask: remove two layers off the images (for avoiding limits)
		for (int x=0; x<nix; x++) for (int y=0; y<niy; y++) for (int z = 0; z<niz; z++) {
			int xyzi = x+nix*y+nix*niy*z;
			mask[xyzi] = false;
			if (x>0 && x<nix-1 && y>0 && y<niy-1 && z>0 && z<niz-1) {
				for (int c=0;c<nc;c++) if (image[c][xyzi]!=0) mask[xyzi] = true;	
			}
		}
	}
	
	public void finalize() {
		segmentation = null;
	}
	
	public final void setBestProbabilities(float[][] bestpb, byte[][] bestlb, int[] objlb) {
		bestproba = bestpb;
		bestlabel = bestlb;
		nbest = Numerics.min(nbest, (byte)bestproba.length);
		objlabel = objlb;
		nobj = objlabel.length;
	}
	
	public final void initBestProbaFromSegmentation(int[] seg, int nlb, float dist) {
		bestlabel = new byte[nbest][nxyzi];
		nobj = nlb;
		objlabel = new int[nobj];
		
		// set the first labels
		byte nset=0;
		for (int xyzi=0; xyzi<nxyzi; xyzi++) if (mask[xyzi]) {
			byte lb=-1;
			for (byte n=0;n<nset && lb==-1;n++) if (seg[xyzi]==objlabel[n]) lb=n;
			if (lb==-1) {
				lb = nset;
				objlabel[nset] = seg[xyzi];
				nset++;
			}
			bestlabel[0][xyzi] = lb;
		}
		
		// propagate the distances MGDM-style
		// ...
	}
	
	public final float[][] getProbabilities() { return bestproba; }
	
	public final byte[][] getLabels() { return bestlabel; }
	
	public final byte[] getSegmentation() { return segmentation; }
    

    public final void diffuseCertainty(int iter, float scale, float factor, int ngbsize, float mincertainty) {
    	
		//mix with the neighbors?
		float[] certainty = new float[nix*niy*niz];   	
		float[][] newproba = new float[nbest][nix*niy*niz];
		byte[][] newlabel = new byte[nbest][nix*niy*niz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nix*niy*niz][26];   	
		byte[][] mapdepth = new byte[nobj][nix*niy*niz];	
		
		// rescale the certainty threshold so that it maps to the computed certainty values
		mincertainty = certaintyFunction(mincertainty,factor);
		
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				imgweight[xyzi][j] = diffusionImageWeightFunction(xyzi,xyzj,scale)/ngbsize;
			}
		}
		for (int m=0;m<nbest;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			newproba[m][xyzi] = bestproba[m][xyzi];
			newlabel[m][xyzi] = bestlabel[m][xyzi];
		}
		for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			for (byte n=0;n<nobj;n++) {
				mapdepth[n][xyzi] = nbest;
				for (byte m=0;m<nbest;m++) {
					if (bestlabel[m][xyzi]==n) {
						mapdepth[n][xyzi] = m;
					}
				}
			}
		}
		float maxdiff = 1.0f;
		float meandiff = 1.0f;
		for (int t=0;t<iter && meandiff>0.001f;t++) {
			BasicInfo.displayMessage("iter "+(t+1));
			maxdiff = 0.0f;
			meandiff = 0.0f;
			float ndiff = 0.0f;
			for (byte n=0;n<nobj;n++) {
				//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
				BasicInfo.displayMessage(".");
				// get the gain ; normalize
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
					if (mapdepth[n][xyzi]<nbest) {
						if (mapdepth[n][xyzi]==0) certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[1][xyzi],factor);
						else certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[mapdepth[n][xyzi]][xyzi],factor);
					} else {
						certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[nbest-1][xyzi],factor);
					}
				}
				// propagate the values : diffusion
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
				//for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
					//int xyzi = x+nix*y+nix*niy*z;
					if (mapdepth[n][xyzi]<nbest && certainty[xyzi]<=mincertainty) {
						float den = certainty[xyzi];
						float num = den*bestproba[mapdepth[n][xyzi]][xyzi];
						float prev = bestproba[mapdepth[n][xyzi]][xyzi];
					
						for (byte j=0;j<26;j++) {
							int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
							ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
						}
						byte[] rank = Numerics.argmax(ngbweight, ngbsize);
						for (int l=0;l<ngbsize;l++) {
							int xyzl = Ngb.neighborIndex(rank[l], xyzi, nix, niy, niz);
							if (mapdepth[n][xyzl]<nbest) num += ngbweight[rank[l]]*bestproba[mapdepth[n][xyzl]][xyzl];
							else num += ngbweight[rank[l]]*bestproba[nbest-1][xyzl];
							den += ngbweight[rank[l]];
						}
						num /= den;
						
						newproba[mapdepth[n][xyzi]][xyzi] = num;
						
						meandiff += Numerics.abs(num-prev);
						ndiff++;
						maxdiff = Numerics.max(maxdiff, Numerics.abs(num-prev));
					}
				}
			}
			meandiff /= ndiff;
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiff+"\n");
			// make a hard copy
			for (int m=0;m<nbest+1;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
				bestproba[m][xyzi] = newproba[m][xyzi];
				bestlabel[m][xyzi] = newlabel[m][xyzi];
			}
		}
		newproba = null;
		newlabel = null;
		
		// re-sort the gain functions
		for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			for (int m=0;m<nbest-1;m++) {
				boolean stop=true;
				for (int l=nbest-1;l>m;l--) {
					if (bestproba[l][xyzi]>bestproba[l-1][xyzi]) {
						float swap = bestproba[l-1][xyzi];
						bestproba[l-1][xyzi] = bestproba[l][xyzi];
						bestproba[l][xyzi] = swap;
						byte swaplb = bestlabel[l-1][xyzi];
						bestlabel[l-1][xyzi] = bestlabel[l][xyzi];
						bestlabel[l][xyzi] = swaplb;
						stop=false;
					}
				}
				if (stop) m=nbest;
			}
		}
    }
    
    private final float diffusionImageWeightFunction(int xyz, int ngb, float scale) {
    	float maxdiff = 0.0f;
    	for (int c=0;c<nc;c++) if (imused[c]) {
    		float diff = Numerics.abs(image[c][xyz]/imscale[c] - image[c][ngb]/imscale[c]);
    		if (diff>maxdiff) maxdiff = diff;
    	}
    	return 1.0f/(1.0f+Numerics.square( maxdiff/scale));
    }
    
    private final float certaintyFunction(float delta, float scale) {
    	return 1.0f - 1.0f/(1.0f+Numerics.square(delta/scale));	
    }
    
}

