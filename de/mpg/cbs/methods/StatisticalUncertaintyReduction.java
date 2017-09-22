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
	//private 	float			rix,riy,riz;   		// image resolutions
	private		float[]			imscale;			// image intensity scaling
	private		float[][]			imvar;			// image intensity scaling
	private		int				nc;					// number of channels
	
	// labeling parameters
	private 	int 			nobj;    			// number of shapes
	private 	int[]			objlabel;			// label values in the original image
	private		float[][]		bestproba;			// best probability function
	private		byte[][]		bestlabel;			// corresponding labels
	private static	byte 	  	nbest;				// total number of probability functions
	//private 	byte[] 			segmentation;   	// MGDM's segmentation (object indices 1 to N)
	private		boolean[]		mask;				// masking regions not used in computations
	
	// for debug and display
	private static final boolean		debug=false;
	private static final boolean		verbose=true;
	
	/**
	 *  constructors for different cases: with/out outliers, with/out selective constraints
	 */
	public StatisticalUncertaintyReduction(float[][] img_, byte[] mask_, float[] sca_, int nc_,
									int nix_, int niy_, int niz_, 
									//float rix_, float riy_, float riz_,
									byte nbest_) {
		this(img_, null, mask_, sca_, nc_, nix_, niy_, niz_, nbest_);
	}
	public StatisticalUncertaintyReduction(float[][] img_, float[][] var_, byte[] mask_, float[] sca_, int nc_,
									int nix_, int niy_, int niz_, 
									//float rix_, float riy_, float riz_,
									byte nbest_) {
	
		image = img_;
		imscale = sca_;
		imvar = var_;
		nc = nc_;
		nix = nix_;
		niy = niy_;
		niz = niz_;
		nxyzi = nix*niy*niz;
		//rix = rix_;
		//riy = riy_;
		//riz = riz_;
		nbest = nbest_;
		
		// init all the arrays in atlas space
		try {
			//segmentation = new byte[nix*niy*niz];	
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
	
	//public final byte[] getSegmentation() { return segmentation; }
    
	public final void estimateMeanImageNoise() {
		imvar = new float[nc][nxyzi];   	
		
		for (int xyzi=0; xyzi<nxyzi; xyzi++) if (mask[xyzi]) {
			for (int c=0;c<nc;c++) {
				float num = 0.0f;
				float den = 0.0f;
			
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
					if (mask[xyzj]) {
						num += Numerics.abs(image[c][xyzi]-image[c][xyzj]);
						den ++;
					}
				}
				// average of the difference is not quite the standard deviation: need corrective factor
				if (den>0) imvar[c][xyzi] = num/den;
				else imvar[c][xyzi] = 1.0f;
			}
		}
		return;
	}

	public final void estimateMedianImageNoise() {
		imvar = new float[nc][nxyzi];   	
		
		double[] sample = new double[26];
		Percentile measure = new Percentile();
		//double rfactor = Erf.erf(1.0/FastMath.sqrt(2.0));
		double rfactor = FastMath.sqrt(2.0)*Erf.erfInv(0.5);
		for (int xyzi=0; xyzi<nxyzi; xyzi++) if (mask[xyzi]) {
			for (int c=0;c<nc;c++) {
				int ns=0;
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
					if (mask[xyzj]) {
						sample[ns] = Numerics.abs(image[c][xyzi]-image[c][xyzj]);
						ns ++;
					}
				}
				// med = sigma x sqrt(2) x erf-1(1/2); sigma of difference is twice the original
				if (ns>0) imvar[c][xyzi] = (float)(0.5*measure.evaluate(sample, 0, ns, 50.0)/rfactor);
				else imvar[c][xyzi] = 1.0f;
			}
		}
		return;
	}

	public final void estimateImageRankScale(int ngbsize) {
		imvar = new float[nc][nxyzi];   	
		
		float[] sample = new float[26];
		//Percentile measure = new Percentile();
		//double ratio = 100.0*ngbsize/26.0;
		for (int xyzi=0; xyzi<nxyzi; xyzi++) if (mask[xyzi]) {
			for (int c=0;c<nc;c++) {
				int ns=0;
				for (byte j=0;j<26;j++) {
					int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
					if (mask[xyzj]) {
						sample[ns] = Numerics.abs(image[c][xyzi]-image[c][xyzj]);
						ns ++;
					}
				}
				// sets up an arbitrary weighting such that w(diff_ngbsize)=1/2
				//if (ns>0) imvar[c][xyzi] = (float)measure.evaluate(sample, 0, ns, ratio);
				if (ns>=ngbsize) {
					byte[] rank = Numerics.argmin(sample, ns, ngbsize);
					imvar[c][xyzi] = (float)sample[rank[ngbsize-1]];
				} else {
					imvar[c][xyzi] = 1.0f;
				}
			}
		}
		return;
	}
	
	public final float[] computeMaxImageWeight(float scale) {
		float[] imgweight = new float[nix*niy*niz];   	
		
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			imgweight[xyzi] = 0.0f;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				imgweight[xyzi] = Numerics.max(imgweight[xyzi], diffusionImageWeightFunction(xyzi,xyzj,scale));
			}
		}
		return imgweight;
	}
	public final float[] computeMinImageWeight(float scale) {
		float[] imgweight = new float[nix*niy*niz];   	
		
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			imgweight[xyzi] = 1.0f;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				imgweight[xyzi] = Numerics.min(imgweight[xyzi], diffusionImageWeightFunction(xyzi,xyzj,scale));
			}
		}
		return imgweight;
	}
	public final float[][] computeAllImageWeight(float scale) {
		float[][] imgweight = new float[26][nix*niy*niz];   	
		
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				imgweight[j][xyzi] = diffusionImageWeightFunction(xyzi,xyzj,scale);
			}
		}
		return imgweight;
	}
	public final float[][] computeBestImageWeight(float scale, int ngbsize) {
		float[][] imgweight = new float[ngbsize][nix*niy*niz];   	
		
		float[] w0 = new float[26];
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				w0[j] = diffusionImageWeightFunction(xyzi,xyzj,scale);
			}
			byte[] rank = Numerics.argmax(w0, ngbsize);
			for (byte n=0;n<ngbsize;n++) {
				imgweight[n][xyzi] = w0[rank[n]];
			}
		}
		return imgweight;
	}
	public final float[] computeMaxCertainty(float factor) {
		float[] certainty = new float[nix*niy*niz];   	
		
		for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
			certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[1][xyzi],factor);
		}
		return certainty;	
	}
	
    public final void diffuseCertainty(int iter, float scale, float factor, int ngbsize, float mincertainty, boolean computeDistribution, float diffratio) {
    	
		//mix with the neighbors?
		float[] certainty = new float[nix*niy*niz];   	
		float[][] newproba = new float[nbest][nix*niy*niz];
		byte[][] newlabel = new byte[nbest][nix*niy*niz];
		float[] ngbweight = new float[26];
		float[][] imgweight = new float[nix*niy*niz][26];   	
		//byte[][] mapdepth = new byte[nobj][nix*niy*niz];	
		byte[] mapdepth = new byte[nix*niy*niz];	
		
		float[][] objmean = new float[nobj][nc];
		float[][] objvar = new float[nobj][nc];
		float[] objcount = new float[nobj];
		
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
		
		for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
			int xyzi = x+nix*y+nix*niy*z;
			for (byte j=0;j<26;j++) {
				int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
				imgweight[xyzi][j] = sorfactor*diffusionImageWeightFunction(xyzi,xyzj,scale)/ngbsize;
			}
		}
		for (int m=0;m<nbest;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
			newproba[m][xyzi] = bestproba[m][xyzi];
			newlabel[m][xyzi] = bestlabel[m][xyzi];
		}
		/* redo every time?? 
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
		*/
		
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
			/*
			// re-compute depth
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
			*/
			
			/* not a good idea
			// estimate a gaussian intensity distribution for each region, and modulate the corresponding labels
			if (computeDistribution) { // should be done at first, then in last step of the loop
				for (int n=0;n<nobj;n++) {
					for (int c=0;c<nc;c++) {
						objmean[n][c] = 0.0f;
						objvar[n][c] = 0.0f;
					}
					objcount[n] = 0.0f;
				}
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
					int n = bestlabel[0][xyzi];
					if (n>-1) {
						for (int c=0;c<nc;c++) {
							objmean[n][c] += bestproba[0][xyzi]*image[c][xyzi];
						}
						objcount[n]+=bestproba[0][xyzi];
					}
				}
				for (int n=0;n<nobj;n++) if (objcount[n]>0) {
					for (int c=0;c<nc;c++) {
						objmean[n][c] /= objcount[n];
					}
				}
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
					int n = bestlabel[0][xyzi];
					if (n>-1) {
						for (int c=0;c<nc;c++) {
							objvar[n][c] += bestproba[0][xyzi]*(image[c][xyzi]-objmean[n][c])*(image[c][xyzi]-objmean[n][c]);
						}
					}
				}
				for (int n=0;n<nobj;n++) if (objcount[n]>0) {
					for (int c=0;c<nc;c++) {
						objvar[n][c] /= objcount[n];
					}
				}
				// for debug
				for (int c=0;c<nc;c++) {
					BasicInfo.displayMessage("\n contrast "+(c+1));
					BasicInfo.displayMessage("\n mean |");
					for (int n=0;n<nobj;n++) BasicInfo.displayMessage(objmean[n][c]+"| ");
					BasicInfo.displayMessage("\n stdev |");
					for (int n=0;n<nobj;n++) BasicInfo.displayMessage((float)FastMath.sqrt(objvar[n][c])+"| ");
					BasicInfo.displayMessage("\n count |");
					for (int n=0;n<nobj;n++) BasicInfo.displayMessage(objcount[n]+"| ");
				}				
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
					for (int m=0;m<nbest;m++) {
						int n = bestlabel[m][xyzi];
						if (n>-1) {
							//float proba = (float)(1.0/FastMath.sqrt(2.0*FastMath.PI*objvar[n][0])*FastMath.exp(-0.5*(image[0][xyzi]-objmean[n][0])*(image[0][xyzi]-objmean[n][0])/objvar[n][0]));
							float proba = (float)FastMath.exp(-0.5*(image[0][xyzi]-objmean[n][0])*(image[0][xyzi]-objmean[n][0])/objvar[n][0]);
							for (int c=1;c<nc;c++) {
								//float probac = (float)(1.0/FastMath.sqrt(2.0*FastMath.PI*objvar[n][c])*FastMath.exp(-0.5*(image[c][xyzi]-objmean[n][c])*(image[c][xyzi]-objmean[n][c])/objvar[n][c]));
								float probac = (float)FastMath.exp(-0.5*(image[c][xyzi]-objmean[n][c])*(image[c][xyzi]-objmean[n][c])/objvar[n][c]);
								if (probac<proba) proba = probac;
							}
							bestproba[m][xyzi] *= proba;
						}
					}
				}
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
			*/
				
			// main loop: label-per-label
			for (byte n=0;n<nobj;n++) {
				
				// re-compute depth
				// mapdepth only needed for obj n: single map, recomputed every time?
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
					mapdepth[xyzi] = nbest;
					for (byte m=0;m<nbest;m++) {
						if (bestlabel[m][xyzi]==n) {
							mapdepth[xyzi] = m;
						}
					}
				}
				
				//BasicInfo.displayMessage("propagate gain for label "+n+"\n");
				BasicInfo.displayMessage(".");

				// get the gain ; normalize
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
					if (mapdepth[xyzi]<nbest) {
						if (mapdepth[xyzi]==0) certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[1][xyzi],certaintyfactor);
						else certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[mapdepth[xyzi]][xyzi],certaintyfactor);
					} else {
						certainty[xyzi] = certaintyFunction(bestproba[0][xyzi]-bestproba[nbest-1][xyzi],certaintyfactor);
					}
				}
				// propagate the values : diffusion
				for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
				//for (int x=1;x<nix-1;x++) for (int y=1;y<niy-1;y++) for (int z=1;z<niz-1;z++) {
					//int xyzi = x+nix*y+nix*niy*z;
					if (mapdepth[xyzi]<nbest && certainty[xyzi]<=mincertainty) {
						nproc++;
						
						float den = certainty[xyzi];
						float num = den*bestproba[mapdepth[xyzi]][xyzi];
						float prev = bestproba[mapdepth[xyzi]][xyzi];
					
						for (byte j=0;j<26;j++) {
							int xyzj = Ngb.neighborIndex(j, xyzi, nix, niy, niz);
							ngbweight[j] = certainty[xyzj]*imgweight[xyzi][j];
						}
						byte[] rank = Numerics.argmax(ngbweight, ngbsize);
						for (int l=0;l<ngbsize;l++) {
							int xyzl = Ngb.neighborIndex(rank[l], xyzi, nix, niy, niz);
							if (mapdepth[xyzl]<nbest) num += ngbweight[rank[l]]*bestproba[mapdepth[xyzl]][xyzl];
							else num += ngbweight[rank[l]]*bestproba[nbest-1][xyzl];
							den += ngbweight[rank[l]];
						}
						// other forms of regularization??
						if (computeDistribution) {
							//num /= 0.5f*(1.0f+0.5f*1.95f); // about 1.0
						} else {
							if (den>1e-9f) num /= den;
						}
						newproba[mapdepth[xyzi]][xyzi] = num;
						newlabel[mapdepth[xyzi]][xyzi] = n;
						
						meandiff += Numerics.abs(num-prev);
						ndiff++;
						maxdiff = Numerics.max(maxdiff, Numerics.abs(num-prev));
						if (prev<0.5f && num>0.5f) nflip++;
						if (prev>0.5f && num<0.5f) nflip++;
					}
				}
				// maybe not really useful..
				/*
				if (computeDistribution) {
					float maxproba = 0.0f;
					for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
						if (mapdepth[xyzi]<nbest && certainty[xyzi]<=mincertainty) {
							if (newproba[mapdepth[xyzi]][xyzi]>maxproba) maxproba = newproba[mapdepth[xyzi]][xyzi];
						}
					}
					meandiff = 0.0f;
					maxdiff = 0.0f;
					for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) if (mask[xyzi]) {
						if (mapdepth[xyzi]<nbest && certainty[xyzi]<=mincertainty) {
							newproba[mapdepth[xyzi]][xyzi] /= maxproba;
							meandiff += Numerics.abs(newproba[mapdepth[xyzi]][xyzi]-bestproba[mapdepth[xyzi]][xyzi]);
							maxdiff = Numerics.max(maxdiff, Numerics.abs(newproba[mapdepth[xyzi]][xyzi]-bestproba[mapdepth[xyzi]][xyzi]));
						}
					}
				}
				*/
			}
			if (ndiff>0) meandiff /= ndiff;
			if (t==0) t0diff = meandiff;
			
			// make a hard copy
			for (int m=0;m<nbest;m++) for (int xyzi=0;xyzi<nix*niy*niz;xyzi++) {
				bestproba[m][xyzi] = newproba[m][xyzi];
				bestlabel[m][xyzi] = newlabel[m][xyzi];
			}
			float nresort=0.0f;
			float nreseg=0.0f;
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
							nresort++;
							if (l==1) nreseg++;
						}
					}
					//if (stop) m=nbest;
				}
			}
			BasicInfo.displayMessage("mean diff. "+meandiff+", max diff. "+maxdiff+", changed. "+nreseg+"\n");
			//BasicInfo.displayMessage("n processed "+nproc+", n flipped "+nflip+"\n");
			
			//BasicInfo.displayMessage("n resorted"+nresort+"\n");
			
		}
		newproba = null;
		newlabel = null;
		/* every time?
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
		*/
    }
       
    private final float diffusionImageWeightFunction(int xyz, int ngb, float scale) {
    	float maxdiff = 0.0f;
    	if (imvar==null) {
		for (int c=0;c<nc;c++) {
			float diff = Numerics.abs((image[c][xyz] - image[c][ngb])/imscale[c]);
			// argmax
			if (diff>maxdiff) maxdiff = diff;
			// arg max dist from 1
			/*
			if (c==0) maxdiff = diff;
			else if (maxdiff>1 && diff>maxdiff) maxdiff = diff;
			else if (maxdiff<1 && diff<maxdiff) maxdiff = diff;
			else if (maxdiff>1 && diff<1.0f/maxdiff) maxdiff = diff;
			else if (maxdiff<1 && diff>1.0f/maxdiff) maxdiff = diff;
			*/
		}
	} else {
		for (int c=0;c<nc;c++) {
			float diff = Numerics.abs((image[c][xyz] - image[c][ngb])/imvar[c][xyz]);
			// argmax
			if (diff>maxdiff) maxdiff = diff;
			// arg max dist from 1
			/*
			if (c==0) maxdiff = diff;
			else if (maxdiff>1 && diff>maxdiff) maxdiff = diff;
			else if (maxdiff<1 && diff<maxdiff) maxdiff = diff;
			else if (maxdiff>1 && diff<1.0f/maxdiff) maxdiff = diff;
			else if (maxdiff<1 && diff>1.0f/maxdiff) maxdiff = diff;
			*/
		}
	}
    	return 1.0f/(1.0f+Numerics.square( maxdiff/scale));
    }
    
    private final float certaintyFunction(float delta, float scale) {
    	//return 1.0f - 1.0f/(1.0f+Numerics.square(delta/scale));	
    	return (float)FastMath.pow(1.0+Numerics.abs(delta), scale);	
    	//return (1.0f+scale)*(1.0f-1.0f/(1.0f+Numerics.abs(delta/scale)));
    	//return (1.0f+Numerics.quad(scale))*(1.0f - 1.0f/(1.0f+Numerics.quad(delta/scale)));	
    	
    }
    
}

