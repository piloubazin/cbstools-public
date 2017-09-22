package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm handles basic manipulations of stacked probability maps
 *
 *	@version    July 2016
 *	@author     Pierre-Louis Bazin 
 *		
 *
 */
 
public class MaxProbaRepresentation {
	
	
	// data and membership buffers
	private 	float[][] 		maxproba;  	// maximum probabilities
	private 	byte[][] 		maxlabel;   // corresponding labels
	private static	byte 	  	nmax;					// total number of values to keep
	private 	byte 			nobj;    	// number of objects
	private 	int				nx,ny,nz,nxyz;   		// image dimensions
	
	public MaxProbaRepresentation(byte nmax_, byte nobj_, int nx_, int ny_, int nz_) {
		nmax = nmax_;
		nobj = nobj_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		nxyz = nx*ny*nz;
	}
	
	public final void setMaxProba(float[][] mp_, byte[][] ml_) { maxproba = mp_; maxlabel = ml_; }
	
	public final float[][] getMaxProba() { return maxproba; }
	public final byte[][] getMaxLabel() { return maxlabel; }
	
	// slightly different if background is included or not, if probabilities are normalized or not
	public final void buildFromCompleteProbabilities(float[][] proba) {
	   	maxproba = new float[nmax][nxyz];
    	maxlabel = new byte[nmax][nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
    		boolean zero=true;
    		for (byte n=0;n<nobj && zero;n++) if (proba[n][xyz]>0) zero=false;
    		if (!zero) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[m][xyz]>proba[nbest][xyz]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[nbest][xyz];
					proba[nbest][xyz] = -1.0f;
				}
			} else {
				for (int n=0;n<nmax;n++) {
					maxproba[n][xyz] = -1.0f;
					maxlabel[n][xyz] = -1;
				}
			}
     	}
		proba = null;
	}
	public final void buildFromCompleteProbabilities(float[] proba) {
	   	maxproba = new float[nmax][nxyz];
    		maxlabel = new byte[nmax][nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
    		boolean zero=true;
    		for (byte n=0;n<nobj && zero;n++) if (proba[xyz+nxyz*n]>0) zero=false;
    		if (!zero) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[xyz+nxyz*m]>proba[xyz+nxyz*nbest]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[xyz+nxyz*nbest];
					proba[xyz+nxyz*nbest] = -1.0f;
				}
			} else {
				for (int n=0;n<nmax;n++) {
					maxproba[n][xyz] = -1.0f;
					maxlabel[n][xyz] = -1;
				}
			}
     	}
		proba = null;
	}
	public final void buildFromNormalizedProbabilitiesAndBackground(float[][] proba) {
	   	maxproba = new float[nmax][nx*ny*nz];
    	maxlabel = new byte[nmax][nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			float sum=0.0f;
			for (byte n=0;n<nobj;n++) sum += proba[n][xyz];
			if (sum>0) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[m][xyz]>proba[nbest][xyz]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[nbest][xyz];
					proba[nbest][xyz] = -1.0f;
				}
				for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<1.0f-sum) {
					for (byte m=(byte)(nmax-1);m>n;m--) {
						maxproba[m][xyz] = maxproba[m-1][xyz];
						maxlabel[m][xyz] = maxlabel[m-1][xyz];
					}
					maxproba[n][xyz] = 1.0f-sum;
					maxlabel[n][xyz] = nobj;
					n=nmax;
				}
			} else {
				// background is first, all others are zero
				maxproba[0][xyz] = 1.0f;
				maxlabel[0][xyz] = nobj;
				for (int n=1;n<nmax;n++) {
					maxproba[n][xyz] = 0.0f;
					maxlabel[n][xyz] = -1;
				}
			}
     		}
		// now the background is its own label
		nobj++;
		proba = null;
	}
	public final void buildFromNormalizedProbabilitiesAndBackground(float[] proba) {
	   	maxproba = new float[nmax][nx*ny*nz];
		maxlabel = new byte[nmax][nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) {
    		float sum=0.0f;
    		for (byte n=0;n<nobj;n++) sum += proba[xyz+nxyz*n];
    		if (sum>0) {
			for (byte n=0;n<nmax;n++) {
				byte nbest=0;
				
				for (byte m=1;m<nobj;m++) if (proba[xyz+nxyz*m]>proba[xyz+nxyz*nbest]) {
					nbest = m;
				}
				maxlabel[n][xyz] = nbest;
				maxproba[n][xyz] = proba[xyz+nxyz*nbest];
				proba[xyz+nxyz*nbest] = -1.0f;
			}
			for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<1.0f-sum) {
				for (byte m=(byte)(nmax-1);m>n;m--) {
					maxproba[m][xyz] = maxproba[m-1][xyz];
					maxlabel[m][xyz] = maxlabel[m-1][xyz];
				}
				maxproba[n][xyz] = 1.0f-sum;
				maxlabel[n][xyz] = nobj;
				n=nmax;
			}
		} else {
			// background is first, all others are zero
			maxproba[0][xyz] = 1.0f;
			maxlabel[0][xyz] = nobj;
			for (int n=1;n<nmax;n++) {
				maxproba[n][xyz] = 0.0f;
				maxlabel[n][xyz] = -1;
			}
		}
	}
     	// now the background is its own label
     	nobj++;
		proba = null;
	}
	public final void buildFromCompetingProbabilitiesAndBackground(float[][] proba) {
	   	maxproba = new float[nmax][nx*ny*nz];
    		maxlabel = new byte[nmax][nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			float max=0.0f;
			for (byte n=0;n<nobj;n++) max = Numerics.max(max, proba[n][xyz]);
			if (max>0) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[m][xyz]>proba[nbest][xyz]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[nbest][xyz];
					proba[nbest][xyz] = -1.0f;
				}
				for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<1.0f-max) {
					for (byte m=(byte)(nmax-1);m>n;m--) {
						maxproba[m][xyz] = maxproba[m-1][xyz];
						maxlabel[m][xyz] = maxlabel[m-1][xyz];
					}
					maxproba[n][xyz] = 1.0f-max;
					maxlabel[n][xyz] = nobj;
					n=nmax;
				}
			} else {
				
				// background is first, all others are zero
				maxproba[0][xyz] = 1.0f;
				maxlabel[0][xyz] = nobj;
				for (byte n=1;n<nmax;n++) {
					maxproba[n][xyz] = 0.0f;
					maxlabel[n][xyz] = n;
				}
				/*
				for (int n=0;n<nmax;n++) {
					maxproba[n][xyz] = -1.0f;
					maxlabel[n][xyz] = -1;
				}
				*/
			}
		}
		// now the background is its own label
		nobj++;
		proba = null;
	}
	public final void buildFromCompetingProbabilitiesAndBackground(float[] proba) {
	   	maxproba = new float[nmax][nx*ny*nz];
    		maxlabel = new byte[nmax][nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			float max=0.0f;
			for (byte n=0;n<nobj;n++) max = Numerics.max(max, proba[xyz+nxyz*n]);
			if (max>0) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[xyz+nxyz*m]>proba[xyz+nxyz*nbest]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[xyz+nxyz*nbest];
					proba[xyz+nxyz*nbest] = -1.0f;
				}
				for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<1.0f-max) {
					for (byte m=(byte)(nmax-1);m>n;m--) {
						maxproba[m][xyz] = maxproba[m-1][xyz];
						maxlabel[m][xyz] = maxlabel[m-1][xyz];
					}
					maxproba[n][xyz] = 1.0f-max;
					maxlabel[n][xyz] = nobj;
					n=nmax;
				}
			} else {
				
				// background is first, all others are zero
				maxproba[0][xyz] = 1.0f;
				maxlabel[0][xyz] = nobj;
				for (byte n=1;n<nmax;n++) {
					maxproba[n][xyz] = 0.0f;
					maxlabel[n][xyz] = n;
				}
				/*
				for (int n=0;n<nmax;n++) {
					maxproba[n][xyz] = -1.0f;
					maxlabel[n][xyz] = -1;
				}
				*/
			}
		}
		// now the background is its own label
		nobj++;
		proba = null;
	}
	public final void buildFromCompetingProbabilitiesAndConstantBackground(float[] proba, float bgproba) {
	   	maxproba = new float[nmax][nx*ny*nz];
    		maxlabel = new byte[nmax][nx*ny*nz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			float max=0.0f;
			for (byte n=0;n<nobj;n++) max = Numerics.max(max, proba[xyz+nxyz*n]);
			if (max>0) {
				for (byte n=0;n<nmax;n++) {
					byte nbest=0;
					
					for (byte m=1;m<nobj;m++) if (proba[xyz+nxyz*m]>proba[xyz+nxyz*nbest]) {
						nbest = m;
					}
					maxlabel[n][xyz] = nbest;
					maxproba[n][xyz] = proba[xyz+nxyz*nbest];
					proba[xyz+nxyz*nbest] = -1.0f;
				}
				for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<bgproba) {
					for (byte m=(byte)(nmax-1);m>n;m--) {
						maxproba[m][xyz] = maxproba[m-1][xyz];
						maxlabel[m][xyz] = maxlabel[m-1][xyz];
					}
					maxproba[n][xyz] = bgproba;
					maxlabel[n][xyz] = nobj;
					n=nmax;
				}
			} else {
				
				// background is first, all others are zero
				maxproba[0][xyz] = bgproba;
				maxlabel[0][xyz] = nobj;
				for (byte n=1;n<nmax;n++) {
					maxproba[n][xyz] = 0.0f;
					maxlabel[n][xyz] = n;
				}
				/*
				for (int n=0;n<nmax;n++) {
					maxproba[n][xyz] = -1.0f;
					maxlabel[n][xyz] = -1;
				}
				*/
			}
		}
		// now the background is its own label
		nobj++;
		proba = null;
	}
	
	public final void insertNewLabel(int xyz, float val, byte lb) {
		if (val<maxproba[nmax-1][xyz]) return;
		
		for (byte n=0; n<nmax; n++) if (maxproba[n][xyz]<val) {
			for (byte m=(byte)(nmax-1);m>n;m--) {
				maxproba[m][xyz] = maxproba[m-1][xyz];
				maxlabel[m][xyz] = maxlabel[m-1][xyz];
			}
			maxproba[n][xyz] = val;
			maxlabel[n][xyz] = lb;
			n=nmax;
		}
	}
	
	public final void updateLabel(int xyz, float val, byte lb) {
		if (val<maxproba[nmax-1][xyz]) return;
		
		for (byte n=0; n<nmax; n++) if (maxlabel[n][xyz]==lb) {
			if (val>maxproba[n][xyz]) {
				maxproba[n][xyz] = val;
				for (byte d=0;d<n;d++) if (maxproba[d][xyz]<val) {
					for (byte m=d;m<n;m++) {
						maxproba[m+1][xyz] = maxproba[m][xyz];
						maxlabel[m+1][xyz] = maxlabel[m][xyz];
					}
					maxproba[d][xyz] = val;
					maxlabel[d][xyz] = lb;
					return;
				}
			} else if (val<maxproba[n][xyz]) {
				maxproba[n][xyz] = val;
				for (byte d=(byte)(nmax-1);d>n;d++) if (maxproba[d][xyz]>val) {
					for (byte m=d;m>n;m--) {
						maxproba[m-1][xyz] = maxproba[m][xyz];
						maxlabel[m-1][xyz] = maxlabel[m][xyz];
					}
					maxproba[d][xyz] = val;
					maxlabel[d][xyz] = lb;
					return;
				}
			} 
		}
		insertNewLabel(xyz, val, lb);
		return;
	}
}

