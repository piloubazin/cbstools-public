package de.mpg.cbs.libraries;


import java.io.*;
import java.util.*;
import java.lang.*;

import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *   Centralized generic functions for spatial transformation of coordinates
 *
 *
 *	@version    Feb 2009
 *	@author     Pierre-Louis Bazin
 *		
 *
*/
public class ParametricTransform {
	
	int type;
	private static final	int   	NONE = 0;
	private static final	int   	RIGID = 11;
	private static final	int   	SINGLE_SCALE = 12;
	private static final	int		SCALED_RIGID = 13;
	private static final	int   	FULLY_AFFINE = 31;
	
	float 	x0i, y0i, z0i;	// image center
	float 	rix, riy, riz;	// image resolutions
	int 	nix, niy, niz;	// image dimensions
	
	float 	x0t, y0t, z0t;	// template center
	float 	rtx, rty, rtz;	// template resolutions
	int 	ntx, nty, ntz;	// template dimensions
	
	int 	Nt;				// transform dimensions
	
	float Iscale, Iscale2;
	
	boolean debug = true;
	
	public ParametricTransform(String type_, 
								float x0i_, float y0i_, float z0i_,
								float rix_, float riy_, float riz_,
								int nix_, int niy_, int niz_,
								float x0t_, float y0t_, float z0t_,
								float rtx_, float rty_, float rtz_,
								int ntx_, int nty_, int ntz_) {
		
		if (type_.equals("rigid")) { type = RIGID; Nt = 6; } 
		else if (type_.equals("single_scale")) { type = SINGLE_SCALE; Nt = 7; }
		else if (type_.equals("scaled_rigid")) { type = SCALED_RIGID; Nt = 9; }
		else if (type_.equals("fully_affine")) { type = FULLY_AFFINE; Nt = 12; }
		else { type = NONE; Nt = 0; }
		
		if (debug) BasicInfo.displayMessage("transform type: "+type_+"("+type+")\n");
		
		x0i = x0i_; y0i = y0i_; z0i = z0i_;
		rix = rix_; riy = riy_; riz = riz_;
		nix = nix_; niy = niy_; niz = niz_;
		
		x0t = x0t_; y0t = y0t_; z0t = z0t_;
		rtx = rtx_; rty = rty_; rtz = rtz_;
		ntx = ntx_; nty = nty_; ntz = ntz_;
		
		Iscale2 = 0.25f*( nix*rix*nix*rix + niy*riy*niy*riy + niz*riz*niz*riz );
		Iscale = (float)Math.sqrt(Iscale2);	
	}
	
	public final void finalize() {
	}
	
	public final int getDimension() { return Nt; }
	
	public final String getTransformType() {
		if (type==NONE) 						return "none";
		else if (type==RIGID) 					return "rigid";
		else if (type==SINGLE_SCALE) 			return "single_scale";
		else if (type==SCALED_RIGID) 			return "scaled_rigid";
		else if (type==FULLY_AFFINE) 			return "fully_affine";
		else 									return "not found!";
	}		
	
	/** check if the transform is linear in the image coordinates (and thus pre-computable) */
	public final boolean isLinear() {
		if (type==NONE || type==RIGID || type==SINGLE_SCALE || type==SCALED_RIGID || type==FULLY_AFFINE) return true;
		else return false;
	}
	
	/** check if the transform uses a rotation matrix */
	public final boolean useRotation() {
		if (type==RIGID || type==SINGLE_SCALE || type==SCALED_RIGID) return true;
		else return false;
	}
	
	/** 
	 *	computes the transformed coordinates from image to template space
	 *	with scaling (2: half res, 4: quarter res, etc)
	 */
	public final void imageToTemplate(float[] X, float x,float y,float z, float[] trans, float[][] rot, float scale) {
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		if (type==NONE) {
			X[0] = ( Xi/rtx + x0t)/scale;
			X[1] = ( Yi/rty + y0t)/scale;
			X[2] = ( Zi/rtz + z0t)/scale;
		} else if (type==RIGID) {
			X[0] = ( (rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( (rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( (rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			
			X[0] = ( factor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( factor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( factor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		} else if (type==SCALED_RIGID) {
			X[0] = ( (1.0f+trans[6]/Iscale)*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx + x0t)/scale;
			X[1] = ( (1.0f+trans[7]/Iscale)*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty + y0t)/scale;
			X[2] = ( (1.0f+trans[8]/Iscale)*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz + z0t)/scale;
		
		} else if (type==FULLY_AFFINE) {
			X[0] = ( ((1+trans[0])*Xi+   trans[1] *Yi+   trans[2] *Zi+trans[ 9])/rtx + x0t)/scale;
			X[1] = ( (   trans[3] *Xi+(1+trans[4])*Yi+   trans[5] *Zi+trans[10])/rty + y0t)/scale;
			X[2] = ( (   trans[6] *Xi+   trans[7] *Yi+(1+trans[8])*Zi+trans[11])/rtz + z0t)/scale;
		}
		return;
	}

	/** 
	 *	computes the transformed coordinates from image to template space
	 *	with scaling (2: half res, 4: quarter res, etc)
	 */
	public final void templateToImage(float[] X, float x,float y,float z, float[] trans, float[][] rot, float scale) {
		float Xt = (x*scale-x0t)*rtx;
		float Yt = (y*scale-y0t)*rty;
		float Zt = (z*scale-z0t)*rtz;
		if (type==NONE) {
			X[0] = ( Xt/rix + x0i)/scale;
			X[1] = ( Yt/riy + y0i)/scale;
			X[2] = ( Zt/riz + z0i)/scale;
		} else if (type==RIGID) {
			X[0] = ( (rot[0][0]*(Xt-trans[3])+rot[1][0]*(Yt-trans[4])+rot[2][0]*(Zt-trans[5]) )/rix + x0i)/scale;
			X[1] = ( (rot[0][1]*(Xt-trans[3])+rot[1][1]*(Yt-trans[4])+rot[2][1]*(Zt-trans[5]) )/riy + y0i)/scale;
			X[2] = ( (rot[0][2]*(Xt-trans[3])+rot[1][2]*(Yt-trans[4])+rot[2][2]*(Zt-trans[5]) )/riz + z0i)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			
			X[0] = ( (rot[0][0]*(Xt/factor-trans[3])+rot[1][0]*(Yt/factor-trans[4])+rot[2][0]*(Zt/factor-trans[5]) )/rix + x0i)/scale;
			X[1] = ( (rot[0][1]*(Xt/factor-trans[3])+rot[1][1]*(Yt/factor-trans[4])+rot[2][1]*(Zt/factor-trans[5]) )/riy + y0i)/scale;
			X[2] = ( (rot[0][2]*(Xt/factor-trans[3])+rot[1][2]*(Yt/factor-trans[4])+rot[2][2]*(Zt/factor-trans[5]) )/riz + z0i)/scale;
		} else if (type==SCALED_RIGID) {
			float fx = 1.0f+trans[6]/Iscale;
			float fy = 1.0f+trans[7]/Iscale;
			float fz = 1.0f+trans[8]/Iscale;
			
			X[0] = ( (rot[0][0]*(Xt/fx-trans[3])+rot[1][0]*(Yt/fy-trans[4])+rot[2][0]*(Zt/fz-trans[5]) )/rix + x0i)/scale;
			X[1] = ( (rot[0][1]*(Xt/fx-trans[3])+rot[1][1]*(Yt/fy-trans[4])+rot[2][1]*(Zt/fz-trans[5]) )/riy + y0i)/scale;
			X[2] = ( (rot[0][2]*(Xt/fx-trans[3])+rot[1][2]*(Yt/fy-trans[4])+rot[2][2]*(Zt/fz-trans[5]) )/riz + z0i)/scale;
		} else if (type==FULLY_AFFINE) {
			float[][] mat = new float[3][4];
			mat[0][0] = (1+trans[0]);
			mat[0][1] = trans[1];
			mat[0][2] = trans[2];
			mat[0][3] = trans[9];
			mat[1][0] = trans[3];
			mat[1][1] = (1+trans[4]);
			mat[1][2] = trans[5];
			mat[1][3] = trans[10];
			mat[2][0] = trans[6];
			mat[2][1] = trans[7];
			mat[2][2] = (1+trans[8]);
			mat[2][3] = trans[11];
			
			mat = Matrix3D.invertAffineTransform(mat);
			
			X[0] = ( (mat[0][0]*Xt + mat[0][1]*Yt + mat[0][2]*Zt + mat[0][3])/rix + x0i)/scale;
			X[1] = ( (mat[1][0]*Xt + mat[1][1]*Yt + mat[1][2]*Zt + mat[1][3])/riy + y0i)/scale;
			X[2] = ( (mat[2][0]*Xt + mat[2][1]*Yt + mat[2][2]*Zt + mat[2][3])/riz + z0i)/scale;
		}
		return;
	}

	public final void imageToTemplateDerivatives(float[][] dX, int x,int y,int z,float[] trans, float[][] rot, float[][] dRa, float[][] dRb, float[][] dRc, float scale) {
		for (int i=0;i<3;i++) 
			for (int t=0;t<Nt;t++)
				dX[i][t] = 0.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
			return;
		} else if (type==RIGID) {
			dX[0][0] = ( (dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( (dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( (dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = 1.0f/rtx/scale;
			
			dX[1][0] = ( (dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( (dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( (dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = 1.0f/rty/scale;
			
			dX[2][0] = ( (dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( (dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( (dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = 1.0f/rtz/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dX[0][0] = ( factor*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( factor*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( factor*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = factor/rtx/scale;
			dX[0][6] = ( dfactor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][0] = ( factor*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( factor*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( factor*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = factor/rty/scale;
			dX[1][6] = ( dfactor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][0] = ( factor*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( factor*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( factor*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = factor/rtz/scale;
			dX[2][6] = ( dfactor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==SCALED_RIGID) {
			float factorx = 1.0f + trans[6]/Iscale;
			float factory = 1.0f + trans[7]/Iscale;
			float factorz = 1.0f + trans[8]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dX[0][0] = ( factorx*(dRa[0][0]*Xi+dRa[0][1]*Yi+dRa[0][2]*Zi)/rtx )/scale;
			dX[0][1] = ( factorx*(dRb[0][0]*Xi+dRb[0][1]*Yi+dRb[0][2]*Zi)/rtx )/scale;
			dX[0][2] = ( factorx*(dRc[0][0]*Xi+dRc[0][1]*Yi+dRc[0][2]*Zi)/rtx )/scale;
			dX[0][3] = factorx/rtx/scale;
			dX[0][6] = ( dfactor*(rot[0][0]*Xi+rot[0][1]*Yi+rot[0][2]*Zi+trans[3])/rtx )/scale;
			
			dX[1][0] = ( factory*(dRa[1][0]*Xi+dRa[1][1]*Yi+dRa[1][2]*Zi)/rty )/scale;
			dX[1][1] = ( factory*(dRb[1][0]*Xi+dRb[1][1]*Yi+dRb[1][2]*Zi)/rty )/scale;
			dX[1][2] = ( factory*(dRc[1][0]*Xi+dRc[1][1]*Yi+dRc[1][2]*Zi)/rty )/scale;
			dX[1][4] = factory/rty/scale;
			dX[1][7] = ( dfactor*(rot[1][0]*Xi+rot[1][1]*Yi+rot[1][2]*Zi+trans[4])/rty )/scale;
			
			dX[2][0] = ( factorz*(dRa[2][0]*Xi+dRa[2][1]*Yi+dRa[2][2]*Zi)/rtz )/scale;
			dX[2][1] = ( factorz*(dRb[2][0]*Xi+dRb[2][1]*Yi+dRb[2][2]*Zi)/rtz )/scale;
			dX[2][2] = ( factorz*(dRc[2][0]*Xi+dRc[2][1]*Yi+dRc[2][2]*Zi)/rtz )/scale;
			dX[2][5] = factorz/rtz/scale;
			dX[2][8] = ( dfactor*(rot[2][0]*Xi+rot[2][1]*Yi+rot[2][2]*Zi+trans[5])/rtz )/scale;
		} else if (type==FULLY_AFFINE) {
			dX[0][0] = ( Xi/rtx )/scale;
			dX[0][1] = ( Yi/rtx )/scale;
			dX[0][2] = ( Zi/rtx )/scale;
			dX[0][9] = 1/rtx/scale;
			
			dX[1][3] = ( Xi/rty )/scale;
			dX[1][4] = ( Yi/rty )/scale;
			dX[1][5] = ( Zi/rty )/scale;
			dX[1][10] = 1/rty/scale;
			
			dX[2][6] = ( Xi/rtz )/scale;
			dX[2][7] = ( Yi/rtz )/scale;
			dX[2][8] = ( Zi/rtz )/scale;
			dX[2][11] = 1/rtz/scale;
		}
		return;
	}

	/** 
	 *	computes the transformed coordinates from image to template space
	 *	with scaling (2: half res, 4: quarter res, etc)
	 */
	public final float normScalingFactor(int x,int y,int z, float[] trans, float[][] rot, float scale) {
		float scaling = 1.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
			scaling = 1.0f;
		} else if (type==RIGID) {
			scaling = 1.0f;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			scaling = factor*factor*factor;
		} else if (type==SCALED_RIGID) {
			scaling = (1.0f+trans[6]/Iscale)*(1.0f+trans[7]/Iscale)*(1.0f+trans[8]/Iscale);
		} else if (type==FULLY_AFFINE) {
			scaling = (1+trans[0])*(1+trans[4])*(1+trans[8])
						+ trans[1]*trans[5]*trans[6]
						+ trans[2]*trans[3]*trans[7]
						- (1+trans[0])*trans[5]*trans[7]
						- trans[2]*(1+trans[4])*trans[6]
						- trans[1]*trans[3]*(1+trans[8]);
		}
		return scaling;
	}

	public final float[] normScalingFactorDerivatives(int x,int y,int z,float[] trans, float[][] rot, float[][] dRa, float[][] dRb, float[][] dRc, float scale) {
		float[] dF = new float[Nt];
		
		for (int t=0;t<Nt;t++)
			dF[t] = 0.0f;
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		
		if (type==NONE) {
		} else if (type==RIGID) {
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f + trans[6]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dF[6] = 3.0f*dfactor*factor*factor;
		} else if (type==SCALED_RIGID) {
			float factorx = 1.0f + trans[6]/Iscale;
			float factory = 1.0f + trans[7]/Iscale;
			float factorz = 1.0f + trans[8]/Iscale;
			float dfactor = 1.0f/Iscale;
			
			dF[6] = dfactor*factory*factorz;
			dF[7] = dfactor*factorz*factorx;
			dF[8] = dfactor*factorx*factory;
		} else if (type==FULLY_AFFINE) {
			dF[0] = (1+trans[4])*(1+trans[8]) - trans[5]*trans[7];
			dF[1] = trans[5]*trans[6] - trans[3]*(1+trans[8]);
			dF[2] = trans[3]*trans[7] - (1+trans[4])*trans[6];
			dF[3] = trans[2]*trans[7] - trans[1]*(1+trans[8]);
			dF[4] = (1+trans[0])*(1+trans[8]) - trans[2]*trans[6];
			dF[5] = trans[1]*trans[6] - (1+trans[0])*trans[7];
			dF[6] = trans[1]*trans[5] - trans[2]*(1+trans[4]);
			dF[7] = trans[2]*trans[3] - (1+trans[0])*trans[5];
			dF[8] = (1+trans[0])*(1+trans[4]) - trans[1]*trans[3];
		}
		return dF;
	}

	public final void precomputeImageToTemplateMatrix(float[][] transformMatrix, float[] trans, float[][] rot, float scale) {
		
		if (type==NONE) {
			transformMatrix[0][0] = rix/rtx;
			transformMatrix[0][1] = 0;
			transformMatrix[0][2] = 0;
			transformMatrix[1][0] = 0;
			transformMatrix[1][1] = riy/rty;
			transformMatrix[1][2] = 0;
			transformMatrix[2][0] = 0;
			transformMatrix[2][1] = 0;
			transformMatrix[2][2] = riz/rtz;
			
			transformMatrix[0][3] = ((-x0i)*rix/rtx + x0t)/scale;
			transformMatrix[1][3] = ((-y0i)*riy/rty + y0t)/scale;
			transformMatrix[2][3] = ((-z0i)*riz/rtz + z0t)/scale;		
		} else if (type==RIGID) {
			transformMatrix[0][0] = rot[0][0]*rix/rtx;
			transformMatrix[0][1] = rot[0][1]*riy/rtx;
			transformMatrix[0][2] = rot[0][2]*riz/rtx;
			transformMatrix[1][0] = rot[1][0]*rix/rty;
			transformMatrix[1][1] = rot[1][1]*riy/rty;
			transformMatrix[1][2] = rot[1][2]*riz/rty;
			transformMatrix[2][0] = rot[2][0]*rix/rtz;
			transformMatrix[2][1] = rot[2][1]*riy/rtz;
			transformMatrix[2][2] = rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = ((rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = ((rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = ((rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f+trans[6]/Iscale;
			
			transformMatrix[0][0] = factor*rot[0][0]*rix/rtx;
			transformMatrix[0][1] = factor*rot[0][1]*riy/rtx;
			transformMatrix[0][2] = factor*rot[0][2]*riz/rtx;
			transformMatrix[1][0] = factor*rot[1][0]*rix/rty;
			transformMatrix[1][1] = factor*rot[1][1]*riy/rty;
			transformMatrix[1][2] = factor*rot[1][2]*riz/rty;
			transformMatrix[2][0] = factor*rot[2][0]*rix/rtz;
			transformMatrix[2][1] = factor*rot[2][1]*riy/rtz;
			transformMatrix[2][2] = factor*rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = (factor*(rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = (factor*(rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = (factor*(rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		}  else if (type==SCALED_RIGID) {
			float factorx = 1.0f+trans[6]/Iscale;
			float factory = 1.0f+trans[7]/Iscale;
			float factorz = 1.0f+trans[8]/Iscale;
			
			transformMatrix[0][0] = factorx*rot[0][0]*rix/rtx;
			transformMatrix[0][1] = factorx*rot[0][1]*riy/rtx;
			transformMatrix[0][2] = factorx*rot[0][2]*riz/rtx;
			transformMatrix[1][0] = factory*rot[1][0]*rix/rty;
			transformMatrix[1][1] = factory*rot[1][1]*riy/rty;
			transformMatrix[1][2] = factory*rot[1][2]*riz/rty;
			transformMatrix[2][0] = factorz*rot[2][0]*rix/rtz;
			transformMatrix[2][1] = factorz*rot[2][1]*riy/rtz;
			transformMatrix[2][2] = factorz*rot[2][2]*riz/rtz;
			
			transformMatrix[0][3] = (factorx*(rot[0][0]*(-x0i)*rix+rot[0][1]*(-y0i)*riy+rot[0][2]*(-z0i)*riz+trans[3])/rtx + x0t)/scale;
			transformMatrix[1][3] = (factory*(rot[1][0]*(-x0i)*rix+rot[1][1]*(-y0i)*riy+rot[1][2]*(-z0i)*riz+trans[4])/rty + y0t)/scale;
			transformMatrix[2][3] = (factorz*(rot[2][0]*(-x0i)*rix+rot[2][1]*(-y0i)*riy+rot[2][2]*(-z0i)*riz+trans[5])/rtz + z0t)/scale;
		} else if (type==FULLY_AFFINE) {
			transformMatrix[0][0] = (1+trans[0])*rix/rtx;
			transformMatrix[0][1] = (0+trans[1])*riy/rtx;
			transformMatrix[0][2] = (0+trans[2])*riz/rtx;
			transformMatrix[1][0] = (0+trans[3])*rix/rty;
			transformMatrix[1][1] = (1+trans[4])*riy/rty;
			transformMatrix[1][2] = (0+trans[5])*riz/rty;
			transformMatrix[2][0] = (0+trans[6])*rix/rtz;
			transformMatrix[2][1] = (0+trans[7])*riy/rtz;
			transformMatrix[2][2] = (1+trans[8])*riz/rtz;
			
			transformMatrix[0][3] = ((((1+trans[0])*(-x0i)*rix+   trans[1] *(-y0i)*riy+   trans[2] *(-z0i)*riz)+trans[ 9])/rtx + x0t)/scale;
			transformMatrix[1][3] = (((   trans[3] *(-x0i)*rix+(1+trans[4])*(-y0i)*riy+   trans[5] *(-z0i)*riz)+trans[10])/rty + y0t)/scale;
			transformMatrix[2][3] = (((   trans[6] *(-x0i)*rix+   trans[7] *(-y0i)*riy+(1+trans[8])*(-z0i)*riz)+trans[11])/rtz + z0t)/scale;
		} 	
		
		return;
	}
	
	public final void precomputeTemplateToImageMatrix(float[][] transformMatrix, float[] trans, float[][] rot, float scale) {
		
		if (type==NONE) {
			transformMatrix[0][0] = rtx/rix;
			transformMatrix[0][1] = 0;
			transformMatrix[0][2] = 0;
			transformMatrix[1][0] = 0;
			transformMatrix[1][1] = rty/riy;
			transformMatrix[1][2] = 0;
			transformMatrix[2][0] = 0;
			transformMatrix[2][1] = 0;
			transformMatrix[2][2] = rtz/riz;
			
			transformMatrix[0][3] = ((-x0t)*rtx/rix + x0i)/scale;
			transformMatrix[1][3] = ((-y0t)*rty/riy + y0i)/scale;
			transformMatrix[2][3] = ((-z0t)*rtz/riz + z0i)/scale;		
		} else if (type==RIGID) {
			transformMatrix[0][0] = rot[0][0]*rtx/rix;
			transformMatrix[0][1] = rot[1][0]*rty/rix;
			transformMatrix[0][2] = rot[2][0]*rtz/rix;
			transformMatrix[1][0] = rot[0][1]*rtx/riy;
			transformMatrix[1][1] = rot[1][1]*rty/riy;
			transformMatrix[1][2] = rot[2][1]*rtz/riy;
			transformMatrix[2][0] = rot[0][2]*rtx/riz;
			transformMatrix[2][1] = rot[1][2]*rty/riz;
			transformMatrix[2][2] = rot[2][2]*rtz/riz;
			
			transformMatrix[0][3] = ((rot[0][0]*((-x0t)*rtx-trans[3])+rot[1][0]*((-y0t)*rty-trans[4])+rot[2][0]*((-z0t)*rtz-trans[5]))/rix + x0i)/scale;
			transformMatrix[1][3] = ((rot[0][1]*((-x0t)*rtx-trans[3])+rot[1][1]*((-y0t)*rty-trans[4])+rot[2][1]*((-z0t)*rtz-trans[5]))/riy + y0i)/scale;
			transformMatrix[2][3] = ((rot[0][2]*((-x0t)*rtx-trans[3])+rot[1][2]*((-y0t)*rty-trans[4])+rot[2][2]*((-z0t)*rtz-trans[5]))/riz + z0i)/scale;
		} else if (type==SINGLE_SCALE) {
			float factor = 1.0f+trans[6]/Iscale;
			
			transformMatrix[0][0] = rot[0][0]*rtx/rix/factor;
			transformMatrix[0][1] = rot[1][0]*rty/rix/factor;
			transformMatrix[0][2] = rot[2][0]*rtz/rix/factor;
			transformMatrix[1][0] = rot[0][1]*rtx/riy/factor;
			transformMatrix[1][1] = rot[1][1]*rty/riy/factor;
			transformMatrix[1][2] = rot[2][1]*rtz/riy/factor;
			transformMatrix[2][0] = rot[0][2]*rtx/riz/factor;
			transformMatrix[2][1] = rot[1][2]*rty/riz/factor;
			transformMatrix[2][2] = rot[2][2]*rtz/riz/factor;
			
			transformMatrix[0][3] = ((rot[0][0]*((-x0t/factor)*rtx-trans[3])+rot[1][0]*((-y0t/factor)*rty-trans[4])+rot[2][0]*((-z0t/factor)*rtz-trans[5]))/rix + x0i)/scale;
			transformMatrix[1][3] = ((rot[0][1]*((-x0t/factor)*rtx-trans[3])+rot[1][1]*((-y0t/factor)*rty-trans[4])+rot[2][1]*((-z0t/factor)*rtz-trans[5]))/riy + y0i)/scale;
			transformMatrix[2][3] = ((rot[0][2]*((-x0t/factor)*rtx-trans[3])+rot[1][2]*((-y0t/factor)*rty-trans[4])+rot[2][2]*((-z0t/factor)*rtz-trans[5]))/riz + z0i)/scale;
		}  else if (type==SCALED_RIGID) {
			float fx = 1.0f+trans[6]/Iscale;
			float fy = 1.0f+trans[7]/Iscale;
			float fz = 1.0f+trans[8]/Iscale;
			
			transformMatrix[0][0] = rot[0][0]*rtx/rix/fx;
			transformMatrix[0][1] = rot[1][0]*rty/rix/fy;
			transformMatrix[0][2] = rot[2][0]*rtz/rix/fz;
			transformMatrix[1][0] = rot[0][1]*rtx/riy/fx;
			transformMatrix[1][1] = rot[1][1]*rty/riy/fy;
			transformMatrix[1][2] = rot[2][1]*rtz/riy/fz;
			transformMatrix[2][0] = rot[0][2]*rtx/riz/fx;
			transformMatrix[2][1] = rot[1][2]*rty/riz/fy;
			transformMatrix[2][2] = rot[2][2]*rtz/riz/fz;
			
			transformMatrix[0][3] = ((rot[0][0]*((-x0t/fx)*rtx-trans[3])+rot[1][0]*((-y0t/fy)*rty-trans[4])+rot[2][0]*((-z0t/fz)*rtz-trans[5]))/rix + x0i)/scale;
			transformMatrix[1][3] = ((rot[0][1]*((-x0t/fx)*rtx-trans[3])+rot[1][1]*((-y0t/fy)*rty-trans[4])+rot[2][1]*((-z0t/fz)*rtz-trans[5]))/riy + y0i)/scale;
			transformMatrix[2][3] = ((rot[0][2]*((-x0t/fx)*rtx-trans[3])+rot[1][2]*((-y0t/fy)*rty-trans[4])+rot[2][2]*((-z0t/fz)*rtz-trans[5]))/riz + z0i)/scale;
		} else if (type==FULLY_AFFINE) {
			float[][] mat = new float[3][4];
			mat[0][0] = (1+trans[0]);
			mat[0][1] = trans[1];
			mat[0][2] = trans[2];
			mat[0][3] = trans[9];
			mat[1][0] = trans[3];
			mat[1][1] = (1+trans[4]);
			mat[1][2] = trans[5];
			mat[1][3] = trans[10];
			mat[2][0] = trans[6];
			mat[2][1] = trans[7];
			mat[2][2] = (1+trans[8]);
			mat[2][3] = trans[11];

			mat = Matrix3D.invertAffineTransform(mat);
			
			transformMatrix[0][0] = mat[0][0]*rtx/rix;
			transformMatrix[0][1] = mat[0][1]*rty/rix;
			transformMatrix[0][2] = mat[0][2]*rtz/rix;
			transformMatrix[1][0] = mat[1][0]*rtx/riy;
			transformMatrix[1][1] = mat[1][1]*rty/riy;
			transformMatrix[1][2] = mat[1][2]*rtz/riy;
			transformMatrix[2][0] = mat[2][0]*rtx/riz;
			transformMatrix[2][1] = mat[2][1]*rty/riz;
			transformMatrix[2][2] = mat[2][2]*rtz/riz;
			
			transformMatrix[0][3] = (((mat[0][0]*(-x0t)*rtx + mat[0][1]*(-y0t)*rty + mat[0][2]*(-z0t)*rtz) + mat[0][3])/rix + x0i)/scale;
			transformMatrix[1][3] = (((mat[1][0]*(-x0t)*rtx + mat[1][1]*(-y0t)*rty + mat[1][2]*(-z0t)*rtz) + mat[1][3])/riy + y0i)/scale;
			transformMatrix[2][3] = (((mat[2][0]*(-x0t)*rtx + mat[2][1]*(-y0t)*rty + mat[2][2]*(-z0t)*rtz) + mat[2][3])/riz + z0i)/scale;
		} 	
		
		return;
	}
	
	public final RotationMatrix computeRotationMatrix(float[] trans) {
		RotationMatrix R = new RotationMatrix();
		R.setParameters(trans[0],trans[1],trans[2]);
		return R;
	}
	
	public static final float[][] computeRotation(float[] trans) {
		RotationMatrix R = new RotationMatrix();
		R.setParameters(trans[0],trans[1],trans[2]);
		return R.getMatrix();
	}
	
	public static final float[] changeTransformType(float[] trans, String oldType_, String newType_) {
		int oldType, newType;
		int No, Nn;
		float[] newtrans=null;
		
			 if (oldType_.equals("rigid")) { oldType = RIGID; No = 6; } 
		else if (oldType_.equals("single_scale")) { oldType = SINGLE_SCALE; No = 7; }
		else if (oldType_.equals("scaled_rigid")) { oldType = SCALED_RIGID; No = 9; }
		else if (oldType_.equals("fully_affine")) { oldType = FULLY_AFFINE; No = 12; }
		else { oldType = NONE; No = 0; }

			 if (newType_.equals("rigid")) { newType = RIGID; Nn = 6; } 
		else if (newType_.equals("single_scale")) { newType = SINGLE_SCALE; Nn = 7; }
		else if (newType_.equals("scaled_rigid")) { newType = SCALED_RIGID; Nn = 9; }
		else if (newType_.equals("fully_affine")) { newType = FULLY_AFFINE; Nn = 12; }
		else { newType = NONE; Nn = 0; }

		if (oldType==newType) return trans;
		else if (oldType==NONE) {
			newtrans = new float[Nn];
			for (int n=0;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==RIGID && (newType==SINGLE_SCALE || newType==SCALED_RIGID) ) {
			newtrans = new float[Nn];
			for (int n=0;n<No;n++) newtrans[n] = trans[n];
			for (int n=No;n<Nn;n++) newtrans[n] = 0.0f;
		}else if (oldType==SINGLE_SCALE && newType==SCALED_RIGID ) {
			newtrans = new float[Nn];
			for (int n=0;n<No;n++) newtrans[n] = trans[n];
			newtrans[No] = trans[No-1];
			newtrans[No+1] = trans[No-1];
			for (int n=No+2;n<Nn;n++) newtrans[n] = 0.0f;
		} else if (oldType==RIGID && newType==FULLY_AFFINE ) {
			float[][] rot = computeRotation(trans);
			newtrans = new float[Nn];
			for (int n=0;n<3;n++) for (int m=0;m<3;m++) {
				if (n==m) newtrans[n+3*m] = rot[m][n]-1.0f;
				else newtrans[n+3*m] = rot[m][n];
			}
			for (int n=0;n<3;n++) newtrans[9+n] = trans[3+n];
			for (int n=12;n<Nn;n++) newtrans[n] = 0.0f;
		} else {
			BasicInfo.displayMessage("change from "+oldType_+" to "+newType_+" currently not supported\n");
			return trans;
		}
		BasicInfo.displayMessage("change from "+oldType_+" to "+newType_+" done ("+newtrans.length+")\n");
		return newtrans;
	}
	
	/** 
	 *	computes transformed directions from template to image space
	 */
	public final void templateToImageDirection(float[] dir, int x,int y,int z, float[] trans, float[][] rot, float scale) {
		float norm = (float)Math.sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
		
		if (norm<1e-30f) return;
		
		float[] V = new float[3];
		
		float Xi = (x*scale-x0i)*rix;
		float Yi = (y*scale-y0i)*riy;
		float Zi = (z*scale-z0i)*riz;
		if (type==NONE) {
			return;
		} else if (type==RIGID) {
			V[0] = (rot[0][0]*dir[0]+rot[1][0]*dir[1]+rot[2][0]*dir[2]);
			V[1] = (rot[0][1]*dir[0]+rot[1][1]*dir[1]+rot[2][1]*dir[2]);
			V[2] = (rot[0][2]*dir[0]+rot[1][2]*dir[1]+rot[2][2]*dir[2]);
		} else if (type==SINGLE_SCALE) {
			V[0] = (rot[0][0]*dir[0]+rot[1][0]*dir[1]+rot[2][0]*dir[2]);
			V[1] = (rot[0][1]*dir[0]+rot[1][1]*dir[1]+rot[2][1]*dir[2]);
			V[2] = (rot[0][2]*dir[0]+rot[1][2]*dir[1]+rot[2][2]*dir[2]);
		} else if (type==SCALED_RIGID) {
			float sx = 1.0f/(1.0f+trans[6]/Iscale);
			float sy = 1.0f/(1.0f+trans[7]/Iscale);
			float sz = 1.0f/(1.0f+trans[8]/Iscale);
			
			V[0] = (rot[0][0]*dir[0]/sx+rot[1][0]*dir[1]/sy+rot[2][0]*dir[2]/sz);
			V[1] = (rot[0][1]*dir[0]/sx+rot[1][1]*dir[1]/sy+rot[2][1]*dir[2]/sz);
			V[2] = (rot[0][2]*dir[0]/sx+rot[1][2]*dir[1]/sy+rot[2][2]*dir[2]/sz);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		} else if (type==FULLY_AFFINE) {
			float[][] M = new float[3][3];
			M[0][0] = (1+trans[0]);
			M[0][1] = trans[1];
			M[0][2] = trans[2];
			M[1][0] = trans[3];
			M[1][1] = (1+trans[4]);
			M[1][2] = trans[5];
			M[2][0] = trans[6];
			M[2][1] = trans[7];
			M[2][2] = (1+trans[8]);
			Matrix3D.invert(M);
			
			V[0] = (M[0][0]*dir[0]+M[0][1]*dir[1]+M[0][2]*dir[2]);
			V[1] = (M[1][0]*dir[0]+M[1][1]*dir[1]+M[1][2]*dir[2]);
			V[2] = (M[2][0]*dir[0]+M[2][1]*dir[1]+M[2][2]*dir[2]);
			
			float nV = (float)Math.sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			V[0] *= norm/nV;
			V[1] *= norm/nV;
			V[2] *= norm/nV;
		}
		dir[0] = V[0];
		dir[1] = V[1];
		dir[2] = V[2];
		
		return;
	}

	public final String displayTransform(float[] trans) {
		String info = "transform: (";
		for (int n=0;n<Nt-1;n++) info += trans[n]+", ";
		info += trans[Nt-1]+")\n";
		
		return info;
	}
	
}
