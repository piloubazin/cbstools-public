package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensityMp2rageT1Fitting {

	// input parameters
	private		float[] 	inv1;
	private		float[] 	inv2;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private		float		TRinversion = 8.0f;
	private		float		TRexcitation = 0.0062f;
	private		float		TI1 = 0.8f;
	private		float		TI2 = 2.7f;
	private		float		angle1 = 7.0f;
	private		float		angle2 = 5.0f;
	private 	int			Nexcitations = 160;
	private		float		inversionEfficiency = 0.96f;
	
	private		float		angleScale = (float)(3142.0f/FastMath.PI);
	private		float		intensityScale = 4000.0f;
	private		float		t1mapThreshold = 4000.0f;
	private		int			lutSamples = 4000;
	
	// output parameters
	private		float[] uni = null;
	private		float[] t1map = null;
	private		float[] snr = null;
	
	// set inputs
	public final void setFirstInversionImage(float[] in) { inv1 = in; }
	public final void setSecondInversionImage(float[] in) { inv2 = in; }
	
	public final void setInversionRepetitionTime(float in) { TRinversion = in; }
	public final void setExcitationRepetitionTime(float in) { TRexcitation = in; }
	public final void setNumberExcitations(int in) { Nexcitations = in; }
	public final void setFirstInversionTime(float in) { TI1 = in; }
	public final void setSecondInversionTime(float in) { TI2 = in; }
	public final void setFirstFlipAngle(float in) { angle1 = in; }
	public final void setSecondFlipAngle(float in) { angle2 = in; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; } 
	public final String getLabel() { return "MP2RAGE T1 Fitting"; }
	public final String getName() { return "MP2RAGE T1 Fitting"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimation of uniform T1-weighted and quantitative T1 parameter maps from the MP2RAGE sequence, following the model of Marques et al., 2010."; }
		
	public final String getVersion() { return "3.1"; }

	// get outputs
	public float[] getUniformT1weightedImage() { return uni; }
	public float[] getQuantitativeT1mapImage() { return t1map; }
	public float[] getRelativeSnrImage() { return snr; }
	
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// we assume 4D images: dim 0 magnitude, dim 1 phase
		
		uni = new float[nxyz];
		snr = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			double prod = inv1[xyz]*inv2[xyz]*FastMath.cos((inv1[xyz+nxyz]-inv2[xyz+nxyz])/angleScale);
			double norm = inv1[xyz]*inv1[xyz] + inv2[xyz]*inv2[xyz];
			
			double diff2 = (inv1[xyz]*inv1[xyz] - inv2[xyz]*inv2[xyz])*(inv1[xyz]*inv1[xyz] - inv2[xyz]*inv2[xyz]);
			double norm3 = norm*norm*norm;
			
			if (intensityScale*norm>prod) uni[xyz] = (float)(prod/norm);
			else uni[xyz] = intensityScale;
			
			if (intensityScale*intensityScale*norm3>diff2) snr[xyz] = (float)FastMath.sqrt(diff2/norm3);
			else snr[xyz] = intensityScale;
		}
		
		// for the T1 map, more complicated: first estimate the inv1, inv2 you would get for a range of t1 values, then find the fit?
		double[] T1lookup = new double[lutSamples];
		
		double TA = TI1 - Nexcitations/2.0*TRexcitation;
		double TB = TI2 - 3.0*Nexcitations/2.0*TRexcitation - TA;
		double TC = TRinversion - 2.0*Nexcitations*TRexcitation - TA - TB;
		
		double a1rad = angle1/180.0*FastMath.PI;
		double a2rad = angle2/180.0*FastMath.PI;
		
		double cosA1cosA2n = FastMath.pow( FastMath.cos(a1rad)*FastMath.cos(a2rad), Nexcitations );
		
		for (int t=0;t<lutSamples;t++) {
			double qt1 = t*t1mapThreshold/lutSamples;
			
			double E1 = FastMath.exp(-TRexcitation/qt1);
			double EI = FastMath.exp(-TRinversion/qt1);
			double EA = FastMath.exp(-TA/qt1);
			double EB = FastMath.exp(-TB/qt1);
			double EC = FastMath.exp(-TC/qt1);
			double E1cosA1 = FastMath.cos(a1rad)*E1;
			double E1cosA2 = FastMath.cos(a2rad)*E1;
			double E1cosA1n = FastMath.pow( (FastMath.cos(a1rad)*E1), Nexcitations);
			double E1cosA2n = FastMath.pow( (FastMath.cos(a2rad)*E1), Nexcitations);
			double E1cosA1n21 = FastMath.pow( (FastMath.cos(a1rad)*E1), Nexcitations/2.0-1.0);
			double E1cosA2n2 = FastMath.pow( (FastMath.cos(a2rad)*E1), Nexcitations/2.0);
			double mza = (1.0-EA)*E1cosA1n + (1.0-E1)*(1.0-E1cosA1n)/(1.0-E1cosA1);
			double mzb = ( mza*EB + (1.0-EB) )*E1cosA2n + (1.0-E1)*(1.0-E1cosA2n)/(1.0-E1cosA2);
			double mzc = ( mzb*EC + (1.0-EC) )/(1.0 + inversionEfficiency*cosA1cosA2n*EI);
			double mzss = ( ( ( ( (1.0-EA)*E1cosA1n + (1.0-E1)*(1.0-E1cosA1n)/(1.0-E1cosA1) )*EB + (1.0-EB) )*E1cosA2n + (1.0-E1)*(1.0-E1cosA2n)/(1.0-E1cosA2) )*EC + (1.0-EC) )/(1.0 + inversionEfficiency*cosA1cosA2n*EI);
			
			double gre1 = FastMath.sin(a1rad)*( (-inversionEfficiency*mzss*EA + 1.0-EA)*E1cosA1n21 + (1.0-E1)*(1.0-E1cosA1n21)/(1.0-E1cosA1) );
			double gre2 = FastMath.sin(a2rad)*( (mzss-(1.0-EC))/(EC*E1cosA2n2) - (1.0-E1)*(1.0/E1cosA2n2-1.0)/(1.0-E1cosA2) );
		
			T1lookup[t] = gre1*gre2/(gre1*gre1+gre2*gre2);
		}
		// invert the LUT for speed
		double[] unilookup = new double[lutSamples];
		for (int t=0;t<lutSamples;t++) {
			double val = (t-lutSamples/2.0)/lutSamples;
			
			int n=0;
			while (n<lutSamples && T1lookup[n]<val) n++;
			
			if (n>0) {
				double ratio = (val-T1lookup[n-1])/(T1lookup[n]-T1lookup[n-1]);
				unilookup[t] = ( ratio*n+(1.0-ratio)*(n-1) )*t1mapThreshold/lutSamples;
			} else {
				unilookup[t] = 0.0;	
			}
			
		}
		
		t1map = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			int t = Numerics.round(uni[xyz]*lutSamples + lutSamples/2.0);
			if (t>=0 && t<lutSamples)
				t1map[xyz] = (float)unilookup[t];
			else if (t<0)
				t1map[xyz] = 0.0f;
			else if (t>=lutSamples)
				t1map[xyz] = t1mapThreshold;
		}
			
		return;
	}


}
