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
	private		float[] 	inv1 = null;
	private		float[] 	inv2 = null;
	private		float[]		b1map = null;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private		float		TRinversion = 5.50f;
	private		float		TRexcitation1 = 0.0062f;
	private		float		TRexcitation2 = 0.0062f;
	private		float		TI1 = 0.8f;
	private		float		TI2 = 2.7f;
	private		float		angle1 = 7.0f;
	private		float		angle2 = 5.0f;
	private 	int			Nexcitations = 160;
	private		float		inversionEfficiency = 0.96f;
	
	private     boolean     scalePhase = true;
	private		float		intensityScale = 4000.0f;
	private		float		t1Max = 4.0f;
	private		float		t1Min = 0.05f;
	private		int			t1Samples = 3950;
	private		int			uniSamples = 4000;
	
	private		boolean		useB1correction = true;
	private		float 		b1Scaling = 1.0f;
	private		float		b1min = 0.05f;
	private		float		b1max = 2.0f;
	private 	int			b1Samples = 1950;
	
	// output parameters
	private		float[] uni = null;
	private		float[] t1map = null;
	private		float[] r1map = null;
	private		float[] snr = null;
	
	private		double[][] T1lookup = null;
	private 	double[][] unilookup = null;
	
	// set inputs
	public final void setFirstInversionImage(float[] in) { inv1 = in; }
	public final void setSecondInversionImage(float[] in) { inv2 = in; }
	public final void setB1mapImage(float[] in) { b1map = in; }
	
	public final void setFirstInversionMagnitude(float[] in) {
	    if (inv1==null) {
	        System.out.print("build combined image\n");
	        inv1 = new float[2*nxyz];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) inv1[xyz] = in[xyz];
	}
	public final void setFirstInversionPhase(float[] in) {
	    if (inv1==null) {
	        System.out.print("build combined image\n");
	        inv1 = new float[2*nxyz];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) inv1[xyz+nxyz] = in[xyz];
	}
	public final void setSecondInversionMagnitude(float[] in) {
	    if (inv2==null) {
	        System.out.print("build combined image\n");
	        inv2 = new float[2*nxyz];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) inv2[xyz] = in[xyz];
	}
	public final void setSecondInversionPhase(float[] in) {
	    if (inv2==null) {
	        System.out.print("build combined image\n");
	        inv2 = new float[2*nxyz];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) inv2[xyz+nxyz] = in[xyz];
	}
	
	public final void setInversionRepetitionTime(float in) { TRinversion = in; }
	public final void setExcitationRepetitionTime(float in) { TRexcitation1 = in; TRexcitation2 = in; }
	public final void setFirstExcitationRepetitionTime(float in) { TRexcitation1 = in; }
	public final void setSecondExcitationRepetitionTime(float in) { TRexcitation2 = in; }
	public final void setNumberExcitations(int in) { Nexcitations = in; }
	public final void setFirstInversionTime(float in) { TI1 = in; }
	public final void setSecondInversionTime(float in) { TI2 = in; }
	public final void setFirstFlipAngle(float in) { angle1 = in; }
	public final void setSecondFlipAngle(float in) { angle2 = in; }
	public final void setInversionEfficiency(float in) { inversionEfficiency = in; }
	public final void setCorrectB1inhomogeneities(boolean in) { useB1correction = in; }
	public final void setB1mapScaling(float in) { b1Scaling = in;}
	public final void setScalePhaseImagess(boolean in) { scalePhase = in; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.MP2RAGE"; } 
	public final String getLabel() { return "MP2RAGE T1 Fitting"; }
	public final String getName() { return "MP2RAGE T1 Fitting"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimation of uniform T1-weighted and quantitative T1 parameter maps from the MP2RAGE sequence, following the model of Marques et al., 2010."; }
		
	public final String getVersion() { return "3.1f"; }

	// get outputs
	public float[] getUniformT1weightedImage() { return uni; }
	public float[] getQuantitativeT1mapImage() { return t1map; }
	public float[] getQuantitativeR1mapImage() { return r1map; }
	public float[] getRelativeSnrImage() { return snr; }
	
	// for debug
	public float[] generateT1LookupImage() { 
		/*
		float[] img = new float[lutSamples*lutSamples];
		for (int x=0;x<lutSamples;x++) {
			int y = Numerics.bounded(Numerics.round(T1lookup[x]/intensityScale*lutSamples),0,lutSamples);
			img[x+nx*y]++;
		}
		*/
		float[] img = new float[t1Samples*uniSamples*b1Samples];
		for (int x=0;x<t1Samples;x++) {
			for (int y=0;y<uniSamples;y++) {
				for (int z=0;z<b1Samples;z++) {
					img[x+t1Samples*y+t1Samples*uniSamples*z] = (float)T1lookup[z][x];
				}
			}
		}
		return img;
	}
	
	public float[] generateUniLookupImage() {  
		/*
		float[] img = new float[lutSamples*lutSamples];
		for (int x=0;x<lutSamples;x++) {
			int y = Numerics.bounded(Numerics.round(unilookup[x]/t1mapThreshold*lutSamples),0,lutSamples);
			img[x+nx*y]++;
		}*/
		float[] img = new float[uniSamples*t1Samples*b1Samples];
		for (int x=0;x<uniSamples;x++) {
			for (int y=0;y<t1Samples;y++) {
				for (int z=0;z<b1Samples;z++) {
					img[x+uniSamples*y+uniSamples*t1Samples*z] = (float)unilookup[z][x];
				}
			}
		}
		return img;
	}
	
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm
		System.out.print("Loading data\n");
		
		// we assume 4D images: dim 0 magnitude, dim 1 phase
		
		// estimate angular scale (range should be [-PI +PI]
		double phscale1 = 1.0;
		double phscale2 = 1.0;
		if (scalePhase) {
            float phsmin1 = inv1[nxyz];
            float phsmax1 = inv1[nxyz];
            float phsmin2 = inv1[nxyz];
            float phsmax2 = inv1[nxyz];
            for (int xyz=0;xyz<nxyz;xyz++) {
                if (inv1[xyz+nxyz]<phsmin1) phsmin1 = inv1[xyz+nxyz];
                if (inv1[xyz+nxyz]>phsmax1) phsmax1 = inv1[xyz+nxyz];
                if (inv2[xyz+nxyz]<phsmin2) phsmin2 = inv2[xyz+nxyz];
                if (inv2[xyz+nxyz]>phsmax2) phsmax2 = inv2[xyz+nxyz];
            }
            phscale1 = (phsmax1-phsmin1)/(2.0*FastMath.PI);
            phscale2 = (phsmax2-phsmin2)/(2.0*FastMath.PI);
            System.out.print("Phase scaling: ["+phsmin1+", "+phsmax1+"] -> [-PI, PI]\n");
            System.out.print("Phase scaling: ["+phsmin2+", "+phsmax2+"] -> [-PI, PI]\n");
		}
		
		uni = new float[nxyz];
		snr = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			double prod = inv1[xyz]*inv2[xyz]*FastMath.cos(inv1[xyz+nxyz]/phscale1-inv2[xyz+nxyz]/phscale2);
			double norm = inv1[xyz]*inv1[xyz] + inv2[xyz]*inv2[xyz];
			
			double diff2 = (inv1[xyz]*inv1[xyz] - inv2[xyz]*inv2[xyz])*(inv1[xyz]*inv1[xyz] - inv2[xyz]*inv2[xyz]);
			double norm3 = norm*norm*norm;
			
			if (intensityScale*norm>prod) uni[xyz] = (float)(prod/norm);
			else uni[xyz] = 0.0f;
			
			if (intensityScale*intensityScale*norm3>diff2) snr[xyz] = (float)FastMath.sqrt(diff2/norm3);
			else snr[xyz] = 1.0f;
		}
		
		if (!useB1correction) b1Samples = 1;	
		
		// for the T1 map, more complicated: first estimate the inv1, inv2 you would get for a range of t1 values, then find the fit?
		T1lookup = new double[b1Samples][t1Samples];
		
		double TA = TI1 - Nexcitations/2.0f*TRexcitation1;
		double TB = TI2 - TI1 - Nexcitations/2.0f*TRexcitation1 - Nexcitations/2.0f*TRexcitation2;
		double TC = TRinversion - TI2 - Nexcitations/2.0f*TRexcitation2;
		
		System.out.print("timings: TA = "+TA+", TB = "+TB+", TC = "+TC+"\n"); 
		
		double a1rad = angle1/180.0*FastMath.PI;
		double a2rad = angle2/180.0*FastMath.PI;
		
		for (int b=0;b<b1Samples;b++) {
		    System.out.print(".");
		
			double b1 = b1min + b*(b1max-b1min)/(b1Samples-1.0);
			if (!useB1correction) b1 = 1.0;
			double cosA1cosA2n = FastMath.pow( FastMath.cos(b1*a1rad)*FastMath.cos(b1*a2rad), Nexcitations );
		
			T1lookup[b][0] = 0.5;
			for (int t=0;t<t1Samples;t++) {
				double qt1 = t1Min + t*(t1Max-t1Min)/(t1Samples-1.0);
			
                double E1 = FastMath.exp(-TRexcitation1/qt1);
                double E2 = FastMath.exp(-TRexcitation2/qt1);
            
                double EI = FastMath.exp(-TRinversion/qt1);
            
                double EA = FastMath.exp(-TA/qt1);
                double EB = FastMath.exp(-TB/qt1);
                double EC = FastMath.exp(-TC/qt1);
            
                double E1cosA1 = FastMath.cos(b1*a1rad)*E1;
                double E2cosA2 = FastMath.cos(b1*a2rad)*E2;
            
                double E1cosA1n = FastMath.pow( (FastMath.cos(b1*a1rad)*E1), Nexcitations);
                double E2cosA2n = FastMath.pow( (FastMath.cos(b1*a2rad)*E2), Nexcitations);
            
                double E1cosA1n2 = FastMath.pow( (FastMath.cos(b1*a1rad)*E1), Nexcitations/2.0f);
                double E2cosA2n2 = FastMath.pow( (FastMath.cos(b1*a2rad)*E2), Nexcitations/2.0f);
            
                double mza = (1.0f-EA)*E1cosA1n + (1.0f-E1)*(1.0f-E1cosA1n)/(1.0f-E1cosA1);
                double mzb = ( mza*EB + (1.0f-EB) )*E2cosA2n + (1.0f-E2)*(1.0f-E2cosA2n)/(1.0f-E2cosA2);
                double mzc = ( mzb*EC + (1.0f-EC) )/(1.0f + inversionEfficiency*cosA1cosA2n*EI);
                double mzss = ( ( ( ( (1.0f-EA)*E1cosA1n + (1.0f-E1)*(1.0f-E1cosA1n)/(1.0f-E1cosA1) )*EB 
                                + (1.0f-EB) )*E2cosA2n + (1.0f-E2)*(1.0f-E2cosA2n)/(1.0f-E2cosA2) )*EC 
                                    + (1.0f-EC) )/(1.0f + inversionEfficiency*cosA1cosA2n*EI);
            
                //double gre1 = FastMath.sin(b1*a1rad)*( (-inversionEfficiency*mzss*EA + 1.0f-EA)*E1cosA1n21 + (1.0f-E1)*(1.0f-E1cosA1n21)/(1.0f-E1cosA1) );
                //double gre2 = FastMath.sin(b1*a2rad)*( (mzss-(1.0f-EC))/(EC*E1cosA2n2) - (1.0f-E1)*(1.0f/E1cosA2n2-1.0f)/(1.0f-E1cosA2) );
            
                double factora = ( (-inversionEfficiency*mzss*EA + 1.0-EA)*E1cosA1n2 + (1.0-E1)*(1.0-E1cosA1n2)/(1.0-E1cosA1) );
                double factorb = factora*E1cosA1n2 + (1.0-E1)*(1.0-E1cosA1n2)/(1.0-E1cosA1n2);
                double factorc = (factorb*EB + 1.0-EB)*E2cosA2n2 + (1.0-E2)*(1.0-E2cosA2n2)/(1.0-E2cosA2);
            
                double gre1 = FastMath.sin(b1*a1rad)*factora;
                double gre2 = FastMath.sin(b1*a2rad)*factorc;
            
                //T1lookup[t] = intensityScale/2.0f + intensityScale*gre1*gre2/(gre1*gre1+gre2*gre2);
                T1lookup[b][t] = gre1*gre2/(gre1*gre1+gre2*gre2);
			}
			/*
			if (Numerics.abs(gre2)>1 || Numerics.abs(gre1)>1) {
				System.out.println("gre1 = "+gre1+", gre2 = "+gre2+", mzss = "+mzss+", mzc = "+mzc+", E1 = "+E1);
				T1lookup[t] = 0.5f;
			} else {
				//T1lookup[t] = intensityScale/2.0f + intensityScale*gre1*gre2/(gre1*gre1+gre2*gre2);
				T1lookup[t] = gre1*gre2/(gre1*gre1+gre2*gre2);
			}*/
		}
		// invert the LUT for speed
		System.out.print("\n invert lookup table\n"); 
		unilookup = new double[b1Samples][uniSamples];
		for (int b=0;b<b1Samples;b++) {
			for (int u=0;u<uniSamples;u++) {
				//double val = intensityScale*t/(double)lutSamples;
				double val = -0.5 + 1.0*u/(uniSamples-1.0);
			
				int t=0;
				while (t<t1Samples && T1lookup[b][t]>val) t++;
			
				if (t>0 && t<t1Samples) {
				    // make a linear interpolation between previous and current value
					double ratio = Numerics.bounded( (T1lookup[b][t-1]-val)/(T1lookup[b][t-1]-T1lookup[b][t]), 0.0, 1.0);
					unilookup[b][u] = t1Min + ( ratio*t+(1.0-ratio)*(t-1) )*(t1Max-t1Min)/(t1Samples-1.0);
				} else if (t>=t1Samples) {
					unilookup[b][u] = t1Max;
				} else {
					unilookup[b][u] = t1Min;	
				}
			}
		}
		
		System.out.print(" compute final maps\n"); 
		t1map = new float[nxyz];
		r1map = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
		    double u = (uni[xyz]+0.5)*(uniSamples-1.0);
			int u0 = Numerics.bounded(Numerics.floor(u), 0, uniSamples-1);
			int u1 = Numerics.bounded(Numerics.ceil(u), 0, uniSamples-1);
			double uratio = Numerics.bounded( (u-u0), 0.0, 1.0);
                
			double b = 0;
			int b0 = 0, b1 = 0;
			double bratio = 0.0;
            if (useB1correction) {
			    b = (b1map[xyz]/b1Scaling-b1min)/(b1max-b1min)*(b1Samples-1.0);
			    b0 = Numerics.bounded(Numerics.floor(b), 0, b1Samples-1);
                b1 = Numerics.bounded(Numerics.ceil(b), 0, b1Samples-1);
                bratio = Numerics.bounded( (b-b0), 0.0, 1.0);
            }
			// linear interpolation here too
			t1map[xyz] = (float)( (1.0-uratio)*(1.0-bratio)*unilookup[b0][u0]
								 +(1.0-uratio)*bratio*unilookup[b1][u0]
								 +uratio*(1.0-bratio)*unilookup[b0][u1]
								 +uratio*bratio*unilookup[b1][u1]);
			
			// r1 map is just inverse
			r1map[xyz] = 1.0f/t1map[xyz];
			
			// re-correct the uni with the B1 information
			if (useB1correction) {
			    double t = (t1map[xyz]-t1Min)/(t1Max-t1Min)*(t1Samples-1.0);
			    int t0 = Numerics.bounded(Numerics.floor(t), 0, t1Samples-1);
			    int t1 = Numerics.bounded(Numerics.ceil(t), 0, t1Samples-1);
			    double tratio = Numerics.bounded( (t-t0), 0.0, 1.0);
			    
			    // assume B1 = 1 values this time
			    b = (1.0-b1min)/(b1max-b1min)*(b1Samples-1.0);
			    b0 = Numerics.bounded(Numerics.floor(b), 0, b1Samples-1);
                b1 = Numerics.bounded(Numerics.ceil(b), 0, b1Samples-1);
                bratio = Numerics.bounded( (b-b0), 0.0, 1.0);
                
			    uni[xyz] = (float)( (1.0-tratio)*(1.0-bratio)*T1lookup[b0][t0]
			    				   +(1.0-tratio)*bratio*T1lookup[b1][t0]
			    				   +tratio*(1.0-bratio)*T1lookup[b0][t1]
			    				   +tratio*bratio*T1lookup[b1][t1]);			    
			}
		}
		System.out.print(" scale outputs\n"); 
		
		// rescale the UNI to [0 4000]
		for (int xyz=0;xyz<nxyz;xyz++) {
			uni[xyz] = Numerics.bounded(intensityScale*uni[xyz]+intensityScale/2.0f,0.0f,intensityScale);
		}

		// rescale the T1 map to milliseconds
		for (int xyz=0;xyz<nxyz;xyz++) {
			t1map[xyz] = Numerics.bounded(1000.0f*t1map[xyz],1000.0f*t1Min,1000.0f*t1Max);
		}
			
		// rescale the R1 map to mHz
		for (int xyz=0;xyz<nxyz;xyz++) {
			r1map[xyz] = Numerics.bounded(1000.0f*r1map[xyz],1000.0f/t1Max,1000.0f/t1Min);
		}
			
		return;
	}

}
