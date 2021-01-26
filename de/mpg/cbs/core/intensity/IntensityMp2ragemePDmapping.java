package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class IntensityMp2ragemePDmapping {

	// input parameters
	private		float[] 	inv1 = null;
	private		float[] 	inv2 = null;
	private		float[] 	uni = null;
	private		float[]		b1map = null;
	private		float[] 	t1map = null;
	private		float[] 	r2smap = null;
	private		float[] 	s0img = null;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private		float		TRinversion = 5.50f;
	private		float		TEcho1 = 0.0030f;
	private		float		TRexcitation1 = 0.0062f;
	private		float		TRexcitation2 = 0.0062f;
	private		float		TI1 = 0.8f;
	private		float		TI2 = 2.7f;
	private		float		angle1 = 7.0f;
	private		float		angle2 = 5.0f;
	private 	int			Nexcitations = 160;
	private		float		inversionEfficiency = 0.96f;
	
	private     boolean     usePhase = true;
	private     boolean     scalePhase = true;
	
	private		boolean		useB1correction = true;
	private		float 		b1Scaling = 1.0f;
	
	// output parameters
	private		float[] pd1map = null;
	private		float[] pd2map = null;
	private		float[] pd0map = null;
	private		float[] pdmap = null;
	
	// set inputs
	public final void setFirstInversionImage(float[] in) { inv1 = in; }
	public final void setSecondInversionImage(float[] in) { inv2 = in; }

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
	
    public final void setB1mapImage(float[] in) { b1map = in; }
	public final void setT1mapImage(float[] in) { t1map = in; }
	public final void setR2smapImage(float[] in) { r2smap = in; }
	public final void setS0Image(float[] in) { s0img = in; }
	public final void setUniformImage(float[] in) { uni = in; }
		
	public final void setInversionRepetitionTime(float in) { TRinversion = in; }
	public final void setFirstEchoTime(float in) { TEcho1 = in; }
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
	public final String getCategory() { return "Intensity.MP2RAGEME"; } 
	public final String getLabel() { return "MP2RAGEME PD Mapping"; }
	public final String getName() { return "MP2RAGEME PD Mappng"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Integrative Model-based Cognitive Neuroscience unit, Universiteit van Amsterdam | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimation of approximate proton density maps from the MP2RAGEME sequence, based on the model of Marques et al., 2010."; }
		
	public final String getVersion() { return "3.2"; }

	// get outputs
	public float[] getProtonDensityImage1() { return pd1map; }
	public float[] getProtonDensityImage2() { return pd2map; }
	public float[] getProtonDensityImage0() { return pd0map; }
	public float[] getProtonDensityImage() { return pdmap; }
	
	
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm
		System.out.print("Loading data\n");
				
		// we assume 4D images: dim 0 magnitude, dim 1 phase, if used
		if (uni!=null) {
		    usePhase = false;
		}
		
		// estimate angular scale (range should be [-PI +PI]
		double phscale1 = 1.0;
		double phscale2 = 1.0;
		if (usePhase && scalePhase) {
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
		double TA = TI1 - Nexcitations/2.0f*TRexcitation1;
		double TB = TI2 - TI1 - Nexcitations/2.0f*TRexcitation1 - Nexcitations/2.0f*TRexcitation2;
		double TC = TRinversion - TI2 - Nexcitations/2.0f*TRexcitation2;
		
		System.out.print("timings: TA = "+TA+", TB = "+TB+", TC = "+TC+"\n"); 
		
		double a1rad = angle1/180.0*FastMath.PI;
		double a2rad = angle2/180.0*FastMath.PI;
		
		pd1map = new float[nxyz];
		pd2map = new float[nxyz];
		if (s0img!=null) pd0map = new float[nxyz];
 		pdmap = new float[nxyz];
		
		for (int xyz=0;xyz<nxyz;xyz++) {
		    
		    // assuming T1 in seconds
		    //double qt1 = t1map[xyz]/1000.0;
		    double qt1 = t1map[xyz];
		    
		    // assuming R2* in Hz
		    double expr2s = 1.0;
		    //if (r2smap!=null) expr2s = FastMath.exp(-r2smap[xyz]*TEcho1*1000.0);
		    if (r2smap!=null) expr2s = FastMath.exp(-r2smap[xyz]*TEcho1);
		
			double b1 = 1.0;
			if (useB1correction) b1 = b1map[xyz];
			
			double cosA1cosA2n = FastMath.pow( FastMath.cos(b1*a1rad)*FastMath.cos(b1*a2rad), Nexcitations );
		
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
            
            double prod;
            if (usePhase)
                prod = inv1[xyz]*inv2[xyz]*FastMath.cos(inv1[xyz+nxyz]/phscale1-inv2[xyz+nxyz]/phscale2);
            else
                prod = uni[xyz]*(inv1[xyz]*inv1[xyz]+inv2[xyz]*inv2[xyz]);
            
            pdmap[xyz] = (float)FastMath.sqrt(Numerics.max(0.0, (prod/(gre1*gre2*expr2s*expr2s)) ) );
            
            //pd1map[xyz] = (float)(inv1[xyz]*FastMath.cos(inv1[xyz+nxyz]/phscale1)/gre1/expr2s);
            //pd2map[xyz] = (float)(inv2[xyz]*FastMath.cos(inv2[xyz+nxyz]/phscale2)/gre2/expr2s);
            pd1map[xyz] = (float)(inv1[xyz]/gre1/expr2s);
            pd2map[xyz] = (float)(inv2[xyz]/gre2/expr2s);
            
            if (s0img!=null) pd0map[xyz] = (float)(s0img[xyz]/gre2);
            /*
            pdmap[xyz] = pd1map[xyz]+pd2map[xyz];
            if (s0img!=null) {
                pdmap[xyz] += pd0map[xyz];
                pdmap[xyz] /= 3.0f;
            } else {
                pdmap[xyz] /= 2.0f;
            }*/
 		}
 			
		return;
	}

}
