package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class IntensityFlashT2sFitting {
	float[][] image = null;
	float[] te;
	
	int nx, ny, nz, nxyz;
	float rx, ry, rz;
	
	int nimg = 4;
	float imax = 10000.0f;
	
	float[] s0img;
	float[] t2img;
	float[] r2img;
	float[] errimg;
	
	// set inputs
	public final void setNumberOfEchoes(int val) { 
	    nimg = val;
	    image = new float[nimg][];
	    te = new float[nimg];
	}
	public final void setEchoImageAt(int num, float[] val) { image[num] = val; }
	public final void setEchoTimeAt(int num, float val) { te[num] = val; }

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.devel"; } 
	public final String getLabel() { return "Flash T2* Fitting"; }
	public final String getName() { return "FlashT2sFitting"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Estimation of T2* maps from linear least squares in log space."; }
		
	public final String getVersion() { return "3.1.3"; }

	// outputs
	public final float[] getS0Image() { return s0img; }
	public final float[] getT2sImage() { return t2img; }
	public final float[] getR2sImage() { return r2img; }
	public final float[] getResidualImage() { return errimg; }

	public void execute() {
	    // fit the exponential
		double Sx, Sx2, Sy, Sxy, delta;
		Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
		for (int i=0;i<nimg;i++) {
			Sx += -te[i];
			Sx2 += Numerics.square(-te[i]);
		}
		delta = nimg*Sx2 - Sx*Sx;
		
		System.out.print("Loading data\n");
		
		System.out.print("Echo times: [ ");
		for (int i=0;i<nimg;i++) System.out.print(te[i]+" ");
		System.out.print("]\n");
		
		s0img = new float[nxyz];
		r2img = new float[nxyz];
		t2img = new float[nxyz];
		errimg = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			boolean process=true;
			for (int i=0;i<nimg;i++) if (image[i][xyz]<=0) process=false;
			if (process) {
			    //System.out.print(".");
				Sy = 0.0f;
				Sxy = 0.0f;
				for (int i=0;i<nimg;i++) {
					double ydata = FastMath.log(image[i][xyz]);
					Sy += ydata;
					Sxy += -ydata*te[i];
				}
				//s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
				s0img[xyz] = (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta );
				r2img[xyz] = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1.0f);
				t2img[xyz] = 1.0f/r2img[xyz];
				/*
				errimg[xyz] = 0.0f;
				for (int i=0;i<nimg;i++) {
					errimg[xyz] += Numerics.square( FastMath.log(image[i][xyz]) 
													 - FastMath.log(s0img[xyz]) + te[i]*r2img[xyz] )/nimg;
				}
				errimg[xyz] = (float)FastMath.sqrt(errimg[xyz]);
				*/
				double residual = 0.0;
                double mean = 0.0;
                double variance = 0.0;
                for (int i=0;i<nimg;i++) mean += FastMath.log(image[i][xyz]);
                mean /= nimg;
                for (int i=0;i<nimg;i++) {
                    double expected = FastMath.log(s0img[xyz]) - te[i]*r2img[xyz];
                    variance += Numerics.square(mean-FastMath.log(image[i][xyz]));
                    residual += Numerics.square(expected-FastMath.log(image[i][xyz]));
                }
                double rsquare = 1.0;
                if (variance>0) rsquare = Numerics.max(1.0 - (residual/variance), 0.0);
                errimg[xyz] = (float)rsquare;        
			}
		}
		System.out.print("Done\n");
		
	}
}
