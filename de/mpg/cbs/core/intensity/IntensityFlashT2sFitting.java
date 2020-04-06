package de.mpg.cbs.core.intensity;

import de.mpg.cbs.libraries.*;
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
	float r2max = 1000.0f;
	
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

	public final void setMaxR2s(float val) { r2max = val; }
	
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
				r2img[xyz] = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
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
	
	public void minEchoEstimation() {
		
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
                // fit the exponential
                double Sx, Sx2, Sy, Sxy, delta;
                Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                for (int i=0;i<nimg;i++) {
                    Sx += -te[i];
                    Sx2 += Numerics.square(-te[i]);
                }
                delta = nimg*Sx2 - Sx*Sx;
				Sy = 0.0f;
				Sxy = 0.0f;
				for (int i=0;i<nimg;i++) {
					double ydata = FastMath.log(image[i][xyz]);
					Sy += ydata;
					Sxy += -ydata*te[i];
				}
				//s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
				s0img[xyz] = (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta );
				r2img[xyz] = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
				t2img[xyz] = 1.0f/r2img[xyz];
				
				// check if over the boundary: if so, reduce the echo train
				int necho = nimg-1;
				while (necho>1) {
				    Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                    for (int i=0;i<necho;i++) {
                        Sx += -te[i];
                        Sx2 += Numerics.square(-te[i]);
                    }
                    delta = necho*Sx2 - Sx*Sx;
                    Sy = 0.0f;
                    Sxy = 0.0f;
                    for (int i=0;i<necho;i++) {
                        double ydata = FastMath.log(image[i][xyz]);
                        Sy += ydata;
                        Sxy += -ydata*te[i];
                    }
                    float r2val = Numerics.bounded( (float)( (necho*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
                    if (r2val<r2img[xyz]) {
                        r2img[xyz] = r2val;
                        //s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
                        s0img[xyz] = (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta );
                        t2img[xyz] = 1.0f/r2img[xyz];
                    }
				    necho = necho-1;
				}   
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

	public void variableEchoEstimation() {
		
		System.out.print("Loading data\n");
		
		System.out.print("Echo times: [ ");
		for (int i=0;i<nimg;i++) System.out.print(te[i]+" ");
		System.out.print("]\n");
		
		// step 1: compute the difference between full R2* fit and minimum
		float[] diff = new float[nxyz];
		boolean[] mask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			mask[xyz] = false;
		    boolean process=true;
			for (int i=0;i<nimg;i++) if (image[i][xyz]<=0) process=false;
			if (process) {
			    //System.out.print(".");
                // fit the exponential
                double Sx, Sx2, Sy, Sxy, delta;
                Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                for (int i=0;i<nimg;i++) {
                    Sx += -te[i];
                    Sx2 += Numerics.square(-te[i]);
                }
                delta = nimg*Sx2 - Sx*Sx;
				Sy = 0.0f;
				Sxy = 0.0f;
				for (int i=0;i<nimg;i++) {
					double ydata = FastMath.log(image[i][xyz]);
					Sy += ydata;
					Sxy += -ydata*te[i];
				}
				//s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
				float r2full = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
				float r2min = r2full;
				
				// mask very high values by default
				if (r2full>2.0*r2max) mask[xyz] = true;
				
				// check if over the boundary: if so, reduce the echo train
				int necho = nimg-1;
				while (necho>1) {
				    Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                    for (int i=0;i<necho;i++) {
                        Sx += -te[i];
                        Sx2 += Numerics.square(-te[i]);
                    }
                    delta = necho*Sx2 - Sx*Sx;
                    Sy = 0.0f;
                    Sxy = 0.0f;
                    for (int i=0;i<necho;i++) {
                        double ydata = FastMath.log(image[i][xyz]);
                        Sy += ydata;
                        Sxy += -ydata*te[i];
                    }
                    float r2val = Numerics.bounded( (float)( (necho*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
                    if (r2val<r2min) r2min = r2val;
				    necho = necho-1;
				}
				diff[xyz] = r2full-r2min;
			}
		}
		// step 2: look for neighborhoods
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    if (diff[xyz]>r2max) {
                int neighbors=0;
                for (byte n=0;n<26;n++) {
                    int ngb = Ngb.neighborIndex(n, xyz, nx, ny, nz);
                    if (diff[ngb]>r2max) neighbors++;
                }
                if (neighbors>8) mask[xyz] = true;
            }   
		}
		// step 3: find nearby voxels with higher R2* value
		
		s0img = new float[nxyz];
		r2img = new float[nxyz];
		t2img = new float[nxyz];
		errimg = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			boolean process=true;
			for (int i=0;i<nimg;i++) if (image[i][xyz]<=0) process=false;
			if (process) {
			    //System.out.print(".");
                // fit the exponential
                double Sx, Sx2, Sy, Sxy, delta;
                Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                for (int i=0;i<nimg;i++) {
                    Sx += -te[i];
                    Sx2 += Numerics.square(-te[i]);
                }
                delta = nimg*Sx2 - Sx*Sx;
				Sy = 0.0f;
				Sxy = 0.0f;
				for (int i=0;i<nimg;i++) {
					double ydata = FastMath.log(image[i][xyz]);
					Sy += ydata;
					Sxy += -ydata*te[i];
				}
				//s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
				float r2full = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
				float s0full = (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta );
				r2img[xyz] = r2full;
				t2img[xyz] = 1.0f/r2full;
				s0img[xyz] = s0full;
			}
		}
		// expand mask on neighbors with high values
		int nmask = ObjectStatistics.volume(mask, nx,ny,nz);
        int added = nmask;
		while (added>0.01*nmask) {
		    System.out.print(".");
            boolean[] expand = new boolean[nxyz];
            for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                int xyz = x+nx*y+nx*ny*z;
                if (mask[xyz]) {
                    int neighbors=0;
                    for (byte n=0;n<6;n++) {
                        int ngb = Ngb.neighborIndex(n, xyz, nx, ny, nz);
                        if (r2img[ngb]>r2img[xyz]) expand[ngb] = true;
                    }
                }   
            }
            added = 0;
            for (int xyz=0;xyz<nxyz;xyz++) {
                if (!mask[xyz] && expand[xyz]) added++;
                mask[xyz] = (mask[xyz] || expand[xyz]);
            }
        }		
		// find minimum again
		for (int xyz=0;xyz<nxyz;xyz++) {
			boolean process=true;
			for (int i=0;i<nimg;i++) if (image[i][xyz]<=0) process=false;
			if (process) {
				if (mask[xyz]) {
                    // check if over the boundary: if so, reduce the echo train
                    int necho = nimg-1;
                    while (necho>1) {
                        double Sx, Sx2, Sy, Sxy, delta;
                        Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
                        for (int i=0;i<necho;i++) {
                            Sx += -te[i];
                            Sx2 += Numerics.square(-te[i]);
                        }
                        delta = necho*Sx2 - Sx*Sx;
                        Sy = 0.0f;
                        Sxy = 0.0f;
                        for (int i=0;i<necho;i++) {
                            double ydata = FastMath.log(image[i][xyz]);
                            Sy += ydata;
                            Sxy += -ydata*te[i];
                        }
                        float r2val = Numerics.bounded( (float)( (necho*Sxy-Sx*Sy)/delta ), 0.001f, 1000.0f);
                        if (r2val<r2img[xyz]) {
                            r2img[xyz] = r2val;
                            //s0img[xyz] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, imax);
                            s0img[xyz] = (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta );
                            t2img[xyz] = 1.0f/r2img[xyz];
                        }
                        necho = necho-1;
                    }   
                }
            }
        }
        // map residuals
		for (int xyz=0;xyz<nxyz;xyz++) {
			boolean process=true;
			for (int i=0;i<nimg;i++) if (image[i][xyz]<=0) process=false;
			if (process) {
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
