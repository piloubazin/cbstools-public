package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.*;

import Jama.*;

/*
 * @author Pierre-Louis Bazin
 */
public class IntensityMp2ragemePCADenoising {

	// input parameters
	private		float[] 	invmag = null;
	private		float[] 	invphs = null;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private     int         nimg = 5;
	private		float		stdevCutoff = 1.0f;
	private     int         minDimension = 2;
	private     int         maxDimension = -1;
	private     int         ngbSize = 5;
	private     boolean     separate = false;
	private     boolean     tvmag = false;
	private     boolean     tvphs = false;
	
	public static final String[]     thresholdingTypes = {"Eigenvalues","Global noise","Second half"};
	private     String      thresholding = "Second half";
		
	// output parameters
	private		float[] invmagden = null;
	private		float[] invphsden = null;
	private		float[] pcadim = null;
	private		float[] errmap = null;
	
	// internal variables (for debug)
	private float[][] images;
	private float[][] denoised;
	private float[][] eigvec;
	private float[][] eigval;
	
	// set inputs
	public final void setMagnitudeImages(float[] in) { invmag = in; }
	public final void setPhaseImages(float[] in) { invphs = in; }

	public final void setFirstInversionImage(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
            invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz] = in[xyz];
	        invphs[xyz] = in[xyz+nxyz];
	    }
	}

	public final void setSecondInversionEcho1Image(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+1*nxyz] = in[xyz];
	        invphs[xyz+1*nxyz] = in[xyz+nxyz];
	    }
	}

	public final void setSecondInversionEcho2Image(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+2*nxyz] = in[xyz];
	        invphs[xyz+2*nxyz] = in[xyz+nxyz];
	    }
	}

	public final void setSecondInversionEcho3Image(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+3*nxyz] = in[xyz];
	        invphs[xyz+3*nxyz] = in[xyz+nxyz];
	    }
	}

	public final void setSecondInversionEcho4Image(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+4*nxyz] = in[xyz];
	        invphs[xyz+4*nxyz] = in[xyz+nxyz];
	    }
	}

	public final void setSecondInversionImage(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/8;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+1*nxyz] = in[xyz];
	        invmag[xyz+2*nxyz] = in[xyz+nxyz];
	        invmag[xyz+3*nxyz] = in[xyz+2*nxyz];
	        invmag[xyz+4*nxyz] = in[xyz+3*nxyz];
	        
	        invphs[xyz+1*nxyz] = in[xyz+4*nxyz];
	        invphs[xyz+2*nxyz] = in[xyz+5*nxyz];
	        invphs[xyz+3*nxyz] = in[xyz+6*nxyz];
	        invphs[xyz+4*nxyz] = in[xyz+7*nxyz];
	    }
	}

	public final void setMagnitudeAndPhaseImage(int num, float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        invmag[xyz+(num-1)*nxyz] = in[xyz];
	        invphs[xyz+(num-1)*nxyz] = in[xyz+nxyz];
	    }
	}
	
	public final void setImageNumber(int in) { nimg = in; }
	public final void setStdevCutoff(float in) { stdevCutoff = in; }
	public final void setMinimumDimension(int in) { minDimension = in; }
	public final void setMaximumDimension(int in) { maxDimension = in; }
	public final void setPatchSize(int in) { ngbSize = in; }
	public final void setSeparateProcessing(boolean in) { separate = in; }
	
	public final void setMagnitudeTVSubtraction(boolean in) { tvmag = in; }
	public final void setPhaseTVSubtraction(boolean in) { tvphs = in; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.devel"; } 
	public final String getLabel() { return "MP2RAGEME PCA denoising"; }
	public final String getName() { return "MP2RAGEME PCA denoising"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Denoise the data with a PCA-based method (adapted from Manjon et al., Plos One 2013."; }
		
	public final String getVersion() { return "3.1.3"; }

	// get outputs
	public float[] getDenoisedMagnitudeImages() { return invmag; }
	public float[] getDenoisedPhaseImages() { return invphs; }
	public float[] getLocalDimensionImage() { return pcadim; }
	public float[] getNoiseMagnitudeImage() { return errmap; }
	
	
	public float[] getFirstInversionImage() {
	    float[] out = new float[2*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz] = invmag[xyz];
	        out[xyz+nxyz] = invphs[xyz];
	    }
	    return out;
	}
	
	public float[] getSecondInversionEcho1Image() {
	    float[] out = new float[2*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz] = invmag[xyz+1*nxyz];
	        out[xyz+nxyz] = invphs[xyz+1*nxyz];
	    }
	    return out;
	}
	
	public float[] getSecondInversionEcho2Image() {
	    float[] out = new float[2*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz] = invmag[xyz+2*nxyz];
	        out[xyz+nxyz] = invphs[xyz+2*nxyz];
	    }
	    return out;
	}
	
	public float[] getSecondInversionEcho3Image() {
	    float[] out = new float[2*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz] = invmag[xyz+3*nxyz];
	        out[xyz+nxyz] = invphs[xyz+3*nxyz];
	    }
	    return out;
	}
	
	public float[] getSecondInversionEcho4Image() {
	    float[] out = new float[2*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz] = invmag[xyz+4*nxyz];
	        out[xyz+nxyz] = invphs[xyz+4*nxyz];
	    }
	    return out;
	}
	
	public float[] getSecondInversionImage() {
	    float[] out = new float[8*nxyz];
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz+0*nxyz] = invmag[xyz+1*nxyz];
	        out[xyz+1*nxyz] = invmag[xyz+2*nxyz];
	        out[xyz+2*nxyz] = invmag[xyz+3*nxyz];
	        out[xyz+3*nxyz] = invmag[xyz+4*nxyz];
	        out[xyz+4*nxyz] = invphs[xyz+1*nxyz];
	        out[xyz+5*nxyz] = invphs[xyz+2*nxyz];
	        out[xyz+6*nxyz] = invphs[xyz+3*nxyz];
	        out[xyz+7*nxyz] = invphs[xyz+4*nxyz];
	    }
	    return out;
	}
	
	public float[] getRawComplexImage() {
	    float[] out = new float[2*nimg*nxyz];
	    for (int i=0;i<2*nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz+i*nxyz] = images[i][xyz];
	    }
	    return out;
	}

	public float[] getDenComplexImage() {
	    float[] out = new float[2*nimg*nxyz];
	    for (int i=0;i<2*nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz+i*nxyz] = denoised[i][xyz];
	    }
	    return out;
	}

	public float[] getEigenvectorImage() {
	    float[] out = new float[2*nimg*nxyz];
	    for (int i=0;i<2*nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz+i*nxyz] = eigvec[i][xyz];
	    }
	    return out;
	}

	public float[] getEigenvalueImage() {
	    float[] out = new float[2*nimg*nxyz];
	    for (int i=0;i<2*nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
	        out[xyz+i*nxyz] = eigval[i][xyz];
	    }
	    return out;
	}
	
	public void execute() {
	    //if (separate) denoiseSeparately();
	    //else denoiseJointly();
	//}

	//public void denoiseJointly() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// we assume 4D images of size nimg
		if (invmag==null || invphs==null) System.out.println("data stacks not properly initialized!");
		
		// renormalize phase
		float phsmin = invphs[0];
		float phsmax = invphs[0];
        for (int xyz=0;xyz<nimg*nxyz;xyz++) {
			if (invphs[xyz]<phsmin) phsmin = invphs[xyz];
			if (invphs[xyz]>phsmax) phsmax = invphs[xyz];
		}
		double phsscale = (phsmax-phsmin)/(2.0*FastMath.PI);
		
		// opt. remove TV global variations
		float[][] tvimgmag = null;
		if (tvmag) {
		    tvimgmag = new float[nimg][];
		    float[] tmp = new float[nxyz];
		    for (int i=0;i<nimg;i++) {
                System.out.println("global variations removal magnitude "+(i+1));
                for (int xyz=0;xyz<nxyz;xyz++) tmp[xyz] = invmag[xyz+i*nxyz];
                // separate the magnitude images
                TotalVariation1D algo = new TotalVariation1D(tmp,null,nx,ny,nz, 0.33f, 0.125f, 0.001f, 100);
                algo.solve();
                tvimgmag[i] = algo.exportResult();
                for (int xyz=0;xyz<nxyz;xyz++) invmag[xyz+i*nxyz] -= tvimgmag[i][xyz];
            }
        }
        float[][] tvimgphs = null;
		if (tvphs) {
		    tvimgphs = new float[nimg][];
		    float[] tmp = new float[nxyz];
		    for (int i=0;i<nimg;i++) {
               System.out.println("global variations removal phase "+(i+1));
                 for (int xyz=0;xyz<nxyz;xyz++) tmp[xyz] = invphs[xyz+i*nxyz];
                // separate the magnitude images
                TotalVariation1D algo = new TotalVariation1D(tmp,null,nx,ny,nz, 0.5f, 0.125f, 0.001f, 100);
                algo.solveWrapped();
                tvimgphs[i] = algo.exportResultWrapped();
                for (int xyz=0;xyz<nxyz;xyz++) invphs[xyz+i*nxyz] -= tvimgphs[i][xyz];
            }
        }
		
		// 1. create all the sin, cos images
		images = new float[2*nimg][nxyz];
		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                images[2*i][xyz] = (float)(invmag[xyz+i*nxyz]*FastMath.cos(invphs[xyz+i*nxyz]/phsscale));
                images[2*i+1][xyz] = (float)(invmag[xyz+i*nxyz]*FastMath.sin(invphs[xyz+i*nxyz]/phsscale));
            }
        }
        invmag = null;
        invphs = null;
		
		// 2. estimate PCA in slabs of NxNxN size
		int ngb = ngbSize;
		int nstep = Numerics.floor(ngb/2.0);
		int nimg2 = 2*nimg;
		denoised = new float[nimg2][nxyz];
		eigvec = new float[nimg2][nxyz];
		eigval = new float[nimg2][nxyz];
		float[] weights = new float[nxyz];
		pcadim = new float[nxyz];
		errmap = new float[nxyz];
		// border issues should be cleaned-up, ignored so far
		for (int x=0;x<nx;x+=nstep) for (int y=0;y<ny;y+=nstep) for (int z=0;z<nz;z+=nstep) {
		    int ngbx = Numerics.min(ngb, nx-x);
		    int ngby = Numerics.min(ngb, ny-y);
		    int ngbz = Numerics.min(ngb, nz-z);
		    int ngb3 = ngbx*ngby*ngbz;
		    if (ngb3<nimg2) {
		        System.out.println("!patch is too small!");
		    } else {
                double[][] patch = new double[ngb3][nimg2];
                for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) for (int i=0;i<nimg2;i++) {
                    patch[dx+ngbx*dy+ngbx*ngby*dz][i] = images[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)];
                }
                // mean over samples
                double[] mean = new double[nimg2];
                for (int i=0;i<nimg2;i++) {
                   for (int n=0;n<ngb3;n++) mean[i] += patch[n][i];
                   mean[i] /= (double)ngb3;
                   for (int n=0;n<ngb3;n++) patch[n][i] -= mean[i];
                }
                // PCA from SVD X = USVt
                //System.out.println("perform SVD");
                Matrix M = new Matrix(patch);
                SingularValueDecomposition svd = M.svd();
            
                // estimate noise
                // simple version: compute the standard deviation of the patch
                double sigma = 0.0;
                for (int n=0;n<ngb3;n++) for (int i=0;i<nimg2;i++) {
                    sigma += patch[n][i]*patch[n][i];
                }
                sigma /= nimg2*ngb3;
                sigma = FastMath.sqrt(sigma);
                
                // cutoff
                //System.out.println("eigenvalues: ");
                double[] eig = new double[nimg2];
                boolean[] used = new boolean[nimg2];
                int nzero=0;
                double eigsum = 0.0;
                for (int n=0;n<nimg2;n++) {
                    eig[n] = svd.getSingularValues()[n];
                    eigsum += Numerics.abs(eig[n]);
                }
                // fit second half linear decay model
                double[] loc = new double[nimg];
                double[][] fit = new double[nimg][1];
                for (int n=nimg;n<nimg2;n++) {
                    loc[n-nimg] = (n-nimg)/(double)nimg;
                    fit[n-nimg][0] = Numerics.abs(eig[n]);
                }
                double[][] poly = new double[nimg][2];
                for (int n=0;n<nimg;n++) {
                    poly[n][0] = 1.0;
                    poly[n][1] = loc[n];
                }
                // invert the linear model
                Matrix mtx = new Matrix(poly);
                Matrix smp = new Matrix(fit);
                Matrix val = mtx.solve(smp);
		
                /*
                // use robust measurement?
                // Teil-Shen estimator ? not working properly...
                double[] slopes = new double[nimg*(nimg-1)];
                int s=0;
                for (int n=nimg;n<nimg2;n++) for (int m=n+1;m<nimg2;m++) {
                    slopes[s] = (Numerics.abs(eig[n]) - Numerics.abs(eig[m]))/(double)(n-m);
                    s++;
                }
                Percentile measure = new Percentile();
				double slope = measure.evaluate(slopes, 50.0);
				double[] intercepts = new double[nimg];
				for (int n=nimg;n<nimg2;n++) {
				    intercepts[n-nimg] = Numerics.abs(eig[n]) - n*slope;
				}
				double intercept = measure.evaluate(intercepts, 50.0);
				*/
                // compute the expected value:
                double[] expected = new double[nimg2];
                for (int n=0;n<nimg2;n++) {
                    double n0 = (n-nimg)/(double)nimg;
                    // linear coeffs,
                    expected[n] = (val.get(0,0) + n0*val.get(1,0));
                    //expected[n] = n*slope + intercept;
                }
                
				
                for (int n=0;n<nimg2;n++) {
                    //System.out.print(" "+(eig[n]/sigma));
                    //if (n>=minDimension && nimg2*Numerics.abs(eig[n]) < stdevCutoff*sigma) {
                    //if (n>=minDimension && nimg2*Numerics.abs(eig[n]) < stdevCutoff*eigsum) {
                    if (n>=minDimension && Numerics.abs(eig[n]) < stdevCutoff*expected[n]) {
                        used[n] = false;
                        nzero++;
                        //System.out.print("(-),");
                    } else  if (maxDimension>0 && n>=maxDimension) {
                        used[n] = false;
                        nzero++;
                        //System.out.print("(|),");
                    } else {
                        used[n] = true;
                        //System.out.print("(+),");
                    }
                }
                // reconstruct
                Matrix U = svd.getU();
                Matrix V = svd.getV();
                for (int n=0;n<ngb3;n++) for (int i=0;i<nimg2;i++) {
                    // Sum_s>0 s_kU_kV_kt
                    patch[n][i] = mean[i];
                    for (int j=0;j<nimg2;j++) if (used[j]) {
                        patch[n][i] += U.get(n,j)*eig[j]*V.get(i,j);
                    }
                }
                // add to the denoised image
                double wpatch = (1.0/(1.0 + nimg2 - nzero));
                for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) {
                    for (int i=0;i<nimg2;i++) {
                        denoised[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*patch[dx+ngbx*dy+ngbx*ngby*dz][i]);
                        eigval[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*eig[i]);
                        //eigvec[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*U.get(dx+ngbx*dy+ngbx*ngby*dz,i));
                        eigvec[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*expected[i]);
                    }
                    weights[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                    pcadim[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(nimg2-nzero));
                    errmap[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*val.get(1,0));
                    //errmap[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*slope);
                }
            }
        }
		//errmap = new float[nxyz];
        for (int xyz=0;xyz<nxyz;xyz++) {
            double err = 0.0;
            for (int i=0;i<nimg2;i++) {
                denoised[i][xyz] /= weights[xyz];
                eigval[i][xyz] /= weights[xyz];
                eigvec[i][xyz] /= weights[xyz];
                err += (denoised[i][xyz]-images[i][xyz])*(denoised[i][xyz]-images[i][xyz]);
            }
            pcadim[xyz] /= weights[xyz];
            errmap[xyz] /= weights[xyz];
           // errmap[xyz] = (float)FastMath.sqrt(err/nimg2);
        }
        //images = null;
          
        // 3. rebuild magnitude and phase images
        invmag = new float[nimg*nxyz];
        invphs = new float[nimg*nxyz];
  		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                invmag[xyz+i*nxyz] = (float)FastMath.sqrt(denoised[2*i][xyz]*denoised[2*i][xyz]+denoised[2*i+1][xyz]*denoised[2*i+1][xyz]);
                invphs[xyz+i*nxyz] = (float)(FastMath.atan2(denoised[2*i+1][xyz],denoised[2*i][xyz])*phsscale);
             }
        }
        
        // opt. add back the TV estimate, if needed
        if (tvmag) {
            for (int i=0;i<nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
                invmag[xyz+i*nxyz] += tvimgmag[i][xyz];
            }
        }
        if (tvphs) {
            for (int i=0;i<nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
                invphs[xyz+i*nxyz] += tvimgphs[i][xyz];
            }
        }
		return;
	}

}
