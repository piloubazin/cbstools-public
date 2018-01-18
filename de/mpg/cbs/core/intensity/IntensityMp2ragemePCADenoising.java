package de.mpg.cbs.core.intensity;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

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
	private		float		stdevCutoff = 2.3f;
		
	// output parameters
	private		float[] invmagden = null;
	private		float[] invphsden = null;
	private		float[] pcadim = null;
	
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

	public final void setImageNumber(int in) { nimg = in; }
	public final void setStdevCutoff(float in) { stdevCutoff = in; }
	
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
	
	public void execute() {
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
		
		// 1. create all the sin, cos images
		float[][] images = new float[2*nimg][nxyz];
		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                images[2*i][xyz] = (float)(invmag[xyz+i*nxyz]*FastMath.cos(invphs[xyz+i*nxyz]/phsscale));
                images[2*i+1][xyz] = (float)(invmag[xyz+i*nxyz]*FastMath.sin(invphs[xyz+i*nxyz]/phsscale));
            }
        }
        invmag = null;
        invphs = null;
		
		// 2. estimate PCA in slabs of NxNxN size
		int ngb = 4;
		int nstep = 2;
		int ngb3 = ngb*ngb*ngb;
		int nimg2 = 2*nimg;
		float[][] denoised = new float[2*nimg][nxyz];
		float[] weights = new float[nxyz];
		pcadim = new float[nxyz];
		// border issues should be cleaned-up, ignored so far
		for (int x=0;x<nx-ngb;x+=nstep) for (int y=0;y<ny-ngb;y+=nstep) for (int z=0;z<nz-ngb;z+=nstep) {
		    double[][] patch = new double[ngb3][nimg2];
		    for (int dx=0;dx<ngb;dx++) for (int dy=0;dy<ngb;dy++) for (int dz=0;dz<ngb;dz++) for (int i=0;i<nimg2;i++) {
		        patch[dx+ngb*dy+ngb*ngb*dz][i] = images[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)];
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
            int nzero=0;
            for (int n=0;n<nimg2;n++) {
                eig[n] = svd.getSingularValues()[n];
                //System.out.print(" "+(eig[n]/sigma));
                if (Numerics.abs(eig[n]) < stdevCutoff*sigma) {
                    eig[n] = 0.0;
                    nzero++;
                    //System.out.print("(-),");
                } else {
                    //System.out.print("(+),");
                }
            }
            // reconstruct
            Matrix U = svd.getU();
            Matrix V = svd.getV();
            for (int n=0;n<ngb3;n++) for (int i=0;i<nimg2;i++) {
                // Sum_s>0 s_kU_kV_kt
                patch[n][i] = mean[i];
                for (int j=0;j<nimg2;j++) if (eig[j]!=0) {
                    patch[n][i] += U.get(n,j)*eig[j]*V.get(i,j);
                }
            }
            // add to the denoised image
            double wpatch = (1.0/(1.0 + nimg2 - nzero));
            for (int dx=0;dx<ngb;dx++) for (int dy=0;dy<ngb;dy++) for (int dz=0;dz<ngb;dz++) {
                for (int i=0;i<nimg2;i++) {
                    denoised[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*patch[dx+ngb*dy+ngb*ngb*dz][i]);
                }
                weights[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                pcadim[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(nimg2-nzero));
		    }
        }
        images = null;
        for (int xyz=0;xyz<nxyz;xyz++) {
            for (int i=0;i<nimg;i++) {
                denoised[i][xyz] /= weights[xyz];
            }
            pcadim[xyz] /= weights[xyz];
        }
           
        // 3. rebuild magnitude and phase images
        invmag = new float[nimg*nxyz];
        invphs = new float[nimg*nxyz];
  		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                invmag[xyz+i*nxyz] = (float)FastMath.sqrt(denoised[2*i][xyz]*denoised[2*i][xyz]+denoised[2*i+1][xyz]*denoised[2*i+1][xyz]);
                invphs[xyz+i*nxyz] = (float)(FastMath.atan2(denoised[2*i+1][xyz],denoised[2*i][xyz])*phsscale);
             }
        }
              
		return;
	}

}
