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
public class IntensityGenericPCADenoising {

	// input parameters
	private		float[][] 	images = null;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private     int         nimg = 5;
	private		float		stdevCutoff = 2.3f;
	private     int         minDimension = 1;
	private     int         ngbSize = 4;
		
	// output parameters
	private		float[][] denoised = null;
	private		float[] pcadim = null;
	private		float[] errmap = null;
	
	// set inputs
	public final void setImageAt(int n, float[] in) {
	    if (images==null) {
	        if (nxyz==0) nxyz = in.length;
            images = new float[nimg][nxyz];
	    }
	    for (int xyz=0;xyz<nxyz;xyz++) {
	        images[n][xyz] = in[xyz];
	    }
	}
	
	public final void setImageNumber(int in) { nimg = in; }
	public final void setStdevCutoff(float in) { stdevCutoff = in; }
	public final void setMinimumDimension(int in) { minDimension = in; }
	public final void setPatchSize(int in) { ngbSize = in; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.devel"; } 
	public final String getLabel() { return "Generic PCA denoising"; }
	public final String getName() { return "Generic PCA denoising"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Denoise the data with a PCA-based method (adapted from Manjon et al., Plos One 2013."; }
		
	public final String getVersion() { return "3.1.3"; }

	// get outputs
	public float[] getDenoisedImageAt(int n, int m) { 
	    if (m-n==1) return denoised[n];
	    else {
            float[] out = new float[(m-n)*nxyz];
            for (int c=n;c<m;c++) for (int xyz=0;xyz<nxyz;xyz++) {
                out[xyz+c*nxyz] = denoised[n+c][xyz];
            }
            return out;
        }
	}
	public float[] getLocalDimensionImage() { return pcadim; }
	public float[] getNoiseMagnitudeImage() { return errmap; }
		
	public void execute() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// we assume 4D images of size nimg
		if (images==null) System.out.println("data stacks not properly initialized!");
				
		// 2. estimate PCA in slabs of NxNxN size
		int ngb = ngbSize;
		int nstep = Numerics.floor(ngb/2.0);
		denoised = new float[nimg][nxyz];
		float[] weights = new float[nxyz];
		pcadim = new float[nxyz];
		// border issues should be cleaned-up, ignored so far
		for (int x=0;x<nx;x+=nstep) for (int y=0;y<ny;y+=nstep) for (int z=0;z<nz;z+=nstep) {
		    int ngbx = Numerics.min(ngb, nx-x);
		    int ngby = Numerics.min(ngb, ny-y);
		    int ngbz = Numerics.min(ngb, nz-z);
		    int ngb3 = ngbx*ngby*ngbz;
		    if (ngb3<nimg) {
		        System.out.println("!patch is too small!");
		    } else {
                double[][] patch = new double[ngb3][nimg];
                for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) for (int i=0;i<nimg;i++) {
                    patch[dx+ngbx*dy+ngbx*ngby*dz][i] = images[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)];
                }
                // mean over samples
                double[] mean = new double[nimg];
                for (int i=0;i<nimg;i++) {
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
                for (int n=0;n<ngb3;n++) for (int i=0;i<nimg;i++) {
                    sigma += patch[n][i]*patch[n][i];
                }
                sigma /= nimg*ngb3;
                sigma = FastMath.sqrt(sigma);
                
                // cutoff
                //System.out.println("eigenvalues: ");
                double[] eig = new double[nimg];
                int nzero=0;
                for (int n=0;n<nimg;n++) {
                    eig[n] = svd.getSingularValues()[n];
                    //System.out.print(" "+(eig[n]/sigma));
                    if (n>=minDimension && Numerics.abs(eig[n]) < stdevCutoff*sigma) {
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
                for (int n=0;n<ngb3;n++) for (int i=0;i<nimg;i++) {
                    // Sum_s>0 s_kU_kV_kt
                    patch[n][i] = mean[i];
                    for (int j=0;j<nimg;j++) if (eig[j]!=0) {
                        patch[n][i] += U.get(n,j)*eig[j]*V.get(i,j);
                    }
                }
                // add to the denoised image
                double wpatch = (1.0/(1.0 + nimg - nzero));
                for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) {
                    for (int i=0;i<nimg;i++) {
                        denoised[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*patch[dx+ngbx*dy+ngbx*ngby*dz][i]);
                    }
                    weights[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                    pcadim[x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(nimg-nzero));
                }
            }
        }
		errmap = new float[nxyz];
        for (int xyz=0;xyz<nxyz;xyz++) {
            double err = 0.0;
            for (int i=0;i<nimg;i++) {
                denoised[i][xyz] /= weights[xyz];
                err += (denoised[i][xyz]-images[i][xyz])*(denoised[i][xyz]-images[i][xyz]);
            }
            pcadim[xyz] /= weights[xyz];
            errmap[xyz] = (float)FastMath.sqrt(err/nimg);
        }
        images = null;
              
		return;
	}

}
