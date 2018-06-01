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
public class IntensityComplexPCADenoising {

	// input parameters
	private		float[] 	invmag = null;
	private		float[] 	invphs = null;
	private		int			nx, ny, nz, nxyz;
	private 	float 		rx, ry, rz;

	private     int         nimg = 5;
	private		float		stdevCutoff = 1.1f;
	private     int         minDimension = 2;
	private     int         maxDimension = -1;
	private     int         ngbSize = 5;
	private     int         winSize = 5;
	//private     boolean     separate = false;
	//private     boolean     tvmag = false;
	//private     boolean     tvphs = false;
	
	//public static final String[]     thresholdingTypes = {"Eigenvalues","Global noise","Second half"};
	//private     String      thresholding = "Second half";
		
	// output parameters
	private		float[] globalpcadim = null;
	private		float[] globalerrmap = null;
	
	// internal variables (for debug)
	private float[][] images;
	private float[][] denoised;
	private float[][] eigvec;
	private float[][] eigval;
	private	float[][] pcadim;
	private	float[][] errmap;
	
	// set inputs
	public final void setMagnitudeImages(float[] in) { invmag = in; }
	public final void setPhaseImages(float[] in) { invphs = in; }

	public final void setMagnitudeAndPhaseImage(float[] in) {
	    if (invmag==null || invphs==null) {
	        if (nxyz==0) nxyz = in.length/2;
	        invmag = new float[nxyz*nimg];
	        invphs = new float[nxyz*nimg];
	    }
	    for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                invmag[xyz+i*nxyz] = in[xyz+i*nxyz];
                invphs[xyz+i*nxyz] = in[xyz+i*nxyz+nimg*nxyz];
            }
        }
	}
	
	public final void setImageNumber(int in) { nimg = in; }
	public final void setStdevCutoff(float in) { stdevCutoff = in; }
	public final void setMinimumDimension(int in) { minDimension = in; }
	public final void setMaximumDimension(int in) { maxDimension = in; }
	public final void setPatchSize(int in) { ngbSize = in; }
	public final void setWindowSize(int in) { winSize = in; }
	//public final void setSeparateProcessing(boolean in) { separate = in; }
	
	//public final void setMagnitudeTVSubtraction(boolean in) { tvmag = in; }
	//public final void setPhaseTVSubtraction(boolean in) { tvphs = in; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity.devel"; } 
	public final String getLabel() { return "Complex PCA denoising"; }
	public final String getName() { return "Complex PCA denoising"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for  Neuroimaging | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Denoise the data with a PCA-based method (adapted from Manjon et al., Plos One 2013."; }
		
	public final String getVersion() { return "3.1.3"; }

	// get outputs
	public float[] getDenoisedMagnitudeImages() { return invmag; }
	public float[] getDenoisedPhaseImages() { return invphs; }
	public float[] getDenoisedMagnitudeAndPhaseImage() {
	    float[] combi = new float[2*nxyz*nimg];
	    for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                combi[xyz+i*nxyz] = invmag[xyz+i*nxyz];
                combi[xyz+i*nxyz+nimg*nxyz] = invphs[xyz+i*nxyz];
            }
	    }
	    return combi;
	}
	
	public float[] getLocalDimensionImage() { return globalpcadim; }
	public float[] getNoiseFitImage() { return globalerrmap; }
	
	/*
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
	*/
	
	public void execute() {
	    if (invphs==null) executeMagnitudeDenoising();
	    else executeComplexDenoising();
	}
	
	private void executeComplexDenoising() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// we assume 4D images of size nimg
		if (invmag==null || invphs==null) System.out.print("data stacks not properly initialized!\n");
		
		// renormalize phase
		float phsmin = invphs[0];
		float phsmax = invphs[0];
        for (int xyz=0;xyz<nimg*nxyz;xyz++) {
			if (invphs[xyz]<phsmin) phsmin = invphs[xyz];
			if (invphs[xyz]>phsmax) phsmax = invphs[xyz];
		}
		double phsscale = (phsmax-phsmin)/(2.0*FastMath.PI);
		
		// unwrap phase and remove TV global variations
        float[][] tvimgphs = null;
		//if (tvphs) {
        tvimgphs = new float[nimg][];
        float[] phs = new float[nxyz];
        for (int i=0;i<nimg;i++) {
            System.out.print("global variations removal phase "+(i+1)+"\n");
            for (int xyz=0;xyz<nxyz;xyz++) phs[xyz] = invphs[xyz+i*nxyz];
             // unwrap phase images
            IntensityFastMarchingUnwrapping unwrap = new IntensityFastMarchingUnwrapping();
            unwrap.setPhaseImage(phs);
            unwrap.setDimensions(nx,ny,nz);
            unwrap.setResolutions(rx,ry,rz);
            unwrap.setTVScale(0.33f);
            unwrap.setTVPostProcessing("TV-approximation");
            unwrap.execute();
            tvimgphs[i] = unwrap.getCorrectedImage();
            for (int xyz=0;xyz<nxyz;xyz++) invphs[xyz+i*nxyz] -= tvimgphs[i][xyz];
        }
        //}
		
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
		
		// 2. estimate PCA in slabs of NxNxN size xT windows
		int ngb = ngbSize;
		int nstep = Numerics.floor(ngb/2.0);
		int nimg2 = 2*nimg;
		int ntime = Numerics.min(winSize,nimg);
		if (ntime<0) ntime = nimg;
		int ntime2 = 2*ntime;
		int tstep = Numerics.floor(ntime/2.0);
		int nsample = Numerics.ceil(nimg/tstep);
		denoised = new float[nimg2][nxyz];
		//eigvec = new float[nimg2][nxyz];
		//eigval = new float[nimg2][nxyz];
		float[][] weights = new float[nimg][nxyz];
		pcadim = new float[nimg][nxyz];
		errmap = new float[nimg][nxyz];
		// border issues should be cleaned-up, ignored so far
        for (int t=0;t<nimg;t+=tstep) {
            boolean last = false;
		    boolean skip = false;
		    System.out.print("step "+(t/tstep));
    		if (t+ntime>nimg) {
		        if (ntime<nimg) {
                    // shift windows at the end of the time domain if needed
                    t = nimg-ntime;
                    System.out.print(":");
                    last = true;
                } else {
                    // skip this round entirely
                    System.out.print("x\n");
                    skip = true;
                }
            }
            if (!skip) {
                System.out.print("...\n");
                for (int x=0;x<nx;x+=nstep) for (int y=0;y<ny;y+=nstep) for (int z=0;z<nz;z+=nstep) {
                    int ngbx = Numerics.min(ngb, nx-x);
                    int ngby = Numerics.min(ngb, ny-y);
                    int ngbz = Numerics.min(ngb, nz-z);
                    int ngb3 = ngbx*ngby*ngbz;
                    boolean process = false;
                    if (ngb3<ntime2) {
                        System.out.print("!patch is too small!\n");
                        process = false;
                    } else {
                        process = true;
                    }
                    if (process) {
                        double[][] patch = new double[ngb3][ntime2];
                        for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) for (int i=2*t;i<2*t+ntime2;i++) {
                            patch[dx+ngbx*dy+ngbx*ngby*dz][i-2*t] = images[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)];
                        }
                        // mean over samples
                        double[] mean = new double[ntime2];
                        for (int i=0;i<ntime2;i++) {
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
                        for (int n=0;n<ngb3;n++) for (int i=0;i<ntime2;i++) {
                            sigma += patch[n][i]*patch[n][i];
                        }
                        sigma /= ntime2*ngb3;
                        sigma = FastMath.sqrt(sigma);
                        
                        // cutoff
                        //System.out.println("eigenvalues: ");
                        double[] eig = new double[ntime2];
                        boolean[] used = new boolean[ntime2];
                        int nzero=0;
                        double eigsum = 0.0;
                        for (int n=0;n<ntime2;n++) {
                            eig[n] = svd.getSingularValues()[n];
                            eigsum += Numerics.abs(eig[n]);
                        }
                        // fit second half linear decay model
                        double[] loc = new double[ntime];
                        double[][] fit = new double[ntime][1];
                        for (int n=ntime;n<ntime2;n++) {
                            loc[n-ntime] = (n-ntime)/(double)ntime;
                            fit[n-ntime][0] = Numerics.abs(eig[n]);
                        }
                        double[][] poly = new double[ntime][2];
                        for (int n=0;n<ntime;n++) {
                            poly[n][0] = 1.0;
                            poly[n][1] = loc[n];
                        }
                        // invert the linear model
                        Matrix mtx = new Matrix(poly);
                        Matrix smp = new Matrix(fit);
                        Matrix val = mtx.solve(smp);
                
                        // compute the expected value:
                        double[] expected = new double[ntime2];
                        for (int n=0;n<ntime2;n++) {
                            double n0 = (n-ntime)/(double)ntime;
                            // linear coeffs,
                            expected[n] = (val.get(0,0) + n0*val.get(1,0));
                            //expected[n] = n*slope + intercept;
                        }
                        double residual = 0.0;
                        double meaneig = 0.0;
                        double variance = 0.0;
                        for (int n=ntime;n<ntime;n++) meaneig += Numerics.abs(eig[n]);
                        meaneig /= ntime;
                        for (int n=ntime;n<ntime2;n++) {
                            variance += (meaneig-Numerics.abs(eig[n]))*(meaneig-Numerics.abs(eig[n]));
                            residual += (expected[n]-Numerics.abs(eig[n]))*(expected[n]-Numerics.abs(eig[n]));
                        }
                        double rsquare = 1.0;
                        if (variance>0) rsquare = Numerics.max(1.0 - (residual/variance), 0.0);
                        
                        for (int n=0;n<ntime2;n++) {
                            //System.out.print(" "+(eig[n]/sigma));
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
                        for (int n=0;n<ngb3;n++) for (int i=0;i<ntime2;i++) {
                            // Sum_s>0 s_kU_kV_kt
                            patch[n][i] = mean[i];
                            for (int j=0;j<ntime2;j++) if (used[j]) {
                                patch[n][i] += U.get(n,j)*eig[j]*V.get(i,j);
                            }
                        }
                        // add to the denoised image
                        double wpatch = (1.0/(1.0 + ntime2 - nzero));
                        for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) {
                            for (int i=0;i<ntime2;i++) {
                                denoised[2*t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*patch[dx+ngbx*dy+ngbx*ngby*dz][i]);
                                //eigval[2*t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*eig[i]);
                                //eigvec[2*t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*U.get(dx+ngbx*dy+ngbx*ngby*dz,i));
                            }
                            for (int i=0;i<ntime;i++) {
                                weights[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                                pcadim[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(ntime2-nzero));
                                errmap[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*rsquare);
                            }
                        }
                    }
                }
                if (last) t=nimg;
            }
        }
        for (int xyz=0;xyz<nxyz;xyz++) {
            double err = 0.0;
            for (int i=0;i<nimg2;i++) {
                int t = Numerics.floor(i/2.0f);
                denoised[i][xyz] /= weights[t][xyz];
                //eigval[i][xyz] /= weights[t][xyz];
                //eigvec[i][xyz] /= weights[t][xyz];
            }
            for (int i=0;i<nimg;i++) {
                pcadim[i][xyz] /= weights[i][xyz];
                errmap[i][xyz] /= weights[i][xyz];
            }
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
        globalpcadim = new float[nimg*nxyz];
        globalerrmap = new float[nimg*nxyz];
  		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                globalpcadim[xyz+i*nxyz] = pcadim[i][xyz];
                globalerrmap[xyz+i*nxyz] = errmap[i][xyz];
            }
        }
        
        // opt. add back the TV estimate, if needed
        //if (tvphs) {
        for (int i=0;i<nimg;i++) for (int xyz=0;xyz<nxyz;xyz++) {
            invphs[xyz+i*nxyz] += tvimgphs[i][xyz];
            // wrap around phase values?
            invphs[xyz+i*nxyz] = Numerics.modulo(invphs[xyz+i*nxyz], phsmax-phsmin);
        }
        //}
		return;
	}

	private void executeMagnitudeDenoising() {
		// this assumes all the inputs are already set
		
		// main algorithm
		
		// we assume 4D images of size nimg
		if (invmag==null) System.out.print("data stacks not properly initialized!\n");
		
		// 1. pass directly the magnitude signal
		images = new float[nimg][nxyz];
		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                images[i][xyz] = invmag[xyz+i*nxyz];
            }
        }
        invmag = null;
        
		// 2. estimate PCA in slabs of NxNxN size xT windows
		int ngb = ngbSize;
		int nstep = Numerics.floor(ngb/2.0);
		int ntime = Numerics.min(winSize,nimg);
		if (ntime<0) ntime = nimg;
		int tstep = Numerics.floor(ntime/2.0);
		int nsample = Numerics.ceil(nimg/tstep);
        System.out.print("patch dimensions ["+ngb+" x "+ntime+"] shifting by ["+nstep+" x "+tstep+"]\n");
		System.out.print("time steps: "+nsample+" (over "+nimg+" time points)\n");
		 
		denoised = new float[nimg][nxyz];
		//eigvec = new float[nimg][nxyz];
		//eigval = new float[nimg][nxyz];
		float[][] weights = new float[nimg][nxyz];
		pcadim = new float[nimg][nxyz];
		errmap = new float[nimg][nxyz];
		// border issues should be cleaned-up, ignored so far
		for (int t=0;t<nimg;t+=tstep) {
		    boolean last = false;
		    boolean skip = false;
		    System.out.print("step "+(t/tstep));
    		if (t+ntime>nimg) {
		        if (ntime<nimg) {
                    // shift windows at the end of the time domain if needed
                    t = nimg-ntime;
                    System.out.print(":");
                    last = true;
                } else {
                    // skip this round entirely
                    System.out.print("x\n");
                    skip = true;
                }
            }
            if (!skip) {
                System.out.print("...\n");
                for (int x=0;x<nx;x+=nstep) for (int y=0;y<ny;y+=nstep) for (int z=0;z<nz;z+=nstep) {
                    int ngbx = Numerics.min(ngb, nx-x);
                    int ngby = Numerics.min(ngb, ny-y);
                    int ngbz = Numerics.min(ngb, nz-z);
                    int ngb3 = ngbx*ngby*ngbz;
                    boolean process = false;
                    if (ngb3<ntime) {
                        System.out.print("!patch is too small!\n");
                        process = false;
                    } else {
                        process = true;
                    }
                    if (process) {
                        double[][] patch = new double[ngb3][ntime];
                        for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) for (int i=t;i<t+ntime;i++) {
                            patch[dx+ngbx*dy+ngbx*ngby*dz][i-t] = images[i][x+dx+nx*(y+dy)+nx*ny*(z+dz)];
                        }
                        // mean over samples
                        double[] mean = new double[ntime];
                        for (int i=0;i<ntime;i++) {
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
                        for (int n=0;n<ngb3;n++) for (int i=0;i<ntime;i++) {
                            sigma += patch[n][i]*patch[n][i];
                        }
                        sigma /= ntime*ngb3;
                        sigma = FastMath.sqrt(sigma);
                        
                        // cutoff
                        //System.out.println("eigenvalues: ");
                        double[] eig = new double[ntime];
                        boolean[] used = new boolean[ntime];
                        int nzero=0;
                        double eigsum = 0.0;
                        for (int n=0;n<ntime;n++) {
                            eig[n] = svd.getSingularValues()[n];
                            eigsum += Numerics.abs(eig[n]);
                        }
                        // fit second half linear decay model
                        int ntimeh = Numerics.ceil(ntime/2.0f);
                        double[] loc = new double[ntimeh];
                        double[][] fit = new double[ntimeh][1];
                        for (int n=ntime-ntimeh;n<ntime;n++) {
                            loc[n-ntime+ntimeh] = (n-ntime+ntimeh)/(double)ntimeh;
                            fit[n-ntime+ntimeh][0] = Numerics.abs(eig[n]);
                        }
                        double[][] poly = new double[ntimeh][2];
                        for (int n=0;n<ntimeh;n++) {
                            poly[n][0] = 1.0;
                            poly[n][1] = loc[n];
                        }
                        // invert the linear model
                        Matrix mtx = new Matrix(poly);
                        Matrix smp = new Matrix(fit);
                        Matrix val = mtx.solve(smp);
                
                        // compute the expected value:
                        double[] expected = new double[ntime];
                        for (int n=0;n<ntime;n++) {
                            double n0 = (n-ntime+ntimeh)/(double)ntimeh;
                            // linear coeffs,
                            expected[n] = (val.get(0,0) + n0*val.get(1,0));
                            //expected[n] = n*slope + intercept;
                        }
                        double residual = 0.0;
                        double meaneig = 0.0;
                        double variance = 0.0;
                        for (int n=ntimeh;n<ntime;n++) meaneig += Numerics.abs(eig[n]);
                        meaneig /= (ntime-ntimeh);
                        for (int n=ntimeh;n<ntime;n++) {
                            variance += (meaneig-Numerics.abs(eig[n]))*(meaneig-Numerics.abs(eig[n]));
                            residual += (expected[n]-Numerics.abs(eig[n]))*(expected[n]-Numerics.abs(eig[n]));
                        }
                        double rsquare = 1.0;
                        if (variance>0) rsquare = Numerics.max(1.0 - (residual/variance), 0.0);
                        
                        for (int n=0;n<ntime;n++) {
                            //System.out.print(" "+(eig[n]/sigma));
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
                        for (int n=0;n<ngb3;n++) for (int i=0;i<ntime;i++) {
                            // Sum_s>0 s_kU_kV_kt
                            patch[n][i] = mean[i];
                            for (int j=0;j<ntime;j++) if (used[j]) {
                                patch[n][i] += U.get(n,j)*eig[j]*V.get(i,j);
                            }
                        }
                        // add to the denoised image
                        double wpatch = (1.0/(1.0 + ntime - nzero));
                        for (int dx=0;dx<ngbx;dx++) for (int dy=0;dy<ngby;dy++) for (int dz=0;dz<ngbz;dz++) {
                            for (int i=0;i<ntime;i++) {
                                denoised[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*patch[dx+ngbx*dy+ngbx*ngby*dz][i]);
                                //eigval[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*eig[i]);
                                //eigvec[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*U.get(dx+ngbx*dy+ngbx*ngby*dz,i));
                                weights[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                                pcadim[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(ntime-nzero));
                                errmap[t+i][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*rsquare);
                            }
                            //weights[(t+i)/tstep][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)wpatch;
                            //pcadim[(t+i)/tstep][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*(ntime-nzero));
                            //errmap[(t+i)/tstep][x+dx+nx*(y+dy)+nx*ny*(z+dz)] += (float)(wpatch*rsquare);
                        }
                    }
                }
                if (last) t=nimg;
            }
        }
        for (int xyz=0;xyz<nxyz;xyz++) {
            double err = 0.0;
            for (int i=0;i<nimg;i++) {
                denoised[i][xyz] /= weights[i][xyz];
                //eigval[i][xyz] /= weights[i][xyz];
                //eigvec[i][xyz] /= weights[i][xyz];
                pcadim[i][xyz] /= weights[i][xyz];
                errmap[i][xyz] /= weights[i][xyz];
            }
        }
        //images = null;
          
        // 3. rebuild magnitude and phase images
        invmag = new float[nimg*nxyz];
        for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                invmag[xyz+i*nxyz] = (float)denoised[i][xyz];
            }
        }
        globalpcadim = new float[nimg*nxyz];
        globalerrmap = new float[nimg*nxyz];
  		for (int i=0;i<nimg;i++) {
            for (int xyz=0;xyz<nxyz;xyz++) {
                globalpcadim[xyz+i*nxyz] = pcadim[i][xyz];
                globalerrmap[xyz+i*nxyz] = errmap[i][xyz];
            }
        }
  		return;
	}

}
