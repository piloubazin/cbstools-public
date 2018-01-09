package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.io.*;
import java.util.*;
import org.apache.commons.math3.util.FastMath;


public class LaminarIterativeSmoothing {
	private float[] intensityImage;
	private float[] layersImage;
	private int[] maskImage = null;
	
	private float    fwhmParam = 5.0f;
	private int      nlayers = 1;  
	
	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;

	private float[] smoothIntensityImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final float HASQ3 = (float)(FastMath.sqrt(3.0f)/2.0f);


	public final void setIntensityImage(float[] val) { intensityImage = val; }
	public final void setProfileSurfaceImage(float[] val) { layersImage = val; }
	public final void setROIMask(int[] val) { maskImage = val; }
	public final void setFWHMmm(float val) { fwhmParam = val; }
	
	public final void setLayers(int val) { nlayers=val; }
	public final void set4thDimension(int val) { nt=val; }
		
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Laminar Analysis.devel"; }
	public final String getLabel() { return "Laminar Iterative Smoothing"; }
	public final String getName() { return "LaminarIterativeSmoothing"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Smoothes data image along layer surfaces with an iterative heat kernel approach."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.1"; };
	
	// create outputs
    public final float[] getSmoothedIntensityImage() { return smoothIntensityImage; }
		
	public void execute() {
				
		// create a set of ROI masks for all the regions and also outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layersImage[xyz]>=0.0f && layersImage[xyz+nlayers*nxyz]<=0.0f);
		}
		
		if (maskImage!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) if (maskImage[xyz]==0) ctxmask[xyz] = false;
		}
		
		// remove image boundaries
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			if (x<=1 || x>=nx-2 || y<=1 || y>=ny-2 || z<=1 || z>=nz-2) {
				ctxmask[xyz] = false;
			}
		}
		
		// mask size
		int nctx=0;
		for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) nctx++;
		System.out.println("cortex mask size: "+nctx+" voxels");
		System.out.println("layers: "+nlayers);
		
		// get estimates for partial voluming of each layer
		System.out.println("Define partial volume coefficients");
		float[][] pvol = new float[nlayers+1][nxyz];
		for (int l=0;l<=nlayers;l++) {
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (ctxmask[xyz]) {
				//System.out.print(".");
					//pvol[l][xyz] = fastApproxPartialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
					pvol[l][xyz] =  Numerics.bounded(0.5f-layersImage[xyz+l*nxyz],0.0f,1.0f);
					//else pvol[l][xyz] = partialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
				}	}
		}
		
		// iterations: number of voxels needed to reach the boundary, assuming that N times G_sigma0 approx. G_sigma1 if sigma1^2 = N sigma0^2 (ok if weight <<1)
		double sigma0 = 0.5f;
		double sigma1 = fwhmParam/rx/(2.0f*FastMath.sqrt(FastMath.log(4.0f)));
		int iterations = Numerics.ceil( (sigma1*sigma1)/(sigma0*sigma0) );
		// re-compute sigma0 to be exact
		sigma0 = sigma1/FastMath.sqrt(iterations);
		
		double weight = FastMath.exp(-1.0f/(2.0f*sigma0*sigma0));
		System.out.println("Standard deviation (voxels): "+sigma1);
		System.out.println("Number of iterations (sigma: "+sigma0+"): "+iterations);
		System.out.println("Neighbour weight = "+weight);
		
		// no black & white iterations: the layers are generally too thin
		float[] sintensityImage = new float[nxyz*nt];
		double[] layerval = new double[nt];
		double layersum = 0.0f;
	
		for (int itr=0; itr<iterations; itr++) {
			//System.out.println("iteration "+(itr+1));
		
			// here we assume a linear combination across neighbors and layers
			// other options could be to keep only the layer with largest pv
			// note that the smoothing happens only parallel to the layers here
			for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
				double sumweight = 0.0f;
				for (int t=0;t<nt;t++) sintensityImage[xyz+nxyz*t] = 0.0f;
				for (int l=0;l<nlayers;l++) {
					double pvweight = Numerics.max(pvol[l+1][xyz]-pvol[l][xyz],0.0f);
					if (pvweight>0) {
						for (int t=0;t<nt;t++) layerval[t] = pvweight*intensityImage[xyz+nxyz*t];
						layersum = pvweight;
						for (byte k=0; k<26; k++) {
							int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
							float dw = 1.0f/Ngb.neighborDistance(k);
							if (ctxmask[xyzn]) {
								double pv = Numerics.max(pvol[l+1][xyzn]-pvol[l][xyzn],0.0f);
								for (int t=0;t<nt;t++) layerval[t] += pv*weight*dw*intensityImage[xyzn+nxyz*t];
								layersum += pv*weight*dw;
							}
						}
						for (int t=0;t<nt;t++) layerval[t] /= layersum;
					}
					for (int t=0;t<nt;t++) sintensityImage[xyz+nxyz*t] += (float)(pvweight*layerval[t]);
					sumweight += pvweight;
				}
				for (int t=0;t<nt;t++) sintensityImage[xyz+nxyz*t] /= (float)sumweight;
			}
			// copy back to original data image
			for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
				for (int t=0;t<nt;t++) intensityImage[xyz+nxyz*t] = sintensityImage[xyz+nxyz*t];
			}
		}				
		smoothIntensityImage = sintensityImage;
		
		System.out.println("Done");
		
		return;
	}
	/*
	float partialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int offset, int nx, int ny, int nz) {
		int xyz0 = x0+nx*y0+nx*ny*z0;
		double phi0 = levelset[xyz0+offset];
		
		if (phi0<-HASQ3) return 1.0f;
		else if (phi0>HASQ3) return 0.0f;
		
		// use direction, distance to create a continuous plane model V(X-X0)+phi(X0) = d(X,P)
		double vx = 0.5f*(levelset[xyz0+1+offset]-levelset[xyz0-1+offset]);
		double vy = 0.5f*(levelset[xyz0+nx+offset]-levelset[xyz0-nx+offset]);
		double vz = 0.5f*(levelset[xyz0+nx*ny+offset]-levelset[xyz0-nx*ny+offset]);
		double norm = FastMath.sqrt(vx*vx+vy*vy+vz*vz);
		
		// integrate the voxel cube over that model
		
		// first shot: numerically (other options: Heaviside integration, case-by-case geometric solution?)
		double volume = 0.0f;
		for (double x=-0.4f5; x<=0.4f5; x+=0.1f) for (double y=-0.4f5; y<=0.4f5; y+=0.1f) for (double z=-0.4f5; z<=0.4f5; z+=0.1f) {
			if (vx*x+vy*y+vz*z+norm*phi0<=0) volume+=0.001f;
		}
		
		// second option: with explicit definition of the intersection (~marching cubes!)
		// requires intersecting all 12 edges from the voxel cube, then defining either case-by-case formulas or 
		// deriving the formula for an integral of the Heaviside function
		
		return (float)volume;
	}

	float fastApproxPartialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int offset, int nx, int ny, int nz) {
		return Numerics.bounded(0.5f-levelset[x0+nx*y0+nx*ny*z0+offset],0.0f,1.0f);
	}
	*/
}
