package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class LaminarSmoothContrastMapping {

	private float[] layersImage;
	private float[] contrastImage;
	private float[] mappingImage = null;
	private int[] maskImage = null;
	
	private String measureParam="mean";
	public static final String[] measureTypes = {"mean", "median", "min", "max", "stdev"};
	private String interpParam="NN";	
	public static final String[] interpTypes = {"NN","linear"};
	private float distwmParam = 0.0f;
	private float distcsfParam = 0.0f;
	private float fwhmParam = 1.0f;
	private boolean smoothFirstParam = false;
	
	private int nx, ny, nz, nt, nxyz;
	private int nix, niy, niz, nit, nixyz;
	private float rx, ry, rz;
	private float rix, riy, riz;

	private float[] mappedContrastImage;
	private int[] mappedMaskImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final byte LINEAR = 1;
	private static final byte NEAREST = 2;
		
	private static final byte MEAN = 1;
	private static final byte STDEV = 2;
	private static final byte MIN = 3;
	private static final byte MAX = 4;
	private static final byte MEDIAN = 5;
		
	private static final float MASKVAL = 1e13f;
	
	// create inputs
	public final void setLayerSurfaceImage(float[] val) { layersImage = val; }
	public final void setIntensityImage(float[] val) { contrastImage = val; }
	public final void setIntensityMaskImage(int[] val) { maskImage = val; }
	public final void setContrastToLayersMapping(float[] val) { mappingImage = val; }
	
	public final void setMeasureAlongProfile(String val) { measureParam = val; }
	public final void setInterpolation(String val) { interpParam = val; }
	public final void setDistanceRatioToWM(float val) { distwmParam = val; }
	public final void setDistanceRatioToCSF(float val) { distcsfParam = val; }
	public final void setSmoothingFWHMmm(float val) { fwhmParam = val; }
	public final void setSmoothBeforeMeasure(boolean val) {smoothFirstParam = val; }
	
	public final void setLayersDimensions(int x, int y, int z, int t) { nx=x; ny=y; nz=z; nt=t; nxyz=nx*ny*nz; }
	public final void setLayersDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nt=dim[3]; nxyz=nx*ny*nz; }
	
	public final void setLayersResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setLayersResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final void setContrastDimensions(int x, int y, int z, int t) { nix=x; niy=y; niz=z; nit=t; nixyz=nix*niy*niz; }
	public final void setContrastDimensions(int[] dim) { nix=dim[0]; niy=dim[1]; niz=dim[2]; nit=dim[3]; nixyz=nix*niy*niz; }
	
	public final void setContrastResolutions(float x, float y, float z) { rix=x; riy=y; riz=z; }
	public final void setContrastResolutions(float[] res) { rix=res[0]; riy=res[1]; riz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Laminar Analysis"; }
	public final String getLabel() { return "Smooth Contrast Mapping"; }
	public final String getName() { return "SmoothContrastMapping"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for Neuroimaging, Netherlands Institute for Neuroscience, Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Sample contrast intensities summarized over cortical profiles and smoothed."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
			
	// create outputs
	public final float[] getCortexMappedImage() { return mappedContrastImage; }
	public final int[] getCortexMappedMask() { return mappedMaskImage; }
	
	public void execute(){
		
		int nlayers = nt-1;
		
		float[][] layers = new float[nlayers+1][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = layersImage[xyz+nxyz*l];
		}
		layersImage = null;
		
		// create a mask for all the regions outside of the area where WM ratio > param and csf ratio > param
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0);
			// measure relative depth? TODO
			if (ctxmask[xyz]) {
				float distwm = 0.0f;
				for (int n=0;n<nlayers;n++) {
					if (layers[n][xyz]>=0 && layers[n+1][xyz]<=0) distwm = (float)n/(float)nlayers + layers[n][xyz]/(layers[n][xyz]-layers[n+1][xyz])/(float)nlayers;
				}
				if (distwm<distwmParam) ctxmask[xyz] = false;
				if (1.0f-distwm<distcsfParam) ctxmask[xyz] = false;
			}
		}
		
		// create self-mapping if not set
		if (mappingImage==null) {
			mappingImage = new float[3*nxyz];
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				int xyz = x + nx*y + nx*ny*z;
				mappingImage[xyz+X*nxyz] = x;
				mappingImage[xyz+Y*nxyz] = y;
				mappingImage[xyz+Z*nxyz] = z;
			}
		}
		
		// create intensity mask
		boolean[] contrastMask = new boolean[nixyz];
		if (maskImage==null) {
			for (int xyz=0;xyz<nixyz;xyz++) {
				contrastMask[xyz] = false;
				for (int t=0;t<nit;t++) {
					if (contrastImage[xyz+t*nixyz]!=0) contrastMask[xyz] = true;
				}
			}
		} else {
			for (int xyz=0;xyz<nixyz;xyz++) {
				contrastMask[xyz] = (maskImage[xyz]>0);
			}
		}
		byte interp = NEAREST;
		if (interpParam.equals("linear")) interp = LINEAR;
					
		// smooth before sample, if selected
		boolean[] mappedcontrastmask = new boolean[nxyz];
		if (smoothFirstParam) {
			contrastImage = mapContrastToLayers(contrastImage, contrastMask, ctxmask, mappedcontrastmask, interp);
			if (fwhmParam>0) {
				directLaminarSmoothing(contrastImage, nit, layers, nlayers, mappedcontrastmask, ctxmask);
			}
		}
		
		// sample profile measure of choice
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		byte measure = MEAN;
		if (measureParam.equals("stdev")) measure = STDEV;
		else if (measureParam.equals("min")) measure = MIN;
		else if (measureParam.equals("max")) measure = MAX;
		else if (measureParam.equals("median")) measure = MEDIAN;
		
		float[] contrastMeasure = new float[nxyz*nit];
		float[][] samples = new float[nlayers+1][nit];
		boolean[] sampled = new boolean[nlayers+1];
		float maskval = MASKVAL;
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				profile.computeTrajectory(layers, x, y, z);
				
				for (int l=0;l<=nlayers;l++) {
					
					// interpolate the contrast
					if (smoothFirstParam) {
						// get the coordinates in layer space
						float xl = profile.getPt(l)[X];
						float yl = profile.getPt(l)[Y];
						float zl = profile.getPt(l)[Z];
						// data already in layer space
						if (interp==NEAREST) {
							for (int t=0;t<nit;t++) 
								samples[l][t] = ImageInterpolation.nearestNeighborInterpolation(contrastImage, t*nxyz, mappedcontrastmask, maskval, xl,yl,zl, nx,ny,nz);
						} else {
							for (int t=0;t<nit;t++) 
								samples[l][t] = ImageInterpolation.linearInterpolation(contrastImage, t*nxyz, mappedcontrastmask, maskval, xl,yl,zl, nx,ny,nz);
						}						
					} else {
						// find the coordinates in contrasts space
						float xi = ImageInterpolation.linearInterpolation(mappingImage, X*nxyz, ctxmask, maskval, 
																			profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																			nx, ny, nz);
						float yi = ImageInterpolation.linearInterpolation(mappingImage, Y*nxyz, ctxmask, maskval, 
																			profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																			nx, ny, nz);
						float zi = ImageInterpolation.linearInterpolation(mappingImage, Z*nxyz, ctxmask, maskval, 
																			profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																			nx, ny, nz);
						// data in contrast space
						if (interp==NEAREST) {
							for (int t=0;t<nit;t++) 
								samples[l][t] = ImageInterpolation.nearestNeighborInterpolation(contrastImage, t*nixyz, contrastMask, maskval, xi,yi,zi, nix,niy,niz);
						} else {
							for (int t=0;t<nit;t++) 
								samples[l][t] = ImageInterpolation.linearInterpolation(contrastImage, t*nixyz, contrastMask, maskval, xi,yi,zi, nix,niy,niz);
						}
					}
					sampled[l] = (samples[l][0]!=maskval);
				}
				
				// compute measure of interest
				if (measure==MEAN) {
					for (int t=0;t<nit;t++) {
						double sum = 0.0;
						double den = 0.0;
						for (int l=0;l<=nlayers;l++) if (sampled[l]) {
							sum += samples[l][t];
							den++;
						}
						if (den>0) {
						    contrastMeasure[xyz+t*nxyz] = (float)(sum/den);
						    mappedcontrastmask[xyz] = true;
						} else {
						    contrastMeasure[xyz+t*nxyz] = maskval;
						    mappedcontrastmask[xyz] = false;
						}
					}
				} else {
					System.out.println("not implemented yet");
				}
			}
		}
		
		// post smoothing
		if (!smoothFirstParam && fwhmParam>0) {
			directLaminarSmoothing(contrastMeasure, nit, layers, nlayers, mappedcontrastmask, ctxmask);
		}
		
		// output
		mappedContrastImage = contrastMeasure;
		mappedMaskImage = new int[nxyz];
		for (int xyz=0;xyz<nixyz;xyz++) {
			if (ctxmask[xyz] && mappedcontrastmask[xyz]) mappedMaskImage[xyz] = 1;
			else mappedMaskImage[xyz] = 0;
		}
	}

	private final void directLaminarSmoothing(float[] data, int nd, float[][] layers, int nlayers, boolean[] datamask, boolean[] ctxmask) {
				
		// get estimates for partial voluming of each layer
		System.out.println("Define partial volume coefficients");
		float[][] pvol = new float[nlayers+1][nxyz];
		for (int l=0;l<=nlayers;l++) {
			for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
				int xyz = x + nx*y + nx*ny*z;
				if (ctxmask[xyz]) {
					//pvol[l][xyz] = fastApproxPartialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
					pvol[l][xyz] =  Numerics.bounded(0.5f-layers[l][xyz],0.0f,1.0f);
					//else pvol[l][xyz] = partialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
				}	
			}
		}
		
		// iterations: number of voxels needed to reach the boundary, assuming that N times G_sigma0 approx. G_sigma1 if sigma1^2 = N sigma0^2 (ok if weight <<1)
		double sigma0 = 0.5f;
		double sigma1 = fwhmParam/Numerics.min(rx,ry,rz)/(2.0f*FastMath.sqrt(FastMath.log(4.0f)));
		int iterations = Numerics.ceil( (sigma1*sigma1)/(sigma0*sigma0) );
		// re-compute sigma0 to be exact
		sigma0 = sigma1/FastMath.sqrt(iterations);
		
		double weight = FastMath.exp(-1.0f/(2.0f*sigma0*sigma0));
		System.out.println("Standard deviation (voxels): "+sigma1);
		System.out.println("Number of iterations (sigma: "+sigma0+"): "+iterations);
		System.out.println("Neighbour weight = "+weight);
		
		// no black & white iterations: the layers are generally too thin
		float[] sdata = new float[nxyz*nd];
		boolean[] sdatamask = new boolean[nxyz];
		double[] layerval = new double[nd];
		double layersum = 0.0f;
	
		System.out.println("Smooth laminar data");
		for (int itr=0; itr<iterations; itr++) {
			//System.out.println("iteration "+(itr+1));
		
			// here we assume a linear combination across neighbors and layers
			// other options could be to keep only the layer with largest pv
			// note that the smoothing happens only parallel to the layers here
			for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
				double sumweight = 0.0f;
				for (int t=0;t<nd;t++) sdata[xyz+nxyz*t] = 0.0f;
				for (int l=0;l<nlayers;l++) {
					double pvweight = Numerics.max(pvol[l+1][xyz]-pvol[l][xyz],0.0f);
					if (pvweight>0) {
					    if (datamask[xyz]) {
                            for (int t=0;t<nd;t++) layerval[t] = pvweight*data[xyz+nxyz*t];
                            layersum = pvweight;
                        } else {
                            for (int t=0;t<nd;t++) layerval[t] = 0.0;
                            layersum = 0.0;
                        }
						for (byte k=0; k<26; k++) {
							int xyzn = Ngb.neighborIndex(k, xyz, nx, ny, nz);
							float dw = 1.0f/Ngb.neighborDistance(k);
							if (ctxmask[xyzn] && datamask[xyzn]) {
								double pv = Numerics.max(pvol[l+1][xyzn]-pvol[l][xyzn],0.0f);
								for (int t=0;t<nd;t++) layerval[t] += pv*weight*dw*data[xyzn+nxyz*t];
								layersum += pv*weight*dw;
							}
						}
						if (layersum>0) {
                            for (int t=0;t<nd;t++) layerval[t] /= layersum;
                        } else {
                            pvweight = 0.0;
                        }
					}
					// if no data, pvweight is zero
					for (int t=0;t<nd;t++) sdata[xyz+nxyz*t] += (float)(pvweight*layerval[t]);
					sumweight += pvweight;
				}
				if (sumweight>0) {
                    for (int t=0;t<nd;t++) sdata[xyz+nxyz*t] /= (float)sumweight;
                    sdatamask[xyz] = true;
                } else {
                    sdatamask[xyz] = false;
                }
			}
			// copy back to original data image
			for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
				for (int t=0;t<nd;t++) data[xyz+nxyz*t] = sdata[xyz+nxyz*t];
				// update mask: values at boundaries with actual data get filled in
				datamask[xyz] = sdatamask[xyz];
			}
		}
		
		return;
	}

	private final float[] mapContrastToLayers(float[] input, boolean[] inputmask, boolean[] ctxmask, boolean[] mappedmask, byte interp) {
				
		float[] data = new float[nxyz*nit];
		float maskval = MASKVAL;
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				// find the coordinates in contrasts space
				float xi = ImageInterpolation.linearInterpolation(mappingImage, X*nxyz, ctxmask, maskval, x,y,z, nx,ny,nz);
				float yi = ImageInterpolation.linearInterpolation(mappingImage, Y*nxyz, ctxmask, maskval, x,y,z, nx,ny,nz);
				float zi = ImageInterpolation.linearInterpolation(mappingImage, Z*nxyz, ctxmask, maskval, x,y,z, nx,ny,nz);
						
				// interpolate the contrast
				if (interp==NEAREST) {
					for (int t=0;t<nit;t++) 
						data[xyz+nxyz*t] = ImageInterpolation.nearestNeighborInterpolation(input, t*nixyz, inputmask, maskval, xi,yi,zi, nix,niy,niz);
				} else {
					for (int t=0;t<nit;t++) 
						data[xyz+nxyz*t] = ImageInterpolation.linearInterpolation(input, t*nixyz, inputmask, maskval, xi,yi,zi, nix,niy,niz);
				}
				// set the output mask
				mappedmask[xyz] = true;
				for (int t=0;t<nit;t++) 
					if (data[xyz+nxyz*t]==maskval)
						mappedmask[xyz] = false;
			}
		}
		return data;
	}

}
