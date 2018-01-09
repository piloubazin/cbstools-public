package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class LaminarProfileGeometry {

	// containers
	private float[] 	layersImage;

	// parameters
	private float		extraParam = 0.0f;
	private	float		smoothingParam = 0.3f;
	private String		calcParam;
	private String		regParam;
	
	public static final String[] calcTypes = {"thickness", 
						    "curvedness", 
						    "shape_index",
						    "mean_curvature", 
						    "gauss_curvature",
						    "profile_length", 
						    "profile_curvature", 
						    "profile_torsion",
						    "profile_direction"};
						    
	public static final String[] regTypes = {"none","Gaussian"};
	
	private float[] calcImage;
	private boolean isCalc4D = false;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;
	
	// create inputs
	public final void setLayerSurfaceImage(float[] val) { layersImage = val; }
	public final void setComputedMeasure(String val) { calcParam = val; }
	public final void setRegularization(String val) { regParam = val; }
	public final void setSmoothing(float val) { smoothingParam = val; }
	public final void setOutsideExtension_mm(float val) { extraParam = val; }
	
	public final void setLayersDimensions(int x, int y, int z, int t) { nx=x; ny=y; nz=z; nt=t; nxyz=nx*ny*nz; }
	public final void setLayersDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nt=dim[3]; nxyz=nx*ny*nz; }
	
	public final void setLayersResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setLayersResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Laminar Analysis"; }
	public final String getLabel() { return "Profile Geometry"; }
	public final String getName() { return "ProfileGeometry"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Spinoza Centre for Neuroimaging, Netherlands Institute for Neuroscience, Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Compute various geometric quantities for a cortical layers."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };

	// create outputs
	public final float[] getResultImage() { return calcImage; }
	public final boolean isResult4D() { return isCalc4D; }

	public final void execute() {
		
		int nlayers = nt-1;
		
		// note: we assume isotropic data
		float[][] layers = new float[nlayers+1][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = layersImage[xyz+nxyz*l];
		}
		layersImage = null;
		
		// create a mask for regions of value 0 on center image (kinda arbitrary..)
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layers[0][xyz] >= -extraParam/rx && layers[nlayers][xyz] <= extraParam/rx);
		}
		
		// main algorithm: thickness & curve parameters
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		float[] calc = null;
		if (calcParam.equals("profile_direction")) {
		    calc = new float[nxyz*3];
		    isCalc4D = true;
		} else {
		    calc = new float[nxyz];
		}
		// if surface processing: compute the smoothing
		/* in all cases
		if (calcParam.equals("mean_curvature") 
			|| calcParam.equals("gauss_curvature") 
			|| calcParam.equals("shape_index") 
			|| calcParam.equals("curvedness") ) {
		*/
		if (regParam.equals("Gaussian")) {
			float[][] kernel = ImageFilters.separableGaussianKernel(smoothingParam, smoothingParam, smoothingParam);
			int ks = (kernel[0].length-1)/2;
			
			for (int l=0;l<=nlayers;l++) {
				layers[l] = ImageFilters.separableConvolution(layers[l], nx, ny, nz, kernel, ks, ks, ks);
			}
		}
		//}
		
		for (int x=1; x<nx-1; x++) for (int y=1; y<ny-1; y++) for (int z=1; z<nz-1; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				if (calcParam.equals("thickness")) {
					// basic thickness definition
					calc[xyz] = (layers[0][xyz]-layers[nlayers][xyz])*rx;
				} else
				if (calcParam.equals("mean_curvature") 
					|| calcParam.equals("gauss_curvature") 
					|| calcParam.equals("shape_index") 
					|| calcParam.equals("curvedness")  
					|| calcParam.equals("profile_direction") ) {
				
					// search for closest surface to the point
					int lbest = 0;
					for (int l=1;l<=nlayers;l++) {
						if (Numerics.abs(layers[l][xyz])<Numerics.abs(layers[lbest][xyz]) ) {
							lbest = l;
						}
					}
					// compute the corresponding measure
					if (calcParam.equals("mean_curvature")) {
						calc[xyz] = ImageGeometry.meanCurvatureValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.equals("gauss_curvature")) {
						calc[xyz] = ImageGeometry.gaussCurvatureValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.equals("shape_index")) {
						calc[xyz] = ImageGeometry.shapeIndexValue(layers[lbest], xyz, nx, ny, nz);
					} else if (calcParam.equals("curvedness")) {
						calc[xyz] = ImageGeometry.curvednessValue(layers[lbest], xyz, nx, ny, nz);
					}else if (calcParam.equals("profile_direction")) {
						calc[xyz+X*nxyz] = 0.5f*(layers[lbest][xyz+1]-layers[lbest][xyz-1]);
						calc[xyz+Y*nxyz] = 0.5f*(layers[lbest][xyz+nx]-layers[lbest][xyz-nx]);
						calc[xyz+Z*nxyz] = 0.5f*(layers[lbest][xyz+nx*ny]-layers[lbest][xyz-nx*ny]);
						float norm = (float)FastMath.sqrt(calc[xyz+X*nxyz]*calc[xyz+X*nxyz]
						                                 +calc[xyz+Y*nxyz]*calc[xyz+Y*nxyz]
						                                 +calc[xyz+Z*nxyz]*calc[xyz+Z*nxyz]);
						if (norm>0) {
							calc[xyz+X*nxyz] /= norm; 
							calc[xyz+Y*nxyz] /= norm; 
							calc[xyz+Z*nxyz] /= norm; 
						}
					}
				} else {
					// first compute the profile	
					float err = profile.computeTrajectory(layers, x, y, z);
					if (calcParam.equals("profile_length")) {
						calc[xyz] = profile.computeLength();
					} else if (calcParam.equals("profile_curvature")) {
						calc[xyz] = profile.computeCurvature();
					} else if (calcParam.equals("profile_torsion")) {
						calc[xyz] = profile.computeTorsion();
					}
				}
			}
		}

		// output
		calcImage = calc;
				
		return;			
	}

}
