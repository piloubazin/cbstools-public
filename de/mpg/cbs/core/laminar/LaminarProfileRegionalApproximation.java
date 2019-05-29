package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class LaminarProfileRegionalApproximation {

	private float[] layersImage;
	//private float[] depthImage;
	private float[] intensityImage;
	private int[] maskImage=null;
	private String interpParam="linear";	
	
	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;
	
	private int    nbins = 20;
	private boolean optimize = true;

	private float[] weightImage;
	private float[] medianProfile;
	private float[] perc25Profile;
	private float[] perc75Profile;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	// create inputs
	public final void setProfileSurfaceImage(float[] val) { layersImage = val; }
	//public final void setDepthImage(float[] val) { depthImage = val; }
	public final void setIntensityImage(float[] val) { intensityImage = val; }
	public final void setRoiMask(int[] val) { maskImage = val; }
	public final void setInterpolation(String val) { interpParam = val; }
	
	public final void setDimensions(int x, int y, int z, int t) { nx=x; ny=y; nz=z; nt=t; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nt=dim[3]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Laminar Analysis"; }
	public final String getLabel() { return "Profile Sampling"; }
	public final String getName() { return "ProfileSampling"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Integrated Model-based Cognitive Neuroscience Reseaerch Unit, Universiteit van Amsterdam | Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Average some intensity cortical profiles for a given region of interest."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.2"; };
			
	// create outputs
	public final float[] getMedianProfile() { return medianProfile; }
	public final float[] getIqrProfile() { return iqrProfile; }
	public final float[] getProfileWeights() { return weightImage; }
	
	public void execute(){
		
		int nlayers = nt-1;
		
		float[][] layers = new float[nlayers+1][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<=nlayers;l++) {
			int xyz = x+nx*y+nx*ny*z;
			layers[l][xyz] = layersImage[xyz+nxyz*l];
		}
		layersImage = null;
		
		float[] intensity = intensityImage;
		
		// create a mask for all the regions outside of the area where layer 1 is > 0 and layer 2 is < 0
		boolean[] ctxmask = new boolean[nxyz];
        for (int xyz=0;xyz<nxyz;xyz++) {
            ctxmask[xyz] = (layers[0][xyz]>=0.0 && layers[nlayers][xyz]<=0.0);
        }
        // include only from the provided mask, if any
		if (maskImage!=null) {
			for (int xyz=0;xyz<nxyz;xyz++) {
				if (maskImage[xyz]==0) ctxmask[xyz] = false;
			}
		}
				
		// main algorithm
		CorticalProfile profile = new CorticalProfile(nlayers, nx, ny, nz, rx, ry, rz);
		
		byte LINEAR = 1;
		byte NEAREST = 2;
		byte interp = LINEAR;
		if (interpParam.equals("nearest")) interp = NEAREST;
		
		// get the profiles, remove any that has bad values
		float maskval = 1e13f;
		float[][] mapping = new float[nxyz][nlayers+1];
		boolean[] sampled = new boolean[nxyz];
		int nsample = 0;
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				profile.computeTrajectory(layers, x, y, z);
				sampled[xyz] = true;
				
				for (int l=0;l<=nlayers;l++) {
					// interpolate the contrast
					if (interp==NEAREST) {
						mapping[xyz][l] = ImageInterpolation.nearestNeighborInterpolation(intensity, ctxmask, maskval, 
																					profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																					nx, ny, nz);
					} else {
						mapping[xyz][l] = ImageInterpolation.linearInterpolation(intensity, ctxmask, maskval, 
																					profile.getPt(l)[X], profile.getPt(l)[Y], profile.getPt(l)[Z], 
																					nx, ny, nz);
					}
					if (mapping[xyz][l]==maskval) {
						sampled[xyz] = false;
					}
				}
				if (sampled[xyz]) nsample++;
			}
		}
		layers = null;
		intensity = null;
		ctxmask = null;
		
		// mapping to save time
		int[] index = new int[nsample];
		int id=0;
		for (int xyz=0;xyz<nxyz;xyz++) if (sampled[xyz]) {
		    index[id] = xyz;
		    id++;
		}
		
		// Chebyshev approximation
		float[][] mapcoefs = new float[nsample][degree+1];
		float[] mapres = new float[nsample];
		int[] mapdeg = new int[nsample];
		degree = Numerics.min(nlayers, 5);
		double[][] chebcoeff = new double[nlayers+1][degree+1];
		double[][] data = new double[nlayers+][1];
		for (int n=0;n<nsample;n++) {
		    for (int l=0;l<=layers;l++) {
                computeChebyshevCoefficients(chebcoeff[l], (double)l/nlayers, degree);
                data[l][0] = mapping[index[n]][l];
			}
            // invert the linear model
            Matrix mtx = new Matrix(chebcoeff);
            Matrix smp = new Matrix(data);
            Matrix val = mtx.solve(smp);
            for (int d=0;d<=degree;d++) {
                mapcoefs[n][d] = (float)val.get(d,0);
            }
            Matrix res = mtx.times(val).minus(smp);
            mapres[n] = (float)res.normInf();
            // optimise the degree?
            mapdeg[n] = (byte)degree;
				
            if (otpimize) {
                double mean = 0.0;
                double variance = 0.0;
                for (int l=0;l<nlayers+1;l++) mean += mapping[index[n]][l]/(nlayers+1.0);
                for (int l=0;l<nlayers+1;l++) variance += Numerics.square(mapping[index[n]][l]-mean)/nlayers;
                double evidence = 0.0;
                double previous = 0.0;
                for (int d=0;d<=degree;d++) {
                    previous = evidence;
					
                    for (int l=0;l<nlayers+1;l++) { 
                        double curve = 0.0;
                        for (int u=0;u<=d;u++) {
                            curve += mapcoefs[n][u]*chebcoeff[l][n];
                        }
                        evidence += Numerics.square(curve-mapping[index[n]][l])/variance;
                    }
                    evidence += (d+1);
                    if (d>0) {
                        // first minimum evidence: stop at the frst time it increases again
                        if (evidence>previous) {
                            mapdeg[n] = d-1;
                            d = degree+1;
						}
					}
				}
			}
		}
		
		// estimate the median, iqr over Chebyshev samples
		double[] median = new double[ndegree+1];
		double[] perc25 = new double[ndegree+1];
		double[] perc75 = new double[ndegree+1];
		for (int d=0;d<ndegree+1;d++) {
		    double[] val = new double[nsample];
		    for (int n=0;n<nsample;n++) val[n] = mapcoef[n][d];
		    Histogram dist = new Histogram(val, nbins, nsample);
		    median[d] = dist.percentage(0.50f);
		    perc75[d] = dist.percentage(0.75f);
		    perc25[d] = dist.percentage(0.25f);
		}

		// compute Gaussian weights
		double[] sampleWeights = new double[nsample];
		
		for (int n=0;n<nsample;n++) {
		    double dist = 0.0;
		    for (int d=0;d<ndegree+1;d++) dist += Numerics.square( (mapcoef[n][d]-median[d])/(1.349*(perc75[d]-perc25[d])) );
		    
		    sampleWeights[n] = FastMath.exp(-0.5*dist);
		}

		// output: weight map, median and iqr profile
		weightImage = new float[nxyz];
		for (int n=0;n<nsample;n++) {
		    weightImage[index[n]] = (float)sampleWeights[n];
		}
		resImage = new float[nxyz];
		for (int n=0;n<nsample;n++) {
		    resImage[index[n]] = (float)sampleWeights[n];
		}
		degImage = new int[nxyz];
		for (int n=0;n<nsample;n++) {
		    degImage[index[n]] = (float)sampleWeights[n];
		}
		medianProfile = new float[nlayers+1];
		iqrProfile = new float[nlayers+1];
		for (int l=0;l<nlayers+1;l++) {
		    medianProfile[l] = (float)median[l];
		    iqrProfile[l] = (float)iqr[l];
		}
		return;
	}

	
	private final void computeChebyshevCoefficients(double[] pol, float depth, int deg) {
		for (int n=0;n<=deg;n++) {
			pol[n] = 0.0;
		}
		double r = (2.0*(double)depth-1.0);
		
		if (deg>=0) pol[0] = 1.0;
		if (deg>=1) pol[1] = r;
		if (deg>=2) pol[2] = 2.0*r*r - 1.0;
		if (deg>=3) pol[3] = 4.0*r*r*r - 3.0*r;
		if (deg>=4) pol[4] = 8.0*r*r*r*r - 8.0*r*r + 1.0;
		if (deg>=5) pol[5] = 16.0*r*r*r*r*r - 20.0*r*r*r + 5.0*r;
		if (deg>=6) pol[6] = 32.0*r*r*r*r*r*r - 48.0*r*r*r*r + 18.0*r*r - 1.0;	
		if (deg==7) pol[7] = 64.0*r*r*r*r*r*r*r - 112.0*r*r*r*r*r + 56.0*r*r*r - 7.0*r;
		return;
	}

}
