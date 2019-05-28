package de.mpg.cbs.core.laminar;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class LaminarProfileAveraging {

	private float[] layersImage;
	private float[] intensityImage;
	private int[] maskImage=null;
	private String interpParam="linear";	
	
	private int nx, ny, nz, nt, nxyz;
	private float rx, ry, rz;
	
	private int    nbins = 20;
	private float  mincorr = 0.0f;

	private float[] weightImage;
	private float[] sampleProfile;
	private float[] medianProfile;
	private float[] iqrProfile;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	// create inputs
	public final void setProfileSurfaceImage(float[] val) { layersImage = val; }
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
	public final float[] getSampleProfile() { return sampleProfile; }
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
		
		// estimate the median, iqr over samples
		double[] median = new double[nlayers+1];
		double[] iqr = new double[nlayers+1];
		for (int l=0;l<nlayers+1;l++) {
		    double[] val = new double[nsample];
		    for (int n=0;n<nsample;n++) val[n] = mapping[index[n]][l];
		    Histogram dist = new Histogram(val, nbins, nsample);
		    median[l] = dist.percentage(0.5f);
		    iqr[l] = dist.percentage(0.75f) - dist.percentage(0.25f);
		}

		// compute correlation weights
		double[] sampleWeights = new double[nsample];
		
		// here we use weighted averages over the profile length
		double[] profileWeights = new double[nlayers+1];
		double profileWeightSum = 0.0;
		for (int l=0;l<nlayers+1;l++) {
		    // constant values are just discarded
		    if (iqr[l]>0) profileWeights[l] = 1.0/iqr[l];   
		    else profileWeights[l] = 0.0;
		    profileWeightSum += profileWeights[l];
		}
		for (int l=0;l<nlayers+1;l++) profileWeights[l] /= profileWeightSum;
		
		double avgmed = 0.0;
		for (int l=0;l<nlayers+1;l++) avgmed += profileWeights[l]*median[l];
		double varmed = 0.0;
		for (int l=0;l<nlayers+1;l++) varmed += profileWeights[l]*(median[l]-avgmed)*(median[l]-avgmed);
		
		// correlation to find the best profile
		double sampleWeightMin = mincorr;
		double sampleWeightMax = 0.0;
		int nbest = -1;
		for (int n=0;n<nsample;n++) {
		    // use the profile average, or the median one (so value differences are important)
		    sampleWeights[n] = 0.0;
		    for (int l=0;l<nlayers+1;l++) sampleWeights[n] += profileWeights[l]*(mapping[index[n]][l]-avgmed)*(median[l]-avgmed)/varmed;
		    if (sampleWeights[n]>sampleWeightMax) {
		        sampleWeightMax = sampleWeights[n];
		        nbest = n;
		    }
		}
		double sampleWeightAvg = 0.0;
		for (int n=0;n<nsample;n++) {
		    sampleWeights[n] = Numerics.max(sampleWeights[n]-sampleWeightMin, 0.0f);
		    sampleWeightAvg += sampleWeights[n]/nsample;
		}
		System.out.println(" avg. correlation score: "+sampleWeightAvg);
		/*
		// correlate and recompute a weighted median, iqr, loop
		for (int t=0;t<10;t++) {
		    System.out.print("\n iteration "+t);
		    double best = sampleWeightAvg;
		    
		    // recompute median
		    for (int l=0;l<nlayers+1;l++) {
                double[] val = new double[nsample];
                for (int n=0;n<nsample;n++) val[n] = mapping[index[n]][l];
                Histogram dist = new Histogram(val, sampleWeights, nbins, nsample);
                median[l] = dist.percentage(0.5f);
                iqr[l] = dist.percentage(0.75f) - dist.percentage(0.25f);
            }
		    // recompute metric weight
		    profileWeightSum = 0.0;
            for (int l=0;l<nlayers+1;l++) {
                // constant values are just discarded
                if (iqr[l]>0) profileWeights[l] = 1.0/iqr[l];   
                else profileWeights[l] = 0.0;
                profileWeightSum += profileWeights[l];
            }
            for (int l=0;l<nlayers+1;l++) profileWeights[l] /= profileWeightSum;
		
            avgmed = 0.0;
            for (int l=0;l<nlayers+1;l++) avgmed += profileWeights[l]*median[l];

            // correlation to find the best profile
            sampleWeightMin = 0.0;
            sampleWeightMax = 0.0;
            sampleWeightAvg = 0.0;
            nbest = -1;
            for (int n=0;n<nsample;n++) {
                //double avgmap = 0.0;
                //for (int l=0;l<nlayers+1;l++) avgmap += profileWeights[l]*mapping[index[n]][l];
                
                sampleWeights[n] = 0.0;
                for (int l=0;l<nlayers+1;l++) sampleWeights[n] += profileWeights[l]*(mapping[index[n]][l]-avgmed)*(median[l]-avgmed);
                sampleWeights[n] = Numerics.max(sampleWeights[n],0.0);
                
                sampleWeightAvg += sampleWeights[n]/nsample;
                //if (sampleWeights[n]<sampleWeightMin) sampleWeightMin = sampleWeights[n];
                if (sampleWeights[n]>sampleWeightMax) {
                    sampleWeightMax = sampleWeights[n];
                    nbest = n;
                }
            }
            for (int n=0;n<nsample;n++) sampleWeights[n] = (sampleWeights[n]-sampleWeightMin)/(sampleWeightMax-sampleWeightMin);
            sampleWeightAvg = (sampleWeightAvg-sampleWeightMin)/(sampleWeightMax-sampleWeightMin);
            
            System.out.print(" avg. correlation score: "+sampleWeightAvg);
            if (sampleWeightAvg>best) {
                if (t>2 && sampleWeightAvg-best<1e-6) t = 1000;
                best = sampleWeightAvg;
            }
		}
		*/
		// random sampling
		int best = -1;
		double bestscore = 0.0;
		for (int t=0;t<1000;t++) {
		    int idx = Numerics.round(FastMath.random()*(nsample-1));
		    if (sampleWeights[idx]>0) {
                //double avg = 0.0;
                //for (int l=0;l<nlayers+1;l++) avg += profileWeights[l]*mapping[index[idx]][l];
                //double var = 0.0;
                //for (int l=0;l<nlayers+1;l++) var += profileWeights[l]*(mapping[index[idx]][l]-avg)*(mapping[index[idx]][l]-avg);
		        double score = 0.0;
		        for (int n=0;n<nsample;n++) {
		            double corr = 0.0;
		            for (int l=0;l<nlayers+1;l++) {
		                corr += profileWeights[l]*(mapping[index[n]][l]-avgmed)*(mapping[index[idx]][l]-avgmed)/varmed;
		            }
		            score += Numerics.max(corr-mincorr, 0.0)/nsample;
		        }
		        System.out.println(" sample correlation score: "+score);
		        if (score>bestscore) {
		            best = idx;
		            bestscore = score;
		        }
		    }
		}
		// recalculateeverything for the bestone
		//double avg = 0.0;
        //for (int l=0;l<nlayers+1;l++) avg += profileWeights[l]*mapping[index[best]][l];
        //double var = 0.0;
        //for (int l=0;l<nlayers+1;l++) var += profileWeights[l]*(mapping[index[best]][l]-avg)*(mapping[index[best]][l]-avg);
		System.out.println(" best sample correlation score: "+bestscore);
        double score = 0.0;
        for (int n=0;n<nsample;n++) {
            double corr = 0.0;
            for (int l=0;l<nlayers+1;l++) {
                corr += profileWeights[l]*(mapping[index[n]][l]-avgmed)*(mapping[index[best]][l]-avgmed)/varmed;
            }
            sampleWeights[n] = Numerics.max(corr-mincorr, 0.0);
        }
        
        // recompute median
        for (int l=0;l<nlayers+1;l++) {
            double[] val = new double[nsample];
            for (int n=0;n<nsample;n++) val[n] = mapping[index[n]][l];
            Histogram dist = new Histogram(val, sampleWeights, nbins, nsample);
            median[l] = dist.percentage(0.5f);
            iqr[l] = dist.percentage(0.75f) - dist.percentage(0.25f);
        }

        profileWeightSum = 0.0;
		for (int l=0;l<nlayers+1;l++) {
		    // constant values are just discarded
		    if (iqr[l]>0) profileWeights[l] = 1.0/iqr[l];   
		    else profileWeights[l] = 0.0;
		    profileWeightSum += profileWeights[l];
		}
		for (int l=0;l<nlayers+1;l++) profileWeights[l] /= profileWeightSum;
		
		avgmed = 0.0;
		for (int l=0;l<nlayers+1;l++) avgmed += profileWeights[l]*median[l];
		varmed = 0.0;
		for (int l=0;l<nlayers+1;l++) varmed += profileWeights[l]*(median[l]-avgmed)*(median[l]-avgmed);

        score = 0.0;
        for (int n=0;n<nsample;n++) {
            double corr = 0.0;
            for (int l=0;l<nlayers+1;l++) {
                corr += profileWeights[l]*(mapping[index[n]][l]-avgmed)*(median[l]-avgmed)/varmed;
            }
            score += Numerics.max(corr-mincorr, 0.0)/nsample;
        }
        System.out.println(" best median correlation score: "+score);

        // output: weight map, median and iqr profile
		weightImage = new float[nxyz];
		for (int n=0;n<nsample;n++) {
		    weightImage[index[n]] = (float)sampleWeights[n];
		}
		sampleProfile = new float[nlayers+1];
		medianProfile = new float[nlayers+1];
		iqrProfile = new float[nlayers+1];
		for (int l=0;l<nlayers+1;l++) {
		    sampleProfile[l] = (float)mapping[index[best]][l];
		    medianProfile[l] = (float)median[l];
		    iqrProfile[l] = (float)iqr[l];
		}
		return;
	}

}
