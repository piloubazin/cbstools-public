package de.mpg.cbs.core.segmentation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationSureSeg {

	private float[] input1Image = null;
	private float[] input2Image = null;
	private float[] input3Image = null;
	
	private float[] var1Image = null;
	private float[] var2Image = null;
	private float[] var3Image = null;
	
	private float scale1Param =1.0f;
	private float scale2Param =1.0f;
	private float scale3Param =1.0f;	
	
	private byte[] 	maskImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private float[] probaImage = null;
	private int[] 	labelImage = null;
	private	boolean	maxProba = false;
	private byte 	nbestParam	=	4;
	private boolean	includeBgParam = false;
	private boolean	rescaleProbaParam = false;
	
	private int 	iterationParam	=	100;
	private float 	imgscaleParam		=	1.0f;
	private	float 	certainscaleParam		= 	1.0f;
	private float 	mincertaintyParam		=	1.0f;
	private int 	neighborParam	=	4;
	private float 	diffratioParam		=	0.01f;
	
	private boolean	computeDistributionParam = false;
	private boolean	computeNoiseParam = false;
	
	// outputs
	private int[] segmentImage;
	private float[] maxprobaImage;
	private byte[] maxidImage;
	private float[] debugImage;
	
	// create inputs
	public final void setContrastImage1(float[] val) { input1Image = val; }
	public final void setContrastImage2(float[] val) { input2Image = val; }
	public final void setContrastImage3(float[] val) { input3Image = val; }
	
	public final void setNoiseImage1(float[] val) { var1Image = val; }
	public final void setNoiseImage2(float[] val) { var2Image = val; }
	public final void setNoiseImage3(float[] val) { var3Image = val; }
	
	public final void setContrastScale1(float val) { scale1Param = val; }
	public final void setContrastScale2(float val) { scale2Param = val; }
	public final void setContrastScale3(float val) { scale3Param = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setImageMask(byte[] val) { maskImage = val; }
	
	public final void setLabel4DProbabilities(float[] val) { probaImage = val; maxProba = false; }
	public final void setLabelMaxProbability(float[] val) { probaImage = val; maxProba = true; }
	public final void setLabelSegmentation(int[] val) { labelImage = val; }
	public final void setLabelDepth(byte val) { nbestParam = val; }
	public final void setBackgroundIncluded(boolean val) { includeBgParam = val; }
	public final void setRescaleIndividualProbabilities(boolean val) { rescaleProbaParam = val; }

	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setImageScale(float val) { imgscaleParam = val; }
	public final void setCertaintyScale(float val) { certainscaleParam = val; }
	public final void setMinCertainty(float val) { mincertaintyParam = val; }
	public final void setNeighborhoodSize(int val) { neighborParam = val; }
	public final void setMinDifferenceRatio(float val) { diffratioParam = val; }
	
	public final void setReestimateIntensityDistributions(boolean val) { computeDistributionParam = val; }
	public final void setEstimateNoise(boolean val) { computeNoiseParam = val; }

	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Segmentation.devel"; }
	public static final String getLabel() { return "Spatial Uncertainty Relaxation (SUR)"; }
	public static final String getName() { return "SUR"; }

	public static final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public static final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public static final String getDescription() { return "Sharpen segmentations by propagating information from high certainty to high uncertainty regions."; }
	public static final String getLongDescription() { return getDescription(); }
		
	public static final String getVersion() { return "3.1.0"; };

	// create outputs
	public final int[] getSegmentationImage() { return segmentImage; }
	public final float[] getMaxProbabilityImage() { return maxprobaImage; }
	public final byte[] getSegmentedIdsImage() { return maxidImage; }
	public final float[] getDebugImage() { return debugImage; }
	
	public void execute(){
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image != null) nimg++;
		if (input3Image != null) nimg++;
		
		float[][] image = new float[nimg][nxyz];
		float[] scaling = new float[nimg];
		float[][] noise = null;
		if (var1Image!=null) noise = new float[nimg][nxyz];
		int n = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[n][xyz] = input1Image[xyz];
			if (var1Image!=null) noise[n][xyz] = var1Image[xyz];
		}
		scaling[n] = scale1Param;
		input1Image = null;
		var1Image = null;
		if (input2Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input2Image[xyz];
				if (var2Image!=null) noise[n][xyz] = var2Image[xyz];
			}
			scaling[n] = scale2Param;
		}	
		input2Image = null;		
		var2Image = null;		
		if (input3Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input3Image[xyz];
				if (var3Image!=null) noise[n][xyz] = var3Image[xyz];
			}
			scaling[n] = scale3Param;
		}			
		input3Image = null;		
		var3Image = null;		
		
		// two input options: a max proba + labels or a set of probabilities (incl. background or not)
		byte nlabels = -1;
		int[] objlb = null;
		float[][] maxproba = null;
		byte[][] maxlabel = null;
		
		if (maxProba) {
			BasicInfo.displayMessage("build from max labeling and max probability...\n");
			objlb = ObjectLabeling.listOrderedLabels(labelImage,nx,ny,nz);
			nlabels = (byte)objlb.length;
			
			BasicInfo.displayMessage("found "+nlabels+" labels\n");
			// create a distance-based probability map
			float[] boundary = new float[nxyz];
			if (!includeBgParam) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (labelImage[xyz]==0) boundary[xyz] = 1.0f-probaImage[xyz];
					else boundary[xyz] = probaImage[xyz];
				}
			} else {
				for (int xyz=0;xyz<nxyz;xyz++) {
					boundary[xyz] = probaImage[xyz];
				}
			}
			byte[] lbs = new byte[nlabels];
			for (int l=0;l<nlabels;l++) lbs[l] = (byte)objlb[l];
			BasicInfo.displayMessage("distance-based MGDM representation...\n");
			MgdmRepresentation mgdm = new MgdmRepresentation(labelImage, boundary, nx,ny,nz, rx,ry,rz, lbs, nlabels, nbestParam-1, false, 9.0f);
			
			//MaxProbaRepresentation maxprep = new MaxProbaRepresentation(nbestParam, nlabels, nx,ny,nz);
			maxproba = new float[nlabels][nxyz];
			maxlabel = mgdm.getLabels();
			float[] distproba = new float[nbestParam];
			for (int xyz=0;xyz<nxyz;xyz++) {
				maxproba[0][xyz] = boundary[xyz];
				distproba[nbestParam-1] = 0.0f;
				for (int b=0;b<nbestParam-1;b++) {
					distproba[b] = 1.0f/(1.0f+mgdm.getFunctions()[b][xyz]);
					distproba[nbestParam-1] += distproba[b];
				}
				for (int b=0;b<nbestParam-1;b++) {
					maxproba[b+1][xyz] = Numerics.min( boundary[xyz], (1.0f-boundary[xyz])*(distproba[b]/distproba[nbestParam-1]) );
				}
			}			
		} else {
			BasicInfo.displayMessage("build from probabilities...\n");
			// create max proba, labels (assuming a background value with label zero)
			nlabels = (byte)(probaImage.length/nxyz);
			
			if (rescaleProbaParam) {
				float[] minpb = new float[nlabels];
				float[] maxpb = new float[nlabels];
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int l=0;l<nlabels;l++)  {
						if (probaImage[xyz+l*nxyz]<minpb[l]) minpb[l] = probaImage[xyz+l*nxyz];
						if (probaImage[xyz+l*nxyz]>maxpb[l]) maxpb[l] = probaImage[xyz+l*nxyz];
					}
				}
				BasicInfo.displayMessage("probability ranges ("+nlabels+" ):\n");
				for (int l=0;l<nlabels;l++) BasicInfo.displayMessage("["+minpb[l]+", "+maxpb[l]+"] ");
				BasicInfo.displayMessage("\n");
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int l=0;l<nlabels;l++)  {
						probaImage[xyz+l*nxyz] = (probaImage[xyz+l*nxyz]-minpb[l])/(maxpb[l]-minpb[l]);
					}
				}
			}
			
			MaxProbaRepresentation maxprep = new MaxProbaRepresentation(nbestParam, nlabels, nx,ny,nz);
			if (includeBgParam) maxprep.buildFromCompleteProbabilities(probaImage);
			else maxprep.buildFromCompetingProbabilitiesAndBackground(probaImage);
			//else maxprep.buildFromCompetingProbabilitiesAndConstantBackground(probaImage,0.5f);
			maxproba = maxprep.getMaxProba();
			maxlabel = maxprep.getMaxLabel();
			
			if (!includeBgParam) nlabels++;
			
			// get labels from input if given
			if (labelImage==null) {
				objlb = new int[nlabels];
				for (int l=0;l<nlabels;l++) objlb[l] = l+1;
			} else {
				objlb = ObjectLabeling.listOrderedLabels(labelImage,nx,ny,nz);
			}
			BasicInfo.displayMessage("found "+nlabels+" labels\n");

		}
		/*
		// rescale highest values to 1 / lowest value to 0 for each label
		if (rescaleProbaParam) {
			float[] minpb = new float[nlabels];
			float[] maxpb = new float[nlabels];
			for (int xyz=0;xyz<nxyz;xyz++) {
				for (int b=0;b<nbestParam;b++) {
					if (maxproba[b][xyz]<minpb[maxlabel[b][xyz]]) minpb[maxlabel[b][xyz]] = maxproba[b][xyz];
					if (maxproba[b][xyz]>maxpb[maxlabel[b][xyz]]) maxpb[maxlabel[b][xyz]] = maxproba[b][xyz];
				}
			}
			BasicInfo.displayMessage("probability ranges ("+nlabels+" ):\n");
			for (int l=0;l<nlabels;l++) BasicInfo.displayMessage("["+minpb[l]+", "+maxpb[l]+"] ");
			BasicInfo.displayMessage("\n");
			for (int xyz=0;xyz<nxyz;xyz++) {
				for (int b=0;b<nbestParam;b++) {
					maxproba[b][xyz] = (maxproba[b][xyz]-minpb[maxlabel[b][xyz]])/(maxpb[maxlabel[b][xyz]]-minpb[maxlabel[b][xyz]]);
				}
			}
		}
		*/
		
		// main algorithm
		
		StatisticalUncertaintyReduction sur = new StatisticalUncertaintyReduction(image, noise, maskImage, scaling, nimg, nx, ny, nz,  nbestParam);
		sur.setBestProbabilities(maxproba, maxlabel, objlb);
		//if (computeNoiseParam) sur.estimateMeanImageNoise();
		//if (computeNoiseParam) sur.estimateMedianImageNoise();
		if (computeNoiseParam) {
			BasicInfo.displayMessage("estimation of rank-based noise level...\n");
			sur.estimateImageRankScale(neighborParam);
		}
		
		BasicInfo.displayMessage("main certainty diffusion...\n");
		sur.diffuseCertainty(iterationParam, imgscaleParam, certainscaleParam, neighborParam, mincertaintyParam, computeDistributionParam, diffratioParam);
				
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
		
		// map to 1D arrays
		segmentImage = new int[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (sur.getLabels()[0][xyz]>=0) {
				if (!includeBgParam && objlb[sur.getLabels()[0][xyz]]==nlabels) {
					segmentImage[xyz] = 0;
				} else {
					segmentImage[xyz] = objlb[sur.getLabels()[0][xyz]];
				}
			} else {
				segmentImage[xyz] = 0;
			}
		}
		maxprobaImage = new float[nxyz*nbestParam];
		maxidImage = new byte[nxyz*nbestParam];
		for (int l=0;l<nbestParam;l++) for (int xyz=0;xyz<nxyz;xyz++) {
			maxprobaImage[xyz+nxyz*l] = sur.getProbabilities()[l][xyz];
			maxidImage[xyz+nxyz*l] = sur.getLabels()[l][xyz];
		}
	
		// for debug
		debugImage = new float[nxyz*(2*neighborParam+1)];
		float[] tmp = sur.computeMaxCertainty(certainscaleParam);
		for (int xyz=0;xyz<nxyz;xyz++) debugImage[xyz] = tmp[xyz];
		/*
		tmp = sur.computeMaxImageWeight(imgscaleParam);
		for (int xyz=0;xyz<nxyz;xyz++) debugImage[xyz+nxyz] = tmp[xyz];
		tmp = sur.computeMinImageWeight(imgscaleParam);
		for (int xyz=0;xyz<nxyz;xyz++) debugImage[xyz+2*nxyz] = tmp[xyz];
		/*
		float[][] imw = sur.computeAllImageWeight(imgscaleParam);
		for (int xyz=0;xyz<nxyz;xyz++) for (int j=0;j<26;j++) debugImage[xyz+j*nxyz] = imw[j][xyz];
		*/
		float[][] imw = sur.computeBestImageWeight(imgscaleParam,2*neighborParam);
		for (int xyz=0;xyz<nxyz;xyz++) for (int j=0;j<2*neighborParam;j++) debugImage[xyz+nxyz+j*nxyz] = imw[j][xyz];
		
		return;
	}

}
