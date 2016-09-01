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
	
	private float scale1Param =1.0f;
	private float scale2Param =1.0f;
	private float scale3Param =1.0f;	
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private float[] probaImage = null;
	private int[] 	labelImage = null;
	private byte 	nbestParam	=	4;
	private boolean	includeBgParam = false;
	private boolean	rescaleProbaParam = false;
	
	private int 	iterationParam	=	500;
	private float 	imgscaleParam		=	0.05f;
	private	float 	certainscaleParam		= 	0.1f;
	private float 	mincertaintyParam		=	0.25f;
	private int 	neighborParam	=	6;
	
	// outputs
	private int[] segmentImage;
	private float[] maxprobaImage;
	private byte[] maxidImage;

	// create inputs
	public final void setContrastImage1(float[] val) { input1Image = val; }
	public final void setContrastImage2(float[] val) { input2Image = val; }
	public final void setContrastImage3(float[] val) { input3Image = val; }
	
	public final void setContrastScale1(float val) { scale1Param = val; }
	public final void setContrastScale2(float val) { scale2Param = val; }
	public final void setContrastScale3(float val) { scale3Param = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setLabelProbabilities(float[] val) { probaImage = val; }
	public final void setLabelSegmentation(int[] val) { labelImage = val; }
	public final void setLabelDepth(byte val) { nbestParam = val; }
	public final void setBackgroundIncluded(boolean val) { includeBgParam = val; }

	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setImageScale(float val) { imgscaleParam = val; }
	public final void setCertaintyScale(float val) { certainscaleParam = val; }
	public final void setMinCertainty(float val) { mincertaintyParam = val; }
	public final void setNeighborhoodSize(int val) { neighborParam = val; }
	
	
	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Segmentation.devel"; }
	public static final String getLabel() { return "Statistical Uncertainty Reduction (SUR)"; }
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

	public void execute(){
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image != null) nimg++;
		if (input3Image != null) nimg++;
		
		float[][] image = new float[nimg][nxyz];
		float[] scaling = new float[nimg];
		int n = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[n][xyz] = input1Image[xyz];
		}
		scaling[n] = scale1Param;
		input1Image = null;
		if (input2Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input2Image[xyz];
			}
			scaling[n] = scale2Param;
		}	
		input2Image = null;		
		if (input3Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input3Image[xyz];
			}
			scaling[n] = scale3Param;
		}			
		input3Image = null;		
		
		// two input options: a max proba + labels or a set of probabilities (incl. background or not)
		byte nlabels = -1;
		int[] objlb = null;
		float[][] maxproba = null;
		byte[][] maxlabel = null;
		
		if (labelImage==null) {
			BasicInfo.displayMessage("build from probabilities...\n");
			// create max proba, labels (assuming a background value with label zero)
			nlabels = (byte)(probaImage.length/nxyz);
			
			MaxProbaRepresentation maxprep = new MaxProbaRepresentation(nbestParam, nlabels, nx,ny,nz);
			if (includeBgParam) maxprep.buildFromCompleteProbabilities(probaImage);
			else maxprep.buildFromCompetingProbabilitiesAndBackground(probaImage);
			maxproba = maxprep.getMaxProba();
			maxlabel = maxprep.getMaxLabel();
			
			if (!includeBgParam) nlabels++;
			objlb = new int[nlabels];
			for (int l=0;l<nlabels;l++) objlb[l] = l+1;

		} else {
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
		}
		
		// main algorithm
		boolean[] used = new boolean[nimg];
		for (int i=0;i<nimg;i++) used[i] = true;
		
		StatisticalUncertaintyReduction sur = new StatisticalUncertaintyReduction(image, used, scaling, nimg, nx, ny, nz,  nbestParam);
		sur.setBestProbabilities(maxproba, maxlabel, objlb);
		
		sur.diffuseCertainty(iterationParam, imgscaleParam, certainscaleParam, neighborParam, mincertaintyParam);
				
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
		
		return;
	}

}
