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
	private byte 	nlabelsParam = 0;
	private byte 	nbestParam	=	4;
	
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
	public final void setLabelNumber(byte val) { nlabelsParam = val; }
	public final void setLabelDepth(byte val) { nbestParam = val; }

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
		
		// create max proba, labels (assuming a background value with label zero)
		MaxProbaRepresentation maxproba = new MaxProbaRepresentation(nbestParam, nlabelsParam, nx,ny,nz);
		maxproba.buildFromCompetingProbabilitiesAndBackground(probaImage);
		
		// main algorithm
		boolean[] used = new boolean[nimg];
		for (int i=0;i<nimg;i++) used[i] = true;
		
		StatisticalUncertaintyReduction sur = new StatisticalUncertaintyReduction(image, used, scaling, nimg, nx, ny, nz,  nbestParam);
		
		int[] objlb = new int[nlabelsParam];
		for (int l=0;l<nlabelsParam;l++) objlb[l] = l+1;
		sur.setBestProbabilities(maxproba.getMaxProba(), maxproba.getMaxLabel(), objlb);
		
		sur.diffuseCertainty(iterationParam, imgscaleParam, certainscaleParam, neighborParam, mincertaintyParam);
					
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
		
		// map to 1D arrays
		segmentImage = sur.getSegmentation();
		maxprobaImage = sur.getProbabilities();
		maxidImage = sur.getLabels();
		
		return;
	}

}
