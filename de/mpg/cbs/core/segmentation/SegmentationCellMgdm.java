package de.mpg.cbs.core.brain;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationCellMgdm {

	private float[] input1Image = null;
	private float[] input2Image = null;
	private float[] input3Image = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private String type1Param = "none";
	private String type2Param = "none";
	private String type3Param = "none";
	private String type4Param = "none";
	public static final String[] inputTypes = {"centroid-proba", "foreground-proba"};
	
	private String dimParam = "none";
	public static final String[] dimTypes = {"2D", "3D"};
													
	private int 	iterationParam	=	500;
	private float 	changeParam		=	0.001f;
	private	float 	forceParam		= 	0.1f;
	private float 	curvParam		=	0.4f;
		
	private String 	topologyParam	=	"wcs";
	public static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	private String	lutdir = null;
	
	// outputs
	private int[] segmentImage;
	private float[] mgdmImage;
	
	private float[] membershipImage;
	private int[] labelImage;
	private int output4Dlength;
	
	// create inputs
	public final void setContrastImage1(float[] val) { input1Image = val; }
	public final void setContrastType1(String val) { type1Param = val; }
	public final void setContrastImage2(float[] val) { input2Image = val; }
	public final void setContrastType2(String val) { type2Param = val; }
	public final void setContrastImage3(float[] val) { input3Image = val; }
	public final void setContrastType3(String val) { type3Param = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setDataWeight(float val) {forceParam = val; }
	public final void setCurvatureWeight(float val) { curvParam = val; }
	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setMinChange(float val) { changeParam = val; }
	public final void setSteps(int val) { stepParam = val; }
	
	public final void setTopology(String val) { topologyParam = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }
	
	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Segmentation.devel"; }
	public static final String getLabel() { return "MGDM Cell Segmentation"; }
	public static final String getName() { return "MGDMCellSegmentation"; }

	public static final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public static final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public static final String getDescription() { return "MGDM estimation of cells from deep learning inputs."; }
	public static final String getLongDescription() { return getDescription(); }
		
	public static final String getVersion() { return "3.1.2"; }

	// create outputs
	public final int[] getSegmentedBrainImage() { return segmentImage; }
	public final float[] getLevelsetBoundaryImage() { return mgdmImage; }
	public final int getOutput4Dlength() { return output4Dlength; }
	public final float[] getPosteriorMaximumMemberships4D() { return membershipImage; }
	public final byte[] getPosteriorMaximumLabels4D() { return labelImage; }

	public void execute(){
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image != null) nimg++;
		if (input3Image != null) nimg++;
		
		String[] modality = new String[nimg];
		float[][] image = new float[nimg][];
		int n=0;
		modality[n] = type1Param;
		image[n] = input1Image;
		if (input2Image != null) { n++; modality[n] = type2Param; image[n] = input2Image; }
		if (input3Image != null) { n++; modality[n] = type3Param; image[n] = input3Image; }
		
		float[] centroids = null;
		float[] fgproba = null;
		for (int n=0;n<nimg;n++) {
		    if (modality[n].equals("centroid-proba")) centroids = image[n];
		    else if (modality[n].equals("foreground-proba")) fgproba = image[n];
		}
		
		// different streams for 2D or 3D
		byte DIM2D = 2;
		byte DIM3D = 3;
		byte dimension = DIM2D;
		if (dimParam.equals("3D")) dimension = DIM3D;
		
		// 1. Get the initial segmentation from cell centroid detection
		int[] initialization = null;
		if (dimension==DIM2D) {
		    // independet Z stacks
		    for (int z=0;z<nz;z++) {
		        int[] labels = ObjectLabeling.connected6Object3D(Object.objectFromImage(centroids, nx,ny,1, 0.0f, ObjectExtraction.SUPERIOR), nx,ny,1);
		        //...
		    }
		}
		
		// 2. Get the forces from foreground proba
		
		// 3. Run MGDM!
		
		float[][] image = new float[nimg][nxyz];
		n = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[n][xyz] = input1Image[xyz];
		}
		input1Image = null;
		if (input2Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input2Image[xyz];
			}
		}	
		input2Image = null;		
		if (input3Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input3Image[xyz];
			}
		}			
		input3Image = null;		
		if (input4Image != null) {
			n++;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = input4Image[xyz];
			}
		}	
		input4Image = null;		
		
		// main algorithm
		BasicInfo.displayMessage("Load atlas\n");

		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam);
		
		//adjust modalities to canonical names
		BasicInfo.displayMessage("Image contrasts:\n");
		for (n=0;n<nimg;n++) {
			modality[n] = atlas.displayContrastName(atlas.contrastId(modality[n]));
			BasicInfo.displayMessage(modality[n]+"\n");
		}
		
		atlas.setImageInfo(nx, ny, nz, rx, ry, rz, orient, orx, ory, orz);
		atlas.adjustAtlasScale(image, nimg);
		atlas.setQuantitativeNormalization(normalizeQuantitative);
		
		atlas.initShapeMapping();
		
		BasicInfo.displayMessage("Compute tissue classification\n");

		ShapeAtlasClassification2 classif = new ShapeAtlasClassification2(image, modality, 
																		  nimg, nx,ny,nz, rx,ry,rz, atlas);
		classif.initialAtlasTissueCentroids();
		classif.computeMemberships();
		
		BasicInfo.displayMessage("First alignment\n");
		
		atlas.alignObjectCenter(classif.getMemberships()[ShapeAtlasClassification.WM], "wm");
		atlas.refreshShapeMapping();
		BasicInfo.displayMessage("transform: "+atlas.displayTransform(atlas.getTransform()));
			
		classif.computeMemberships();
		
		if (adjustIntensPriors) {
			float diff = 1.0f;
			for (int t=0;t<20 && diff>0.01f;t++) {
				diff = classif.computeCentroids();
				BasicInfo.displayMessage("iteration "+t+", max diff: "+diff+"\n");
				classif.computeMemberships();
			}
		}
		
		BasicInfo.displayMessage("Rigid alignment\n");
		
		BasicRigidRegistration rigid = new BasicRigidRegistration(atlas.generateObjectImage("wm"), 
																	classif.getMemberships()[ShapeAtlasClassification.WM], 
																	atlas.getShapeDim()[0],atlas.getShapeDim()[1],atlas.getShapeDim()[2],
																	atlas.getShapeRes()[0],atlas.getShapeRes()[1],atlas.getShapeRes()[2],
																	50, 0.0f, 1);
		
		rigid.register();
		atlas.updateRigidTransform(rigid.getTransform());
		atlas.refreshShapeMapping();

		// clean-up
		rigid.finalize();
		rigid = null;
		
		classif.computeMemberships();

		if (adjustIntensPriors) {
			float diff = 1.0f;
			for (int t=0;t<20 && diff>0.01f;t++) {
				diff = classif.computeCentroids();
				BasicInfo.displayMessage("iteration "+t+", max diff: "+diff+"\n");
				classif.computeMemberships();
			}
		}
		classif.estimateIntensityTransform();
		
		byte[][][] tissueseg = classif.exportClassificationByte();
		
		// first order warping
		BasicInfo.displayMessage("NL warping\n");
		
		BasicDemonsWarping warp1 = new BasicDemonsWarping(atlas.generateDifferentialObjectSegmentation("wm","mask"), 
															classif.getDifferentialSegmentation(classif.WM, classif.BG),
															atlas.getShapeDim()[0],atlas.getShapeDim()[1],atlas.getShapeDim()[2],
															atlas.getShapeRes()[0],atlas.getShapeRes()[1],atlas.getShapeRes()[2],
															1.0f, 4.0f, 1.0f, BasicDemonsWarping.DIFFUSION,
															null);
							
		warp1.initializeTransform();
		
		for (int t=0;t<50;t++) {
			warp1.registerImageToTarget();
		}
		
		atlas.updateNonRigidTransform(warp1);
		
			
		// record corresponding segmentation
		 byte[] target = atlas.generateTransformedClassification();
		
		// clean-up
		warp1.finalize();
		warp1 = null;
		
		BasicInfo.displayMessage("level set segmentation\n");

		/*
		// minimum scale?
		float factor = (float)Math.sqrt(3.0)*Numerics.max(rx/atlas.getShapeRes()[0],ry/atlas.getShapeRes()[1],rz/atlas.getShapeRes()[2]);
		BasicInfo.displayMessage("minimum scale: "+factor+"\n");
		
		// or full scale?
		factor = 1.0f;
		*/
		
		int nmgdm = 5;
		int ngain = 5;
		
		// scale the force with regard to the final level
		float factor = 1.0f/(Numerics.max(rx/atlas.getShapeRes()[0],ry/atlas.getShapeRes()[1],rz/atlas.getShapeRes()[2]));
		BasicInfo.displayMessage("atlas scale: "+factor+"\n");
			
		// constant scale for the distance (=> multiplied by scaling factor)
		float distanceScale = scaleParam/rx;
			
		// count the number of pre-processing MGDMs already done: skip the lowest smoothness steps in the following ones
		int nprocessed = 0;
			
		MgdmFastAtlasSegmentation2 mgdma = new MgdmFastAtlasSegmentation2(image, modality, classif.getImageRanges(), nimg,
																		nx,ny,nz, rx,ry,rz, 
																		atlas,
																		nmgdm, ngain,
																		forceParam*factor,
																		curvParam,
																		0.0f, 0.0f,
																		distanceScale,
																		"wcs", lutdir);
				
		BasicInfo.displayMessage("gain...\n");
		
		mgdma.computeAtlasBestGainFunction();
		
		if (diffuseProbabilities)
			mgdma.diffuseBestGainFunctions(20, 0.5f, diffuseParam);
		
		if (stepParam>0) {
			BasicInfo.displayMessage("atlas-space levelset evolution...\n");
			mgdma.evolveNarrowBand(iterationParam,changeParam);
			
			nprocessed++;
		}
		
		atlas.setTemplate(mgdma.getLabeledSegmentation());
		
		// minimum scale?
		factor = 1.0f/((float)Math.sqrt(3.0)*Numerics.max(rx/atlas.getShapeRes()[0],ry/atlas.getShapeRes()[1],rz/atlas.getShapeRes()[2]));
		BasicInfo.displayMessage("minimum scale: "+factor+"\n");
		
		//short[] counter = mgdma.exportFrozenPointCounter();
		byte[] atlasseg = mgdma.getSegmentation();
		byte[] initseg = null;
		
		MgdmFastScaledSegmentation2 mgdms = null;
		// make a progressive scale increase/decrease
		for (int nscale=0;nscale<10 && factor>1.75f;nscale++) {
			mgdms = new MgdmFastScaledSegmentation2(image, modality, classif.getImageRanges(), nimg,
																			nx,ny,nz, rx,ry,rz, 
																			atlas, atlasseg, initseg,
																			nmgdm, ngain, factor,
																			forceParam*factor, 
																			curvParam,
																			0.0f,
																			distanceScale,
																			"wcs", lutdir);
					
			BasicInfo.displayMessage("gain...\n");
			
			mgdms.importBestGainFunction(mgdma.getBestGainFunctionHD(), mgdma.getBestGainLabelHD());
			
			//BasicInfo.displayMessage("scaled-space levelset evolution...\n");
			if (nprocessed<stepParam) {
				BasicInfo.displayMessage("scaled-space levelset evolution...\n");
				mgdms.evolveNarrowBand(iterationParam,changeParam);

				nprocessed++;
			}
				
			initseg = mgdms.exportScaledSegmentation();
			// scaling factor must be lower than 1/sqrt(3) to preserve topology
			factor *= 0.575f;
			BasicInfo.displayMessage("additional scaling? (new factor: "+factor+")\n");
		}

		MgdmFastSegmentation2 mgdm = new MgdmFastSegmentation2(image, modality, classif.getImageRanges(), nimg,
																nx,ny,nz, rx,ry,rz, 
																atlas, atlasseg, initseg,
																nmgdm, ngain,
																forceParam, 
																curvParam,
																0.0f,
																distanceScale,
																topologyParam, lutdir);
		
		// clean-up
		classif.finalize();
		classif = null;
		
		BasicInfo.displayMessage("gain...\n");
			
		mgdm.importBestGainFunctions(mgdma.getBestGainFunctionHD(), mgdma.getBestGainLabelHD());
			
		//float[] initgain = mgdm.exportBestGainSegmentation();
		//float[][][] initmems = mgdm.exportBestGainFunction();
			
		if (nprocessed<stepParam) {
			BasicInfo.displayMessage("full scale levelset evolution...\n");
			mgdm.evolveNarrowBand(iterationParam,changeParam);
		}
		
		//float[][][] evolmems = mgdm.exportBestGainFunction();
		
		BasicInfo.displayMessage("partial volume estimates...\n");
			
		if (computePosteriors) {
			mgdm.computeApproxPartialVolumes(distanceScale, false);
		}
					
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
			
		segmentImage = mgdm.labelSegmentation();
		BasicInfo.displayMessage("segmentation");
		
		mgdmImage = mgdm.getFunctions()[0];
		BasicInfo.displayMessage(".. boundaries");
		
		idImage = mgdm.getSegmentation();
		BasicInfo.displayMessage("segmentation ids");
		
		if (outputParam.equals("label_memberships")) {		
			output4Dlength = ngain+1;
			
			membershipImage = mgdm.exportBestGainFunctions1D(0, ngain, !computePosteriors);
			BasicInfo.displayMessage(".. memberships(4d)");
			
			labelImage = mgdm.exportBestGainLabelsByte1D(0, ngain);
			BasicInfo.displayMessage(".. labels(4d)");
		} else if (outputParam.equals("raw_memberships")) {		
			output4Dlength = ngain+1;
			
			membershipImage = mgdm.exportBestGainFunctions1D(0, ngain, !computePosteriors);
			BasicInfo.displayMessage(".. memberships(4d)");
			
			labelImage = mgdm.exportBestGainsByte1D(0, ngain);
			BasicInfo.displayMessage(".. labels(4d)");
		} else if (outputParam.equals("segmentation")) {		
			output4Dlength = 1;
			
			// get the best label map and probabilities by default
			membershipImage = mgdm.exportBestGainFunctions1D(0, 0, !computePosteriors);
			BasicInfo.displayMessage(".. best membership");
			
			labelImage = mgdm.exportBestGainsByte1D(0, 0);
			BasicInfo.displayMessage(".. best label");
		} else if (outputParam.equals("debug")) {			
			output4Dlength = 3;
			
			membershipImage = mgdm.exportBestGainFunctions1D(0, 3, !computePosteriors);
			BasicInfo.displayMessage(".. memberships(4d)");
			
			//byte[][][][] byte4d = mgdm.exportBestGainLabelsByte(0, ngain);
			labelImage = new byte[nx*ny*nz*3];
			byte[][][] gain = mgdm.exportBestGainByte();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				labelImage[xyz] = tissueseg[x][y][z];
				labelImage[xyz+nxyz] = target[xyz];
				labelImage[xyz+2*nxyz] = gain[x][y][z];
			}
			BasicInfo.displayMessage(".. debug(4d)");
		}
		return;
	}

}
