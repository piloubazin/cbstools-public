package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistBrainMgdmMultiSegmentationUpdate2 extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;
	private ParamVolume input4Image;
	
	private ParamVolume priorImage;
	
	private ParamOption type1Param;
	private ParamOption type2Param;
	private ParamOption type3Param;
	private ParamOption type4Param;
	private static final String[] inputTypes = {"-- 3T --", "MPRAGE3T", "T1MAP3T", "MP2RAGE3T", 
													"HCPT1w", "HCPT2w", "NormMPRAGE", "FLAIR3T",
													"-- 7T --", "T1MAP7T", "MP2RAGE7T", "T2SW7T", "QSM7T", 
													"-- 9.4fT --", "T1MAP9T", "MP2RAGE9T",
													"-- misc --", "Filters", "PVDURA", "none"};
	
	private static final String inputTypeInfo = "Currently available contrasts:\n"
			+"T1MAP7T: a 7T quantitative T1 map, \n"
			+"MP2RAGE7T: a T1-weighted image from 7T MP2RAGE (UNI), \n"
			+"T2SW7T: a 7T T2*-weighted image (devel), \n"
			+"SQSM7T: a 7T quantitative susceptibility map (devel), \n"
			+"T1MAP9T: a 9.4fT quantitative T1 map (devel), \n"
			+"MP2RAGE9T: a T1-weighted image from 9.4fT MP2RAGE (UNI) (devel), \n"
			+"MPRAGE3T: a 3T T1-weighted MPRAGE image, \n"
			+"T1MAP3T: a 3T quantitative T1 map, \n"
			+"MP2RAGE3T: a T1-weighted image from MP2RAGE (UNI), \n"
			+"HCPT1w: a T1-weighted image using the HCP sequence, \n"
			+"HCPT2w: a T2-weighted image using the HCP sequence, \n"
			+"NormMPRAGE: a 3T MPRAGE normalised for B1 shading, \n"
			+"FLAIR3T: a 3T FLAIR image, \n"
			+"Filters: a composite image of outputs from dura, pv and arteries pre-processing, \n"
			+"PVDURA: a composite image of outputs from dura and pv (obsolete).\n";												
	
	private ParamFile atlasParam;

	private ParamInteger 	iterationParam;
	private ParamInteger 	stepParam;
	private ParamFloat 	changeParam;
	private	ParamFloat 	forceParam;
	private ParamFloat 	curvParam;
	private ParamFloat 	scaleParam;
	
	private ParamBoolean	computePosteriors;
	private ParamBoolean	adjustIntensPriors;
	private ParamBoolean	diffuseProbabilities;
	private ParamFloat		diffuseParam;
	
	private ParamOption 	topologyParam;
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamOption 	outputParam;
	private static final String[] outputTypes = {"segmentation","label_memberships","raw_memberships","debug"};
	//private static final String[] outputTypes = {"segmentation","memberships","cortex"};
	
	private ParamVolume segmentImage;
	private ParamVolume mgdmImage;
	private ParamVolume idImage;
	
	private ParamVolume membershipImage;
	private ParamVolume labelImage;
	
	/*
	private ParamVolume initlImage;
	private ParamVolume initrImage;
	
	private ParamVolume wmlImage;
	private ParamVolume wmrImage;
	private ParamVolume gmlImage;
	private ParamVolume gmrImage;
	private ParamVolume csflImage;
	private ParamVolume csfrImage;
	
	private ParamVolume cwmImage;
	private ParamVolume cgmImage;
	private ParamVolume scsfImage;
	*/
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Data");
		
		imageParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		imageParams.add(type1Param = new ParamOption("Contrast Type 1",inputTypes));
		type1Param.setValue("T1MAP7T");
		type1Param.setDescription("Currently available contrasts:\n"
			+"T1MAP7T: a 7T quantitative T1 map\n"
			+"MP2RAGE7T: a T1-weighted image from MP2RAGE (UNI)\n"
			+"T2SW7T: a 7T T2*-weighted image (devel)\n"
			+"SQSM7T: a 7T quantitative susceptibility map (devel)\n"
			+"MPRAGE3T: a 3T T1-weighted MPRAGE image\n"
			+"T1MAP3T: a 3T quantitative T1 map\n"
			+"MP2RAGE3T: a T1-weighted image from MP2RAGE (UNI)\n"
			+"HCPT1w: a T1-weighted image using the HCP sequence\n"
			+"HCPT2w: a T2-weighted image using the HCP sequence\n"
			+"NormMPRAGE: a 3T MPRAGE normalised for B1 shading\n"
			+"FLAIR3T: a 3T FLAIR image\n"
			+"Filters: a composite image of outputs from dura, pv and arteries pre-processing\n"
			+"PVDURA: a composite image of outputs from dura and pv (obsolete)\n"
			);
		
		imageParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		imageParams.add(type2Param = new ParamOption("Contrast Type 2",inputTypes));
		input2Image.setMandatory(false);
		type2Param.setValue("Filters");
		
		imageParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		imageParams.add(type3Param = new ParamOption("Contrast Type 3",inputTypes));
		input3Image.setMandatory(false);
		type3Param.setValue("MP2RAGE7T");
		
		imageParams.add(input4Image = new ParamVolume("Contrast Image 4 (opt)"));
		imageParams.add(type4Param = new ParamOption("Contrast Type 4",inputTypes));
		input4Image.setMandatory(false);
		type4Param.setValue("T2SW7T");
		
		imageParams.add(priorImage = new ParamVolume("Prior Segmentation"));
		
		imageParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		
		mainParams.add(forceParam = new ParamFloat("Data weight", -1E10f, 1E10f, 0.1f));
		mainParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.4f));
		mainParams.add(scaleParam = new ParamFloat("Posterior scale (mm)", -1E10f, 1E10f, 5.0f));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(changeParam = new ParamFloat("Min change", 0f, 1f, 0.001f));
		mainParams.add(stepParam = new ParamInteger("Steps", 0, 100000, 5));
		
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(computePosteriors = new ParamBoolean("Compute posteriors", true));
		mainParams.add(adjustIntensPriors = new ParamBoolean("Adjust intensity priors", true));
		mainParams.add(diffuseProbabilities = new ParamBoolean("Diffuse probabilities", true));
		mainParams.add(diffuseParam = new ParamFloat("Diffusion factor", 0f, 10f, 0.5f));
		
		mainParams.add(outputParam = new ParamOption("Output images", outputTypes));
		outputParam.setValue("label_memberships");

		inputParams.add(mainParams);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing.devel");
		inputParams.setLabel("MGDM Segmentation Update 2");
		inputParams.setName("MGDMSegmentationUpdate2");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimate brain structures from an atlas for a MRI dataset (multiple input combinations are possible), starting from a given segmentation.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(segmentImage = new ParamVolume("Segmented Brain Image",VoxelType.BYTE));
		outputParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image",VoxelType.BYTE));
		outputParams.add(idImage = new ParamVolume("Segmentation Ids Image",VoxelType.BYTE));
		
		outputParams.add(membershipImage = new ParamVolume("Posterior Maximum Memberships (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(labelImage = new ParamVolume("Posterior Maximum Labels (4D)",VoxelType.INT,-1,-1,-1,-1));
		membershipImage.setMandatory(false);
		labelImage.setMandatory(false);

		outputParams.setName("brain segmentation images");
		outputParams.setLabel("brain segmentation images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat in1Img = new ImageDataFloat(input1Image.getImageData());
		int nx = in1Img.getRows();
		int ny = in1Img.getCols();
		int nz = in1Img.getSlices();
		int nxyz = nx*ny*nz;
		float rx = in1Img.getHeader().getDimResolutions()[0];
		float ry = in1Img.getHeader().getDimResolutions()[1];
		float rz = in1Img.getHeader().getDimResolutions()[2];
		
		int orient = in1Img.getHeader().getImageOrientation().ordinal();
		int orx = in1Img.getHeader().getAxisOrientation()[0].ordinal();
		int ory = in1Img.getHeader().getAxisOrientation()[1].ordinal();
		int orz = in1Img.getHeader().getAxisOrientation()[2].ordinal();
		
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image.getImageData() != null) nimg++;
		if (input3Image.getImageData() != null) nimg++;
		if (input4Image.getImageData() != null) nimg++;
		
		String[] modality = new String[nimg];
		int n=0;
		modality[n] = type1Param.getValue();
		if (input2Image.getImageData() != null) { n++; modality[n] = type2Param.getValue(); }
		if (input3Image.getImageData() != null) { n++; modality[n] = type3Param.getValue(); }
		if (input4Image.getImageData() != null) { n++; modality[n] = type4Param.getValue(); }
		
		float[][] image = new float[nimg][nxyz];
		float[][][] buffer;
		n = 0;
		buffer = in1Img.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			image[n][xyz] = buffer[x][y][z];
		}
		buffer = null;
		
		if (input2Image.getImageData() != null) {
			n++;
			buffer = (new ImageDataFloat(input2Image.getImageData())).toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = buffer[x][y][z];
			}
			buffer = null;
		}			
		if (input3Image.getImageData() != null) {
			n++;
			buffer = (new ImageDataFloat(input3Image.getImageData())).toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = buffer[x][y][z];
			}
			buffer = null;
		}			
		if (input4Image.getImageData() != null) {
			n++;
			buffer = (new ImageDataFloat(input4Image.getImageData())).toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				image[n][xyz] = buffer[x][y][z];
			}
			buffer = null;
		}		
		
		int[] priorseg = new int[nx*ny*nz];
		int[][][] intbuffer = (new ImageDataInt(priorImage.getImageData())).toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			priorseg[xyz] = intbuffer[x][y][z];
		}
		intbuffer = null;
		
		// main algorithm
		BasicInfo.displayMessage("Load atlas\n");

		SimpleShapeAtlas2 atlas = new SimpleShapeAtlas2(atlasParam.getValue().getAbsolutePath());
		
		//adjust modalities to canonical names
		BasicInfo.displayMessage("Image contrasts:\n");
		for (n=0;n<nimg;n++) {
			modality[n] = atlas.displayContrastName(atlas.contrastId(modality[n]));
			BasicInfo.displayMessage(modality[n]+"\n");
		}
		
		atlas.setImageInfo(nx, ny, nz, rx, ry, rz, orient, orx, ory, orz);
		atlas.adjustAtlasScale(image, nimg);
			
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
		
		if (adjustIntensPriors.getValue().booleanValue()) {
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

		if (adjustIntensPriors.getValue().booleanValue()) {
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
		float factor = (float)Math.sqrt(3.0f)*Numerics.max(rx/atlas.getShapeRes()[0],ry/atlas.getShapeRes()[1],rz/atlas.getShapeRes()[2]);
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
		float distanceScale = scaleParam.getValue().floatValue()/rx;
			
		// count the number of pre-processing MGDMs already done: skip the lowest smoothness steps in the following ones
		int nprocessed = 0;
			
		MgdmFastAtlasSegmentation2 mgdma = new MgdmFastAtlasSegmentation2(image, modality, classif.getImageRanges(), nimg,
																		nx,ny,nz, rx,ry,rz, 
																		atlas,
																		nmgdm, ngain,
																		forceParam.getValue().floatValue()*factor,
																		curvParam.getValue().floatValue(),
																		0.0f, 0.0f,
																		distanceScale,
																		"wcs",null);
				
		BasicInfo.displayMessage("gain...\n");
		
		mgdma.computeAtlasBestGainFunction();
		
		if (diffuseProbabilities.getValue().booleanValue())
			mgdma.diffuseBestGainFunctions(20, 0.5f, diffuseParam.getValue().floatValue());
		
		/* all this is now bypassed by the prior segmentation
		if (stepParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("atlas-space levelset evolution...\n");
			mgdma.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());
			
			nprocessed++;
		}
		
		atlas.setTemplate(mgdma.getLabeledSegmentation());
		
		// minimum scale?
		factor = 1.0f/((float)Math.sqrt(3.0f)*Numerics.max(rx/atlas.getShapeRes()[0],ry/atlas.getShapeRes()[1],rz/atlas.getShapeRes()[2]));
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
																			forceParam.getValue().floatValue()*factor, 
																			curvParam.getValue().floatValue(),
																			0.0f,
																			distanceScale,
																			"wcs");
					
			BasicInfo.displayMessage("gain...\n");
			
			mgdms.importBestGainFunction(mgdma.getBestGainFunctionHD(), mgdma.getBestGainLabelHD());
			
			//BasicInfo.displayMessage("scaled-space levelset evolution...\n");
			if (nprocessed<stepParam.getValue().intValue()) {
				BasicInfo.displayMessage("scaled-space levelset evolution...\n");
				mgdms.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());

				nprocessed++;
			}
				
			initseg = mgdms.exportScaledSegmentation();
			// scaling factor must be lower than 1/sqrt(3) to preserve topology
			factor *= 0.5f75f;
			BasicInfo.displayMessage("additional scaling? (new factor: "+factor+")\n");
		}
		*/
		// prior segmentation : from labels to seg
		int maxlb = 0;
		for (int xyz=0;xyz<nx*ny*nz;xyz++) if (priorseg[xyz]>maxlb) maxlb = priorseg[xyz];
		byte[] invertlb = new byte[maxlb+1];
		for (byte l=0; l<atlas.getNumber(); l++) invertlb[atlas.getLabels()[l]] = l;
			
		byte[] initseg = new byte[nx*ny*nz];
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			initseg[xyz] = invertlb[priorseg[xyz]];
		}
		priorseg = null;
		invertlb = null;

		MgdmFastSegmentation2 mgdm = new MgdmFastSegmentation2(image, modality, classif.getImageRanges(), nimg,
																nx,ny,nz, rx,ry,rz, 
																atlas, null, initseg,
																nmgdm, ngain,
																forceParam.getValue().floatValue(), 
																curvParam.getValue().floatValue(),
																0.0f,
																distanceScale,
																topologyParam.getValue(),null);
		
		// clean-up
		classif.finalize();
		classif = null;
		
		BasicInfo.displayMessage("gain...\n");
			
		mgdm.importBestGainFunctions(mgdma.getBestGainFunctionHD(), mgdma.getBestGainLabelHD());
			
		//float[] initgain = mgdm.exportBestGainSegmentation();
		//float[][][] initmems = mgdm.exportBestGainFunction();
			
		if (nprocessed<stepParam.getValue().intValue()) {
			BasicInfo.displayMessage("full scale levelset evolution...\n");
			mgdm.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());
		}
		
		//float[][][] evolmems = mgdm.exportBestGainFunction();
		
		BasicInfo.displayMessage("partial volume estimates...\n");
			
		if (computePosteriors.getValue().booleanValue()) {
			mgdm.computeApproxPartialVolumes(distanceScale, false);
		}
					
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
			
		String imgname = in1Img.getName();
		
		byte[][][] seg = new byte[nx][ny][nz];
		float[] img = mgdm.exportMap(mgdm.labelSegmentation());
		//float[] img = mgdm.exportMap(mgdm.getSegmentation());
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = (byte)img[xyz];
		}
		img = null;
		ImageDataByte segData = new ImageDataByte(seg);	
		seg = null;
		segData.setHeader(in1Img.getHeader());
		segData.setName(imgname+"_seg");
		segmentImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("segmentation");
		
		float[][][] lvl = new float[nx][ny][nz];
		float[] fcn = mgdm.getFunctions()[0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lvl[x][y][z] = fcn[xyz];
		}
		img = null;
		ImageDataFloat lvlData = new ImageDataFloat(lvl);	
		lvl = null;
		lvlData.setHeader(in1Img.getHeader());
		lvlData.setName(imgname+"_mgdm");
		mgdmImage.setValue(lvlData);
		lvlData = null;
		BasicInfo.displayMessage(".. boundaries");
		
		byte[][][] ids = new byte[nx][ny][nz];
		//float[] img = mgdm.exportMap(mgdm.labelSegmentation());
		img = mgdm.exportMap(mgdm.getSegmentation());
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			ids[x][y][z] = (byte)img[xyz];
		}
		img = null;
		ImageDataByte idsData = new ImageDataByte(ids);	
		ids = null;
		idsData.setHeader(in1Img.getHeader());
		idsData.setName(imgname+"_segid");
		idImage.setValue(idsData);
		idsData = null;
		BasicInfo.displayMessage("segmentation ids");
		
		if (outputParam.getValue().equals("label_memberships")) {		
			/* debug */
			float[][][][] buffer4d = mgdm.exportBestGainFunctions(0, ngain, !computePosteriors.getValue().booleanValue());
			ImageDataFloat bufferData = new ImageDataFloat(buffer4d);
			buffer4d = null;
			bufferData.setHeader(in1Img.getHeader());
			bufferData.setName(imgname+"_mems");
			membershipImage.setValue(bufferData);
			bufferData = null;
			BasicInfo.displayMessage(".. memberships(4d)");
			
			byte[][][][] byte4d = mgdm.exportBestGainLabelsByte(0, ngain);
			ImageDataByte byteData = new ImageDataByte(byte4d);
			byte4d = null;
			byteData.setHeader(in1Img.getHeader());
			byteData.setName(imgname+"_lbls");
			labelImage.setValue(byteData);
			byteData = null;			
			BasicInfo.displayMessage(".. labels(4d)");
		} else if (outputParam.getValue().equals("raw_memberships")) {		
			/* debug */
			float[][][][] buffer4d = mgdm.exportBestGainFunctions(0, ngain, !computePosteriors.getValue().booleanValue());
			ImageDataFloat bufferData = new ImageDataFloat(buffer4d);
			buffer4d = null;
			bufferData.setHeader(in1Img.getHeader());
			bufferData.setName(imgname+"_mems");
			membershipImage.setValue(bufferData);
			bufferData = null;
			BasicInfo.displayMessage(".. memberships(4d)");
			
			byte[][][][] byte4d = mgdm.exportBestGainsByte(0, ngain);
			ImageDataByte byteData = new ImageDataByte(byte4d);
			byte4d = null;
			byteData.setHeader(in1Img.getHeader());
			byteData.setName(imgname+"_lbls");
			labelImage.setValue(byteData);
			byteData = null;			
			BasicInfo.displayMessage(".. labels(4d)");
		} else if (outputParam.getValue().equals("segmentation")) {		
			// get the best label map and probabilities by default
			float[][][] mem = mgdm.exportBestGainFunction();
			ImageDataFloat bufferData = new ImageDataFloat(mem);
			mem = null;
			bufferData.setHeader(in1Img.getHeader());
			bufferData.setName(imgname+"_mem1");
			membershipImage.setValue(bufferData);
			bufferData = null;
			BasicInfo.displayMessage(".. best membership");
			
			//byte[][][] lbl = mgdm.exportBestGainLabelByte();
			byte[][][] lbl = mgdm.exportBestGainByte();
			ImageDataByte byteData = new ImageDataByte(lbl);
			lbl = null;
			byteData.setHeader(in1Img.getHeader());
			byteData.setName(imgname+"_lbl1");
			labelImage.setValue(byteData);
			byteData = null;			
			BasicInfo.displayMessage(".. best label");
		} else if (outputParam.getValue().equals("debug")) {			
			/* debug */
			float[][][][] buffer4d = mgdm.exportBestGainFunctions(0, ngain, !computePosteriors.getValue().booleanValue());
			ImageDataFloat bufferData = new ImageDataFloat(buffer4d);
			buffer4d = null;
			bufferData.setHeader(in1Img.getHeader());
			bufferData.setName(imgname+"_mems");
			membershipImage.setValue(bufferData);
			bufferData = null;
			BasicInfo.displayMessage(".. memberships(4d)");
			
			//byte[][][][] byte4d = mgdm.exportBestGainLabelsByte(0, ngain);
			byte[][][][] byte4d = new byte[nx][ny][nz][3];
			byte[][][] gain = mgdm.exportBestGainByte();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				byte4d[x][y][z][0] = tissueseg[x][y][z];
				byte4d[x][y][z][1] = target[xyz];
				byte4d[x][y][z][2] = gain[x][y][z];
			}
			ImageDataByte byteData = new ImageDataByte(byte4d);
			byte4d = null;
			byteData.setHeader(in1Img.getHeader());
			byteData.setName(imgname+"_debug");
			labelImage.setValue(byteData);
			byteData = null;			
			BasicInfo.displayMessage(".. debug(4d)");
		}
		return;
	}

}
