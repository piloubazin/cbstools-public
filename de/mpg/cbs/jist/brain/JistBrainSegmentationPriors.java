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
public class JistBrainSegmentationPriors extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume input1Image;
	private ParamVolume input2Image;
	private ParamVolume input3Image;
	private ParamVolume input4Image;
	
	private ParamOption type1Param;
	private ParamOption type2Param;
	private ParamOption type3Param;
	private ParamOption type4Param;
	private static final String[] inputTypes = {"-- 3T --", "MPRAGE3T", "T1MAP3T", "MP2RAGE3T", 
													"HCPT1w", "HCPT2w", "NormMPRAGE", "FLAIR3T",
													"DWIFA3T", "DWIMD3T",
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
			+"dwiFA: a DWI fractional anisotropy image, \n"
			+"dwiMD: a DWI mean diffusivity image, \n"
			+"Filters: a composite image of outputs from dura, pv and arteries pre-processing, \n"
			+"PVDURA: a composite image of outputs from dura and pv (obsolete).\n";												
													
	private ParamFile atlasParam;
	private ParamBoolean adjustIntensPriors;

	private ParamVolume labelImage;
	private ParamVolume regionImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(input1Image = new ParamVolume("Contrast Image 1"));
		inputParams.add(type1Param = new ParamOption("Contrast Type 1",inputTypes));
		type1Param.setValue("none");
		type1Param.setDescription(inputTypeInfo);
		
		inputParams.add(input2Image = new ParamVolume("Contrast Image 2 (opt)"));
		inputParams.add(type2Param = new ParamOption("Contrast Type 2",inputTypes));
		input2Image.setMandatory(false);
		type2Param.setValue("none");
		type2Param.setDescription(inputTypeInfo);
		
		inputParams.add(input3Image = new ParamVolume("Contrast Image 3 (opt)"));
		inputParams.add(type3Param = new ParamOption("Contrast Type 3",inputTypes));
		input3Image.setMandatory(false);
		type3Param.setValue("none");
		type3Param.setDescription(inputTypeInfo);
		
		inputParams.add(input4Image = new ParamVolume("Contrast Image 4 (opt)"));
		inputParams.add(type4Param = new ParamOption("Contrast Type 4",inputTypes));
		input4Image.setMandatory(false);
		type4Param.setValue("none");
		type4Param.setDescription(inputTypeInfo);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		inputParams.add(adjustIntensPriors = new ParamBoolean("Adjust intensity priors", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing.devel");
		inputParams.setLabel("Brain Segmentation Priors");
		inputParams.setName("BrainSegmentationPriors");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Generates priors from an atlas for a MRI dataset.");
		info.setLongDescription("Generates priors from an atlas for a MRI dataset.\n"+inputTypeInfo);
		
		info.setVersion("3.1.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(labelImage = new ParamVolume("Prior Region Label",VoxelType.BYTE));
		outputParams.add(regionImage = new ParamVolume("Region Boundary Image",VoxelType.FLOAT));
		
		outputParams.setName("brain segmentation priors");
		outputParams.setLabel("brain segmentation priors");
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
		
		// clean-up
		warp1.finalize();
		warp1 = null;
		
		int nobj = atlas.getNumber();
		float[] proba = new float[nobj];
		int nax = atlas.getShapeDim()[0];
		int nay = atlas.getShapeDim()[1];
		int naz = atlas.getShapeDim()[2];
		float[] region = new float[nax*nay*naz];
		byte[] label = new byte[nax*nay*naz];
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			int xyz = x + nax*y + nax*nay*z;
			for (int nb=0;nb<nobj;nb++) {
				proba[nb] = atlas.getShape(nb)[xyz];
			}
			byte[] idx = Numerics.argmax(proba, 2);	
			region[xyz] = proba[idx[0]]/Numerics.max(0.001f, proba[idx[0]]+proba[idx[1]]);
			label[xyz] = atlas.getLabels()[idx[0]];
		}
		// warp into subject space
		float[][][] subjectregion = new float[nx][ny][nz];
		byte[][][] subjectlabel = new byte[nx][ny][nz];
		float[] XA = new float[3];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			atlas.imageToShapeCoordinates(XA, x, y, z);
			
			subjectregion[x][y][z] = ImageInterpolation.linearInterpolation(region, 0.0f, XA[0],XA[1],XA[2],nax,nay,naz);
			subjectlabel[x][y][z] = ImageInterpolation.nearestNeighborInterpolation(label, (byte)0, XA[0],XA[1],XA[2],nax,nay,naz);
		}
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
			
		String imgname = in1Img.getName();
		
		ImageDataByte segData = new ImageDataByte(subjectlabel);	
		subjectlabel = null;
		segData.setHeader(in1Img.getHeader());
		segData.setName(imgname+"_lbl");
		labelImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("labels");
		
		ImageDataFloat regionData = new ImageDataFloat(subjectregion);	
		subjectregion = null;
		regionData.setHeader(in1Img.getHeader());
		regionData.setName(imgname+"_reg");
		regionImage.setValue(regionData);
		regionData = null;
		BasicInfo.displayMessage(".. boundaries");

		return;
	}

}
