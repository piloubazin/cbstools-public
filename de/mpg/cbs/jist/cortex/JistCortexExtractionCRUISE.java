package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistCortexExtractionCRUISE extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume initImage;
	private ParamVolume wmImage;
	private ParamVolume gmImage;
	private ParamVolume csfImage;
	
	//private ParamFloat 	edgeParam;
	private	ParamFloat 	balloonParam;
	private ParamFloat 	curvParam;
	//private ParamFloat 	divParam;
	private ParamInteger 	iterationParam;
	//private ParamOption 	gwImageParam;
	//private ParamOption 	cgImageParam;
	private ParamOption 	topologyParam;
	private ParamBoolean	normalizeParam;
	
	//private static final String[] opts = {"Iso image", "T1 map", "Probabilities", "Iso+Proba", "T1map+Proba", "All"};
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume cortexImage;
	private ParamVolume gwbImage;
	private ParamVolume cgbImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(initImage = new ParamVolume("Initial WM Segmentation Image"));
		imageParams.add(wmImage = new ParamVolume("Filled WM Probability Image"));
		imageParams.add(gmImage = new ParamVolume("GM Probability Image"));
		imageParams.add(csfImage = new ParamVolume("CSF and BG Probability Image"));
		
		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		//mainParams.add(edgeParam = new ParamFloat("Edge weight", -1E10f, 1E10f, 0.1f25));
		mainParams.add(balloonParam = new ParamFloat("Data weight (balloon force)", -1E10f, 1E10f, 0.9f));
		mainParams.add(curvParam = new ParamFloat("Regularization weight (curvature)", -1E10f, 1E10f, 0.1f));
		//mainParams.add(divParam = new ParamFloat("Divergence weight", -1E10f, 1E10f, 0.0f5));
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
			
		//mainParams.add(gwImageParam = new ParamOption("Image used for WM / GM", opts));
		//gwImageParam.setValue("Iso+Proba");

		//mainParams.add(cgImageParam = new ParamOption("Image used for GM / CSF", opts));
		//cgImageParam.setValue("Iso+Proba");

		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		mainParams.add(normalizeParam = new ParamBoolean("Normalize the probabilities", true));

		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing.prev");
		inputParams.setLabel("CRUISE Cortex Extraction");
		inputParams.setName("CruiseCortexExtraction");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Xiao Han", "","http://www.iacl.ece.jhu.edu/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new Citation("X. Han, D.L. Pham, D. Tosun, M.E. Rettmann, C. Xu, and J. L. Prince, "
								+"CRUISE: Cortical Reconstruction Using Implicit Surface Evolution, "
								+"NeuroImage, vol. 23, pp. 997--1012, 2004."));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Segments the cortex from a whole brain segmented data set with the CRUISE method.");
		
		info.setVersion("2.0f");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(cortexImage = new ParamVolume("Cortex Mask",VoxelType.UBYTE));
		outputParams.add(gwbImage = new ParamVolume("WM-GM Levelset",VoxelType.FLOAT));
		outputParams.add(cgbImage = new ParamVolume("GM-CSF Levelset",VoxelType.FLOAT));
		
		outputParams.setName("cortical images");
		outputParams.setLabel("cortical images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat wmImg = new ImageDataFloat(wmImage.getImageData());
		int nx = wmImg.getRows();
		int ny = wmImg.getCols();
		int nz = wmImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = wmImg.getHeader().getDimResolutions()[0];
		float ry = wmImg.getHeader().getDimResolutions()[1];
		float rz = wmImg.getHeader().getDimResolutions()[2];
		float[][][] buffer = wmImg.toArray3d();
		float[] wm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			wm[xyz] = buffer[x][y][z];
		}
		wmImg = null;
		buffer = null;
		
		ImageDataFloat gmImg = new ImageDataFloat(gmImage.getImageData());
		buffer = gmImg.toArray3d();
		float[] gm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			gm[xyz] = buffer[x][y][z];
		}
		gmImg = null;
		buffer = null;
				
		ImageDataFloat csfImg = new ImageDataFloat(csfImage.getImageData());
		buffer = csfImg.toArray3d();
		float[] csf = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			csf[xyz] = buffer[x][y][z];
		}
		csfImg = null;
		buffer = null;
		
		ImageDataUByte initImg = new ImageDataUByte(initImage.getImageData());
		byte[][][] bytebuffer = initImg.toArray3d();
		byte[] init = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			init[xyz] = bytebuffer[x][y][z];
		}
		initImg = null;
		bytebuffer = null;
				
		// normalize the probabilities into memberships
		if (normalizeParam.getValue().booleanValue()) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				float sum = wm[xyz]+gm[xyz]+csf[xyz];
				if (sum>0.0001f) {
					wm[xyz] /= sum;
					gm[xyz] /= sum;
					csf[xyz] /= sum;
				}
			}
		}		
				
		// main algorithm
		boolean useProbas = true;
		boolean gwbflag = true;
		
		// levelset evolution: wm /gm		
		Mp2rageCortexGdm gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
														wm, gm, csf, null,
														null, null,
														0.0f, balloonParam.getValue().floatValue(), 
														curvParam.getValue().floatValue(), 0.0f,
														topologyParam.getValue(),null,useProbas, gwbflag);
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		}
		
		// output
		byte[][][] seg = new byte[nx][ny][nz];
		float[][][] gwb = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = gdm.getSegmentation()[xyz];
			gwb[x][y][z] = gdm.getLevelSet()[xyz];
		}
		
		// second step

		// use the wm result as the filler mask => guarantees the boundaries never cross
		init = gdm.getSegmentation();

		boolean cgbflag = !gwbflag;
		
		// levelset evolution: gm / csf
		gdm = new Mp2rageCortexGdm(init, null, nx, ny, nz, rx, ry, rz,
									wm, gm, csf, null,
									null, null,
									0.0f, balloonParam.getValue().floatValue(), 
									curvParam.getValue().floatValue(), 0.0f,
									topologyParam.getValue(),null,useProbas, cgbflag);
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			
			gdm.evolveNarrowBand(iterationParam.getValue().intValue(), 0.001f);
		}
		
		// output
		String imgname = gmImage.getImageData().getName();
		
		float[][][] cgb = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = (byte)(gdm.getSegmentation()[xyz]+seg[x][y][z]);
			cgb[x][y][z] = gdm.getLevelSet()[xyz];
		}
		
		ImageDataUByte cortexData = new ImageDataUByte(seg);		
		cortexData.setHeader(wmImage.getImageData().getHeader());
		cortexData.setName(imgname+"_cortex");
		cortexImage.setValue(cortexData);
		cortexData = null;
		seg = null;
		
		ImageDataFloat gwbData = new ImageDataFloat(gwb);		
		gwbData.setHeader(wmImage.getImageData().getHeader());
		gwbData.setName(imgname+"_gwb");
		gwbImage.setValue(gwbData);
		gwbData = null;
		gwb = null;
		
		ImageDataFloat cgbData = new ImageDataFloat(cgb);		
		cgbData.setHeader(wmImage.getImageData().getHeader());
		cgbData.setName(imgname+"_cgb");
		cgbImage.setValue(cgbData);
		cgbData = null;
		cgb = null;
		
	}


}
