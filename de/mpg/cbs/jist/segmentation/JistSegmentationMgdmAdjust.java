package de.mpg.cbs.jist.segmentation;

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
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
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
public class JistSegmentationMgdmAdjust extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume labelImage;
	private ParamVolume mgdmImage;
	private ParamVolume correctImage;
	
	private ParamVolume memsImage;
	private ParamVolume lblsImage;
	
	private	ParamFloat 	balloonParam;
	private	ParamFloat 	edgeParam;
	private ParamFloat 	curvParam;
	private ParamFloat 	scaleParam;
	
	private ParamInteger 	iterationParam;
	private ParamFloat 	changeParam;
	
	private ParamOption 	topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamVolume newlabelImage;
	private ParamVolume newmgdmImage;
	private ParamVolume newmemsImage;
	private ParamVolume newlblsImage;
	private ParamVolume maskImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		imageParams=new ParamCollection("Images");
		imageParams.add(correctImage = new ParamVolume("New Segmentation Image"));
		
		imageParams.add(labelImage = new ParamVolume("Initial Segmentation Image"));
		imageParams.add(mgdmImage = new ParamVolume("Levelset Boundary Image"));
		
		imageParams.add(memsImage = new ParamVolume("Maximum Memberships (4D)"));
		imageParams.add(lblsImage = new ParamVolume("Maximum Labels (4D)"));
		memsImage.setMandatory(false);
		lblsImage.setMandatory(false);

		inputParams.add(imageParams);
		
		mainParams=new ParamCollection("Parameters");
		mainParams.add(balloonParam = new ParamFloat("Adjustment weight", -1E10f, 1E10f, 0.5f));
		mainParams.add(curvParam = new ParamFloat("Curvature weight", -1E10f, 1E10f, 0.2f));
		
		mainParams.add(scaleParam = new ParamFloat("Membership scale (mm)", -1E10f, 1E10f, 2.0f));
		
		mainParams.add(iterationParam = new ParamInteger("Max iterations", 0, 100000, 500));
		mainParams.add(changeParam = new ParamFloat("Min change", 0, 1, 0.001f));
			
		mainParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.add(mainParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Adjust Mgdm Segmentation");
		inputParams.setName("AdjustMgdmSegmentation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Adjusts a MGDM-based segmentation to a (manually or automatically) edited segmentation");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(newlabelImage = new ParamVolume("Adjusted Segmentation",VoxelType.BYTE));
		outputParams.add(newmgdmImage = new ParamVolume("Adjusted Boundaries",VoxelType.FLOAT));
		
		outputParams.add(newmemsImage = new ParamVolume("Adjusted Maximum Memberships (4D)",VoxelType.FLOAT,-1,-1,-1,-1));
		outputParams.add(newlblsImage = new ParamVolume("Adjusted Maximum Labels (4D)",VoxelType.BYTE,-1,-1,-1,-1));
		newmemsImage.setMandatory(false);
		newlblsImage.setMandatory(false);

		outputParams.add(maskImage = new ParamVolume("Computation mask",VoxelType.BYTE));
		maskImage.setMandatory(false);
		
		outputParams.setName("adjusted images");
		outputParams.setLabel("adjusted images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataByte labelImg = new ImageDataByte(labelImage.getImageData());
		int nx = labelImg.getRows();
		int ny = labelImg.getCols();
		int nz = labelImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = labelImg.getHeader().getDimResolutions()[0];
		float ry = labelImg.getHeader().getDimResolutions()[1];
		float rz = labelImg.getHeader().getDimResolutions()[2];
		byte[][][] bytebuffer = labelImg.toArray3d();
		byte[] segment = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			segment[xyz] = bytebuffer[x][y][z];
		}
		bytebuffer = null;
		
		ImageDataFloat mgdmImg = new ImageDataFloat(mgdmImage.getImageData());
		float[][][] buffer = mgdmImg.toArray3d();
		float[] mgdm = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			mgdm[xyz] = buffer[x][y][z];
		}
		buffer = null;
		
		ImageDataByte correctImg = new ImageDataByte(correctImage.getImageData());
		bytebuffer = correctImg.toArray3d();
		byte[] correct = new byte[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			correct[xyz] = bytebuffer[x][y][z];
		}
		bytebuffer = null;

		// adjust the posteriors
		int nlb = 0;
		byte[][] labels = null;
		float[][] mems = null;
		if (memsImage.getImageData()!=null && lblsImage.getImageData()!=null) {
			ImageDataByte lbImg = new ImageDataByte(lblsImage.getImageData());
			nlb = lbImg.getComponents();
			if (nlb==1) {
				bytebuffer = lbImg.toArray3d();
				lbImg = null;
				labels = new byte[nlb][nxyz];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					labels[0][xyz] = bytebuffer[x][y][z];
				}
				bytebuffer = null;
			} else {
				byte[][][][] bytebuffer4d = lbImg.toArray4d();
				lbImg = null;
				labels = new byte[nlb][nxyz];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					for (int n=0;n<nlb;n++) {
						labels[n][xyz] = bytebuffer4d[x][y][z][n];
					}
				}
				bytebuffer4d = null;
			}
			ImageDataFloat memImg = new ImageDataFloat(memsImage.getImageData());
			if (nlb==1) {
				buffer = memImg.toArray3d();
				memImg = null;
				mems = new float[nlb][nxyz];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					mems[0][xyz] = buffer[x][y][z];
				}
				buffer = null;
			} else {
				float[][][][] buffer4d = memImg.toArray4d();
				memImg = null;
				mems = new float[nlb][nxyz];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					for (int n=0;n<nlb;n++) {
						mems[n][xyz] = buffer4d[x][y][z][n];
					}
				}
				buffer4d = null;
			}
		}
		
		// change the membership values to match the corrected segmentation
		float scale = scaleParam.getValue().floatValue()/rx;
		if (mems!=null && labels!=null) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (segment[xyz]!=correct[xyz]) {
					// 1. find location of corrected label
					int ncorr=-1;
					float memcorr = 0;
					for (int n=0;n<nlb;n++) {
						if (labels[n][xyz]==correct[xyz]) {
							ncorr = n;
							memcorr = mems[n][xyz];
						}
					}
					
					// 2. move all the higher label values down
					if (ncorr>-1) {
						for (int n=ncorr;n>0;n--) {
							labels[n][xyz] = labels[n-1][xyz];
							mems[n][xyz] = mems[n-1][xyz];
						}
					} else {
						for (int n=nlb-1;n>0;n--) {
							labels[n][xyz] = labels[n-1][xyz];
							mems[n][xyz] = mems[n-1][xyz];
						}
					}
					
					// 3. artificially increase value to be above the rest
					labels[0][xyz] = correct[xyz];
					mems[0][xyz] = Numerics.max(memcorr, (0.5f+2.0f*mgdm[xyz]/scale)/(1.0f+2.0f*mgdm[xyz]/scale));
					if (nlb>1) mems[0][xyz] = Numerics.max(mems[0][xyz],Numerics.min(1.0f, mems[1][xyz]+0.1f));
				}
			}
		} else {
			// compute a basic label, membership pair from segmentation	
			labels = new byte[1][nxyz];
			mems = new float[1][nxyz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				labels[0][xyz] = correct[xyz];
				mems[0][xyz] = (0.5f+2.0f*mgdm[xyz]/scale)/(1.0f+2.0f*mgdm[xyz]/scale);
			}
		}
			
		// MGDM evolution		
		MgdmSegmentationCorrection refine = new MgdmSegmentationCorrection(correct, nx, ny, nz, rx, ry, rz,
																			segment, mgdm, 3,
																			labels, mems, nlb,
																			balloonParam.getValue().floatValue(), 
																			curvParam.getValue().floatValue(), 
																			topologyParam.getValue());
		
		if (iterationParam.getValue().intValue()>0) {
			BasicInfo.displayMessage("level set segmentation...\n");
			refine.evolveNarrowBand(iterationParam.getValue().intValue(),changeParam.getValue().floatValue());
		}
		
		byte[][][] seg = new byte[nx][ny][nz];
		float[] img = refine.exportMap(refine.labelSegmentation());
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			seg[x][y][z] = (byte)img[xyz];
		}
		img = null;
		ImageDataByte segData = new ImageDataByte(seg);	
		seg = null;
		segData.setHeader(mgdmImg.getHeader());
		segData.setName(labelImg.getName()+"_correct");
		newlabelImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("segmentation");
		
		float[][][] lvl = new float[nx][ny][nz];
		float[] fcn = refine.getFunctions()[0];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			lvl[x][y][z] = fcn[xyz];
		}
		ImageDataFloat lvlData = new ImageDataFloat(lvl);	
		lvl = null;
		lvlData.setHeader(mgdmImg.getHeader());
		lvlData.setName(mgdmImg.getName()+"_correct");
		newmgdmImage.setValue(lvlData);
		lvlData = null;
		BasicInfo.displayMessage(".. boundaries");
		
		mems = refine.exportBestGainFunctions();
		float[][][][] buffer4d = new float[nx][ny][nz][nlb];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			for (int n=0;n<nlb;n++) {
				buffer4d[x][y][z][n] = mems[n][xyz];
			}
		}
		mems = null;
		ImageDataFloat bufferData = new ImageDataFloat(buffer4d);
		buffer4d = null;
		bufferData.setHeader(mgdmImg.getHeader());
		bufferData.setName(labelImg.getName()+"_cxmems");
		newmemsImage.setValue(bufferData);
		bufferData = null;
		BasicInfo.displayMessage(".. memberships(4d)");
	
		float[][] lbls = new float[nlb][];
		for (int n=0;n<nlb;n++) {
			lbls[n] = refine.exportBestGainLabels(n);
		}
		byte[][][][] byte4d = new byte[nx][ny][nz][nlb];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			for (int n=0;n<nlb;n++) {
				byte4d[x][y][z][n] = (byte)lbls[n][xyz];
			}
		}
		lbls = null;
		ImageDataByte byteData = new ImageDataByte(byte4d);
		byte4d = null;
		byteData.setHeader(mgdmImg.getHeader());
		byteData.setName(labelImg.getName()+"_cxlbls");
		newlblsImage.setValue(byteData);
		byteData = null;			
		BasicInfo.displayMessage(".. labels(4d)");

		byte[][][] msk = new byte[nx][ny][nz];
		boolean[] mask = refine.getMask();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (mask[xyz])	msk[x][y][z] = 1;
			else msk[x][y][z] = 0;
		}
		mask = null;
		ImageDataByte mskData = new ImageDataByte(msk);	
		msk = null;
		mskData.setHeader(mgdmImg.getHeader());
		mskData.setName(labelImg.getName()+"_cxmask");
		maskImage.setValue(mskData);
		mskData = null;
		BasicInfo.displayMessage("mask");
		
	}

}
