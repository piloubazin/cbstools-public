package de.mpg.cbs.jist.registration;

import java.io.*;

import Jama.Matrix;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamMatrix;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipavWrapper;

import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.MatrixHolder;
import gov.nih.mipav.model.structures.TransMatrix;
import gov.nih.mipav.model.algorithms.AlgorithmWSinc;
import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.file.FileNIFTI;
import gov.nih.mipav.view.Preferences;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.libraries.*;

import javax.vecmath.Point3i;

/*
 * Copy Header
 * @author Pilou Bazin
 *
 */
public class JistRegistrationScannerTransform extends ProcessingAlgorithm{
	private ParamVolume srcVol;
	private ParamVolume trgVol;
	private ParamBoolean notransBoolean;
	private ParamBoolean scannerBoolean;
	//private ParamBoolean offsetBoolean;
	//private ParamPointInteger offsrcPoint;
	//private ParamPointInteger offtrgPoint;
	
	private ParamVolume alignedSrcVol;
	private ParamVolume alignedTrgVol;
	public ParamMatrix outSrcTransform;
	public ParamMatrix outTrgTransform;

	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(srcVol = new ParamVolume("Source Image", null,-1,-1,-1,-1));
		inputParams.add(trgVol = new ParamVolume("Target Image", null,-1,-1,-1,-1));
		inputParams.add(notransBoolean = new ParamBoolean("Ignore translation",false));
		inputParams.add(scannerBoolean = new ParamBoolean("Output in scanner space",true));
		//inputParams.add(offsetBoolean = new ParamBoolean("Compensate for offset",false));
		//inputParams.add(offsrcPoint = new ParamPointInteger("Offset source point", new Point3i(0,0,0)));
		//inputParams.add(offtrgPoint = new ParamPointInteger("Offset target point", new Point3i(0,0,0)));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration.devel");
		inputParams.setLabel("Scanner Transform");
		inputParams.setName("ScannerTransform");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Transform two images to match via scanner coordinates");
		info.setVersion("3.0.3");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(alignedSrcVol=new ParamVolume("Aligned Source", null,-1,-1,-1,-1));
		outputParams.add(alignedTrgVol=new ParamVolume("Aligned Target", null,-1,-1,-1,-1));
		outputParams.add(outSrcTransform=new ParamMatrix("Source Transformation",4,4));
		outputParams.add(outTrgTransform=new ParamMatrix("Target Transformation",4,4));
	}


	protected void execute(CalculationMonitor monitor) {
		
		// must force to uncheck the 'flip NIFTI orientation' before running this!!
		String flip = Preferences.getProperty(Preferences.PREF_FLIP_NIFTI_READ);
		Preferences.setProperty(Preferences.PREF_FLIP_NIFTI_READ,"false");
		
		ImageDataMipavWrapper src = (ImageDataMipavWrapper)srcVol.getImageData();
		TransMatrix src2scan = src.getModelImageDirect().getMatrixHolder().getNIFTICompositeMatrices()[0];
		// pull the translation info directly from NIFTI header
		String srcfile = src.getModelImageDirect().getImageFileName();
		String srcdir = src.getModelImageDirect().getImageDirectory();
		boolean srccompressed = (srcfile.endsWith("zip") || srcfile.endsWith("gz") || srcfile.endsWith("bz2"));
		NiftiInterface niftisrc = new NiftiInterface(srcfile, srcdir);
		try { 
			niftisrc.readHeader(srcfile, srcdir, srccompressed);
		} catch (IOException e) {
			Interface.displayError(e.getMessage());
		}
		float[] osrc = new float[3];
		osrc[0] = niftisrc.qoffset_x;
		osrc[1] = niftisrc.qoffset_y;
		osrc[2] = niftisrc.qoffset_z;
		niftisrc.finalize(); niftisrc = null;
		
		ImageDataMipavWrapper trg = (ImageDataMipavWrapper)trgVol.getImageData();
		TransMatrix trg2scan = trg.getModelImageDirect().getMatrixHolder().getNIFTICompositeMatrices()[0];
		// pull the translation info directly from NIFTI header
		String trgfile = trg.getModelImageDirect().getImageFileName();
		String trgdir = trg.getModelImageDirect().getImageDirectory();
		boolean trgcompressed = (trgfile.endsWith("zip") || trgfile.endsWith("gz") || trgfile.endsWith("bz2"));
		NiftiInterface niftitrg = new NiftiInterface(trgfile, trgdir);
		try { 
			niftitrg.readHeader(trgfile, trgdir, trgcompressed);
		} catch (IOException e) {
			Interface.displayError(e.getMessage());
		}
		float[] otrg = new float[3];
		otrg[0] = niftitrg.qoffset_x;
		otrg[1] = niftitrg.qoffset_y;
		otrg[2] = niftitrg.qoffset_z;
		niftitrg.finalize(); niftitrg = null;
		
		// build the two transform matrices
		Matrix srcmat = Matrix.identity(4,4);
		Matrix trgmat = Matrix.identity(4,4);
	
		// assume the rotation are good
		for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
			srcmat.set(i,j, src2scan.get(i,j));
			trgmat.set(i,j, trg2scan.get(i,j));
		}
		/*
		// rotation part
		for (int j=0;j<4;j++) {
			srcmat.set(0,j, -src2scan.get(0,j));
			srcmat.set(1,j, -src2scan.get(1,j));
			srcmat.set(2,j, src2scan.get(2,j));
			
			trgmat.set(0,j, -trg2scan.get(0,j));
			trgmat.set(1,j, -trg2scan.get(1,j));
			trgmat.set(2,j, trg2scan.get(2,j));
		}		
		
		// check rotation determinant
		Matrix srcR = srcmat.getMatrix(0,2, 0,2);
		srcR = mat33_polar(srcR);
		if (srcR.det() < 0) {
			// qfac = -1
			for (int j=0;j<3;j++) srcmat.set(j,2, -srcmat.get(j,2));
        }
		Matrix trgR = trgmat.getMatrix(0,2, 0,2);
		trgR = mat33_polar(trgR);
		if (trgR.det() < 0) {
			// qfac = -1
			for (int j=0;j<3;j++) trgmat.set(j,2, -trgmat.get(j,2));
        }
        */
        /* replaced by pulling the q form directly from the NIFTI header
        // load translation from origin (weird conversion of data)
        float[] srco = new float[3];
        float[] trgo = new float[3];
        for (int j=0;j<3;j++) {
        	srco[j] = src.getModelImageDirect().getFileInfo()[0].getOrigin()[j];
        	trgo[j] = trg.getModelImageDirect().getFileInfo()[0].getOrigin()[j];
        }
		// special cases for translation
		/*
		int[] srcOrientation = getAxisOrientation(src2scan);
		int[] trgOrientation = getAxisOrientation(trg2scan);
		*/
		/*
		int[] srcOrientation = src.getModelImageDirect().getFileInfo()[0].getAxisOrientation();
		int[] trgOrientation = trg.getModelImageDirect().getFileInfo()[0].getAxisOrientation();
		for (int j = 0; j < 3; j++) {
			if (srcOrientation[j] == FileInfoBase.ORI_L2R_TYPE) {
				srcmat.set(0,3, -Math.abs(srco[j]));
			} else if (srcOrientation[j] == FileInfoBase.ORI_R2L_TYPE) {
				srcmat.set(0,3, Math.abs(srco[j]));
			} else if (srcOrientation[j] == FileInfoBase.ORI_P2A_TYPE) {
				srcmat.set(1,3, -Math.abs(srco[j]));
			} else if (srcOrientation[j] == FileInfoBase.ORI_A2P_TYPE) {
				srcmat.set(1,3, Math.abs(srco[j]));
			} else if (srcOrientation[j] == FileInfoBase.ORI_I2S_TYPE) {
				srcmat.set(2,3, -Math.abs(srco[j]));
			} else if (srcOrientation[j] == FileInfoBase.ORI_S2I_TYPE) {
				srcmat.set(2,3, Math.abs(srco[j]));
			}
			if (trgOrientation[j] == FileInfoBase.ORI_L2R_TYPE) {
				 trgmat.set(0,3, -Math.abs(trgo[j]));
			} else if (trgOrientation[j] == FileInfoBase.ORI_R2L_TYPE) {
				trgmat.set(0,3, Math.abs(trgo[j]));
			} else if (trgOrientation[j] == FileInfoBase.ORI_P2A_TYPE) {
				trgmat.set(1,3, -Math.abs(trgo[j]));
			} else if (trgOrientation[j] == FileInfoBase.ORI_A2P_TYPE) {
				trgmat.set(1,3, Math.abs(trgo[j]));
			} else if (trgOrientation[j] == FileInfoBase.ORI_I2S_TYPE) {
				trgmat.set(2,3, -Math.abs(trgo[j]));
			} else if (trgOrientation[j] == FileInfoBase.ORI_S2I_TYPE) {
				trgmat.set(2,3, Math.abs(trgo[j]));
			}
		}	
        // times -1,-1,+1 ?? it depends :(
        // not a clean systemeatic approach !!!
        srcmat.set(0,3, -srcmat.get(0,3));
        srcmat.set(1,3, -srcmat.get(1,3));
        trgmat.set(0,3, -trgmat.get(0,3));
        trgmat.set(1,3, -trgmat.get(1,3));
        trgmat.set(2,3, -trgmat.get(2,3));
        */
        for (int j=0;j<3;j++) {
        	srcmat.set(j,3, osrc[j]);
        	trgmat.set(j,3, otrg[j]);
        }
        // x,y are turned to negative
        srcmat.set(0,3, -srcmat.get(0,3));
        srcmat.set(1,3, -srcmat.get(1,3));
        trgmat.set(0,3, -trgmat.get(0,3));
        trgmat.set(1,3, -trgmat.get(1,3));
        
        /*
		// offset: subtract to equate the points in translation
		if (offsetBoolean.getValue()) {
			srcmat.set(0,3, srcmat.get(0,3)-offsrcPoint.getValue().x);
			srcmat.set(1,3, srcmat.get(1,3)-offsrcPoint.getValue().y);
			srcmat.set(2,3, srcmat.get(2,3)-offsrcPoint.getValue().z);
			
			trgmat.set(0,3, trgmat.get(0,3)-offtrgPoint.getValue().x);
			trgmat.set(1,3, trgmat.get(1,3)-offtrgPoint.getValue().y);
			trgmat.set(2,3, trgmat.get(2,3)-offtrgPoint.getValue().z);
		}
		*/
		
		// auxiliary matrices
		Matrix srcrot = Matrix.identity(4,4);
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) srcrot.set(i,j, srcmat.get(i,j));
		Matrix srcinv = srcmat.inverse();
		Matrix trgrot = Matrix.identity(4,4);
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) trgrot.set(i,j, trgmat.get(i,j));
		Matrix trginv = trgmat.inverse();
		
		// combine matrices (different options)
		Matrix src2trg = Matrix.identity(4,4);
		Matrix trg2src = Matrix.identity(4,4);
		if (notransBoolean.getValue()) {
			if (scannerBoolean.getValue()) {
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) src2trg.set(i,j, srcrot.get(i,j));
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) trg2src.set(i,j, trgrot.get(i,j));
			} else {
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
					src2trg.set(i,j, 0.0f);
					// this uses the inverse rotation (R^-1 = R^T)
					for (int k=0;k<4;k++) src2trg.set(i,j, src2trg.get(i,j)+trgrot.get(k,i)*srcrot.get(k,j));
					trg2src.set(i,j, 0.0f);
					// this uses the inverse rotation (R^-1 = R^T)
					for (int k=0;k<4;k++) trg2src.set(i,j, trg2src.get(i,j)+srcrot.get(k,i)*trgrot.get(k,j));
				}
			}
		} else {
			if (scannerBoolean.getValue()) {
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) src2trg.set(i,j, srcmat.get(i,j));
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) trg2src.set(i,j, trgmat.get(i,j));
			} else {
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
					src2trg.set(i,j, 0.0f);
					for (int k=0;k<4;k++) src2trg.set(i,j, src2trg.get(i,j)+trginv.get(i,k)*srcmat.get(k,j));
				}
				for (int i=0;i<4;i++) for (int j=0;j<4;j++) {
					trg2src.set(i,j, 0.0f);
					for (int k=0;k<4;k++) trg2src.set(i,j, trg2src.get(i,j)+srcinv.get(i,k)*trgmat.get(k,j));
				}
			}
		}
		Matrix invsrc2trg = src2trg.inverse();
		Matrix invtrg2src = trg2src.inverse();
				
		// transform images (in scanner space, make the field of view large enough for both)
		ImageDataFloat srcImg = new ImageDataFloat(src);
		ImageDataFloat trgImg = new ImageDataFloat(trg);
		
		int nsx = srcImg.getRows();
		int nsy = srcImg.getCols();
		int nsz = srcImg.getSlices();
		int nst = srcImg.getComponents();
		
		int ntx = trgImg.getRows();
		int nty = trgImg.getCols();
		int ntz = trgImg.getSlices();
		int ntt = trgImg.getComponents();
		
		float[][][] source;
		float[][][] target;
		if (nst>1) {
			float[][][][] tmp = srcImg.toArray4d();
			source = new float[nsx][nsy][nsz];
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) source[x][y][z] = tmp[x][y][z][0];
		} else {
			source = srcImg.toArray3d();
		}
		if (ntt>1) {
			float[][][][] tmp = trgImg.toArray4d();
			target = new float[ntx][nty][ntz];
			for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) target[x][y][z] = tmp[x][y][z][0];
		} else {
			target = trgImg.toArray3d();
		}
		
		int x0 = 0, y0 = 0, z0 = 0;
		if (scannerBoolean.getValue()) {
			// 1. find image extents from deformed ones	
			Matrix spt = new Matrix(4,1, 0.0f); spt.set(3,0, 1);
			Matrix tpt = new Matrix(4,1, 0.0f); tpt.set(3,0, 1);
			Matrix pt = null;
			
			// source bounding box
			pt = src2trg.times(spt);
			x0 = Numerics.floor(pt.get(0,0)); int xN = Numerics.ceil(pt.get(0,0));
			y0 = Numerics.floor(pt.get(1,0)); int yN = Numerics.ceil(pt.get(1,0));
			z0 = Numerics.floor(pt.get(2,0)); int zN = Numerics.ceil(pt.get(2,0));

			spt.set(0,0, nsx); spt.set(1,0, 0); spt.set(2,0, 0); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, 0); spt.set(1,0, nsy); spt.set(2,0, 0); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, 0); spt.set(1,0, 0); spt.set(2,0, nsz); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, nsx); spt.set(1,0, nsy); spt.set(2,0, 0); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, 0); spt.set(1,0, nsy); spt.set(2,0, nsz); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, nsx); spt.set(1,0, 0); spt.set(2,0, nsz); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			spt.set(0,0, nsx); spt.set(1,0, nsy); spt.set(2,0, nsz); spt.set(3,0, 1);
			pt = src2trg.times(spt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));

			// target bounding box
			tpt.set(0,0, 0); tpt.set(1,0, 0); tpt.set(2,0, 0); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, ntx); tpt.set(1,0, 0); tpt.set(2,0, 0); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, 0); tpt.set(1,0, nty); tpt.set(2,0, 0); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, 0); tpt.set(1,0, 0); tpt.set(2,0, ntz); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, ntx); tpt.set(1,0, nty); tpt.set(2,0, 0); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, 0); tpt.set(1,0, nty); tpt.set(2,0, ntz); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, ntx); tpt.set(1,0, 0); tpt.set(2,0, ntz); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			tpt.set(0,0, ntx); tpt.set(1,0, nty); tpt.set(2,0, ntz); tpt.set(3,0, 1);
			pt = trg2src.times(tpt);
			x0 = Numerics.min(x0, Numerics.floor(pt.get(0,0))); xN = Numerics.max(xN, Numerics.ceil(pt.get(0,0)));
			y0 = Numerics.min(y0, Numerics.floor(pt.get(1,0))); yN = Numerics.max(yN, Numerics.ceil(pt.get(1,0)));
			z0 = Numerics.min(z0, Numerics.floor(pt.get(2,0))); zN = Numerics.max(zN, Numerics.ceil(pt.get(2,0)));
			
			// find extents
			System.out.println("Scanner coordinates: ("+x0+", "+y0+", "+z0+") - ("+xN+", "+yN+", "+zN+")");
			
			int nx = xN-x0+1;
			int ny = yN-y0+1;
			int nz = zN-z0+1;
			
			float[][][] scansrc = new float[nx][ny][nz];
			float[][][] scantrg = new float[nx][ny][nz];
			
			// 2. map everything
			for (int x=x0;x<=xN;x++) for (int y=y0;y<=yN;y++) for (int z=z0;z<=zN;z++) {
				pt.set(0,0, x);
				pt.set(1,0, y);
				pt.set(2,0, z);
				pt.set(3,0, 1);
				// use inverse transforms here
				spt = invsrc2trg.times(pt);
				tpt = invtrg2src.times(pt);
				
				scansrc[x-x0][y-y0][z-z0] = ImageInterpolation.linearInterpolation(source, 0.0f, (float)spt.get(0,0), (float)spt.get(1,0), (float)spt.get(2,0), nsx, nsy, nsz);
				scantrg[x-x0][y-y0][z-z0] = ImageInterpolation.linearInterpolation(target, 0.0f, (float)tpt.get(0,0), (float)tpt.get(1,0), (float)tpt.get(2,0), ntx, nty, ntz);
			}
			
			ImageDataFloat srcData = new ImageDataFloat(scansrc);		
			//srcData.setHeader(srcImg.getHeader());
			srcData.setName(srcImg.getName()+"_2scan");
			alignedSrcVol.setValue(srcData);
			
			ImageDataFloat trgData = new ImageDataFloat(scantrg);		
			//trgData.setHeader(trgImg.getHeader());
			trgData.setName(trgImg.getName()+"_2scan");
			alignedTrgVol.setValue(trgData);
				
		} else {
			float[][][] mapdsrc = new float[ntx][nty][ntz];
			float[][][] mapdtrg = new float[nsx][nsy][nsz];
			
			// 2. map everything
			Matrix pt = new Matrix(4,1, 0.0f);
			for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
				pt.set(0,0, x);
				pt.set(1,0, y);
				pt.set(2,0, z);
				pt.set(3,0, 1);
				// use inverse transforms here
				Matrix spt = invsrc2trg.times(pt);
				
				mapdsrc[x][y][z] = ImageInterpolation.linearInterpolation(source, 0.0f, (float)spt.get(0,0), (float)spt.get(1,0), (float)spt.get(2,0), nsx, nsy, nsz);
			}
			for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
				pt.set(0,0, x);
				pt.set(1,0, y);
				pt.set(2,0, z);
				pt.set(3,0, 1);
				// use inverse transforms here
				Matrix tpt = invtrg2src.times(pt);
				
				mapdtrg[x][y][z] = ImageInterpolation.linearInterpolation(target, 0.0f, (float)tpt.get(0,0), (float)tpt.get(1,0), (float)tpt.get(2,0), ntx, nty, ntz);
			}
			
			ImageDataFloat srcData = new ImageDataFloat(mapdsrc);		
			srcData.setHeader(trgImg.getHeader());
			srcData.setName(srcImg.getName()+"_2trg");
			alignedSrcVol.setValue(srcData);
			
			ImageDataFloat trgData = new ImageDataFloat(mapdtrg);		
			trgData.setHeader(srcImg.getHeader());
			trgData.setName(trgImg.getName()+"_2src");
			alignedTrgVol.setValue(trgData);
		}
		
		
		// output transfor matrices rescaled to remove resolution impact (for use with MIPAV's transform module)
		float[] sres = srcImg.getHeader().getDimResolutions();
		float[] tres = trgImg.getHeader().getDimResolutions();
		if (scannerBoolean.getValue()) {
			// modify the inverse, then re-invert
			
			// convolve the x0,y0,z0 vector before changing anything
			Matrix pt = new Matrix(4,1, 0.0f); pt.set(3,0, 1);
			Matrix spt = new Matrix(4,1, 0.0f); spt.set(3,0, 1);
			Matrix tpt = new Matrix(4,1, 0.0f); tpt.set(3,0, 1);
			
			pt.set(0,0, x0);
			pt.set(1,0, y0);
			pt.set(2,0, z0);
			pt.set(3,0, 0); // use only the rotation part
			
			// use inverse transforms here
			spt = invsrc2trg.times(pt);
			tpt = invtrg2src.times(pt);
			
			// update rotation
			for (int i=0;i<3;i++) for (int j=0;j<3;j++) invsrc2trg.set(i,j, invsrc2trg.get(i,j)*sres[i]);
			for (int i=0;i<3;i++) for (int j=0;j<3;j++) invtrg2src.set(i,j, invtrg2src.get(i,j)*tres[i]);
			// update translation
			for (int i=0;i<3;i++) invsrc2trg.set(i,3, invsrc2trg.get(i,3)*sres[i]+spt.get(i,0)*sres[i]);
			for (int i=0;i<3;i++) invtrg2src.set(i,3, invtrg2src.get(i,3)*tres[i]+tpt.get(i,0)*tres[i]);
			
			src2trg = invsrc2trg.inverse();
			trg2src = invtrg2src.inverse();
		} else {
			// rotation
			//for (int i=0;i<3;i++) for (int j=0;j<3;j++) invsrc2trg.set(i,j, invsrc2trg.get(i,j)/tres[i]*sres[j]);
			//for (int i=0;i<3;i++) for (int j=0;j<3;j++) invtrg2src.set(i,j, invtrg2src.get(i,j)/sres[i]*tres[j]);
			for (int i=0;i<3;i++) for (int j=0;j<3;j++) invsrc2trg.set(i,j, invsrc2trg.get(i,j)/tres[j]*sres[i]);
			for (int i=0;i<3;i++) for (int j=0;j<3;j++) invtrg2src.set(i,j, invtrg2src.get(i,j)/sres[j]*tres[i]);
			// translation
			for (int i=0;i<3;i++) invsrc2trg.set(i,3, invsrc2trg.get(i,3)*sres[i]);
			for (int i=0;i<3;i++) invtrg2src.set(i,3, invtrg2src.get(i,3)*tres[i]);
			
			src2trg = invsrc2trg.inverse();
			trg2src = invtrg2src.inverse();
		}		
		outSrcTransform.setValue(src2trg);
		outTrgTransform.setValue(trg2src);
		
		// reset preferences to previous setting
		Preferences.setProperty(Preferences.PREF_FLIP_NIFTI_READ, flip);
		
		return;
	}


	// from FileNIFTI.class
	private int[] getAxisOrientation(TransMatrix mat) {
        int[] axisOrientation = new int[3];
        //double[][] array;
        double xi, xj, xk, yi, yj, yk, zi, zj, zk, val, detQ, detP;
        Matrix P, Q, M;
        int i, j, k = 0, p, q, r, ibest, jbest, kbest, pbest, qbest, rbest;
        double vbest;

        /* load column vectors for each (i,j,k) direction from matrix */

        /*-- i axis --*/
        /*-- j axis --*/
        /*-- k axis --*/
        //array = mat.getMatrix(0, 2, 0, 2).getArray();
        xi = mat.get(0, 0);
        xj = mat.get(0, 1);
        xk = mat.get(0, 2);
        yi = mat.get(1, 0);
        yj = mat.get(1, 1);
        yk = mat.get(1, 2);
        zi = mat.get(2, 0);
        zj = mat.get(2, 1);
        zk = mat.get(2, 2);

        /* normalize column vectors to get unit vectors along each ijk-axis */

        /* normalize i axis */
        axisOrientation[0] = FileInfoBase.ORI_UNKNOWN_TYPE;
        val = Math.sqrt((xi * xi) + (yi * yi) + (zi * zi));

        if (val == 0.0f) {
            return axisOrientation; /* stupid input */
        }

        xi /= val;
        yi /= val;
        zi /= val;

        /* normalize j axis */

        val = Math.sqrt((xj * xj) + (yj * yj) + (zj * zj));

        if (val == 0.0f) {
            return axisOrientation; /* stupid input */
        }

        xj /= val;
        yj /= val;
        zj /= val;

        /* orthogonalize j axis to i axis, if needed */

        val = (xi * xj) + (yi * yj) + (zi * zj); /* dot product between i and j */

        if (Math.abs(val) > 1.e-4) {
            xj -= val * xi;
            yj -= val * yi;
            zj -= val * zi;
            val = Math.sqrt((xj * xj) + (yj * yj) + (zj * zj)); /* must renormalize */

            if (val == 0.0f) {
                return axisOrientation; /* j was parallel to i? */
            }

            xj /= val;
            yj /= val;
            zj /= val;
        }

        /* normalize k axis; if it is zero, make it the cross product i x j */

        val = Math.sqrt((xk * xk) + (yk * yk) + (zk * zk));

        if (val == 0.0f) {
            xk = (yi * zj) - (zi * yj);
            yk = (zi * xj) - (zj * xi);
            zk = (xi * yj) - (yi * xj);
        } else {
            xk /= val;
            yk /= val;
            zk /= val;
        }

        /* orthogonalize k to i */

        val = (xi * xk) + (yi * yk) + (zi * zk); /* dot product between i and k */

        if (Math.abs(val) > 1.e-4) {
            xk -= val * xi;
            yk -= val * yi;
            zk -= val * zi;
            val = Math.sqrt((xk * xk) + (yk * yk) + (zk * zk));

            if (val == 0.0f) {
                return axisOrientation; /* bad */
            }

            xk /= val;
            yk /= val;
            zk /= val;
        }

        /* orthogonalize k to j */

        val = (xj * xk) + (yj * yk) + (zj * zk); /* dot product between j and k */

        if (Math.abs(val) > 1.e-4) {
            xk -= val * xj;
            yk -= val * yj;
            zk -= val * zj;
            val = Math.sqrt((xk * xk) + (yk * yk) + (zk * zk));

            if (val == 0.0f) {
                return axisOrientation; /* bad */
            }

            xk /= val;
            yk /= val;
            zk /= val;
        }

        Q = new Matrix(3, 3);

        Q.set(0, 0, xi);
        Q.set(0, 1, xj);
        Q.set(0, 2, xk);
        Q.set(1, 0, yi);
        Q.set(1, 1, yj);
        Q.set(1, 2, yk);
        Q.set(2, 0, zi);
        Q.set(2, 1, zj);
        Q.set(2, 2, zk);

        /* at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes */

        detQ = Q.det();

        if (detQ == 0.0f) {
            Interface.displayError("detQ == 0.0f in getAxisOrientation");

            return axisOrientation;
        }

        P = new Matrix(3, 3);
        /* Build and test all possible +1/-1 coordinate permutation matrices P;
         * then find the P such that the rotation matrix M=PQ is closest to theidentity, in the sense of M having the
         * smallest total rotation angle. */

        /* Despite the formidable looking 6 nested loops, there are
         *only 3*3*3*2*2*2 = 216 passes, which will run very quickly. */

        vbest = -666.0f;
        ibest = pbest = qbest = rbest = 1;
        jbest = 2;
        kbest = 3;

        for (i = 1; i <= 3; i++) { /* i = column number to use for row #1 */

            for (j = 1; j <= 3; j++) { /* j = column number to use for row #2 */

                if (i == j) {
                    continue;
                }

                for (k = 1; k <= 3; k++) { /* k = column number to use for row #3 */

                    if ((i == k) || (j == k)) {
                        continue;
                    }

                    P.set(0, 0, 0.0f);
                    P.set(0, 1, 0.0f);
                    P.set(0, 2, 0.0f);
                    P.set(1, 0, 0.0f);
                    P.set(1, 1, 0.0f);
                    P.set(1, 2, 0.0f);
                    P.set(2, 0, 0.0f);
                    P.set(2, 1, 0.0f);
                    P.set(2, 2, 0.0f);

                    for (p = -1; p <= 1; p += 2) { /* p,q,r are -1 or +1      */

                        for (q = -1; q <= 1; q += 2) { /* and go into rows #1,2,3 */

                            for (r = -1; r <= 1; r += 2) {
                                P.set(0, i - 1, p);
                                P.set(1, j - 1, q);
                                P.set(2, k - 1, r);
                                detP = P.det(); /* sign of permutation */

                                if ((detP * detQ) <= 0.0f) {
                                    continue; /* doesn't match sign of Q */
                                }

                                M = P.times(Q);

                                /* angle of M rotation = 2.0f*acos(0.5f*sqrt(1.0f+trace(M)))       */
                                /* we want largest trace(M) == smallest angle == M nearest to I */

                                val = M.get(0, 0) + M.get(1, 1) + M.get(2, 2); /* trace */

                                if (val > vbest) {
                                    vbest = val;
                                    ibest = i;
                                    jbest = j;
                                    kbest = k;
                                    pbest = p;
                                    qbest = q;
                                    rbest = r;
                                }
                            }
                        }
                    }
                }
            }
        }

        /* At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.
         *
         * The matrix P that corresponds is the best permutation approximation to Q-inverse; that is, P (approximately)
         * takes (x,y,z) coordinates to the (i,j,k) axes.
         *
         * For example, the first row of P (which contains pbest in column ibest) determines the way the i axis points
         * relative to the anatomical (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis, which is
         * direction P2A (if pbest > 0) or A2P (if pbest < 0).
         *
         * So, using ibest and pbest, we can assign the output code forthe i axis.  Mutatis mutandis for the j and k axes,
         * of course. */

        switch (ibest * pbest) {

            case 1:
                i = FileInfoBase.ORI_R2L_TYPE;
                break;

            case -1:
                i = FileInfoBase.ORI_L2R_TYPE;
                break;

            case 2:
                i = FileInfoBase.ORI_A2P_TYPE;
                break;

            case -2:
                i = FileInfoBase.ORI_P2A_TYPE;
                break;

            case 3:
                i = FileInfoBase.ORI_I2S_TYPE;
                break;

            case -3:
                i = FileInfoBase.ORI_S2I_TYPE;
                break;
        }

        switch (jbest * qbest) {

            case 1:
                j = FileInfoBase.ORI_R2L_TYPE;
                break;

            case -1:
                j = FileInfoBase.ORI_L2R_TYPE;
                break;

            case 2:
                j = FileInfoBase.ORI_A2P_TYPE;
                break;

            case -2:
                j = FileInfoBase.ORI_P2A_TYPE;
                break;

            case 3:
                j = FileInfoBase.ORI_I2S_TYPE;
                break;

            case -3:
                j = FileInfoBase.ORI_S2I_TYPE;
                break;

            default:
                j = 1;
        }

        switch (kbest * rbest) {

            case 1:
                k = FileInfoBase.ORI_R2L_TYPE;
                break;

            case -1:
                k = FileInfoBase.ORI_L2R_TYPE;
                break;

            case 2:
                k = FileInfoBase.ORI_A2P_TYPE;
                break;

            case -2:
                k = FileInfoBase.ORI_P2A_TYPE;
                break;

            case 3:
                k = FileInfoBase.ORI_I2S_TYPE;
                break;

            case -3:
                k = FileInfoBase.ORI_S2I_TYPE;
                break;
        }

        axisOrientation[0] = i;
        axisOrientation[1] = j;
        axisOrientation[2] = k;

        return axisOrientation;
    }

    /**
     * Polar decomposition of a 3x3 matrix: finds the closest orthogonal matrix to input A (in both Frobenius and L2
     * norms). Algorithm is that from NJ Higham, SIAM JSci Stat Comput, 7:1160-1174.
     *
     * @param   A  DOCUMENT ME!
     *
     * @return  DOCUMENT ME!
     */
    private Matrix mat33_polar(Matrix A) {
        Matrix X, Y, Z;
        double alp, bet, gam, gmi;
        double dif = 1.0f;
        int k = 0;
        double val;

        X = A.copy();
        Z = new Matrix(3, 3);

        // force matrix to be nonsingular
        gam = X.det();

        while (gam == 0.0f) { // perturb matrix
            gam = 0.00001 * (0.001f + mat33_rownorm(X));
            val = X.get(0, 0);
            X.set(0, 0, val + gam);
            val = X.get(1, 1);
            X.set(1, 1, val + gam);
            val = X.get(2, 2);
            X.set(2, 2, val + gam);
            gam = X.det();
        }

        while (true) {
            Y = X.inverse();

            if (dif > 0.3f) { // far from convergence
                alp = Math.sqrt(mat33_rownorm(X) * mat33_colnorm(X));
                bet = Math.sqrt(mat33_rownorm(Y) * mat33_colnorm(Y));
                gam = Math.sqrt(bet / alp);
                gmi = 1.0f / gam;
            } else {
                gam = gmi = 1.0f; // close to convergence
            }

            Z.set(0, 0, 0.5f * ((gam * X.get(0, 0)) + (gmi * Y.get(0, 0))));
            Z.set(0, 1, 0.5f * ((gam * X.get(0, 1)) + (gmi * Y.get(1, 0))));
            Z.set(0, 2, 0.5f * ((gam * X.get(0, 2)) + (gmi * Y.get(2, 0))));
            Z.set(1, 0, 0.5f * ((gam * X.get(1, 0)) + (gmi * Y.get(0, 1))));
            Z.set(1, 1, 0.5f * ((gam * X.get(1, 1)) + (gmi * Y.get(1, 1))));
            Z.set(1, 2, 0.5f * ((gam * X.get(1, 2)) + (gmi * Y.get(2, 1))));
            Z.set(2, 0, 0.5f * ((gam * X.get(2, 0)) + (gmi * Y.get(0, 2))));
            Z.set(2, 1, 0.5f * ((gam * X.get(2, 1)) + (gmi * Y.get(1, 2))));
            Z.set(2, 2, 0.5f * ((gam * X.get(2, 2)) + (gmi * Y.get(2, 2))));

            dif = Math.abs(Z.get(0, 0) - X.get(0, 0)) + Math.abs(Z.get(0, 1) - X.get(0, 1)) +
                  Math.abs(Z.get(0, 2) - X.get(0, 2)) + Math.abs(Z.get(1, 0) - X.get(1, 0)) +
                  Math.abs(Z.get(1, 1) - X.get(1, 1)) + Math.abs(Z.get(1, 2) - X.get(1, 2)) +
                  Math.abs(Z.get(2, 0) - X.get(2, 0)) + Math.abs(Z.get(2, 1) - X.get(2, 1)) +
                  Math.abs(Z.get(2, 2) - X.get(2, 2));

            k = k + 1;

            if ((k > 100) || (dif < 3.0e-6f)) {
                break; // convergence or exhaustion
            }

            X = Z.copy();
        }

        return Z;
    }

    /**
     * max column norm of 3x3 matrix.
     *
     * @param   A  DOCUMENT ME!
     *
     * @return  DOCUMENT ME!
     */
    private double mat33_colnorm(Matrix A) {
        double r1, r2, r3;
        r1 = Math.abs(A.get(0, 0)) + Math.abs(A.get(1, 0)) + Math.abs(A.get(2, 0));
        r2 = Math.abs(A.get(0, 1)) + Math.abs(A.get(1, 1)) + Math.abs(A.get(2, 1));
        r3 = Math.abs(A.get(0, 2)) + Math.abs(A.get(1, 2)) + Math.abs(A.get(2, 2));

        if (r1 < r2) {
            r1 = r2;
        }

        if (r1 < r3) {
            r1 = r3;
        }

        return r1;
    }

    /**
     * max row norm of 3x3 matrix.
     *
     * @param   A  DOCUMENT ME!
     *
     * @return  DOCUMENT ME!
     */
    private double mat33_rownorm(Matrix A) {
        double r1, r2, r3;
        r1 = Math.abs(A.get(0, 0)) + Math.abs(A.get(0, 1)) + Math.abs(A.get(0, 2));
        r2 = Math.abs(A.get(1, 0)) + Math.abs(A.get(1, 1)) + Math.abs(A.get(1, 2));
        r3 = Math.abs(A.get(2, 0)) + Math.abs(A.get(2, 1)) + Math.abs(A.get(2, 2));

        if (r1 < r2) {
            r1 = r2;
        }

        if (r1 < r3) {
            r1 = r3;
        }

        return r1;
    }
}

