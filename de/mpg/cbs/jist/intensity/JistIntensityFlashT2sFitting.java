package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensityFlashT2sFitting extends ProcessingAlgorithm{
	ParamVolume vol1Param;
	ParamVolume vol2Param;
	ParamVolume vol3Param;
	ParamVolume vol4Param;
	ParamVolume vol5Param;
	ParamVolume vol6Param;
	ParamVolume vol7Param;
	ParamVolume vol8Param;
	ParamFloat te1Param;
	ParamFloat te2Param;
	ParamFloat te3Param;
	ParamFloat te4Param;
	ParamFloat te5Param;
	ParamFloat te6Param;
	ParamFloat te7Param;
	ParamFloat te8Param;
	//ParamBoolean noiseParam;
	
	ParamVolume s0Param;
	ParamVolume t2Param;
	ParamVolume r2Param;
	ParamVolume errParam;
	
	private static final String cvsversion = "$Revision: 1.1f0 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Fits multiple echos of a T2* sequence with an exponential to get T2* estimates";
	private static final String longDescription = "";


	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(vol1Param=new ParamVolume("Echo 1"));
		inputParams.add(te1Param=new ParamFloat("TE 1"));
		inputParams.add(vol2Param=new ParamVolume("Echo 2"));
		inputParams.add(te2Param=new ParamFloat("TE 2"));
		inputParams.add(vol3Param=new ParamVolume("Echo 3"));
		inputParams.add(te3Param=new ParamFloat("TE 3"));
		vol3Param.setMandatory(false);
		te3Param.setMandatory(false);
		inputParams.add(vol4Param=new ParamVolume("Echo 4"));
		inputParams.add(te4Param=new ParamFloat("TE 4"));
		vol4Param.setMandatory(false);
		te4Param.setMandatory(false);
		inputParams.add(vol5Param=new ParamVolume("Echo 5"));
		inputParams.add(te5Param=new ParamFloat("TE 5"));
		vol5Param.setMandatory(false);
		te5Param.setMandatory(false);
		inputParams.add(vol6Param=new ParamVolume("Echo 6"));
		inputParams.add(te6Param=new ParamFloat("TE 6"));
		vol6Param.setMandatory(false);
		te6Param.setMandatory(false);
		inputParams.add(vol7Param=new ParamVolume("Echo 7"));
		inputParams.add(te7Param=new ParamFloat("TE 7"));
		vol7Param.setMandatory(false);
		te7Param.setMandatory(false);
		inputParams.add(vol8Param=new ParamVolume("Echo 8"));
		inputParams.add(te8Param=new ParamFloat("TE 8"));
		vol8Param.setMandatory(false);
		te8Param.setMandatory(false);
		//inputParams.add(noiseParam=new ParamBoolean("Estimate background noise"));
		

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.devel");
		inputParams.setLabel("Flash T2* Fitting");
		inputParams.setName("FlashT2sFitting");


		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion("3.1.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(s0Param=new ParamVolume("S0 Volume"));
		outputParams.add(t2Param=new ParamVolume("T2 Volume"));
		outputParams.add(r2Param=new ParamVolume("R2 Volume"));
		outputParams.add(errParam=new ParamVolume("Residuals Volume"));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageDataFloat vol1Img = new ImageDataFloat(vol1Param.getImageData());
		float[][][] vol1 = vol1Img.toArray3d();
		int nx = vol1Img.getRows();
		int ny = vol1Img.getCols();
		int nz = vol1Img.getSlices();
		ImageDataFloat vol2Img = new ImageDataFloat(vol2Param.getImageData());
		float[][][] vol2 = vol2Img.toArray3d();
		int nimg = 2;
		float[][][] vol3=null, vol4=null, vol5=null, vol6=null, vol7=null, vol8=null;
		if (vol3Param.getImageData()!=null) {
			ImageDataFloat vol3Img = new ImageDataFloat(vol3Param.getImageData());
			vol3 = vol3Img.toArray3d();
			nimg++;
		}
		if (vol4Param.getImageData()!=null) {
			ImageDataFloat vol4Img = new ImageDataFloat(vol4Param.getImageData());
			vol4 = vol4Img.toArray3d();
			nimg++;
		}
		if (vol5Param.getImageData()!=null) {
			ImageDataFloat vol5Img = new ImageDataFloat(vol5Param.getImageData());
			vol5 = vol5Img.toArray3d();
			nimg++;
		}
		if (vol6Param.getImageData()!=null) {
			ImageDataFloat vol6Img = new ImageDataFloat(vol6Param.getImageData());
			vol6 = vol6Img.toArray3d();
			nimg++;
		}
		if (vol7Param.getImageData()!=null) {
			ImageDataFloat vol7Img = new ImageDataFloat(vol7Param.getImageData());
			vol7 = vol7Img.toArray3d();
			nimg++;
		}
		if (vol8Param.getImageData()!=null) {
			ImageDataFloat vol8Img = new ImageDataFloat(vol8Param.getImageData());
			vol8 = vol8Img.toArray3d();
			nimg++;
		}
		float[][][][] vol = new float[nimg][][][];
		vol[0] = vol1;
		vol[1] = vol2;
		if (nimg>2) vol[2] = vol3;
		if (nimg>3) vol[3] = vol4;
		if (nimg>4) vol[4] = vol5;
		if (nimg>5) vol[5] = vol6;
		if (nimg>6) vol[6] = vol7;
		if (nimg>7) vol[7] = vol8;
		
		// retrieve the corresponding TE
		float[] TE = new float[nimg];
		TE[0] = te1Param.getValue().floatValue();
		TE[1] = te2Param.getValue().floatValue();
		if (nimg>2) TE[2] = te3Param.getValue().floatValue();
		if (nimg>3) TE[3] = te4Param.getValue().floatValue();
		if (nimg>4) TE[4] = te5Param.getValue().floatValue();
		if (nimg>5) TE[5] = te6Param.getValue().floatValue();
		if (nimg>6) TE[6] = te7Param.getValue().floatValue();
		if (nimg>7) TE[7] = te8Param.getValue().floatValue();
		
		// fit the exponential
		double Sx, Sx2, Sy, Sxy, delta;
		Sx = 0.0f; Sx2 = 0.0f; delta = 0.0f;
		for (int i=0;i<nimg;i++) {
			Sx += -TE[i];
			Sx2 += Numerics.square(-TE[i]);
		}
		delta = nimg*Sx2 - Sx*Sx;
		
		float[][][] s0 = new float[nx][ny][nz];
		float[][][] r2 = new float[nx][ny][nz];
		float[][][] t2 = new float[nx][ny][nz];
		float[][][] err = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boolean process=true;
			for (int i=0;i<nimg;i++) if (vol[i][x][y][z]<=0) process=false;
			if (process) {
				Sy = 0.0f;
				Sxy = 0.0f;
				for (int i=0;i<nimg;i++) {
					double ydata = FastMath.log(vol[i][x][y][z]);
					Sy += ydata;
					Sxy += -ydata*TE[i];
				}
				s0[x][y][z] = Numerics.bounded( (float)FastMath.exp( (Sx2*Sy-Sxy*Sx)/delta ), 0, 4000);
				r2[x][y][z] = Numerics.bounded( (float)( (nimg*Sxy-Sx*Sy)/delta ), 0.001f, 1.0f);
				t2[x][y][z] = 1.0f/r2[x][y][z];
				
				err[x][y][z] = 0.0f;
				for (int i=0;i<nimg;i++) {
					err[x][y][z] += Numerics.square( FastMath.log(vol[i][x][y][z]) 
													 - FastMath.log(s0[x][y][z]) + TE[i]*r2[x][y][z] );
				}
				err[x][y][z] = (float)FastMath.sqrt(err[x][y][z]);
			} else {
				
			}
		}
		
		ImageDataFloat s0Data = new ImageDataFloat(s0);
		s0Data.setHeader(vol1Param.getImageData().getHeader());
		s0Data.setName(vol1Param.getImageData().getName()+"_s0");
		s0Param.setValue(s0Data);
		s0Data = null;
		s0 = null;
		
		ImageDataFloat r2Data = new ImageDataFloat(r2);
		r2Data.setHeader(vol1Param.getImageData().getHeader());
		r2Data.setName(vol1Param.getImageData().getName()+"_r2");
		r2Param.setValue(r2Data);
		r2Data = null;
		r2 = null;
		
		ImageDataFloat t2Data = new ImageDataFloat(t2);
		t2Data.setHeader(vol1Param.getImageData().getHeader());
		t2Data.setName(vol1Param.getImageData().getName()+"_t2");
		t2Param.setValue(t2Data);
		t2Data = null;
		t2 = null;
		
		ImageDataFloat errData = new ImageDataFloat(err);
		errData.setHeader(vol1Param.getImageData().getHeader());
		errData.setName(vol1Param.getImageData().getName()+"_err");
		errParam.setValue(errData);
		errData = null;
		err = null;
	}
}
