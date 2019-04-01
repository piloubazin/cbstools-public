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
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;
import de.mpg.cbs.core.intensity.IntensityFlashT2sFitting;

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
	
	IntensityFlashT2sFitting algorithm;
	
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
		
		algorithm = new IntensityFlashT2sFitting();
		inputParams.setPackage(algorithm.getPackage());
		inputParams.setCategory(algorithm.getCategory());
		inputParams.setLabel(algorithm.getLabel());
		inputParams.setName(algorithm.getName());

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(References.getAuthor(algorithm.getAlgorithmAuthors()[0]));
		info.setAffiliation(algorithm.getAffiliation());
		info.setDescription(algorithm.getDescription());
		
		info.setVersion(algorithm.getVersion());
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(s0Param=new ParamVolume("S0 Volume",null,-1,-1,-1,-1));
		outputParams.add(t2Param=new ParamVolume("T2 Volume",null,-1,-1,-1,-1));
		outputParams.add(r2Param=new ParamVolume("R2 Volume",null,-1,-1,-1,-1));
		outputParams.add(errParam=new ParamVolume("Residuals Volume",null,-1,-1,-1,-1));
		
		outputParams.setName("t2s_fitted_images");
		outputParams.setLabel("T2s fitted Images");
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageDataFloat vol1Img = new ImageDataFloat(vol1Param.getImageData());
		int nx = vol1Img.getRows();
		int ny = vol1Img.getCols();
		int nz = vol1Img.getSlices();
		int nxyz = nx*ny*nz;
		ImageHeader header = vol1Img.getHeader();
		String name = vol1Param.getImageData().getName();
        float[] res = Interface.getResolutions(vol1Param);
		
		float[][][] vol1;
		if (Interface.isImage4D(vol1Param)) {
		    float[][][][] tmp = vol1Img.toArray4d();
		    vol1 = new float[nx][ny][nz];
		    for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		        vol1[x][y][z] = tmp[x][y][z][0];
		    }
		} else vol1 = vol1Img.toArray3d();
		
		ImageDataFloat vol2Img = new ImageDataFloat(vol2Param.getImageData());
		float[][][] vol2;
		if (Interface.isImage4D(vol2Param)) {
		    float[][][][] tmp = vol2Img.toArray4d();
		    vol2 = new float[nx][ny][nz];
		    for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		        vol2[x][y][z] = tmp[x][y][z][0];
		    }
		} else vol2 = vol2Img.toArray3d();
		
		int nimg = 2;
		float[][][] vol3=null, vol4=null, vol5=null, vol6=null, vol7=null, vol8=null;
		
		if (vol3Param.getImageData()!=null) {
			ImageDataFloat vol3Img = new ImageDataFloat(vol3Param.getImageData());
			if (Interface.isImage4D(vol3Param)) {
                float[][][][] tmp = vol3Img.toArray4d();
                vol3 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol3[x][y][z] = tmp[x][y][z][0];
                }
            } else vol3 = vol3Img.toArray3d();
		    nimg++;
		}
		if (vol4Param.getImageData()!=null) {
			ImageDataFloat vol4Img = new ImageDataFloat(vol4Param.getImageData());
			if (Interface.isImage4D(vol4Param)) {
                float[][][][] tmp = vol4Img.toArray4d();
                vol4 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol4[x][y][z] = tmp[x][y][z][0];
                }
            } else vol4 = vol4Img.toArray3d();
			nimg++;
		}
		if (vol5Param.getImageData()!=null) {
			ImageDataFloat vol5Img = new ImageDataFloat(vol5Param.getImageData());
			if (Interface.isImage4D(vol5Param)) {
                float[][][][] tmp = vol5Img.toArray4d();
                vol5 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol5[x][y][z] = tmp[x][y][z][0];
                }
            } else vol5 = vol5Img.toArray3d();
			nimg++;
		}
		if (vol6Param.getImageData()!=null) {
			ImageDataFloat vol6Img = new ImageDataFloat(vol6Param.getImageData());
			if (Interface.isImage4D(vol6Param)) {
                float[][][][] tmp = vol6Img.toArray4d();
                vol6 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol6[x][y][z] = tmp[x][y][z][0];
                }
            } else vol6 = vol6Img.toArray3d();
			nimg++;
		}
		if (vol7Param.getImageData()!=null) {
			ImageDataFloat vol7Img = new ImageDataFloat(vol7Param.getImageData());
			if (Interface.isImage4D(vol7Param)) {
                float[][][][] tmp = vol7Img.toArray4d();
                vol7 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol7[x][y][z] = tmp[x][y][z][0];
                }
            } else vol7 = vol7Img.toArray3d();
			nimg++;
		}
		if (vol8Param.getImageData()!=null) {
			ImageDataFloat vol8Img = new ImageDataFloat(vol8Param.getImageData());
			if (Interface.isImage4D(vol8Param)) {
                float[][][][] tmp = vol8Img.toArray4d();
                vol8 = new float[nx][ny][nz];
                for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                    vol8[x][y][z] = tmp[x][y][z][0];
                }
            } else vol8 = vol8Img.toArray3d();
			nimg++;
		}
		float[][] vol = new float[nimg][nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
            vol[0][xyz] = vol1[x][y][z];
            vol[1][xyz] = vol2[x][y][z];
            if (nimg>2) vol[2][xyz] = vol3[x][y][z];
            if (nimg>3) vol[3][xyz] = vol4[x][y][z];
            if (nimg>4) vol[4][xyz] = vol5[x][y][z];
            if (nimg>5) vol[5][xyz] = vol6[x][y][z];
            if (nimg>6) vol[6][xyz] = vol7[x][y][z];
            if (nimg>7) vol[7][xyz] = vol8[x][y][z];
        }		
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

        algorithm = new IntensityFlashT2sFitting();
        
        algorithm.setNumberOfEchoes(nimg);
        algorithm.setDimensions(nx, ny, nz);
        algorithm.setResolutions(res);
        for (int  i=0;i<nimg;i++) {
            algorithm.setEchoImageAt(i, vol[i]);
            algorithm.setEchoTimeAt(i, TE[i]);
        }
        
		System.out.println("data loaded...");
		
		algorithm.execute();

		System.out.println("T2* fitting..");
		
		float[][][] tmp = new float[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    tmp[x][y][z] = algorithm.getS0Image()[xyz];
		}
		ImageDataFloat s0Data = new ImageDataFloat(tmp);
		//s0Data.setHeader(header);
		s0Data.setName(name+"_s0");
		s0Param.setValue(s0Data);
		s0Data = null;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    tmp[x][y][z] = algorithm.getR2sImage()[xyz];
		}
		ImageDataFloat r2Data = new ImageDataFloat(tmp);
		//r2Data.setHeader(header);
		r2Data.setName(name+"_r2");
		r2Param.setValue(r2Data);
		r2Data = null;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    tmp[x][y][z] = algorithm.getT2sImage()[xyz];
		}
		ImageDataFloat t2Data = new ImageDataFloat(tmp);
		//t2Data.setHeader(header);
		t2Data.setName(name+"_t2");
		t2Param.setValue(t2Data);
		t2Data = null;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    int xyz = x+nx*y+nx*ny*z;
		    tmp[x][y][z] = algorithm.getResidualImage()[xyz];
		}
		ImageDataFloat errData = new ImageDataFloat(tmp);
		//errData.setHeader(header);
		errData.setName(name+"_err");
		errParam.setValue(errData);
		errData = null;
		
		System.out.println("Done.");
	}
}
