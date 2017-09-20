package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamMatrix;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import Jama.Matrix;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistRegistrationGenerateLinearMapping extends ProcessingAlgorithm{
	ParamVolume volParam;
	ParamVolume srcParam;
	ParamMatrix	 transParam;
	ParamBoolean invertParam;
	ParamVolume coordParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolume("Reference Image Volume",null,-1,-1,-1,-1));
		
		Matrix identityMtx = new Matrix(4,4);
		for(int i = 0; i < 4; i++) identityMtx.set(i,i,1);
		inputParams.add(transParam=new ParamMatrix("Matrix (opt)",identityMtx));
		inputParams.add(invertParam=new ParamBoolean("Invert matrix?",false));
		inputParams.add(srcParam=new ParamVolume("Source Image Volume (opt)",null,-1,-1,-1,-1));
		srcParam.setMandatory(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		inputParams.setLabel("Generate Linear Mapping");
		inputParams.setName("GenerateLinearMapping");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("generate a coordinate mapping for this image and linear transformation (for applying linear and non-linear deformations together)");
		info.setVersion("3.0.3");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(coordParam=new ParamVolume("Coordinate Mapping",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageDataFloat vol = new ImageDataFloat(volParam.getImageData());
		int nx=vol.getRows();
		int ny=vol.getCols();
		int nz=vol.getSlices();
		
		int X=0, Y=1, Z=2, T=3;
		float rsx = 1, rsy = 1, rsz = 1;
		float rtx = 1, rty = 1, rtz = 1;
		rtx = volParam.getImageData().getHeader().getDimResolutions()[X];
		rty = volParam.getImageData().getHeader().getDimResolutions()[Y];
		rtz = volParam.getImageData().getHeader().getDimResolutions()[Z];
		// source image resolution, if different from target
		if (srcParam.getImageData()!=null) {
			rsx = srcParam.getImageData().getHeader().getDimResolutions()[X];
			rsy = srcParam.getImageData().getHeader().getDimResolutions()[Y];
			rsz = srcParam.getImageData().getHeader().getDimResolutions()[Z];
		} else {
			rsx = rtx;
			rsy = rty;
			rsz = rtz;
		}
		
		// main algorithm
		float[][][][] result = new float[nx][ny][nz][3];
		ImageDataFloat resultData;
		
		Matrix transMatrix = transParam.getValue();
		// invert by default
		if (!invertParam.getValue().booleanValue()) transMatrix = transMatrix.inverse();
			
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			double xt = transMatrix.get(X,X)*x*rtx + transMatrix.get(X,Y)*y*rty + transMatrix.get(X,Z)*z*rtz + transMatrix.get(X,T);
			double yt = transMatrix.get(Y,X)*x*rtx + transMatrix.get(Y,Y)*y*rty + transMatrix.get(Y,Z)*z*rtz + transMatrix.get(Y,T);
			double zt = transMatrix.get(Z,X)*x*rtx + transMatrix.get(Z,Y)*y*rty + transMatrix.get(Z,Z)*z*rtz + transMatrix.get(Z,T);
			double wt = transMatrix.get(T,X)*x*rtx + transMatrix.get(T,Y)*y*rty + transMatrix.get(T,Z)*z*rtz + transMatrix.get(T,T);
			
			/* or transposed?
			double xt = transMatrix.get(X,X)*x*rtx + transMatrix.get(Y,X)*y*rty + transMatrix.get(Z,X)*z*rtz + transMatrix.get(T,X);
			double yt = transMatrix.get(X,Y)*x*rtx + transMatrix.get(Y,Y)*y*rty + transMatrix.get(Z,Y)*z*rtz + transMatrix.get(T,Y);
			double zt = transMatrix.get(X,Z)*x*rtx + transMatrix.get(Y,Z)*y*rty + transMatrix.get(Z,Z)*z*rtz + transMatrix.get(T,Z);
			double wt = transMatrix.get(X,T)*x*rtx + transMatrix.get(Y,T)*y*rty + transMatrix.get(Z,T)*z*rtz + transMatrix.get(T,T);
			*/
			
			/* no normalization by last term
			result[x][y][z][X] = (float)(xt/wt/rsx);
			result[x][y][z][Y] = (float)(yt/wt/rsy);
			result[x][y][z][Z] = (float)(zt/wt/rsz);
			*/
			result[x][y][z][X] = (float)(xt/rsx);
			result[x][y][z][Y] = (float)(yt/rsy);
			result[x][y][z][Z] = (float)(zt/rsz);
		}
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(vol.getHeader());
		resultData.setName(vol.getName()+"_mapping");
		coordParam.setValue(resultData);
		resultData = null;
		result = null;
	}
}
