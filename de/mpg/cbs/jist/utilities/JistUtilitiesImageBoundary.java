package de.mpg.cbs.jist.utilities;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

import org.apache.commons.math3.util.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistUtilitiesImageBoundary extends ProcessingAlgorithm{
	ParamVolume volParam;
	ParamVolume resultVolParam;
	ParamOption operationParam;
	ParamInteger distParam;
	ParamBoolean addboundaryParam;
	
	private static final String cvsversion = "$Revision: 1.10 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Sets the values at the image boundary to specific values";
	private static final String longDescription = "(we can use zero, min, max)";

	private static final byte ZERO = 0;
	private static final byte MIN = 1;
	private static final byte MAX = 2;
	private static final byte NOISE = 3;

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolume("Image Volume"));
		inputParams.add(operationParam=new ParamOption("Image boundary value",new String[]{"zero","min","max","noise"}));
		inputParams.add(distParam=new ParamInteger("Image boundary size (voxels)", 1, 100, 1));
		inputParams.add(addboundaryParam=new ParamBoolean("resize image to add boundary", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Set Image Boundary");
		inputParams.setName("SetImageBoundary");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion(revnum);
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultVolParam=new ParamVolume("Result Volume",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageDataFloat vol = new ImageDataFloat(volParam.getImageData());
		int nx=vol.getRows();
		int ny=vol.getCols();
		int nz=vol.getSlices();
		int nt=vol.getComponents();
		
		ImageDataFloat resultData;
		
		if (nt==1) {
			float[][][] image = vol.toArray3d();
			
			// main algorithm
			float Imin = ImageStatistics.minimum(image,nx,ny,nz);
			float Imax = ImageStatistics.maximum(image,nx,ny,nz);
			int d = distParam.getValue().intValue();
			
			int nbx,nby,nbz,b;
			if (addboundaryParam.getValue().booleanValue()) {
				b = d;
				nbx = nx+2*b;
				nby = ny+2*b;
				nbz = nz+2*b;
			} else {
				b = 0;
				nbx = nx;
				nby = ny;
				nbz = nz;
			}
			float[][][] result = new float[nbx][nby][nbz];
			
					
			byte op = ZERO;
			if (operationParam.getValue().equals("min")) op = MIN;
			else if (operationParam.getValue().equals("max")) op = MAX;
			else if (operationParam.getValue().equals("noise")) op = NOISE;
			
			for (int x=0;x<nbx;x++) for (int y=0;y<nby;y++) for (int z=0;z<nbz;z++) {
				if (x<d || x>=nbx-d || y<d || y>=nby-d || z<d || z>=nbz-d) {
					if (op==ZERO) result[x][y][z] = 0.0f;
					else if (op==MIN) result[x][y][z] = Imin;
					else if (op==MAX) result[x][y][z] = Imax;
					else if (op==NOISE) result[x][y][z] = (float)(0.01f*FastMath.random()*(Imax-Imin)+Imin);
				} else {
					result[x][y][z] = image[x-b][y-b][z-b];
				}
			}
			resultData = new ImageDataFloat(result);	
		} else {
			float[][][][] image = vol.toArray4d();
			
			// main algorithm	
			float Imin = ImageStatistics.minimum(image,nx,ny,nz,nt);
			float Imax = ImageStatistics.maximum(image,nx,ny,nz,nt);
			int d = distParam.getValue().intValue();
			
			int nbx,nby,nbz,b;
			if (addboundaryParam.getValue().booleanValue()) {
				b = d;
				nbx = nx+2*b;
				nby = ny+2*b;
				nbz = nz+2*b;
			} else {
				b = 0;
				nbx = nx;
				nby = ny;
				nbz = nz;
			}
			float[][][][] result = new float[nbx][nby][nbz][nt];
			
			byte op = ZERO;
			if (operationParam.getValue().equals("min")) op = MIN;
			else if (operationParam.getValue().equals("max")) op = MAX;
			
			for (int x=0;x<nbx;x++) for (int y=0;y<nby;y++) for (int z=0;z<nbz;z++) for (int t=0;t<nt;t++) {
				if (x<d || x>=nbx-d || y<d || y>=nby-d || z<d || z>=nbz-d) {
					if (op==ZERO) result[x][y][z][t] = 0.0f;
					else if (op==MIN) result[x][y][z][t] = Imin;
					else if (op==MAX) result[x][y][z][t] = Imax;
				} else {
					result[x][y][z][t] = image[x-b][y-b][z-b][t];
				}
			}
			resultData = new ImageDataFloat(result);	
		}
		resultData.setHeader(vol.getHeader());
		resultData.setName(vol.getName()+"_imb");
		resultVolParam.setValue(resultData);
		resultData = null;
	}
}
