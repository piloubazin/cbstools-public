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
public class JistUtilitiesMirrorImage extends ProcessingAlgorithm{
	ParamVolume volParam;
	ParamVolume resultVolParam;
	ParamOption directionParam;
	ParamInteger distParam;
	ParamBoolean addboundaryParam;
	
	private static final String cvsversion = "$Revision: 1.10 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Sets the values at the image boundary to specific values";
	private static final String longDescription = "(we can use zero, min, max)";

	private static final String[] directions = {"-X","+X","-Y","+Y","-Z","+Z"};
	
	private static final byte ZERO = 0;
	private static final byte MIN = 1;
	private static final byte MAX = 2;
	private static final byte NOISE = 3;

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolume("Image Volume"));
		inputParams.add(directionParam=new ParamOption("Mirroring direction",directions));
		inputParams.add(distParam=new ParamInteger("Added boundary size (voxels)", -100, 100, 0));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Build Mirror Image");
		inputParams.setName("BuildMirrorImage");

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
			int d = distParam.getValue().intValue();
			
			int nbx,nby,nbz;
			nbx = nx;
			nby = ny;
			nbz = nz;
			if (directionParam.getValue().endsWith("X")) {
				nbx = 2*nx+2*d;
			} else if (directionParam.getValue().endsWith("Y")) {
				nby = 2*ny+2*d;
			} else if (directionParam.getValue().endsWith("Z")) {
				nbz = 2*nz+2*d;
			}
			float[][][] result = new float[nbx][nby][nbz];
								
			if (directionParam.getValue().equals("-X")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					result[x][y][z] = image[nx-1-x][y][z];
				}
				for (int x=nx+2*d;x<nbx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					result[x][y][z] = Numerics.max(result[x][y][z], image[x-nx-2*d][y][z]);
				}
			}
			resultData = new ImageDataFloat(result);	
		} else {
			float[][][][] image = vol.toArray4d();
			
			
			// main algorithm
			int d = distParam.getValue().intValue();
			
			int nbx,nby,nbz;
			nbx = nx;
			nby = ny;
			nbz = nz;
			if (directionParam.getValue().endsWith("X")) {
				nbx = 2*nx+2*d;
			} else if (directionParam.getValue().endsWith("Y")) {
				nby = 2*ny+2*d;
			} else if (directionParam.getValue().endsWith("Z")) {
				nbz = 2*nz+2*d;
			}
			float[][][][] result = new float[nbx][nby][nbz][nt];
								
			if (directionParam.getValue().equals("-X")) {
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
					result[x][y][z][t] = image[nx-1-x][y][z][t];
				}
				for (int x=nx+2*d;x<nbx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
					result[x][y][z][t] = Numerics.max(result[x][y][z][t], image[x-nx-2*d][y][z][t]);
				}
			}
			resultData = new ImageDataFloat(result);	
		}
		resultData.setHeader(vol.getHeader());
		resultData.setName(vol.getName()+"_mirr");
		resultVolParam.setValue(resultData);
		resultData = null;
	}
}
