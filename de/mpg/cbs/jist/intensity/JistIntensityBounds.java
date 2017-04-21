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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensityBounds extends ProcessingAlgorithm{
	ParamVolume volParam;
	ParamVolume resultVolParam;
	ParamOption operationParam;
	ParamBoolean normalizeParam;
	ParamFloat minParam, maxParam;
	
	private static final String cvsversion = "$Revision: 1.1f0 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Perform simple intensity thresholding into a [min, max] interval.";
	private static final String longDescription = "(options allow to use raw intensities, normalized intensities, and robust normalization)";


	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolume("Image Volume"));
		inputParams.add(operationParam=new ParamOption("Bounds definition",new String[]{"raw_intensity","normalized","robust"}));
		inputParams.add(minParam=new ParamFloat("Min", -1e10f, 1e10f, 0.0f));
		inputParams.add(maxParam=new ParamFloat("Max", -1e10f, 1e10f, 1.0f));
		inputParams.add(normalizeParam=new ParamBoolean("Normalize output", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.prev");
		inputParams.setLabel("Intensity Bounds");
		inputParams.setName("Intensity_Bounds");

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
		float[][][] image = vol.toArray3d();
		
		// main algorithm
		float[][][] result = new float[nx][ny][nz];
		ImageDataFloat resultData;
		
		float Imin, Imax;
		float minVal = minParam.getValue().floatValue();
		float maxVal = maxParam.getValue().floatValue();
		
		if (operationParam.getValue().equals("robust")) {
			if (minVal<0 || minVal>1 || maxVal<0 || maxVal>1) {
				System.err.println("IntensityBound: min, max parameters must be in [0,1]!!");
				return;
			}
			Imin = ImageStatistics.robustMinimum(image, minVal, 4, nx, ny, nz);
			Imax = ImageStatistics.robustMinimum(image, maxVal, 4, nx, ny, nz);
		} else if (operationParam.getValue().equals("normalized")) {
			float I1 = ImageStatistics.minimum(image,nx,ny,nz);
			float I2 = ImageStatistics.maximum(image,nx,ny,nz);
			Imin = I1 + minVal*(I2-I1);
			Imax = I1 + maxVal*(I2-I1);
		} else {
			Imin = minVal;
			Imax = maxVal;
		}
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			result[x][y][z] = Numerics.bounded(image[x][y][z], Imin, Imax);
			
			if (normalizeParam.getValue().booleanValue()) {
				result[x][y][z] = (result[x][y][z]-Imin)/(Imax-Imin);
			}
		}
		resultData = new ImageDataFloat(result);		
		resultData.setHeader(vol.getHeader());
		resultData.setName(vol.getName()+"_bound");
		resultVolParam.setValue(resultData);
		resultData = null;
		result = null;
	}
}
