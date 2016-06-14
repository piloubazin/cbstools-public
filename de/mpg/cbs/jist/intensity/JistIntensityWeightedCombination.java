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
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensityWeightedCombination extends ProcessingAlgorithm{
	ParamVolume vol1Param;
	ParamVolume vol2Param;
	ParamVolume weight1Param;
	ParamVolume weight2Param;
	ParamVolume resultVolParam;
	ParamOption operation;

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(vol1Param=new ParamVolume("Volume 1"));
		inputParams.add(weight1Param=new ParamVolume("Weight 1"));
		inputParams.add(vol2Param=new ParamVolume("Volume 2"));
		inputParams.add(weight2Param=new ParamVolume("Weight 2"));
		inputParams.add(operation=new ParamOption("Operation",new String[]{"weighted_average","maximum_or_avg","preferred_max1"}));


		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("Weighted Image Combination");
		inputParams.setName("WeightedImageCombination");


		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("combines two images based on weights (different options are possible: simple average, use the maximum value weight unless equal, use preferentially the first image).");
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultVolParam=new ParamVolume("Result Volume",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageData vol1=vol1Param.getImageData();
		ImageData vol2=vol2Param.getImageData();
		ImageData weight1=weight1Param.getImageData();
		ImageData weight2=weight2Param.getImageData();
		int rows=vol1.getRows();
		int cols=vol1.getCols();
		int slices=vol1.getSlices();
		int comps = vol1.getComponents();
		
		
		String resultName = vol1.getName();
		if (operation.getValue().equals("weighted_average")) {
			resultName+= "_wavg";
		} else
		if (operation.getValue().equals("maximum_or_avg")) {
			resultName+= "_max";
		} else
		if (operation.getValue().equals("preferred_max1")) {
			resultName+= "_max1";
		}
		ImageDataFloat resultVol=new ImageDataFloat(resultName,rows,cols,slices,comps);
		
		double tmp1,tmp2;
		double w1,w2,den;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					for (int l = 0; l < comps; l++) {
						tmp1=vol1.getDouble(i, j, k, l);
						tmp2=vol2.getDouble(i, j, k, l);
						w1=weight1.getDouble(i, j, k, l);
						w2=weight2.getDouble(i, j, k, l);
						den = w1+w2;

						if (operation.getValue().equals("weighted_average")) {
							if (den>0) resultVol.set(i, j, k, l, (w1*tmp1+w2*tmp2)/den);
							else resultVol.set(i, j, k, l, 0.0);
						} else if (operation.getValue().equals("maximum_or_avg")) {
							if (w1>w2) resultVol.set(i, j, k, l, tmp1);
							else if (w1<w2) resultVol.set(i, j, k, l, tmp2);
							else resultVol.set(i, j, k, l, 0.5*tmp1+0.5*tmp2);
						} else if (operation.getValue().equals("preferred_max1")) {
							if (w1>=w2) resultVol.set(i, j, k, l, tmp1);
							else resultVol.set(i, j, k, l, (w1/w2*tmp1+(w2-w1)/w2*tmp2));
						}
					}
				}
			}
		}
		resultVol.setHeader(vol1.getHeader());
		resultVolParam.setValue(resultVol);
	}
}
