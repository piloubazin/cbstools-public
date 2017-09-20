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
 * @author Blake Lucas (bclucas@jhu.edu)
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensityCombinations extends ProcessingAlgorithm{
	ParamVolume vol1Param;
	ParamVolume vol2Param;
	ParamVolume resultVolParam;
	ParamOption operation;

	private static final String cvsversion = "$Revision: 1.10 $";
	private static final String revnum = cvsversion.replace("Revision: ", "").replace("$", "").replace(" ", "");
	private static final String shortDescription = "Perform simple image combination operations on two images. The operations include 'Min', 'Max', "
													+"\n 'SignedMax' (=I1 if I1>I2, =-I2 else), 'Stack' (=I1 if I1>I2, =max(I1)+I2+1 else)";
	

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(vol1Param=new ParamVolume("Volume 1"));
		inputParams.add(vol2Param=new ParamVolume("Volume 2"));
		inputParams.add(operation=new ParamOption("Operation",new String[]{"Min","Max","SignedMax","Stack"}));


		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.prev");
		inputParams.setLabel("Image Combinations");
		inputParams.setName("Image_Combinations");


		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription);
		info.setVersion(revnum);
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultVolParam=new ParamVolume("Result Volume",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageData vol1=vol1Param.getImageData();
		ImageData vol2=vol2Param.getImageData();
		int rows=vol1.getRows();
		int cols=vol1.getCols();
		int slices=vol1.getSlices();
		int comps = vol1.getComponents();
		
		String resultName = vol1.getName();
		if (operation.getValue().equals("Min")) {
			resultName+= "_min";
		} else
		if (operation.getValue().equals("Max")) {
			resultName+= "_max";
		} else
		if (operation.getValue().equals("SignedMax")) {
			resultName+= "_sgm";
		} else
		if (operation.getValue().equals("Stack")) {
			resultName+= "_stc";
		}
		ImageDataFloat resultVol=new ImageDataFloat(resultName,rows,cols,slices,comps);
		
		double max1 = 0.0;
		if (operation.getValue().equals("Stack")) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					for (int k = 0; k < slices; k++) {
						for (int l = 0; l < comps; l++) {
							max1=Numerics.max(max1,vol1.getDouble(i, j, k, l));
						}
					}
				}
			}
		}
		
		double tmp1,tmp2;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					for (int l = 0; l < comps; l++) {
						tmp1=vol1.getDouble(i, j, k, l);
						tmp2=vol2.getDouble(i, j, k, l);

						if (operation.getValue().equals("Min")) {
							resultVol.set(i, j, k, l, Numerics.min(tmp1,tmp2));
						} else
						if (operation.getValue().equals("Max")) {
							resultVol.set(i, j, k, l, Numerics.max(tmp1,tmp2));
						} else
						if (operation.getValue().equals("SignedMax")) {
							if (tmp1>=tmp2) resultVol.set(i, j, k, l, tmp1);
							else resultVol.set(i, j, k, l, -tmp2);
						} else 
						if (operation.getValue().equals("Stack")) {
							if (tmp1>=tmp2) resultVol.set(i, j, k, l, tmp1);
							else resultVol.set(i, j, k, l, max1+tmp2+1.0);
						}
					}
				}
			}
		}
		resultVol.setHeader(vol1.getHeader());
		resultVolParam.setValue(resultVol);
	}
}
