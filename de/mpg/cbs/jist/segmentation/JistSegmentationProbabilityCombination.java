package de.mpg.cbs.jist.segmentation;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistSegmentationProbabilityCombination extends ProcessingAlgorithm{
	ParamVolumeCollection volParam;
	ParamInteger			sizeParam;

	ParamVolume maxVolParam;
	ParamVolume lbVolParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolumeCollection("Probability Images"));
		inputParams.add(sizeParam=new ParamInteger("Combined result size",1,10,3));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Probability Combination");
		inputParams.setName("ProbabilityCombination");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("combine multiple probability maps into a (truncated) maximum probability image");
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(maxVolParam=new ParamVolume("Maximum Probability",null,-1,-1,-1,-1));
		outputParams.add(lbVolParam=new ParamVolume("Maximum Label",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		int rows=volParam.getImageDataList().get(0).getRows();
		int cols=volParam.getImageDataList().get(0).getCols();
		int slices=volParam.getImageDataList().get(0).getSlices();
		
		int N = volParam.getImageDataList().size();
		int M = Numerics.min(N, sizeParam.getValue().intValue());
		
		ImageDataFloat maxVol=new ImageDataFloat(rows,cols,slices,M);
		ImageDataUByte lbVol=new ImageDataUByte(rows,cols,slices,M);
		
		float[] val = new float[N];
		float[] best = new float[M];
		byte[] arg = new byte[M];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					
					for(byte l=0; l<N; l++) 
						val[l]=volParam.getImageDataList().get(l).getFloat(i, j, k, l);
					
					Numerics.argmax(arg, best, val, M);
					for(byte l=0; l<M; l++){
						maxVol.set(i, j, k, l, best[l]);
						lbVol.set(i, j, k, l, arg[l]);
					}
				}
			}
		}
		lbVol.setName(volParam.getImageDataList().get(0).getName()+"_lb");
		maxVol.setName(volParam.getImageDataList().get(0).getName()+"_max");
		
		lbVol.setHeader(volParam.getImageDataList().get(0).getHeader());
		maxVol.setHeader(volParam.getImageDataList().get(0).getHeader());
		
		lbVolParam.setValue(lbVol);
		maxVolParam.setValue(maxVol);
	}
}
