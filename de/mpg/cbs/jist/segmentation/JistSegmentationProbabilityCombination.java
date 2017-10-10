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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
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
	ParamVolumeCollection 	volParam;
	ParamInteger			sizeParam;
	ParamBoolean			includebgParam;
	ParamBoolean			rescaleParam;

	ParamVolume maxVolParam;
	ParamVolume lbVolParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolumeCollection("Probability Images"));
		inputParams.add(sizeParam=new ParamInteger("Combined result size",1,10,3));
		inputParams.add(includebgParam=new ParamBoolean("Background included",true));
		inputParams.add(rescaleParam=new ParamBoolean("Rescale probabilities",false));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Probability Combination");
		inputParams.setName("ProbabilityCombination");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("combine multiple probability maps into a (truncated) maximum probability image");
		info.setVersion("3.0.9");
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
		int NP = N;
		boolean addbg = (!includebgParam.getValue().booleanValue());
		if (addbg) NP++;
		
		ImageDataFloat maxVol=new ImageDataFloat(rows,cols,slices,M);
		ImageDataUByte lbVol=new ImageDataUByte(rows,cols,slices,M);
		
		float[] max = new float[N];
		float[] min = new float[N];
		for(byte l=0; l<N; l++) {
			max[l] = 1.0f;
			min[l] = 0.0f;
		}
		if (rescaleParam.getValue().booleanValue()) {
			for(byte l=0; l<N; l++) {
				float vmin = 1e15f;
				float vmax = -1e15f;
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						for (int k = 0; k < slices; k++) {
							float val=volParam.getImageDataList().get(l).getFloat(i, j, k, l);
							if (val<vmin) vmin = val;
							if (val>vmax) vmax = val;
						}
					}
				}
				min[l] = vmin;
				max[l] = vmax;
			}			
		}
		
		float[] val = new float[NP];
		float[] best = new float[M];
		byte[] arg = new byte[M];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					
					for(byte l=0; l<N; l++) {
						if (min[l]<max[l]) val[l]=(volParam.getImageDataList().get(l).getFloat(i, j, k, l)-min[l])/(max[l]-min[l]);
						else val[l] = 0.0f;
					}
						
					if (addbg) {
					    /* use the sum or the max?
						val[N] = 1.0f;
						for(byte l=0; l<N; l++) val[N] -= val[l];
						val[N] = Numerics.max(val[N],0.0f);
						*/
						val[N] = 1.0f;
						for(byte l=0; l<N; l++) val[N] = Numerics.min(1.0f-val[l],val[N]);
						val[N] = Numerics.max(val[N],0.0f);
					}
					Numerics.argmax(arg, best, val, M);
					for(byte l=0; l<M; l++){
						maxVol.set(i, j, k, l, best[l]);
						if (arg[l]<N) lbVol.set(i, j, k, l, arg[l]+1);
						else lbVol.set(i, j, k, l, 0);
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
