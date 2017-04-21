package de.mpg.cbs.jist.brain;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;

import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistBrainFilterStacking extends ProcessingAlgorithm{
	ParamVolume 		vol1Param;
	ParamVolume 		vol2Param;
	ParamVolume 		vol3Param;
	ParamFloat			bgParam;
	
	ParamVolume stackVolParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(vol1Param=new ParamVolume("Dura Filter Image (opt)"));
		inputParams.add(vol2Param=new ParamVolume("Partial Volume Image (opt)"));
		inputParams.add(vol3Param=new ParamVolume("Arteries Filter Image (opt)"));
		vol1Param.setMandatory(false);
		vol2Param.setMandatory(false);
		vol3Param.setMandatory(false);
		inputParams.add(bgParam=new ParamFloat("Background threshold [0-1]", 0.0f, 1.0f, 0.001f));
		bgParam.setDescription("maps low values from any filter to zero, for visualization purposes only.");
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Brain Processing");
		inputParams.setLabel("MP2RAGE Filter Stacking");
		inputParams.setName("FilterStacking");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("combines filter results into a single 'filters' image to input to the brain segmentation module.");
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(stackVolParam=new ParamVolume("Stacked Filter Probabilities"));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		int rows=0,cols=0,slices=0;
		String name=null;
		ImageHeader header=null;
		if (vol1Param.getImageData()!=null) {
			rows=vol1Param.getImageData().getRows();
			cols=vol1Param.getImageData().getCols();
			slices=vol1Param.getImageData().getSlices();
			name = vol1Param.getImageData().getName();
			header = vol1Param.getImageData().getHeader();
		} else if (vol2Param.getImageData()!=null) {
			rows=vol2Param.getImageData().getRows();
			cols=vol2Param.getImageData().getCols();
			slices=vol2Param.getImageData().getSlices();
			name = vol2Param.getImageData().getName();
			header = vol2Param.getImageData().getHeader();
		} else if (vol3Param.getImageData()!=null) {
			rows=vol3Param.getImageData().getRows();
			cols=vol3Param.getImageData().getCols();
			slices=vol3Param.getImageData().getSlices();
			name = vol3Param.getImageData().getName();
			header = vol3Param.getImageData().getHeader();
		}	
		ImageDataFloat stackVol=new ImageDataFloat(rows,cols,slices);
		
		float bg = bgParam.getValue().floatValue();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					
					float val=0.0f;
					if (vol1Param.getImageData()!=null) 
						if (vol1Param.getImageData().getFloat(i, j, k)>bg 
							&& vol1Param.getImageData().getFloat(i, j, k)>val) 
							val = vol1Param.getImageData().getFloat(i, j, k);
					if (vol2Param.getImageData()!=null) 
						if (vol2Param.getImageData().getFloat(i, j, k)>bg 
							&& vol2Param.getImageData().getFloat(i, j, k)>val) 
							val = vol2Param.getImageData().getFloat(i, j, k)+2.0f;
					if (vol3Param.getImageData()!=null) 
						if (vol3Param.getImageData().getFloat(i, j, k)>bg 
							&& ( (val>2.0f && vol3Param.getImageData().getFloat(i, j, k)>val-2.0f)
								|| (vol3Param.getImageData().getFloat(i, j, k)>val) ) ) 
							val = vol3Param.getImageData().getFloat(i, j, k)+4.0f;
					
					stackVol.set(i, j, k, val);
				}
			}
		}
		stackVol.setName(name+"_filters");
		stackVol.setHeader(header);
		stackVolParam.setValue(stackVol);
	}
}
