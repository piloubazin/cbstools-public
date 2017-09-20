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
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistSegmentationProbabilityExtraction extends ProcessingAlgorithm{
	ParamVolume maxVolParam;
	ParamVolume lbVolParam;
	ParamInteger indexParam;
	
	ParamVolumeCollection volParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(maxVolParam=new ParamVolume("Maximum Probability",null,-1,-1,-1,-1));
		maxVolParam.setLoadAndSaveOnValidate(false);
		inputParams.add(lbVolParam=new ParamVolume("Maximum Label",null,-1,-1,-1,-1));
		lbVolParam.setLoadAndSaveOnValidate(false);
		inputParams.add(indexParam=new ParamInteger("Probability to extract (-1 for all)",-1,100,-1));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Probability Extraction");
		inputParams.setName("ProbabilityExtraction");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("extract probabilities from a combined probability map");
		info.setVersion("3.0");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(volParam=new ParamVolumeCollection("Extracted Volumes"));
		volParam.setLoadAndSaveOnValidate(false);
	}

	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		
		ImageDataUByte lbData = new ImageDataUByte(lbVolParam.getImageData());
		ImageDataFloat probData = new ImageDataFloat(maxVolParam.getImageData());
		
		int nx=probData.getRows();
		int ny=probData.getCols();
		int nz=probData.getSlices();
		
		int nc = probData.getComponents();
		
		String probName = probData.getName();
		ImageHeader probHeader = probData.getHeader();
		
		byte[][][][] labels; 
		float[][][][] probas;
		byte[][][] lb0;
		if (nc>1) {
			labels = lbData.toArray4d();
			probas = probData.toArray4d();
			lb0  = new byte[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				lb0[x][y][z] = labels[x][y][z][0];	
			}
		} else {
			lb0 = lbData.toArray3d();
			float[][][] tmp = probData.toArray3d();
			labels  = new byte[nx][ny][nz][1];
			probas  = new float[nx][ny][nz][1];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				labels[x][y][z][0] = lb0[x][y][z];
				probas[x][y][z][0] = tmp[x][y][z];
			}
			tmp = null;
		}
		probData.dispose();
		maxVolParam.dispose();
		lbData.dispose();
		lbVolParam.dispose();
		
		byte[] lblist = ObjectLabeling.listOrderedLabels(lb0, nx, ny, nz);
		int nlb = lblist.length;
		lb0 = null;
		
		int idx = indexParam.getValue().intValue();
		
		for (int lb=0;lb<nlb;lb++) if (idx==-1 || idx==lblist[lb]) {
			float[][][] proba = new float[nx][ny][nz];
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				proba[x][y][z] = 0.0f;
				for (int c=0;c<nc;c++) {
					if (labels[x][y][z][c]==lblist[lb]) proba[x][y][z] = probas[x][y][z][c];
				}
			}
			ImageDataFloat res = new ImageDataFloat(proba);
			res.setHeader(maxVolParam.getImageData().getHeader());
			res.setName(maxVolParam.getImageData().getName()+"_p"+lblist[lb]);
			
			volParam.add(res);
			volParam.writeAndFreeNow(this);
		}
	}
}
