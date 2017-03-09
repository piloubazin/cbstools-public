package de.mpg.cbs.jist.fmri;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.special.*;

import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.ModelStorageBase;

/*
 * @author Pierre-Louis Bazin
 */
public class JistFmriNoiseEstimate extends ProcessingAlgorithm {

	// jist containers
	
	// input 
	private ParamVolume 	dataImage;
	
	private static final byte X=0;
	private static final byte Y=1;
	private static final byte Z=2;
	
	// ouput
	private ParamVolume noiseImage;
	private ParamFloat noiseParam;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(dataImage = new ParamVolume("fMRI Data Image",VoxelType.FLOAT,-1,-1,-1,-1));
		dataImage.setLoadAndSaveOnValidate(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("fMRI");
		inputParams.setLabel("fMRI Noise Estimation");
		inputParams.setName("fMRINoiseEstimation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Estimates locally-varying noise as the median of the difference between consecutive time points.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(noiseImage = new ParamVolume("Noise Map",VoxelType.FLOAT));
		outputParams.add(noiseParam = new ParamFloat("Median Noise Level"));
		
		outputParams.setLabel("fMRI Noise Estimation");
		outputParams.setName("fMRINoiseEstimation");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat	dataImg = new ImageDataFloat(dataImage.getImageData());
		int nx = dataImg.getRows();
		int ny = dataImg.getCols();
		int nz = dataImg.getSlices();
		int nd = dataImg.getComponents();
		
		String dataName = dataImg.getName();
		ImageHeader dataHeader = dataImg.getHeader();
		
		float[][][][] data = dataImg.toArray4d();
		dataImg.dispose();
		dataImage.dispose();
		
		System.out.println("fmri data loaded");

		// local noise estimate: median of consecutive time point distance
		float[][][] noise = new float[nx][ny][nz];
		double[] sample = new double[nd-1];
		// median (half-normal) = sigma * erfinv(1/2) [0.476936] * sqrt(2) [1.4142136] : ~1.3 (sigma is x2)
		double factor = 1.0/Erf.erfInv(0.5)/FastMath.sqrt(2.0)/2.0;
		int ndata=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (data[x][y][z][0]>0) {
			 for (int d=1;d<nd;d++) sample[d-1] = Numerics.abs(data[x][y][z][d]-data[x][y][z][d-1]);
			 Percentile measure = new Percentile();
			 measure.setData(sample);
			 noise[x][y][z] = (float)(factor*measure.evaluate(50.0));
			 ndata++;
		}
		
		// global estimate: median
		double[] med = new double[ndata];
		int n=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (data[x][y][z][0]>0) {
			med[n] = noise[x][y][z];
			n++;
		}
		Percentile measure = new Percentile();
		measure.setData(med);
		float mednoise = (float)measure.evaluate(50.0);
		
		System.out.println("output..");
		
		// output
		ImageDataFloat noiseImg = new ImageDataFloat(noise);		
		noiseImg.setHeader(dataHeader);
		noiseImg.setName(dataName + "_noise");
		noiseImage.setValue(noiseImg);
		
		noiseParam.setValue(mednoise);
	}

}
