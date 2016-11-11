package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
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
import org.apache.commons.math3.random.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistIntensityRobustNoiseEstimation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	private ParamVolume maskImage;
	private ParamOption ngbParam;
	private ParamBoolean outlierParam;
	//private ParamOption distParam;
	private ParamOption modelParam;
	private ParamOption samplingParam;
	
	private ParamVolume edgeImage;
	private ParamVolume outImage;
	private ParamVolume countImage;
	private ParamFloat  medianParam;
	
	// parameters
	private		static final String[]	ngbTypes = {"6C","18C","26C"};
	private		String		ngbType = "6C";
	
	//private		static final String[]	distTypes = {"distance", "difference"};
	//private		String		distType = "distance";
	
	private		static final String[]	samplingTypes = {"all", "median", "random"};
	private		String		samplingType = "median";
	
	private		static final String[]	modelTypes = {"Gaussian", "Exponential"};
	private		String		modelType = "Gaussian";
	
	private static final byte[] ngbx = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	private static final byte[] ngby = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	private static final byte[] ngbz = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(maskImage = new ParamVolume("Mask Image (opt)"));
		maskImage.setMandatory(false);
		
		inputParams.add(modelParam = new ParamOption("Noise model", modelTypes));
		modelParam.setValue(modelType);
		
		//inputParams.add(distParam = new ParamOption("Measure", distTypes));
		//distParam.setValue(distType);
		
		inputParams.add(ngbParam = new ParamOption("Neighborhood connectivity", ngbTypes));
		ngbParam.setValue(ngbType);
		
		inputParams.add(samplingParam = new ParamOption("Neighborhood sampling", samplingTypes));
		ngbParam.setValue(samplingType);
		
		inputParams.add(outlierParam = new ParamBoolean("Estimate outliers", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.devel");
		inputParams.setLabel("Robust Noise Estimation");
		inputParams.setName("RobustNoiseEstimation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Noise estimation in images based on median statistics with outliers.");
		
		info.setVersion("3.0.6");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outImage = new ParamVolume("Outlier Image",VoxelType.FLOAT));
		outputParams.add(medianParam = new ParamFloat(" Estimated noise parameter"));
		
		outputParams.setName("robust noise images");
		outputParams.setLabel("robust noise images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataFloat inImg = new ImageDataFloat(inputImage.getImageData());
		float[][][] image = inImg.toArray3d();
		int nx = inImg.getRows();
		int ny = inImg.getCols();
		int nz = inImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = inImg.getHeader().getDimResolutions()[0];
		float ry = inImg.getHeader().getDimResolutions()[1];
		float rz = inImg.getHeader().getDimResolutions()[2];
		inImg = null;
		
		// basic mask for zero values
		boolean[][][] mask = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]==0) mask[x][y][z] = false;
			else mask[x][y][z]=true;
		}
		
		// find min,max for scale
		float INF = 1e9f;
		float Imin = INF, Imax = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
		}
		
		// main algorithm
		
		// first estimate the noise level
		int ngb = 6; // 6C neighborhood
		if (ngbParam.getValue().equals("18C")) ngb = 18;
		else if (ngbParam.getValue().equals("26C")) ngb = 26;
		
		// for median sampling
		double[] neighbors = new double[ngb];
		Percentile measure = new Percentile();
		// for random sampling
		Well19937c random = new Well19937c();
		
		// 1. generate the data
		int nb = 0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) nb++;
		int nngb = 1;
		if (samplingParam.getValue().equals("all")) nngb = ngb;
		
		double data[] = new double[nngb*nb];
		int n = 0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			if (samplingParam.getValue().equals("all")) {
				// all data
				for (int d=0;d<nngb;d++) {
					data[n] = Numerics.abs(image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]);
				 	n++;
				}
			} else if (samplingParam.getValue().equals("median")) {
				// sample the median of the neighbors
				for (int d=0;d<nngb;d++) {
					neighbors[d] = Numerics.abs(image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]);
				}
				data[n] = (float)measure.evaluate(neighbors,50.0);
				n++;
			} else if (samplingParam.getValue().equals("random")) {
				// sample a random neighbor
				int d = random.nextInt(ngb);
				data[n] = Numerics.abs(image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]);
				n++;
			}
		}
		
		float sigma = 0.0f;
		float beta = 0.0f;
		if (modelParam.getValue().equals("Gaussian")) {
			sigma = ImageStatistics.robustHalfGaussianFit(data, outlierParam.getValue().booleanValue(), nngb*nb);
		} else if (modelParam.getValue().equals("Exponential")) {
			beta = ImageStatistics.robustExponentialFit(data, outlierParam.getValue().booleanValue(), nngb*nb);
		}
		float[][][] outliers = new float[nx][ny][nz];
		double outratio = 0.0;
		double sum = 0.0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			outliers[x][y][z] = 0.0f;
			for (int d=0;d<ngb;d++) {
				 if (mask[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]) {
				 	float diff = image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]];
				 	
				 	float proba = 0.0f;
				 	if (modelParam.getValue().equals("Gaussian")) {
				 		proba = (float)(FastMath.exp( -0.5*diff*diff/(sigma*sigma))*FastMath.sqrt(2.0/(FastMath.PI*sigma*sigma)) );
				 	} else if (modelParam.getValue().equals("Exponential")) {
				 		proba = (float)(FastMath.exp( -diff/beta)/beta );
				 	}				 		
				 	float ratio = 1.0f/(Imax-Imin)/(1.0f/(Imax-Imin)+proba);
				 	outratio += ratio;
				 	sum++;
				 	if (ratio>outliers[x][y][z]) outliers[x][y][z] = ratio;
				}
			}
		}
		outratio /= sum;
		
		System.out.println("outlier ratio: "+outratio);
		
		// output
		if (modelParam.getValue().equals("Gaussian")) {
			medianParam.setValue(sigma);
		} else if (modelParam.getValue().equals("Exponential")) {
			medianParam.setValue(beta);
		}
		
		ImageDataFloat bufferData = new ImageDataFloat(outliers);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_out");
		outImage.setValue(bufferData);
		bufferData = null;
		outliers = null;
	}


}
