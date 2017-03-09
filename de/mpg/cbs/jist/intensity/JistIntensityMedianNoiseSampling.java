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
public class JistIntensityMedianNoiseSampling extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume inputImage;
	private ParamOption ngbParam;
	private ParamFloat ratioParam;
	private ParamInteger sampleParam;
	private ParamInteger iterParam;
	private ParamBoolean adjustParam;
	private ParamOption distParam;
	
	private ParamVolume edgeImage;
	private ParamVolume outImage;
	private ParamVolume countImage;
	private ParamFloat  medianParam;
	
	// parameters
	private		static final String[]	ngbTypes = {"6C","18C","26C"};
	private		String		ngbType = "6C";
	
	private		static final String[]	distTypes = {"distance", "difference"};
	private		String		distType = "distance";
	
	private static final byte[] ngbx = {+1,  0,  0, -1,  0,  0, +1,  0, +1, +1,  0, -1, -1,  0, +1, -1,  0, -1, +1, +1, +1, +1, -1, -1, -1, -1};
	private static final byte[] ngby = { 0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1,  0, +1, -1, -1, +1, +1, -1, -1, +1};
	private static final byte[] ngbz = { 0,  0, +1,  0,  0, -1,  0, +1, +1,  0, -1, +1,  0, +1, -1,  0, -1, -1, +1, -1, +1, -1, +1, -1, +1, -1};
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		
		inputParams.add(sampleParam = new ParamInteger("Sample size", 10, 1000000, 200));
		inputParams.add(iterParam = new ParamInteger("Max iterations", 0, 1000000, 200));
		inputParams.add(distParam = new ParamOption("Measure", distTypes));
		distParam.setValue(distType);
		
		inputParams.add(ngbParam = new ParamOption("Neighborhood connectivity", ngbTypes));
		ngbParam.setValue(ngbType);
		
		//inputParams.add(ratioParam = new ParamFloat("Scaling ratio", 0.0f, 1.0f, 0.0f));
		inputParams.add(adjustParam = new ParamBoolean("Two-level estimation", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity.devel");
		inputParams.setLabel("Median Noise Sampling");
		inputParams.setName("MedianNoiseSampling");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Noise estimation in images based on median sampling.");
		
		info.setVersion("3.0.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(edgeImage = new ParamVolume("Edge Image",VoxelType.FLOAT));
		outputParams.add(outImage = new ParamVolume("Outlier Image",VoxelType.FLOAT));
		outputParams.add(countImage = new ParamVolume("Sampling Image",VoxelType.FLOAT));
		outputParams.add(medianParam = new ParamFloat("Median noise"));
		
		outputParams.setName("median noise images");
		outputParams.setLabel("median noise images");
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
		
		// find min,max for scaling
		float INF = 1e9f;
		float Imin = INF, Imax = -INF;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]<Imin) Imin = image[x][y][z];
			if (image[x][y][z]>Imax) Imax = image[x][y][z];
		}
		/*
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			image[x][y][z] = (image[x][y][z]-Imin)/(Imax-Imin);
		}
		*/
		
		// main algorithm
		
		// first estimate the noise level
		int ngb = 6; // 6C neighborhood
		if (ngbParam.getValue().equals("18C")) ngb = 18;
		else if (ngbParam.getValue().equals("26C")) ngb = 26;
		
		int Nsample = sampleParam.getValue().intValue();
		int Niter = iterParam.getValue().intValue();
		
		float[] median = new float[Niter];
		float[] range = new float[Niter];
		Well19937c rand = new Well19937c();
		
		boolean abs = false;
		if (distParam.getValue().equals("distance")) abs = true;
		
		// sample uniformly location and difference direction
		float dist = 1e9f;
		float meanmed = 0.0f;
		float meanrng = 0.0f;
		float nsample = 0.0f;
		float[][][] sampling = new float[nx][ny][nz];
		int maxtrial = (nx*ny+ny*nz+nz*nx)*ngb;
		float mindist = 0.01f;
		boolean stop = false;
		int nstop = 0;
		for (int t=0;t<Niter && !stop;t++) {
			// build a sample array
			double[] sample = new double[Nsample];
			int count=0;
			int trial=0;
			while (count<Nsample && trial<maxtrial) {
				trial++;
				// generate pseudo-random locations
				int x = 1+rand.nextInt(nx-1);		
				int y = 1+rand.nextInt(ny-1);		
				int z = 1+rand.nextInt(nz-1);		
				if (mask[x][y][z]) {
					// pseudo-random direction
					int d = rand.nextInt(ngb);
					if (mask[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]) {
						if (abs) sample[count] = Numerics.abs(image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]);
						else sample[count] = image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]];
						count++;
						sampling[x][y][z]++;
					}
				}
			}
			if (trial>=maxtrial) System.out.println("too many failed trials!! (succeeded: "+count+")");
			
			// estimate mean, variance robustly
			Percentile measure = new Percentile();
			measure.setData(sample);
			median[t] = (float)measure.evaluate(50.0);
			//range[t] = (float)(measure.evaluate(84.1) - measure.evaluate(15.9));
			range[t] = (float)(measure.evaluate(69.15) - measure.evaluate(30.85));
			
			// check for convergence with the mean (rather than median that has non-zero proba to be zero randomly)
			if (t>FastMath.sqrt(Niter)) {
				dist = Numerics.abs(meanmed/(float)t - (meanmed+median[t])/(float)(t+1));
				dist = Numerics.max(dist,Numerics.abs(meanrng/(float)t - (meanrng+range[t])/(float)(t+1)));
				if (dist<mindist*0.01f*(Imax-Imin)) {
					nstop++;
					if (nstop>=5) stop = true;
				} else {
					nstop=0;
				}
			}
			meanmed += median[t];
			meanrng += range[t];
			nsample++;
			
			//System.out.println((t+1)+". sample (median, range): "+median[t]+", "+range[t]+" --> estimate: "+(meanmed/nsample)+", "+(meanrng/nsample));
		}
		// also check the median^2?
		double[] med = new double[(int)nsample];
		double[] rng = new double[(int)nsample];
		for (int t=0;t<nsample;t++) {
			med[t] = median[t];
			rng[t] = range[t];
		}
		Percentile measure = new Percentile();
		measure.setData(med);
		float medmed = (float)measure.evaluate(50.0);
		measure.setData(rng);
		float medrng = (float)measure.evaluate(50.0);
		System.out.println("median (median, range): "+medmed+", "+medrng);
		
		float[][][] outliers = new float[nx][ny][nz];
		double outratio = 0.0;
		double sum = 0.0;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
			outliers[x][y][z] = 0.0f;
			for (int d=0;d<ngb;d++) {
				 if (mask[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]) {
				 	float diff = image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]];
				 	if (abs) diff = Numerics.abs(diff);
				 	//float proba = (float)(FastMath.exp( -0.5*Numerics.square((diff-medmed)/(0.5*medrng)))
				 	//					  /FastMath.sqrt(2.0*FastMath.PI*0.25*medrng*medrng) );
				 	float proba = (float)(FastMath.exp( -0.5*Numerics.square((diff-medmed)/(medrng)))
				 						  /FastMath.sqrt(2.0*FastMath.PI*medrng*medrng) );
				 	float ratio = 1.0f/(Imax-Imin)/(1.0f/(Imax-Imin)+proba);
				 	outratio += ratio;
				 	sum++;
				 	if (ratio>outliers[x][y][z]) outliers[x][y][z] = ratio;
				}
			}
		}
		outratio /= sum;
		
		System.out.println("outlier ratio: "+outratio);
		
		// second pass?
		float[][][] edges = null;
		if (adjustParam.getValue().booleanValue()) {

			for (int n=0;n<20;n++) {
				double prevratio = outratio;
				double volundercurve = (0.1915*(1-outratio) + outratio*medrng/(Imax-Imin))*100.0;
				//double volundercurve = (0.5*(0.682*(1-outratio) + outratio*medrng/(Imax-Imin)))*100.0;
				
				System.out.println("half volume under the [mu-sigma,mu+sigma] curves (%): "+volundercurve);
		
				// sample uniformly location and difference direction
				dist = 1e9f;
				meanmed = 0.0f;
				meanrng = 0.0f;
				nsample = 0.0f;
				maxtrial = (nx*ny+ny*nz+nz*nx)*ngb;
				mindist = 0.01f;
				stop = false;
				nstop = 0;
				for (int t=0;t<Niter && !stop;t++) {
					// build a sample array
					double[] sample = new double[Nsample];
					int count=0;
					int trial=0;
					while (count<Nsample && trial<maxtrial) {
						trial++;
						// generate pseudo-random locations
						int x = 1+rand.nextInt(nx-1);		
						int y = 1+rand.nextInt(ny-1);		
						int z = 1+rand.nextInt(nz-1);		
						if (mask[x][y][z]) {
							// pseudo-random direction
							int d = rand.nextInt(ngb);
							if (mask[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]) {
								if (abs) sample[count] = Numerics.abs(image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]);
								else sample[count] = image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]];
								count++;
								sampling[x][y][z]++;
							}
						}
					}
					if (trial>=maxtrial) System.out.println("too many failed trials!! (succeeded: "+count+")");
					
					// estimate mean, variance robustly
					measure = new Percentile();
					measure.setData(sample);
					median[t] = (float)measure.evaluate(50.0);
					range[t] = (float)(measure.evaluate(50.0+volundercurve) - measure.evaluate(50.0-volundercurve));
					
					// check for convergence with the mean (rather than median that has non-zero proba to be zero randomly)
					if (t>FastMath.sqrt(Niter)) {
						dist = Numerics.abs(meanmed/(float)t - (meanmed+median[t])/(float)(t+1));
						dist = Numerics.max(dist,Numerics.abs(meanrng/(float)t - (meanrng+range[t])/(float)(t+1)));
						if (dist<mindist*0.01f*(Imax-Imin)) {
							nstop++;
							if (nstop>=5) stop = true;
						} else {
							nstop=0;
						}
					}
					meanmed += median[t];
					meanrng += range[t];
					nsample++;
					
					//System.out.println((t+1)+". sample (median, range): "+median[t]+", "+range[t]+" --> estimate: "+(meanmed/nsample)+", "+(meanrng/nsample));
				}
				// also check the median^2?
				med = new double[(int)nsample];
				rng = new double[(int)nsample];
				for (int t=0;t<nsample;t++) {
					med[t] = median[t];
					rng[t] = range[t];
				}
				measure = new Percentile();
				measure.setData(med);
				medmed = (float)measure.evaluate(50.0);
				measure.setData(rng);
				medrng = (float)measure.evaluate(50.0);
				System.out.println("median (median, range): "+medmed+", "+medrng);
			
				edges = new float[nx][ny][nz];
				outratio = 0.0f;
				sum = 0.0f;
				for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (mask[x][y][z]) {
					edges[x][y][z] = 0.0f;
					for (int d=0;d<ngb;d++) {
						 if (mask[x+ngbx[d]][y+ngby[d]][z+ngbz[d]]) {
							float diff = image[x][y][z]-image[x+ngbx[d]][y+ngby[d]][z+ngbz[d]];
							if (abs) diff = Numerics.abs(diff);
							//float proba = (float)(FastMath.exp( -0.5*Numerics.square((diff-medmed)/(0.5*medrng)))
							//					  /FastMath.sqrt(2.0*FastMath.PI*0.25*medrng*medrng) );
							float proba = (float)(FastMath.exp( -0.5*Numerics.square((diff-medmed)/(medrng)))
												  /FastMath.sqrt(2.0*FastMath.PI*medrng*medrng) );
							float ratio = 1.0f/(Imax-Imin)/(1.0f/(Imax-Imin)+proba);
							outratio += ratio;
							sum++;
							if (ratio>edges[x][y][z])
								edges[x][y][z] = ratio;
						}
					}
				}
				outratio /= sum;
				
				System.out.println("outlier ratio: "+outratio);
				
				if (Numerics.abs(prevratio-outratio)<0.001) n = 10000;
			}
		}
		
		// output
		medianParam.setValue(medmed);
		
		ImageDataFloat bufferData = new ImageDataFloat(edges);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_edges");
		edgeImage.setValue(bufferData);
		bufferData = null;
		edges = null;
				
		bufferData = new ImageDataFloat(outliers);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_out");
		outImage.setValue(bufferData);
		bufferData = null;
		outliers = null;
		
		bufferData = new ImageDataFloat(sampling);
		bufferData.setHeader(inputImage.getImageData().getHeader());
		bufferData.setName(inputImage.getImageData().getName()+"_count");
		countImage.setValue(bufferData);
		bufferData = null;
		sampling = null;
	}


}
