package de.mpg.cbs.jist.filter;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.special.Erf;


/*
 * @author Pierre-Louis Bazin
 */
public class JistFilterNormalizeResponse extends ProcessingAlgorithm {

	// parameters
	private		static final String[]	probatypes = {"Gaussian","Exponential","Exp+Log-normal"};
	private		String		probatype = "Gaussian";
	private		static final String[]	measuretypes = {"median","mean","mean-iterated","median-iterated"};
	private		String		measuretype = "median";
	
	// jist containers
	private ParamVolume inputImage;
	private ParamVolume resultImage;
	private ParamOption probaParam;
	private ParamOption measureParam;
	//private ParamBoolean skipParam;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputImage = new ParamVolume("Input Image"));
		inputParams.add(probaParam = new ParamOption("Null distribution", probatypes));
		probaParam.setDescription("Probabilistic model for the null response of the filter.");
		probaParam.setValue(probatype);
		inputParams.add(measureParam = new ParamOption("Measured statistic", measuretypes));
		measureParam.setDescription("Statistic and algorithm used to estimate the distribution");
		measureParam.setValue(measuretype);
		//inputParams.add(skipParam = new ParamBoolean("skip zero values", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Filter");
		inputParams.setLabel("Normalize Response");
		inputParams.setName("NormalizeResponse");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Maps a filter response to a probability based on a mixture of null distribution and uniform outliers.");
		
		info.setVersion("3.0.5");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultImage = new ParamVolume("Probability Image",VoxelType.FLOAT));
		outputParams.setName("proba_image");
		outputParams.setLabel("Proba Image");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat in = new ImageDataFloat(inputImage.getImageData());
		float[][][] image = in.toArray3d();
		int nx = in.getRows();
		int ny= in.getCols();
		int nz = in.getSlices();
		
		String model = probaParam.getValue();
		String metric = measureParam.getValue();

		// main algorithm
		
		// first pass: non-robust statistics
		double mean = 0.0;
		double sqmean = 0.0;
		double max = 0.0;
		int nb = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
			
			mean += image[x][y][z];
			sqmean += image[x][y][z]*image[x][y][z];
			max = Numerics.max(max,image[x][y][z]);
			nb++;
		}
		mean /= nb;
		sqmean /= nb;
		
		BasicInfo.displayMessage("parameter estimates: mean "+mean+", sq mean "+FastMath.sqrt(sqmean)+"\n");
		
		// robust measures
		double[] response = new double[nb];
		int n=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
			response[n] = image[x][y][z];
			n++;
		}
		Percentile measure = new Percentile();
		double median = measure.evaluate(response, 50.0);
		measure = null;
		response = null;
		
		BasicInfo.displayMessage("parameter estimates: median "+median+"\n");
		
		// exponential fit or gaussian fit ? gaussian seems more sensitive, too much for the strict filter
		
		float[][][] result = new float[nx][ny][nz];

		// exponential estimates		
		double beta = mean;
		if (metric.startsWith("median")) beta = median/FastMath.log(2.0);
		
		BasicInfo.displayMessage("parameter estimates: beta "+beta+"\n");
		
		// half-normal estimates
		double sigma2 = sqmean;
		// med = sigma x sqrt(2) x erf-1(1/2)
		if (metric.startsWith("median")) sigma2 = median*median/0.45493642;
		double normg = FastMath.sqrt(2.0/(sigma2*FastMath.PI));
		
		BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
		
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
			if (model.startsWith("Exp")) {
				double pe = FastMath.exp( -image[x][y][z]/beta)/beta;
				result[x][y][z] = (float)(1.0/max/( 1.0/max+pe));
			} else if (model.equals("Gaussian")) {
				double pg = normg*FastMath.exp(-image[x][y][z]*image[x][y][z]/(2.0*sigma2));
				result[x][y][z] = (float)(1.0/max/( 1.0/max+pg));
			}
		}
		
		if (metric.endsWith("iterated")) {
			double betaprev;
			double sigma2prev;
			for (int t=0;t<20;t++) {
				BasicInfo.displayMessage("iteration "+(t+1)+"\n");

				betaprev = beta;
				sigma2prev = sigma2;
				
				mean = 0.0;
				sqmean = 0.0;
				nb = 0;
				double den = 0.0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
					
					mean += result[x][y][z]*image[x][y][z];
					sqmean += result[x][y][z]*image[x][y][z]*image[x][y][z];
					den += result[x][y][z];
					if (result[x][y][z]<0.5) nb++;
				}
				mean /= den;
				sqmean /= den;
				
				BasicInfo.displayMessage("parameter estimates: mean "+mean+", sq mean "+FastMath.sqrt(sqmean)+"\n");
				
				// robust measures
				response = new double[nb];
				n=0;
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
					if (result[x][y][z]<0.5) {
						response[n] = image[x][y][z];
						n++;
					}
				}
				measure = new Percentile();
				median = measure.evaluate(response, 50.0);
				measure = null;
				response = null;
				
				BasicInfo.displayMessage("parameter estimates: median "+median+"\n");

				// exponential estimates		
				beta = mean;
				if (metric.startsWith("median")) beta = median/FastMath.log(2.0);
				
				BasicInfo.displayMessage("parameter estimates: beta "+beta+"\n");
				
				// half-normal estimates
				sigma2 = sqmean;
				// med = sigma x sqrt(2) x erf-1(1/2)
				if (metric.startsWith("median")) sigma2 = median*median/0.45493642;
				normg = FastMath.sqrt(2.0/(sigma2*FastMath.PI));
				
				BasicInfo.displayMessage("parameter estimates: sigma "+FastMath.sqrt(sigma2)+"\n");
				
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
					if (model.startsWith("Exp")) {
						double pe = FastMath.exp( -image[x][y][z]/beta)/beta;
						result[x][y][z] = (float)(1.0/max/( 1.0/max+pe));
					} else if (model.equals("Gaussian")) {
						double pg = normg*FastMath.exp(-image[x][y][z]*image[x][y][z]/(2.0*sigma2));
						result[x][y][z] = (float)(1.0/max/( 1.0/max+pg));
					}
				}
				if (model.startsWith("Exp") && 2.0*Numerics.abs(beta-betaprev)/(beta+betaprev)<0.01) t=1000;
				if (model.equals("Gaussian") && 2.0*Numerics.abs(FastMath.sqrt(sigma2)-FastMath.sqrt(sigma2prev))
												/(FastMath.sqrt(sigma2)+FastMath.sqrt(sigma2prev))<0.01) t=1000;
			} 

		}
		
		if (model.equals("Exp+Log-normal")) {
			
			// model the filter response as something more interesting, e.g. log-normal (removing the bg samples)
			nb=0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
				nb++;
			}
			response = new double[nb];
			n=0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (image[x][y][z]>0) {
				response[n] = image[x][y][z];
				n++;
			}
			double[] weights = new double[nb];
			for (int b=0;b<nb;b++) { 
				weights[b] = (1.0-FastMath.exp( -response[b]/beta));
				response[b] = FastMath.log(response[b]);
			}
			
			double fmean = ImageStatistics.weightedPercentile(response,weights,50.0,nb);
			
			// stdev: 50% +/- 1/2*erf(1/sqrt(2)) (~0.341344746..)
			double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
			double fdev = 0.5*(ImageStatistics.weightedPercentile(response,weights,50.0+dev,nb) - ImageStatistics.weightedPercentile(response,weights,50.0-dev,nb));
			
			BasicInfo.displayMessage("Log-normal parameter estimates: mean = "+FastMath.exp(fmean)+", stdev = "+FastMath.exp(fdev)+",\n");
			
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++){
				int xyz = x + nx*y + nx*ny*z;
				if (image[x][y][z]>0) {
					double pe = FastMath.exp( -image[x][y][z]/beta)/beta;
					double plg = FastMath.exp(-Numerics.square(FastMath.log(image[x][y][z])-fmean)/(2.0*fdev*fdev))/FastMath.sqrt(2.0*FastMath.PI*fdev*fdev);
					result[x][y][z] = (float)(plg/(plg+pe));
				}
			}
			
		}
		
		
		String suffix = "_";
		if (model.startsWith("Exp")) suffix += "e";
		else suffix += "g";
		if (model.endsWith("Log-normal")) suffix += "l";
		if (metric.startsWith("median")) suffix += "r";
		else suffix += "m";
		if (metric.endsWith("iterated")) suffix += "i";
		else suffix += "s";
		
		ImageDataFloat resultData = new ImageDataFloat(result);		
		resultData.setHeader(in.getHeader());
		resultData.setName(in.getName()+suffix);
		resultImage.setValue(resultData);
		resultData = null;
		result = null;
		
	}

}
