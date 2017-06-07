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

import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.special.Erf;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensity4DCombination extends ProcessingAlgorithm{
	ParamVolume volParam;
	ParamVolume weightParam;
	ParamVolume resultVolParam;
	ParamOption operation;

	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(volParam=new ParamVolume("4D Volume"));
		inputParams.add(weightParam=new ParamVolume("4D Weight (opt)"));
		weightParam.setMandatory(false);
		inputParams.add(operation=new ParamOption("Operation",new String[]{"average","median","maximum","minimum","weighted_avg","max_weight","stdev","robust_dev"}));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("4D Image Combination");
		inputParams.setName("4DImageCombination");


		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription("combines the last dimension of a 4D image stack (different options are possible: simple average, median, maximum or minimum value, weighted average or value for the maximum weight).");
		info.setVersion("3.0.9");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultVolParam=new ParamVolume("Result Volume",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageData vol=volParam.getImageData();
		ImageData weight=weightParam.getImageData();
		int rows=vol.getRows();
		int cols=vol.getCols();
		int slices=vol.getSlices();
		int comps = vol.getComponents();
		
		double dev = 100.0*0.5*Erf.erf(1.0/FastMath.sqrt(2.0));
		
		String resultName = vol.getName();
		if (operation.getValue().equals("average")) {
			resultName+= "_avg";
		}  else
		if (operation.getValue().equals("median")) {
			resultName+= "_med";
		}  else
		if (operation.getValue().equals("maximum")) {
			resultName+= "_max";
		} else
		if (operation.getValue().equals("minimum")) {
			resultName+= "_min";
		} else
		if (operation.getValue().equals("weighted_avg")) {
			resultName+= "_wavg";
		} else
		if (operation.getValue().equals("max_weight")) {
			resultName+= "_maxw";
		} else
		if (operation.getValue().equals("stdev")) {
			resultName+= "_std";
		} else
		if (operation.getValue().equals("robust_dev")) {
			resultName+= "_rstd";
		}
		ImageDataFloat resultVol=new ImageDataFloat(resultName,rows,cols,slices);
		
		double tmp;
		double sum,min,max;
		double w,wsum,wden,wmax;
		int maxwid;
		double[] vals = new double[comps];
		Percentile measure = new Percentile();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < slices; k++) {
					sum=0.0; wsum=0.0; wden=0.0;
					min=1e12; max=-1e12;
					wmax=-1e12;maxwid=-1;
					for (int l = 0; l < comps; l++) {
						tmp=vol.getDouble(i, j, k, l);
						vals[l] = tmp;
						sum+=tmp;
						if (tmp>max) max = tmp;
						if (tmp<min) min = tmp;
						
						if (weight!=null) {
							w=weight.getDouble(i, j, k, l);
							wsum += w*sum;
							wden += w;
							if (w>wmax) {
								wmax = w;
								maxwid = l;
							}
						}
					}
					if (operation.getValue().equals("average")) {
						resultVol.set(i, j, k, sum/comps);
					} else if (operation.getValue().equals("median")) {
						resultVol.set(i, j, k, measure.evaluate(vals, 0, comps, 50.0));
					} else if (operation.getValue().equals("maximum")) {
						resultVol.set(i, j, k, max);
					} else if (operation.getValue().equals("minimum")) {
						resultVol.set(i, j, k, min);
					} else if (operation.getValue().equals("weighted_avg")) {
						resultVol.set(i, j, k, wsum/wden);
					} else if (operation.getValue().equals("max_weight")) {
						resultVol.set(i, j, k, vals[maxwid]);
					} else if (operation.getValue().equals("stdev")) {
						double var=0.0;
						double mean=sum/comps;
						for (int l = 0; l < comps; l++) {
							tmp=vol.getDouble(i, j, k, l);
							var+=(tmp-mean)*(tmp-mean);
						}
						resultVol.set(i,j,k, FastMath.sqrt(var/comps));
					}  else if (operation.getValue().equals("robust_dev")) {
						double fdev = 0.5*(measure.evaluate(vals, 0, comps, 50.0+dev) - measure.evaluate(vals, 0, comps, 50.0-dev));
						resultVol.set(i,j,k, fdev);
					}
				}
			}
		}
		resultVol.setHeader(vol.getHeader());
		resultVolParam.setValue(resultVol);
	}
}
