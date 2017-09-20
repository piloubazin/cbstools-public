package de.mpg.cbs.jist.statistics;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile.DialogType;

import java.net.*;
import java.io.*;
import java.util.*;
import java.lang.*;
import java.text.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistStatisticsSurfaceTimeSeries extends ProcessingAlgorithm {

	// jist containers
	private ParamSurface 	dataSurface;
	private ParamSurface 	labelSurface;
	private ParamFile 		atlasParam;
	private ParamBoolean	removeFirstParam;
	
	private ParamFile 		statsParam;
	
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private static final String[] statTypes = {"Mean", "Stdev", "none"};
	
	private ParamFile 		outputParam;
		
	private String delim = ",";
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(dataSurface = new ParamSurface("Surface with Time series data"));
		inputParams.add(labelSurface = new ParamSurface("Surface with ROI labels"));
		
		inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		atlasParam.setMandatory(false);
		inputParams.add(removeFirstParam = new ParamBoolean("skip first label",true));
		
		inputParams.add(statsParam = new ParamFile("Spreadsheet file directory", DialogType.DIRECTORY));
		
		inputParams.add(stat1Param = new ParamOption("Statistic 1", statTypes));
		inputParams.add(stat2Param = new ParamOption("Statistic 2", statTypes));
		inputParams.add(stat3Param = new ParamOption("Statistic 3", statTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics");
		inputParams.setLabel("Time Series Statistics (surface)");
		inputParams.setName("SurfaceTimeSeriesStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various statistics from time series and volumetric labels and outputs to a spreadsheet.");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputParam=new ParamFile("statistics spreadsheet"));
		
		outputParams.setName("spreadsheet result");
		outputParams.setLabel("spreadsheet result");
	}

	@Override
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surfdata = dataSurface.getSurface();
		int nv = surfdata.getVertexCount();
		String name = surfdata.getName();
		double[][] data = surfdata.getVertexData();
		int nt = surfdata.getVertexData(0).length;
		surfdata = null;
		
		BasicInfo.displayMessage("n vertex = "+nv+"\n");
		BasicInfo.displayMessage("n time series = "+nt+"\n");
		
		// get label list from atlas or label map
		int nlabels;
		String[] lbname;
		int[] lbid;
		
		int[] labelmap = new int[nv];
		double[][] surflabels = labelSurface.getSurface().getVertexData();
		for(int v=0; v<nv; v++) {
			labelmap[v] = (int)surflabels[v][0];
		}
		surflabels = null;
		
		if (atlasParam.getValue()!=null && atlasParam.getValue().length()>0) {
			BasicInfo.displayMessage("Load atlas\n");
	
			SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam.getValue().getAbsolutePath());
			
			nlabels = atlas.getNumber();
			lbname = atlas.getNames();
			byte[] tmp = atlas.getLabels();
			lbid = new int[nlabels];
			for (int n=0;n<nlabels;n++) lbid[n] = tmp[n];
			tmp = null;
		} else {
			// find the number of labels
			lbid = ObjectLabeling.listOrderedLabels(labelmap, nv, 1, 1);
			nlabels = lbid.length;
			lbname = new String[nlabels];
			for (int n=0;n<nlabels;n++) lbname[n] = new String("Label_"+lbid[n]);
		}
		
		
		// Main algorithm
		
		String line;
		
		// standard computations for the object
		ArrayList<String> output = new ArrayList();

		// remove the zero label, if needed
		if (removeFirstParam.getValue().booleanValue()) {
			int[] tmp = new int[nlabels-1];
			String[] tmpname = new String[nlabels-1];
			for (int n=1;n<nlabels;n++) {
				tmp[n-1] = lbid[n];
				tmpname[n-1] = lbname[n];
			}
			nlabels--;
			lbid = tmp;
			lbname = tmpname;
		}
		
		
		ArrayList<String> statistics = new ArrayList();
		statistics.add(stat1Param.getValue());
		statistics.add(stat2Param.getValue());
		statistics.add(stat3Param.getValue());
		
		// pre-compute redundant measures ?
		
		// compute the statistics
		for (int s=0; s<statistics.size(); s++) {
			System.out.print("Statistic: "+statistics.get(s)+"\n");
			if (statistics.get(s).equals("Mean")) {
				float[][] mean = new float[nlabels][nt];
				float[] den = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					for (int t=0;t<nt;t++) mean[n][t] = 0.0f;
					den[n] = 0.0f;
				}
				for (int v=0;v<nv;v++) {
					for (int n=0;n<nlabels;n++) if (labelmap[v]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += data[v][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int n=0;n<nlabels;n++) {
					line = "Mean"+delim+lbname[n];
					for (int t=0;t<nt;t++) line+=(delim+mean[n][t]);
					line+=("\n");
					output.add(line);
				}
			}
			if (statistics.get(s).equals("Stdev")) {
				float[][] mean = new float[nlabels][nt];
				float[][] std = new float[nlabels][nt];
				float[] den = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					for (int t=0;t<nt;t++) mean[n][t] = 0.0f;
					for (int t=0;t<nt;t++) std[n][t] = 0.0f;
					den[n] = 0.0f;
				}
				for (int v=0;v<nv;v++)  {
					for (int n=0;n<nlabels;n++) if (labelmap[v]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += data[v][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int v=0;v<nv;v++) {
					for (int n=0;n<nlabels;n++) if (labelmap[v]==lbid[n]) {
						for (int t=0;t<nt;t++) std[n][t] += (data[v][t]-mean[n][t])*(data[v][t]-mean[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) std[n][t] = (float)FastMath.sqrt(std[n][t]/den[n]);
				}
				for (int n=0;n<nlabels;n++) {
					line = "Stdev"+delim+lbname[n];
					for (int t=0;t<nt;t++) line+=(delim+std[n][t]);
					line+=("\n");
					output.add(line);
				}
			}
		}	
		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam.getURI()+"/"+name+"_statistics.csv");
		String filename = statsParam.getValue().getAbsolutePath()+"/"+name+"_statistics.csv";
		
		// write the output to file
		writeStatisticsFile(output, filename);
		
		// make a copy in output?
		
		outputParam.setValue(new File(filename));
		
		return;
	}

	private final void writeStatisticsFile(ArrayList<String> list, String name) {
		try {
            File f = new File(name);
            FileWriter fw = new FileWriter(f);
            PrintWriter pw = new PrintWriter( fw );
			pw.write("CBSTools Surface Time Series Statistics File\n");
            for (int n=0;n<list.size();n++) {
				pw.write(list.get(n));
				System.out.print(list.get(n));
			}
			pw.close();
            fw.close();
       }
        catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
        catch (IOException e) {
            System.out.println(e.getMessage());
        } 
		catch (OutOfMemoryError e){
			System.out.println(e.getMessage());
		}
		catch (Exception e) {
			System.out.println(e.getMessage());
        }
	}	
}
