package de.mpg.cbs.jist.statistics;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
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
public class JistStatisticsTimeSeries extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	segImage;
	private ParamVolume 	intensImage;
	private ParamFile 		atlasParam;
	private ParamBoolean	removeFirstParam;
	
	private ParamFile 		statsParam;
	
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private ParamOption 	stat4Param;
	private static final String[] statTypes = {"Mean", "Stdev", "Skewness", "Kurtosis", "none"};
	
	private ParamFile 		outputParam;
		
	private String delim = ",";
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensImage = new ParamVolume("Time series Image"));
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		
		inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		atlasParam.setMandatory(false);
		inputParams.add(removeFirstParam = new ParamBoolean("skip first label",true));
		
		inputParams.add(statsParam = new ParamFile("Spreadsheet file directory", DialogType.DIRECTORY));
		
		inputParams.add(stat1Param = new ParamOption("Statistic 1", statTypes));
		inputParams.add(stat2Param = new ParamOption("Statistic 2", statTypes));
		inputParams.add(stat3Param = new ParamOption("Statistic 3", statTypes));
		inputParams.add(stat4Param = new ParamOption("Statistic 4", statTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics");
		inputParams.setLabel("Time Series Statistics");
		inputParams.setName("TimeSeriesStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various statistics from time series and volumetric labels and outputs to a spreadsheet.");
		
		info.setVersion("3.0");
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
		
		// import the segmentation data into 1D arrays
		ImageDataInt segImg = new ImageDataInt(segImage.getImageData());
		int nx = segImg.getRows();
		int ny = segImg.getCols();
		int nz = segImg.getSlices();
		int nc = segImg.getComponents();
		int nxyz = nx*ny*nz;
		BasicInfo.displayMessage("Image dims: "+nx+", "+ny+", "+nz+"; "+nc+"\n");
		float rx = segImg.getHeader().getDimResolutions()[0];
		float ry = segImg.getHeader().getDimResolutions()[1];
		float rz = segImg.getHeader().getDimResolutions()[2];
		
		int orient = segImg.getHeader().getImageOrientation().ordinal();
		int orx = segImg.getHeader().getAxisOrientation()[0].ordinal();
		int ory = segImg.getHeader().getAxisOrientation()[1].ordinal();
		int orz = segImg.getHeader().getAxisOrientation()[2].ordinal();
		
		String segname = segImg.getName();
		
		// get label list from atlas or segmentation
		
		int nlabels;
		String[] lbname;
		int[] lbid;
		int[][][] segmentation;

		segmentation = segImg.toArray3d();
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
			lbid = ObjectLabeling.listOrderedLabels(segmentation, nx, ny, nz);
			nlabels = lbid.length;
			lbname = new String[nlabels];
			for (int n=0;n<nlabels;n++) lbname[n] = new String("Label_"+lbid[n]);
		}
		
		float[][][][] intensity = null;
		ImageDataFloat intensImg = new ImageDataFloat(intensImage.getImageData());
		intensity = intensImg.toArray4d();	 
		int nt = intensImg.getComponents();
		String name = intensImg.getName();
		intensImg = null;
		
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
		statistics.add(stat4Param.getValue());
		
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
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += intensity[x][y][z][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int n=0;n<nlabels;n++) {
					line = "Mean"+delim+segname+delim+lbname[n];
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
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += intensity[x][y][z][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) std[n][t] += (intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>1) for (int t=0;t<nt;t++) std[n][t] = (float)FastMath.sqrt(std[n][t]/(den[n]-1));
				}
				for (int n=0;n<nlabels;n++) {
					line = "Stdev"+delim+segname+delim+lbname[n];
					for (int t=0;t<nt;t++) line+=(delim+std[n][t]);
					line+=("\n");
					output.add(line);
				}
			} else if (statistics.get(s).equals("Skewness")) {
				float[][] mean = new float[nlabels][nt];
				float[][] std = new float[nlabels][nt];
				float[][] skew = new float[nlabels][nt];
				float[] den = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					for (int t=0;t<nt;t++) mean[n][t] = 0.0f;
					for (int t=0;t<nt;t++) std[n][t] = 0.0f;
					for (int t=0;t<nt;t++) skew[n][t] = 0.0f;
					den[n] = 0.0f;
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += intensity[x][y][z][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) std[n][t] += (intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) std[n][t] = (float)FastMath.sqrt(std[n][t]/den[n]);
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) skew[n][t] += (intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t])/(std[n][t]*std[n][t]*std[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) skew[n][t] /= den[n];
				}
				for (int n=0;n<nlabels;n++) {
					line = "Skewness"+delim+segname+delim+lbname[n];
					for (int t=0;t<nt;t++) line+=(delim+skew[n][t]);
					line+=("\n");
					output.add(line);
				}
			} else if (statistics.get(s).equals("Kurtosis")) {
				float[][] mean = new float[nlabels][nt];
				float[][] mo2 = new float[nlabels][nt];
				float[][] kurt = new float[nlabels][nt];
				float[] den = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					for (int t=0;t<nt;t++) mean[n][t] = 0.0f;
					for (int t=0;t<nt;t++) mo2[n][t] = 0.0f;
					for (int t=0;t<nt;t++) kurt[n][t] = 0.0f;
					den[n] = 0.0f;
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) mean[n][t] += intensity[x][y][z][t];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mean[n][t] /= den[n];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) mo2[n][t] += (intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) mo2[n][t] = mo2[n][t]/den[n];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n]) {
						for (int t=0;t<nt;t++) kurt[n][t] += (intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t])
														*(intensity[x][y][z][t]-mean[n][t])*(intensity[x][y][z][t]-mean[n][t])
														/(mo2[n][t]*mo2[n][t]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) for (int t=0;t<nt;t++) kurt[n][t] = kurt[n][t]/den[n]-3.0f;
				}
				for (int n=0;n<nlabels;n++) {
					line = "Kurtosis"+delim+segname+delim+lbname[n];
					for (int t=0;t<nt;t++) line+=(delim+kurt[n][t]);
					line+=("\n");
					output.add(line);
				}
			}		}	
		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam.getURI()+"/"+name+"_statistics.csv");
		String filename = statsParam.getValue().getAbsolutePath()+"/"+name+"_statistics.csv";
		
		// write the output to file
		//writeStatisticsFile(output, filename);
		addStatisticsToFile(filename, output, "");
			
		// make a copy in output?
		
		outputParam.setValue(new File(filename));
		
		return;
	}

	private final void addStatisticsToFile(String name, ArrayList<String> output, String lbline) {
		
		// open the file
		ArrayList<String> 	previous = loadStatisticsFile(name);
		
		// merge the output
		appendStatistics(previous, output, lbline);
		
		// save the result
		writeStatisticsFile(previous, name);
	}
	
	private	final void appendStatistics(ArrayList<String> main, ArrayList<String> added, String lbline) {
		for (int n=0;n<added.size();n++) {
			// extract statistics type
			String type = added.get(n).substring(0,added.get(n).indexOf(delim,added.get(n).indexOf(delim)+1));
			System.out.println(added.get(n));
			System.out.println(type);
			
			// find the last line with this type
			int last=-1;
			for (int m=0;m<main.size();m++) {
				if (main.get(m).indexOf(delim)>-1)
					if (main.get(m).substring(0,main.get(m).indexOf(delim,main.get(m).indexOf(delim)+1)).equals(type)) last = m;
			}
			if (last>-1) {
				main.add(last+1, added.get(n));
			} else {
				// add a space after each different statistic as well as labels before
				main.add(lbline);
				main.add(added.get(n));
				main.add(" \n");
			}
		}
	}
	
	private final ArrayList<String> loadStatisticsFile(String name) {
		ArrayList<String> list = new ArrayList();
		try {
            System.out.println("reading previous statistic file: "+name);
            File f = new File(name);
            FileReader fr = new FileReader(f);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
			
			// Exact corresponding template for first line ?
            if (!line.startsWith("CBSTools Time Series Statistics File")) {
                System.out.println("not a proper statistics file");
                br.close();
                fr.close();
                return null;
            }
			line = br.readLine();
			while (line!=null) {
				list.add(line+"\n");
				line = br.readLine();
				System.out.println(line);
			}
			br.close();
            fr.close();
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
		return list;	
	}

	private final void writeStatisticsFile(ArrayList<String> list, String name) {
		try {
            File f = new File(name);
            FileWriter fw = new FileWriter(f);
            PrintWriter pw = new PrintWriter( fw );
			pw.write("CBSTools Time Series Statistics File\n");
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
