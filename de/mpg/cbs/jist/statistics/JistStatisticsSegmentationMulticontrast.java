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
public class JistStatisticsSegmentationMulticontrast extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	segImage;
	private ParamVolume 	intens1Image;
	private ParamVolume 	intens2Image;
	private ParamVolume 	intens3Image;
	private ParamFile 		atlasParam;
	private ParamBoolean	removeFirstParam;
	private ParamBoolean	ignoreZeroParam;
	
	private ParamFile 		statsParam;
	
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private static final String[] statTypes = {"--- single contrast ---", "Mean_intensity", "Std_intensity",
												"Median_intensity","Inter_quantile_range",
												"--- comparisons ---", 
												"Mutual_information",
												"--- joint modeling ---",
												"Region_PCA", "Boundary_PCA"};
	
	private ParamFile 		outputParam;
		
	private String delim = ",";
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(intens1Image = new ParamVolume("Intensity Image 1"));
		inputParams.add(intens2Image = new ParamVolume("Intensity Image 2 (opt)"));
		inputParams.add(intens3Image = new ParamVolume("Intensity Image 3 (opt)"));
		intens2Image.setMandatory(false);
		intens3Image.setMandatory(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		atlasParam.setMandatory(false);
		inputParams.add(removeFirstParam = new ParamBoolean("skip first label",true));
		inputParams.add(ignoreZeroParam = new ParamBoolean("ignore zero intensities",true));
		
		inputParams.add(statsParam = new ParamFile("Spreadsheet file directory", DialogType.DIRECTORY));
		
		inputParams.add(stat1Param = new ParamOption("Statistic 1", statTypes));
		inputParams.add(stat2Param = new ParamOption("Statistic 2", statTypes));
		inputParams.add(stat3Param = new ParamOption("Statistic 3", statTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics");
		inputParams.setLabel("Multi-contrast Statistics");
		inputParams.setName("MultiContrastStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various statistics for intensity images and outputs the results into an organized spreadsheet.");
		
		info.setVersion("3.1.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputParam=new ParamFile());
		
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
		
		int nint=0;
		if (intens1Image.getImageData() != null) nint++;
		if (intens2Image.getImageData() != null) nint++;
		if (intens3Image.getImageData() != null) nint++;
		
		float[][][][] intensity = new float[nint][][][];
		if (intens1Image.getImageData() != null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens1Image.getImageData());
			intensity[0] = intensImg.toArray3d();	 
			intensImg = null;
		}
		if (intens2Image.getImageData() != null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens2Image.getImageData());
			intensity[1] = intensImg.toArray3d();	 
			intensImg = null;
		}
		if (intens3Image.getImageData() != null) {
			ImageDataFloat intensImg = new ImageDataFloat(intens3Image.getImageData());
			intensity[2] = intensImg.toArray3d();	 
			intensImg = null;
		}
		
		// Main algorithm
		
		String line;
		
		String imgtag,notag,reftag;
		String[] inttag = new String[nint];
		imgtag = delim+segImage.getImageData().getName();
		notag =  delim+"-";
		if (intens1Image.getImageData()!=null) inttag[0] = delim+intens1Image.getImageData().getName();
		if (intens2Image.getImageData()!=null) inttag[1] = delim+intens2Image.getImageData().getName();
		if (intens3Image.getImageData()!=null) inttag[2] = delim+intens3Image.getImageData().getName();
			
		// standard computations for the object
		ArrayList<String> output = new ArrayList();

		// remove the zero label, if needed
		if (removeFirstParam.getValue().booleanValue()) {
			int[] tmp = new int[nlabels-1];
			String[] txt = new String[nlabels-1];
			for (int n=1;n<nlabels;n++) {
				tmp[n-1] = lbid[n];
				txt[n-1] = lbname[n];
			}
			nlabels--;
			lbid = tmp;
			lbname = txt;
		}
		
		// create label name list
		String lbline = "labels"+notag+notag+notag+notag;
		for (int n=0;n<nlabels;n++) lbline += (delim+lbname[n]);
		lbline +=("\n");
		System.out.print(lbline);
		
		ArrayList<String> statistics = new ArrayList();
		statistics.add(stat1Param.getValue());
		statistics.add(stat2Param.getValue());
		statistics.add(stat3Param.getValue());
		
		// pre-compute redundant measures ?
		
		// compute the statistics
		for (int s=0; s<statistics.size(); s++) {
			System.out.print("Statistic: "+statistics.get(s)+"\n");
			if (statistics.get(s).equals("Mean_intensity")) {
				for (int i=0;i<nint;i++) {
					float[] mean = new float[nlabels];
					float[] den = new float[nlabels];
					boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
					for (int n=0;n<nlabels;n++) {
						mean[n] = 0.0f;
						den[n] = 0.0f;
					}
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							mean[n] += intensity[i][x][y][z];
							den[n] += 1.0f;
						}
					}
					for (int n=0;n<nlabels;n++) {
						if (den[n]>0) mean[n] = mean[n]/den[n];
					}
					line = "Mean_intensity"+imgtag+inttag[i]+notag+notag;
					for (int n=0;n<nlabels;n++) line+=(delim+mean[n]);
					line+=("\n");
					output.add(line);
				}
			}
			if (statistics.get(s).equals("Std_intensity")) {
				for (int i=0;i<nint;i++) {
					float[] mean = new float[nlabels];
					float[] std = new float[nlabels];
					float[] den = new float[nlabels];
					boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
					for (int n=0;n<nlabels;n++) {
						mean[n] = 0.0f;
						std[n] = 0.0f;
						den[n] = 0.0f;
					}
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							mean[n] += intensity[i][x][y][z];
							den[n] += 1.0f;
						}
					}
					for (int n=0;n<nlabels;n++) {
						if (den[n]>0) mean[n] = mean[n]/den[n];
					}
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							std[n] += (intensity[i][x][y][z]-mean[n])*(intensity[i][x][y][z]-mean[n]);
						}
					}
					for (int n=0;n<nlabels;n++) {
						if (den[n]>0) std[n] = (float)Math.sqrt(std[n]/den[n]);
					}
					line = "Std_intensity"+imgtag+inttag[i]+notag+notag;
					for (int n=0;n<nlabels;n++) line+=(delim+std[n]);
					line+=("\n");
					output.add(line);
				}
			}
			if (statistics.get(s).equals("Median_intensity")) {
				for (int i=0;i<nint;i++) {
					float Imin = intensity[i][0][0][0];
					float Imax = intensity[i][0][0][0];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						if (intensity[i][x][y][z]> Imax) Imax = intensity[i][x][y][z];
						if (intensity[i][x][y][z]< Imin) Imin = intensity[i][x][y][z];
					}
					int Nbins = 100;
					float[][] hist = new float[nlabels][Nbins+1];
					boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							int bin = Numerics.floor((intensity[i][x][y][z]-Imin)/(Imax-Imin)*Nbins);
							if (bin<0) bin = 0;
							if (bin>=Nbins) bin = Nbins-1;
							hist[n][bin]++;
							hist[n][Nbins]++;	// total number
						}
					}
					float[] per = new float[nlabels];
					for (int n=0;n<nlabels;n++) {
						float count = 0.0f;
						int bin = 0;
						while (count<0.5f*hist[n][Nbins]) {
							count += hist[n][bin];
							bin++;
						}
						per[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
					}
					line = "Median_intensity"+imgtag+inttag[i]+notag+notag;
					for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
					line+=("\n");
					output.add(line);
				}
			}
			if (statistics.get(s).equals("Median_intensity")) {
				for (int i=0;i<nint;i++) {
					float Imin = intensity[i][0][0][0];
					float Imax = intensity[i][0][0][0];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						if (intensity[i][x][y][z]> Imax) Imax = intensity[i][x][y][z];
						if (intensity[i][x][y][z]< Imin) Imin = intensity[i][x][y][z];
					}
					int Nbins = 100;
					float[][] hist = new float[nlabels][Nbins+1];
					boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							int bin = Numerics.floor((intensity[i][x][y][z]-Imin)/(Imax-Imin)*Nbins);
							if (bin<0) bin = 0;
							if (bin>=Nbins) bin = Nbins-1;
							hist[n][bin]++;
							hist[n][Nbins]++;	// total number
						}
					}
					float[] per25 = new float[nlabels];
					float[] per75 = new float[nlabels];
					for (int n=0;n<nlabels;n++) {
						float count = 0.0f;
						int bin = 0;
						while (count<0.25f*hist[n][Nbins]) {
							count += hist[n][bin];
							bin++;
						}
						per25[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
						
						count = 0.0f;
						bin = 0;
						while (count<0.75f*hist[n][Nbins]) {
							count += hist[n][bin];
							bin++;
						}
						per75[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
					}
					line = "Inter_quantile_range"+imgtag+inttag[i]+notag+notag;
					for (int n=0;n<nlabels;n++) line+=(delim+(per75[n]-per25[n]));
					line+=("\n");
					output.add(line);
				}
			}
			if (statistics.get(s).equals("Mutual_information")) {
				int Nbins = 100;
				int njoint = 0;
				if (nint==2) njoint = 1;
				else if (nint==3) njoint = 3;
				float[][][] hist = new float[nint][nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				float[] Imin = new float[nint];
				float[] Imax = new float[nint];
				for (int i=0;i<nint;i++) {
					Imin[i] = intensity[i][0][0][0];
					Imax[i] = intensity[i][0][0][0];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						if (intensity[i][x][y][z]> Imax[i]) Imax[i] = intensity[i][x][y][z];
						if (intensity[i][x][y][z]< Imin[i]) Imin[i] = intensity[i][x][y][z];
					}
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[i][x][y][z]!=0) ) {
							int bin = Numerics.floor((intensity[i][x][y][z]-Imin[i])/(Imax[i]-Imin[i])*Nbins);
							if (bin<0) bin = 0;
							if (bin>=Nbins) bin = Nbins-1;
							hist[i][n][bin]++;
							hist[i][n][Nbins]++;	// total number
						}
					}
				}
				float[][][][] jointhist = new float[njoint][nlabels][Nbins+1][Nbins+1];
				for (int i=0;i<njoint;i++) {
					int i1=0,i2=0;
					if (i==0) { i1 = 0; i2 = 1; }
					else if (i==1) {i1 = 1; i2 = 2; }
					else if (i==2) {i1 = 2; i2 = 0; }
					
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || (intensity[i1][x][y][z]!=0 && intensity[i2][x][y][z]!=0) ) ) {
							int bin1 = Numerics.floor((intensity[i1][x][y][z]-Imin[i1])/(Imax[i1]-Imin[i1])*Nbins);
							int bin2 = Numerics.floor((intensity[i2][x][y][z]-Imin[i2])/(Imax[i2]-Imin[i2])*Nbins);
							if (bin1<0) bin1 = 0;
							if (bin1>=Nbins) bin1 = Nbins-1;
							if (bin2<0) bin2 = 0;
							if (bin2>=Nbins) bin2 = Nbins-1;
							jointhist[i][n][bin1][bin2]++;
							jointhist[i][n][Nbins][Nbins]++;	// total number
						}
					}
					double[] mi = new double[nlabels];
					for (int n=0;n<nlabels;n++) {
						for (int a=0;a<Nbins;a++) for (int b=0;b<Nbins;b++) if (jointhist[i][n][a][b]>0) {
							mi[n] += jointhist[i][n][a][b]/jointhist[i][n][Nbins][Nbins]*FastMath.log(jointhist[i][n][a][b]/jointhist[i][n][Nbins][Nbins]
																										/hist[i1][n][a]*hist[i1][n][Nbins]/hist[i2][n][b]*hist[i2][n][Nbins]);
						}
					}
					line = "Mutual_information"+imgtag+inttag[i1]+inttag[i2]+notag;
					for (int n=0;n<nlabels;n++) line+=(delim+mi[n]);
					line+=("\n");
					output.add(line);
				}
			}
		}	
		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam.getURI()+"/volume_statistics.csv");
		String filename = statsParam.getValue().getAbsolutePath()+"/volume_statistics.csv";
		
		// write the output to file
		addStatisticsToFile(filename, output, lbline);
		
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
			String type = added.get(n).substring(0,added.get(n).indexOf(delim));
			System.out.println(added.get(n));
			System.out.println(type);
			
			// find the last line with this type
			int last=-1;
			for (int m=0;m<main.size();m++) {
				if (main.get(m).indexOf(delim)>-1)
					if (main.get(m).substring(0,main.get(m).indexOf(delim)).equals(type)) last = m;
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
            if (!line.startsWith("MIPAV Volumetric Statistics File")) {
                System.out.println("not a proper MIPAV statistics file");
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
			pw.write("MIPAV Volumetric Statistics File\n");
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
