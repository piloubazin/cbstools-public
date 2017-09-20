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
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
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


/*
 * @author Pierre-Louis Bazin
 */
public class JistStatisticsLesions extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	lesionImage;
	private ParamVolume 	ventricleImage;
	private ParamVolume 	cortexImage;
	private ParamVolume 	intensImage;
	
	private ParamFile 		statsParam;
	
	/*
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private static final String[] statTypes = {"Lesion_volume", "Lesion_wm_ratio", "Periventricular_ratio", "Pericortical_ratio", 
												"Lesion_count", "Mean_size", "Max_size",
												"Mean_intensity", "Std_intensity", "none"};
	*/
	private ParamFile 		outputParam;
	private ParamVolume		lesionlbImage;
	private ParamVolume		locationImage;
		
	private String delim = ",";
		
	private static final byte GENERAL 			= 1;
	private static final byte PERIVENTRICULAR 	= 2;
	private static final byte DEEPWM 			= 3;
	private static final byte PERICORTICAL 	= 4;
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(lesionImage = new ParamVolume("Lesion mask"));
		inputParams.add(ventricleImage = new ParamVolume("Ventricle mask (opt)"));
		inputParams.add(cortexImage = new ParamVolume("Cortex mask (opt)"));
		inputParams.add(intensImage = new ParamVolume("Intensity Image (opt)"));
		ventricleImage.setMandatory(false);
		cortexImage.setMandatory(false);
		intensImage.setMandatory(false);
		
		inputParams.add(statsParam = new ParamFile("Spreadsheet file directory", DialogType.DIRECTORY));
		/*
		inputParams.add(stat1Param = new ParamOption("Statistic 1", statTypes));
		inputParams.add(stat2Param = new ParamOption("Statistic 2", statTypes));
		inputParams.add(stat3Param = new ParamOption("Statistic 3", statTypes));
		*/
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics.devel");
		inputParams.setLabel("Lesions Statistics");
		inputParams.setName("LesionsStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute sub-maps and statistics of detected lesions and outputs the results into an organized spreadsheet.");
		
		info.setVersion("3.0.1");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputParam=new ParamFile("Statistics File"));
		outputParams.add(lesionlbImage=new ParamVolume("Individual Lesions",VoxelType.INT));
		outputParams.add(locationImage=new ParamVolume("Lesions location",VoxelType.UBYTE));

		outputParams.setName("spreadsheet result");
		outputParams.setLabel("spreadsheet result");
	}

	@Override
	protected void execute(CalculationMonitor monitor) {
		
		// import the segmentation data into 1D arrays
		ImageDataUByte lesionImg = new ImageDataUByte(lesionImage.getImageData());
		int nx = lesionImg.getRows();
		int ny = lesionImg.getCols();
		int nz = lesionImg.getSlices();
		int nc = lesionImg.getComponents();
		int nxyz = nx*ny*nz;
		BasicInfo.displayMessage("Image dims: "+nx+", "+ny+", "+nz+"; "+nc+"\n");
		float rx = lesionImg.getHeader().getDimResolutions()[0];
		float ry = lesionImg.getHeader().getDimResolutions()[1];
		float rz = lesionImg.getHeader().getDimResolutions()[2];
		
		int orient = lesionImg.getHeader().getImageOrientation().ordinal();
		int orx = lesionImg.getHeader().getAxisOrientation()[0].ordinal();
		int ory = lesionImg.getHeader().getAxisOrientation()[1].ordinal();
		int orz = lesionImg.getHeader().getAxisOrientation()[2].ordinal();
		
		// get label list from atlas or segmentation
		byte[][][] lesions = lesionImg.toArray3d();

		byte[][][] ventricles = null;
		if (ventricleImage.getImageData() != null) {
			ImageDataUByte	ventImg = new ImageDataUByte(ventricleImage.getImageData());
			ventricles = ventImg.toArray3d();	 
			ventImg = null;
		}
		
		byte[][][] cortex = null;
		if (cortexImage.getImageData() != null) {
			ImageDataUByte	ctxImg = new ImageDataUByte(cortexImage.getImageData());
			cortex = ctxImg.toArray3d();	 
			ctxImg = null;
		}
		
		float[][][] intensity = null;
		if (intensImage.getImageData() != null) {
			ImageDataFloat intensImg = new ImageDataFloat(intensImage.getImageData());
			intensity = intensImg.toArray3d();	 
			intensImg = null;
		}
		
		// Main algorithm
		
		// 1. Find separate lesions
		boolean[][][] obj = new boolean[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (lesions[x][y][z]>0) obj[x][y][z] = true;
			else obj[x][y][z] = false;
		}
		int[][][] lesionlb = ObjectLabeling.connected26Object3D(obj,nx,ny,nz);
		int[] lblist = ObjectLabeling.listOrderedLabels(lesionlb, nx,ny,nz);
		int nlb = lblist.length;
		
		// 2. find peri-ventricular and peri-cortical lesions
		byte[][][] location = new byte[nx][ny][nz];
		int[] vcount = new int[nlb];
		int[] ccount = new int[nlb];
		int[] gcount = new int[nlb];
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) if (obj[x][y][z]) {
			if (ventricles==null && cortex==null) {
				location[x][y][z] = GENERAL;
			} else {
				for (int n=0;n<nlb;n++) if (lblist[n]==lesionlb[x][y][z]) {
					gcount[n]++;
					if (ventricles!=null && maxNeighbor26(ventricles,x,y,z)>0) vcount[n]++;
					if (cortex!=null && maxNeighbor26(cortex,x,y,z)>0) ccount[n]++;
				}
			}
		}
		byte[] type = new byte[nlb];
		for (int n=1;n<nlb;n++) {
			if (ventricles==null && cortex==null) {
				type[n] = GENERAL;
			} else {
				if (vcount[n]>ccount[n]) type[n] = PERIVENTRICULAR;
				else if (ccount[n]>vcount[n]) type[n] = PERICORTICAL;
				else type[n] = DEEPWM;
			}
		}
			
		if (ventricles!=null || cortex!=null) {
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (obj[x][y][z]) {
				for (int n=0;n<nlb;n++) if (lblist[n]==lesionlb[x][y][z]) {
					location[x][y][z] = type[n];
				}
			}
		}
		
		// Statistics
		
		String line;
		
		String imgtag,notag,inttag;
		imgtag = delim+lesionImage.getImageData().getName();
		notag =  delim+"-";
		if (intensImage.getImageData()!=null) inttag = delim+intensImage.getImageData().getName();
		else inttag = delim+"-";
			
		// standard computations for the object
		ArrayList<String> output = new ArrayList();

		// create label name list
		String lbline = "labels"+notag+notag;
		lbline += (delim+"general");
		lbline += (delim+"periventricular");
		lbline += (delim+"pericortical");
		lbline += (delim+"deep_wm");
		lbline +=("\n");
		System.out.print(lbline);
		
		int nstat = 4;
		byte[] statlb = new byte[nstat];
		statlb[0] = GENERAL;
		statlb[1] = PERIVENTRICULAR;
		statlb[2] = PERICORTICAL;
		statlb[3] = DEEPWM;
		
		// compute the statistics : 	"Lesion_volume",
		//								"Lesion_count", "Mean_size", "Max_size",
		//								"Mean_intensity", "Std_intensity"

		float[] volume = new float[nstat];
		for (int n=0;n<nstat;n++) {
			obj = ObjectExtraction.objectFromLabelImage(location,nx,ny,nz,statlb[n],ObjectLabeling.EQUAL);
			volume[n] = ObjectStatistics.volume(obj,nx,ny,nz);
		}
		line = "Lesion_volumes"+imgtag+notag;
		for (int n=0;n<nstat;n++) line += (delim+volume[n]*rx*ry*rz);
		line +="\n";
		output.add(line);

		float[] meansize = new float[nstat];
		float[] maxsize = new float[nstat];
		float[] count = new float[nstat];
		for (int l=1;l<nlb;l++) {
			for (int n=0;n<nstat;n++) if (type[l]==statlb[n]) {
				obj = ObjectExtraction.objectFromLabelImage(lesionlb,nx,ny,nz,lblist[l],ObjectLabeling.EQUAL);
				float size = ObjectStatistics.volume(obj,nx,ny,nz);
				meansize[n] += size;
				count[n]++;
				if (size>maxsize[n]) maxsize[n] = size;
			}
		}
		for (int n=0;n<nstat;n++) if (count[n]>0) meansize[n] /= count[n];
				
		line = "Lesion_count"+imgtag+notag;
		for (int n=0;n<nstat;n++) line += (delim+count[n]);
		line +="\n";
		output.add(line);

		line = "Mean_size"+imgtag+notag;
		for (int n=0;n<nstat;n++) line += (delim+meansize[n]*rx*ry*rz);
		line +="\n";
		output.add(line);

		line = "Max_size"+imgtag+notag;
		for (int n=0;n<nstat;n++) line += (delim+maxsize[n]*rx*ry*rz);
		line +="\n";
		output.add(line);

		if (intensity!=null) {
			float[] mean = new float[nstat];
			float[] den = new float[nstat];
			float[] std = new float[nstat];
			for (int n=0;n<nstat;n++) {
				mean[n] = 0.0f;
				std[n] = 0.0f;
					den[n] = 0.0f;
			}
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				for (int n=0;n<nstat;n++) if (location[x][y][z]==lblist[n]) {
					mean[n] += intensity[x][y][z];
					den[n] += 1.0f;
				}
			}
			for (int n=0;n<nstat;n++) {
				if (den[n]>0) mean[n] = mean[n]/den[n];
			}
			line = "Mean_intensity"+imgtag+inttag;
			for (int n=0;n<nstat;n++) line+=(delim+mean[n]);
			line+=("\n");
			output.add(line);

			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				for (int n=0;n<nstat;n++) if (location[x][y][z]==lblist[n]) {
					std[n] += (intensity[x][y][z]-mean[n])*(intensity[x][y][z]-mean[n]);
				}
			}
			for (int n=0;n<nstat;n++) {
				if (den[n]>0) std[n] = (float)Math.sqrt(std[n]/den[n]);
			}
			line = "Std_intensity"+imgtag+inttag;
			for (int n=0;n<nstat;n++) line+=(delim+std[n]);
			line+=("\n");
			output.add(line);
		}

		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam.getURI()+"/lesions_statistics.csv");
		String filename = statsParam.getValue().getAbsolutePath()+"/lesions_statistics.csv";
		
		// write the output to file
		addStatisticsToFile(filename, output, lbline);
		
		// make a copy in output
		outputParam.setValue(new File(filename));
		
		// write the images
		ImageDataInt bufferData = new ImageDataInt(lesionlb);
		bufferData.setHeader(lesionImage.getImageData().getHeader());
		bufferData.setName(lesionImage.getImageData().getName()+"_lb");
		lesionlbImage.setValue(bufferData);
		bufferData = null;
		
		ImageDataUByte byteData = new ImageDataUByte(location);
		byteData.setHeader(lesionImage.getImageData().getHeader());
		byteData.setName(lesionImage.getImageData().getName()+"_loc");
		locationImage.setValue(byteData);
		byteData = null;
		
		return;
	}
	
	private byte maxNeighbor26(byte[][][] img, int x, int y, int z) {
		byte max = 0;
		for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int k=-1;k<=1;k++) if (i*i+j*j+k*k>0) {
			if (img[x+i][y+j][z+k]>max) max = img[x+i][y+j][z+k];
		}
		return max;
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
