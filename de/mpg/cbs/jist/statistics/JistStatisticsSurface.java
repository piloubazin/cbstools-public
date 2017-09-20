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

import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
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


/*
 * @author Pierre-Louis Bazin
 */
public class JistStatisticsSurface extends ProcessingAlgorithm {

	// jist containers
	private ParamSurface 	segSurface;
	private ParamSurface 	refSurface;
	private ParamFile 		atlasParam;
	private ParamBoolean	removeFirstParam;
	
	private ParamFile 		statsParam;
	
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private static final String[] statTypes = {"Dice_overlap", "Jaccard_overlap", "Surface_area", "Surface_ratio", "none"};
	
	private ParamFile 		outputParam;
		
	private String delim = ",";
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segSurface=new ParamSurface("Segmented Surface"));
		inputParams.add(refSurface=new ParamSurface("Reference Surface (opt)"));
		segSurface.setLoadAndSaveOnValidate(false);
		refSurface.setLoadAndSaveOnValidate(false);
		refSurface.setMandatory(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		atlasParam.setMandatory(false);
		inputParams.add(removeFirstParam = new ParamBoolean("skip first label",true));
		
		inputParams.add(statsParam = new ParamFile("Spreadsheet file directory", DialogType.DIRECTORY));
		
		inputParams.add(stat1Param = new ParamOption("Statistic 1", statTypes));
		inputParams.add(stat2Param = new ParamOption("Statistic 2", statTypes));
		inputParams.add(stat3Param = new ParamOption("Statistic 3", statTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics");
		inputParams.setLabel("Segmentation Statistics (surface)");
		inputParams.setName("SurfaceSegmentationStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various statistics from surface labels and outputs to a spreadsheet (assumes both surfaces have same geometry or at least matched vertices).");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.BETA);
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
		EmbeddedSurface seg = segSurface.getSurface();
		double[][] segData = seg.getVertexData();
		
		EmbeddedSurface ref = null;
		double[][] refData = null;
		if (refSurface.getSurface()!=null) {
			ref = refSurface.getSurface();
			refData = ref.getVertexData();
		}
		
		String imgtag,notag,reftag,inttag;
		imgtag = delim+seg.getName();
		reftag =  delim+"-";
		notag =  delim+"-";
		inttag = delim+"-";
		if (ref!=null) reftag = delim+ref.getName();
		
		// turn into label lists
		int nlist = segData.length;
		BasicInfo.displayMessage("Mesh size: "+nlist+"\n");
		
		int[] segList = new int[nlist];
		for (int n=0;n<nlist;n++) segList[n] = (int)segData[n][0];
		
		int[] refList = null;
		if (ref!=null) {
			refList = new int[nlist];
			for (int n=0;n<nlist;n++) refList[n] = (int)refData[n][0];
		}
		segData = null; refData = null;
		segSurface.dispose();
		refSurface.dispose();
		
		// get label list from atlas or segmentation
		
		int nlabels;
		String[] lbname;
		int[] lbid;
		
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
			// find the number of labels from the segmented image
			lbid = ObjectLabeling.listOrderedLabels(segList, nlist, 1, 1);
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
		String lbline = "labels"+notag+notag+notag;
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
			if (statistics.get(s).equals("Dice_overlap")) {
				// compare: requires same labels (use the reference as basis)
				float[] dice = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segList,nlist,1,1,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nlist,1,1);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(refList,nlist,1,1,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nlist,1,1);
					boolean[] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.XOR,nlist,1,1);
					float diff = ObjectStatistics.volume(xor,nlist,1,1);
									
					dice[n] = ((vol1+vol2-diff)/(vol1+vol2));
				}
				line = "Dice_overlap"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+dice[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Jaccard_overlap")) {
				// compare: requires same labels (use the template as basis)
				float[] jaccard = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segList,nlist,1,1,lbid[n],ObjectLabeling.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nlist,1,1);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(refList,nlist,1,1,lbid[n],ObjectLabeling.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nlist,1,1);
					boolean[] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectLabeling.XOR,nlist,1,1);
					float diff = ObjectStatistics.volume(xor,nlist,1,1);
									
					jaccard[n] = ((vol1+vol2-diff)/(vol1+vol2+diff));
				}
				line = "Jaccard_overlap"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+jaccard[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Surface_area")) {
				// compare: requires same labels (use the reference as basis)
				double[] surface = new double[nlabels+1];
				for (int f=0;f<seg.getFaceCount();f++) {
					float area = seg.getFaceArea(f);
					
					int p1 = seg.getCoordinateIndex(3*f + 0);
					int p2 = seg.getCoordinateIndex(3*f + 1);
					int p3 = seg.getCoordinateIndex(3*f + 2);
		
					for (int n=0;n<nlabels;n++) {
						if (segList[p1]==lbid[n]) {
							surface[n] += area/3.0f;
						}
						if (segList[p2]==lbid[n]) {
							surface[n] += area/3.0f;
						}
						if (segList[p3]==lbid[n]) {
							surface[n] += area/3.0f;
						}
					}
					/* does not include borders
					if (segList[p1]==segList[p2] && segList[p2]==segList[p3]) {
						for (int n=0;n<nlabels;n++) {
							if (segList[p1]==lbid[n]) {
								surface[n] += area;
							}
						}
					}
					*/
				}
				line = "Surface_area"+imgtag+notag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+surface[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Surface_ratio")) {
				// compare: requires same labels (use the reference as basis)
				double[] surface = new double[nlabels+1];
				for (int f=0;f<seg.getFaceCount();f++) {
					float area = seg.getFaceArea(f);
					
					int p1 = seg.getCoordinateIndex(3*f + 0);
					int p2 = seg.getCoordinateIndex(3*f + 1);
					int p3 = seg.getCoordinateIndex(3*f + 2);
		
					for (int n=0;n<nlabels;n++) {
						if (segList[p1]==lbid[n]) {
							surface[n] += area/3.0f;
						}
						if (segList[p2]==lbid[n]) {
							surface[n] += area/3.0f;
						}
						if (segList[p3]==lbid[n]) {
							surface[n] += area/3.0f;
						}
					}
					/* does not include borders
					if (segList[p1]==segList[p2] && segList[p2]==segList[p3]) {
						for (int n=0;n<nlabels;n++) {
							if (segList[p1]==lbid[n]) {
								surface[n] += area;
							}
						}
					}
					*/
					// global area measurement (includes all triangles, even at the boundaries)
					surface[nlabels] += area;
				}
				line = "Surface_ratio"+imgtag+notag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+surface[n]/surface[nlabels]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
		}	
		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam.getURI()+"/surface_statistics.csv");
		String filename = statsParam.getValue().getAbsolutePath()+"/surface_statistics.csv";
		
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
