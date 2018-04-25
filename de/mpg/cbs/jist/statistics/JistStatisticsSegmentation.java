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
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataByte;
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
public class JistStatisticsSegmentation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	segImage;
	private ParamVolume 	intensImage;
	private ParamVolume 	tempImage;
	private ParamFile 		atlasParam;
	private ParamBoolean	removeFirstParam;
	private ParamBoolean	ignoreZeroParam;
	
	private ParamFile 		statsParam;
	
	private ParamOption 	stat1Param;
	private ParamOption 	stat2Param;
	private ParamOption 	stat3Param;
	private static final String[] statTypes = {"none", "Voxels", "Volume", "--- intensity ---", "Mean_intensity", "Std_intensity",
												"10_intensity","25_intensity","50_intensity","75_intensity","90_intensity",
												"--- overlap ---", 
												"Volumes", "Dice_overlap", "Jaccard_overlap", "Volume_difference",
												"False_positives","False_negatives",
												"Dilated_Dice_overlap","Dilated_false_positive","Dilated_false_negative",
												"Dilated_false_negative_volume","Dilated_false_positive_volume",
												"--- Clusters ---",
												"Detected_clusters", "False_detections",
												"Cluster_numbers", "Mean_cluster_sizes", "Cluster_maps",
												"--- boundaries ---",
												"Average_surface_distance", "Average_surface_difference", 
												"Average_squared_surface_distance", "Hausdorff_distance"};
	
	private ParamFile 		outputParam;
    private ParamVolume 	outputImage;
		
	private String delim = ",";
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(intensImage = new ParamVolume("Intensity Image (opt)"));
		inputParams.add(tempImage = new ParamVolume("Template Image (opt)"));
		intensImage.setMandatory(false);
		tempImage.setMandatory(false);
		
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
		inputParams.setLabel("Segmentation Statistics");
		inputParams.setName("SegmentationStatistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute various statistics (geometric, intensity or compared to template labeling) for volumetric segmentation labels and outputs the results into an organized spreadsheet.");
		
		info.setVersion("3.0.7");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputParam=new ParamFile());
		outputParams.add(outputImage=new ParamVolume("Ouput Image",null,-1,-1,-1,-1));
		outputImage.setMandatory(false);
		
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
		
		String name = Interface.getName(segImage);
		ImageHeader header = Interface.getHeader(segImage);
		
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
		
		float[][][] intensity = null;
		if (intensImage.getImageData() != null) {
			ImageDataFloat intensImg = new ImageDataFloat(intensImage.getImageData());
			intensity = intensImg.toArray3d();	 
			intensImg = null;
		}
		
		int[][][] template = null;
		if (tempImage.getImageData() != null) {
			ImageDataInt	tempImg = new ImageDataInt(tempImage.getImageData());
			template = tempImg.toArray3d();	 
			tempImg = null;
		}
		
		// Main algorithm
		
		String line;
		
		String imgtag,notag,reftag,inttag;
		imgtag = delim+segImage.getImageData().getName();
		notag =  delim+"-";
		if (intensImage.getImageData()!=null) inttag = delim+intensImage.getImageData().getName();
		else inttag = delim+"-";
		if (tempImage.getImageData() != null) reftag = delim+tempImage.getImageData().getName();
		else reftag = delim+"-";
			
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
			if (statistics.get(s).equals("Voxels")) {
				// compute the volumes
				float[] volume = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					volume[n] = ObjectStatistics.volume(obj,nx,ny,nz);
				}
				line = "Voxels"+imgtag+notag+notag;
				for (int n=0;n<nlabels;n++) line += (delim+(int)volume[n]);
				line +="\n";
				output.add(line);
			}
			if (statistics.get(s).equals("Volume")) {
				// compute the volumes
				float[] volume = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					volume[n] = ObjectStatistics.volume(obj,nx,ny,nz);
				}
				line = "Volume"+imgtag+notag+notag;
				for (int n=0;n<nlabels;n++) line += (delim+volume[n]*rx*ry*rz);
				line +="\n";
				output.add(line);
			}
			if (statistics.get(s).equals("Volumes")) {
				// compute the volumes
				float[][] volumes = new float[2][nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					volumes[0][n] = ObjectStatistics.volume(obj1,nx,ny,nz);
					volumes[1][n] = ObjectStatistics.volume(obj2,nx,ny,nz);
				}
				line = "SegVolume"+imgtag+notag+notag;
				for (int n=0;n<nlabels;n++) line += (delim+volumes[0][n]*rx*ry*rz);
				line +="\n";
				output.add(line);
				System.out.println(line);
				line = "RefVolume"+notag+reftag+notag;
				for (int n=0;n<nlabels;n++) line += (delim+volumes[1][n]*rx*ry*rz);
				line +="\n";
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dice_overlap")) {
				// compare: requires same labels (use the reference as basis)
				float[] dice = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.XOR,nx,ny,nz);
					float diff = ObjectStatistics.volume(xor,nx,ny,nz);
									
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
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectLabeling.SUPEQUAL,ObjectLabeling.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectLabeling.XOR,nx,ny,nz);
					float diff = ObjectStatistics.volume(xor,nx,ny,nz);
									
					jaccard[n] = ((vol1+vol2-diff)/(vol1+vol2+diff));
				}
				line = "Jaccard_overlap"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+jaccard[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Volume_difference")) {
				// compare: requires same labels (use the reference as basis)
				float[] voldiff = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
									
					voldiff[n] = Numerics.abs(vol1/vol2-1.0f);
				}
				line = "Volume_difference"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+voldiff[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("False_positives")) {
				// compare: requires same labels (use the reference as basis)
				float[] fpos = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fpos[n] = (vol1-diff)/vol1;
				}
				line = "False_positives"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fpos[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("False_negatives")) {
				// compare: requires same labels (use the reference as basis)
				float[] fneg = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fneg[n] = (vol2-diff)/vol2;
				}
				line = "False_negatives"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fneg[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dilated_Dice_overlap")) {
				// compare: requires same labels (use the reference as basis)
				float[] dice = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] objd1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					boolean[][][] and1 = ObjectGeometry.binaryOperation(objd1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[][][] objd2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					boolean[][][] and2 = ObjectGeometry.binaryOperation(obj1,objd2,ObjectExtraction.AND,nx,ny,nz);
					float inter = ObjectStatistics.volume(and1,nx,ny,nz) + ObjectStatistics.volume(and2,nx,ny,nz);
									
					dice[n] = (inter/(vol1+vol2));
				}
				line = "Dilated_Dice_overlap"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+dice[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dilated_false_positive")) {
				// compare: requires same labels (use the reference as basis)
				float[] fpos = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					obj2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fpos[n] = (vol1-diff)/vol1;
				}
				line = "Dilated_false_positive"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fpos[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dilated_false_negative")) {
				// compare: requires same labels (use the reference as basis)
				float[] fneg = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					obj1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fneg[n] = (vol2-diff)/vol2;
				}
				line = "Dilated_false_negative"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fneg[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dilated_false_positive_volume")) {
				// compare: requires same labels (use the reference as basis)
				float[] fpos = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					obj2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fpos[n] = (vol1-diff)*rx*ry*rz;
				}
				line = "Dilated_false_positive_volume"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fpos[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Dilated_false_negative_volume")) {
				// compare: requires same labels (use the reference as basis)
				float[] fneg = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					obj1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					float diff = ObjectStatistics.volume(and,nx,ny,nz);
									
					fneg[n] = (vol2-diff)*rx*ry*rz;
				}
				line = "Dilated_false_negative_volume"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fneg[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Detected_clusters")) {
				// compare: requires same labels (use the reference as basis)
				float[] detc = new float[nlabels];
				float[] nclu = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					// not needed here... int[][][] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label2 = ObjectLabeling.connected18Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[] detected = new boolean[list2.length];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (and[x][y][z]) {
					    for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) detected[l] = true;    	
					}
					// count the final numbers
					detc[n] = 0;
					for (int l=1;l<list2.length;l++) if (detected[l]) detc[n]++;
					nclu[n] = list2.length-1;
				}
				line = "Detected_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+detc[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Total_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+nclu[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("False_detections")) {
				// compare: requires same labels (use the reference as basis)
				float[] fdet = new float[nlabels];
				float[] ndet = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[] detected = new boolean[list1.length];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) if (and[x][y][z]) {
					    for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l]) detected[l] = true;    	
					}
					// count the final numbers
					fdet[n] = 0;
					for (int l=1;l<list1.length;l++) if (detected[l]) fdet[n]++;
					ndet[n] = list1.length-1;
					fdet[n] = ndet[n]-fdet[n];
				}
				line = "False_detections"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fdet[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Total_detections"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+ndet[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Cluster_numbers")) {
				// compare: requires same labels (use the reference as basis)
				float[] trupos = new float[nlabels];
				float[] falpos = new float[nlabels];
				float[] falneg = new float[nlabels];
				float[] totseg = new float[nlabels];
				float[] tottpl = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					    if (and[x][y][z]) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) true_positive.set(l, true);    
                            for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l]) false_positive.set(l, false);    
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) false_negative.set(l, false);   
                       }
                    }
					// count the final numbers
					trupos[n] = true_positive.cardinality();
					falpos[n] = false_positive.cardinality();
					falneg[n] = false_negative.cardinality();
					totseg[n] = list1.length-1;
					tottpl[n] = list2.length-1;
				}
				line = "True_positive_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+trupos[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "False_positive_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+falpos[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "False_negative_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+falneg[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Segmented_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+totseg[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Template_clusters"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+tottpl[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Mean_cluster_sizes")) {
				// compare: requires same labels (use the reference as basis)
				float[] tp_size = new float[nlabels];
				float[] fp_size = new float[nlabels];
				float[] fn_size = new float[nlabels];
				float[] seg_size = new float[nlabels];
				float[] tpl_size = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					// not needed here... int[][][] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					    if (and[x][y][z]) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) true_positive.set(l, true);    
                            for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l]) false_positive.set(l, false);    
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) false_negative.set(l, false);   
                       }
                    }
					// count the labels for true positive, false positive, false negative, total detected, total truth
					tp_size[n] = 0.0f;
					fp_size[n] = 0.0f;
					fn_size[n] = 0.0f;
					seg_size[n] = 0.0f;
					tpl_size[n] = 0.0f;
					for (int l=1;l<list1.length;l++) {
					    float vol = ObjectStatistics.volume(label1, list1[l], nx,ny,nz);
					    seg_size[n] += vol;
					    if (false_positive.get(l)) fp_size[n] += vol;
					}   
					for (int l=1;l<list2.length;l++) {
					    float vol = ObjectStatistics.volume(label2, list2[l], nx,ny,nz);
					    tpl_size[n] += vol;
					    if (false_negative.get(l)) fn_size[n] += vol;
					    if (true_positive.get(l)) tp_size[n] += vol;
					}
					// normalize
					tp_size[n] /= true_positive.length();
					fp_size[n] /= false_positive.length();
					fn_size[n] /= false_negative.length();
					seg_size[n] /= list1.length;
					tpl_size[n] /= list2.length;
				}
				line = "True_positive_mean_cluster_size"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+tp_size[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "False_positive_mean_cluster_size"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fp_size[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "False_negative_mean_cluster_size"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+fn_size[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Segmented_mean_cluster_size"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+seg_size[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
				line = "Template_mean_cluster_size"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+tpl_size[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Cluster_maps")) {
				// compare: requires same labels (use the reference as basis)
				float[] trupos = new float[nlabels];
				float[] falpos = new float[nlabels];
				float[] falneg = new float[nlabels];
				int[][][] clustermap = new int[nx][ny][nz];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					    if (and[x][y][z]) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) true_positive.set(l, true);    
                            for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l]) false_positive.set(l, false);    
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) false_negative.set(l, false);   
                       }
                    }
                    for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
                        if (label2[x][y][z]>0) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l] && true_positive.get(l))
                                clustermap[x][y][z] = 3*n+3;
                        }
                        if (label1[x][y][z]>0) {
                            for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l] && false_positive.get(l))
                                clustermap[x][y][z] = 3*n+2;
                        }
                        if (label2[x][y][z]>0) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l] && false_negative.get(l))
                                clustermap[x][y][z] = 3*n+1;
                        }
                    }
				}
                ImageDataInt bufferData = new ImageDataInt(clustermap);
                bufferData.setHeader(header);
                bufferData.setName(name+"_clustermap");
                outputImage.setValue(bufferData);
			}
			if (statistics.get(s).equals("Average_surface_distance")) {
				// compare: requires same labels (use the reference as basis)
				float[] dist = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					if (vol1>0 && vol2>0)
						dist[n] = ObjectStatistics.averageSurfaceDistance(obj1,obj2,rx,ry,rz,nx,ny,nz);
					else
						dist[n] = -1.0f;
					//System.out.println("asd "+n+": "+dist[n]);
				}
				line = "Average_surface_distance"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+dist[n]);
				line+=("\n");
				output.add(line);
				System.out.println(line);
			}
			if (statistics.get(s).equals("Average_surface_difference")) {
				// compare: requires same labels (use the reference as basis)
				float[] dist = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					if (vol1>0 && vol2>0)
						dist[n] = ObjectStatistics.averageSurfaceDifference(obj1,obj2,rx,ry,rz,nx,ny,nz);
					else
						dist[n] = -1.0f;
					//System.out.println("asd "+n+": "+dist[n]);
				}
				line = "Average_surface_difference"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+dist[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Average_squared_surface_distance")) {
				// compare: requires same labels (use the reference as basis)
				float[] dist = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					if (vol1>0 && vol2>0)
						dist[n] = ObjectStatistics.averageSquaredSurfaceDistance(obj1,obj2,rx,ry,rz,nx,ny,nz);
					else
						dist[n] = -1.0f;
					//System.out.println("assd "+n+": "+dist[n]);
				}
				line = "Average_squared_surface_distance"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+dist[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Hausdorff_distance")) {
				// compare: requires same labels (use the reference as basis)
				float[] hausdorff = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					if (vol1>0 && vol2>0)
						hausdorff[n] = ObjectStatistics.hausdorffDistance(obj1,obj2,rx,ry,rz,nx,ny,nz);
					else
						hausdorff[n] = -1.0f;
					//System.out.println("hd "+n+": "+hausdorff[n]);
				}
				line = "Hausdorff_distance"+imgtag+reftag+notag;
				for (int n=0;n<nlabels;n++) line+=(delim+hausdorff[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Mean_intensity")) {
				float[] mean = new float[nlabels];
				float[] den = new float[nlabels];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int n=0;n<nlabels;n++) {
					mean[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						mean[n] += intensity[x][y][z];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) mean[n] = mean[n]/den[n];
				}
				line = "Mean_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+mean[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Std_intensity")) {
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
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						mean[n] += intensity[x][y][z];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) mean[n] = mean[n]/den[n];
				}
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						std[n] += (intensity[x][y][z]-mean[n])*(intensity[x][y][z]-mean[n]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) std[n] = (float)Math.sqrt(std[n]/den[n]);
				}
				line = "Std_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+std[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("50_intensity")) {
				float Imin = intensity[0][0][0];
				float Imax = intensity[0][0][0];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
					if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
				}
				int Nbins = 100;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
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
				line = "50_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("90_intensity")) {
				float Imin = intensity[0][0][0];
				float Imax = intensity[0][0][0];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
					if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
				}
				int Nbins = 100;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
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
					while (count<0.9f*hist[n][Nbins]) {
						count += hist[n][bin];
						bin++;
					}
					per[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
				}
				line = "90_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("10_intensity")) {
				float Imin = intensity[0][0][0];
				float Imax = intensity[0][0][0];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
					if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
				}
				int Nbins = 100;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
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
					while (count<0.1f*hist[n][Nbins]) {
						count += hist[n][bin];
						bin++;
					}
					per[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
				}
				line = "10_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("75_intensity")) {
				float Imin = intensity[0][0][0];
				float Imax = intensity[0][0][0];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
					if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
				}
				int Nbins = 100;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
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
					while (count<0.75f*hist[n][Nbins]) {
						count += hist[n][bin];
						bin++;
					}
					per[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
				}
				line = "75_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("25_intensity")) {
				float Imin = intensity[0][0][0];
				float Imax = intensity[0][0][0];
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					if (intensity[x][y][z]> Imax) Imax = intensity[x][y][z];
					if (intensity[x][y][z]< Imin) Imin = intensity[x][y][z];
				}
				int Nbins = 100;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam.getValue().booleanValue();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					for (int n=0;n<nlabels;n++) if (segmentation[x][y][z]==lbid[n] && (!ignoreZero || intensity[x][y][z]!=0) ) {
						int bin = Numerics.floor((intensity[x][y][z]-Imin)/(Imax-Imin)*Nbins);
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
					while (count<0.25f*hist[n][Nbins]) {
						count += hist[n][bin];
						bin++;
					}
					per[n] = Imin + (bin-1)/(float)Nbins*(Imax-Imin);
				}
				line = "25_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Cluster_maps")) {
			    byte[][][][] map = new byte[nx][ny][nz][nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[][][] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					// not needed here... int[][][] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[][][] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],lbid[n],ObjectExtraction.SUPEQUAL,ObjectExtraction.INFEQUAL);
					int[][][] label2 = ObjectLabeling.connected18Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[][][] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					    if (and[x][y][z]) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) map[x][y][z][n] = 3;   	
                        } else if (obj1[x][y][z]) {
                            for (int l=1;l<list1.length;l++) if (label1[x][y][z]==list1[l]) map[x][y][z][n] = 2;  
                        } else if (obj2[x][y][z]) {
                            for (int l=1;l<list2.length;l++) if (label2[x][y][z]==list2[l]) map[x][y][z][n] = 1;
                        }
					}
					ImageDataUByte outputData = new ImageDataUByte(map);		
                    outputData.setHeader(segImg.getHeader());
                    outputData.setName(segImg.getName()+"_cluster_map");
                    outputImage.setValue(outputData);
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
