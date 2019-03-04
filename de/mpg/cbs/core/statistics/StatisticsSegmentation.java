package de.mpg.cbs.core.statistics;

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
public class StatisticsSegmentation {

	// jist containers
	private int[] 	segImage;
	private float[] intensImage = null;
	private int[] 	tempImage = null;
	
	private String segName = "segmentation";
	private String intensName = "intensity";
	private String tempName = "template";
	
	private String 	atlasParam = null;
	private boolean	removeFirstParam = false;
	private boolean	ignoreZeroParam = true;
	
	int nx, ny, nz, nc, nxyz;
	float rx, ry, rz;
	
	private String 		statsParam;
	
	private int nstat = 3;
	private String[] 	statParam = new String[nstat];
	private static final String[] statTypes = {"none", "Voxels", "Volume", "--- intensity ---", "Mean_intensity", "Std_intensity",
												"10_intensity","25_intensity","50_intensity","75_intensity","90_intensity",
												"Median_intensity","IQR_intensity","SNR_intensity","rSNR_intensity",
												"Boundary_gradient","Boundary_magnitude",
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
	
	private String 		outputParam;
    private float[] 	outputImage = null;
		
	private String delim = ",";
		
	
	// input parameters
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setIntensityImage(float[] val) { intensImage = val; }
	public final void setTemplateImage(int[] val) { tempImage = val; }
	public final void setSegmentationName(String val) { segName = val; }
	public final void setIntensityName(String val) { intensName = val; }
	public final void setTemplateName(String val) { tempName = val; }
	
	public final void setAtlasFile(String val) { atlasParam = val; }
	public final void setSkipFirstLabel(boolean val) { removeFirstParam = val; }
	public final void setIgnoreZeroIntensities(boolean val) { ignoreZeroParam = val; }
		
	public final void setSpreadsheetDirectory(String val) { statsParam = val+"/volume_statistics.csv"; }
	public final void setSpreadsheetFile(String val) { statsParam = val; }
		
	public final void setStatisticAt(int n, String val) {statParam[n] = val; }
	public final void setStatisticNumber(int  val) {nstat = val; statParam = new String[nstat]; }
	public final void setStatistic1(String val) {statParam[0] = val; }
	public final void setStatistic2(String val) {statParam[1] = val; }
	public final void setStatistic3(String val) {statParam[2] = val; }

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nc=1; nxyz=nx*ny*nz; }
	public final void setDimensions(int x, int y, int z, int c) { nx=x; ny=y; nz=z; nc=c; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { 
	    nx=dim[0]; ny=dim[1]; nz=dim[2]; 
	    if (dim.length>3) nc=dim[3]; else nc=1;
	    nxyz=nx*ny*nz; 
	}
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Statistics"; } 
	public final String getLabel() { return "Segmentation Statistics"; }
	public final String getName() { return "SegmentationStatistics"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Compute various statistics (geometric, intensity or compared to template labeling) for volumetric segmentation labels and outputs the results into an organized spreadsheet."; }

	public final String getVersion() { return "3.0.7"; }

	// outputs
    public final String getOutputFile() { return outputParam; }
    public final float[] getOutputImage() { return outputImage; }    
        
	public void execute() {
		
		// get label list from atlas or segmentation
		
		int nlabels;
		String[] lbname;
		int[] lbid;
		int[] segmentation = segImage;

		if (atlasParam!=null && atlasParam.length()>0) {
			BasicInfo.displayMessage("Load atlas\n");
	
			SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam);
			
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
		
		float[] intensity = intensImage;
		int[] template = tempImage;
		
		// Main algorithm
		
		String line;
		
		String imgtag,notag,reftag,inttag;
		imgtag = delim+segName;
		notag =  delim+"-";
		if (intensImage!=null) inttag = delim+intensName;
		else inttag = delim+"-";
		if (tempImage!=null) reftag = delim+tempName;
		else reftag = delim+"-";
			
		// standard computations for the object
		ArrayList<String> output = new ArrayList();

		// remove the zero label, if needed
		if (removeFirstParam) {
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
		//String lbline = "labels"+notag+notag+notag;
		//for (int n=0;n<nlabels;n++) lbline += (delim+lbname[n]);
		//lbline +=("\n");
		String lbline = "Measure"+delim+"Segmentation"+delim+"Template"+delim+"Intensity";
		for (int n=0;n<nlabels;n++) lbline += (delim+lbname[n]);
		lbline +=("\n");
		System.out.print(lbline);
		
		ArrayList<String> statistics = new ArrayList();
		for (int n=0;n<nstat;n++) 
            if (statParam[n]!=null) statistics.add(statParam[n]);
		
		// pre-compute redundant measures ?
		
		// compute the statistics
		for (int s=0; s<statistics.size(); s++) {
			System.out.print("Statistic: "+statistics.get(s)+"\n");
			if (statistics.get(s).equals("Voxels")) {
				// compute the volumes
				float[] volume = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[] obj = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
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
					boolean[] obj = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.XOR,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectLabeling.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] xor = ObjectGeometry.binaryOperation(obj1,obj2,ObjectLabeling.XOR,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] objd1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					boolean[] and1 = ObjectGeometry.binaryOperation(objd1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[] objd2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					boolean[] and2 = ObjectGeometry.binaryOperation(obj1,objd2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					obj2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					obj1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					obj2 = Morphology.dilateObject(obj2, nx,ny,nz, 1,1,1);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					obj1 = Morphology.dilateObject(obj1, nx,ny,nz, 1,1,1);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol2 = ObjectStatistics.volume(obj2,nx,ny,nz);
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					// not needed here... int[] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label2 = ObjectLabeling.connected18Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[] detected = new boolean[list2.length];
					for (int xyz=0;xyz<nxyz;xyz++) if (and[xyz]) {
					    for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) detected[l] = true;    	
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					boolean[] detected = new boolean[list1.length];
					for (int xyz=0;xyz<nxyz;xyz++) if (and[xyz]) {
					    for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l]) detected[l] = true;    	
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int xyz=0;xyz<nxyz;xyz++) if (and[xyz]) {
                        for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) true_positive.set(l, true);    
                        for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l]) false_positive.set(l, false);    
                        for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) false_negative.set(l, false);   
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					// not needed here... int[] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int xyz=0;xyz<nxyz;xyz++) if (and[xyz]) {
                        for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) true_positive.set(l, true);    
                        for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l]) false_positive.set(l, false);    
                        for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) false_negative.set(l, false);   
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
				float[] clustermap = new float[nxyz];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label1 = ObjectLabeling.connected26Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label2 = ObjectLabeling.connected26Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					BitSet true_positive = new BitSet(list2.length);
					BitSet false_positive = new BitSet(list1.length);
					BitSet false_negative = new BitSet(list2.length);
					true_positive.set(0, list2.length, false);
					false_positive.set(1, list1.length, true);
					false_negative.set(1, list2.length, true);
					for (int xyz=0;xyz<nxyz;xyz++) {
					    if (and[xyz]) {
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) true_positive.set(l, true);    
                            for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l]) false_positive.set(l, false);    
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) false_negative.set(l, false);   
                       }
                    }
                    for (int xyz=0;xyz<nxyz;xyz++) {
					    if (label2[xyz]>0) {
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l] && true_positive.get(l))
                                clustermap[xyz] = 3*n+3;
                        }
                        if (label1[xyz]>0) {
                            for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l] && false_positive.get(l))
                                clustermap[xyz] = 3*n+2;
                        }
                        if (label2[xyz]>0) {
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l] && false_negative.get(l))
                                clustermap[xyz] = 3*n+1;
                        }
                    }
				}
                outputImage = clustermap;
			}
			if (statistics.get(s).equals("Average_surface_distance")) {
				// compare: requires same labels (use the reference as basis)
				float[] dist = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
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
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					float vol1 = ObjectStatistics.volume(obj1,nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
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
				boolean ignoreZero = ignoreZeroParam;
				for (int n=0;n<nlabels;n++) {
					mean[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						mean[n] += intensity[xyz];
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
				boolean ignoreZero = ignoreZeroParam;
				for (int n=0;n<nlabels;n++) {
					mean[n] = 0.0f;
					std[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						mean[n] += intensity[xyz];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) mean[n] = mean[n]/den[n];
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						std[n] += (intensity[xyz]-mean[n])*(intensity[xyz]-mean[n]);
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
			if (statistics.get(s).equals("SNR_intensity")) {
				float[] mean = new float[nlabels];
				float[] std = new float[nlabels];
				float[] den = new float[nlabels];
				boolean ignoreZero = ignoreZeroParam;
				for (int n=0;n<nlabels;n++) {
					mean[n] = 0.0f;
					std[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						mean[n] += intensity[xyz];
						den[n] += 1.0f;
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) mean[n] = mean[n]/den[n];
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						std[n] += (intensity[xyz]-mean[n])*(intensity[xyz]-mean[n]);
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) std[n] = (float)Math.sqrt(std[n]/den[n]);
				}
				line = "SNR_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+mean[n]/std[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("50_intensity")) {
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
			if (statistics.get(s).equals("Median_intensity")) {
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
				line = "Median_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("90_intensity")) {
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
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
			if (statistics.get(s).equals("IQR_intensity")) {
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
						if (bin<0) bin = 0;
						if (bin>=Nbins) bin = Nbins-1;
						hist[n][bin]++;
						hist[n][Nbins]++;	// total number
					}
				}
				float[] per = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					float count = 0.0f;
					int binmin = 0;
					while (count<0.25f*hist[n][Nbins]) {
						count += hist[n][binmin];
						binmin++;
					}
					int binmax = binmin;
					while (count<0.75f*hist[n][Nbins]) {
						count += hist[n][binmax];
						binmax++;
					}
					per[n] = (binmax-binmin)/(float)Nbins*(Imax-Imin);
				}
				line = "IQR_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("rSNR_intensity")) {
				float Imin = intensity[0];
				float Imax = intensity[0];
				for (int xyz=0;xyz<nxyz;xyz++) {
					if (intensity[xyz]> Imax) Imax = intensity[xyz];
					if (intensity[xyz]< Imin) Imin = intensity[xyz];
				}
				int Nbins = 1000;
				float[][] hist = new float[nlabels][Nbins+1];
				boolean ignoreZero = ignoreZeroParam;
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
						int bin = Numerics.floor((intensity[xyz]-Imin)/(Imax-Imin)*Nbins);
						if (bin<0) bin = 0;
						if (bin>=Nbins) bin = Nbins-1;
						hist[n][bin]++;
						hist[n][Nbins]++;	// total number
					}
				}
				float[] per = new float[nlabels];
				for (int n=0;n<nlabels;n++) {
					float count = 0.0f;
					int binmin = 0;
					while (count<0.25f*hist[n][Nbins]) {
						count += hist[n][binmin];
						binmin++;
					}
					int binmed = binmin;
					while (count<0.5f*hist[n][Nbins]) {
						count += hist[n][binmed];
						binmed++;
					}
					int binmax = binmed;
					while (count<0.75f*hist[n][Nbins]) {
						count += hist[n][binmax];
						binmax++;
					}
					per[n] = (Imin + (binmed-1)/(float)Nbins*(Imax-Imin)) / ((binmax-binmin)/(float)Nbins*(Imax-Imin));
				}
				line = "rSNR_intensity"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+per[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Boundary_gradient")) {
				float[] grad = new float[nlabels];
				float[] den = new float[nlabels];
				boolean ignoreZero = ignoreZeroParam;
				for (int n=0;n<nlabels;n++) {
					grad[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
					    // look for neighbors
					    for (byte k=0;k<26;k++) {
					        int ngb = Ngb.neighborIndex(k, xyz, nx, ny, nz);
					        if ( (!ignoreZero || intensity[ngb]!=0) && segmentation[ngb]!=lbid[n]) {
					            grad[n] += (intensity[xyz]-intensity[ngb]);
                                den[n] += 1.0f;
                            }
                        }
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) grad[n] = grad[n]/den[n];
				}
				line = "Boundary_gradient"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+grad[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Boundary_magnitude")) {
				float[] mag = new float[nlabels];
				float[] den = new float[nlabels];
				boolean ignoreZero = ignoreZeroParam;
				for (int n=0;n<nlabels;n++) {
					mag[n] = 0.0f;
					den[n] = 0.0f;
				}
				for (int xyz=0;xyz<nxyz;xyz++) {
					for (int n=0;n<nlabels;n++) if (segmentation[xyz]==lbid[n] && (!ignoreZero || intensity[xyz]!=0) ) {
					    // look for neighbors
					    for (byte k=0;k<26;k++) {
					        int ngb = Ngb.neighborIndex(k, xyz, nx, ny, nz);
					        if ( (!ignoreZero || intensity[ngb]!=0) && segmentation[ngb]!=lbid[n]) {
					            mag[n] += Numerics.abs(intensity[xyz]-intensity[ngb]);
                                den[n] += 1.0f;
                            }
                        }
					}
				}
				for (int n=0;n<nlabels;n++) {
					if (den[n]>0) mag[n] = mag[n]/den[n];
				}
				line = "Boundary_magnitude"+imgtag+notag+inttag;
				for (int n=0;n<nlabels;n++) line+=(delim+mag[n]);
				line+=("\n");
				output.add(line);
			}
			if (statistics.get(s).equals("Cluster_maps")) {
			    float[] map = new float[nxyz*nlabels];
				for (int n=0;n<nlabels;n++) {
					// get the connected components
					boolean[] obj1 = ObjectExtraction.objectFromLabelImage(segmentation,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					int[] list1 = ObjectLabeling.listOrderedLabels(label1, nx,ny,nz);
					// not needed here... int[] label1 = ObjectLabeling.connected18Object3D(obj1, nx,ny,nz);
					boolean[] obj2 = ObjectExtraction.objectFromLabelImage(template,nx,ny,nz,lbid[n],ObjectExtraction.EQUAL);
					int[] label2 = ObjectLabeling.connected18Object3D(obj2, nx,ny,nz);
					int[] list2 = ObjectLabeling.listOrderedLabels(label2, nx,ny,nz);
					// use the intersection to find which ones overlap
					boolean[] and = ObjectGeometry.binaryOperation(obj1,obj2,ObjectExtraction.AND,nx,ny,nz);
					for (int xyz=0;xyz<nxyz;xyz++) {
					    if (and[xyz]) {
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) map[xyz+n*nxyz] = 3;   	
                        } else if (obj1[xyz]) {
                            for (int l=1;l<list1.length;l++) if (label1[xyz]==list1[l]) map[xyz+n*nxyz] = 2;  
                        } else if (obj2[xyz]) {
                            for (int l=1;l<list2.length;l++) if (label2[xyz]==list2[l]) map[xyz+n*nxyz] = 1;
                        }
					}
					outputImage = map;
				}
			}
		}	
		System.out.println("Result:");
		System.out.println(output);
		
		System.out.println("write image data to: "+statsParam);
		String filename = statsParam;
		
		// write the output to file
		addStatisticsToFile(filename, output, lbline);
				
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
	    if (main.size()<=0) main.add(lbline);
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
				// removed for cleaner .csv files
				// main.add(lbline);
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
            //if (!line.startsWith("CBSTools Volumetric Statistics File")) {
            //    System.out.println("not a proper CBSTools statistics file");
            //    br.close();
            //    fr.close();
            //    return null;
            //}
			//line = br.readLine();
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
			//pw.write("CBSTools Volumetric Statistics File\n");
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
