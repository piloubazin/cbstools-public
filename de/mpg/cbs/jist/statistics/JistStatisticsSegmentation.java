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
import de.mpg.cbs.core.statistics.*;

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
		
	private StatisticsSegmentation algorithm;
	
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
		
		algorithm = new StatisticsSegmentation();
		
		inputParams.setPackage(algorithm.getPackage());
		inputParams.setCategory(algorithm.getCategory());
		inputParams.setLabel(algorithm.getLabel());
		inputParams.setName(algorithm.getName());

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(References.getAuthor(algorithm.getAlgorithmAuthors()[0]));
		info.setAffiliation(algorithm.getAffiliation());
		info.setDescription(algorithm.getDescription());
		
		info.setVersion(algorithm.getVersion());
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
		
		// main algorithm
		algorithm = new StatisticsSegmentation();

		// i/o variables
		String name = Interface.getName(segImage);
		ImageHeader header = Interface.getHeader(segImage);
		
		int[] dims = Interface.getDimensions(segImage);
		algorithm.setDimensions(dims);
		algorithm.setResolutions(Interface.getResolutions(segImage));
		
		algorithm.setSegmentationImage(Interface.getIntegerImage3D(segImage));
		algorithm.setSegmentationName(name);
       
		if (Interface.isValid(intensImage)) {
            algorithm.setIntensityImage(Interface.getFloatImage3D(intensImage));
            algorithm.setIntensityName(Interface.getName(intensImage));
		}
        if (Interface.isValid(tempImage)) {
            algorithm.setTemplateImage(Interface.getIntegerImage3D(tempImage));
            algorithm.setTemplateName(Interface.getName(tempImage));
        }
        
        algorithm.setAtlasFile(atlasParam.getValue().getAbsolutePath());
        algorithm.setSkipFirstLabel(removeFirstParam.getValue().booleanValue());
		algorithm.setIgnoreZeroIntensities(ignoreZeroParam.getValue().booleanValue());
		
		algorithm.setSpreadsheetDirectory(statsParam.getValue().getAbsolutePath());
		
		algorithm.setStatistic1(stat1Param.getValue());
		algorithm.setStatistic2(stat2Param.getValue());
		algorithm.setStatistic3(stat3Param.getValue());
		
        algorithm.execute();
		
		// outputs
		if (algorithm.getOutputImage()!=null)
            Interface.setFloatImage3D(algorithm.getOutputImage(), dims, outputImage, name+"_stat-map", header);
		outputParam.setValue(new File(algorithm.getOutputFile()));
		
		return;
	}

}
