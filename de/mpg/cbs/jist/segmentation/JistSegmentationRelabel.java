package de.mpg.cbs.jist.segmentation;

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

import gov.nih.mipav.view.*;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSegmentationRelabel extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	segImage;
	private ParamFile 		labelsParam;
	
	private ParamVolume		outputImage;
		
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(segImage = new ParamVolume("Segmentation Image"));
		inputParams.add(labelsParam = new ParamFile("Re-labeling list",new FileExtensionFilter(new String[]{"txt"})));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Segmentation");
		inputParams.setLabel("Relabel Segmentation");
		inputParams.setName("RelabelSegmentation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Relabel a segmentation image based on a list of label correspondences.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputImage = new ParamVolume("Relabeled Segmentation",VoxelType.INT));
		
		outputParams.setName("relabeling result");
		outputParams.setLabel("relabeling result");
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
		
		int[][][] segmentation = segImg.toArray3d();
		
		// get label list from atlas or segmentation
		
		int nlb = 0;
		String[] lbnames = null;
		int[] lbid = null;
		int[] lbnew = null;
		
		BasicInfo.displayMessage("Load labeling file\n");
	
		try {
            File f = new File(labelsParam.getValue().getAbsolutePath());
			String dir = f.getParent();
            FileReader fr = new FileReader(f);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
			StringTokenizer st;
			String imageFile, objLabelFile;
            // Exact corresponding template
            if (!line.equals("Relabeling File (edit at your own risks)")) {
                System.out.println("not a proper Relabeling file");
                br.close();
                fr.close();
                return;
            }
			line = br.readLine();
			line = br.readLine();
			if (line.startsWith("Structures")) {
				// Syntax:	lbid	name	lbnew
				st = new StringTokenizer(line, "	");
				st.nextToken();
				nlb = MipavUtil.getInt(st);
				lbnames = new String[nlb];
				lbid = new int[nlb];
				lbnew = new int[nlb];
				for (int n=0;n<nlb;n++) {
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					lbid[n] = MipavUtil.getInt(st);
					lbnames[n] = st.nextToken();
					lbnew[n] = MipavUtil.getInt(st);
					System.out.println("structure "+lbid[n]+": "+lbnames[n]+" -> "+lbnew[n]);
				}
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

        // change the labels
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
        	boolean found = false;
        	for (int n=0;n<nlb && !found;n++) {
        		if (segmentation[x][y][z]==lbid[n]) {
        			segmentation[x][y][z] = lbnew[n];
        			found = true;
        		}
        	}
        	if (!found)	segmentation[x][y][z] = 0;
		}	
		
		ImageDataInt segData = new ImageDataInt(segmentation);	
		segmentation = null;
		segData.setHeader(segImg.getHeader());
		segData.setName(segImg.getName()+"_relb");
		outputImage.setValue(segData);
		segData = null;
		
		return;
	}

}
