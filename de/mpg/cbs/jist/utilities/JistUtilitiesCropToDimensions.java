package de.mpg.cbs.jist.utilities;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamPointInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import javax.vecmath.Point3i;

import de.mpg.cbs.utilities.CropParameters;
import de.mpg.cbs.utilities.CubicVolumeCropper;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

/*
 * @author Christine Lucas Tardif (ctardif@cbs.mpg.de)
 */
public class JistUtilitiesCropToDimensions extends ProcessingAlgorithm{
	
	private ParamVolume inputVol;
	
	private ParamVolume outputVol;

	private ParamInteger xmin;
	private ParamInteger xmax;
	private ParamInteger ymin;
	private ParamInteger ymax;
	private ParamInteger zmin;
	private ParamInteger zmax;
	
	private ParamString roiName;
	
	private ParamPointInteger offset;
	private ParamPointInteger dim;
	private ParamPointInteger dimNew;


	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputVol = new ParamVolume("Input Volume"));
		
		inputParams.add(xmin = new ParamInteger("X min", 0, Integer.MAX_VALUE, 3));
		inputParams.add(xmax = new ParamInteger("X max", 0, Integer.MAX_VALUE, 3));
		inputParams.add(ymin = new ParamInteger("Y min", 0, Integer.MAX_VALUE, 3));
		inputParams.add(ymax = new ParamInteger("Y max", 0, Integer.MAX_VALUE, 3));
		inputParams.add(zmin = new ParamInteger("Z min", 0, Integer.MAX_VALUE, 3));
		inputParams.add(zmax = new ParamInteger("Z max", 0, Integer.MAX_VALUE, 3));
		
		inputParams.add(roiName = new ParamString("ROI Name"));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Utilities");
		inputParams.setLabel("Crop To Dimensions");
		inputParams.setName("CropToDimensions");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Crops a 3D image to the entered voxel dimensions");
		
		info.setVersion("3.0.6");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}


	protected void createOutputParameters(ParamCollection outputParams) {

		outputParams.add(outputVol = new ParamVolume("Cropped Image",null,-1,-1,-1,-1));
	
		offset=new ParamPointInteger("Offset");
		dim=new ParamPointInteger("Uncropped Dimensions");
		dimNew=new ParamPointInteger("Cropped Dimensions");
		
		outputParams.add(offset);
		outputParams.add(dim);
		outputParams.add(dimNew);

		outputParams.setName("CroppedImage");
		outputParams.setLabel("Cropped Image");
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		
		ImageDataFloat inVol = new ImageDataFloat(inputVol.getImageData());
		
		int nx = inVol.getRows();
		int ny = inVol.getCols();
		int nz = inVol.getSlices();
		int nt = inVol.getComponents();
		int nxyz = nx*ny*nz;
		float rx = inVol.getHeader().getDimResolutions()[0];
		float ry = inVol.getHeader().getDimResolutions()[1];
		float rz = inVol.getHeader().getDimResolutions()[2];
		
		if (nt==1) {
			float[][][] in = inVol.toArray3d();
			
			int newnx = xmax.getInt()-xmin.getInt();
			int newny = ymax.getInt()-ymin.getInt();
			int newnz = zmax.getInt()-zmin.getInt();	
			
			float[][][] out = new float[newnx][newny][newnz];
			
			for (int x=0;x<newnx;x++) for (int y=0;y<newny;y++) for (int z=0;z<newnz;z++) {
				out[x][y][z] = in[x+xmin.getInt()][y+ymin.getInt()][z+zmin.getInt()];
			}
			in = null;
			inVol = null;
			
			offset.setValue(new Point3i(xmin.getInt(),ymin.getInt(),zmin.getInt()));
			dimNew.setValue(new Point3i(newnx,newny,newnz));
			dim.setValue(new Point3i(nx,ny,nz));
	
			
			ImageDataFloat outVol = new ImageDataFloat(out);		
			outVol.setHeader(inputVol.getImageData().getHeader());
			outVol.setName(inputVol.getImageData().getName()+"_"+roiName.getValue());
			outputVol.setValue(outVol);
			out = null;
			outVol = null;
		} else {
			float[][][][] in = inVol.toArray4d();
			
			int newnx = xmax.getInt()-xmin.getInt();
			int newny = ymax.getInt()-ymin.getInt();
			int newnz = zmax.getInt()-zmin.getInt();	
			
			float[][][][] out = new float[newnx][newny][newnz][nt];
			
			for (int x=0;x<newnx;x++) for (int y=0;y<newny;y++) for (int z=0;z<newnz;z++) for (int t=0;t<nt;t++) {
				out[x][y][z][t] = in[x+xmin.getInt()][y+ymin.getInt()][z+zmin.getInt()][t];
			}
			in = null;
			inVol = null;
			
			offset.setValue(new Point3i(xmin.getInt(),ymin.getInt(),zmin.getInt()));
			dimNew.setValue(new Point3i(newnx,newny,newnz));
			dim.setValue(new Point3i(nx,ny,nz));
	
			ImageDataFloat outVol = new ImageDataFloat(out);		
			outVol.setHeader(inputVol.getImageData().getHeader());
			outVol.setName(inputVol.getImageData().getName()+"_"+roiName.getValue());
			outputVol.setValue(outVol);
			out = null;
			outVol = null;
		}			
	}
}
