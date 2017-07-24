package de.mpg.cbs.jist.surface;

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
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.io.FileExtensionFilter;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceMgdmRepresentation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume labelImage;
	private ParamVolume functionImage;
	private ParamFile atlasParam;
	private ParamInteger 	mgdmParam;
	private ParamBoolean	zeroParam;
	private ParamFloat 	distParam;
	//private ParamOption 	typeParam;
	//private static final String[] actionTypes = {"rebuild"};
	
	private ParamVolume mgdmlabelImage;
	private ParamVolume mgdmfunctImage;
		
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(labelImage = new ParamVolume("Label Image"));
		inputParams.add(functionImage = new ParamVolume("Function Image (opt)"));
		functionImage.setMandatory(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		atlasParam.setMandatory(false);
        
		inputParams.add(mgdmParam = new ParamInteger("MGDM degree", 1, 10, 3));
		inputParams.add(distParam = new ParamFloat("Max. distance (voxels)", -1, 10000, -1));
		distParam.setDescription("Maximum distance computed from the boundary (-1 if unlimited).");
		inputParams.add(zeroParam = new ParamBoolean("Skip zero label", true));
		
        //inputParams.add(typeParam = new ParamOption("Data type", functionTypes));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setLabel("Mgdm Representation");
		inputParams.setName("MgdmRepresentation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Builds a multi-object geometric deformable model (MGDM) representation for a collection of surfaces.");
		
		info.setVersion("3.0.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mgdmlabelImage = new ParamVolume("MGDM Labels Image",VoxelType.BYTE,-1,-1,-1,-1));
		outputParams.add(mgdmfunctImage = new ParamVolume("MGDM Function Image",VoxelType.FLOAT,-1,-1,-1,-1));
		
		outputParams.setName("MGDM representation");
		outputParams.setLabel("MGDM representation");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		ImageDataByte	labelImg = new ImageDataByte(labelImage.getImageData());
		int nx = labelImg.getRows();
		int ny = labelImg.getCols();
		int nz = labelImg.getSlices();
		int nxyz = nx*ny*nz;
		BasicInfo.displayMessage("Image dims: "+nx+", "+ny+", "+nz+"\n");
		float rx = labelImg.getHeader().getDimResolutions()[0];
		float ry = labelImg.getHeader().getDimResolutions()[1];
		float rz = labelImg.getHeader().getDimResolutions()[2];
		
		int orient = labelImg.getHeader().getImageOrientation().ordinal();
		int orx = labelImg.getHeader().getAxisOrientation()[0].ordinal();
		int ory = labelImg.getHeader().getAxisOrientation()[1].ordinal();
		int orz = labelImg.getHeader().getAxisOrientation()[2].ordinal();
		
		byte[] label = new byte[nx*ny*nz];
		byte[][][] buffer;
		buffer = labelImg.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			label[xyz] = buffer[x][y][z];
		}
		buffer = null;

		float[] function = null;
		if (functionImage.getImageData()!=null) {
			ImageDataFloat functImg = new ImageDataFloat(functionImage.getImageData());
			function = new float[nx*ny*nz];
			float[][][] fbuffer;
			fbuffer = functImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				function[xyz] = fbuffer[x][y][z];
			}
			buffer = null;
		}
		
		byte[] objlabel = null;
		int nobj = -1;
		if (atlasParam.getValue()!=null) {
			BasicInfo.displayMessage("Load atlas\n");
	
			SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam.getValue().getAbsolutePath());
		
			objlabel = atlas.getLabels();	
		} else {
			// buld from label image
			objlabel = ObjectLabeling.listOrderedLabels(label, nx, ny, nz);
		}
		nobj = objlabel.length;
		
		int nmgdm = mgdmParam.getValue().intValue();
		boolean skip = zeroParam.getValue().booleanValue();
		float dist = distParam.getValue().floatValue();
		
		if (skip && objlabel[0]==0) {
			byte[] tmp = new byte[nobj-1];
			for (int n=0;n<nobj-1;n++) tmp[n] = objlabel[n+1];
			objlabel = tmp;
			nobj--;
		}
		
		MgdmRepresentation mgdm = new MgdmRepresentation(label, function, nx, ny, nz, rx, ry, rz,
															objlabel, nobj, nmgdm, skip, dist); 
			
		// outputs
		String imgname;
		imgname = labelImg.getName();
		
		byte[][][][] lbl = mgdm.exportLabels();
		ImageDataByte lblData = new ImageDataByte(lbl);	
		lbl = null;
		lblData.setHeader(labelImg.getHeader());
		lblData.setName(labelImg.getName()+"_lbls");
		mgdmlabelImage.setValue(lblData);
		lblData = null;
		BasicInfo.displayMessage("labels");
		
		float[][][][] fun = mgdm.exportFunctions();
		ImageDataFloat funData = new ImageDataFloat(fun);	
		fun = null;
		funData.setHeader(labelImg.getHeader());
		funData.setName(labelImg.getName()+"_func");
		mgdmfunctImage.setValue(funData);
		funData = null;
		BasicInfo.displayMessage("functions");

		return;
	}
	
}
