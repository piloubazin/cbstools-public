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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurfaceCollection;
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

import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.graphics.smooth.SurfaceInflate;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;

import javax.vecmath.Point3f;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceMgdmToMeshes extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume labelImage;
	private ParamVolume functionImage;
	//private ParamFile atlasParam;
	//private ParamInteger 	mgdmParam;
	private ParamBoolean	zeroParam;
	//private ParamFloat 	distParam;
	ParamOption		topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6"};
	//private ParamOption 	typeParam;
	//private static final String[] actionTypes = {"rebuild"};
	
	//private ParamVolume mgdmlabelImage;
	//private ParamVolume mgdmfunctImage;
	
	private ParamSurfaceCollection meshSurfaces;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(labelImage = new ParamVolume("Label Image"));
		inputParams.add(functionImage = new ParamVolume("Function Image"));
		//functionImage.setMandatory(false);
		
		//inputParams.add(atlasParam = new ParamFile("Atlas file (opt)",new FileExtensionFilter(new String[]{"txt"})));
		//atlasParam.setMandatory(false);
        
		inputParams.add(zeroParam = new ParamBoolean("Skip background label", true));
		
        //inputParams.add(typeParam = new ParamOption("Data type", functionTypes));
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("6/18");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setLabel("Mgdm to mesh collection");
		inputParams.setName("Mgdm to mesh collection");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Builds a multi-object geometric deformable model (MGDM) representation for a collection of surfaces.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(meshSurfaces=new ParamSurfaceCollection("Surface Meshes"));
		
		outputParams.setName("MGDM to Meshes");
		outputParams.setLabel("MGDM to Meshes");
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
		//if (functionImage.getImageData()!=null) {
			ImageDataFloat functImg = new ImageDataFloat(functionImage.getImageData());
			function = new float[nx*ny*nz];
			float[][][] fbuffer;
			fbuffer = functImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				function[xyz] = fbuffer[x][y][z];
				if (function[xyz]<0) function[xyz] *= -1.0f;
			}
			buffer = null;
		//}
		
		byte[] objlabel = null;
		int nobj = -1;
		/*
		if (atlasParam.getValue()!=null) {
			BasicInfo.displayMessage("Load atlas\n");
	
			SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam.getValue().getAbsolutePath());
		
			objlabel = atlas.getLabels();	
		} else {
			// buld from label image
			objlabel = ObjectLabeling.listOrderedLabels(label, nx, ny, nz);
		}
		*/
		objlabel = ObjectLabeling.listOrderedLabels(label, nx, ny, nz);
		nobj = objlabel.length;
		
		//int nmgdm = mgdmParam.getValue().intValue();
		boolean skip = zeroParam.getValue().booleanValue();
		//float dist = distParam.getValue().floatValue();
		
		if (skip) {
		//if (objlabel[0]==0) {
			byte[] tmp = new byte[nobj-1];
			for (int n=0;n<nobj-1;n++) tmp[n] = objlabel[n+1];
			objlabel = tmp;
			nobj--;
		}
		/*
		MgdmRepresentation mgdm = new MgdmRepresentation(label, function, nx, ny, nz, rx, ry, rz,
															objlabel, nobj, nmgdm, skip, dist); 
		*/
		
		// surface extraction
		int conn=0, cObj=0, cBg=0;
			 if (topologyParam.getValue().equals("6/18"))	{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }
		else if (topologyParam.getValue().equals("6/26"))	{ conn=ConnectivityRule.CONNECT_6_26; cObj = 6; cBg = 26; }
		else if (topologyParam.getValue().equals("18/6"))	{ conn=ConnectivityRule.CONNECT_18_6; cObj = 18; cBg = 6; }
		else if (topologyParam.getValue().equals("26/6"))	{ conn=ConnectivityRule.CONNECT_26_6; cObj = 26; cBg = 6; }
		else 												{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }

		// for each label, create the levelset then generate a mesh
		float[][][] levelset = new float[nx][ny][nz];
		for (int n=0;n<nobj;n++) {
			BasicInfo.displayMessage("extract object: "+objlabel[n]+"\n");
			// extract corresponding region	
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				if (label[xyz]==objlabel[n]) {
					levelset[x][y][z] = -function[xyz];
				} else {
					levelset[x][y][z] = function[xyz];
				}
			}
			ImageDataFloat levelsetImg = new ImageDataFloat(levelset);
			levelsetImg.setName(labelImg.getName()+"_lb"+objlabel[n]);
			
			IsoSurfaceOnGrid surfGen=new IsoSurfaceOnGrid();
			monitor.observe(surfGen);		
			EmbeddedSurface surf=surfGen.solveOriginal(levelsetImg,conn,-1e-9f,true);
			surf.scaleVertices(new float[]{rx, ry, rz});
			
			surf.scaleVertices(new float[]{1.0f, -1.0f, -1.0f});
			Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
			surf.translate(pt);
				
			surf.computeNormals();
			
			// weird discrepancy with the inflation module??
			//SurfaceInflate inflate=new SurfaceInflate(surf);
			
			// outputs
			meshSurfaces.add(surf);
		}
		
		return;
	}
	
}
