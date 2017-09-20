package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.graphics.smooth.SurfaceInflate;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import javax.vecmath.Point3f;

/*
 * @author Pierre-Louis Bazin and Christine L Tardif
 */
public class JistSurfaceLevelsetToMesh extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume levelsetImage;
	ParamBoolean 			mipavTransform;
	
	ParamOption		topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6"};
	
	private ParamSurface origSurface;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImage = new ParamVolume("Surface Level Set Image"));
		//inputParams.add(mipavTransform=new ParamBoolean("Align mesh to MIPAV image space", true));		
				
			
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("6/18");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setLabel("Level Set To Mesh");
		inputParams.setName("LevelSetToMesh");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Extracts a surface mesh from a level set function.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(origSurface=new ParamSurface("Surface Mesh"));

		outputParams.setName("extracted surface mesh");
		outputParams.setLabel("extracted surface mesh");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// import the image data into 1D arrays
		
		ImageDataFloat	levelsetImg = new ImageDataFloat(levelsetImage.getImageData());
		
		int nx = levelsetImg.getRows();
		int ny = levelsetImg.getCols();
		int nz = levelsetImg.getSlices();
		int nxyz = nx*ny*nz;
		float rx = levelsetImg.getHeader().getDimResolutions()[0];
		float ry = levelsetImg.getHeader().getDimResolutions()[1];
		float rz = levelsetImg.getHeader().getDimResolutions()[2];
		//float ox = levelsetImg.getHeader().getOrigin()[0];
		//float oy = levelsetImg.getHeader().getOrigin()[1];
		//float oz = levelsetImg.getHeader().getOrigin()[2];
		
		// surface extraction?
		int conn=0, cObj=0, cBg=0;
			 if (topologyParam.getValue().equals("6/18"))	{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }
		else if (topologyParam.getValue().equals("6/26"))	{ conn=ConnectivityRule.CONNECT_6_26; cObj = 6; cBg = 26; }
		else if (topologyParam.getValue().equals("18/6"))	{ conn=ConnectivityRule.CONNECT_18_6; cObj = 18; cBg = 6; }
		else if (topologyParam.getValue().equals("26/6"))	{ conn=ConnectivityRule.CONNECT_26_6; cObj = 26; cBg = 6; }
		else 												{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }

		IsoSurfaceOnGrid surfGen=new IsoSurfaceOnGrid();
		monitor.observe(surfGen);		
		EmbeddedSurface surf=surfGen.solveOriginal(levelsetImg,conn,-1e-9f,true);
		surf.scaleVertices(new float[]{rx, ry, rz});
		
		//if (mipavTransform.getValue().booleanValue()) {
		surf.scaleVertices(new float[]{1.0f, -1.0f, -1.0f});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		surf.translate(pt);
		//}
		
		surf.computeNormals();
		
		// weird discrepancy with the inflation module??
		SurfaceInflate inflate=new SurfaceInflate(surf);
		
		// outputs
		origSurface.setValue(surf);
	}


}
