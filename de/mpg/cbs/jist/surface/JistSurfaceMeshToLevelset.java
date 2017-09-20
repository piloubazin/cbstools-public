package de.mpg.cbs.jist.surface;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;
import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.ResampleLevelSet;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.topology.TopologyCorrection;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.cs.cisst.algorithms.geometry.surface.*;
import edu.jhu.cs.cisst.algorithms.segmentation.gac.DistanceField3D;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;


public class JistSurfaceMeshToLevelset  extends ProcessingAlgorithm{
	ParamVolume origVol;
	ParamSurface origSurf;
	ParamSurface resultSurf;
	ParamVolume result;
	//ParamOption connectivity;
	ParamBoolean mipavTransform;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(origSurf=new ParamSurface("Surface Mesh"));
		inputParams.add(origVol=new ParamVolume("Reference Volume"));
		origVol.setDescription("Reference volume to use for level set representation dimensions.");
		//inputParams.add(connectivity=new ParamOption("Connectivity (Foreground,Background)",new String[]{"(18,6)","(6,18)","(26,6)","(6,26)"}));		
		//inputParams.add(mipavTransform=new ParamBoolean("Align to MIPAV image space", true));
		
		inputParams.setName("MeshToLevelSet");
		inputParams.setLabel("Mesh to Level Set");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		
		AlgorithmInformation info=getAlgorithmInformation();
		info.setWebsite("");
		info.setVersion(SurfaceToMask.getVersion());
		info.setDescription("Finds a volumetric level-set representation of an input surface.");
		
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(PrinceGroupAuthors.blakeLucas);
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(result=new ParamVolume("Level Set Image"));
		outputParams.add(resultSurf=new ParamSurface("Aligned Surface Mesh (opt)"));
		resultSurf.setMandatory(false);
	}
	protected void execute(CalculationMonitor monitor) {
		
		ImageDataFloat	origVolData = new ImageDataFloat(origVol.getImageData());
		int nx = origVolData.getRows();
		int ny = origVolData.getCols();
		int nz = origVolData.getSlices();
		float rx = origVolData.getHeader().getDimResolutions()[0];
		float ry = origVolData.getHeader().getDimResolutions()[1];
		float rz = origVolData.getHeader().getDimResolutions()[2];
		float ox = origVolData.getHeader().getOrigin()[0];
		float oy = origVolData.getHeader().getOrigin()[1];
		float oz = origVolData.getHeader().getOrigin()[2];
		
		EmbeddedSurface surf=origSurf.getSurface();
		
		/*
		// maps the surface to the voxel space
		if (mipavTransform.getValue().booleanValue()) {
			surf.scaleVertices(new float[]{-1.0f, 1.0f, -1.0f});
			Point3f pt0 = new Point3f((nx-1)*rx/2.0f, (ny-1)*ry/2.0f, (nz-1)*rz/2.0f);
			surf.translate(pt0);
		}
		*/
		EmbeddedSurface rsurf = surf.clone();
		
		// scale the surface data into voxel space
		rsurf.scaleVertices(new float[]{1.0f/rx, -1.0f/ry, -1.0f/rz});
		Point3f pt1 = new Point3f(0, (ny-1), (nz-1));
		rsurf.translate(pt1);
		
		/*
		SurfaceToMask surf2vol=new SurfaceToMask();
		monitor.observe(surf2vol);
		ImageDataFloat maskVol=surf2vol.solve(origVolData, rsurf);
		*/
		
		// from Blake Lucas's new code
		TriangleMesh mesh = new TriangleMesh(rsurf);
		MeshDistanceHash hash;
		hash = new MeshDistanceHash(mesh, 1.0,	new Point3f(0, 0, 0), new Point3f(nx, ny, nz));
		monitor.observe(hash);
		ImageDataFloat df = hash.getDistanceField();
		DistanceField3D dfs = new DistanceField3D();
		df = dfs.solve(df, 10.0f);
		result.setValue(df);
		//isoSurfParam.setValue(hash.getIsoSurface());
		//reshapeMatrixParam.setValue(hash.getTransform());

		df.setHeader(origVolData.getHeader());
		df.setName(surf.getName()+"_ls");
		result.setValue(df);
		
		if (mipavTransform.getValue().booleanValue()) {
			surf.setName(surf.getName()+"_al");
			resultSurf.setValue(surf);
		}
	}

}
