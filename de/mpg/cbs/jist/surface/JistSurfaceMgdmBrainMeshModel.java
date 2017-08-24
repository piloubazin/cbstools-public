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
import java.util.BitSet;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceMgdmBrainMeshModel extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume labelImage;
	private ParamVolume functionImage;
	private ParamFile atlasParam;
	//private ParamInteger 	mgdmParam;
	//private ParamBoolean	zeroParam;
	ParamOption		topologyParam;
	private ParamFloat 	distParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6"};
	//private ParamOption 	typeParam;
	//private static final String[] actionTypes = {"rebuild"};
	private ParamBoolean	simpleParam;
	
	//private ParamVolume mgdmlabelImage;
	//private ParamVolume mgdmfunctImage;
	
	private ParamSurfaceCollection meshSurfaces;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(labelImage = new ParamVolume("Label Image"));
		inputParams.add(functionImage = new ParamVolume("Function Image"));
		//functionImage.setMandatory(false);
		
		inputParams.add(atlasParam = new ParamFile("Atlas file",new FileExtensionFilter(new String[]{"txt"})));
		//atlasParam.setMandatory(false);
        
		//inputParams.add(zeroParam = new ParamBoolean("Skip zero label", true));
		
        //inputParams.add(typeParam = new ParamOption("Data type", functionTypes));
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("6/18");

		inputParams.add(distParam=new ParamFloat("Min edge distance",0,0.5f,0.2f));
		inputParams.add(simpleParam = new ParamBoolean("Simple mesh", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setLabel("Mgdm to brain mesh model");
		inputParams.setName("Mgdm to brain mesh model");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Turns a multi-object geometric deformable model (MGDM) representation into a surface-based brain model for simulations.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(meshSurfaces=new ParamSurfaceCollection("Surface Meshes"));
		
		outputParams.setName("MGDM to Brain Mesh Model");
		outputParams.setLabel("MGDM to Brain Mesh Model");
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
		
		// generate regions from brain atlas
		BasicInfo.displayMessage("Load atlas\n");
	
		SimpleShapeAtlas atlas = new SimpleShapeAtlas(atlasParam.getValue().getAbsolutePath());

		int maxlb = 0;
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			if (atlas.getLabels()[nobj]>maxlb) maxlb = atlas.getLabels()[nobj];
		}
		
		System.out.println("Extracting regions\n");
		
		BitSet isOuterShell = new BitSet(maxlb);
		BitSet isOuterBrainSurface = new BitSet(maxlb);
		BitSet isCerebrumSurface = new BitSet(maxlb);
		BitSet isCerebellumSurface = new BitSet(maxlb);
		BitSet isCRWhiteMatterSurface = new BitSet(maxlb);
		BitSet isCBWhiteMatterSurface = new BitSet(maxlb);
		BitSet isBrainstemSurface = new BitSet(maxlb);
		BitSet isSubcorticalSurface = new BitSet(maxlb);
		BitSet isVentricleSurface = new BitSet(maxlb);
		
		// use negative region definitions
		for (int nobj=0;nobj<atlas.getNumber();nobj++) {
			// Outer shell of CSF
				 if (atlas.getNames()[nobj].equals("Background")) isOuterShell.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Dura")) isOuterShell.set(atlas.getLabels()[nobj], false);
			//else if (atlas.getNames()[nobj].equals("Sinuses")) isOuterShell.set(atlas.getLabels()[nobj]);
			else isOuterShell.set(atlas.getLabels()[nobj], true);
			
			// Outer brain surface
				 if (!isOuterShell.get(atlas.getLabels()[nobj])) isOuterBrainSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Sulcal-CSF")) isOuterBrainSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Arteries")) isOuterBrainSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Sinuses")) isOuterBrainSurface.set(atlas.getLabels()[nobj], false);
			else isOuterBrainSurface.set(atlas.getLabels()[nobj], true);
			
			// cerebrum surface
				 if (!isOuterBrainSurface.get(atlas.getLabels()[nobj])) isCerebrumSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isCerebrumSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isCerebrumSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Brainstem")) isCerebrumSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Ventricle4")) isCerebrumSurface.set(atlas.getLabels()[nobj], false);
			else isCerebrumSurface.set(atlas.getLabels()[nobj], true);
			
			// WM surface
				 if (!isOuterBrainSurface.get(atlas.getLabels()[nobj])) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-GML")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-GMR")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-GM")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("CerebralGM")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Brainstem")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Ventricle4")) isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			else isCRWhiteMatterSurface.set(atlas.getLabels()[nobj], true);

			// cerebellum surface
				 if (atlas.getNames()[nobj].equals("Cerebellum-GM")) isCerebellumSurface.set(atlas.getLabels()[nobj], true);
			else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isCerebellumSurface.set(atlas.getLabels()[nobj], true);
			else isCerebellumSurface.set(atlas.getLabels()[nobj], false);
			
			// cerebellum WM surface
				 if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isCBWhiteMatterSurface.set(atlas.getLabels()[nobj], true);
			else isCBWhiteMatterSurface.set(atlas.getLabels()[nobj], false);
			
			// brainstem WM surface
				 if (atlas.getNames()[nobj].equals("Brainstem")) isBrainstemSurface.set(atlas.getLabels()[nobj], true);
			else isBrainstemSurface.set(atlas.getLabels()[nobj], false);
			
			// Subcortical surface (incl. ventricles)
				 if (!isCRWhiteMatterSurface.get(atlas.getLabels()[nobj])) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-WML")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-WMR")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("CerebralWM")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebrum-WM")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Cerebellum-WM")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("Brainstem")) isSubcorticalSurface.set(atlas.getLabels()[nobj], false);
			else isSubcorticalSurface.set(atlas.getLabels()[nobj], true);
			
			// Ventricles
				 if (!isSubcorticalSurface.get(atlas.getLabels()[nobj])) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("AmygdalaL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("CaudateL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("HippocampusL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("PutamenL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("ThalamusL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("GlobusPallidusL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("StriatumL")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("AmygdalaR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("CaudateR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("HippocampusR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("PutamenR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("ThalamusR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("GlobusPallidusR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("StriatumR")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("SubcorticalGM")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else if (atlas.getNames()[nobj].equals("SubCorticalGM")) isVentricleSurface.set(atlas.getLabels()[nobj], false);
			else isVentricleSurface.set(atlas.getLabels()[nobj], true);			
		}
		
		// surface extraction
		int conn=0, cObj=0, cBg=0;
			 if (topologyParam.getValue().equals("6/18"))	{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }
		else if (topologyParam.getValue().equals("6/26"))	{ conn=ConnectivityRule.CONNECT_6_26; cObj = 6; cBg = 26; }
		else if (topologyParam.getValue().equals("18/6"))	{ conn=ConnectivityRule.CONNECT_18_6; cObj = 18; cBg = 6; }
		else if (topologyParam.getValue().equals("26/6"))	{ conn=ConnectivityRule.CONNECT_26_6; cObj = 26; cBg = 6; }
		else 											{ conn=ConnectivityRule.CONNECT_6_18; cObj = 6; cBg = 18; }

		float dist = distParam.getFloat();
		
		if (simpleParam.getValue().booleanValue()) {
			// for each label, create the levelset then generate a mesh
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isOuterShell, labelImg.getName()+"_"+"shell",conn,dist));
			//meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isOuterBrainSurface, labelImg.getName()+"_"+"brain",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCerebrumSurface, labelImg.getName()+"_"+"cerebrum_gm",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCerebellumSurface, labelImg.getName()+"_"+"cerebellum_gm",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCRWhiteMatterSurface, labelImg.getName()+"_"+"cerebrum_wm",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCBWhiteMatterSurface, labelImg.getName()+"_"+"cerebellum_wm",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isBrainstemSurface, labelImg.getName()+"_"+"brainstem",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isSubcorticalSurface, labelImg.getName()+"_"+"subcortex",conn,dist));
			meshSurfaces.add(extractSimpleMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isVentricleSurface, labelImg.getName()+"_"+"ventricles",conn,dist));
			
		} else {
			// for each label, create the levelset then generate a mesh
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isOuterShell, labelImg.getName()+"_"+"shell",conn,dist));
			//meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isOuterBrainSurface, labelImg.getName()+"_"+"brain",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCerebrumSurface, labelImg.getName()+"_"+"cerebrum_gm",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCerebellumSurface, labelImg.getName()+"_"+"cerebellum_gm",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCRWhiteMatterSurface, labelImg.getName()+"_"+"cerebrum_wm",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isCBWhiteMatterSurface, labelImg.getName()+"_"+"cerebellum_wm",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isBrainstemSurface, labelImg.getName()+"_"+"brainstem",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isSubcorticalSurface, labelImg.getName()+"_"+"subcortex",conn,dist));
			meshSurfaces.add(extractMgdmSurface(function,label,nx,ny,nz,rx,ry,rz, isVentricleSurface, labelImg.getName()+"_"+"ventricles",conn,dist));
		}			
		return;
	}
	
	EmbeddedSurface extractMgdmSurface(float[] function, byte[] label, int nx, int ny, int nz, float rx, float ry, float rz, BitSet includedLabels, String name, int conn, float mindist) {
		float[][][] levelset = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			// check for distances
			if (function[xyz]<mindist) function[xyz] = mindist; 

			// build level set 	
			if (includedLabels.get(label[xyz]))  levelset[x][y][z] = -function[xyz];
			else levelset[x][y][z] = function[xyz];
		}
		ImageDataFloat levelsetImg = new ImageDataFloat(levelset);
		levelsetImg.setName(name);
			
		IsoSurfaceOnGrid surfGen=new IsoSurfaceOnGrid();
		//monitor.observe(surfGen);		
		EmbeddedSurface surf=surfGen.solveOriginal(levelsetImg,conn,-1e-9f,true);
		surf.scaleVertices(new float[]{rx, ry, rz});
			
		surf.scaleVertices(new float[]{1.0f, -1.0f, -1.0f});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		surf.translate(pt);
				
		surf.computeNormals();
			
		// weird discrepancy with the inflation module??
		SurfaceInflate inflate=new SurfaceInflate(surf);
		
		return surf;
	}
	
	EmbeddedSurface extractSimpleMgdmSurface(float[] function, byte[] label, int nx, int ny, int nz, float rx, float ry, float rz, BitSet includedLabels, String name, int conn, float mindist) {
	
		float[][][] levelset = new float[nx][ny][nz];

		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			// check for distances
			if (function[xyz]<mindist) function[xyz] = mindist; 

			// build level set 	
			if (includedLabels.get(label[xyz]))  levelset[x][y][z] = -function[xyz];
			else levelset[x][y][z] = function[xyz];
		}
		EmbeddedSurface surf = extractSimpleSurface(levelset, nx, ny, nz, rx, ry, rz);
		surf.setName(name);
		/*
		surf.scaleVertices(new float[]{rx, ry, rz});
		surf.scaleVertices(new float[]{1.0f, -1.0f, -1.0f});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		surf.translate(pt);
		*/
		surf.computeNormals();

		return surf;
	}
	
	EmbeddedSurface extractSimpleSurface(float[][][] levelset, int nx, int ny, int nz, float rx, float ry, float rz) {
		
		Point3f[] points;
		int[] indexes;
		
		int Nv = 0;
		int Ne = 0;
		int Nf = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (levelset[x][y][z]<0) {
				if (levelset[x+1][y][z]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
				if (levelset[x-1][y][z]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
				if (levelset[x][y+1][z]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
				if (levelset[x][y-1][z]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
				if (levelset[x][y][z+1]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
				if (levelset[x][y][z-1]>=0) {
					Nv+=4; Ne+=5; Nf+=2;
				}
			}
		}
		points = new Point3f[Nv];
		for (int v=0;v<Nv;v++) points[v] = new Point3f();
		indexes = new int[Nf*3];
		
		Point3f[] newpoints = new Point3f[4];
		for (int v=0;v<4;v++) newpoints[v] = new Point3f();
		int v=0;
		int e=0;
		int f=0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (levelset[x][y][z]<0) {
				// object: check for boundary
				// check the six directions
				// +X
				if (levelset[x+1][y][z]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x+0.5f; newpoints[0].y = y-0.5f; newpoints[0].z = z-0.5f;
					newpoints[1].x = x+0.5f; newpoints[1].y = y+0.5f; newpoints[1].z = z-0.5f;
					newpoints[2].x = x+0.5f; newpoints[2].y = y+0.5f; newpoints[2].z = z+0.5f;
					newpoints[3].x = x+0.5f; newpoints[3].y = y-0.5f; newpoints[3].z = z+0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
				// -X
				if (levelset[x-1][y][z]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x-0.5f; newpoints[0].y = y-0.5f; newpoints[0].z = z-0.5f;
					newpoints[1].x = x-0.5f; newpoints[1].y = y+0.5f; newpoints[1].z = z-0.5f;
					newpoints[2].x = x-0.5f; newpoints[2].y = y+0.5f; newpoints[2].z = z+0.5f;
					newpoints[3].x = x-0.5f; newpoints[3].y = y-0.5f; newpoints[3].z = z+0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
				// +Y
				if (levelset[x][y+1][z]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x-0.5f; newpoints[0].y = y+0.5f; newpoints[0].z = z-0.5f;
					newpoints[1].x = x+0.5f; newpoints[1].y = y+0.5f; newpoints[1].z = z-0.5f;
					newpoints[2].x = x+0.5f; newpoints[2].y = y+0.5f; newpoints[2].z = z+0.5f;
					newpoints[3].x = x-0.5f; newpoints[3].y = y+0.5f; newpoints[3].z = z+0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
				// -Y
				if (levelset[x][y-1][z]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x-0.5f; newpoints[0].y = y-0.5f; newpoints[0].z = z-0.5f;
					newpoints[1].x = x+0.5f; newpoints[1].y = y-0.5f; newpoints[1].z = z-0.5f;
					newpoints[2].x = x+0.5f; newpoints[2].y = y-0.5f; newpoints[2].z = z+0.5f;
					newpoints[3].x = x-0.5f; newpoints[3].y = y-0.5f; newpoints[3].z = z+0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
				// +Z
				if (levelset[x][y][z+1]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x-0.5f; newpoints[0].y = y-0.5f; newpoints[0].z = z+0.5f;
					newpoints[1].x = x+0.5f; newpoints[1].y = y-0.5f; newpoints[1].z = z+0.5f;
					newpoints[2].x = x+0.5f; newpoints[2].y = y+0.5f; newpoints[2].z = z+0.5f;
					newpoints[3].x = x-0.5f; newpoints[3].y = y+0.5f; newpoints[3].z = z+0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
				// -Z
				if (levelset[x][y][z-1]>=0) {
					// boundary: create face, vertices, edges
					// vertices
					newpoints[0].x = x-0.5f; newpoints[0].y = y-0.5f; newpoints[0].z = z-0.5f;
					newpoints[1].x = x+0.5f; newpoints[1].y = y-0.5f; newpoints[1].z = z-0.5f;
					newpoints[2].x = x+0.5f; newpoints[2].y = y+0.5f; newpoints[2].z = z-0.5f;
					newpoints[3].x = x-0.5f; newpoints[3].y = y+0.5f; newpoints[3].z = z-0.5f;
					// faces
					// check for already exisiting point, otherwise new
					indexes[f+0] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[0].x && points[vp].y==newpoints[0].y && points[vp].z==newpoints[0].z) {
							indexes[f+0] = vp;
							vp = v;
						}
					}
					if (indexes[f+0]==-1) {
						points[v].x = newpoints[0].x;
						points[v].y = newpoints[0].y;
						points[v].z = newpoints[0].z;
						indexes[f+0] = v;
						v++;
					}
					indexes[f+5] = indexes[f+0];
					
					indexes[f+4] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[1].x && points[vp].y==newpoints[1].y && points[vp].z==newpoints[1].z) {
							indexes[f+4] = vp;
							vp = v;
						}
					}
					if (indexes[f+4]==-1) {
						points[v].x = newpoints[1].x;
						points[v].y = newpoints[1].y;
						points[v].z = newpoints[1].z;
						indexes[f+4] = v;
						v++;
					}
					
					indexes[f+2] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[2].x && points[vp].y==newpoints[2].y && points[vp].z==newpoints[2].z) {
							indexes[f+2] = vp;
							vp = v;
						}
					}
					if (indexes[f+2]==-1) {
						points[v].x = newpoints[2].x;
						points[v].y = newpoints[2].y;
						points[v].z = newpoints[2].z;
						indexes[f+2] = v;
						v++;
					}
					indexes[f+3] = indexes[f+2];
					
					indexes[f+1] = -1;
					for (int vp=0;vp<v;vp++) {
						if (points[vp].x==newpoints[3].x && points[vp].y==newpoints[3].y && points[vp].z==newpoints[3].z) {
							indexes[f+1] = vp;
							vp = v;
						}
					}
					if (indexes[f+1]==-1) {
						points[v].x = newpoints[3].x;
						points[v].y = newpoints[3].y;
						points[v].z = newpoints[3].z;
						indexes[f+1] = v;
						v++;
					}					
					// update sizes
					f+=6;
				}
			}
		}
		Point3f[] finalpts = new Point3f[v];
		for (int vp=0;vp<v;vp++) {
			finalpts[vp] = new Point3f();
			finalpts[vp].x = points[vp].x;
			finalpts[vp].y = points[vp].y;
			finalpts[vp].z = points[vp].z;
		}			
		return new EmbeddedSurface(finalpts, indexes);
	}
}
