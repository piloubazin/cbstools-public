package de.mpg.cbs.jist.cortex;

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
 * @author Pierre-Louis Bazin
 * @author Duygu Tosun
 */
public class JistCortexSurfaceMeshInflation extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume levelsetImage;
	
	ParamFloat sor;
	ParamFloat err;
	ParamInteger maxiters;
	ParamInteger step;
	ParamSurface inSurf;
	ParamSurface outSurf;
	ParamBoolean lorentzian;
	ParamOption		topologyParam;
	
	private static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	
	private ParamSurface inflatedSurface;
	private ParamSurface origSurface;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImage = new ParamVolume("Levelset Image"));
				
		inputParams.add(sor=new ParamFloat("SOR Parameter",0f,1f,0.75f));

		inputParams.add(err=new ParamFloat("Mean Curvature Threshold",0,100,8));
		inputParams.add(step=new ParamInteger("Step Size",0,100,50));
		inputParams.add(maxiters=new ParamInteger("Max Iterations",0,10000000,10000));
		inputParams.add(lorentzian=new ParamBoolean("Lorentzian Norm",true));
			
		inputParams.add(topologyParam = new ParamOption("Topology", topoTypes));
		topologyParam.setValue("wcs");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing");
		inputParams.setLabel("Surface Mesh Inflation");
		inputParams.setName("SurfaceMeshInflation");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Duygu Tosun", "",""));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new Citation("D. Tosun, M. E. Rettmann, X. Han, X. Tao, C. Xu, S. M. Resnick, D. Pham, and J. L. Prince, "
								+"Cortical Surface Segmentation and Mapping, "
								+"NeuroImage, vol. 23, pp. S108--S118, 2004."));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Inflates a cortical surface mesh.");
		
		info.setVersion("3.0.2");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(origSurface=new ParamSurface("Original Surface"));
		outputParams.add(inflatedSurface=new ParamSurface("Inflated Surface"));

		outputParams.setName("inflated images");
		outputParams.setLabel("inflated images");
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
		surf.scaleVertices(new float[]{rx, -ry, -rz});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		surf.translate(pt);
		surf.computeNormals();
		
		// surface inflation
		SurfaceInflate inflate=new SurfaceInflate(surf);
		monitor.observe(inflate);
		EmbeddedSurface insurf=inflate.inflate(sor.getFloat(), step.getInt(),err.getFloat(), maxiters.getInt(),lorentzian.getValue());
		
		// compute the gyri / sulci map from the intersection?? not needed!
		/* too fancy
		SurfaceToMask surf2vol=new SurfaceToMask();
		monitor.observe(surf2vol);
		ImageDataFloat maskVol=surf2vol.solve(new ImageDataFloat(levelsetImg), insurf);
		
		double[][] data = surf.getVertexData();

		for(int i=0; i<data.length; i++){
			Point3f p = surf.getVertex(i);
			
			// mapping from mesh to voxel space
			int x = Numerics.round(p.x/rx);
			int y = Numerics.round((ny-1)-p.y/ry);
			int z = Numerics.round((nz-1)-p.z/rz);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				data[i] = new double[]{maskVol.getFloat(x,y,z,0)};
			}
		}
		surf.setVertexData(data);
		insurf.setVertexData(data);
		*/
		// just use the original levelset!
		double[][] data = insurf.getVertexData();

		for(int i=0; i<data.length; i++){
			Point3f p = insurf.getVertex(i);
			
			// mapping from mesh to voxel space
			int x = Numerics.round(p.x/rx);
			int y = Numerics.round((ny-1)-p.y/ry);
			int z = Numerics.round((nz-1)-p.z/rz);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				if (levelsetImg.getFloat(x,y,z,0)>=0) data[i] = new double[]{-1};
				else data[i] = new double[]{+1};
			}
		}
		insurf.setVertexData(data);
		surf.setVertexData(data);
		
		// outputs
		origSurface.setValue(surf);
		inflatedSurface.setValue(insurf);
	}


}
