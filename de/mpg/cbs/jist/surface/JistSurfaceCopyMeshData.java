package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import org.apache.commons.math3.util.FastMath;

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


public class JistSurfaceCopyMeshData  extends ProcessingAlgorithm{
	ParamSurface 	surfaceWithData;
	ParamSurface 	surfaceWithoutData;
	
	ParamInteger 	dataIndexParam;
	ParamBoolean 	allDataParam;
	ParamBoolean 	appendDataParam;
	
	ParamSurface 	outputSurface;
	
	private static final float	UNKNOWN = -1.0f;
       
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(surfaceWithData=new ParamSurface("Surface Mesh with Data"));
		inputParams.add(surfaceWithoutData=new ParamSurface("Surface Mesh without Data"));
		inputParams.add(dataIndexParam = new ParamInteger("Data Index", 0, 100, 0));
		inputParams.add(allDataParam=new ParamBoolean("Include all data", false));
		inputParams.add(appendDataParam=new ParamBoolean("Append data", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setName("CopyMeshData");
		inputParams.setLabel("Copy Mesh Data");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Copies data embedded in one mesh onto another mesh with matching vertices.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputSurface=new ParamSurface("Output Surface"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surfwdata = surfaceWithData.getSurface();
		EmbeddedSurface newsurf = surfaceWithoutData.getSurface().clone();
		
		int idx = dataIndexParam.getValue().intValue();
			
		int N;
		if (appendDataParam.getValue().booleanValue()) {
			N = newsurf.getVertexData(0).length;
		} else {
			N = 0;
		}
		
		int M;
		if (allDataParam.getValue().booleanValue()) {
			M = surfwdata.getVertexData(0).length;
		} else {
			M = 1;
		}
		
		double[][] data = new double[newsurf.getVertexData().length][N+M];
				
		if (appendDataParam.getValue().booleanValue()) {
			for(int i=0; i<newsurf.getVertexData().length; i++) {
				for(int n=0; n<N; n++) {
					data[i][n] = newsurf.getVertexData(i)[n];
				}
			}
		}
		
		if (allDataParam.getValue().booleanValue()) {
			for(int i=0; i<newsurf.getVertexData().length; i++) {
				for(int m=0; m<M; m++) {
					data[i][N+m] = surfwdata.getVertexData(i)[m];
				}
			}
		} else {
			for(int i=0; i<newsurf.getVertexData().length; i++) {
				data[i][N] = surfwdata.getVertexData(i)[idx];
			}
		}
		
		newsurf.setVertexData(data);
		
		newsurf.setName(newsurf.getName()+"_data");
		outputSurface.setValue(newsurf);
		
	}


}
