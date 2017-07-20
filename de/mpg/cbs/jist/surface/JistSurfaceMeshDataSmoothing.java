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


public class JistSurfaceMeshDataSmoothing  extends ProcessingAlgorithm{
	ParamSurface 	inputSurface;
	ParamFloat		fwhmParam;
	
	ParamSurface 	outputSurface;
	    
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(inputSurface = new ParamSurface("Input Surface"));
		inputParams.add(fwhmParam = new ParamFloat("FWHM", 0.0f, 20.0f, 5.0f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		inputParams.setName("SmoothSurfaceMeshData");
		inputParams.setLabel("Smooth Surface Mesh Data");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Smooths the data mapped onto a surface mesh using an iterative smoothing algorithm that takes the average of a vertex and its neighbours. Hagler et al. NeuroImage 2006");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outputSurface=new ParamSurface("Output Surface"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surf = inputSurface.getSurface();
		EmbeddedSurface smoothsurf = inputSurface.getSurface();

		float fwhm = fwhmParam.getValue().floatValue();
		int nbIterations = FastMath.round((fwhm/1.25f)*(fwhm/1.25f));
		
		BasicInfo.displayMessage("Number of iterations = "+nbIterations+"\n");
		
		surf.buildAllTables();
		int[][] neighbors = surf.getNeighborVertexVertexTable();
		
		for (int i=0; i<nbIterations; i++) {
			BasicInfo.displayMessage("iter "+(i+1)+"\n");
			double[][] sdata = surf.getVertexData();
			
			for(int v=0; v<surf.getVertexCount(); v++) {
			 
				int N = neighbors[v].length;
				for (int d=0; d<sdata[0].length; d++) {
					for (int n=0; n<N; n++) {			
						sdata[v][d] += surf.getVertexData(neighbors[v][n])[d];
					}
					sdata[v][d] /= (N+1);
				}
								
				smoothsurf.setVertexData(sdata);
			}
			surf = smoothsurf.clone();
		}
		
		smoothsurf.setName(surf.getName()+"_smoothdata");
		outputSurface.setValue(smoothsurf);
		
	}


}
