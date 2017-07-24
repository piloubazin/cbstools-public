package de.mpg.cbs.jist.statistics;

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


public class JistStatisticsSurfaceData  extends ProcessingAlgorithm{
	
	ParamSurface	groupDataSurface;
	ParamSurface 	statisticsSurface;     
	
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(groupDataSurface = new ParamSurface("Surface with Group Data"));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Statistics");
		inputParams.setName("SurfaceGroupStatistics");
		inputParams.setLabel("Surface Group Statistics");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine Lucas Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Performs basic statistics (mean, stdev) on aligned group surface data embedded on a single mesh.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(statisticsSurface = new ParamSurface("Output Surface"));
		
		outputParams.setName("surface statistics");
		outputParams.setLabel("surface statistics");
		
	}
	
	@Override
	protected void execute(CalculationMonitor monitor) {
		
		EmbeddedSurface surfdata = groupDataSurface.getSurface();
		EmbeddedSurface surfstats = groupDataSurface.getSurface().clone();
			
		int N = surfdata.getVertexData(0).length;
		
		BasicInfo.displayMessage("N = "+N+"\n");
		
		double[][] data = surfdata.getVertexData();
		double[][] stats = new double[data.length][2];
		
		for(int i=0; i<data.length; i++) {
			stats[i][0] = 0.0f;
			for(int n=0; n<N; n++) {
				stats[i][0] += data[i][n];
			}
			stats[i][0] /= N;
			stats[i][1] = 0.0f;
			for(int n=0; n<N; n++) {
				stats[i][1] += (data[i][n]-stats[i][0])*(data[i][n]-stats[i][0]);
			}
			stats[i][1] = FastMath.sqrt(stats[i][1]/N);
		}
		
		surfstats.setVertexData(stats);
		
		surfstats.setName(surfdata.getName()+"_stats");
		statisticsSurface.setValue(surfstats);
		
	}


}
