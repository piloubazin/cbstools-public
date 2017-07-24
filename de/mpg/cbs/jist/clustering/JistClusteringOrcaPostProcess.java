package de.mpg.cbs.jist.clustering;

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
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
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

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class JistClusteringOrcaPostProcess extends ProcessingAlgorithm {
													
	private ParamVolume orcaClusterImage;
	private ParamVolume orcaMeanImage;
	private ParamVolume orcaStdevImage;
	private ParamVolume orcaNptImage;
	
	private ParamInteger 	nbestParam;
	private ParamFloat		coverageParam;
	
	private ParamVolume nbestImage;
	private ParamVolume coverageImage;
	private ParamVolume boundaryImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(orcaClusterImage = new ParamVolume("Orca Clusters"));
		inputParams.add(orcaMeanImage = new ParamVolume("Orca Mean"));
		inputParams.add(orcaStdevImage = new ParamVolume("Orca Stdev"));
		inputParams.add(orcaNptImage = new ParamVolume("Orca Npt"));
		
		inputParams.add(nbestParam = new ParamInteger("Number of clusters to keep", 0, 100000, 20));
		inputParams.add(coverageParam = new ParamFloat("Minimum coverage percentage", 0.0f, 100.0f, 95.0f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Clustering.devel");
		inputParams.setLabel("Orca Post-processing");
		inputParams.setName("OrcaPostProcessing");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Extracts meaningful quantities from raw ORCA clustering.");
		
		info.setVersion("3.1.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(nbestImage = new ParamVolume("N-best Labels",VoxelType.INT));
		outputParams.add(coverageImage = new ParamVolume("Coverage Labels",VoxelType.INT));
		//outputParams.add(boundaryImage = new ParamVolume("boundary Image",VoxelType.FLOAT));
		
		outputParams.setName("orca post process");
		outputParams.setLabel("orca-post-process");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		// orca input
		ImageData input = orcaClusterImage.getImageData();
		int nx = input.getRows();
		int ny = input.getCols();
		int nz = input.getSlices();
		int nxyz = nx*ny*nz;
		int[] orcaclusters = new int[nxyz];
		float[] orcamean = new float[nxyz];
		float[] orcastdev = new float[nxyz];
		float[] orcalognpt = new float[nxyz];
		int[][][] intbuffer = (new ImageDataInt(orcaClusterImage.getImageData())).toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			orcaclusters[xyz] = intbuffer[x][y][z];
		}
		intbuffer = null;
		float[][][] buffer = (new ImageDataFloat(orcaMeanImage.getImageData())).toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			orcamean[xyz] = buffer[x][y][z];
		}
		buffer = null;	
		buffer = (new ImageDataFloat(orcaStdevImage.getImageData())).toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			orcastdev[xyz] = buffer[x][y][z];
		}
		buffer = null;	
		buffer = (new ImageDataFloat(orcaNptImage.getImageData())).toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			orcalognpt[xyz] = buffer[x][y][z];
		}
		buffer = null;	
		
		// main algorithm: make ORCA cluster lists
		BasicInfo.displayMessage("extracting ORCA labels...\n");
		int[] clusterlist = ObjectLabeling.listLabels(orcaclusters, nx, ny, nz);
		float[] clusteravg = new float[clusterlist.length];
		float[] clusterstd = new float[clusterlist.length];
		int[] clustersize= new int[clusterlist.length];
		int ncluster = clusterlist.length;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (orcaclusters[xyz]>0) {
				clusteravg[orcaclusters[xyz]-1] = orcamean[xyz];
				clusterstd[orcaclusters[xyz]-1] = orcastdev[xyz];
				clustersize[orcaclusters[xyz]-1] = (int)FastMath.exp(orcalognpt[xyz]-1);
			}
		}

		// get the N best / max. coverage
		BasicInfo.displayMessage("finding the best clusters...\n");
		float coverage = coverageParam.getValue().floatValue()/100.0f;
		int nbest = nbestParam.getValue().intValue();
		float ratio=0.0f;
		int fullsize = 0;
		for (int n=0;n<ncluster;n++) fullsize += clustersize[n];
		coverage *= fullsize;
			
		int nkept=0;
		int nratio=0;
		int[] keptlabel = new int[ncluster];
		while (ratio<coverage || nkept<nbest) {
			BasicInfo.displayMessage(".");
			int id=0;
			for (int n=1;n<ncluster;n++) {
				if (clustersize[n]>clustersize[id]) id = n;
			}
			nkept++;
			if (ratio<coverage) nratio++;
			ratio += clustersize[id];
			
			keptlabel[id] = nkept;
			
			clustersize[id] *= -1;
		}
		BasicInfo.displayMessage("\n re-mapping the clusters");
		
		// warp into subject space, and compare with orca clusters
		int[][][] nbestlabel = new int[nx][ny][nz];
		int[][][] coverlabel = new int[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (orcaclusters[xyz]>0) {
				if (keptlabel[orcaclusters[xyz]-1]>0 && keptlabel[orcaclusters[xyz]-1]<=nbest)
					nbestlabel[x][y][z] = keptlabel[orcaclusters[xyz]-1];
				if (keptlabel[orcaclusters[xyz]-1]>0 && keptlabel[orcaclusters[xyz]-1]<=nratio)
					coverlabel[x][y][z] = keptlabel[orcaclusters[xyz]-1];
			}
		}
		
		// outputs
		BasicInfo.displayMessage("generating outputs...\n");
			
		String imgname = input.getName();
		
		ImageDataInt segData = new ImageDataInt(nbestlabel);	
		nbestlabel = null;
		segData.setHeader(input.getHeader());
		segData.setName(imgname+"_nbest");
		nbestImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("n best labels");
		
		segData = new ImageDataInt(coverlabel);	
		nbestlabel = null;
		segData.setHeader(input.getHeader());
		segData.setName(imgname+"cover");
		coverageImage.setValue(segData);
		segData = null;
		BasicInfo.displayMessage("coverage labels");
		
		return;
	}

}
