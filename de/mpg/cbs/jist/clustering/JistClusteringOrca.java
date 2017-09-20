package de.mpg.cbs.jist.clustering;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamNumberCollection;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class JistClusteringOrca extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume 	intensityImage;
	private ParamVolume 	maskImage;
	private ParamVolume 	stdevImage;
	
	private ParamFloat		globaldevParam;
	private ParamFloat		thresholdParam;
	private ParamInteger		maxscaleParam;
	private ParamInteger		minscaleParam;
	private ParamFloat		sizeParam;
	private ParamFloat		recomputeParam;
	private ParamBoolean	skipParam;
	private ParamFloat		devfactorParam;
	private ParamBoolean	propagParam;
	private ParamOption		metricParam;
	private ParamOption		varupdateParam;
	private ParamInteger		clustersizeParam;
	private ParamInteger		connectParam;
	private ParamBoolean	cleanclusterParam;
	
	private		static final String[]	metricTypes = {"joint Jensen-Shannon","max Jensen-Shannon","exact Jensen-Shannon","table Jensen-Shannon"};
	//private		String		metricType = "max Jensen-Shannon";
	private		String		metricType = "table Jensen-Shannon";
	
	private		static final String[]	varupdateTypes = {"Maximum","Linear","Zero","Constant"};
	//private		String		varupdateType = "Maximum";
	private		String		varupdateType = "Zero";
	
	private ParamVolume clusterImage;
	private ParamVolume scaleImage;
	private ParamVolume meanImage;
	private ParamVolume varImage;
	private ParamVolume probaImage;
	private ParamVolume octreeImage;
	private ParamVolume mergingImage;
	private ParamVolume growingImage;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensityImage = new ParamVolume("Intensity Image (4D)",null,-1,-1,-1,-1));
		inputParams.add(maskImage = new ParamVolume("Mask Image (opt)"));
		maskImage.setMandatory(false);
		inputParams.add(stdevImage = new ParamVolume("Standard Deviation Image (opt)"));
		stdevImage.setMandatory(false);
		
		inputParams.add(thresholdParam=new ParamFloat("Threshold parameter", 0.0f, 1e9f, 0.5f));
		inputParams.add(globaldevParam=new ParamFloat("Noise standard deviation", 0.0f, 1e9f, 100.0f));
		inputParams.add(maxscaleParam=new ParamInteger("Max octree scale", 0, 20, 10));
		inputParams.add(minscaleParam=new ParamInteger("Min clustering scale", 0, 20, 0));
		inputParams.add(sizeParam=new ParamFloat("Size ratio for relation pruning", 0f, 1f, 0.125f));
		inputParams.add(recomputeParam=new ParamFloat("Score ratio for recomputing", 0f, 1f, 0.9f));
		inputParams.add(skipParam=new ParamBoolean("Skip clustering", false));
		inputParams.add(devfactorParam=new ParamFloat("Standard deviation factor", 0f, 1e9f, 1.0f));
		inputParams.add(propagParam=new ParamBoolean("Propagate octree labels", true));
		inputParams.add(metricParam = new ParamOption("Clustering metric", metricTypes));
		metricParam.setValue(metricType);
		inputParams.add(varupdateParam = new ParamOption("Variance update", varupdateTypes));
		varupdateParam.setValue(varupdateType);
		inputParams.add(clustersizeParam=new ParamInteger("Minimum cluster size to keep", 0, 10000, 50));
		inputParams.add(connectParam=new ParamInteger("Connectivity", 6, 26, 6));
		inputParams.add(cleanclusterParam=new ParamBoolean("Automatic cluster threshold", false));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Clustering.devel");
		inputParams.setLabel("ORCA Clustering");
		inputParams.setName("OrcaClustering");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Clusters 4D data with octree region competition aggregative clustering.");
		
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(clusterImage = new ParamVolume("Clusters", VoxelType.INT));
		outputParams.add(scaleImage = new ParamVolume("Scale", VoxelType.FLOAT));
		outputParams.add(meanImage = new ParamVolume("Mean",null,-1,-1,-1,-1));
		outputParams.add(varImage = new ParamVolume("Stdev",null,-1,-1,-1,-1));
		outputParams.add(probaImage = new ParamVolume("Proba",null,-1,-1,-1,-1));
		outputParams.add(octreeImage = new ParamVolume("Octree", VoxelType.UBYTE));
		outputParams.add(mergingImage = new ParamVolume("Merging", null,-1,-1,-1,-1));
		outputParams.add(growingImage = new ParamVolume("Growing", null,-1,-1,-1,-1));
		
		outputParams.setName("orca clustering");
		outputParams.setLabel("orca clustering");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		System.out.println("ORCA clustering process");
		System.out.flush();
		
		// import the image data into 1D arrays
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		
		int nx = intensImg.getRows();
		int ny = intensImg.getCols();
		int nz = intensImg.getSlices();
		int nt = intensImg.getComponents();
		int nxyz = nx*ny*nz;
		float rx = intensImg.getHeader().getDimResolutions()[0];
		float ry = intensImg.getHeader().getDimResolutions()[1];
		float rz = intensImg.getHeader().getDimResolutions()[2];
		String imgname = intensImg.getName();
		
		float[][] intensity = new float[nxyz][nt];
		if (nt>1) {
			float[][][][] buffer4 = intensImg.toArray4d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[xyz][t] = buffer4[x][y][z][t];
			}
			buffer4 = null;
		} else {
			float[][][] buffer3 = intensImg.toArray3d();
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				intensity[xyz][0] = buffer3[x][y][z];
			}
			buffer3 = null;			
		}
		intensImg = null;
		
		// mask for regions with no data
		System.out.print("import mask...");
		System.out.flush();
		boolean[] mask = new boolean[nxyz];
		if (maskImage.getImageData()!=null) {
			ImageDataUByte	maskImg = new ImageDataUByte(maskImage.getImageData());
			byte[][][] bytebuffer = maskImg.toArray3d();
			
			int nvol=0;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				//boolean same=true;
				//for (int l=0;l<nlayers && same;l++) if (intensity[l+1][xyz]!=intensity[l][xyz]) same = false;
				//ctxmask[xyz] = (!same);
				mask[xyz] = (bytebuffer[x][y][z]>0);
				if (mask[xyz]) nvol++;
			}
			bytebuffer = null;
			maskImg = null;
			System.out.println("volume: "+nvol);
		} 
		// remove the image boundary by default (to avoid neighborhood checks)
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
			int z=0;
			int xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;		
			z = nz-1;
			xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;
		}
		for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int x=0;
			int xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;
			x=nx-1;
			xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;
		}
		for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
			int y=0;
			int xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;
			y=ny-1;
			xyz = x+nx*y+nx*ny*z;
			mask[xyz] = false;
		}
		
		float[][] stdev = null;
		if (stdevImage.getImageData()!=null) {
			stdev = new float[nxyz][nt];
			ImageDataFloat	stdevImg = new ImageDataFloat(stdevImage.getImageData());
			if (nt>1) {
				float[][][][] buffer4 = stdevImg.toArray4d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int t=0;t<nt;t++) {
					int xyz = x+nx*y+nx*ny*z;
					stdev[xyz][t] = buffer4[x][y][z][t];
				}
				buffer4 = null;
			} else {
				float[][][] buffer3 = stdevImg.toArray3d();
				for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
					int xyz = x+nx*y+nx*ny*z;
					stdev[xyz][0] = buffer3[x][y][z];
				}
				buffer3 = null;			
			}
			stdevImg = null;
		}

		// main algorithm
		float[] globaldev = new float[nt];
		for (int t=0;t<nt;t++) {
			globaldev[t] = globaldevParam.getValue().floatValue();
		}
		float threshold = thresholdParam.getValue().floatValue();
		float sizeratio = sizeParam.getValue().floatValue();
		float recomputeratio = recomputeParam.getValue().floatValue();
		float devfactor = devfactorParam.getValue().floatValue();
		
		System.out.println("parameters: th = "+threshold+", size ratio = "+sizeratio+", recompute ratio = "+recomputeratio+", noise scaling = "+devfactor);
		System.out.print("noise estimate: "); for (int t=0;t<nt;t++) System.out.print(globaldev[t]+" "); System.out.print("\n");
		long start_time = System.currentTimeMillis();

		System.out.println("init octree based tree");
		System.out.flush();
		
		String metric = metricParam.getValue();
		String varupdate = varupdateParam.getValue();
		int clustersize = clustersizeParam.getValue().intValue();
		int clusterthreshold = clustersize;
		int nclusters = 0;
				
		int maxlevel = maxscaleParam.getValue().intValue();
		
		int connect = connectParam.getValue().intValue();
		
		OctreeMultiClusterSimplification octree = new OctreeMultiClusterSimplification(intensity, mask, stdev, nx, ny, nz, nt, globaldev, devfactor, threshold, metric, varupdate, clustersize, connect);
		
		System.out.println("prune");
		System.out.flush();
			
		maxlevel = octree.mergeToJensenShannonDistance(maxlevel);
		
		// get the shape of the octree
		byte[][][] octimg = octree.exportOctreeStructure(maxlevel);
		ImageDataUByte octData = new ImageDataUByte(octimg);		
		octData.setHeader(intensityImage.getImageData().getHeader());
		octData.setName(imgname+"_octree");
		octreeImage.setValue(octData);
		octData = null;
		octimg = null;
		
		int startlevel = Numerics.min(maxlevel, maxscaleParam.getValue().intValue());
		int stoplevel = Numerics.max(0, minscaleParam.getValue().intValue());
		
		int[][][][] merging = new int[startlevel-stoplevel+1][][][];
		int[][][][] growing = new int[startlevel-stoplevel][][][];
		
		if (!skipParam.getValue().booleanValue()) {
			float[] var = new float[nt];
			for (int t=0;t<nt;t++) var[t] = devfactor*devfactor*globaldev[t]*globaldev[t];
		
			AggregativeOctreeMultiClustering clustering = new AggregativeOctreeMultiClustering(octree, var, threshold, metric, varupdate, clustersize, connect);
		
			// debug
			int[] current;
		
			for (int level=startlevel;level>=stoplevel;level--) {
				System.out.println("create aggregative clustering for level "+level);
				System.out.flush();
			
				if (level==startlevel) {
					clustering.initFirstClusters(level);
				} else {
					clustering.initNextClusters(level, nclusters, propagParam.getValue().booleanValue(), recomputeratio);
					growing[level] = octree.exportLabelsFrom(level);
				}
				nclusters = clustering.getInitClusterNumber();
				
				System.out.println("("+nclusters+" initial clusters)");
				System.out.flush();
				
				System.out.println("compute aggregative clusters");
				System.out.flush();
		
				nclusters = clustering.hierarchicalClustering(0, (recomputeratio>0), recomputeratio, sizeratio);
			
				System.out.println(" best = "+nclusters);
			
				if (level==stoplevel) {
					if (cleanclusterParam.getValue().booleanValue()) clusterthreshold = clustering.smallClusterThreshold(nclusters);
					else clusterthreshold = clustersize;
				}
				// this operation is destructive, make sure to do all the clustering analysis before this
				nclusters = clustering.exportClusteringResults(octree.getLbl(level), octree.getSum(level), octree.getSqd(level), octree.getNpt(level), nclusters);
				System.out.println(" final = "+nclusters);
				
				// output the result in a 4D array
				merging[level-stoplevel] = octree.exportLabelsFrom(level);
			}
		} else {
			// only the propagation of labels / region competition	
			for (int level=startlevel;level>=stoplevel;level--) {
				System.out.println("create octree clustering for level "+level);
				System.out.flush();
			
				if (level==startlevel) {
					nclusters = octree.generateLabelsAt(level);
				} else {
					octree.generateNextLevelMaps(level);

					// labeling: import from previous scale and add new scale elements
					if (propagParam.getValue().booleanValue()) 
						nclusters = octree.competitionWithJensenShannonDistance(level, nclusters, recomputeratio);
					else 
						nclusters = octree.generateNextLabelsAt(level, nclusters);
				}
				
				System.out.println("("+nclusters+" initial clusters)");
				System.out.flush();
			}
		
		}
		// clean up the smaller clusters
		//if (clusterthreshold>1) nclusters = octree.removeSmallClusters(stoplevel, clusterthreshold, nclusters); 
		//System.out.println("filtered clusters (size >= "+clusterthreshold+") = "+nclusters);
		
		System.out.println("extract output maps...");
		System.out.flush();
		
		
		ImageDataInt cl1Data = new ImageDataInt(octree.exportLblAtLevel(stoplevel));		
		cl1Data.setHeader(intensityImage.getImageData().getHeader());
		cl1Data.setName(imgname+"_clusters");
		clusterImage.setValue(cl1Data);
		cl1Data = null;
		
		ImageDataFloat cl2Data = new ImageDataFloat(octree.exportLogNptAtLevel(stoplevel));		
		cl2Data.setHeader(intensityImage.getImageData().getHeader());
		cl2Data.setName(imgname+"_npt");
		scaleImage.setValue(cl2Data);
		cl2Data = null;
		
		ImageDataFloat meanData = new ImageDataFloat(octree.exportAvgAtLevel(stoplevel));		
		meanData.setHeader(intensityImage.getImageData().getHeader());
		meanData.setName(imgname+"_mean");
		meanImage.setValue(meanData);
		meanData = null;
		
		ImageDataFloat varData = new ImageDataFloat(octree.exportStdAtLevel(stoplevel));		
		varData.setHeader(intensityImage.getImageData().getHeader());
		varData.setName(imgname+"_stdev");
		varImage.setValue(varData);
		varData = null;
		
		ImageDataFloat probaData = new ImageDataFloat(octree.exportProbaAtLevel(stoplevel, intensity));		
		probaData.setHeader(intensityImage.getImageData().getHeader());
		probaData.setName(imgname+"_proba");
		probaImage.setValue(probaData);
		probaData = null;
		
		int[][][][] tmp = new int[nx][ny][nz][startlevel-stoplevel+1];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<startlevel-stoplevel+1;l++) {
			tmp[x][y][z][l] = merging[l][x][y][z];
		}
		ImageDataInt mergingData = new ImageDataInt(tmp);		
		mergingData.setHeader(intensityImage.getImageData().getHeader());
		mergingData.setName(imgname+"_merging");
		mergingImage.setValue(mergingData);
		mergingData = null;
		
		tmp = new int[nx][ny][nz][startlevel-stoplevel+1];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) for (int l=0;l<startlevel-stoplevel;l++) {
			tmp[x][y][z][l+1] = growing[l][x][y][z];
		}
		ImageDataInt growData = new ImageDataInt(tmp);		
		growData.setHeader(intensityImage.getImageData().getHeader());
		growData.setName(imgname+"_growing");
		growingImage.setValue(growData);
		growData = null;

		System.out.println("processing time (milliseconds): " + (System.currentTimeMillis()-start_time));
	}
	

}
