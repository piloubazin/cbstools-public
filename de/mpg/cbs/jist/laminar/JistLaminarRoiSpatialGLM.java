package de.mpg.cbs.jist.laminar;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
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

import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.apache.commons.math3.util.FastMath;
import Jama.Matrix;
import java.util.*;

/*
 * @author Pierre-Louis Bazin
 */
public class JistLaminarRoiSpatialGLM extends ProcessingAlgorithm {

	// jist containers
	private ParamCollection imageParams;
	private ParamCollection mainParams;
	
	private ParamVolume layersImage;
	private ParamVolume dataImage;
	private ParamVolume roiImage;
	private ParamOption 	glmtypeParam;
	private static final String[] glmTypes = {"local_pv_full", 
											  "local_pv_fast",
											  "basic_bold_full",
											  "basic_bold_fast"};
	private ParamFloat boldratioParam;	
											  
	private ParamFile	textFile;
	private ParamVolume glmImage;
	private ParamVolume resImage;
	private ParamVolume pvolImage;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final float HASQ3 = (float)(FastMath.sqrt(3.0)/2.0);

	
	protected void createInputParameters(ParamCollection inputParams) {
		
		//imageParams=new ParamCollection("Images");
		inputParams.add(dataImage = new ParamVolume("Data Image",null,-1,-1,-1,-1));
		inputParams.add(layersImage = new ParamVolume("Profile Surface Image",null,-1,-1,-1,-1));
		inputParams.add(roiImage = new ParamVolume("Region(s) of interest",null,-1,-1,-1,-1));
		roiImage.setDescription("a 3D or 4D image of regions of interest");
		inputParams.add(glmtypeParam = new ParamOption("GLM estimation method", glmTypes));
		glmtypeParam.setDescription("local pv: only taking into account the local layer contribution; basic bold: assumes downstream mixing; full: high precision pve; fast: approximate pve");
		inputParams.add(boldratioParam = new ParamFloat("Downstream mixing ratio", 0.0f, 1.0f, 0.2f));
		boldratioParam.setDescription("ratio of mixing in the basic bold model");
		
		//inputParams.add(imageParams);
			
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Laminar Analysis.devel");
		inputParams.setLabel("ROI Spatial GLM");
		inputParams.setName("ROISpatialGLM");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Spinoza Centre | Netherlands Institute for Neuroscience | Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Sample some data image across layer surfaces assuming a signal mixture (spatial GLM) in a given set of ROIs.");
		
		info.setVersion("3.1.1");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(textFile = new ParamFile("ROI Estimates"));

		outputParams.add(glmImage = new ParamVolume("GLM-based data",null,-1,-1,-1,-1));
		outputParams.add(resImage = new ParamVolume("GLM residuals",null,-1,-1,-1,-1));
		outputParams.add(pvolImage = new ParamVolume("partial volume estimates",null,-1,-1,-1,-1));
		
		outputParams.setName("roi GLM images");
		outputParams.setLabel("roi GLM images");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		String intensname = Interface.getName(dataImage);
		ImageHeader header = Interface.getHeader(dataImage);
		
		int[] dims = Interface.getDimensions(layersImage);
		float[] res = Interface.getResolutions(layersImage);
		int nlayers = Interface.getComponents(layersImage)-1;
		int nxyz = dims[X]*dims[Y]*dims[Z];
		int nx = dims[X]; int ny = dims[Y]; int nz = dims[Z]; 
		
		// at least two layer surfaces (gwb and cgb)
		float[] layers = Interface.getFloatImage4D(layersImage);
		
		float[] data = null;
		int nt = Interface.getComponents(dataImage);
		if (Interface.isImage4D(dataImage)) {
			data = Interface.getFloatImage4D(dataImage);
		} else {
			data = Interface.getFloatImage3D(dataImage);
		}
		
		// create a set of ROI masks for all the regions and also outside of the area where layer 1 is > 0 and layer N is < 0 (including boundaries: sqrt(3)/2 distance)
		boolean[] ctxmask = new boolean[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) {
			ctxmask[xyz] = (layers[xyz]>=-HASQ3 && layers[xyz+nlayers*nxyz]<=HASQ3);
		}
		
		String roiname = Interface.getName(roiImage);
		boolean[][] rois = null;
		int[] roilist = null;
		int nroi = 0;
		if (Interface.isImage4D(roiImage)) {
			int[] roilabels = Interface.getIntegerImage4D(roiImage);
			nroi = Interface.getComponents(roiImage);
			rois = new boolean[nroi][nxyz];
			for (int r=0;r<nroi;r++) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					rois[r][xyz] = (ctxmask[xyz] && roilabels[xyz+nxyz*r]>0);
				}
			}
			roilist = new int[nroi];
			for (int r=0;r<nroi;r++) roilist[r] = r+1;
			roilabels = null;
		} else {
			int[] roilabels = Interface.getIntegerImage3D(roiImage);
			roilist = ObjectLabeling.listLabels(roilabels, nx, ny, nz);
			for (int r=0;r<roilist.length;r++) if (roilist[r]>0) nroi++;
			System.out.println("number of ROIs found: "+nroi+" ("+roilist.length+")\n");
			rois = new boolean[nroi][nxyz];
			int rp=0;
			for (int r=0;r<roilist.length;r++) if (roilist[r]>0) {
				for (int xyz=0;xyz<nxyz;xyz++) {
					rois[rp][xyz] = (ctxmask[xyz] && roilabels[xyz]==roilist[r]);
				}
				rp++;
			}
			roilabels = null;
		}
		// speed-up: get rid of any voxel outside of the ROIs
		for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
			ctxmask[xyz] = false;
			for (int r=0;r<nroi && !ctxmask[xyz];r++) {
				if (rois[r][xyz]) ctxmask[xyz] = true;
			}
		}
		
		// main algorithm
		boolean fast = false;
		if (glmtypeParam.getValue().contains("fast")) fast = true;
		boolean bold = false;
		if (glmtypeParam.getValue().contains("basic_bold")) bold = true;
		float boldratio = boldratioParam.getValue().floatValue();
		
		// 1. define partial voume for each layer, each voxel
		float[][] pvol = new float[nlayers+1][nxyz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z = 0; z<nz; z++) {
			int xyz = x + nx*y + nx*ny*z;
			if (ctxmask[xyz]) {
				for (int l=0;l<=nlayers;l++) {
					if (fast) pvol[l][xyz] = fastApproxPartialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
					else pvol[l][xyz] = partialVolumeFromSurface(x, y, z, layers, l*nxyz, nx, ny, nz);
				}
			}
		}
		
		// 2. build and invert the GLM for each ROI
		float[][][] estimates = new float[nroi][nlayers][nt];
		float[] residuals = new float[nroi];
		for (int r=0;r<nroi;r++) {
			int nsample = 0;
			for (int xyz=0;xyz<nxyz;xyz++) if (rois[r][xyz]) nsample++;
			System.out.println("ROI "+(r+1)+": "+nsample+" samples\n");
			
			// GLM*VAL = SAMPLES
			// nsamples x layers * nlayers x nt = nsamples x nt
			double[][] glm = new double[nsample][nlayers];
			double[][] samples = new double[nsample][nt];
			
			int n=0;
			for (int xyz=0;xyz<nxyz;xyz++) if (rois[r][xyz]) {
				// different version for bold integration
				for (int l=0;l<nlayers;l++) glm[n][l] = Numerics.max(pvol[l+1][xyz]-pvol[l][xyz],0.0f);
				for (int t=0;t<nt;t++) samples[n][t] = data[xyz+nxyz*t];
				n++;
			}
			// invert the linear model
			Matrix mtx = new Matrix(glm);
			Matrix smp = new Matrix(samples);
			Matrix val = mtx.solve(smp);
			for (int l=0;l<nlayers;l++) for (int t=0;t<nt;t++) {
				estimates[r][l][t] = (float)val.get(l,t);
			}
			Matrix err = mtx.times(val).minus(smp);
			residuals[r] = (float)err.normInf();
		}
		
		// create an data volume from the estimated GLM and compute the distance (optional, but nice for checking)
		float[] glmdata = new float[nxyz*nt];
		float[] count = new float[nxyz];
		float[] resmap = new float[nxyz];
		float[] pvolmap = new float[nxyz*nlayers];
		for (int xyz=0;xyz<nxyz;xyz++) if (ctxmask[xyz]) {
			for (int r=0;r<nroi;r++) if (rois[r][xyz]) {
				for (int l=0;l<nlayers;l++) {
					pvolmap[xyz+nxyz*l] =  Numerics.max(pvol[l+1][xyz]-pvol[l][xyz], 0.0f);
					for (int t=0;t<nt;t++) {
						glmdata[xyz+t*nxyz] += pvolmap[xyz+nxyz*l]*estimates[r][l][t];
					}
					count[xyz] += pvolmap[xyz+nxyz*l];
				}
			}
			if (count[xyz]>0) {
				for (int t=0;t<nt;t++) {
					glmdata[xyz+t*nxyz] /= count[xyz];
					resmap[xyz] += Numerics.square(data[xyz+t*nxyz] - glmdata[xyz+t*nxyz]);
				}
				resmap[xyz] = (float)FastMath.sqrt(resmap[xyz]/nt);
			}
		}
				
		// output
		String filename = intensname+"_roiglm_"+roiname+".txt";
		textFile.setValue(writeROIestimateFile(filename, estimates, residuals, roilist, nroi, nlayers, nt));
		
		Interface.setFloatImage4D(glmdata, dims, nt, glmImage, intensname+"_roiglm_approx", header);
		
		Interface.setFloatImage3D(resmap, dims, resImage, intensname+"_roiglm_err", header);
		
		Interface.setFloatImage4D(pvolmap, dims, nlayers, pvolImage, intensname+"_roiglm_pvol", header);
	}
	
	final File writeROIestimateFile(String filename, float[][][] estimates, float[] residuals, int[] roilist, int nroi, int nlayers, int nt) {
		File resultFile = null;
		File destdir = null;
		try{
			destdir = new File(this.getOutputDirectory().getCanonicalFile()+File.separator+this.getAlgorithmName());
			if(!destdir.isDirectory()){
				(new File(destdir.getCanonicalPath())).mkdir();
			}
			resultFile = new File(destdir+File.separator+filename);
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		FileWriter fw = null;
		try {
			fw = new FileWriter(resultFile);
		} catch (IOException e) {
			e.printStackTrace();
		}
		PrintWriter pw = new PrintWriter( fw );
		pw.write("GLM-inverted data per layer per roi : residuals per roi");
		pw.printf("\n");
		// open the file to write
		try {
			for (int r=0;r<nroi;r++) {
				pw.printf("ROI "+roilist[r]+"\n");
				for (int l=0;l<nlayers;l++) {
					pw.printf("layer "+(l+1)+", ");
					for (int t=0;t<nt;t++) pw.printf(" %.20f,", estimates[r][l][t]);
					pw.printf("\n");
				}
				pw.printf("residuals, "+residuals[r]+"\n \n");
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}	

		pw.close();
		try {
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		fw = null;
		pw = null;

		 return resultFile;
	}

	

	// brute force
	void findProfileNeighborhood(BitSet sampled, int x, int y, int z, int delta, float[][] layers, CorticalProfile profile, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
		sampled.clear();
		
		// compute profiles around a certain size
		for (int dx=-delta;dx<=delta;dx++) for (int dy=-delta;dy<=delta;dy++) for (int dz=-delta;dz<=delta;dz++) {
			profile.computeTrajectory(layers, x+dx, y+dy, z+dz);
		
			// add the underlying voxels to the estimation matrix
			for (int l=0;l<=nlayers;l++) {
				int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
				if (ctxmask[xyzs]) sampled.set(xyzs);
			}
		}
	}
	
	// should be faster and more precise
	void findFastMarchingProfileNeighborhood(BitSet sampled, float[] dist, int x, int y, int z, float delta, float[][] layers, CorticalProfile profile, boolean[] ctxmask, int nx, int ny, int nz, int nlayers) {
		sampled.clear();
		
		// compute profile centered on point of interest
		profile.computeTrajectory(layers, x, y, z);
		
		// add the underlying voxels to the estimation matrix
		for (int l=0;l<=nlayers;l++) {
			int xyzs = Numerics.round(profile.getPt(l)[X]) + nx*Numerics.round(profile.getPt(l)[Y]) + nx*ny*Numerics.round(profile.getPt(l)[Z]);
			if (ctxmask[xyzs]) sampled.set(xyzs);
		}
		
		// find enighbors through fast marching
		BinaryHeap2D heap = new BinaryHeap2D(Numerics.ceil(2*delta)*sampled.size(), BinaryHeap2D.MINTREE);
		int idx = 0;
		for (int n=0;n<sampled.size();n++) {
			// get the next non-zero value
			idx = sampled.nextSetBit(idx);
			heap.addValue(0.0f,idx,(byte)1);
		}
		// reset the sampled values: the original ones will be re-added from the heap
		sampled.clear();
		// main loop
		while (heap.isNotEmpty()) {
        	// extract point with minimum distance
        	float curdist = heap.getFirst();
        	int xyz = heap.getFirstId();
        	heap.removeFirst();
        	
        	// for points selected multiple times
        	if (sampled.get(xyz)) continue;
        	
        	// add to sampled values
        	sampled.set(xyz);
        	dist[xyz] = curdist;
        	
        	for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++) {
				int ngb = xyz+dx+nx*dy+nx*ny*dz;
				if (!sampled.get(ngb) && ctxmask[ngb]) {
					// add to the tree if close enough from first profile
					float newdist = curdist+(float)FastMath.sqrt(dx*dx+dy*dy+dz*dz);
					if (newdist<=delta) heap.addValue(newdist, ngb, (byte)1);
				}
			}
		}
	}

	float partialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int offset, int nx, int ny, int nz) {
		int xyz0 = x0+nx*y0+nx*ny*z0;
		double phi0 = levelset[xyz0+offset];
		
		if (phi0<-HASQ3) return 1.0f;
		else if (phi0>HASQ3) return 0.0f;
		
		// use direction, distance to create a continuous plane model V(X-X0)+phi(X0) = d(X,P)
		double vx = 0.5*(levelset[xyz0+1+offset]-levelset[xyz0-1+offset]);
		double vy = 0.5*(levelset[xyz0+nx+offset]-levelset[xyz0-nx+offset]);
		double vz = 0.5*(levelset[xyz0+nx*ny+offset]-levelset[xyz0-nx*ny+offset]);
		double norm = FastMath.sqrt(vx*vx+vy*vy+vz*vz);
		
		// integrate the voxel cube over that model
		
		// first shot: numerically (other options: Heaviside integration, case-by-case geometric solution?)
		double volume = 0.0;
		for (double x=-0.45; x<=0.45; x+=0.1) for (double y=-0.45; y<=0.45; y+=0.1) for (double z=-0.45; z<=0.45; z+=0.1) {
			if (vx*x+vy*y+vz*z+norm*phi0<=0) volume+=0.001;
		}
		
		// second option: with explicit definition of the intersection (~marching cubes!)
		// requires intersecting all 12 edges from the voxel cube, then defining either case-by-case formulas or 
		// deriving the formula for an integral of the Heaviside function
		
		return (float)volume;
	}

	float fastApproxPartialVolumeFromSurface(int x0, int y0, int z0, float[] levelset, int offset, int nx, int ny, int nz) {
		return Numerics.bounded(0.5f-levelset[x0+nx*y0+nx*ny*z0+offset],0.0f,1.0f);
	}	
}
