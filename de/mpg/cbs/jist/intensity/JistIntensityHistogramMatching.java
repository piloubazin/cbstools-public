package de.mpg.cbs.jist.intensity;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipav;

import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
public class JistIntensityHistogramMatching extends ProcessingAlgorithm{
	ParamVolume subParam;
	ParamVolume tmpParam;
	ParamVolume resultVolParam;
	ParamBoolean zeroParam;
	ParamBoolean projectParam;

	private static final String shortDescription = "Normalizes the intensity histogram of an image to a template histogram";
	private static final String longDescription = "";


	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(subParam=new ParamVolume("Image to normalize"));
		inputParams.add(tmpParam=new ParamVolume("Template to match"));
		inputParams.add(zeroParam=new ParamBoolean("Skip zero values?", true));
		inputParams.add(projectParam=new ParamBoolean("Output in template range?", true));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Intensity");
		inputParams.setLabel("Histogram Matching");
		inputParams.setName("HistogramMatching");

		AlgorithmInformation info = getAlgorithmInformation();
		info.setWebsite("http://www.cbs.mpg.de/");
		info.setDescription(shortDescription);
		info.setLongDescription(shortDescription + longDescription);
		info.setVersion("3.0.7");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
	}


	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(resultVolParam=new ParamVolume("Matched Image",null,-1,-1,-1,-1));
	}


	protected void execute(CalculationMonitor monitor) throws AlgorithmRuntimeException {
		ImageDataFloat	tmpImg = new ImageDataFloat(tmpParam.getImageData());
		int ntx = tmpImg.getRows();
		int nty = tmpImg.getCols();
		int ntz = tmpImg.getSlices();
		float[][][] tmp = tmpImg.toArray3d();

		ImageDataFloat	subImg = new ImageDataFloat(subParam.getImageData());
		int nsx = subImg.getRows();
		int nsy = subImg.getCols();
		int nsz = subImg.getSlices();
		float[][][] sub = subImg.toArray3d();
		
		// build a cumulative mapping
		boolean skip = zeroParam.getValue().booleanValue();
		
		int ntsample=0;
		float tmin = 1e12f;
		float tmax = -1e12f;
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) if (!skip || (skip && tmp[x][y][z]!=0)) {
			ntsample++;
			if (tmp[x][y][z]>tmax) tmax = tmp[x][y][z];
			if (tmp[x][y][z]<tmin) tmin = tmp[x][y][z];
		}
		int nssample=0;
		float smin = 1e12f;
		float smax = -1e12f;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) if (!skip || (skip && sub[x][y][z]!=0)) {
			nssample++;
			if (sub[x][y][z]>smax) smax = sub[x][y][z];
			if (sub[x][y][z]<smin) smin = sub[x][y][z];
		}
		
		// build cumulative histogram for mapping
		int sres = 1000;
		int tres = 2000;
		float[] subhist = new float[sres+1];
		float[] tmphist = new float[tres+1];
		float[] mapping = new float[sres+1];
		for (int n=0;n<sres;n++) subhist[n] = 0;
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) if (!skip || (skip && sub[x][y][z]!=0)) {
			int ns = Numerics.floor( (sub[x][y][z]-smin)/(smax-smin)*sres);
			subhist[ns]+= 1.0f/nssample;
		}
		for (int n=1;n<sres;n++) subhist[n] += subhist[n-1];
		
		for (int n=0;n<tres;n++) tmphist[n] = 0;
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) if (!skip || (skip && tmp[x][y][z]!=0)) {
			int nt = Numerics.floor( (tmp[x][y][z]-tmin)/(tmax-tmin)*tres);
			tmphist[nt]+= 1.0f/ntsample;
		}
		for (int n=1;n<tres;n++) tmphist[n] += tmphist[n-1];
		
		// find the correspondences
		boolean project = projectParam.getValue().booleanValue();
		
		int nt=0;
		for (int n=0;n<sres;n++) {
			while (nt<tres && tmphist[nt] < subhist[n]) nt++;
			// mapping: use a weighted sum (linear approximation)
			float ratio = 0.5f;
			if (nt==0) ratio=0.0f;
			else if (tmphist[nt]>tmphist[nt-1]) ratio = (tmphist[nt]-subhist[n])/(tmphist[nt]-tmphist[nt-1]);
			if (project) mapping[n] = tmin + (ratio*(nt-1) + (1-ratio)*nt)/(float)(tres)*(tmax-tmin);
			else mapping[n] = smin + (ratio*(nt-1) + (1-ratio)*nt)/(float)(tres)*(smax-smin);
		}
		if (project) mapping[sres] = tmax;
		else mapping[sres] = smax;
		
		float[][][] res = new float[nsx][nsy][nsz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			if (!skip || (skip && sub[x][y][z]!=0)) {
				float ratio = (sub[x][y][z]-smin)/(smax-smin)*sres;
				int ns = Numerics.floor(ratio);
				ratio -= ns;
				// linear mapping here too
				if (ns<sres) res[x][y][z] = (1.0f-ratio)*mapping[ns]+ratio*mapping[ns+1];
				else res[x][y][z] = mapping[ns];
			} else {
				res[x][y][z] = 0.0f;	
			}
		}
		ImageDataFloat resultVol = new ImageDataFloat(res);		
		resultVol.setHeader(subImg.getHeader());
		resultVol.setName(subImg.getName()+"_histm");
		resultVolParam.setValue(resultVol);
	}
}
