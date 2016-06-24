package de.mpg.cbs.core.intensity;

/**
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
**/

import de.mpg.cbs.utilities.*;

/*
 * @author Pierre-Louis bazin (bazin@cbs.mpg.de)
 *
 */
 

public class IntensityHistogramMatching {

	//inputs
	private float[] subParam;
	private float[] tmpParam;
	private float[] resultVolParam;
	private boolean zeroParam;
	private boolean projectParam;

	private int nsx, nsy, nsz, nsxyz;
	private float rsx, rsy, rsz;
	private int ntx, nty, ntz, ntxyz;
	private float rtx, rty, rtz;
		
	// set inputs
	public final void setImageToNormalize(float[] val) { subParam = val; }
	public final void setTemplateToMatch(float[] val) { tmpParam = val; }
	public final void setSkipZeroValues(boolean val) { zeroParam = val; }
	public final void setOutputInTemplateRange(boolean val) { projectParam = val; }
	
	// set generic inputs	
	public final void setDimensions(int x, int y, int z) { nsx=x; nsy=y; nsz=z; nsxyz=nsx*nsy*nsz; ntx=x; nty=y; ntz=z; ntxyz=ntx*nty*ntz; }
	public final void setDimensions(int[] dim) { nsx=dim[0]; nsy=dim[1]; nsz=dim[2]; nsxyz=nsx*nsy*nsz; ntx=dim[0]; nty=dim[1]; ntz=dim[2]; ntxyz=ntx*nty*ntz; }
	public final void setResolutions(float x, float y, float z) { rsx=x; rsy=y; rsz=z; rtx=x; rty=y; rtz=z; }
	public final void setResolutions(float[] res) { rsx=res[0]; rsy=res[1]; rsz=res[2]; rtx=res[0]; rty=res[1]; rtz=res[2]; }
	
	public final void setImageDimensions(int x, int y, int z) { nsx=x; nsy=y; nsz=z; nsxyz=nsx*nsy*nsz; }
	public final void setImageDimensions(int[] dim) { nsx=dim[0]; nsy=dim[1]; nsz=dim[2]; nsxyz=nsx*nsy*nsz; }
	public final void setImageResolutions(float x, float y, float z) { rsx=x; rsy=y; rsz=z; }
	public final void setImageResolutions(float[] res) { rsx=res[0]; rsy=res[1]; rsz=res[2]; }
	
	public final void setTemplateDimensions(int x, int y, int z) { ntx=x; nty=y; ntz=z; ntxyz=ntx*nty*ntz; }
	public final void setTemplateDimensions(int[] dim) { ntx=dim[0]; nty=dim[1]; ntz=dim[2]; ntxyz=ntx*nty*ntz; }
	public final void setTemplateResolutions(float x, float y, float z) { rtx=x; rty=y; rtz=z; }
	public final void setTemplateResolutions(float[] res) { rtx=res[0]; rty=res[1]; rtz=res[2]; }
	
	//set JIST definitions
	//to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Intensity"; }
	public final String getLabel() { return "Histogram Matching"; }
	public final String getName() { return "HistogramMatching"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Normalizes the intensity histogram of an image to a template histogram"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
		
	// set outputs
	public float[] getMatchedImage() { return resultVolParam;}
		
	public void execute() {
	
		float[][][] tmp = new float[ntx][nty][ntz];
		for (int x=0;x<ntx;x++) for (int y=0;y<nty;y++) for (int z=0;z<ntz;z++) {
			tmp[x][y][z] = tmpParam[x+ntx*y+ntx*nty*z];
		}
		
		float[][][] sub = new float[nsx][nsy][nsz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			sub[x][y][z] = subParam[x+nsx*y+nsx*nsy*z];
		}
		// build a cumulative mapping
		boolean skip = zeroParam;
		
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
		boolean project = projectParam;
		
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
		resultVolParam = new float[nsxyz];
		for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int z=0;z<nsz;z++) {
			resultVolParam[x+nsx*y+nsx*nsy*z] = res[x][y][z];
		}
	}
}
