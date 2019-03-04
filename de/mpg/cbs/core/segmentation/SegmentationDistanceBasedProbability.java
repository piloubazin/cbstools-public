package de.mpg.cbs.core.segmentation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationDistanceBasedProbability {

	// jist containers
	public static final String[] 	mergingTypes = {"product","minimum","sum","none"};
	public static final String 		mergingType = "product";
	
	private float[] probaImage;
	private int[] labelImage;
	private int[] bgmaskImage;
	private float[] mgdmImage;
	
	private float bgscaleParam = 5.0f;
	private float bgprobaParam = 0.5f;
	private float ratioParam = 0.5f;
	private boolean bgincludedParam = true;
	private String mergeParam;
	private boolean maxprobaParam = true;
	
	private int[] segImage = null;
	private float[] priorImage = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;
	
	private int nlabels = 0;
	private int nmax = 4;
	
	// create inputs
	public final void setSegmentationImage(int[] val) { segImage = val; }
	public final void setPriorProbabilityImage(float[] val) { priorImage = val; }
	public final void setBackgroundDistance_mm(float val) { bgscaleParam = val; }
	public final void setBackgroundProbability(float val) { bgprobaParam = val; }
	public final void setDistanceRatio(float val) { ratioParam = val; }
	public final void setBackgroundIncluded(boolean val) { bgincludedParam = val; }
	public final void setProbabilityMerging(String val) { mergeParam = val; }
	public final void setComputeMaxProbability(boolean val) { maxprobaParam = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Segmentation"; }
	public final String getLabel() { return "Distance-based Probability"; }
	public final String getName() { return "DistanceBasedProbability"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Convert a segmentation map into probabilities based on distance to the inside of each object"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	public final float[] getProbabilityImage() { return probaImage; }
	public final int[] getMaxLabelImage() { return labelImage; }
	public final int[] getBackgroundMaskImage() { return bgmaskImage; }
	public final float[] getMgdmImage() { return mgdmImage; }
	public final int getLabelNumber() { return nlabels; }
	public final int getMaxProbaNumber() { return nmax; }
	
	public void execute() {
		
		
		if (segImage==null) {
			BasicInfo.displayMessage("build segmentation from priors...\n");
			int nlb = priorImage.length/nxyz;
			segImage = new int[nxyz];
			for (int xyz=0;xyz<nxyz;xyz++) {
				float maxproba = 0.0f;
				int maxlb = -1;
				if (!bgincludedParam) {
					maxproba = bgprobaParam;
					maxlb = 0;
				}
				for (int l=0;l<nlb;l++) if (priorImage[xyz+l*nxyz]>maxproba) {
					maxproba = priorImage[xyz+l*nxyz];
					maxlb = (l+1);
				}
				segImage[xyz] = maxlb;
			}
		}
		
		BasicInfo.displayMessage("build from segmentation...\n");
		int[] objlb = ObjectLabeling.listOrderedLabels(segImage,nx,ny,nz);
		nlabels = (byte)objlb.length;
			
		BasicInfo.displayMessage("found "+nlabels+" labels\n");
		// create a distance-based probability map
		float[] boundary = new float[nxyz];
		for (int xyz=0;xyz<nxyz;xyz++) boundary[xyz] = 0.5f;
		
		byte[] lbs = new byte[nlabels];
		for (int l=0;l<nlabels;l++) lbs[l] = (byte)objlb[l];
		BasicInfo.displayMessage("distance-based MGDM representation...\n");
		MgdmRepresentation mgdm = new MgdmRepresentation(segImage, boundary, nx,ny,nz, rx,ry,rz, lbs, nlabels, 4, false, 9.0f);
			
		boolean[] mask = mgdm.getMask();
		
		if (maxprobaParam) {
		    probaImage = new float[nmax*nxyz];
		    labelImage = new int[nmax*nxyz];
		} else {
            probaImage = new float[nlabels*nxyz];
            labelImage = null;
        }
		mgdmImage = new float[nxyz];
		float[] maxobjdist = new float[nlabels];
		float[] proba = new float[nlabels];
		for (int xyz=0;xyz<nxyz;xyz++) {
		    mgdmImage[xyz] =  mgdm.getFunctions()[0][xyz];
		    if (mask[xyz]) for (int l=0;l<nlabels;l++) if (segImage[xyz] == objlb[l]) {
				//if (mgdm.getFunctions()[0][xyz]>maxobjdist[l]) maxobjdist[l] = mgdm.getFunctions()[0][xyz];
				if (mgdmImage[xyz]>maxobjdist[l]) maxobjdist[l] = mgdmImage[xyz];
			}
		}
		
		// background: use the largest distance from foreground object instead as basis ? Or just 1-sum?
		for (int xyz=0;xyz<nxyz;xyz++) {
			if (mask[xyz]) {
				// background : go to a given distance? or go to closest neighborś distance ratio?
				float dist;
				if (bgscaleParam>0) {
					dist = mgdm.reconstructedLevelSetAt(xyz, (byte)0)*Numerics.min(rx,ry,rz);
					if (dist<-bgscaleParam) proba[0] = 1.0f;
					else if (dist<0) proba[0] = 0.5f - 0.5f*dist/bgscaleParam;
					else if (dist<bgscaleParam) proba[0] = 0.5f - 0.5f*dist/bgscaleParam;
					else proba[0] = 0.0f;
				} else {
					dist = mgdm.reconstructedLevelSetAt(xyz, (byte)0);
					int l0 = mgdm.getLabels()[0][xyz];
					if (l0==0) {
						l0 = mgdm.getLabels()[1][xyz];
						if (dist<-ratioParam*maxobjdist[l0]) proba[0] = 1.0f;
						else if (dist<0) proba[0] = 0.5f - 0.5f*dist/maxobjdist[l0]/ratioParam;
						else if (dist<ratioParam*maxobjdist[l0]) proba[0] = 0.5f - 0.5f*dist/maxobjdist[l0]/ratioParam;
						else proba[0] = 0.0f;							
					} else {
						if (dist<-ratioParam*maxobjdist[l0]) proba[0] = 0.0f;
						else if (dist<0) proba[0] = 0.5f + 0.5f*dist/maxobjdist[l0]/ratioParam;
						else if (dist<ratioParam*maxobjdist[l0]) proba[0] = 0.5f + 0.5f*dist/maxobjdist[l0]/ratioParam;
						else proba[0] = 1.0f;
					}
				}
				// foreground: go to a ratio of object's size
				float probaFg = 0.0f;
				for (byte l=1;l<nlabels;l++) {
					dist = mgdm.reconstructedLevelSetAt(xyz, l);
					if (dist<-ratioParam*maxobjdist[l]) proba[l] = 1.0f;
					else if (dist<0) proba[l] = 0.5f - 0.5f*dist/maxobjdist[l]/ratioParam;
					else if (dist<ratioParam*maxobjdist[l]) proba[l] = 0.5f - 0.5f*dist/maxobjdist[l]/ratioParam;
					else proba[l] = 0.0f;
					
					probaFg += proba[l];
				}	
				proba[0] = Numerics.min(proba[0], Numerics.max(1.0f-probaFg, 0.0f) );
				if (maxprobaParam) {
				    for (int n=0;n<nmax;n++) {
                        int best=0;
                        for (int l=0;l<nlabels;l++) {
                            if (proba[l]>proba[best]) best = l;
                        }
                        probaImage[n*nxyz+xyz] = proba[best];
                        labelImage[n*nxyz+xyz] = best;
                        
                        proba[best] = 0.0f;
                    }
				} else {
                    for (byte l=0;l<nlabels;l++) probaImage[l*nxyz+xyz] = proba[l];				    
				}
			} else {
				if (maxprobaParam) {
				    probaImage[xyz] = 1.0f;
                    labelImage[xyz] = 0;
                    for (byte n=1;n<nmax;n++) {
                        probaImage[n*nxyz+xyz] = 0.0f;
                        labelImage[n*nxyz+xyz] = n;
                    }
				} else {
                    probaImage[xyz] = 1.0f;
                    for (byte l=1;l<nlabels;l++) probaImage[l*nxyz+xyz] = 0.0f;
                }
			}
		}
		
		if (priorImage!=null) {
			BasicInfo.displayMessage("merging distance-based proba and prior...\n");
			for (int xyz=0;xyz<nxyz;xyz++) {
			    byte l0 = 0;
                if (!bgincludedParam) l0 = 1;
                
                if (maxprobaParam) {
                    for (int l=0;l<nlabels;l++) proba[l] = 0.0f;
                    for (byte n=0;n<nmax;n++) proba[labelImage[n*nxyz+xyz]] = probaImage[n*nxyz+xyz];
                } else {
                    for (int l=0;l<nlabels;l++) proba[l] = probaImage[l*nxyz+xyz];
                }
                if (mergeParam.equals("product")) {
                    for (byte l=l0;l<nlabels;l++) {
                        proba[l] = 	(float)FastMath.sqrt(proba[l]*priorImage[(l-l0)*nxyz+xyz]);
                    }
                } else if (mergeParam.equals("minimum")) {
                    for (byte l=l0;l<nlabels;l++) {
                        proba[l] = 	Numerics.min(proba[l], priorImage[(l-l0)*nxyz+xyz]);
                    }
                } else if (mergeParam.equals("sum")) {
                    for (byte l=l0;l<nlabels;l++) {
                        proba[l] = 	0.5f*(proba[l]+priorImage[(l-l0)*nxyz+xyz]);
                    }
                }
                
                if (maxprobaParam) {
				    for (int n=0;n<nmax;n++) {
                        int best=0;
                        for (int l=0;l<nlabels;l++) {
                            if (proba[l]>proba[best]) best = l;
                        }
                        probaImage[n*nxyz+xyz] = proba[best];
                        labelImage[n*nxyz+xyz] = best;
                        
                        proba[best] = 0.0f;
                    }
				} else {
                    for (byte l=0;l<nlabels;l++) probaImage[l*nxyz+xyz] = proba[l];				    
				}
            }
		}
		
		// compute a mask of non-one background values
		bgmaskImage = new int[nxyz];
		int b = 1;
		for (int x=b;x<nx-b;x++) for (int y=b;y<ny-b;y++) for (int z=b;z<nz-b;z++) {
			int xyz = x+nx*y+nx*ny*z;
			if (maxprobaParam) {
			    if (labelImage[xyz]!=0 || probaImage[xyz]<1) bgmaskImage[xyz] = 1;
			} else {
			    if (probaImage[xyz]<1) bgmaskImage[xyz] = 1;
			}
		}
		
		return;
	}
	
}