package de.mpg.cbs.core.segmentation;

import java.net.URL;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;


/*
 * @author Pierre-Louis Bazin
 */
public class SegmentationCellMgdm {

	private float[] input1Image = null;
	private float[] input2Image = null;
	private float[] input3Image = null;
	
	private int nx, ny, nz, nxyz;
	private float rx, ry, rz;

	private String type1Param = "none";
	private String type2Param = "none";
	private String type3Param = "none";
	public static final String[] inputTypes = {"centroid-proba", "local-maxima", "foreground-proba","image-intensities"};
	
	private String dimParam = "none";
	public static final String[] dimTypes = {"2D", "3D"};
													
	private int 	iterationParam	=	500;
	private float 	changeParam		=	0.001f;
	private	float 	forceParam		= 	0.1f;
	private float 	curvParam		=	0.4f;
	private float   cellthresholdParam = 0.8f;
	private float   centroidthresholdParam = 0.1f;
		
	private String 	topologyParam	=	"wcs";
	public static final String[] topoTypes = {"26/6", "6/26", "18/6", "6/18", "6/6", "wcs", "wco", "no"};
	private String	lutdir = null;
	
	// outputs
	private int[] segmentImage;
	private float[] mgdmImage;
	
	// create inputs
	public final void setContrastImage1(float[] val) { input1Image = val; }
	public final void setContrastType1(String val) { type1Param = val; }
	public final void setContrastImage2(float[] val) { input2Image = val; }
	public final void setContrastType2(String val) { type2Param = val; }
	public final void setContrastImage3(float[] val) { input3Image = val; }
	public final void setContrastType3(String val) { type3Param = val; }
	
	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }
	
	public final void setDataWeight(float val) {forceParam = val; }
	public final void setCurvatureWeight(float val) { curvParam = val; }
	public final void setMaxIterations(int val) { iterationParam = val; }
	public final void setMinChange(float val) { changeParam = val; }
	public final void setCellThreshold(float val) { cellthresholdParam = val; }

	public final void setDataStackDimension(String val) { dimParam = val; }

	public final void setTopology(String val) { topologyParam = val; }
	public final void setTopologyLUTdirectory(String val) { lutdir = val; }
	
	// to be used for JIST definitions, generic info / help
	public static final String getPackage() { return "CBS Tools"; }
	public static final String getCategory() { return "Segmentation.devel"; }
	public static final String getLabel() { return "MGDM Cell Segmentation"; }
	public static final String getName() { return "MGDMCellSegmentation"; }

	public static final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public static final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public static final String getDescription() { return "MGDM estimation of cells from deep learning inputs."; }
	public static final String getLongDescription() { return getDescription(); }
		
	public static final String getVersion() { return "3.1.2"; }

	// create outputs
	public final int[] getSegmentedImage() { return segmentImage; }
	public final float[] getLevelsetBoundaryImage() { return mgdmImage; }

	public void execute(){
		// import the image data into 1D arrays
		int nimg = 1;
		if (input2Image != null) nimg++;
		if (input3Image != null) nimg++;
		
		String[] modality = new String[nimg];
		float[][] image = new float[nimg][];
		int n=0;
		modality[n] = type1Param;
		image[n] = input1Image;
		if (input2Image != null) { n++; modality[n] = type2Param; image[n] = input2Image; }
		if (input3Image != null) { n++; modality[n] = type3Param; image[n] = input3Image; }
		
		float[] centroids = null;
		float[] maxima = null;
		float[] fgproba = null;
		float[] intens = null;
		for (n=0;n<nimg;n++) {
		    if (modality[n].equals("centroid-proba")) centroids = image[n];
		    else if (modality[n].equals("local-maxima")) maxima = image[n];
		    else if (modality[n].equals("foreground-proba")) fgproba = image[n];
		    else if (modality[n].equals("image-intensities")) intens = image[n];
		}
		
		// different streams for 2D or 3D
		byte DIM2D = 2;
		byte DIM3D = 3;
		byte dimension = DIM2D;
		if (dimParam.equals("3D")) dimension = DIM3D;
		
		int nmgdm = 4;
                
		//output images
		segmentImage = new int[nx*ny*nz];
		mgdmImage = new float[nx*ny*nz];
		
		if (dimension==DIM2D) {
		    // independet Z stacks
		    for (int z=0;z<nz;z++) {
		        // 1. Get the initial segmentation from cell centroid detection + local maxima
		        boolean[] seg = new boolean[nx*ny];
		        for (int xy=0;xy<nx*ny;xy++) seg[xy] = (maxima[xy+z*nx*ny]>0);
		        int[] initialization = ObjectLabeling.connected6Object3D(seg, nx,ny,1);
                // simple region growing: could be done more nicely
		        int[] growing = new int[nx*ny];
		        float[] proba = new float[nx*ny];
		        float[] maxproba = new float[nx*ny];
		        for (int xy=0;xy<nx*ny;xy++) proba[xy] = centroids[xy+z*nx*ny];
		        for (int t=0;t<20;t++) {
		            for (int xy=0;xy<nx*ny;xy++) {
		                growing[xy] = initialization[xy];
		                maxproba[xy] = proba[xy];
		            }
		            for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) {
		                int xy = x+nx*y;
		                if (initialization[xy]>0) {
		                    if (growing[xy-1]==0 && proba[xy-1]>centroidthresholdParam*proba[xy]) {
		                        growing[xy-1] = initialization[xy];
		                        maxproba[xy-1] = Numerics.max(maxproba[xy-1], proba[xy-1], proba[xy]);
		                    }
		                    if (growing[xy+1]==0 && proba[xy+1]>centroidthresholdParam*proba[xy]) {
		                        growing[xy+1] = initialization[xy];
		                        maxproba[xy+1] = Numerics.max(maxproba[xy+1], proba[xy+1], proba[xy]);
		                    }
		                    if (growing[xy-nx]==0 && proba[xy-nx]>centroidthresholdParam*proba[xy]) {
		                        growing[xy-nx] = initialization[xy];
		                        maxproba[xy-nx] = Numerics.max(maxproba[xy-nx], proba[xy-nx], proba[xy]);
		                    }
		                    if (growing[xy+nx]==0 && proba[xy+nx]>centroidthresholdParam*proba[xy]) {
		                        growing[xy+nx] = initialization[xy];
		                        maxproba[xy+nx] = Numerics.max(maxproba[xy+nx], proba[xy+nx], proba[xy]);
		                    }
		                }
		            }
		            for (int xy=0;xy<nx*ny;xy++) {
		                initialization[xy] = growing[xy];
		                proba[xy] = maxproba[xy];
		            }
		        }
		            
                int nlb = ObjectLabeling.countLabels(initialization, nx,ny,1);
                System.out.print("slice "+z+", "+nlb+" labels \n");
		        
                // 2. Get the forces from foreground proba
                float[][] forces = new float[nlb][];
                if (intens==null) {
                    float[] fgmap = new float[nx*ny];
                    float[] bgmap = new float[nx*ny];
                    for (int xy=0;xy<nx*ny;xy++) {
                        fgmap[xy] = fgproba[xy+z*nx*ny];
                        bgmap[xy] = 1.0f - fgproba[xy+z*nx*ny];
                    }
                    // trick so that each object has the same proba (foreground)
                    for (int lb=0;lb<nlb;lb++) forces[lb] = fgmap;
               
                    // find background
                    // assume background to be the label with highest cumulative bg proba?
                    float[] bgscore = new float[nlb];
                    for (int lb=0;lb<nlb;lb++)
                        for (int xy=0;xy<nx*ny;xy++)
                            bgscore[lb] += bgmap[xy];
                    int bglb = Numerics.argmax(bgscore);
                    forces[bglb] = bgmap;
                } else {
                    // different model: propagate intensities to 50% of original intensity per cluster
                    float[] initmax = new float[nlb];
                    int bglb = 0;
                    for (int xy=0;xy<nx*ny;xy++) {
                        if (initialization[xy]>0) {
                            initmax[initialization[xy]] = Numerics.max(initmax[initialization[xy]],intens[xy+z*nx*ny]);
                        } else {
                            bglb = initialization[xy];
                        }
                    }
                    // for the background, use the mean (over-estimating)?
                    double sumbg = 0.0, denbg = 0.0;
                    for (int xy=0;xy<nx*ny;xy++) {
                        if (initialization[xy]==0) {
                            sumbg += intens[xy+z*nx*ny];
                            denbg++;
                        }
                    }
                    float initbg = (float)(sumbg/denbg);
                    for (int lb=0;lb<nlb;lb++) {
                        if (lb!=bglb) {
                            forces[lb] = new float[nx*ny];
                            for (int xy=0;xy<nx*ny;xy++) {
                                //forces[lb][xy] = Numerics.bounded( (intens[xy+z*nx*ny]/initmax[lb]-cellthresholdParam)/(1.0f-cellthresholdParam), -1.0f, 1.0f);
                                forces[lb][xy] = Numerics.bounded( ( (1.0f-cellthresholdParam)*initmax[lb] - Numerics.abs(initmax[lb]-intens[xy+z*nx*ny]))
                                                                        /( (1.0f-cellthresholdParam)*initmax[lb]), -1.0f, 1.0f);
                            }
                        } else {
                            forces[lb] = new float[nx*ny];
                            for (int xy=0;xy<nx*ny;xy++) {
                                // not bad..
                                //forces[lb][xy] = Numerics.bounded( (float)((initbg-intens[xy+z*nx*ny])/initbg), -1.0f, 1.0f);
                                forces[lb][xy] = Numerics.bounded( ( (1.0f-cellthresholdParam)*initbg -intens[xy+z*nx*ny] + initbg)
                                                                        /( (1.0f-cellthresholdParam)*initbg), -1.0f, 1.0f);
                            }
                        }
                    }
                }
                    
                // 3. Run MGDM!
                Mgdm2d mgdm = new Mgdm2d(initialization, nx, ny, nlb, nmgdm, rx, ry, null, forces, 
						                0.0f, forceParam, curvParam, 0.0f, 
                                        topologyParam, lutdir);
                
                mgdm.evolveNarrowBand(iterationParam, changeParam);
                
                // 4. copy the results
                for (int xy=0;xy<nx*ny;xy++) {
                    segmentImage[xy+z*nx*ny] = mgdm.getLabels()[0][xy];
                    mgdmImage[xy+z*nx*ny] = mgdm.getFunctions()[0][xy];
                }
		    }
		} else if (dimension==DIM3D) {
            // 1. Get the initial segmentation from cell centroid detection + local maxima
            boolean[] seg = new boolean[nx*ny*nz];
            for (int xyz=0;xyz<nx*ny*nz;xyz++) seg[xyz] = (maxima[xyz]>0);
            int[] initialization = ObjectLabeling.connected6Object3D(seg, nx,ny,nz);
            // simple region growing: could be done more nicely
            int[] growing = new int[nx*ny*nz];
            float[] proba = new float[nx*ny*nz];
            float[] maxproba = new float[nx*ny*nz];
            for (int xyz=0;xyz<nx*ny*nz;xyz++) proba[xyz] = centroids[xyz];
            for (int t=0;t<20;t++) {
                for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                    growing[xyz] = initialization[xyz];
                    maxproba[xyz] = proba[xyz];
                }
                for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
                    int xyz = x+nx*y+nx*ny*z;
                    if (initialization[xyz]>0) {
                        if (growing[xyz-1]==0 && proba[xyz-1]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz-1] = initialization[xyz];
                            maxproba[xyz-1] = Numerics.max(maxproba[xyz-1], proba[xyz-1], proba[xyz]);
                        }
                        if (growing[xyz+1]==0 && proba[xyz+1]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz+1] = initialization[xyz];
                            maxproba[xyz+1] = Numerics.max(maxproba[xyz+1], proba[xyz+1], proba[xyz]);
                        }
                        if (growing[xyz-nx]==0 && proba[xyz-nx]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz-nx] = initialization[xyz];
                            maxproba[xyz-nx] = Numerics.max(maxproba[xyz-nx], proba[xyz-nx], proba[xyz]);
                        }
                        if (growing[xyz+nx]==0 && proba[xyz+nx]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz+nx] = initialization[xyz];
                            maxproba[xyz+nx] = Numerics.max(maxproba[xyz+nx], proba[xyz+nx], proba[xyz]);
                        }
                        if (growing[xyz-nx*ny]==0 && proba[xyz-nx*ny]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz-nx*ny] = initialization[xyz];
                            maxproba[xyz-nx*ny] = Numerics.max(maxproba[xyz-nx*ny], proba[xyz-nx*ny], proba[xyz]);
                        }
                        if (growing[xyz+nx*ny]==0 && proba[xyz+nx*ny]>centroidthresholdParam*proba[xyz]) {
                            growing[xyz+nx*ny] = initialization[xyz];
                            maxproba[xyz+nx*ny] = Numerics.max(maxproba[xyz+nx*ny], proba[xyz+nx*ny], proba[xyz]);
                        }
                    }
                }
                for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                    initialization[xyz] = growing[xyz];
                    proba[xyz] = maxproba[xyz];
                }
            }
                
            int nlb = ObjectLabeling.countLabels(initialization, nx,ny,nz);
            System.out.print("full stack: "+nlb+" labels \n");
            
            // 2. Get the forces from foreground proba
            //float[][] forces = new float[nlb][];
            float[][] forces = new float[1][];
            float[] initmax = new float[nlb];
            int bglb = 0;
            float initbg = 0.0f;
            if (intens==null) {
                float[] fgmap = new float[nx*ny*nz];
                float[] bgmap = new float[nx*ny*nz];
                for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                    fgmap[xyz] = fgproba[xyz];
                    bgmap[xyz] = 1.0f - fgproba[xyz];
                }
                // trick so that each object has the same proba (foreground)
                for (int lb=0;lb<nlb;lb++) forces[lb] = fgmap;
           
                // find background
                // assume background to be the label with highest cumulative bg proba?
                float[] bgscore = new float[nlb];
                for (int lb=0;lb<nlb;lb++)
                    for (int xyz=0;xyz<nx*ny*nz;xyz++)
                        bgscore[lb] += bgmap[xyz];
                bglb = Numerics.argmax(bgscore);
                forces[bglb] = bgmap;
            } else {
                // different model: propagate intensities to 50% of original intensity per cluster
                for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                    if (initialization[xyz]>0) {
                        initmax[initialization[xyz]] = Numerics.max(initmax[initialization[xyz]],intens[xyz]);
                    } else {
                        bglb = initialization[xyz];
                    }
                }
                // for the background, use the mean (over-estimating)?
                double sumbg = 0.0, denbg = 0.0;
                for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                    if (initialization[xyz]==0) {
                        sumbg += intens[xyz];
                        denbg++;
                    }
                }
                initbg = (float)(sumbg/denbg);
                forces[0] = intens;
                /*
                for (int lb=0;lb<nlb;lb++) {
                    if (lb!=bglb) {
                        forces[lb] = new float[nx*ny*nz];
                        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                            //forces[lb][xy] = Numerics.bounded( (intens[xy+z*nx*ny]/initmax[lb]-cellthresholdParam)/(1.0f-cellthresholdParam), -1.0f, 1.0f);
                            forces[lb][xyz] = Numerics.bounded( ( (1.0f-cellthresholdParam)*initmax[lb] - Numerics.abs(initmax[lb]-intens[xyz]))
                                                                    /( (1.0f-cellthresholdParam)*initmax[lb]), -1.0f, 1.0f);
                        }
                    } else {
                        forces[lb] = new float[nx*ny*nz];
                        for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                            // not bad..
                            //forces[lb][xy] = Numerics.bounded( (float)((initbg-intens[xy+z*nx*ny])/initbg), -1.0f, 1.0f);
                            forces[lb][xyz] = Numerics.bounded( ( (1.0f-cellthresholdParam)*initbg -intens[xyz] + initbg)
                                                                    /( (1.0f-cellthresholdParam)*initbg), -1.0f, 1.0f);
                        }
                    }
                }*/
            }
                
            // 3. Run MGDM!
            CellMgdm3d mgdm = new CellMgdm3d(initialization, nx, ny, nz, nlb, nmgdm, rx, ry, rz, null, 
                                    forces, initmax, initbg, bglb,
                                    0.0f, forceParam, curvParam, 0.0f, 
                                    topologyParam, lutdir);
            
            mgdm.evolveNarrowBand(iterationParam, changeParam);
            
            // 4. copy the results
            for (int xyz=0;xyz<nx*ny*nz;xyz++) {
                segmentImage[xyz] = mgdm.getLabels()[0][xyz];
                mgdmImage[xyz] = mgdm.getFunctions()[0][xyz];
            }
        }		    
		return;
	}

}
