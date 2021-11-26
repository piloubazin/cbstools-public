package de.mpg.cbs.core.registration;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin
 */
public class RegistrationApplyAngularDeformations2D {

	// jist containers
	private float[] sourceImage;
	private float[] referenceImage = null;
	private float[] deformation1Image;
	private float[] deformation2Image = null;
	private float[] deformation3Image = null;
	private float[] deformation4Image = null;
	private String type1Option;
	private String type2Option = "none";
	private String type3Option = "none";
	private String type4Option = "none";
	private String interpOption = "nearest";
	private String padOption = "closest";
	
	public static final String[] types = {"none", "deformation(voxels)", "mapping(voxels)", "deformation(mm)", "mapping(mm)"};
	public static final String[] interp = {"nearest", "linear"};
	public static final String[] pads = {"closest", "zero", "min", "max"};
	
	private float[] deformedImage;
	
	private int nsx, nsy, nst, nsxy;
	private float rsx, rsy;
	private int nrx, nry, nrt, nrxy;
	private float rrx, rry;
	
	private int nd1x, nd1y, nd1xy;
	private float rd1x, rd1y;
	private int nd2x, nd2y, nd2xy;
	private float rd2x, rd2y;
	private int nd3x, nd3y, nd3xy;
	private float rd3x, rd3y;
	private int nd4x, nd4y, nd4xy;
	private float rd4x, rd4y;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte T = 2;

	// create inputs
	public final void setImageToDeform(float[] val) { sourceImage = val; }
	//public final void setReferenceImage(float[] val) { referenceImage = val; }
	public final void setDeformationMapping1(float[] val) { deformation1Image = val; }
	public final void setDeformationType1(String val) { type1Option = val; }
	public final void setDeformationMapping2(float[] val) { deformation2Image = val; }
	public final void setDeformationType2(String val) { type2Option = val; }
	public final void setDeformationMapping3(float[] val) { deformation3Image = val; }
	public final void setDeformationType3(String val) { type3Option = val; }
	public final void setDeformationMapping4(float[] val) { deformation4Image = val; }
	public final void setDeformationType4(String val) { type4Option = val; }
	public final void setInterpolationType(String val) { interpOption = val; }
	public final void setImagePadding(String val) { padOption = val; }
		
	
	public final void setImageDimensions(int x, int y) { nsx=x; nsy=y; nst=1; nsxy=nsx*nsy; }
	public final void setImageDimensions(int x, int y, int t) { nsx=x; nsy=y; nst=t; nsxy=nsx*nsy; }
	public final void setImageDimensions(int[] dim) { nsx=dim[0]; nsy=dim[1]; nst=dim[2]; nsxy=nsx*nsy; }
	
	public final void setImageResolutions(float x, float y) { rsx=x; rsy=y; }
	public final void setImageResolutions(float[] res) { rsx=res[0]; rsy=res[1]; }

	//public final void setReferenceDimensions(int x, int y, int z) { nrx=x; nry=y; nrz=z; nrt=t; nrxy=nrx*nry*nrz; }
	//public final void setReferenceDimensions(int[] dim) { nrx=dim[0]; nry=dim[1]; nrz=dim[2]; nrt=dim[3]; nrxy=nrx*nry*nrz; }
	
	//public final void setReferenceResolutions(float x, float y, float z) { rrx=x; rry=y; rrz=z; }
	//public final void setReferenceResolutions(float[] res) { rrx=res[0]; rry=res[1]; rrz=res[2]; }

	public final void setDeformation1Dimensions(int x, int y) { nd1x=x; nd1y=y; nd1xy=nd1x*nd1y; }
	public final void setDeformation1Dimensions(int[] dim) { nd1x=dim[0]; nd1y=dim[1]; nd1xy=nd1x*nd1y; }
	
	public final void setDeformation1Resolutions(float x, float y) { rd1x=x; rd1y=y; }
	public final void setDeformation1Resolutions(float[] res) { rd1x=res[0]; rd1y=res[1]; }

	public final void setDeformation2Dimensions(int x, int y) { nd2x=x; nd2y=y; nd2xy=nd2x*nd2y; }
	public final void setDeformation2Dimensions(int[] dim) { nd2x=dim[0]; nd2y=dim[1]; nd2xy=nd2x*nd2y; }
	
	public final void setDeformation2Resolutions(float x, float y) { rd2x=x; rd2y=y; }
	public final void setDeformation2Resolutions(float[] res) { rd2x=res[0]; rd2y=res[1]; }

	public final void setDeformation3Dimensions(int x, int y) { nd3x=x; nd3y=y; nd3xy=nd3x*nd3y; }
	public final void setDeformation3Dimensions(int[] dim) { nd3x=dim[0]; nd3y=dim[1]; nd3xy=nd3x*nd3y; }
	
	public final void setDeformation3Resolutions(float x, float y) { rd3x=x; rd3y=y; }
	public final void setDeformation3Resolutions(float[] res) { rd3x=res[0]; rd3y=res[1]; }

	public final void setDeformation4Dimensions(int x, int y) { nd4x=x; nd4y=y; nd4xy=nd4x*nd4y; }
	public final void setDeformation4Dimensions(int[] dim) { nd4x=dim[0]; nd4y=dim[1]; nd4xy=nd4x*nd4y; }
	
	public final void setDeformation4Resolutions(float x, float y) { rd4x=x; rd4y=y; }
	public final void setDeformation4Resolutions(float[] res) { rd4x=res[0]; rd4y=res[1]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Registration"; }
	public final String getLabel() { return "Apply Deformations"; }
	public final String getName() { return "ApplyDeformations"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Apply multiple deformations from non-linear transformations. Transformations can be a deformation field or a coordinate mapping."; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.3"; };
	
	public final float[] getDeformedImage() { return deformedImage; }
	
	public void execute() {
				
        // deformation: in reference space
        System.out.println("load deformation 1");
        // scale to voxels if needed
        if (type1Option.endsWith("(mm)")) {
            System.out.println("normalize to resolution ("+rd1x+", "+rd1y+")");
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) {
                int xy = x + nd1x*y;
                deformation1Image[xy + X*nd1xy] /= rd1x; 	
                deformation1Image[xy + Y*nd1xy] /= rd1y; 	
            }
        }
        // turn into a mapping if needed
        if (type1Option.startsWith("deformation")) {
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) {
                int xy = x + nd1x*y;
                deformation1Image[xy + X*nd1xy] += x;
                deformation1Image[xy + Y*nd1xy] += y;
            }
        }
        // check for bad borders
        boolean[] boundary = new boolean[nd1x*nd1y];
        boolean growBoundaries = false;
        for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) {
            int xy = x + nd1x*y;
            if (deformation1Image[xy + X*nd1xy]==0 && deformation1Image[xy + Y*nd1xy]==0) {
                for (byte k=0;k<4;k++) {
                    if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd1x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd1y) {
                        int ngb = Ngb2.neighborIndex(k, xy,nd1x,nd1y,1);
                        if (deformation1Image[ngb + X*nd1xy]!=0 || deformation1Image[ngb + Y*nd1xy]!=0 ) {
                            growBoundaries = true;
                            boundary[ngb] = true;
                            k=4;
                        }
                    }
                }
            }
        }
        while (growBoundaries) {
            boolean[] changed = new boolean[nd1x*nd1y];
            growBoundaries = false;
            for (int x=0;x<nd1x;x++) for (int y=0;y<nd1y;y++) {
                int xy = x + nd1x*y;
                if (boundary[xy]) {
                    for (byte k=0;k<4;k++) {
                        if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd1x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd1y) {
                            int ngb = Ngb2.neighborIndex(k, xy,nd1x,nd1y,1);
                            if (deformation1Image[ngb + X*nd1xy]==0 && deformation1Image[ngb + Y*nd1xy]==0) {
                                deformation1Image[ngb + X*nd1xy] = deformation1Image[xy + X*nd1xy];
                                deformation1Image[ngb + Y*nd1xy] = deformation1Image[xy + Y*nd1xy];
                                growBoundaries = true;
                                changed[ngb] = true;
                            }
                        }
                    }
                }
            }
            boundary = changed;
        }
        
        // add to final deformation
        float[] deformation = deformation1Image;
        
        nrx = nd1x; nry = nd1y;
        rrx = rd1x; rry = rd1y;
        
        // second deformation if given
        if (!type2Option.equals("none") && deformation2Image!=null) {
            System.out.println("load deformation 2");
            // scale to voxels if needed
            if (type2Option.endsWith("(mm)")) {
                System.out.println("normalize to resolution ("+rd2x+", "+rd2y+")");
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) {
                    int xy = x + nd2x*y;
                    deformation2Image[xy + X*nd2xy] /= rd2x; 	
                    deformation2Image[xy + Y*nd2xy] /= rd2y; 	
                }
            }
            // turn into a mapping if needed
            if (type2Option.startsWith("deformation")) {
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) {
                    int xy = x + nd2x*y;
                    deformation2Image[xy + X*nd2xy] += x;
                    deformation2Image[xy + Y*nd2xy] += y;
                }
            }
            // check for bad borders
            boundary = new boolean[nd2x*nd2y];
            growBoundaries = false;
            for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) {
                int xy = x + nd2x*y;
                if (deformation2Image[xy + X*nd2xy]==0 && deformation2Image[xy + Y*nd2xy]==0) {
                    for (byte k=0;k<4;k++) {
                        if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd2x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd2y) {
                            int ngb = Ngb2.neighborIndex(k, xy,nd2x,nd2y,1);
                            if (deformation2Image[ngb + X*nd2xy]!=0 || deformation2Image[ngb + Y*nd2xy]!=0) {
                                growBoundaries = true;
                                boundary[ngb] = true;
                                k=6;
                            }
                        }
                    }
                }
            }
            while (growBoundaries) {
                boolean[] changed = new boolean[nd2x*nd2y];
                growBoundaries = false;
                for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) {
                    int xy = x + nd2x*y;
                    if (boundary[xy]) {
                        for (byte k=0;k<4;k++) {
                            if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd2x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd2y) {
                                int ngb = Ngb2.neighborIndex(k, xy,nd2x,nd2y,1);
                                if (deformation2Image[ngb + X*nd2xy]==0 && deformation2Image[ngb + Y*nd2xy]==0) {
                                    deformation2Image[ngb + X*nd2xy] = deformation2Image[xy + X*nd2xy];
                                    deformation2Image[ngb + Y*nd2xy] = deformation2Image[xy + Y*nd2xy];
                                    growBoundaries = true;
                                    changed[ngb] = true;
                                }
                            }
                        }
                    }
                }
                boundary = changed;
            }
        
            // compose the deformations: X' = def1(def2(X))
            System.out.println("compose deformations");
            float[] composed12 = new float[nd2x*nd2y*2];
            for (int x=0;x<nd2x;x++) for (int y=0;y<nd2y;y++) {
                int xy = x + nd2x*y;
                composed12[xy+X*nd2xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation2Image[xy+X*nd2xy], deformation2Image[xy+Y*nd2xy], X, nd1x, nd1y, 2);
                composed12[xy+Y*nd2xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation2Image[xy+X*nd2xy], deformation2Image[xy+Y*nd2xy], Y, nd1x, nd1y, 2);
            }
            deformation = composed12;
            deformation1Image = null;
            deformation2Image = null;

            nrx = nd2x; nry = nd2y;
            rrx = rd2x; rry = rd2y;

            // third deformation if given
            if (!type3Option.equals("none") && deformation3Image!=null) {
                System.out.println("load deformation 3");
                // scale to voxels if needed
                if (type3Option.endsWith("(mm)")) {
                    System.out.println("normalize to resolution ("+rd3x+", "+rd3y+")");
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) {
                        int xy = x + nd3x*y;
                        deformation3Image[xy + X*nd3xy] /= rd3x; 	
                        deformation3Image[xy + Y*nd3xy] /= rd3y; 	
                    }
                }
                // turn into a mapping if needed
                if (type3Option.startsWith("deformation")) {
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) {
                        int xy = x + nd3x*y;
                        deformation3Image[xy + X*nd3xy] += x;
                        deformation3Image[xy + Y*nd3xy] += y;
                    }
                }
                // check for bad borders
                boundary = new boolean[nd3x*nd3y];
                growBoundaries = false;
                for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) {
                    int xy = x + nd3x*y;
                    if (deformation3Image[xy + X*nd3xy]==0 && deformation3Image[xy + Y*nd3xy]==0) {
                        for (byte k=0;k<4;k++) {
                            if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd3x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd3y) {
                                int ngb = Ngb2.neighborIndex(k, xy,nd3x,nd3y,1);
                                if (deformation3Image[ngb + X*nd3xy]!=0 || deformation3Image[ngb + Y*nd3xy]!=0) {
                                    growBoundaries = true;
                                    boundary[ngb] = true;
                                    k=4;
                                }
                            }
                        }
                    }
                }
                while (growBoundaries) {
                    boolean[] changed = new boolean[nd3x*nd3y];
                    growBoundaries = false;
                    for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) {
                        int xy = x + nd3x*y;
                        if (boundary[xy]) {
                            for (byte k=0;k<4;k++) {
                                if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd3x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd3y) {
                                    int ngb = Ngb2.neighborIndex(k, xy,nd3x,nd3y,1);
                                    if (deformation3Image[ngb + X*nd3xy]==0 && deformation3Image[ngb + Y*nd3xy]==0) {
                                        deformation3Image[ngb + X*nd3xy] = deformation3Image[xy + X*nd3xy];
                                        deformation3Image[ngb + Y*nd3xy] = deformation3Image[xy + Y*nd3xy];
                                        growBoundaries = true;
                                        changed[ngb] = true;
                                    }
                                }
                            }
                        }
                    }
                    boundary = changed;
                }
            
                // compose the deformations: X' = def1(def2(def3(X)))
                System.out.println("compose deformations");
                float[] composed123 = new float[nd3x*nd3y*2];
                for (int x=0;x<nd3x;x++) for (int y=0;y<nd3y;y++) {
                    int xy = x + nd3x*y;
                    composed123[xy+X*nd3xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation3Image[xy+X*nd3xy], deformation3Image[xy+Y*nd3xy], X, nd2x, nd2y, 2);
                    composed123[xy+Y*nd3xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation3Image[xy+X*nd3xy], deformation3Image[xy+Y*nd3xy], Y, nd2x, nd2y, 2);
                }
                deformation = composed123;
                deformation3Image = null;
                composed12 = null;
                
                nrx = nd3x; nry = nd3y;
                rrx = rd3x; rry = rd3y;
                
                if (!type4Option.equals("none") && deformation4Image!=null) {
                    System.out.println("load deformation 4");
                    // scale to voxels if needed
                    if (type4Option.endsWith("(mm)")) {
                        System.out.println("normalize to resolution ("+rd4x+", "+rd4y+")");
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) {
                            int xy = x + nd4x*y;
                            deformation4Image[xy + X*nd4xy] /= rd4x; 	
                            deformation4Image[xy + Y*nd4xy] /= rd4y; 	
                        }
                    }
                    // turn into a mapping if needed
                    if (type4Option.startsWith("deformation")) {
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) {
                            int xy = x + nd4x*y;
                            deformation4Image[xy + X*nd4xy] += x;
                            deformation4Image[xy + Y*nd4xy] += y;
                        }
                    }
                    // check for bad borders
                    boundary = new boolean[nd4x*nd4y];
                    growBoundaries = false;
                    for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) {
                        int xy = x + nd4x*y;
                        if (deformation4Image[xy + X*nd4xy]==0 && deformation4Image[xy + Y*nd4xy]==0) {
                            for (byte k=0;k<4;k++) {
                                if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd4x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd4y) {
                                    int ngb = Ngb2.neighborIndex(k, xy,nd4x,nd4y,1);
                                    if (deformation4Image[ngb + X*nd4xy]!=0 || deformation4Image[ngb + Y*nd4xy]!=0) {
                                        growBoundaries = true;
                                        boundary[ngb] = true;
                                        k=4;
                                    }
                                }
                            }
                        }
                    }
                    while (growBoundaries) {
                        boolean[] changed = new boolean[nd4x*nd4y];
                        growBoundaries = false;
                        for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) {
                            int xy = x + nd4x*y;
                            if (boundary[xy]) {
                                for (byte k=0;k<4;k++) {
                                    if (x+Ngb2.x[k]>=0 && x+Ngb2.x[k]<nd4x && y+Ngb2.y[k]>=0 && y+Ngb2.y[k]<nd4y) {
                                        int ngb = Ngb2.neighborIndex(k, xy,nd4x,nd4y,1);
                                        if (deformation4Image[ngb + X*nd4xy]==0 && deformation4Image[ngb + Y*nd4xy]==0) {
                                            deformation4Image[ngb + X*nd4xy] = deformation4Image[xy + X*nd4xy];
                                            deformation4Image[ngb + Y*nd4xy] = deformation4Image[xy + Y*nd4xy];
                                            growBoundaries = true;
                                            changed[ngb] = true;
                                        }
                                    }
                                }
                            }
                        }
                        boundary = changed;
                    }
                    // compose the deformations: X' = def1(def2(def3(X)))
                    System.out.println("compose deformations");
                    float[] composed1234 = new float[nd4x*nd4y*2];
                    for (int x=0;x<nd4x;x++) for (int y=0;y<nd4y;y++) {
                        int xy = x + nd4x*y;
                        composed1234[xy+X*nd4xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation4Image[xy+X*nd4xy], deformation4Image[xy+Y*nd4xy], X, nd3x, nd3y, 2);
                        composed1234[xy+Y*nd4xy] = ImageInterpolation.linearClosestInterpolation(deformation, deformation4Image[xy+X*nd4xy], deformation4Image[xy+Y*nd4xy], Y, nd3x, nd3y, 2);
                    }
                    deformation = composed1234;
                    deformation4Image = null;
                    composed123 = null;
                    
                    nrx = nd4x; nry = nd4y;
                    rrx = rd4x; rry = rd4y;
                }
            }
        }
        nrxy = nrx*nry;
        
        System.out.println("output dimensions: "+nrx+" x "+nry+"("+nst+")");
        
        // new image
        System.out.println("deform image");
        float min = 1e10f, max = -1e10f;
        if (padOption.equals("min") || padOption.equals("max")) {
            for (int x=0;x<nsx;x++) for (int y=0;y<nsy;y++) for (int t=0;t<nst;t++) {
                int xyt = x + nsx*y + nsx*nsy*t;
                if (sourceImage[xyt]<min) min = sourceImage[xyt];
                if (sourceImage[xyt]>max) max = sourceImage[xyt];
            }
        }
        deformedImage = new float[nrx*nry*nst];
        for (int x=0;x<nrx;x++) for (int y=0;y<nry;y++) {
            int xy = x + nrx*y;
            
            double dx;
            if (x==0) dx = deformation[xy+1+X*nrxy] - deformation[xy+X*nrxy];
            else if (x==nrx-1) dx = deformation[xy+X*nrxy] - deformation[xy-1+X*nrxy];
            else dx = 0.5*(deformation[xy+1+X*nrxy] - deformation[xy-1+X*nrxy]);
            double dy;
            if (y==0) dy = deformation[xy+nrx+Y*nrxy] - deformation[xy+Y*nrxy];
            else if (y==nry-1) dy = deformation[xy+Y*nrxy] - deformation[xy-nrx+Y*nrxy];
            else dy = 0.5*(deformation[xy+nrx+Y*nrxy] - deformation[xy-nrx+Y*nrxy]);
                
            for (int t=0;t<nst;t++) {
                if (interpOption.equals("nearest")) {
                    if (padOption.equals("closest"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearClosestInterpolation2D(sourceImage, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("zero"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, 0.0f, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("min"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, min, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("max"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, max, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);  
                } else if (interpOption.equals("linear")) {
                    if (padOption.equals("closest"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearClosestInterpolation2D(sourceImage, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("zero"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, 0.0f, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("min"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, min, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                    else if (padOption.equals("max"))
                        deformedImage[xy + nrxy*t] = ImageInterpolation.linearInterpolation2D(sourceImage, max, deformation[xy+X*nrxy], deformation[xy+Y*nrxy], t, nsx, nsy, nst);
                }
                double angle = FastMath.atan2(dy*FastMath.sin(deformedImage[xy + nrxy*t]),dx*FastMath.cos(deformedImage[xy + nrxy*t]));
                deformedImage[xy + nrxy*t] = (float)angle;
            }
        }
        sourceImage = null;
        deformation = null;
    }
}
