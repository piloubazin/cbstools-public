package de.mpg.cbs.core.surface;


import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.*;

/**
 * Inflate surface to specified mean curvature using Lorentzian norm as
 * termination criterion. A weighted umbrella operator is used to regularize the
 * surface. This method is useful for generating surfaces that are at multiple
 * levels of geometric complexity, although cusp formation is common with this
 * inflation procedure.
 * 
 * adapted from the surface inflation code
 * from Duygu Tosun and Blake Lucas
 * 
 */
public class SurfaceInflation {

    // image containers
    private float[] pointList;
    private int[] triangleList;
    private float[] inflationScore;
    
    private float   stepSize = 0.75f;       // BETA step size
    private int     step = 10;               // number of iterations at which to rescale surface to preserve volume
    private float   max_curv = 10.0f;       // EPS curvature threshold for termintion
    private int     max_iter = 2000;          // imax maximum number of iterations
    //private boolean useLorentzian = true;   // if Lorentzian scaling should be used
   
    public static final int AREA=1;       // use area to weight neighbor contributions
    public static final int DIST=2;       // use distance to weight neighbor contributions
    public static final int NUMV=3;       // don't weight neighbor contributions
    private int    weighting = AREA;
    private float   regularization = 1e-3f;
    private float   centering = 1e-3f;
    
    public final void setWeightingMethod(int val) { weighting = val; }
    public final void setRegularization(float val) { regularization = val; }
    public final void setCentering(float val) { centering = val; }
    
	public final void setSurfacePoints(float[] val) { pointList = val; }
	public final void setSurfaceTriangles(int[] val) { triangleList = val; }

	public final float[] 	getInflatedSurfacePoints() { return pointList; }
	public final int[] 		getInflatedSurfaceTriangles() { return triangleList; }
	public final float[] 	getInflatedSurfaceValues() { return inflationScore; }

	public final void setStepSize(float val) { stepSize = val; }
	public final void setMaxIter(int val) { max_iter = val; }
	public final void setMaxCurv(float val) { max_curv = val; }
	//public final void setUseLorentzian(boolean val) { useLorentzian = val; }
	
    private int[][] generatePointNeighborTable(int npt, int[] faces) {
        ArrayList[] table = new ArrayList[npt];
        for (int n=0;n<npt;n++) table[n] = new ArrayList<Integer>();
        for (int f=0;f<faces.length;f+=3) {
            table[faces[f+0]].add(faces[f+1]);
            //table[faces[f+0]].add(faces[f+2])
            
            table[faces[f+1]].add(faces[f+2]);
            //table[faces[f+1]].add(faces[f+0])
            
            table[faces[f+2]].add(faces[f+0]);
            //table[faces[f+2]].add(faces[f+1])
        }
        int[][] map = new int[npt][];
        for (int n=0;n<npt;n++) {
            map[n] = new int[table[n].size()];
            for (int m=0;m<table[n].size();m++)
                map[n][m] = (int)table[n].get(m);
        }
        return map;
    }
    
    private int[][] generateFaceNeighborTable(int npt, int[] faces) {
        ArrayList[] table = new ArrayList[npt];
        for (int n=0;n<npt;n++) table[n] = new ArrayList();
        for (int f=0;f<faces.length;f+=3) {
            table[faces[f+0]].add(f/3);
            table[faces[f+1]].add(f/3);
            table[faces[f+2]].add(f/3);
        }
        int[][] map = new int[npt][];
        for (int n=0;n<npt;n++) {
            map[n] = new int[table[n].size()];
            for (int m=0;m<table[n].size();m++)
                map[n][m] = (int)table[n].get(m);
        }
        return map;
    }
    
    private float getTriangleArea(int id, float[] pts, int[] faces) {
        double vx, vy, vz;
        
        float V0x = pts[3*faces[3*id+0]+0];
        float V0y = pts[3*faces[3*id+0]+1];
        float V0z = pts[3*faces[3*id+0]+2];
        
        float V1x = pts[3*faces[3*id+1]+0];
        float V1y = pts[3*faces[3*id+1]+1];
        float V1z = pts[3*faces[3*id+1]+2];
        
        float V2x = pts[3*faces[3*id+2]+0];
        float V2y = pts[3*faces[3*id+2]+1];
        float V2z = pts[3*faces[3*id+2]+2];
        
        vx = ((V0y * V1z) - (V1y * V0z) + (V1y * V2z) - (V2y * V1z) + (V2y * V0z) - (V0y * V2z));
        vy = ((V0z * V1x) - (V1z * V0x) + (V1z * V2x) - (V2z * V1x) + (V2z * V0x) - (V0z * V2x));
        vz = ((V0x * V1y) - (V1x * V0y) + (V1x * V2y) - (V2x * V1y) + (V2x * V0y) - (V0x * V2y));
        
        return (float) (0.5 * Math.sqrt((vx * vx) + (vy * vy) + (vz * vz)));
    }
    
    private float[] getTriangleCenter(int id, float[] pts, int[] faces) {
        float[] center = new float[3];
        center[0] += pts[3*faces[3*id+0]+0]/3.0f;
        center[1] += pts[3*faces[3*id+0]+1]/3.0f;
        center[2] += pts[3*faces[3*id+0]+2]/3.0f;
        
        center[0] += pts[3*faces[3*id+1]+0]/3.0f;
        center[1] += pts[3*faces[3*id+1]+1]/3.0f;
        center[2] += pts[3*faces[3*id+1]+2]/3.0f;
        
        center[0] += pts[3*faces[3*id+2]+0]/3.0f;
        center[1] += pts[3*faces[3*id+2]+1]/3.0f;
        center[2] += pts[3*faces[3*id+2]+2]/3.0f;
        
        return center;
    }
    
    private float[] getTriangleIncenter(int id, float[] pts, int[] faces) {
        float[] length = new float[3];
        float[] incenter = new float[3];
        length[0] = centering+(float)FastMath.sqrt(Numerics.square(pts[3*faces[3*id+1]+0]-pts[3*faces[3*id+2]+0])
                                                       +Numerics.square(pts[3*faces[3*id+1]+1]-pts[3*faces[3*id+2]+1])
                                                       +Numerics.square(pts[3*faces[3*id+1]+2]-pts[3*faces[3*id+2]+2]));

        length[1] = centering+(float)FastMath.sqrt(Numerics.square(pts[3*faces[3*id+2]+0]-pts[3*faces[3*id+0]+0])
                                                       +Numerics.square(pts[3*faces[3*id+2]+1]-pts[3*faces[3*id+0]+1])
                                                       +Numerics.square(pts[3*faces[3*id+2]+2]-pts[3*faces[3*id+0]+2]));

        length[2] = centering+(float)FastMath.sqrt(Numerics.square(pts[3*faces[3*id+0]+0]-pts[3*faces[3*id+1]+0])
                                                       +Numerics.square(pts[3*faces[3*id+0]+1]-pts[3*faces[3*id+1]+1])
                                                       +Numerics.square(pts[3*faces[3*id+0]+2]-pts[3*faces[3*id+1]+2]));
        
        float sum = length[0]+length[1]+length[2];

        incenter[0] += length[0]*pts[3*faces[3*id+0]+0]/sum;
        incenter[1] += length[0]*pts[3*faces[3*id+0]+1]/sum;
        incenter[2] += length[0]*pts[3*faces[3*id+0]+2]/sum;
        
        incenter[0] += length[1]*pts[3*faces[3*id+1]+0]/sum;
        incenter[1] += length[1]*pts[3*faces[3*id+1]+1]/sum;
        incenter[2] += length[1]*pts[3*faces[3*id+1]+2]/sum;
        
        incenter[0] += length[2]*pts[3*faces[3*id+2]+0]/sum;
        incenter[1] += length[2]*pts[3*faces[3*id+2]+1]/sum;
        incenter[2] += length[2]*pts[3*faces[3*id+2]+2]/sum;
        
        return incenter;
    }
    
    private float[] getCenterOfMass(float[] pts) {
        float[] center = new float[3];
        float npt = pts.length/3.0f;
        for (int n=0;n<pts.length;n+=3) {
            center[0] += pts[n+0]/npt;
            center[1] += pts[n+1]/npt;
            center[2] += pts[n+2]/npt;
        }
        return center;
    }
        
    private double getVolume(float[] pts, int[] faces) {
        double volume = 0.0;
        for (int f=0;f<faces.length;f+=3) {
            volume += pts[3*faces[f+0]+0]*(pts[3*faces[f+1]+1]*pts[3*faces[f+2]+2]-pts[3*faces[f+1]+2]*pts[3*faces[f+2]+1]);
            volume += pts[3*faces[f+0]+1]*(pts[3*faces[f+1]+2]*pts[3*faces[f+2]+0]-pts[3*faces[f+1]+0]*pts[3*faces[f+2]+2]);
            volume += pts[3*faces[f+0]+2]*(pts[3*faces[f+1]+0]*pts[3*faces[f+2]+1]-pts[3*faces[f+1]+1]*pts[3*faces[f+2]+0]);
        }
        return volume;
    }
        
    public void execute() {

        double BETA = stepSize;
        float EPS = max_curv;
        int imax = max_iter;
		float scaleL = 0;
							
		int npt = pointList.length/3;
		int nface = triangleList.length/3;
		
		int[][] ngbPoints = generatePointNeighborTable(npt, triangleList);					
		int[][] ngbFaces = generateFaceNeighborTable(npt, triangleList);					
		
		System.out.printf("mean curvature measure threshold = %.2f\n", EPS);
		double areaO = 0;
		for (int i = 0; i < nface; ++i) {
			areaO += getTriangleArea(i, pointList, triangleList);
		}
        scaleL = (float) updateScaleL(pointList, triangleList, ngbPoints, ngbFaces, npt);
		System.out.printf("Lorentzian scale factor = %g\n", scaleL);
		double value = FastMath.sqrt(FastMath.log(1 + 0.5 / (scaleL * scaleL)));
		EPS *= value;
					
		System.out.printf("\n");
		double hCurv, hCurv_0, hCurvOld;
		double origVolume = getVolume(pointList, triangleList);
		// create nbd list for each vertex --- won't be changed
		float[] oldpt = new float[3*npt];
        for (int i = 0; i < 3*npt; i++) {
            oldpt[i]  = pointList[i];
        }
        
		hCurv = meanCurv(pointList, triangleList, ngbPoints, ngbFaces, npt, scaleL);
		hCurv_0 = hCurv;
		hCurvOld = hCurv;
		double hCurvFirst = -1;
		int iter = 0;
		System.out.println("initial mean curvature  = "+hCurv_0);
		do {
		    //float[] newpt = new float[3*npt];
			for (int i = 0; i < npt; i++) {
			    double ptx = 0.0, pty = 0.0, ptz = 0.0;
			    double totalA = 0.0;
				
			    // Regularize mesh by small amount BETA by moving point to
				// center
				for (int k = 0; k < ngbFaces[i].length; k++) {
					int face = ngbFaces[i][k];
					// Calculate center of surrounding region
					//float[] C = getTriangleCenter(face, pointList, triangleList);
					float[] C = getTriangleIncenter(face, pointList, triangleList);
					float weight = 1.0f;
					if (weighting==AREA)
                        weight = regularization+getTriangleArea(face, pointList, triangleList);
                    else if (weighting==DIST)
                        weight = regularization+
                                    (float)FastMath.sqrt( (pointList[3*i+0]-C[0])*(pointList[3*i+0]-C[0])
                                                         +(pointList[3*i+1]-C[1])*(pointList[3*i+1]-C[1])
                                                         +(pointList[3*i+2]-C[2])*(pointList[3*i+2]-C[2]) );
                                               
					totalA += weight;

					ptx += weight * C[0];
					pty += weight * C[1];
					ptz += weight * C[2];
				}
				if (totalA > 0) {
					ptx /= totalA;
					pty /= totalA;
					ptz /= totalA;
				}
				//newpt[3*i+0] = (float) (pointList[3*i+0] * (1 - BETA) + BETA * ptx);
				//newpt[3*i+1] = (float) (pointList[3*i+1] * (1 - BETA) + BETA * pty);
				//newpt[3*i+2] = (float) (pointList[3*i+2] * (1 - BETA) + BETA * ptz);
				pointList[3*i+0] = (float) (pointList[3*i+0] * (1 - BETA) + BETA * ptx);
				pointList[3*i+1] = (float) (pointList[3*i+1] * (1 - BETA) + BETA * pty);
				pointList[3*i+2] = (float) (pointList[3*i+2] * (1 - BETA) + BETA * ptz);
				// System.out.println("POINT "+i+" "+V);
			}
			/*
			for (int i = 0; i < npt; i++) {
			    pointList[3*i+0] = newpt[3*i+0];
			    pointList[3*i+1] = newpt[3*i+1];
			    pointList[3*i+2] = newpt[3*i+2];
			}*/
			if (iter % step == 0) {
			    
				// Rescale vertices to preserve volume
				double scale = FastMath.pow((origVolume / getVolume(pointList,triangleList)), 0.3333333333);
				float[] centroid = getCenterOfMass(pointList);
				for (int i = 0; i < npt; ++i) {
				    pointList[3*i+0] = (float) ((pointList[3*i+0] - centroid[0]) * scale + centroid[0]);
					pointList[3*i+1] = (float) ((pointList[3*i+1] - centroid[1]) * scale + centroid[1]);
					pointList[3*i+2] = (float) ((pointList[3*i+2] - centroid[2]) * scale + centroid[2]);
				}
				
				// Calculate Total Mean Curvature
				hCurv = meanCurv(pointList, triangleList, ngbPoints, ngbFaces, npt, scaleL);
				float max = (float) (Math.abs(hCurv - hCurvOld) / hCurv_0);
				hCurvOld = hCurv;

				//System.out.println("Completed: "+(1 - (hCurv - EPS) / (hCurvFirst - EPS)));
				System.out.println(iter+". Mean Curvature " + hCurv);
			}
			iter++;
		} while ((hCurv > EPS) && (iter <= imax));

		// compute an "inflation score: relative distance change to centroid
		float[] centroid = getCenterOfMass(pointList);
		
		inflationScore = new float[npt];
		for (int i=0;i<npt;i++) {
		    inflationScore[i] = (float)(
		                        FastMath.sqrt((oldpt[3*i+0]-centroid[0])*(oldpt[3*i+0]-centroid[0])
		                                     +(oldpt[3*i+1]-centroid[1])*(oldpt[3*i+1]-centroid[1])
		                                     +(oldpt[3*i+2]-centroid[2])*(oldpt[3*i+2]-centroid[2]))
		                       -FastMath.sqrt((pointList[3*i+0]-centroid[0])*(pointList[3*i+0]-centroid[0])
		                                     +(pointList[3*i+1]-centroid[1])*(pointList[3*i+1]-centroid[1])
		                                     +(pointList[3*i+2]-centroid[2])*(pointList[3*i+2]-centroid[2]))
		                       );
		}
		return;
	}

	/********************************************************/
	private double meanCurv(float[] pointList, int[] faceList, int[][] ngbPoints, int[][] ngbFaces, int npt, double scaleL) {

		double curv = 0;
		double totalA = 0.0;
		for (int i = 0; i < npt; ++i) {
		    // compute the whole area: the 1-ring neighborhood is simply 1/3
		    double area = 0.0;
		    for (int j=0;j<ngbPoints[i].length;j++) {
		        area += getTriangleArea(ngbFaces[i][j], pointList, faceList);
		    }
		    area /= 3.0;
		    totalA += area;
		    
		    double vx = 0.0, vy = 0.0, vz = 0.0;
		    
		    for (int j=0;j<ngbPoints[i].length;j++) {
                // for each vertex, compute the local curvature with the next neighbor
                int ngb = ngbPoints[i][j];
                int ngbFace1 = ngbFaces[i][j];
                // find the corresponding in the other side
                int ngbFace2 = -1;
                for (int k=0;k<ngbPoints[ngb].length;k++) {
                    if (ngbPoints[ngb][k]==i) {
                        ngbFace2 = ngbFaces[ngb][k];
                    }
                }
                if (ngbFace2==-1) System.out.print("!");
                
                // get the cotangent of the outside angles
                int third = -1;
                if (faceList[3*ngbFace1+0]!=i && faceList[3*ngbFace1+0]!=ngb) third = faceList[3*ngbFace1+0];
                else if (faceList[3*ngbFace1+1]!=i && faceList[3*ngbFace1+1]!=ngb) third = faceList[3*ngbFace1+1];
                else if (faceList[3*ngbFace1+2]!=i && faceList[3*ngbFace1+2]!=ngb) third = faceList[3*ngbFace1+2];
                
                float vax = pointList[3*i+0]-pointList[3*third+0];
                float vay = pointList[3*i+1]-pointList[3*third+1];
                float vaz = pointList[3*i+2]-pointList[3*third+2];
                
                float vbx = pointList[3*ngb+0]-pointList[3*third+0];
                float vby = pointList[3*ngb+1]-pointList[3*third+1];
                float vbz = pointList[3*ngb+2]-pointList[3*third+2];
                
                double alpha = FastMath.acos( vax*vbx + vay*vby + vaz*vbz );
                
                third = -1;
                if (faceList[3*ngbFace2+0]!=i && faceList[3*ngbFace2+0]!=ngb) third = faceList[3*ngbFace2+0];
                else if (faceList[3*ngbFace2+1]!=i && faceList[3*ngbFace2+1]!=ngb) third = faceList[3*ngbFace2+1];
                else if (faceList[3*ngbFace2+2]!=i && faceList[3*ngbFace2+2]!=ngb) third = faceList[3*ngbFace2+2];
                
                vax = pointList[3*i+0]-pointList[3*third+0];
                vay = pointList[3*i+1]-pointList[3*third+1];
                vaz = pointList[3*i+2]-pointList[3*third+2];
                
                vbx = pointList[3*ngb+0]-pointList[3*third+0];
                vby = pointList[3*ngb+1]-pointList[3*third+1];
                vbz = pointList[3*ngb+2]-pointList[3*third+2];
                
                double gamma = FastMath.acos( vax*vbx + vay*vby + vaz*vbz );
                
                double factor = 1.0/FastMath.tan(alpha) + 1.0/FastMath.tan(gamma);
                if (Double.isNaN(factor) || Double.isInfinite(factor)) factor = 0.0;
                
                vx += factor*(pointList[3*i+0]-pointList[3*ngb+0]);
                vy += factor*(pointList[3*i+1]-pointList[3*ngb+1]);
                vz += factor*(pointList[3*i+2]-pointList[3*ngb+2]);
            }
                
		    //curv += ( vx*vx + vy*vy + vz*vz )/(2.0*area);
		    curv += FastMath.log(1.0 + 0.5*( vx*vx + vy*vy + vz*vz ) / (scaleL*scaleL))*area;
		    /*
            if (useLorentzian)
                curv += Math.log(1 + 0.5 * ( vx*vx + vy*vy + vz*vz )
                        / (scaleL * scaleL))
                        * area;
            else
                curv += ( vx*vx + vy*vy + vz*vz )/(2.0*area);
            */  
        }
        curv = FastMath.sqrt(curvNorm*curv);
		// curv /= area;
		//curv *= curvNorm;
		//curv = FastMath.sqrt(curv);
		return (curv);
	}

	private static final float curvNorm = (float) (1.0f / (4 * Math.acos(-1)));

	private double updateScaleL(float[] pointList, int[] faceList, int[][] ngbPoints, int[][] ngbFaces, int npt) {
		
	    double[] kcurv = new double[npt];
		
		for (int i = 0; i < npt; ++i) {
		    // compute the whole area: the 1-ring neighborhood is simply 1/3
		    double area = 0.0;
		    for (int j=0;j<ngbPoints[i].length;j++) {
		        area += getTriangleArea(ngbFaces[i][j], pointList, faceList);
		    }
		    area /= 3.0;
		    
		    double vx = 0.0, vy = 0.0, vz = 0.0;
		    
		    for (int j=0;j<ngbPoints[i].length;j++) {
                // for each vertex, compute the local curvature with the next neighbor
                int ngb = ngbPoints[i][j];
                int ngbFace1 = ngbFaces[i][j];
                // find the corresponding in the other side
                int ngbFace2 = -1;
                for (int k=0;k<ngbPoints[ngb].length;k++) {
                    if (ngbPoints[ngb][k]==i) {
                        ngbFace2 = ngbFaces[ngb][k];
                    }
                }
                if (ngbFace2==-1) System.out.print("!");
                
                // get the cotangent of the outside angles
                int third = -1;
                if (faceList[3*ngbFace1+0]!=i && faceList[3*ngbFace1+0]!=ngb) third = faceList[3*ngbFace1+0];
                else if (faceList[3*ngbFace1+1]!=i && faceList[3*ngbFace1+1]!=ngb) third = faceList[3*ngbFace1+1];
                else if (faceList[3*ngbFace1+2]!=i && faceList[3*ngbFace1+2]!=ngb) third = faceList[3*ngbFace1+2];
                
                float vax = pointList[3*i+0]-pointList[3*third+0];
                float vay = pointList[3*i+1]-pointList[3*third+1];
                float vaz = pointList[3*i+2]-pointList[3*third+2];
                
                float vbx = pointList[3*ngb+0]-pointList[3*third+0];
                float vby = pointList[3*ngb+1]-pointList[3*third+1];
                float vbz = pointList[3*ngb+2]-pointList[3*third+2];
                
                double alpha = FastMath.acos( vax*vbx + vay*vby + vaz*vbz );
                
                third = -1;
                if (faceList[3*ngbFace2+0]!=i && faceList[3*ngbFace2+0]!=ngb) third = faceList[3*ngbFace2+0];
                else if (faceList[3*ngbFace2+1]!=i && faceList[3*ngbFace2+1]!=ngb) third = faceList[3*ngbFace2+1];
                else if (faceList[3*ngbFace2+2]!=i && faceList[3*ngbFace2+2]!=ngb) third = faceList[3*ngbFace2+2];
                
                vax = pointList[3*i+0]-pointList[3*third+0];
                vay = pointList[3*i+1]-pointList[3*third+1];
                vaz = pointList[3*i+2]-pointList[3*third+2];
                
                vbx = pointList[3*ngb+0]-pointList[3*third+0];
                vby = pointList[3*ngb+1]-pointList[3*third+1];
                vbz = pointList[3*ngb+2]-pointList[3*third+2];
                
                double gamma = FastMath.acos( vax*vbx + vay*vby + vaz*vbz );
                
                double factor = 1.0/FastMath.tan(alpha) + 1.0/FastMath.tan(gamma);
                if (Double.isNaN(factor) || Double.isInfinite(factor)) factor = 0.0;
                
                vx += factor*(pointList[3*i+0]-pointList[3*ngb+0]);
                vy += factor*(pointList[3*i+1]-pointList[3*ngb+1]);
                vz += factor*(pointList[3*i+2]-pointList[3*ngb+2]);
            }
                
		    kcurv[i] = FastMath.sqrt( vx*vx + vy*vy + vz*vz )/(2.0*area); 
		}
		Percentile measure = new Percentile();
        double median = measure.evaluate(kcurv, 50.0);
		// printf("median = %f ---> ",median);
		for (int i = 0; i < npt; ++i)
			kcurv[i] = FastMath.abs(kcurv[i] - median);
		median = measure.evaluate(kcurv, 50.0);
		// printf("%f\n",median);
		// if (median > 1)
		// Magic number for Median Absolute Deviation MAD
		double scaleL = MedianAbsoluteDeviation * median;
		// else
		// scaleL = 1.4826;
		// printf("scaleL = %f\n",scaleL);
		return scaleL;
	}

	private static final float MedianAbsoluteDeviation = 1.4826f;


}
