package de.mpg.cbs.core.surface;


import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.*;

// Notes: Approximation - distance = min distance to vertices, triangle centroid, edge mid-distance. 
// Could be improved by implementing the projection onto the edges and triangle plane (see JistToolsMeshToLevelset3 and 5)
// Distance function sign calculated using angle weighted pseudo normals, according to the method by: 
// Bærentzen, J A, Aanæs, H, Signed Distance Computation Using the Angle Weighted Pseudonormal, IEEE transactions on visualization and computer graphics 11(3) 243-253 (2005). 

public class SurfaceMeshToLevelsetPseudoNormals {

    // image containers
    private float[] levelsetImage;
	private int nx=-1,ny=-1,nz=-1,nxyz=-1;
	private float rx, ry, rz;

    private float[] pointList;
    private int[] triangleList;
    private float[] inflationScore;
	
	private static final float epsilon = 1e-30f;
	private static final float PADDING = 100.0f;
	private static final int narrowbandDist = 0;
	
	public final void setSurfacePoints(float[] val) { pointList = val; }
	public final void setSurfaceTriangles(int[] val) { triangleList = val; }

	public final void setDimensions(int x, int y, int z) { nx=x; ny=y; nz=z; nxyz=nx*ny*nz; }
	public final void setDimensions(int[] dim) { nx=dim[0]; ny=dim[1]; nz=dim[2]; nxyz=nx*ny*nz; }
	
	public final void setResolutions(float x, float y, float z) { rx=x; ry=y; rz=z; }
	public final void setResolutions(float[] res) { rx=res[0]; ry=res[1]; rz=res[2]; }

	// to be used for JIST definitions, generic info / help
	public final String getPackage() { return "CBS Tools"; }
	public final String getCategory() { return "Surfaces"; }
	public final String getLabel() { return "Mesh To Level Set (with pseudo-normals)"; }
	public final String getName() { return "MeshToLevelSetPseudoNormals"; }

	public final String[] getAlgorithmAuthors() { return new String[]{"Pierre-Louis Bazin","Christine Tardif"}; }
	public final String getAffiliation() { return "Max Planck Institute for Human Cognitive and Brain Sciences"; }
	public final String getDescription() { return "Convert a surface mesh into a levelset"; }
	public final String getLongDescription() { return getDescription(); }
		
	public final String getVersion() { return "3.1.0"; };
	
	public final float[] getLevelset() { return levelsetImage; }
	
	public final void execute() {
		
	    System.out.print("\nSurface Mesh to Levelset Function");
		
		boolean[][][] mask = new boolean[nx][ny][nz];
		float[][][] df = new float[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) {
			df[x][y][z] = PADDING;
			mask[x][y][z] = false;
		}
		
		float[][][] sum = new float[nx][ny][nz];
		
		//System.out.print("Repairing degenerate triangles.\n");
		//surf.repairDegenerateTriangles(0.01f,0.005f);
		int nfaces = triangleList.length/3;
		int npoints = pointList.length/3;
		int nedges = triangleList.length/2;
		
		System.out.print("Number of faces/triangles = "+nfaces+"\n");
		
		System.out.print("Calculate face normals. Looping through triangles...\n");
		float[][] faceNormals = new float[nfaces][3];
		for(int fid=0; fid<nfaces; fid++) {
		    double p0x = pointList[3*triangleList[3*fid+0]+0];
			double p0y = pointList[3*triangleList[3*fid+0]+1];
			double p0z = pointList[3*triangleList[3*fid+0]+2];
			
			double p1x = pointList[3*triangleList[3*fid+1]+0];
			double p1y = pointList[3*triangleList[3*fid+1]+1];
			double p1z = pointList[3*triangleList[3*fid+1]+2];
			
			double p2x = pointList[3*triangleList[3*fid+2]+0];
			double p2y = pointList[3*triangleList[3*fid+2]+1];
			double p2z = pointList[3*triangleList[3*fid+2]+2];
			
			double px = (p1y-p0y)*(p2z-p0z) - (p1z-p0z)*(p2y-p0y);
			double py = (p1z-p0z)*(p2x-p0x) - (p1x-p0x)*(p2z-p0z);
			double pz = (p1x-p0x)*(p2y-p0y) - (p1y-p0y)*(p2x-p0x);
			double norm = FastMath.sqrt( px*px + py*py + pz*pz);
			
			faceNormals[fid][0] = (float)(px/norm);
			faceNormals[fid][1] = (float)(py/norm);
			faceNormals[fid][2] = (float)(pz/norm);
		}
		
		System.out.print("Calculate edge maps...\n");
		int[] edgeList = generateEdgePointList(triangleList);
		int[][] ngbFaceList = generateFaceNeighborTable(npoints, triangleList);
		int[] edgeFaceList = generateEdgeFaceList(edgeList,ngbFaceList);
		int[] faceEdgeList = new int[3*nfaces];
		for (int fid=0;fid<3*nfaces;fid++) faceEdgeList[fid] = -1;
		
		for (int eid=0; eid<nedges; eid++) {
		    if (faceEdgeList[3*edgeFaceList[2*eid+0]+0]==-1) faceEdgeList[3*edgeFaceList[2*eid+0]+0] = eid;
		    else if (faceEdgeList[3*edgeFaceList[2*eid+0]+1]==-1) faceEdgeList[3*edgeFaceList[2*eid+0]+1] = eid;
		    else if (faceEdgeList[3*edgeFaceList[2*eid+0]+2]==-1) faceEdgeList[3*edgeFaceList[2*eid+0]+2] = eid;
		    
		    if (faceEdgeList[3*edgeFaceList[2*eid+1]+0]==-1) faceEdgeList[3*edgeFaceList[2*eid+1]+0] = eid;
		    else if (faceEdgeList[3*edgeFaceList[2*eid+1]+1]==-1) faceEdgeList[3*edgeFaceList[2*eid+1]+1] = eid;
		    else if (faceEdgeList[3*edgeFaceList[2*eid+1]+2]==-1) faceEdgeList[3*edgeFaceList[2*eid+1]+2] = eid;
		}   
		    
		System.out.print("Calculate edge normals. Looping through edges...\n");
		float[][] edgeNormals = new float[nedges][3];
		for(int eid=0; eid<nedges; eid++) {
			int fid1 = edgeFaceList[2*eid+0];
			int fid2 = edgeFaceList[2*eid+1];
			double pex = (faceNormals[fid1][0] + faceNormals[fid2][0]);
			double pey = (faceNormals[fid1][1] + faceNormals[fid2][1]);
			double pez = (faceNormals[fid1][2] + faceNormals[fid2][2]);
			double norm = FastMath.sqrt( pex*pex + pey*pey + pez*pez);
			edgeNormals[eid][0] = (float)(pex/norm);
			edgeNormals[eid][1] = (float)(pey/norm);
			edgeNormals[eid][2] = (float)(pez/norm);
		}
		
		System.out.print("Calculate vertex normals. Looping through vertices...\n");
		float[][] vertexNormals = new float[npoints][3];
		for(int i=0; i<npoints; i++) {
		    double vx = 0.0;
		    double vy = 0.0;
		    double vz = 0.0;
			for(int j=0; j<ngbFaceList[i].length; j++) {
			    int p0 = triangleList[3*ngbFaceList[i][j]+0];
			    int p1 = triangleList[3*ngbFaceList[i][j]+1];
			    int p2 = triangleList[3*ngbFaceList[i][j]+2];
			    double cos = 0.0;
			    if (p0==i) {
			        cos = (pointList[3*p1+0]-pointList[3*p0+0])*(pointList[3*p2+0]-pointList[3*p0+0])
			             +(pointList[3*p1+1]-pointList[3*p0+1])*(pointList[3*p2+1]-pointList[3*p0+1])
			             +(pointList[3*p1+2]-pointList[3*p0+2])*(pointList[3*p2+2]-pointList[3*p0+2]);
			    } else if (p1==i) {
			        cos = (pointList[3*p2+0]-pointList[3*p1+0])*(pointList[3*p0+0]-pointList[3*p1+0])
			             +(pointList[3*p2+1]-pointList[3*p1+1])*(pointList[3*p0+1]-pointList[3*p1+1])
			             +(pointList[3*p2+2]-pointList[3*p1+2])*(pointList[3*p0+2]-pointList[3*p1+2]);
			    } else if (p2==i) {
			        cos = (pointList[3*p0+0]-pointList[3*p2+0])*(pointList[3*p1+0]-pointList[3*p2+0])
			             +(pointList[3*p0+1]-pointList[3*p2+1])*(pointList[3*p1+1]-pointList[3*p2+1])
			             +(pointList[3*p0+2]-pointList[3*p2+2])*(pointList[3*p1+2]-pointList[3*p2+2]);
			    }
			    double angle = FastMath.acos(cos);
			    vx += angle*faceNormals[ngbFaceList[i][j]][0];
				vy += angle*faceNormals[ngbFaceList[i][j]][1];
				vz += angle*faceNormals[ngbFaceList[i][j]][2];
			}
			double norm = FastMath.sqrt(vx*vx + vy*vy + vz*vz);
			vertexNormals[i][0] = (float)(vx/norm);
			vertexNormals[i][1] = (float)(vy/norm);
			vertexNormals[i][2] = (float)(vz/norm);
		}
		
		System.out.print("Calculating distance field. Looping through triangles...\n");
		//For each triangle/face of the mesh, calculate distance of each voxel within a bounding box to plane, vertices and edges
		//fid: face Id		
		
		for(int fid=0; fid<nfaces; fid++) {

			int xmin = Numerics.floor(Numerics.min(pointList[3*triangleList[3*fid+0]+0],pointList[3*triangleList[3*fid+1]+0], pointList[3*triangleList[3*fid+2]+0]))-narrowbandDist;
			int xmax = Numerics.ceil(Numerics.max(pointList[3*triangleList[3*fid+0]+0],pointList[3*triangleList[3*fid+1]+0], pointList[3*triangleList[3*fid+2]+0]))+narrowbandDist;
			int ymin = Numerics.floor(Numerics.min(pointList[3*triangleList[3*fid+0]+1],pointList[3*triangleList[3*fid+1]+1], pointList[3*triangleList[3*fid+2]+1]))-narrowbandDist;
			int ymax = Numerics.ceil(Numerics.max(pointList[3*triangleList[3*fid+0]+1],pointList[3*triangleList[3*fid+1]+1], pointList[3*triangleList[3*fid+2]+1]))+narrowbandDist;
			int zmin = Numerics.floor(Numerics.min(pointList[3*triangleList[3*fid+0]+2],pointList[3*triangleList[3*fid+1]+2], pointList[3*triangleList[3*fid+2]+2]))-narrowbandDist;
			int zmax = Numerics.ceil(Numerics.max(pointList[3*triangleList[3*fid+0]+2],pointList[3*triangleList[3*fid+1]+2], pointList[3*triangleList[3*fid+2]+2]))+narrowbandDist;
			
			//for every voxel in box, calculate min distance to mesh
			//could improve accuracy by calculating the min distance to the face plane instead of distance to face centroid
			for(int x=xmin; x<=xmax; x++) for(int y=ymin; y<=ymax; y++) for(int z=zmin; z<=zmax; z++) {
				//df will be defined at this voxel
				mask[x][y][z] = true;
				
				double dist;
				
				//distance to 3 vertices
				//vertex 0
				int pid = triangleList[3*fid+0];
				dist = FastMath.sqrt((pointList[3*pid+0]-x)*(pointList[3*pid+0]-x)
				                    +(pointList[3*pid+1]-y)*(pointList[3*pid+1]-y)
				                    +(pointList[3*pid+2]-z)*(pointList[3*pid+2]-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (pointList[3*pid+0]-x)*vertexNormals[pid][0]
					              +(pointList[3*pid+1]-y)*vertexNormals[pid][1]
					              +(pointList[3*pid+2]-z)*vertexNormals[pid][2];
					
				}
				
				//vertex 1
				pid = triangleList[3*fid+1];
				dist = FastMath.sqrt((pointList[3*pid+0]-x)*(pointList[3*pid+0]-x)
				                    +(pointList[3*pid+1]-y)*(pointList[3*pid+1]-y)
				                    +(pointList[3*pid+2]-z)*(pointList[3*pid+2]-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (pointList[3*pid+0]-x)*vertexNormals[pid][0]
					              +(pointList[3*pid+1]-y)*vertexNormals[pid][1]
					              +(pointList[3*pid+2]-z)*vertexNormals[pid][2];
					
				}
				
				//vertex 2
				pid = triangleList[3*fid+2];
				dist = FastMath.sqrt((pointList[3*pid+0]-x)*(pointList[3*pid+0]-x)
				                    +(pointList[3*pid+1]-y)*(pointList[3*pid+1]-y)
				                    +(pointList[3*pid+2]-z)*(pointList[3*pid+2]-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (pointList[3*pid+0]-x)*vertexNormals[pid][0]
					              +(pointList[3*pid+1]-y)*vertexNormals[pid][1]
					              +(pointList[3*pid+2]-z)*vertexNormals[pid][2];
					
				}
				
				//distance to face (use face centroid)
				float[] face = getTriangleCenter(fid, pointList, triangleList);
				dist = FastMath.sqrt((face[0]-x)*(face[0]-x)
				                    +(face[1]-y)*(face[1]-y)
				                    +(face[2]-z)*(face[2]-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (face[0]-x)*faceNormals[fid][0]
					              +(face[1]-y)*faceNormals[fid][1]
					              +(face[2]-z)*faceNormals[fid][2];
					
				} 
								
				//distance to 3 edges (use edge midpoint)
				//edge 0
				int eid = faceEdgeList[3*fid+0];
				int pid1 = edgeList[2*eid+0];
				int pid2 = edgeList[2*eid+1];
				dist = FastMath.sqrt((0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*(0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)
				                    +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)
				                    +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*edgeNormals[eid][0]
					              +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*edgeNormals[eid][1]
					              +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*edgeNormals[eid][2];
					
				} 
				
				//edge 1
				eid = faceEdgeList[3*fid+1];
				pid1 = edgeList[2*eid+0];
				pid2 = edgeList[2*eid+1];
				dist = FastMath.sqrt((0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*(0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)
				                    +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)
				                    +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*edgeNormals[eid][0]
					              +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*edgeNormals[eid][1]
					              +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*edgeNormals[eid][2];
					
				} 
				
				//edge 2
				eid = faceEdgeList[3*fid+2];
				pid1 = edgeList[2*eid+0];
				pid2 = edgeList[2*eid+1];
				dist = FastMath.sqrt((0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*(0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)
				                    +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)
				                    +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z));
				
				if (dist<df[x][y][z]) {
					df[x][y][z] = (float)dist;
					sum[x][y][z] = (0.5f*(pointList[3*pid1+0]+pointList[3*pid2+0])-x)*edgeNormals[eid][0]
					              +(0.5f*(pointList[3*pid1+1]+pointList[3*pid2+1])-y)*edgeNormals[eid][1]
					              +(0.5f*(pointList[3*pid1+2]+pointList[3*pid2+2])-z)*edgeNormals[eid][2];
					
				} 
				
			}
		
		}

		System.out.print("done.\n");
		System.out.print("Adding sign to distance field.\n");
		//add sign to distance field
		float[] sdf = new float[nx*ny*nz];
		boolean[] bgmask = new boolean[nx*ny*nz];
		boolean[] datamask = new boolean[nx*ny*nz];
		
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) {
			int xyz = x+nx*y+nx*ny*z;
			bgmask[xyz] = true;
			
			if(mask[x][y][z]) {
				datamask[xyz] = true;
				sdf[xyz] = FastMath.copySign(df[x][y][z], sum[x][y][z]);
			} else {
				datamask[xyz] = false;
			}
		}
		
		sdf = cleanupSign(sdf, datamask, 26, nx, ny, nz);
		
		System.out.print("Padding levelset.\n");
		sdf = padLevelset(sdf, datamask, 6, nx, ny, nz);
		//sdf = padLevelset(sdf, datamask, 26, nx, ny, nz);
		//sdf = padLevelset(sdf, datamask, 18, nx, ny, nz);
		
		for (int x=3; x<nx-3; x++) for (int y=3; y<ny-3; y++) for (int z=3; z<nz-3; z++) {
			int xyz = x+nx*y+nx*ny*z;
			datamask[xyz] = true;
			//df[x][y][z] = sdf[xyz];
		}
		
		System.out.print("Removing sign errors.\n");
		sdf = cleanupPaddingSign(sdf, datamask, 26, nx, ny, nz);
		
		System.out.print("Clean up sign 2.\n");
		sdf = cleanupSign2(sdf, 26, nx, ny, nz);
		
		System.out.print("Evolving narrow band.\n");
		InflateGdm gdm = new InflateGdm(sdf, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		gdm.evolveNarrowBand(0, 1.0f);
		
		System.out.print("Output.\n");
		
		float[][][] lvl = gdm.exportLevelset();		
		levelsetImage = new float[nxyz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
		    levelsetImage[x+nx*y+nx*ny*z] = lvl[x][y][z];
		}
	}
	
	private final float[] cleanupSign(float[] sdf, boolean[] datamask, int connectivity, int nx, int ny, int nz) {
		float[] cleansdf = new float[nx*ny*nz];
		
		int[] xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
		int[] yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
		int[] zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny};
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			float meanSdf = 0.0f;
			int count = 0;
			if (datamask[xyz]) {
				for (int k=0; k<connectivity; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (datamask[xyzn]) {
						meanSdf += sdf[xyzn];
						count++;
					}
				}
				meanSdf/= count;
				 
				if (FastMath.abs(sdf[xyz]-meanSdf) > 1.732) { //1.732
					cleansdf[xyz] = -sdf[xyz];
				} else {
					cleansdf[xyz] = sdf[xyz];
				}
			} else {
				cleansdf[xyz] = sdf[xyz];
			}
		}
		return cleansdf;
			
	}
	
	private final float[] cleanupPaddingSign(float[] sdf, boolean[] datamask, int connectivity, int nx, int ny, int nz) {
		float[] cleansdf = new float[nx*ny*nz];
		
		int[] xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
		int[] yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
		int[] zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny};
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;
			
			cleansdf[xyz] = sdf[xyz];
			if (datamask[xyz]) {
				if (Numerics.abs(sdf[xyz])==PADDING) {
					//System.out.print(".");
					// count the number of padded values with each sign
					int pos=0;
					int neg=0;
					for (int k=0; k<connectivity; k++) {
						int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
						if (datamask[xyzn]) {
							if (sdf[xyzn]>=0) {
								pos++;
							} else if (sdf[xyzn]<0) {
								neg++;
							}
						}
					}
					if (pos>neg) cleansdf[xyz] = PADDING;
					else cleansdf[xyz] = -PADDING;
					if (cleansdf[xyz]!=sdf[xyz]) System.out.print(".");
				}
			}
		}
		return cleansdf;
			
	}
	
	private final float[] cleanupSign2(float[] sdf, int connectivity, int nx, int ny, int nz) {
		float[] cleansdf = new float[nx*ny*nz];
		boolean change = true;
		int itrs = 0;
		int count = 0;
		int[] xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
		int[] yoff = new int[]{0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1};
		int[] zoff = new int[]{0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1,-1, 1, 1, -1,-1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
		
		// if neighbour is 100 (PADDING value), can't have a negative sign
		for (int i=0;i<5;i++) {
			System.out.print("Iteration "+itrs);
			itrs++;
			count = 0;
			change = false;
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				
				if (sdf[xyz]<0) {
					for (int k=0; k<connectivity; k++) {
						if (x+xoff[k]>=0 && x+xoff[k]<nx && y+yoff[k]>=0 && y+yoff[k]<ny && z+zoff[k]>=0 && z+zoff[k]<nz) {
							int xyzn = x+xoff[k] + (y+yoff[k])*nx + (z+zoff[k])*nx*ny;
							if (sdf[xyzn]== PADDING) {
								cleansdf[xyz] = FastMath.abs(sdf[xyz]);
								change = true;
								count++;
							} else {
								cleansdf[xyz] = sdf[xyz];
							}
						}
					}		 
		
				} else {
					cleansdf[xyz] = sdf[xyz];
				}
			}
			System.out.print(" Count = "+count+"\n");
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				sdf[xyz] = cleansdf[xyz];
			}
		}
		
		return cleansdf;
			
	}
		
	private final float[] padLevelset(float[] data, boolean[] datamask, int connectivity, int nx, int ny, int nz) {
		System.out.print("Pad levelset.\n");
		
		float[] dilData = new float[nx*ny*nz];
		
		int[] xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
		int[] yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
		int[] zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny};
		
		BitSet boundary = new BitSet(nx*ny*nz);
        BitSet nextboundary = new BitSet(nx*ny*nz);
        BitSet used = new BitSet(nx*ny*nz);
        BitSet nextused = new BitSet(nx*ny*nz);
        
        System.out.print("\n Set initial processed data");
        //set initial used
        for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			if ((x<2 || x>=nx-2) || (y<2 || y>=ny-2) || (z<2 || z>=nz-2)) {
				dilData[xyz] = PADDING;
				used.set(xyz);
			} else {
				if (datamask[xyz]) {
					used.set(xyz);
					dilData[xyz] = data[xyz];
				}
			}
        }
       	
        System.out.println("\n Set initial boundary");
        //set initial boundary
        for (int x=2;x<nx-2;x++) for (int y=2;y<ny-2;y++) for (int z=2;z<nz-2;z++) {
			int xyz = x+nx*y+nx*ny*z;	
			
			if (used.get(xyz)) { 
				for (int k=0; k<connectivity; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
					if (!used.get(xyzn) && !boundary.get(xyzn)) {
						boundary.set(xyzn);
					}
				}
			}
        }
        
		//dilate boundary
		boolean dilate = true;
		
        nextused = (BitSet)used.clone();
        
    	System.out.println("\n Dilate");
               
        while (dilate) {
        	//System.out.println(".");
        	dilate = false;
			for (int xyz = boundary.nextSetBit(0); xyz >= 0; xyz = boundary.nextSetBit(xyz+1)) {
				// dilate and build the next boundary
				
				for (int k=0; k<connectivity; k++) {
					int xyzn = xyz + xoff[k] + yoff[k] + zoff[k];
										
					if (used.get(xyzn)) {
						if (dilData[xyzn]<0.0f) dilData[xyz] = -PADDING;
						else if (dilData[xyzn]>=0.0f) dilData[xyz] = PADDING;
					} else {
						nextboundary.set(xyzn);
						dilate = true;
					}
				}
			
				// add to the processed list
				nextused.set(xyz);
				
			}
			
			// replace active boundary
			boundary = (BitSet) nextboundary.clone();
			nextboundary.clear();

			// replace active processed list
			used = (BitSet) nextused.clone();
			
        }
        
        return dilData;
	}
	
	private final int[] generateEdgePointList(int[] faceList) {
	    int nedges = faceList.length/2;
	    
	    int[] edgeList = new int[2*nedges];
	    int ne=0;
	    for (int f=0;f<faceList.length;f+=3) {
	        // only add the edges from low to high labels (arbitrary):
	        // that way edges are counted only once (triangles are oriented)
	        if (faceList[f+0]<faceList[f+1]) {
                edgeList[ne+0] =  faceList[f+0];
                edgeList[ne+1] =  faceList[f+1];
                ne += 2;
            }
            
	        if (faceList[f+1]<faceList[f+2]) {
                edgeList[ne+0] =  faceList[f+1];
                edgeList[ne+1] =  faceList[f+2];
                ne += 2;
            }
	        if (faceList[f+2]<faceList[f+0]) {
                edgeList[ne+0] =  faceList[f+2];
                edgeList[ne+1] =  faceList[f+0];
                ne += 2;
            }
        }
        return edgeList;
	}
	
	private final int[] generateEdgeFaceList(int[] edgeList, int[][] pointFaceList) {
	    int nedges = edgeList.length;
	    
        int[] edgeFaceList = new int[2*nedges];
        for (int e=0;e<nedges;e++) {
            int p0 = edgeList[2*e+0];
            int p1 = edgeList[2*e+1];
            
            int de=0;
            for (int f0=0;f0<pointFaceList[p0].length;f0++) {
                for (int f1=0;f1<pointFaceList[p1].length;f1++) {
                    if (pointFaceList[p0][f0]==pointFaceList[p1][f1]) {
                        edgeFaceList[2*e+de] = pointFaceList[p0][f0];
                        de++;
                    }
                }
            }
            if (de!=2) System.out.print("!");
	    }
	    return edgeFaceList;   
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

}
