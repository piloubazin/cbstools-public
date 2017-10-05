package de.mpg.cbs.jist.surface;

import javax.vecmath.Point3f;
import javax.vecmath.Tuple3f;
import javax.vecmath.Vector3f;

import java.io.*;
import java.util.*;
import gov.nih.mipav.view.*;

import gov.nih.mipav.model.structures.jama.*;

import org.apache.commons.math3.util.FastMath;

import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;
import edu.jhu.ece.iacl.algorithms.graphics.intersector.SurfaceIntersector;
import edu.jhu.ece.iacl.algorithms.graphics.isosurf.IsoSurfaceOnGrid;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.ResampleLevelSet;
import edu.jhu.ece.iacl.algorithms.graphics.utilities.SurfaceToMask;
import edu.jhu.ece.iacl.algorithms.topology.ConnectivityRule;
import edu.jhu.ece.iacl.algorithms.topology.TopologyCorrection;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;
import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataInt;
import edu.jhu.ece.iacl.algorithms.volume.DistanceField;

import edu.jhu.cs.cisst.algorithms.geometry.surface.*;
import edu.jhu.cs.cisst.algorithms.segmentation.gac.DistanceField3D;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

// Notes: Approximation - distance = min distance to vertices, triangle centroid, edge mid-distance. 
// Could be improved by implementing the projection onto the edges and triangle plane (see JistToolsMeshToLevelset3 and 5)
// Distance function sign calculated using angle weighted pseudo normals, according to the method by: 
// Bærentzen, J A, Aanæs, H, Signed Distance Computation Using the Angle Weighted Pseudonormal, IEEE transactions on visualization and computer graphics 11(3) 243-253 (2005). 

public class JistSurfaceMeshToLevelsetPseudoNormals  extends ProcessingAlgorithm {
	ParamVolume origVol;
	ParamSurface origSurf;
	ParamVolume result;
	//ParamVolume result2;
	ParamBoolean mipavTransform;
	
	private static final float epsilon = 1e-30f;
	private static final float PADDING = 100.0f;
	private static final int narrowbandDist = 0;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(origSurf=new ParamSurface("Surface"));
		inputParams.add(origVol=new ParamVolume("Reference Volume"));
		origVol.setDescription("Reference volume to use for level set representation dimensions.");
		inputParams.add(mipavTransform=new ParamBoolean("Align to MIPAV image space", true));
		
		inputParams.setName("MeshToLevelSetPseudoNormals");
		inputParams.setLabel("Mesh to Levelset using Pseudonormals");

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces.devel");
		
		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine L Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Computes a level set representation of an input surface mesh using the method of pseudonormals.");
		
		info.setVersion("3.0.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
		

	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(result=new ParamVolume("Level set"));
		//outputParams.add(result2=new ParamVolume("Signed distance function"));
		//outputParams.add(resultSurf=new ParamSurface("Aligned Surface (opt)"));
		//resultSurf.setMandatory(false);
	}
	
	
	protected void execute(CalculationMonitor monitor) {
		
		ImageDataFloat	origVolData = new ImageDataFloat(origVol.getImageData());
		int nx = origVolData.getRows();
		int ny = origVolData.getCols();
		int nz = origVolData.getSlices();
		float rx = origVolData.getHeader().getDimResolutions()[0];
		float ry = origVolData.getHeader().getDimResolutions()[1];
		float rz = origVolData.getHeader().getDimResolutions()[2];
		float ox = origVolData.getHeader().getOrigin()[0];
		float oy = origVolData.getHeader().getOrigin()[1];
		float oz = origVolData.getHeader().getOrigin()[2];
		
		
		EmbeddedSurface surf=origSurf.getSurface();
		
		// maps the surface to the voxel space
		if (mipavTransform.getValue().booleanValue()) {
		surf.scaleVertices(new float[]{-1.0f, 1.0f, -1.0f});
		Point3f pt0 = new Point3f((nx-1)*rx/2.0f, (ny-1)*ry/2.0f, (nz-1)*rz/2.0f);
		surf.translate(pt0);
		}
		
		// scale the surface data into voxel space
		surf.scaleVertices(new float[]{1.0f/rx, -1.0f/ry, -1.0f/rz});
		Point3f pt1 = new Point3f(0, (ny-1), (nz-1));
		surf.translate(pt1);
		
		boolean[][][] mask = new boolean[nx][ny][nz];
		float[][][] df = new float[nx][ny][nz];
		for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) {
			df[x][y][z] = PADDING;
			mask[x][y][z] = false;
		}
		
		float[][][] sum = new float[nx][ny][nz];
		
		//BasicInfo.displayMessage("Repairing degenerate triangles.\n");
		//surf.repairDegenerateTriangles(0.01f,0.005f);
		
		BasicInfo.displayMessage("Number of faces/triangles = "+surf.getFaceCount()+"\n");
		surf.buildAllTables();
		
		BasicInfo.displayMessage("Calculate face normals. Looping through triangles...\n");
		EmbeddedSurface.Face[] faces = surf.getFaces();
		Vector3f[] faceNormals = new Vector3f[surf.getFaceCount()];
		Vector3f edge0 = new Vector3f();  
		Vector3f edge1 = new Vector3f();  
		Vector3f edge2 = new Vector3f();  
		for(int fid=0; fid<surf.getFaceCount(); fid++) {
			Point3f[] p = surf.getFacePoints(fid);
			edge2.sub(p[1], p[0]);
			edge1.sub(p[2], p[0]);
			faceNormals[fid] = new Vector3f();
			faceNormals[fid].cross(edge2, edge1);
			faceNormals[fid].normalize();
		}
		
		BasicInfo.displayMessage("Calculate edge normals. Looping through edges...\n");
		Vector3f[] edgeNormals = new Vector3f[surf.getIndexCount()/2];
		EmbeddedSurface.Face[][] neighborEdgeFaceTable = surf.getNeighborEdgeFaceTable();
		for(int eid=0; eid<surf.getIndexCount()/2; eid++) {
			int fid1 = neighborEdgeFaceTable[eid][0].id;
			int fid2 = neighborEdgeFaceTable[eid][1].id;
			edgeNormals[eid] = new Vector3f();
			edgeNormals[eid].add(faceNormals[fid1], faceNormals[fid2]);
			edgeNormals[eid].scale((float) FastMath.PI);
		}
		
		BasicInfo.displayMessage("Calculate vertex normals. Looping through vertices...\n");
		Vector3f[] vertexNormals = new Vector3f[surf.getVertexCount()];
		EmbeddedSurface.Face[][] neighborVertexFaceTable = surf.getNeighborVertexFaceTable();
		for(int i=0; i<surf.getVertexCount(); i++) {
			vertexNormals[i] = new Vector3f();
			for(int j=0; j<neighborVertexFaceTable[i].length; j++) {
				int fid = neighborVertexFaceTable[i][j].id;
				int[] id = surf.getFaceVertexIds(fid);
				if (i==id[0]) {
					edge0.sub(surf.getVertex(id[1]), surf.getVertex(i));
					edge1.sub(surf.getVertex(id[2]), surf.getVertex(i));
				} else if (i==id[1]) {
					edge0.sub(surf.getVertex(id[0]), surf.getVertex(i));
					edge1.sub(surf.getVertex(id[2]), surf.getVertex(i));
				} else if (i==id[2]) {
					edge0.sub(surf.getVertex(id[1]), surf.getVertex(i));
					edge1.sub(surf.getVertex(id[0]), surf.getVertex(i));
				}
				vertexNormals[i].scaleAdd(edge0.angle(edge1), faceNormals[fid], vertexNormals[i]);
				//vertexNormals[i].scale((float)(0.5/FastMath.PI));
			}
		}
		
		for(int fid=0; fid<surf.getFaceCount(); fid++) {
			faceNormals[fid].scale((float)(2*FastMath.PI));
		}
		
		
		BasicInfo.displayMessage("Calculating distance field. Looping through triangles...\n");
		//For each triangle/face of the mesh, calculate distance of each voxel within a bounding box to plane, vertices and edges
		//fid: face Id		
		
		for(int fid=0; fid<surf.getFaceCount(); fid++) {
			
			//bounding box around triangle
			float a = surf.getMinAngle(fid);
			int[] id = surf.getFaceVertexIds(fid);
			int eid;
			Point3f[] p = surf.getFacePoints(fid);
			Point3f face = surf.getCentroid(fid);
			
			Vector3f normVect = new Vector3f();
			Vector3f n = new Vector3f();
			Vector3f distVect = new Vector3f();
			Point3f[] e = new Point3f[3];
			
			edge2.sub(p[1], p[0]);
			edge1.sub(p[2], p[0]);
			edge0.sub(p[2], p[1]);
			normVect.cross(edge2, edge1);
			normVect.normalize();
		
			//e[2] = new Point3f((p[0].x + p[1].x) * 0.5f, (p[0].y + p[1].y) * 0.5f,(p[0].z + p[1].z) * 0.5f);
			//e[1] = new Point3f((p[0].x + p[2].x) * 0.5f, (p[0].y + p[2].y) * 0.5f,(p[0].z + p[2].z) * 0.5f);
			//e[0] = new Point3f((p[2].x + p[1].x) * 0.5f, (p[2].y + p[1].y) * 0.5f,(p[2].z + p[1].z) * 0.5f);
		
			int xmin = Numerics.floor(Numerics.min(p[0].x, p[1].x, p[2].x))-narrowbandDist;
			int xmax = Numerics.ceil(Numerics.max(p[0].x, p[1].x, p[2].x))+narrowbandDist;
			int ymin = Numerics.floor(Numerics.min(p[0].y, p[1].y, p[2].y))-narrowbandDist;
			int ymax = Numerics.ceil(Numerics.max(p[0].y, p[1].y, p[2].y))+narrowbandDist;
			int zmin = Numerics.floor(Numerics.min(p[0].z, p[1].z, p[2].z))-narrowbandDist;
			int zmax = Numerics.ceil(Numerics.max(p[0].z, p[1].z, p[2].z))+narrowbandDist;
			
			//for every voxel in box, calculate min distance to mesh
			//could improve accuracy by calculating the min distance to the face plane instead of distance to face centroid
			for(int x=xmin; x<=xmax; x++) for(int y=ymin; y<=ymax; y++) for(int z=zmin; z<=zmax; z++) {
				
				//df will be defined at this voxel
				mask[x][y][z] = true;
				
				Point3f p0 = new Point3f(x,y,z);
				float d;
				
				
				//distance to 3 vertices
				//vertex 0
				distVect.sub(p0,p[0]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					n.set(vertexNormals[id[0]]);
					sum[x][y][z] = n.dot(distVect);
					
				}
				
				//vertex 1
				distVect.sub(p0,p[1]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					n.set(vertexNormals[id[1]]);
					sum[x][y][z] = n.dot(distVect);
					
				}
				
				//vertex 2
				distVect.sub(p0,p[2]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					n.set(vertexNormals[id[2]]);
					sum[x][y][z] = n.dot(distVect);
					
				} 
				
				
				//distance to face (use face centroid)
				distVect.sub(p0,face);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					n.set(faceNormals[fid]);
					n.scale((float) (2.0f*FastMath.PI));
					sum[x][y][z] = n.dot(distVect);
					
				} 
				
				
				//distance to 3 edges (use edge midpoint)
				//edge 0
				e[0] = surf.midPoint(faces[fid].edges[0]);
				distVect.sub(p0,e[0]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					eid = faces[fid].edges[0].id;
					n.set(edgeNormals[eid]);
					sum[x][y][z] = n.dot(distVect);
					
				} 
				
				//edge 1
				e[1] = surf.midPoint(faces[fid].edges[1]);
				distVect.sub(p0,e[1]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					eid = faces[fid].edges[1].id;
					n.set(edgeNormals[eid]);
					sum[x][y][z] = n.dot(distVect);
					
				} 
				
				//edge 2
				e[2] = surf.midPoint(faces[fid].edges[2]);
				distVect.sub(p0,e[2]);
				d = distVect.length();
				if (d<df[x][y][z]) {
					df[x][y][z] = d;
					eid = faces[fid].edges[2].id;
					n.set(edgeNormals[eid]);
					sum[x][y][z] = n.dot(distVect);
					
				} 
				
			}
		
		}
		
		BasicInfo.displayMessage("done.\n");
		BasicInfo.displayMessage("Adding sign to distance field.\n");
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
		
		BasicInfo.displayMessage("Padding levelset.\n");
		sdf = padLevelset(sdf, datamask, 6, nx, ny, nz);
		//sdf = padLevelset(sdf, datamask, 26, nx, ny, nz);
		//sdf = padLevelset(sdf, datamask, 18, nx, ny, nz);
		
		for (int x=3; x<nx-3; x++) for (int y=3; y<ny-3; y++) for (int z=3; z<nz-3; z++) {
			int xyz = x+nx*y+nx*ny*z;
			datamask[xyz] = true;
			//df[x][y][z] = sdf[xyz];
		}
		
		BasicInfo.displayMessage("Removing sign errors.\n");
		sdf = cleanupPaddingSign(sdf, datamask, 26, nx, ny, nz);
		
		BasicInfo.displayMessage("Clean up sign 2.\n");
		sdf = cleanupSign2(sdf, 26, nx, ny, nz);
		
		BasicInfo.displayMessage("Evolving narrow band.\n");
		InflateGdm gdm = new InflateGdm(sdf, nx, ny, nz, rx, ry, rz, bgmask, 0.4f, 0.4f, "no",null);
		gdm.evolveNarrowBand(0, 1.0f);
		
		BasicInfo.displayMessage("Output.\n");
		
		ImageDataFloat resultData = new ImageDataFloat(gdm.exportLevelset());		
		//ImageDataFloat resultData = new ImageDataFloat(df);
		resultData.setHeader(origVol.getImageData().getHeader());
		resultData.setName(origSurf.getName()+"_ls");
		result.setValue(resultData);
		resultData = null;
		
		/*
		ImageDataFloat resultData2 = new ImageDataFloat(df);
		resultData2.setHeader(origVol.getImageData().getHeader());
		resultData2.setName(origSurf.getName()+"_sdf");
		result2.setValue(resultData2);
		resultData2 = null;
		*/
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
			BasicInfo.displayMessage("Iteration "+itrs);
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
			BasicInfo.displayMessage(" Count = "+count+"\n");
			for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
				int xyz = x+nx*y+nx*ny*z;
				sdf[xyz] = cleansdf[xyz];
			}
		}
		
		return cleansdf;
			
	}
		
	private final float[] padLevelset(float[] data, boolean[] datamask, int connectivity, int nx, int ny, int nz) {
		BasicInfo.displayMessage("Pad levelset.\n");
		
		float[] dilData = new float[nx*ny*nz];
		
		int[] xoff = new int[]{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};
		int[] yoff = new int[]{0, 0, nx, -nx, 0, 0, nx, -nx, -nx, nx, 0, 0, 0, 0, nx, -nx, nx, -nx, nx, nx, -nx, -nx, nx, nx, -nx, -nx};
		int[] zoff = new int[]{0, 0, 0, 0, nx*ny, -nx*ny, 0, 0, 0, 0, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny,-nx*ny, nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny, nx*ny, -nx*ny};
		
		BitSet boundary = new BitSet(nx*ny*nz);
        BitSet nextboundary = new BitSet(nx*ny*nz);
        BitSet used = new BitSet(nx*ny*nz);
        BitSet nextused = new BitSet(nx*ny*nz);
        
        BasicInfo.displayMessage("\n Set initial processed data");
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
	
}
