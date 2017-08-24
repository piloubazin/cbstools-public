package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolumeCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.libraries.*;
import de.mpg.cbs.methods.*;

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3f;
import org.apache.commons.math3.util.FastMath;

/*
 * @author Pierre-Louis Bazin, Christine Lucas Tardif
 */
public class JistRegistrationSurfaceMeshGroupData extends ProcessingAlgorithm {

	// jist containers
	private ParamVolumeCollection 	levelsetImages;
	private ParamVolumeCollection 	deformImages;
	private ParamSurface 			targetSurface;
	private ParamVolumeCollection 	contrastImages;
	private ParamOption  			mappingOption;
	private static final String[]	mappingTypes = {"projected","raw_average"};
	private	ParamBoolean 			includeTarget;
	//private	ParamBoolean 			computeStdev;
	
	private ParamSurface 			averageSurface;
	private ParamSurface 			outTargetSurface;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private int nx,ny,nz,nxyz;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImages = new ParamVolumeCollection("Source Levelset Surfaces"));
		inputParams.add(deformImages = new ParamVolumeCollection("Deformation Fields"));
		inputParams.add(targetSurface=new ParamSurface("Target Surface Mesh"));
		inputParams.add(mappingOption=new ParamOption("Averaging method", mappingTypes));
		inputParams.add(contrastImages = new ParamVolumeCollection("Source Surface Contrast (opt)"));
		inputParams.add(includeTarget=new ParamBoolean("include target", false));
		levelsetImages.setLoadAndSaveOnValidate(false);
		deformImages.setLoadAndSaveOnValidate(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration.devel");
		inputParams.setLabel("Surface Mesh Group Data");
		inputParams.setName("SurfaceMeshGroupData");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Christine L. Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Creates a geometric average of original cortical surfaces after surface alignment and embbeds the data from each subject.");
		
		info.setVersion("3.0.4");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(outTargetSurface=new ParamSurface("Target Surface with Data"));
		outputParams.add(averageSurface=new ParamSurface("Average Surface with Data"));
		
		outputParams.setName("group surface data");
		outputParams.setLabel("group surface data");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		nx = levelsetImages.getImageDataList().get(0).getRows();
		ny = levelsetImages.getImageDataList().get(0).getCols();
		nz = levelsetImages.getImageDataList().get(0).getSlices();
		nxyz = nx*ny*nz;
		float rx = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[0];
		float ry = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[1];
		float rz = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[2];
		System.out.println("Image dimensions: "+nx+" x "+ny+" x "+nz);
		System.out.println("Image resolutions: "+rx+" x "+ry+" x "+rz);
		
		EmbeddedSurface target = targetSurface.getSurface();
		/*
		target.scaleVertices(new float[]{1.0f/rx, -1.0f/ry, -1.0f/rz});
		Point3f pt1 = new Point3f(0, (ny-1), (nz-1));
		target.translate(pt1);
		*/
		
		int Nsubjects = levelsetImages.getImageDataList().size();

		int Nvertex = target.getVertexCount();
		
		Point3f[] avg = new Point3f[Nvertex];
		int count;
		double[][] data;
		
		if (includeTarget.getValue().booleanValue()) {
			data = new double[Nvertex][Nsubjects+1];
			for(int i=0; i<Nvertex; i++) {
				avg[i] = new Point3f(target.getVertex(i));
				data[i][0] = target.getVertexData()[i][0];
			}
		} else {
			data = new double[Nvertex][Nsubjects];
			for(int i=0; i<Nvertex; i++) {
				avg[i] = new Point3f(new float[]{0.0f,0.0f,0.0f});
			}
		}
		
		BasicInfo.displayMessage("First pass: compute average surface");
		for (int s=0;s<Nsubjects;s++) {
			BasicInfo.displayMessage("Processing subject "+(s+1));
			float[][][] levelset = (new ImageDataFloat(levelsetImages.getImageDataList().get(s))).toArray3d();
			levelsetImages.dispose();
			float[][][][] deform = (new ImageDataFloat(deformImages.getImageDataList().get(s))).toArray4d();
			deformImages.dispose();
			float[][][] contrast = null;
			if (contrastImages.getImageDataList()!=null) {
				contrast = (new ImageDataFloat(contrastImages.getImageDataList().get(s))).toArray3d();
				contrastImages.dispose();
			}
			
			if (includeTarget.getValue().booleanValue()) {
				count = s+1;
			} else {
				count = s;
			}
			for(int i=0; i<Nvertex; i++){
				Point3f p = target.getVertex(i);
				// makes the transition from MIPAV meshes to voxel space here
				int x = Numerics.floor(p.x/rx);
				int y = Numerics.floor((ny-1)-p.y/ry);
				int z = Numerics.floor((nz-1)-p.z/rz);
				if (x>0 && x<nx-2 && y>0 && y<ny-2 && z>0 && z<nz-2) {
					float dx = p.x/rx - x;
					float dy = (ny-1)-p.y/ry - y;
					float dz = (nz-1)-p.z/rz - z;
					float Xs = (1-dx)*(1-dy)*(1-dz)*deform[x][y][z][X]+dx*dy*dz*deform[x+1][y+1][z+1][X]
								+dx*(1-dy)*(1-dz)*deform[x+1][y][z][X]+(1-dx)*dy*(1-dz)*deform[x][y+1][z][X]+(1-dx)*(1-dy)*dz*deform[x][y][z+1][X]
								+dx*dy*(1-dz)*deform[x+1][y+1][z][X]+(1-dx)*dy*dz*deform[x][y+1][z+1][X]+dx*(1-dy)*dz*deform[x+1][y][z+1][X];
								
					float Ys = (1-dx)*(1-dy)*(1-dz)*deform[x][y][z][Y]+dx*dy*dz*deform[x+1][y+1][z+1][Y]
								+dx*(1-dy)*(1-dz)*deform[x+1][y][z][Y]+(1-dx)*dy*(1-dz)*deform[x][y+1][z][Y]+(1-dx)*(1-dy)*dz*deform[x][y][z+1][Y]
								+dx*dy*(1-dz)*deform[x+1][y+1][z][Y]+(1-dx)*dy*dz*deform[x][y+1][z+1][Y]+dx*(1-dy)*dz*deform[x+1][y][z+1][Y];
								
					float Zs = (1-dx)*(1-dy)*(1-dz)*deform[x][y][z][Z]+dx*dy*dz*deform[x+1][y+1][z+1][Z]
								+dx*(1-dy)*(1-dz)*deform[x+1][y][z][Z]+(1-dx)*dy*(1-dz)*deform[x][y+1][z][Z]+(1-dx)*(1-dy)*dz*deform[x][y][z+1][Z]
								+dx*dy*(1-dz)*deform[x+1][y+1][z][Z]+(1-dx)*dy*dz*deform[x][y+1][z+1][Z]+dx*(1-dy)*dz*deform[x+1][y][z+1][Z];
					
					if (mappingOption.getValue().equals("raw_average")) {
						avg[i].x += Xs;	
						avg[i].y += Ys;	
						avg[i].z += Zs;	
						if (contrast!=null) {
							data[i][count] += ImageInterpolation.linearClosestInterpolation(contrast, Xs, Ys, Zs, nx, ny, nz);
						}
						
					} else if (mappingOption.getValue().equals("projected")) {
						float[] proj = new float[3];
						projectToLevelset(levelset, new float[]{Xs,Ys,Zs}, proj);
						avg[i].x += proj[X];	
						avg[i].y += proj[Y];	
						avg[i].z += proj[Z];	
						if (contrast!=null) {
							data[i][count] += ImageInterpolation.linearClosestInterpolation(contrast, proj[X], proj[Y], proj[Z], nx, ny, nz);
						}
						
					}
				}		
			}
		}
		
		for(int i=0; i<Nvertex; i++) { //if (count[i]!=0) {
			avg[i].x /= Nsubjects;
			avg[i].y /= Nsubjects;
			avg[i].z /= Nsubjects;
		}
		
		// copy target surface mesh
		EmbeddedSurface average = target.clone();
		for(int i=0; i<Nvertex; i++) {
			average.setVertex(i, avg[i]);
			average.setVertexData(i, data[i]);
			target.setVertexData(i, data[i]);
		}
		
		average.scaleVertices(new float[]{rx, -ry, -rz});
		Point3f pt = new Point3f(0, (ny-1)*ry, (nz-1)*rz);
		average.translate(pt);
		average.setName(target.getName()+"_avgsurf_groupdata");
		averageSurface.setValue(average);
		
		target.setName(target.getName()+"_trgsurf_groupdata");
		outTargetSurface.setValue(target);
		
	}

	private final float projectToLevelset(float[][][] levelset, float[] pt0, float[] pt) {
		
		// best gradient approximation among various choices (most regular)
		double I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		double Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]-1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]+1.0f, pt0[Y], pt0[Z], nx, ny, nz);
		double Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]-1.0f, pt0[Z], nx, ny, nz);
		double Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]+1.0f, pt0[Z], nx, ny, nz);
		double Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]-1.0f, nx, ny, nz);
		double Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]+1.0f, nx, ny, nz);
		
		//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
		//	System.out.println("NaN: levelset ("+pt0[X]+", "+pt0[Y]+", "+pt0[Z]+")");
		
		double length = I;

		double Dx = 0.5*(Ipx-Imx);
		double Dy = 0.5*(Ipy-Imy);
		double Dz = 0.5*(Ipz-Imz);
		
		double grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		
		//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

		double res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nx, ny, nz);
		
		int t=0;

		if (grad > 0.01) {
			// closed form approximation to the closest point on the N-th layer surface 
			// (accurate if close enough to approximate the surface by a plane)
			pt[X] = (float)(pt0[X] - length*Dx/grad);
			pt[Y] = (float)(pt0[Y] - length*Dy/grad);
			pt[Z] = (float)(pt0[Z] - length*Dz/grad);
			
			//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
			//	System.out.println("NaN: coord ("+pt0[X]+","+pt0[Y]+", "+pt0[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
			
			// correction: go back halfway and recompute gradient, etc: should converge toward the surface
			// approaches using small steps are not very good (unstable)
			// approaches doing the full correction from the new solution are not very stable either
			res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
			while (res*res>0.0001 && t<100 && grad>0.01) {
				t++;
				//System.out.print(".");
				
				// correction: come back halfway, recompute (going further back does not improve the result)
				pt[X] = (float)(pt[X] + 0.5*length*Dx/grad);
				pt[Y] = (float)(pt[Y] + 0.5*length*Dy/grad);
				pt[Z] = (float)(pt[Z] + 0.5*length*Dz/grad);
				
				//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
				//	System.out.println("NaN: coord 0.5 ("+pt[X]+","+pt[Y]+", "+pt[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
			
				// correction: project once more onto the curve (with / without minimization of distance to origin point)
				I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
				Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nx, ny, nz);
				Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nx, ny, nz);
				Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nx, ny, nz);
				Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nx, ny, nz);
				Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nx, ny, nz);
				Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nx, ny, nz);
				
				//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
				//	System.out.println("NaN: levelset ("+pt[X]+", "+pt[Y]+", "+pt[Z]+")");

				length = I;
				Dx = 0.5*(Ipx-Imx);
				Dy = 0.5*(Ipy-Imy);
				Dz = 0.5*(Ipz-Imz);
				
				grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
				//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

				if (grad > 0.01) {
					pt[X] = (float)(pt[X] - length*Dx/grad);
					pt[Y] = (float)(pt[Y] - length*Dy/grad);
					pt[Z] = (float)(pt[Z] - length*Dz/grad);
				
					//if (Float.isNaN(pt[X]) || Float.isNaN(pt[Y]) || Float.isNaN(pt[Z]))
					//	System.out.println("NaN: coord 2 ("+pt[X]+","+pt[Y]+", "+pt[Z]+"| "+length+"x "+Dx+", "+Dy+", "+Dz+" / "+grad);
				}
				res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nx, ny, nz);
			}

		} else {
			//System.err.print(".");
			pt[X] = pt0[X];
			pt[Y] = pt0[Y];
			pt[Z] = pt0[Z];
		}
		return (float)res;
	}

}
