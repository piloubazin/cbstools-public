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
 * @author Pierre-Louis Bazin and Christine L Tardif
 */
public class JistRegistrationSurfaceMeshGroupAverage extends ProcessingAlgorithm {

	// jist containers
	private ParamVolumeCollection levelsetImages;
	private ParamVolumeCollection deformImages;
	private ParamSurface targetSurface;
	private ParamVolumeCollection contrastImages;
	private ParamOption  	mappingOption;
	private static final String[]		mappingTypes = {"raw_average","projected"};
	private ParamOption  	interpOption;
	private static final String[]		interpTypes = {"linear","nearest"};
	private	ParamBoolean includeTarget;
	private	ParamBoolean computeStdev;
	
	private ParamSurface outTargetSurface;
	private ParamSurface averageSurface;
	private ParamSurface stdevSurface;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private int nsx,nsy,nsz,nsxyz;
	private int ntx,nty,ntz,ntxyz;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(levelsetImages = new ParamVolumeCollection("Source Levelset Surfaces"));
		inputParams.add(deformImages = new ParamVolumeCollection("Deformation Fields"));
		inputParams.add(targetSurface=new ParamSurface("Target Surface Mesh"));
		inputParams.add(mappingOption=new ParamOption("Averaging method", mappingTypes));
		inputParams.add(contrastImages = new ParamVolumeCollection("Source Surface Contrast (opt)"));
		inputParams.add(includeTarget=new ParamBoolean("include target", false));
		inputParams.add(computeStdev=new ParamBoolean("compute average distance", true));
		inputParams.add(interpOption=new ParamOption("Interpolation method", interpTypes));
		levelsetImages.setLoadAndSaveOnValidate(false);
		deformImages.setLoadAndSaveOnValidate(false);
		contrastImages.setMandatory(false);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration.devel");
		inputParams.setLabel("Surface Mesh Group Average");
		inputParams.setName("SurfaceMeshGroupAverage");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Christine L. Tardif", "ctardif@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Creates a geometric average of original cortical surfaces after surface alignment.");
		
		info.setVersion("3.0.4");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(averageSurface=new ParamSurface("Average Surface"));
		outputParams.add(stdevSurface=new ParamSurface("Mapped average distance (opt)"));
		outputParams.add(outTargetSurface=new ParamSurface("Target Surface with Average Contrast"));
		stdevSurface.setMandatory(false);
		
		outputParams.setName("average surface");
		outputParams.setLabel("average surface");
	}

	@Override
	protected void execute(CalculationMonitor monitor){

		nsx = levelsetImages.getImageDataList().get(0).getRows();
		nsy = levelsetImages.getImageDataList().get(0).getCols();
		nsz = levelsetImages.getImageDataList().get(0).getSlices();
		nsxyz = nsx*nsy*nsz;
		float rsx = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[0];
		float rsy = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[1];
		float rsz = levelsetImages.getImageDataList().get(0).getHeader().getDimResolutions()[2];
		System.out.println("Source image dimensions: "+nsx+" x "+nsy+" x "+nsz);
		System.out.println("Source Image resolutions: "+rsx+" x "+rsy+" x "+rsz);
		
		ntx = deformImages.getImageDataList().get(0).getRows();
		nty = deformImages.getImageDataList().get(0).getCols();
		ntz = deformImages.getImageDataList().get(0).getSlices();
		ntxyz = ntx*nty*ntz;
		float rtx = deformImages.getImageDataList().get(0).getHeader().getDimResolutions()[0];
		float rty = deformImages.getImageDataList().get(0).getHeader().getDimResolutions()[1];
		float rtz = deformImages.getImageDataList().get(0).getHeader().getDimResolutions()[2];
		System.out.println("Target image dimensions: "+ntx+" x "+nty+" x "+ntz);
		System.out.println("Target Image resolutions: "+rtx+" x "+rty+" x "+rtz);
		
		
		EmbeddedSurface target = targetSurface.getSurface();
		/*
		target.scaleVertices(new float[]{1.0f/rx, -1.0f/ry, -1.0f/rz});
		Point3f pt1 = new Point3f(0, (ny-1), (nz-1));
		target.translate(pt1);
		*/
		
		int Nsubjects = levelsetImages.getImageDataList().size();
		BasicInfo.displayMessage("Number of subjectS: "+Nsubjects);
			
		
		int Nvertex = target.getVertexCount();
		
		Point3f[] avg = new Point3f[Nvertex];
		float[] count = new float[Nvertex];
		double[][] data = new double[Nvertex][1];
		
		if (includeTarget.getValue().booleanValue()) {
			for(int i=0; i<Nvertex; i++) {
				avg[i] = new Point3f(target.getVertex(i));
				// map to voxel space
				avg[i].x = avg[i].x/rtx;
				avg[i].y = (nty-1)-avg[i].y/rty;
				avg[i].z = (ntz-1)-avg[i].z/rtz;
				
				count[i] = 1.0f;
				if (contrastImages.getValue()!=null) {
					data[i][0] = target.getVertexData()[i][0];
				}
			}
		} else {
			for(int i=0; i<Nvertex; i++) {
				avg[i] = new Point3f(new float[]{0.0f,0.0f,0.0f});
				count[i] = 0.0f;
				data[i][0] = 0.0;
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
			
			if (contrastImages.getImageDataList()!=null && contrastImages.size()==Nsubjects) {
				BasicInfo.displayMessage("contrasts "+contrastImages.size());
				if (contrastImages.getImageDataList().get(s)!=null) {
					contrast = (new ImageDataFloat(contrastImages.getImageDataList().get(s))).toArray3d();
					contrastImages.dispose();
				}
			}
			
			for(int i=0; i<Nvertex; i++){
				Point3f p = target.getVertex(i);
				// makes the transition from MIPAV meshes to voxel space here
				int x = Numerics.floor(p.x/rtx);
				int y = Numerics.floor((nty-1)-p.y/rty);
				int z = Numerics.floor((ntz-1)-p.z/rtz);
				if (x>0 && x<ntx-2 && y>0 && y<nty-2 && z>0 && z<ntz-2) {
					float dx = p.x/rtx - x;
					float dy = (nty-1)-p.y/rty - y;
					float dz = (ntz-1)-p.z/rtz - z;
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
						count[i]++;
						if (contrast!=null) {
							if (interpOption.getValue().equals("linear")) { 
								data[i][0] += ImageInterpolation.linearClosestInterpolation(contrast, Xs, Ys, Zs, nsx, nsy, nsz);
							} else {
								data[i][0] += ImageInterpolation.nearestNeighborInterpolation(contrast, Xs, Ys, Zs, nsx, nsy, nsz);
							}
						}
					} else if (mappingOption.getValue().equals("projected")) {
						float[] proj = new float[3];
						projectToLevelset(levelset, new float[]{Xs,Ys,Zs}, proj);
						avg[i].x += proj[X];	
						avg[i].y += proj[Y];	
						avg[i].z += proj[Z];	
						count[i]++;
						if (contrast!=null) {
							if (interpOption.getValue().equals("linear")) { 
								data[i][0] += ImageInterpolation.linearClosestInterpolation(contrast, proj[X], proj[Y], proj[Z], nsx, nsy, nsz);
							} else {
								data[i][0] += ImageInterpolation.nearestNeighborInterpolation(contrast, proj[X], proj[Y], proj[Z], nsx, nsy, nsz);
							}	
						}
					}
				}		
			}
		}
		
		for(int i=0; i<Nvertex; i++) if (count[i]!=0) {
			avg[i].x /= count[i];
			avg[i].y /= count[i];
			avg[i].z /= count[i];
			if (contrastImages.getValue()!=null) data[i][0] /= count[i];
		}
		
		float[] dist;
		// second pass: dist (opt)
		if (computeStdev.getValue().booleanValue()) {
			BasicInfo.displayMessage("Second pass: estimate average distance");
			dist = new float[Nvertex];
			for (int s=0;s<Nsubjects;s++) {
				BasicInfo.displayMessage("Processing subject "+(s+1));
				float[][][] levelset = (new ImageDataFloat(levelsetImages.getImageDataList().get(s))).toArray3d();
				levelsetImages.dispose();
				float[][][][] deform = (new ImageDataFloat(deformImages.getImageDataList().get(s))).toArray4d();
				deformImages.dispose();
				
				for(int i=0; i<Nvertex; i++){
					Point3f p = target.getVertex(i);
					// makes the transition from MIPAV meshes to voxel space here
					int x = Numerics.floor(p.x/rtx);
					int y = Numerics.floor((nty-1)-p.y/rty);
					int z = Numerics.floor((ntz-1)-p.z/rtz);
					if (x>0 && x<ntx-2 && y>0 && y<nty-2 && z>0 && z<ntz-2) {
						float dx = p.x/rtx - x;
						float dy = (nty-1)-p.y/rty - y;
						float dz = (ntz-1)-p.z/rtz - z;
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
							dist[i] += (Xs-avg[i].x)*(Xs-avg[i].x)*rsx*rsx;	
							dist[i] += (Ys-avg[i].y)*(Ys-avg[i].y)*rsy*rsy;	
							dist[i] += (Zs-avg[i].z)*(Zs-avg[i].z)*rsz*rsz;	
						} else if (mappingOption.getValue().equals("projected")) {
							float[] proj = new float[3];
							projectToLevelset(levelset, new float[]{Xs,Ys,Zs}, proj);
							dist[i] += (proj[X]-avg[i].x)*(proj[X]-avg[i].x)*rsx*rsx;	
							dist[i] += (proj[Y]-avg[i].y)*(proj[Y]-avg[i].y)*rsy*rsy;	
							dist[i] += (proj[Z]-avg[i].z)*(proj[Z]-avg[i].z)*rsz*rsz;	
						}
					}		
				}
			}
			
			for(int i=0; i<Nvertex; i++) if (count[i]>1) {
				dist[i] = (float)FastMath.sqrt(dist[i]/count[i]);
			} else {
				dist[i] = 0.0f;
			}
			
			EmbeddedSurface average2 = target.clone();
			for(int i=0; i<Nvertex; i++) {
				average2.setVertex(i, avg[i]);
				average2.setVertexData(i, dist[i]);
			}
			
			average2.scaleVertices(new float[]{rsx, -rsy, -rsz});
			Point3f pt = new Point3f(0, (nsy-1)*rsy, (nsz-1)*rsz);
			average2.translate(pt);
			average2.setName(target.getName()+"_groupavgdist");
			stdevSurface.setValue(average2);
			
		}
		
		// copy target surface mesh
		EmbeddedSurface average = target.clone();
		for(int i=0; i<Nvertex; i++) {
			average.setVertex(i, avg[i]);
		}
		average.setVertexData(data);
		target.setVertexData(data);
		
		average.scaleVertices(new float[]{rsx, -rsy, -rsz});
		Point3f pt = new Point3f(0, (nsy-1)*rsy, (nsz-1)*rsz);
		average.translate(pt);
		average.setName(target.getName()+"_groupavgsurf");
		averageSurface.setValue(average);
		
		target.setName(target.getName()+"_trgsurf_groupavgdata");
		outTargetSurface.setValue(target);
	
	}

	private final float projectToLevelset(float[][][] levelset, float[] pt0, float[] pt) {
		
		// best gradient approximation among various choices (most regular)
		double I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]-1.0f, pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X]+1.0f, pt0[Y], pt0[Z], nsx, nsy, nsz);
		double Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]-1.0f, pt0[Z], nsx, nsy, nsz);
		double Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y]+1.0f, pt0[Z], nsx, nsy, nsz);
		double Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]-1.0f, nsx, nsy, nsz);
		double Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z]+1.0f, nsx, nsy, nsz);
		
		//if (Double.isNaN(I) || Double.isNaN(Imx) || Double.isNaN(Ipx) || Double.isNaN(Imy) || Double.isNaN(Ipy) || Double.isNaN(Imz) || Double.isNaN(Ipz))
		//	System.out.println("NaN: levelset ("+pt0[X]+", "+pt0[Y]+", "+pt0[Z]+")");
		
		double length = I;

		double Dx = 0.5*(Ipx-Imx);
		double Dy = 0.5*(Ipy-Imy);
		double Dz = 0.5*(Ipz-Imz);
		
		double grad = FastMath.sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		
		//if (Double.isNaN(grad)) System.out.println("NaN: gradient ("+Dx+", "+Dy+", "+Dz+")");

		double res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt0[X], pt0[Y], pt0[Z], nsx, nsy, nsz);
		
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
			res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
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
				I = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
				Imx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]-1.0f, pt[Y], pt[Z], nsx, nsy, nsz);
				Ipx = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X]+1.0f, pt[Y], pt[Z], nsx, nsy, nsz);
				Imy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]-1.0f, pt[Z], nsx, nsy, nsz);
				Ipy = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y]+1.0f, pt[Z], nsx, nsy, nsz);
				Imz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]-1.0f, nsx, nsy, nsz);
				Ipz = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z]+1.0f, nsx, nsy, nsz);
				
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
				res = ImageInterpolation.linearInterpolation(levelset, 0.0f, pt[X], pt[Y], pt[Z], nsx, nsy, nsz);
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
