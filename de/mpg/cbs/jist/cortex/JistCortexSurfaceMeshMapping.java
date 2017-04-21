package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamString;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
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

/*
 * @author Pierre-Louis Bazin
 */
public class JistCortexSurfaceMeshMapping extends ProcessingAlgorithm {

	// jist containers
	private ParamVolume intensityImage;
	private ParamString intensityName;
	private ParamSurface inflatedSurface;
	private ParamSurface origSurface;
	private ParamOption  	mappingOption;
	private static final String[]		mappingTypes = {"closest_point","linear_interp","highest_value"};
	
	private ParamSurface mappedInfSurface;
	private ParamSurface mappedOrgSurface;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;

	private static final byte NONE = 10;
	private static final byte CLOSEST = 11;
	private static final byte LINEAR = 12;
	private static final byte HIGHEST = 13;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(intensityImage = new ParamVolume("Intensity Image"));
		inputParams.add(intensityName = new ParamString("Image Label Name"));
		inputParams.add(origSurface=new ParamSurface("Original Surface"));
		inputParams.add(inflatedSurface=new ParamSurface("Inflated Surface (opt)"));
		inflatedSurface.setMandatory(false);
		inputParams.add(mappingOption=new ParamOption("Mapping method", mappingTypes));

		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing");
		inputParams.setLabel("Surface Mesh Mapping");
		inputParams.setName("SurfaceMeshMapping");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Juliane Dinse", "dinse@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.add(new AlgorithmAuthor("Marcel Weiss", "weiss@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Map some intensity image on a surface or pair of original+inflated surfaces (use 'Image Label' as suffix for the mapped data).");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(mappedOrgSurface=new ParamSurface("Mapped Original Surface"));
		outputParams.add(mappedInfSurface=new ParamSurface("Mapped Inflated Surface (opt)"));
		mappedInfSurface.setMandatory(false);

		outputParams.setName("mapped surfaces");
		outputParams.setLabel("mapped surfaces");
	}

	@Override
	protected void execute(CalculationMonitor monitor){
		
		ImageDataFloat	intensImg = new ImageDataFloat(intensityImage.getImageData());
		
		int nx = intensImg.getRows();
		int ny = intensImg.getCols();
		int nz = intensImg.getSlices();
		int nt = intensImg.getComponents();
		int nxyz = nx*ny*nz;
		float rx = intensImg.getHeader().getDimResolutions()[0];
		float ry = intensImg.getHeader().getDimResolutions()[1];
		float rz = intensImg.getHeader().getDimResolutions()[2];
		System.out.println("Image dimensions: "+nx+" x "+ny+" x "+nz);
		System.out.println("Image resolutions: "+rx+" x "+ry+" x "+rz);
		
		float[][][] intensity3d = null;
		float[][][][] intensity4d = null;
		if (nt>1) intensity4d = intensImg.toArray4d();
		else intensity3d = intensImg.toArray3d();
		
		String imgname = intensImg.getName();
		intensImg = null;
		
		EmbeddedSurface inflated = null;
		if (inflatedSurface.getValue()!=null) inflated = inflatedSurface.getSurface();
		EmbeddedSurface orig = origSurface.getSurface();
		
		// check for (robust) data range; increase by a factor of x1000 if within [0,1]
		//float Imin = ImageStatistics.robustMinimum(intensity, 0.01f, 4, nx, ny, nz);
		//float Imax = ImageStatistics.robustMaximum(intensity, 0.01f, 4, nx, ny, nz);
		//System.out.println("Image range : ["+Imin+", "+Imax+"]");
		//float factor = 1.0f;
		//if (Imax-Imin < 10.0f) factor = 1000.0f;
		
		double[][] data = orig.getVertexData();

		byte mapStyle = NONE;
		if (mappingOption.getValue().equals("closest_point")) mapStyle = CLOSEST;
		else if (mappingOption.getValue().equals("linear_interp")) mapStyle = LINEAR;
		else if (mappingOption.getValue().equals("highest_value")) mapStyle = HIGHEST;
		
		for(int i=0; i<data.length; i++){
			Point3f p = orig.getVertex(i);
			
			// mapping from mesh to voxel space: multiple options!
			
			if (mapStyle==CLOSEST) {
				int x = Numerics.round(p.x/rx);
				int y = Numerics.round((ny-1)-p.y/ry);
				int z = Numerics.round((nz-1)-p.z/rz);
				if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
					//data[i] = new double[]{factor*Numerics.bounded(intensity[x][y][z],2*Imin-Imax,2*Imax-Imin)};
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++) data[i][t] = intensity4d[x][y][z][t];
					} else {
						data[i] = new double[]{intensity3d[x][y][z]};
					}
				} else {
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i] = new double[]{0.0f};
					}
				}
			} else if (mapStyle==LINEAR) {
				int x = Numerics.floor(p.x/rx);
				int y = Numerics.floor((ny-1)-p.y/ry);
				int z = Numerics.floor((nz-1)-p.z/rz);
				if (x>0 && x<nx-2 && y>0 && y<ny-2 && z>0 && z<nz-2) {
					float dx = p.x/rx - x;
					float dy = (ny-1)-p.y/ry - y;
					float dz = (nz-1)-p.z/rz - z;
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++) 						
							data[i][t] = (1-dx)*(1-dy)*(1-dz)*intensity4d[x][y][z][t]
										+dx*(1-dy)*(1-dz)*intensity4d[x+1][y][z][t]
										+(1-dx)*dy*(1-dz)*intensity4d[x][y+1][z][t]
										+(1-dx)*(1-dy)*dz*intensity4d[x][y][z+1][t]
										+dx*dy*(1-dz)*intensity4d[x+1][y+1][z][t]
										+(1-dx)*dy*dz*intensity4d[x][y+1][z+1][t]
										+dx*(1-dy)*dz*intensity4d[x+1][y][z+1][t]
												+dx*dy*dz*intensity4d[x+1][y+1][z+1][t];
					} else {
						data[i] = new double[]{(1-dx)*(1-dy)*(1-dz)*intensity3d[x][y][z]
												+dx*(1-dy)*(1-dz)*intensity3d[x+1][y][z]
												+(1-dx)*dy*(1-dz)*intensity3d[x][y+1][z]
												+(1-dx)*(1-dy)*dz*intensity3d[x][y][z+1]
												+dx*dy*(1-dz)*intensity3d[x+1][y+1][z]
												+(1-dx)*dy*dz*intensity3d[x][y+1][z+1]
												+dx*(1-dy)*dz*intensity3d[x+1][y][z+1]
												+dx*dy*dz*intensity3d[x+1][y+1][z+1]};
					}
				} else {
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i] = new double[]{0.0f};
					}
				}				
			} else if (mapStyle==HIGHEST) {
				int x = Numerics.floor(p.x/rx);
				int y = Numerics.floor((ny-1)-p.y/ry);
				int z = Numerics.floor((nz-1)-p.z/rz);
				if (x>0 && x<nx-2 && y>0 && y<ny-2 && z>0 && z<nz-2) {
					float dx = p.x/rx - x;
					float dy = (ny-1)-p.y/ry - y;
					float dz = (nz-1)-p.z/rz - z;
					
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++)
							data[i][t] = Numerics.max(intensity4d[x][y][z][t],
														intensity4d[x+1][y][z][t],
														intensity4d[x][y+1][z][t],
														intensity4d[x][y][z+1][t],
														intensity4d[x+1][y+1][z][t],
														intensity4d[x][y+1][z+1][t],
														intensity4d[x+1][y][z+1][t],
														intensity4d[x+1][y+1][z+1][t]);
					} else {
						data[i] = new double[]{Numerics.max(intensity3d[x][y][z],
															intensity3d[x+1][y][z],
															intensity3d[x][y+1][z],
															intensity3d[x][y][z+1],
															intensity3d[x+1][y+1][z],
															intensity3d[x][y+1][z+1],
															intensity3d[x+1][y][z+1],
															intensity3d[x+1][y+1][z+1])};
					}
				} else {
					if (nt>1) {
						data[i] = new double[nt];
						for (int t=0;t<nt;t++) data[i][t] = 0.0f;
					} else {
						data[i] = new double[]{0.0f};
					}
				}				
				
			}
		}
		// should we clone or not? not needed it seems
		orig.setVertexData(data);
		if (inflated!=null) inflated.setVertexData(data);
		
		orig.setName(orig.getName()+"_"+intensityName.getValue());
		if (inflated!=null) inflated.setName(inflated.getName()+"_"+intensityName.getValue());
		
		mappedOrgSurface.setValue(orig);
		if (inflated!=null) mappedInfSurface.setValue(inflated);		
	}


}
