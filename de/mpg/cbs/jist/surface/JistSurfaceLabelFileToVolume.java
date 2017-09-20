package de.mpg.cbs.jist.surface;

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
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamBoolean;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;

import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;

import gov.nih.mipav.view.MipavUtil;

import de.mpg.cbs.utilities.Numerics;
import de.mpg.cbs.utilities.Interface;
import de.mpg.cbs.utilities.BasicInfo;
import de.mpg.cbs.structures.BinaryHeap4D;

import java.io.*;
import java.util.*;

import javax.vecmath.Point3f;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.structures.geom.EmbeddedSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;


public class JistSurfaceLabelFileToVolume  extends ProcessingAlgorithm{
	ParamVolume originalImage;
	ParamSurface originalSurface;
	ParamFile labelFile;
	ParamSurface dataSurface;
	ParamVolume dataImage;
	ParamVolume maskImage;
	ParamBoolean mipavTransform;
	ParamFloat labelExtension;
	
	protected void createInputParameters(ParamCollection inputParams) {
		inputParams.add(originalSurface=new ParamSurface("Surface Mesh"));
		inputParams.add(originalImage=new ParamVolume("Reference Volume"));
		originalImage.setDescription("Reference volume to use for level set representation dimensions.");
		inputParams.add(labelFile=new ParamFile("Label File"));
		labelFile.setDescription("File of vertex labels (e.g. from FreeSurfer)");
		inputParams.add(mipavTransform=new ParamBoolean("Align to MIPAV image space", true));
		inputParams.add(labelExtension=new ParamFloat("Label Extension Distance", 0.0f, 100.0f, 5.0f));
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Surfaces");
		inputParams.setName("LabelFileToVolume");
		inputParams.setLabel("Mesh Label File to Volume");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Embeds separate labels linked to an input surface mesh (e.g. from Freesurfer labels) into volumetric space.");
		
		info.setVersion("3.0");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(dataSurface=new ParamSurface("Labels Surface"));
		outputParams.add(dataImage=new ParamVolume("Labels Volume"));
		outputParams.add(maskImage=new ParamVolume("Mask Volume"));
	}
	
	protected void execute(CalculationMonitor monitor) {
		
		BasicInfo.displayMessage("load image\n");		

		ImageDataFloat	originalImgData = new ImageDataFloat(originalImage.getImageData());
		int nx = originalImgData.getRows();
		int ny = originalImgData.getCols();
		int nz = originalImgData.getSlices();
		float rx = originalImgData.getHeader().getDimResolutions()[0];
		float ry = originalImgData.getHeader().getDimResolutions()[1];
		float rz = originalImgData.getHeader().getDimResolutions()[2];
		//float ox = originalImgData.getHeader().getOrigin()[0];
		//float oy = originalImgData.getHeader().getOrigin()[1];
		//float oz = originalImgData.getHeader().getOrigin()[2];
		
		BasicInfo.displayMessage("load surface\n");		

		EmbeddedSurface surf = originalSurface.getSurface();
		
		BasicInfo.displayMessage("map surface coordinates\n");		

		// maps the surface to the voxel space
		if (mipavTransform.getValue().booleanValue()) {
			surf.scaleVertices(new float[]{-1.0f, 1.0f, -1.0f});
			Point3f pt0 = new Point3f((nx-1)*rx/2.0f, (ny-1)*ry/2.0f, (nz-1)*rz/2.0f);
			surf.translate(pt0);
		}

		float[][][] labels = new float[nx][ny][nz];
		byte[][][] mask = new byte[nx][ny][nz];
		double[][] data = new double[surf.getVertexCount()][1];

		BasicInfo.displayMessage("load labels\n");		

		float[] list = getLabelsFromFreesurferFile(labelFile.getValue(), surf.getVertexCount()); 
		
		BasicInfo.displayMessage("map surface labels to image\n");		

		for(int i=0; i<surf.getVertexCount(); i++){
			// labels on surface
			data[i][0] = list[i];
			
			// mapping from mesh to voxel space
			Point3f p = surf.getVertex(i);
			
			int x = Numerics.round(p.x/rx);
			int y = Numerics.round((ny-1)-p.y/ry);
			int z = Numerics.round((nz-1)-p.z/rz);
			if (x>0 && x<nx-1 && y>0 && y<ny-1 && z>0 && z<nz-1) {
				labels[x][y][z] = list[i];
				mask[x][y][z] = (byte)1;
			}
		}
		surf.setVertexData(data);
		
		BasicInfo.displayMessage("propagate surface labels\n");		

		// expand labels out for a certain distance
		if (labelExtension.getValue().floatValue()>0)
			fastMarchingPropagation(labels, mask, labelExtension.getValue().floatValue(), nx, ny, nz);
		
		BasicInfo.displayMessage("generate outputs\n");		

		surf.setName(surf.getName()+"_labels");
		dataSurface.setValue(surf);
		
		ImageDataFloat dataImg = new ImageDataFloat(labels);		
		dataImg.setHeader(originalImage.getImageData().getHeader());
		dataImg.setName(surf.getName()+"_labels");
		dataImage.setValue(dataImg);
		dataImg = null;
		labels = null;
		
		ImageDataUByte maskImg = new ImageDataUByte(mask);		
		maskImg.setHeader(originalImage.getImageData().getHeader());
		maskImg.setName(surf.getName()+"_mask");
		maskImage.setValue(maskImg);
		maskImg = null;
		mask = null;
		
	}
	
	private final float[] getLabelsFromFreesurferFile(File input, int npt) {
		try {
            System.out.println("reading label file: "+input.getName());
            FileReader fr = new FileReader(input);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
			
			// Exact corresponding template for first line ?
            if (!line.startsWith("#!ascii label")) {
                System.out.println("not a proper Freesurfer label file");
                br.close();
                fr.close();
                return null;
            }
            // store the results: every vertex point
			float[] labels = new float[npt];
			
			 // label number
			line = br.readLine();
			StringTokenizer st = new StringTokenizer(line, " ");
			int nlabel = MipavUtil.getInt(st);
			for (int n=0;n<nlabel;n++) {
				line = br.readLine();
				st = new StringTokenizer(line, " ");
				// format: id X Y Z lb
				int id = MipavUtil.getInt(st);
				st.nextToken();
				st.nextToken();
				st.nextToken();
				float lb = MipavUtil.getFloat(st);
				//System.out.print("("+id+","+lb+")");
				
				labels[id] = lb;
			}
			br.close();
            fr.close();
            
            return labels;
        }
        catch (FileNotFoundException e) {
            System.out.println(e.getMessage());
        }
        catch (IOException e) {
            System.out.println(e.getMessage());
        } 
		catch (OutOfMemoryError e){
			System.out.println(e.getMessage());
		}
		catch (Exception e) {
			System.out.println(e.getMessage());
        }
        return null;
	}

	public final void fastMarchingPropagation(float[][][] data, byte[][][] mask, float maxDist, int nx, int ny, int nz) {
		float	UNKNOWN = -1.0f;
        
		// computation variables
       float[][][] distance = new float[nx][ny][nz];
		float[] nbdist = new float[6];
		boolean[] nbflag = new boolean[6];
		float curdist,newdist;	
		boolean done, isprocessed;
		
		// compute the neighboring labels and corresponding distance functions (! not the MGDM functions !)
		BasicInfo.displayMessage("fast marching\n");		

		int[] xoff = new int[]{1, -1, 0, 0, 0, 0};
		int[] yoff = new int[]{0, 0, 1, -1, 0, 0};
		int[] zoff = new int[]{0, 0, 0, 0, 1, -1};

		BinaryHeap4D heap = new BinaryHeap4D(nx+ny+nz, BinaryHeap4D.MINTREE);
		
        heap.reset();
		
        BasicInfo.displayMessage("init\n");		
        // initialize the heap from boundaries
        for (int x=0; x<nx; x++) for (int y=0; y<ny; y++) for (int z=0; z<nz; z++) if (mask[x][y][z]>0) {
        	distance[x][y][z] = 0.0f;
        	// search for boundaries
        	for (int n = 0; n<6; n++) {
				int xn = x+xoff[n];
				int yn = y+yoff[n];
				int zn = z+zoff[n];
				if (mask[xn][yn][zn]==0) {
					// add to the heap with previous value
					heap.addValue(1.0f,xn,yn,zn,(int)data[x][y][z]);
				}
            }
        }
		BasicInfo.displayMessage("main loop\n");		

        // grow the labels and functions
        while (heap.isNotEmpty()) {
        	//BasicInfo.displayMessage(".");
        	
        	// extract point with minimum distance
        	curdist = heap.getFirst();
        	int x = heap.getFirstX();
        	int y = heap.getFirstY();
        	int z = heap.getFirstZ();
        	int lb = heap.getFirstK();
        	heap.removeFirst();

			// if more than nlb labels have been found already, this is done
			if (mask[x][y][z]>0)  continue;
			
			// update the distance functions at the current level
			
			data[x][y][z] = lb;
			distance[x][y][z] = curdist;
			mask[x][y][z] = 1; // update the current level
 			
			// find new neighbors
			for (int n = 0; n<6; n++) {
				int xn = x+xoff[n];
				int yn = y+yoff[n];
				int zn = z+zoff[n];
				if (xn>1 && xn<nx-2 && yn>1 && yn<ny-2 && zn>1 && zn<nz-2) {
					// must be in outside the object or its processed neighborhood
					if (mask[xn][yn][zn]==0) {
						// compute new distance based on processed neighbors for the same object
						for (int l=0; l<6; l++) {
							nbdist[l] = UNKNOWN;
							nbflag[l] = false;
							int xb = xn + xoff[l];
							int yb = yn + yoff[l];
							int zb = zn + zoff[l];
							// note that there is at most one value used here
							if (data[xb][yb][zb]==lb && mask[xb][yb][zb]>0) {
								nbdist[l] = distance[xb][yb][zb];
								nbflag[l] = true;
							}			
						}
						newdist = minimumMarchingDistance(nbdist, nbflag);
					
						if (newdist<=maxDist) {
							// add to the heap
							heap.addValue(newdist,xn,yn,zn,lb);
						}
					}
				}
			}
		}

       return;
    }

    
	/**
     * the Fast marching distance computation 
     * (!assumes a 6D array with opposite coordinates stacked one after the other)
     * 
     */
    public final float minimumMarchingDistance(float[] val, boolean[] flag) {
		double s, s2, tmp; 
		int count;
		double dist;
    
        // s = a + b +c; s2 = a*a + b*b +c*c
        s = 0;
        s2 = 0;
        count = 0;

        for (int n=0; n<6; n+=2) {
			if (flag[n] && flag[n+1]) {
				tmp = Numerics.min(val[n], val[n+1]); // Take the smaller one if both are processed
				s += tmp;
				s2 += tmp*tmp;
				count++;
			} else if (flag[n]) {
				s += val[n]; // Else, take the processed one
				s2 += val[n]*val[n];
				count++;
			} else if (flag[n+1]) {
				s += val[n+1];
				s2 += val[n+1]*val[n+1];
				count++;
			}
		}
         // count must be greater than zero since there must be at least one processed pt in the neighbors
        
        tmp = (s+Math.sqrt( (s*s-count*(s2-1.0f))))/count;

        // The larger root
        return (float)tmp;
    }

}
