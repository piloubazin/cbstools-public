package de.mpg.cbs.jist.surface;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamOption;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFloat;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamInteger;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamObject;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamFile;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamSurface;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.geom.CurveCollection;
import edu.jhu.ece.iacl.jist.structures.geom.CurvePath;
import edu.jhu.ece.iacl.jist.io.CurveVtkReaderWriter;
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

import javax.vecmath.Point3f;
import org.apache.commons.math3.util.FastMath;


/*
 * @author Pierre-Louis Bazin
 */
public class JistSurfaceMeshDataTree extends ProcessingAlgorithm {

	// jist containers
	private ParamSurface	inputSurface;

	private ParamOption		outputType;
	private static final String[] outputTypes = {"increasing", "decreasing", "positive", "negative"};

	private ParamInteger stepParam;
	private ParamFloat factorParam;
	
	
	private	ParamObject<CurveCollection>	meshtreeLines;
	
	// global variables
	private static final byte X = 0;
	private static final byte Y = 1;
	private static final byte Z = 2;
	private static final byte N = 3;

	private int	nx, ny, nz;
	
	protected void createInputParameters(ParamCollection inputParams) {
		
		inputParams.add(inputSurface = new ParamSurface("Input Surface"));
		
		inputParams.add(outputType = new ParamOption("value tree", outputTypes));
		
		inputParams.add(stepParam=new ParamInteger("Decimation",1,100000,10));
		inputParams.add(factorParam=new ParamFloat("Scaling factor",0.0f,10.0f,0.5f));
		
		inputParams.setPackage("CBS");
		inputParams.setCategory("Surface");
		inputParams.setLabel("Mesh Data Tree");
		inputParams.setName("MeshDataTree");

		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new AlgorithmAuthor("Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/"));
		info.setAffiliation("Max Planck Institute for Human Cognitive and Brain Sciences");
		info.setDescription("Compute a tree representation of surface values.");
		
		info.setVersion("3.1.3");
		info.setStatus(DevelopmentStatus.RC);
		info.setEditable(false);
	}

	@Override
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(meshtreeLines=new ParamObject<CurveCollection>("Tree mesh (VTK)",new CurveVtkReaderWriter()));		
	}

	@Override
	protected void execute(CalculationMonitor monitor) {	
	    //executeBranchToCenter(monitor);
	    executeCenterToBranch(monitor);
	}
	    
	private void executeBranchToCenter(CalculationMonitor monitor){
		
		EmbeddedSurface surf = inputSurface.getSurface();
		
		int steps = stepParam.getValue().intValue();

		BasicInfo.displayMessage("Building neighbor table\n");
		
		surf.buildAllTables();
		int[][] neighbors = surf.getNeighborVertexVertexTable();
		double[][] data = surf.getVertexData();
		int npts = data.length;
		
		double min = 1e15, max = -1e15;
		for (int i=0; i<npts; i++) {
		    if (data[i][0]<min) min = data[i][0];
		    if (data[i][0]>max) max = data[i][0];
		}
		
        double[] origin = new double[4];
        for  (int p=0;p<npts;p++) {
            Point3f pt = surf.getVertex(p);
            origin[X] += pt.x;
            origin[Y] += pt.y;
            origin[Z] += pt.z;
            origin[N] += 1.0;
        }
        origin[X] /= origin[N];
        origin[Y] /= origin[N];
        origin[Z] /= origin[N];


		BasicInfo.displayMessage("Data range: ["+min+", "+max+"], "+steps+" steps\n");
		
		int[] labels = new int[npts];    
		int[] prev = new int[npts];
		int nblb = 0;
		BinaryHeap2D heap = new BinaryHeap2D(npts, BinaryHeap2D.MAXTREE);
		CurveCollection meshtree = new CurveCollection();
		for (int s=0;s<steps;s++) {
		    double low, high;
		    if (outputType.getValue().equals("increasing")) {
		        high = min + (s+1)*(max-min)/steps;
		        low = min - 0.001*(max-min);
		    } else if (outputType.getValue().equals("decreasing")) {
		        high = max;
		        low = max - (s+1)*(max-min)/steps;
		    } else if (outputType.getValue().equals("negative")) {
		        high = min + (s+1)*(0-min)/steps;
		        low = min - 0.001*(0-min);
		    } else if (outputType.getValue().equals("positive")) {
		        high = max;
		        low = max - (s+1)*(max-0)/steps;
		    } else {
		        high = max;
		        low = 0.999*min;
		    }
		    BasicInfo.displayMessage("step "+s+" range: ["+low+", "+high+"]\n");
		    
		    // find and label regions within this 
		    BasicInfo.displayMessage("labeling ");
		    for (int p=0;p<npts;p++) {
		        if (s>0) prev[p] = labels[p];
		        labels[p] = 0;
		    }
		    heap.reset();
		    int newlb = 1;
		    for (int p=0;p<npts;p++) {
		        if (labels[p]==0 && (data[p][0]>low && data[p][0]<=high)) {
		            
		            BasicInfo.displayMessage(".");
		            // propagate to all valid neighbors
		            heap.addValue(0.0f, p, (byte)1);
		            
		            while (heap.isNotEmpty()) {
		                float dist = heap.getFirst();
		                int pt = heap.getFirstId();
		                heap.removeFirst();

		                labels[pt] = newlb;
		                for (int n=0;n<neighbors[pt].length;n++) {
		                    int ngb = neighbors[pt][n];
		                    if (labels[ngb]==0 && (data[ngb][0]>low && data[ngb][0]<=high)) {
		                        heap.addValue(dist+1.0f, ngb, (byte)1);
		                    }
		                }
		            }
		            newlb++;
		        }
		    }
		    BasicInfo.displayMessage(" done\n");
		    
		    // compute mean location for each cluster
		    BasicInfo.displayMessage("centroids \n");
		    double[][] centroids = new double[newlb][4];
		    for  (int p=0;p<npts;p++) if (labels[p]>0) {
		        Point3f pt = surf.getVertex(p);
                centroids[labels[p]][X] += pt.x;
                centroids[labels[p]][Y] += pt.y;
                centroids[labels[p]][Z] += pt.z;
                centroids[labels[p]][N] += 1.0;
		   }
		   for (int l=1;l<newlb;l++) {
		       centroids[l][X] /= centroids[l][N];
		       centroids[l][Y] /= centroids[l][N];
		       centroids[l][Z] /= centroids[l][N];
		   }
		   // weight centroids by distance to origin
		   double wdist = (double)s/(double)steps;
		   for (int l=1;l<newlb;l++) {
		       centroids[l][X] = (1.0-wdist)*centroids[l][X] + wdist*origin[X];
		       centroids[l][Y] = (1.0-wdist)*centroids[l][Y] + wdist*origin[Y];
		       centroids[l][Z] = (1.0-wdist)*centroids[l][Z] + wdist*origin[Z];
		   }
		   
		   // add to a growing graph structure
		   BasicInfo.displayMessage("add to tree ");
		   if (s==0) {
		       for (int l=1;l<newlb;l++) {
                   CurvePath curve = new CurvePath();
                   curve.add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
                   meshtree.add(curve);
               }
               nblb = newlb;
               BasicInfo.displayMessage((nblb-1)+" branches \n");
		   } else {
		       // find the maximum link in previous: no, we need to build the full connectivity matrix..
		       boolean[][] link = new boolean[nblb][newlb];
		       for (int p=0;p<npts;p++) if (labels[p]>0 && prev[p]>0) {
		           link[prev[p]][labels[p]] = true;
		       }
		       // relabel: old label if connected, new label otherwise
		       int[] relabel = new int[newlb];
		       int nblbnew = nblb;
		       for (int l=1;l<newlb;l++) {
		           boolean isnew = true;
		           for (int n=0;n<nblb;n++) if (link[n][l]) {
		               isnew = false;
		               ((CurvePath)meshtree.get(n-1)).add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
		               relabel[l] = n;
		           }
		           if (isnew) {
		               CurvePath curve = new CurvePath();
                       curve.add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
                       meshtree.add(curve);
                       relabel[l] = nblbnew;
                       nblbnew++;
                   }
		       }
		       nblb = nblbnew;
               BasicInfo.displayMessage((nblbnew-nblb)+" branches \n");

               // update the global number of labels
               BasicInfo.displayMessage("update labels \n");
               for (int p=0;p<npts;p++) if (labels[p]>0) {
                   labels[p] = relabel[labels[p]];
               }
		   }
		   
		}
		   
		// output
		meshtree.setName(surf.getName()+"_mtree");
		meshtreeLines.setObject(meshtree);
		
		return;
	}

	private void executeCenterToBranch(CalculationMonitor monitor){
		
		EmbeddedSurface surf = inputSurface.getSurface();
		
		int steps = stepParam.getValue().intValue();
		
		float factor = factorParam.getValue().floatValue();

		BasicInfo.displayMessage("Building neighbor table\n");
		
		surf.buildAllTables();
		int[][] neighbors = surf.getNeighborVertexVertexTable();
		double[][] data = surf.getVertexData();
		int npts = data.length;
		
		double min = 1e15, max = -1e15;
		for (int i=0; i<npts; i++) {
		    if (data[i][0]<min) min = data[i][0];
		    if (data[i][0]>max) max = data[i][0];
		}
		
        double[] origin = new double[4];
        for  (int p=0;p<npts;p++) {
            Point3f pt = surf.getVertex(p);
            origin[X] += pt.x;
            origin[Y] += pt.y;
            origin[Z] += pt.z;
            origin[N] += 1.0;
        }
        origin[X] /= origin[N];
        origin[Y] /= origin[N];
        origin[Z] /= origin[N];


		BasicInfo.displayMessage("Data range: ["+min+", "+max+"], "+steps+" steps\n");
		
		int[] labels = new int[npts];    
		int[] prev = new int[npts];
		int nblb = 0;
		BinaryHeap2D heap = new BinaryHeap2D(npts, BinaryHeap2D.MAXTREE);
		CurveCollection meshtree = new CurveCollection();
		for (int s=0;s<steps;s++) {
		    double low, high;
		    if (outputType.getValue().equals("increasing")) {
		        high = min + (steps-s)*(max-min)/steps;
		        low = min - 0.001*(max-min);
		    } else if (outputType.getValue().equals("decreasing")) {
		        high = max;
		        low = max - (steps-s)*(max-min)/steps;
		    } else if (outputType.getValue().equals("negative")) {
		        high = min + (steps-s)*(0-min)/steps;
		        low = min - 0.001*(0-min);
		    } else if (outputType.getValue().equals("positive")) {
		        high = max;
		        low = max - (steps-s)*(max-0)/steps;
		    } else {
		        high = max;
		        low = 0.999*min;
		    }
		    BasicInfo.displayMessage("step "+s+" range: ["+low+", "+high+"]\n");
		    
		    // find and label regions within this 
		    BasicInfo.displayMessage("labeling ");
		    for (int p=0;p<npts;p++) {
		        if (s>0) prev[p] = labels[p];
		        labels[p] = 0;
		    }
		    heap.reset();
		    int newlb = 1;
		    for (int p=0;p<npts;p++) {
		        if (labels[p]==0 && (data[p][0]>low && data[p][0]<=high)) {
		            
		            BasicInfo.displayMessage(".");
		            // propagate to all valid neighbors
		            heap.addValue(0.0f, p, (byte)1);
		            
		            while (heap.isNotEmpty()) {
		                float dist = heap.getFirst();
		                int pt = heap.getFirstId();
		                heap.removeFirst();

		                labels[pt] = newlb;
		                for (int n=0;n<neighbors[pt].length;n++) {
		                    int ngb = neighbors[pt][n];
		                    if (labels[ngb]==0 && (data[ngb][0]>low && data[ngb][0]<=high)) {
		                        heap.addValue(dist+1.0f, ngb, (byte)1);
		                    }
		                }
		            }
		            newlb++;
		        }
		    }
		    BasicInfo.displayMessage((newlb-1)+" branches\n");
		    
		    // compute mean location for each cluster
		    BasicInfo.displayMessage("centroids \n");
		    double[][] centroids = new double[newlb][4];
		    for  (int p=0;p<npts;p++) if (labels[p]>0) {
		        Point3f pt = surf.getVertex(p);
                centroids[labels[p]][X] += pt.x;
                centroids[labels[p]][Y] += pt.y;
                centroids[labels[p]][Z] += pt.z;
                centroids[labels[p]][N] += 1.0;
		   }
		   for (int l=1;l<newlb;l++) {
		       centroids[l][X] /= centroids[l][N];
		       centroids[l][Y] /= centroids[l][N];
		       centroids[l][Z] /= centroids[l][N];
		   }
		   /*
		   // weight centroids by distance to origin
		   double wdist = (double)s/(double)steps;
		   for (int l=1;l<newlb;l++) {
		       centroids[l][X] = (wdist)*centroids[l][X] + (1.0-wdist)*origin[X];
		       centroids[l][Y] = (wdist)*centroids[l][Y] + (1.0-wdist)*origin[Y];
		       centroids[l][Z] = (wdist)*centroids[l][Z] + (1.0-wdist)*origin[Z];
		   }
		   */
		   // add to a growing graph structure
		   BasicInfo.displayMessage("add to tree ");
		   if (s==0) {
		       for (int l=1;l<newlb;l++) {
                   CurvePath curve = new CurvePath();
                   curve.add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
                   meshtree.add(curve);
               }
               nblb = newlb;
               BasicInfo.displayMessage((nblb-1)+" branches \n");
		   } else {
		       // find the maximum link in previous: no, we need to build the full connectivity matrix..
		       boolean[][] link = new boolean[nblb][newlb];
		       for (int p=0;p<npts;p++) if (labels[p]>0 && prev[p]>0) {
		           link[prev[p]][labels[p]] = true;
		       }
		       // relabel: old label if first connected, new label otherwise
		       int[] relabel = new int[newlb];
		       int nblbnew = nblb;
		       for (int n=0;n<nblb;n++) {
		           boolean islinked = false;
		           for (int l=1;l<newlb;l++) {
                       if (link[n][l]) {
                           BasicInfo.displayMessage("+");
                           Point3f[] pts = ((CurvePath)meshtree.get(n-1)).getCurve();
                           Point3f pt = pts[pts.length-1];
                           // weight centroids by distance to previous
                           double wdist = FastMath.pow((double)s/(double)steps, factor);
                           //centroids[l][X] = (wdist)*centroids[l][X] + (1.0-wdist)*pt.x;
                           //centroids[l][Y] = (wdist)*centroids[l][Y] + (1.0-wdist)*pt.y;
                           //centroids[l][Z] = (wdist)*centroids[l][Z] + (1.0-wdist)*pt.z;
                           centroids[l][X] = (wdist)*centroids[l][X] + (1.0-wdist)*origin[X];
                           centroids[l][Y] = (wdist)*centroids[l][Y] + (1.0-wdist)*origin[Y];
                           centroids[l][Z] = (wdist)*centroids[l][Z] + (1.0-wdist)*origin[Z];
 
                           // get the previous value?
                           if (!islinked) {
                               BasicInfo.displayMessage("x");
                               ((CurvePath)meshtree.get(n-1)).add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
                               relabel[l] = n;
                               islinked = true;
                           } else {
                               BasicInfo.displayMessage("-");
                               CurvePath curve = new CurvePath();
                               // add previous centroid TODO
                               curve.add(pt);
                               // continue
                               curve.add(new Point3f((float)centroids[l][X], (float)centroids[l][Y], (float)centroids[l][Z]));
                               meshtree.add(curve);
                               relabel[l] = nblbnew;
                               nblbnew++;
                           }
                       }
                   }
		       }
		       BasicInfo.displayMessage((nblbnew-nblb)+" branches \n");
		       nblb = nblbnew;
               
               // update the global number of labels
               BasicInfo.displayMessage("update labels \n");
               for (int p=0;p<npts;p++) if (labels[p]>0) {
                   labels[p] = relabel[labels[p]];
               }
		   }
		   
		}
		   
		// output
		meshtree.setName(surf.getName()+"_mtree");
		meshtreeLines.setObject(meshtree);
		
		return;
	}

}
