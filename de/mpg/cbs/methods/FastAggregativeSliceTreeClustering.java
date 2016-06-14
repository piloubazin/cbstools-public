package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import gov.nih.mipav.view.*;
import gov.nih.mipav.model.structures.jama.*;
import gov.nih.mipav.model.file.FileInfoBase;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.util.FastMath;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.structures.*;
import de.mpg.cbs.utilities.*;

/**
 *
 *  This algorithm builds a hierarchical clustering of connectivity matrix data
 *	The connectivity is assumed to be always >0; a value of -1 makes no connection
 *
 *	The key component is an array structure that contains the association weights for all pairs of connected points.
 *	cluster as well as the degree in the first element of the list. Once constructed, the association array is fully
 *	abstract (no contextual information remains).
 *
 *	@version    September 2012
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class FastAggregativeSliceTreeClustering {
		
	// data buffers
	private 	float[][][]	image;  			// original data
	private 	int			nx,ny,nz;   		// image dimensions
	private     String		slicedir;
	private 	int[] 		slicecount;
	
	private		float		scale;
	private		float		basis;
	private		float 		pvalue;
	private		float		minsize;
	
	private		byte		metric;
	private	static final	byte	DICE = 10;
	private	static final	byte	JACCARD = 11;
	private	static final	byte	OVERLAP_SIZE = 12;
	private	static final	byte	MIX = 13;
	private	static final	byte	RATIO = 14;
	
	private		byte		mode;
	private	static final	byte	DIRECT = 20;
	private	static final	byte	GAUSS = 21;
	private	static final	byte	NP = 22;
	
	private static final double PI2 = Math.PI/2.0;
	private static final double SQPI = Math.sqrt(Math.PI);
	private static final float INF = 1e12f;
	
	// algorithm quantities
	private		int[]		clustering;
	private		int			nlb;
	private		int			key;
	
	//private		ArrayList<Float>	degree;
	//private		ArrayList<ArrayList<Triple>>	assoc;
	private		Triple[][]						assoc;
	private		BitSet							active;
	private		BinaryHeapPair					bintree;
	private		float[]							cost;
	private		float[]							self;
	private		float[]							other;
	private		int[]							latest;
	private		int[][]							clusterPos;
	private		float[]							boundaryScore;
	
	private		Histogram		distribution;
	
    static final boolean		debug				=	true;
	static final boolean		verbose				=	true;
    
	// simple id, value class
	private static class Triple {
		public int id;
		public float weight;
		public float delta;
		public float size;
		
		public Triple(int id) {
			this.id = id;
			this.weight = 0.0f;
			this.delta = 0.0f;
			this.size = 1;
		}
		public Triple(int id, float weight) {
			this.id = id;
			this.weight = weight;
			this.delta = 0.0f;
			this.size = 1;
		}
		public Triple(int id, float weight, float size) {
			this.id = id;
			this.weight = weight;
			this.delta = 0.0f;
			this.size = size;
		}
		public Triple(int id, float weight, float delta, float size) {
			this.id = id;
			this.weight = weight;
			this.delta = delta;
			this.size = size;
		}
		
	}
	
	/* simple variable size array */
	private static class IntArray {
		public int[] val;
		public int	size;
		public int	last;
		
		public IntArray(int size) {
			this.val = new int[size];
			this.size = size;
			this.last = 0;
		}
		public void add(int id) {
			if (last>=val.length) {
				int[] tmp = new int[val.length+size];
				for (int n=0;n<last;n++) tmp[n] = val[n];
				val = tmp;
				if (debug) System.out.print("#");
			}
			val[last] = id;
			last++;
		}
		public void reset() {
			last = 0;
		}
	}
	

	/**
	 *  constructor
	 */
	public FastAggregativeSliceTreeClustering(float[][][] img_, int nx_, int ny_, int nz_, String dir_, 
												float sca_, float bas_, String met_, String mod_) {
		image = img_;
		
		nx = nx_;
		ny = ny_;
		nz = nz_;
		slicedir = dir_;
		
		scale = sca_;
		basis = bas_;
		pvalue = sca_;
		
		if (met_.equals("Dice")) 				metric = DICE;
		else if (met_.equals("Jaccard")) 		metric =JACCARD;
		else if (met_.equals("overlap_size")) 	metric = OVERLAP_SIZE;
		else if (met_.equals("mix")) 			metric = MIX;
		else if (met_.equals("ratio")) 			metric = RATIO;
		
		if (mod_.equals("direct")) 					mode = DIRECT;
		else if (mod_.equals("gaussian")) 			mode =GAUSS;
		else if (mod_.equals("non_parametric")) 	mode = NP;
		
		//estimate the min size of overlap parameter
		minsize = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (image[x][y][z]!=0) minsize++;
		}
		if (slicedir.equals("X")) minsize /= nx;
		if (slicedir.equals("Y")) minsize /= ny;
		if (slicedir.equals("Z")) minsize /= nz;
		// we assume basis% of the average slice to be a good minimum size
		minsize = scale*minsize;
		
		// build a unique list of ids for all labels 
		slicecount = countAndRelabel(slicedir);
	}

	final public void finalize() {
		image = null;
		assoc = null;
		System.gc();
	}
	
   public final int[] countAndRelabel(String dir) {
        int x,y,z;
        int Nlb;
        ArrayList<Float> lb;
        boolean newLabel;
        int[] count = null;
        
        Nlb = 0;
		int prev = 0;
		lb = new ArrayList();
		if (dir.equals("X")) {
			count = new int[nx];
			for (x=0;x<nx;x++) {
        		for (y=0;y<ny;y++) for (z=0;z<nz;z++) {
        			if (image[x][y][z]!=0) {
						newLabel=true;
						for (int n=prev;n<Nlb;n++) {
							if (image[x][y][z]==lb.get(n)) { 
								newLabel=false;
								// replace the image value
								image[x][y][z] = n;
								break; 
							}
						}
						if (newLabel) {
							lb.add(Nlb,image[x][y][z]);
							image[x][y][z] = Nlb;
							Nlb++;
						}
					}
				}
				count[x] = Nlb-prev;
				prev = Nlb;
			}
		}
		if (dir.equals("Y")) {
			count = new int[ny];
			for (y=0;y<ny;y++) {
        		for (z=0;z<nz;z++) for (x=0;x<nx;x++) {
        			if (image[x][y][z]!=0) {
						newLabel=true;
						for (int n=prev;n<Nlb;n++) {
							if (image[x][y][z]==lb.get(n)) { 
								newLabel=false;
								// replace the image value
								image[x][y][z] = n+1;
								break; 
							}
						}
						if (newLabel) {
							lb.add(Nlb,image[x][y][z]);
							image[x][y][z] = Nlb+1;
							Nlb++;
						}
					}
				}
				count[y] = Nlb-prev;
				prev = Nlb;
			}
		}
		if (dir.equals("Z")) {
			count = new int[nz];
			for (z=0;z<nz;z++) {
        		for (x=0;x<nx;x++) for (y=0;y<ny;y++) {
        			if (image[x][y][z]!=0) {
						newLabel=true;
						for (int n=prev;n<Nlb;n++) {
							if (image[x][y][z]==lb.get(n)) { 
								newLabel=false;
								// replace the image value
								image[x][y][z] = n;
								break; 
							}
						}
						if (newLabel) {
							lb.add(Nlb,image[x][y][z]);
							image[x][y][z] = Nlb;
							Nlb++;
						}
					}
				}
				count[z] = Nlb-prev;
				prev = Nlb;
			}
		}
        return count;
    }
    
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final float[] computeWeightMetrics() {
		
		if (debug) System.out.println("-- metric pre-processing --");
		
		nlb=0;
		float nconnect = 0;
		for (int n=0;n<slicecount.length;n++) {
			nlb += slicecount[n];
			if (n>0) nconnect += slicecount[n]*slicecount[n-1];
		}
		if (debug) System.out.println("total number of labels: "+nlb+", max number of connections: "+nconnect);
		
		// too big to store: buid the histogram right away??
		float[] distance = new float[(int)nconnect];
		int[] weight = new int[(int)nconnect];
		int nb = 0;
		
		if (slicedir.equals("X")) {
			int offset = 1+slicecount[0];
			for (int x=1;x<nx-1;x++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[x]>0  && slicecount[x-1]>0) {
					float[] vol0 = new float[slicecount[x]];
					float[] voln = new float[slicecount[x-1]];
					float[][] assoc0n = new float[slicecount[x]][slicecount[x-1]];
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x-1][y][z]!=0) {
							lm = (int)(image[x-1][y][z]-offset+slicecount[x-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x-1][y][z]!=0) assoc0n[l0][lm]++;
					}
					// build all the distances based on overlap
					for (int l0=0;l0<slicecount[x];l0++) {
						// build the distance array
						for (int ln=0;ln<slicecount[x-1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							distance[nb] = dataDistance(assoc0n[l0][ln], vol0[l0], voln[ln]);
							weight[nb] = (int)assoc0n[l0][ln];
							nb++;								
						}
					}
				}
				offset += slicecount[x];
			}
		} else
		if (slicedir.equals("Y")) {
			int offset = 1+slicecount[0];
			for (int y=1;y<ny-1;y++) {
				//System.out.println("slice "+y+" ("+slicecount[y-1]+", "+slicecount[y]+", "+slicecount[y+1]+")");
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[y]>0  && slicecount[y-1]>0) {
					float[] vol0 = new float[slicecount[y]];
					float[] voln = new float[slicecount[y-1]];
					float[][] assoc0n = new float[slicecount[y]][slicecount[y-1]];
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y-1][z]!=0) {
							lm = (int)(image[x][y-1][z]-offset+slicecount[y-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x][y-1][z]!=0) assoc0n[l0][lm]++;
					}
					// build all the weights based on overlap
					for (int l0=0;l0<slicecount[y];l0++) {
						// build the weight array
						for (int ln=0;ln<slicecount[y-1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							distance[nb] = dataDistance(assoc0n[l0][ln], vol0[l0], voln[ln]);
							weight[nb] = (int)assoc0n[l0][ln];
							nb++;								
						}
					}
				}
				offset += slicecount[y];
			}
		} else
		if (slicedir.equals("Z")) {
			int offset = 1+slicecount[0];
			for (int z=1;z<nz-1;z++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[z]>0  && slicecount[z-1]>0) {
					float[] vol0 = new float[slicecount[z]];
					float[] voln = new float[slicecount[z-1]];
					float[][] assoc0n = new float[slicecount[z]][slicecount[z-1]];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y][z-1]!=0) {
							lm = (int)(image[x][y][z-1]-offset+slicecount[z-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x][y][z-1]!=0) assoc0n[l0][lm]++;
					}
					// build all the weights based on overlap
					for (int l0=0;l0<slicecount[z];l0++) {
						// build the weight array
						for (int ln=0;ln<slicecount[z-1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							distance[nb] = dataDistance(assoc0n[l0][ln], vol0[l0], voln[ln]);
							weight[nb] = (int)assoc0n[l0][ln];
							nb++;								
						}
					}
				}
				offset += slicecount[z];
			}
		}
		
		if (debug) System.out.println("statistics");

		Histogram hist = new Histogram(distance, weight, 100, nb);
		System.out.println("Mean : "+hist.mean());
		System.out.println("Stdev : "+hist.stdev());
		System.out.println("Min : "+hist.min());
		System.out.println("Max : "+hist.max());
		System.out.println("Rayleigh : "+hist.rayleighParameter());
		System.out.println("10% : "+hist.percentage(0.10f));
		System.out.println("25% : "+hist.percentage(0.25f));
		System.out.println("50% : "+hist.percentage(0.50f));
		System.out.println("75% : "+hist.percentage(0.75f));
		System.out.println("90% : "+hist.percentage(0.90f));
		
		if (mode==GAUSS) {
			scale = (float)(0.5*hist.mean()*SQPI);
		} 
		distribution = hist;
		distribution.normalizeToMaximum();
		//distribution.normalizeToNonZeroMaximum();
		//distribution.inverseCumulativeHistogram();
		//System.out.println(distribution.printHistogram());
		System.out.println("adjusted scale: "+scale);
		
		return distance;
		
	}
		
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsAverage() {
		
		if (debug) System.out.println("-- weight initialization --");
		
		nlb=0;
		for (int n=0;n<slicecount.length;n++) nlb += slicecount[n];
		if (debug) System.out.println("total number of labels: "+nlb);
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		if (debug) System.out.println("first pass");

		if (slicedir.equals("X")) {
			int offset = 1+slicecount[0];
			for (int x=1;x<nx-1;x++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[x]>0  && slicecount[x-1]+slicecount[x+1]>0) {
					float[] vol0 = new float[slicecount[x]];
					float[] voln = new float[slicecount[x-1]+slicecount[x+1]];
					float[][] assoc0n = new float[slicecount[x]][slicecount[x-1]+slicecount[x+1]];
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x-1][y][z]!=0) {
							lm = (int)(image[x-1][y][z]-offset+slicecount[x-1]);
							voln[lm]++;
						}
						if (image[x+1][y][z]!=0) {
							lp = (int)(image[x+1][y][z]-offset-slicecount[x]);
							voln[lp+slicecount[x-1]]++;
						}
						if (image[x][y][z]!=0 && image[x-1][y][z]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x+1][y][z]!=0) assoc0n[l0][lp+slicecount[x-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[x-1]+slicecount[x+1]];
					for (int l0=0;l0<slicecount[x];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[x-1]+slicecount[x+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[x-1])
								trp[nb] = new Triple(ln+offset-slicecount[x-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							else
								trp[nb] = new Triple(ln-slicecount[x-1]+offset+slicecount[x], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, vol0[l0]);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
				offset += slicecount[x];
			}
		} else
		if (slicedir.equals("Y")) {
			int offset = 1+slicecount[0];
			for (int y=1;y<ny-1;y++) {
				//System.out.println("slice "+y+" ("+slicecount[y-1]+", "+slicecount[y]+", "+slicecount[y+1]+")");
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[y]>0  && slicecount[y-1]+slicecount[y+1]>0) {
					float[] vol0 = new float[slicecount[y]];
					float[] voln = new float[slicecount[y-1]+slicecount[y+1]];
					float[][] assoc0n = new float[slicecount[y]][slicecount[y-1]+slicecount[y+1]];
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y-1][z]!=0) {
							lm = (int)(image[x][y-1][z]-offset+slicecount[y-1]);
							voln[lm]++;
						}
						if (image[x][y+1][z]!=0) {
							lp = (int)(image[x][y+1][z]-offset-slicecount[y]);
							voln[lp+slicecount[y-1]]++;
						}
						if (image[x][y][z]!=0 && image[x][y-1][z]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x][y+1][z]!=0) assoc0n[l0][lp+slicecount[y-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[y-1]+slicecount[y+1]];
					for (int l0=0;l0<slicecount[y];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[y-1]+slicecount[y+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[y-1])
								trp[nb] = new Triple(ln+offset-slicecount[y-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							else
								trp[nb] = new Triple(ln-slicecount[y-1]+offset+slicecount[y], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, vol0[l0]);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
				offset += slicecount[y];
			}
		} else
		if (slicedir.equals("Z")) {
			int offset = 1+slicecount[0];
			for (int z=1;z<nz-1;z++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[z]>0  && slicecount[z-1]+slicecount[z+1]>0) {
					float[] vol0 = new float[slicecount[z]];
					float[] voln = new float[slicecount[z-1]+slicecount[z+1]];
					float[][] assoc0n = new float[slicecount[z]][slicecount[z-1]+slicecount[z+1]];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y][z-1]!=0) {
							lm = (int)(image[x][y][z-1]-offset+slicecount[z-1]);
							voln[lm]++;
						}
						if (image[x][y][z+1]!=0) {
							lp = (int)(image[x][y][z+1]-offset-slicecount[z]);
							voln[lp+slicecount[z-1]]++;
						}
						if (image[x][y][z]!=0 && image[x][y][z-1]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x][y][z+1]!=0) assoc0n[l0][lp+slicecount[z-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[z-1]+slicecount[z+1]];
					for (int l0=0;l0<slicecount[z];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[z-1]+slicecount[z+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[z-1])
								trp[nb] = new Triple(ln+offset-slicecount[z-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							else
								trp[nb] = new Triple(ln-slicecount[z-1]+offset+slicecount[z], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), assoc0n[l0][ln]);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, vol0[l0]);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
				offset += slicecount[z];
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		for (int l=1;l<=nlb;l++) {
			//System.out.println("label "+l);
			//if (debug) System.out.print(".");
			Triple[] node = assoc[l];
			float di = node[0].delta;
			for (int n=1;n<node.length;n++) {
				//float dj = degree.get(node.get(n).id);
				float wij = node[n].weight;
				float wii = node[0].weight;
				float wjj = assoc[node[n].id][0].weight;
				float dj = assoc[node[n].id][0].delta;
				float sij = node[n].size;
				float si = node[0].size;
				float sj = assoc[node[n].id][0].size;
				//float dval = 2.0f*node.get(n).weight/(di + dj);
				// use the formula with self-weights?
				// D = (wii+wjj+2node)/(di+dj) -wii/di -wjj/dj
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj) - wii/di - wjj/dj;
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(di + dj);
				//node.get(n).delta = (wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				//node.get(n).delta = 2.0f*wij/sij*(di + dj)/(si + sj);
				
				// best so far..
				node[n].delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				
				// use the same score as the stopping criterion? slows down the process
				//node.get(n).delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij)
				//					- (1.0f- wij/sij)*wii/si*wjj/sj;
				
				//other[0] += wij/sij;
			}
			cost[0] += node[0].weight/node[0].size/node[0].delta;
			//self[0] += node.get(0).weight/node.get(0).size;
		}
	}
		
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void initAllEdgeWeightsMinimax() {
		
		if (debug) System.out.println("-- weight initialization --");
		
		nlb=0;
		for (int n=0;n<slicecount.length;n++) nlb += slicecount[n];
		
		assoc = new Triple[2*nlb][];
		assoc[0] = new Triple[1];
		assoc[0][0] = new Triple(0);
		
		latest = new int[2*nlb];
		
		if (debug) System.out.println("first pass");

		if (slicedir.equals("X")) {
			int offset = 1+slicecount[0];
			for (int x=1;x<nx-1;x++) {
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[x]>0  && slicecount[x-1]+slicecount[x+1]>0) {
					float[] vol0 = new float[slicecount[x]];
					float[] voln = new float[slicecount[x-1]+slicecount[x+1]];
					float[][] assoc0n = new float[slicecount[x]][slicecount[x-1]+slicecount[x+1]];
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x-1][y][z]!=0) {
							lm = (int)(image[x-1][y][z]-offset+slicecount[x-1]);
							voln[lm]++;
						}
						if (image[x+1][y][z]!=0) {
							lp = (int)(image[x+1][y][z]-offset-slicecount[x]);
							voln[lp+slicecount[x-1]]++;
						}
						if (image[x][y][z]!=0 && image[x-1][y][z]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x+1][y][z]!=0) assoc0n[l0][lp+slicecount[x-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[x-1]+slicecount[x+1]];
					for (int l0=0;l0<slicecount[x];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[x-1]+slicecount[x+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[x-1])
								trp[nb] = new Triple(ln+offset-slicecount[x-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							else
								trp[nb] = new Triple(ln-slicecount[x-1]+offset+slicecount[x], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, 1.0f);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
			}
		} else
		if (slicedir.equals("Y")) {
			int offset = 1+slicecount[0];
			for (int y=1;y<ny-1;y++) {
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[y]>0  && slicecount[y-1]+slicecount[y+1]>0) {
					float[] vol0 = new float[slicecount[y]];
					float[] voln = new float[slicecount[y-1]+slicecount[y+1]];
					float[][] assoc0n = new float[slicecount[y]][slicecount[y-1]+slicecount[y+1]];
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y-1][z]!=0) {
							lm = (int)(image[x][y-1][z]-offset+slicecount[y-1]);
							voln[lm]++;
						}
						if (image[x][y+1][z]!=0) {
							lp = (int)(image[x][y+1][z]-offset-slicecount[y]);
							voln[lp+slicecount[y-1]]++;
						}
						if (image[x][y][z]!=0 && image[x][y-1][z]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x][y+1][z]!=0) assoc0n[l0][lp+slicecount[y-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[y-1]+slicecount[y+1]];
					for (int l0=0;l0<slicecount[y];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[y-1]+slicecount[y+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[y-1])
								trp[nb] = new Triple(ln+offset-slicecount[y-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							else
								trp[nb] = new Triple(ln-slicecount[y-1]+offset+slicecount[y], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, 1.0f);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
			}
		} else
		if (slicedir.equals("Z")) {
			int offset = 1+slicecount[0];
			for (int z=1;z<nz-1;z++) {
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[z]>0  && slicecount[z-1]+slicecount[z+1]>0) {
					float[] vol0 = new float[slicecount[z]];
					float[] voln = new float[slicecount[z-1]+slicecount[z+1]];
					float[][] assoc0n = new float[slicecount[z]][slicecount[z-1]+slicecount[z+1]];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y][z-1]!=0) {
							lm = (int)(image[x][y][z-1]-offset+slicecount[z-1]);
							voln[lm]++;
						}
						if (image[x][y][z+1]!=0) {
							lp = (int)(image[x][y][z+1]-offset-slicecount[z]);
							voln[lp+slicecount[z-1]]++;
						}
						if (image[x][y][z]!=0 && image[x][y][z-1]!=0) assoc0n[l0][lm]++;
						if (image[x][y][z]!=0 && image[x][y][z+1]!=0) assoc0n[l0][lp+slicecount[z-1]]++;
					}
					// build all the weights based on overlap
					Triple[] trp = new Triple[slicecount[z-1]+slicecount[z+1]];
					for (int l0=0;l0<slicecount[z];l0++) {
						// build the weight array
						int nb=0;
						for (int ln=0;ln<slicecount[z-1]+slicecount[z+1];ln++) if (assoc0n[l0][ln]>0) {
							// note: here we can remove unnecessary nodes if desired
							if (ln<slicecount[z-1])
								trp[nb] = new Triple(ln+offset-slicecount[z-1], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							else
								trp[nb] = new Triple(ln-slicecount[z-1]+offset+slicecount[z], associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]), 1.0f);
							
							nb++;								
						}
			
						// link in all conserved directions
						Triple[] node = new Triple[nb+1];
						node[0] = new Triple(l0+offset, 1.0f, 1.0f);
						for (int n=0;n<nb;n++) {
							node[n+1] = trp[n];
						}
						// build the degree array
						float deg = node[0].weight/node[0].size;
						for (int n=0;n<nb;n++) {
							deg += trp[n].weight/trp[n].size/nb;
						}
						//degree.set(labeling[xyz], new Float(deg));
						node[0].delta = deg;
			
						// store the values
						assoc[l0+offset] = node;
			
						// store the latest active index for everything
						latest[l0+offset] = l0+offset;
					}
				}
			}
		}
		
		if (debug) System.out.println("second pass");

		// second pass for delta, cost
		cost = new float[nlb+2];		// store all the cost values
		cost[0] = 0.0f;
		self = new float[nlb+2];		// store all the min self association values
		other = new float[nlb+2];		// store all the max joint association values
		self[0] = 0.0f;
		other[0] = 0.0f;
		for (int l=1;l<=nlb;l++) {
			Triple[] node = assoc[l];
			float di = node[0].delta;
			for (int n=1;n<node.length;n++) {
				//float dj = degree.get(node.get(n).id);
				float wij = node[n].weight;
				float wii = node[0].weight;
				float wjj = assoc[node[n].id][0].weight;
				float dj = assoc[node[n].id][0].delta;
				float sij = node[n].size;
				float si = node[0].size;
				float sj = assoc[node[n].id][0].size;
				
				// best so far..
				//node[n].delta = wij/sij*(wii+wjj+2.0f*wij)/(si + sj + 2.0f*sij);
				//node[n].delta = wij*Numerics.min(wii,wjj,wij);
				//node[n].delta = wij*Numerics.min(wii,wjj,wij) - (1.0f-wij)*(float)Math.sqrt(wii*wjj);
				node[n].delta = Numerics.square(wij*Numerics.min(wii,wjj,wij)) - Numerics.square(1.0f-wij)*wii*wjj;
			}
			cost[0] += node[0].weight/node[0].size/node[0].delta;
		}
	}
	
	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringAverage(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : average (alt)");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// self-weight & degree
					// d(uv) = d(u) + d(v)
					// w(uv,uv) = w(u,u) + w(v,v) + 2*w(uv.uv)
					aNode[0] = new Triple(id, bNode[0].weight + pNode[0].weight + 2.0f*bNode[nbest].weight, 
												bNode[0].delta + pNode[0].delta,  
												bNode[0].size + pNode[0].size + 2.0f*bNode[nbest].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += bNode[n].weight;
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//lbn = latest[lbn];
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								ngbWeight[lbn] += pNode[n].weight;
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						// use the same criterion than for stopping? makes sense, but slows down the process
						// geometric progression:
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/(nlb*nlb)) );
						//float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/((connect+1.0f))) );
						//float wsize = 1.0f;
						aNode[n].delta = wsize*( Numerics.square(aNode[n].weight/aNode[n].size)
													 *Numerics.square( ( aNode[0].weight + wngb[0].weight + 2.0f*aNode[n].weight )
																		/( aNode[0].size + wngb[0].size + 2.0f*aNode[n].size ) )
													 -Numerics.square(1.0f-aNode[n].weight/aNode[n].size)
													  *aNode[0].weight/aNode[0].size
													  *wngb[0].weight/wngb[0].size  );
						
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = 	Numerics.square(bNode[nbest].weight/bNode[nbest].size * aNode[0].weight/aNode[0].size);
					other[iter] = Numerics.square(1.0f-bNode[nbest].weight/bNode[nbest].size) * bNode[0].weight/bNode[0].size
																							* pNode[0].weight/pNode[0].size;
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
					// not a correct stopping criterion								
					if (verbose) if (bcost<0 && first) {
					   System.out.println(iter+" / "+nclusters+": c= "+bcost
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
					}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);

						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}

					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					assoc[lbest][0] = new Triple(id);
					assoc[lpair][0] = new Triple(id);
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
					
					// reduce the number of clusters
					nclusters--;
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
	
	// perform the hierarchical clustering until we reach k0 clusters
	public final void hierarchicalClusteringMinimax(int k0, boolean firststop, int maxlength) {
		
		if (debug) System.out.println("hierarchical clustering : minimax");
		
		// 1. Build a binary tree for the delta values
		bintree = new BinaryHeapPair(nlb+1, Numerics.ceil(0.1f*nlb), BinaryHeapPair.MAXTREE);
		active = new BitSet(nlb+1);
		
		if (debug) System.out.println("initialization");

		for (int lb=1;lb<=nlb;lb++) {
			//if (debug) System.out.print(".");
			Triple[] node = assoc[lb];
			// only store the largest delta (the others will never be selected because it gets relabeled)
			if (node.length>1) {
				int best=1;
				//float bestscore = node.get(1).weight/node.get(1).delta/(node.get(0).delta+assoc.get(node.get(1).id).get(0).delta);
				for (int b=2;b<node.length;b++) {
					if (node[b].delta>node[best].delta) best = b;
				}
				//if (debug) System.out.print(""+node.get(best).delta+","+l+":"+best);
				bintree.addValue(node[best].delta, node[0].id, best);
			} else {
				if (debug) System.out.print("!");	
			}
			// can be active even if there's no labels it is linking to (not symmetric)
			// not the case anymore, no?
			//if (debug) System.out.print(")");
			active.set(lb, true);
			latest[lb] = lb;
		}
		
		if (debug) System.out.println("\n recursive tree building");
		
		// using a large scale storage for neighbors : is it faster or slower ?
		BitSet ngbList = new BitSet(2*nlb);
		float[] ngbWeight = new float[2*nlb];
		float[] ngbSize = new float[2*nlb];
		// init for max() weights
		//for (int n=0;n<2*nlb;n++) ngbWeight[n] = -INF;
		// init for min() weights
		for (int n=0;n<2*nlb;n++) ngbWeight[n] = INF;
		
		// 2. Depile the tree and build recursively the new clusters
		int nclusters = nlb;
		int iter = 0;
		key = nlb;
		boolean first = true;
		IntArray tracks = new IntArray(1000);
		boolean stop = false;
		
		int nstep = 0;
		long newclustertime = 0;
		long depthsearchtime = 0;
		long widthsearchtime = 0;
		long looptime = System.currentTimeMillis();
		
		while (bintree.isNotEmpty() && nclusters>k0 && !stop) {
			//if (debug) System.out.print(".");
			
			// retrive the best delta
			float 	bcost = bintree.getFirst(); 
			int  	lbest = bintree.getFirstId1();
			int 	nbest = bintree.getFirstId2();
			bintree.removeFirst();

			// retrieve corresponding values (if they still exist)
			if (active.get(lbest)) {
				Triple[] bNode = assoc[lbest];
				int	lpair = bNode[nbest].id;
				
				//if (debug) System.out.print("|"+lbest+"-"+lpair);
			
				// update the link label? no, because the weights are now different
				if (active.get(lpair)) {
					Triple[] pNode = assoc[lpair];

					// only count iterations when changing the labels
					iter++;

					// create new cluster
					int id = key+1;
					key++;
					
					// new values
					//ArrayList<Triple> aNode = new ArrayList<Triple>(Numerics.max(bNode.size(),pNode.size())-1);
					Triple[] aNode = new Triple[bNode.length+pNode.length-3];
					// self-weight & degree
					// d(uv) = d(u) + d(v)
					// w(uv,uv) = w(u,u) + w(v,v) + 2*w(uv.uv)
					aNode[0] = new Triple(id, Numerics.min(bNode[0].weight, pNode[0].weight, bNode[nbest].weight), 
												Numerics.min(bNode[0].delta, pNode[0].delta),  
												bNode[0].size + pNode[0].size );
					
					// new mixing weights
					// w(uv,x) = w(u,x) + w(v,x)
					// must check if the links still exist, not duplicate
					
					newclustertime += bNode.length+pNode.length;
					
					// using latest[] allows to preserve the tree structure inside assoc, but skip steps when attributing the labels
					for (int n=1;n<bNode.length;n++) {
						//int lbn = bNode[n].id;
						int lbn = latest[bNode[n].id];
						if (lbn!=lpair) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[bNode[n].id] = lbn;
								//bNode[n].id = lbn;
							}
							if (lbn!=lpair) {
								// add to the label
								ngbList.set(lbn, true);
								// use the min or max of weights along the path?
								// max corresp. to graph-based, not noise prone
								//ngbWeight[lbn] = Numerics.max(bNode[n].weight, ngbWeight[lbn]);
								// min corresp. to most conservative
								ngbWeight[lbn] = Numerics.min(bNode[n].weight, ngbWeight[lbn]);
								ngbSize[lbn] += bNode[n].size;
							}							
							// no need to check on all already created values!!
						}
					}
					for (int n=1;n<pNode.length;n++) {
						//int lbn = pNode[n].id;
						int lbn = latest[pNode[n].id];
						if (lbn!=lbest) {
							if (!active.get(lbn)) {
								// make sure it's the most up-to-date version of the label
								tracks.reset();
								while (!active.get(lbn)) {
									tracks.add(lbn);
									lbn = latest[assoc[lbn][0].id];
									//lbn = assoc[lbn][0].id;
									depthsearchtime++;
								}
								// update the id if not up to date
								for (int lb=0;lb<tracks.last;lb++)
									latest[tracks.val[lb]] = lbn;
									//assoc[tracks.val[lb]][0].id = lbn;
								latest[pNode[n].id] = lbn;
								//pNode[n].id = lbn;
								
							}
							if (lbn!=lbest) {
								// add to the label
								ngbList.set(lbn, true);
								//ngbWeight[lbn] = Numerics.max(pNode[n].weight, ngbWeight[lbn]);
								ngbWeight[lbn] = Numerics.min(pNode[n].weight, ngbWeight[lbn]);
								ngbSize[lbn] += pNode[n].size;
							}
							// no need to check on all already created values!!
						}
					}
					// build the list
					int l=1;
					for (int lbn = ngbList.nextSetBit(0); lbn >= 0; lbn = ngbList.nextSetBit(lbn+1)) {
						// create a new one
						aNode[l] = new Triple(lbn, ngbWeight[lbn], ngbSize[lbn]);
						l++;
						// reset the values
						ngbWeight[lbn] = 0.0f;
						ngbSize[lbn] = 0.0f;
					}
					ngbList.clear();

					// make sure we don't have extra empty values
					//aNode.trimToSize();
					if (l<aNode.length) {
						Triple[] tmp = new Triple[l];
						for (int n=0;n<l;n++) tmp[n] = aNode[n];
						aNode = tmp;
					}
					// new deltas : no need of the u,v values anymore
					// D(uv,x) = ( w(uv,uv) + w(x,x) + 2w(uv,x) )/( d(uv)+d(x) ) - w(uv,uv)/d(uv) - w(x,x)/d(x)
					for (int n=1; n<aNode.length; n++) {
						Triple[] wngb = assoc[aNode[n].id];
						// best method so far
						//aNode[n].delta = aNode[n].weight*Numerics.min(aNode[0].weight, wngb[0].weight, aNode[n].weight);
												
						// use the same criterion than for stopping? makes sense, but slows down the process
						float wsize = 1.0f/(1.0f + Numerics.square( (aNode[0].size+wngb[0].size)/nlb) );
						//float wsize = 1.0f;
						aNode[n].delta = wsize*( Numerics.square(aNode[n].weight *Numerics.min(aNode[0].weight, wngb[0].weight, aNode[n].weight ))
											     -Numerics.square(1.0f-aNode[n].weight)*aNode[0].weight*wngb[0].weight );
					}
					
					/// probably not needed
					// recompute the degree? (for averaged links)
					self[iter] = Numerics.square(bNode[nbest].weight * aNode[0].weight);
					other[iter] = Numerics.square(1.0f-bNode[nbest].weight) * bNode[0].weight * pNode[0].weight;
					
					cost[iter] = bcost;

					if (verbose) if (iter%(nlb/100)==0) {
						long newtime = System.currentTimeMillis();
						nstep++;
						System.out.println("n="+nstep+", t="+(newtime-looptime)+", "+iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length
																	+"| n "+newclustertime+", d "+depthsearchtime+", w "+widthsearchtime+")");
						looptime = newtime;
						newclustertime = 0;
						depthsearchtime = 0;
						widthsearchtime = 0;
					}
																	
				   if (verbose) if (bcost<0 && first) {
				   	   // that's not the originally computed score, is it??
				   	   System.out.println(iter+" / "+nclusters+": c= "+bcost+" | "+(self[iter]-other[iter])
																	+", s= "+self[iter]
																	+", o= "+other[iter]
																	+" ("+bintree.getCurrentSize()+"|"
																	+bNode.length+", "+pNode.length+")");	
					   first=false;
					   if (firststop) stop = true;
				   	}								
					// add the new values to list, binary tree
					assoc[id] = aNode;
					active.set(id, true);
					
					if (aNode.length>1) {
						int best=1;
						for (int b=2;b<aNode.length;b++) {
							if (aNode[b].delta>aNode[best].delta) best = b;
						}
						bintree.addValue(aNode[best].delta, id, best);

						if (debug) if (aNode[best].delta>1) {
							System.out.println(nclusters+": c= "+cost[iter]+" ("
																	+bNode[0].delta+", "
																	+bNode[0].weight+", "
																	+bNode[0].size+" | "
																	+pNode[0].delta+", "
																	+pNode[0].weight+", "
																	+pNode[0].size+" | "
																	+aNode[0].delta+", "
																	+aNode[0].weight+", "
																	+aNode[0].size+")");					   
							
							for (int n=1; n<aNode.length; n++) {
								System.out.print("<"+aNode[n].delta+", "+aNode[n].weight+", "+aNode[n].size+">");
							}
							System.out.print("\n");
						}

					}
					
					// replace the older values with info on what is the new label
					assoc[lbest] = null;
					assoc[lpair] = null;
					assoc[lbest] = new Triple[1];
					assoc[lpair] = new Triple[1];
					assoc[lbest][0] = new Triple(id);
					assoc[lpair][0] = new Triple(id);
					
					// de-activate the labels
					active.set(lbest, false);
					active.set(lpair, false);
					
					// update the label mapping
					latest[lbest] = id;
					latest[lpair] = id;
					latest[id] = id;
										
					// reduce the number of clusters
					nclusters--;
				}
			}
		}
		// done!
		if (debug) System.out.println("completed ("+nclusters+")");
		
	}
			
	public final int firstAssociationLoss(int nt) {
		if (debug) System.out.println("find largest number of clusters with high association");

		int zerolb = -1;
		for (int t=1;t<nlb;t++) {
			//if (self[t]<other[t]) {
			if (cost[t]<0) {
				zerolb = t;
				break;
			}
		}
		if (zerolb==-1) zerolb = nlb-nt-1;
		if (debug) System.out.println("-> "+(nlb-zerolb));
		return (nlb-zerolb);
	}
	
	private final void 	updateLatestLabels(int maxlb) {	
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
	}
	
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final void mergeOverlapForward() {
		
		if (debug) System.out.println("-- merge slices forward--");
		
		int nchanged = 0;
		
		nlb=1;
		for (int n=0;n<slicecount.length;n++) {
			nlb += slicecount[n];
		}
		if (debug) System.out.println("total number of labels: "+nlb);
		
		latest = new int[nlb];
		for (int n=0;n<nlb;n++) latest[n] = n;
		
		if (slicedir.equals("X")) {
			int offset = 1+slicecount[0];
			for (int x=1;x<nx-1;x++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[x]>0  && slicecount[x-1]>0) {
					float[] vol0 = new float[slicecount[x]];
					float[] voln = new float[slicecount[x-1]];
					float[][] assoc0n = new float[slicecount[x]][slicecount[x-1]];
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x-1][y][z]!=0) {
							lm = (int)(image[x-1][y][z]-offset+slicecount[x-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x-1][y][z]!=0) assoc0n[l0][lm]++;
					}
					// build all the distances based on overlap
					for (int l0=0;l0<slicecount[x];l0++) {
						// build the distance array
						for (int ln=0;ln<slicecount[x-1];ln++) if (assoc0n[l0][ln]>0) {
							// decide there and then to change the previous slice labels
							float cost = associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]);		
							if (cost>basis) {
								latest[ln+offset-slicecount[x-1]]= l0+offset;
								nchanged++;
							}
						}
					}
				}
				offset += slicecount[x];
			}
		} else
		if (slicedir.equals("Y")) {
			int offset = 1+slicecount[0];
			for (int y=1;y<ny-1;y++) {
				//System.out.println("slice "+y+" ("+slicecount[y-1]+", "+slicecount[y]+", "+slicecount[y+1]+")");
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[y]>0  && slicecount[y-1]>0) {
					float[] vol0 = new float[slicecount[y]];
					float[] voln = new float[slicecount[y-1]];
					float[][] assoc0n = new float[slicecount[y]][slicecount[y-1]];
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y-1][z]!=0) {
							lm = (int)(image[x][y-1][z]-offset+slicecount[y-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x][y-1][z]!=0) assoc0n[l0][lm]++;
					}
					// build all the weights based on overlap
					for (int l0=0;l0<slicecount[y];l0++) {
						// build the weight array
						for (int ln=0;ln<slicecount[y-1];ln++) if (assoc0n[l0][ln]>0) {
							// decide there and then to change the previous slice labels
							float cost = associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]);		
							if (cost>basis) {
								latest[ln+offset-slicecount[y-1]] = l0+offset;
								nchanged++;
							}
						}
					}
				}
				offset += slicecount[y];
			}
		} else
		if (slicedir.equals("Z")) {
			int offset = 1+slicecount[0];
			for (int z=1;z<nz-1;z++) {
				//if (debug) System.out.print("=");
				// links to preceding and following slices: count the overlaps globally
				if (slicecount[z]>0  && slicecount[z-1]>0) {
					float[] vol0 = new float[slicecount[z]];
					float[] voln = new float[slicecount[z-1]];
					float[][] assoc0n = new float[slicecount[z]][slicecount[z-1]];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l0=-1,lm=-1,lp=-1;
						if (image[x][y][z]!=0) {
							l0 = (int)(image[x][y][z]-offset);
							vol0[l0]++;
						}
						if (image[x][y][z-1]!=0) {
							lm = (int)(image[x][y][z-1]-offset+slicecount[z-1]);
							voln[lm]++;
						}
						if (image[x][y][z]!=0 && image[x][y][z-1]!=0) assoc0n[l0][lm]++;
					}
					// build all the weights based on overlap
					for (int l0=0;l0<slicecount[z];l0++) {
						// build the weight array
						for (int ln=0;ln<slicecount[z-1];ln++) if (assoc0n[l0][ln]>0) {
							// decide there and then to change the previous slice labels
							float cost = associationWeight(assoc0n[l0][ln], vol0[l0], voln[ln]);		
							if (cost>basis) {
								latest[ln+offset-slicecount[z-1]] = l0+offset;
								nchanged++;
							}
						}
					}
				}
				offset += slicecount[z];
			}
		}
		if (debug) System.out.println("number of changed labels: "+nchanged);
		
		return;		
	}
		
	// initial lists: create all the links, thus the images are not needed anymore ?
	public final float[][][] relabelFromSegmentation(int[][][] segmentation) {
		
		if (debug) System.out.println("-- relabel segmentation--");
		
		int nchanged = 0;
		
		nlb=1;
		for (int n=0;n<slicecount.length;n++) {
			nlb += slicecount[n];
		}
		if (debug) System.out.println("total number of labels: "+nlb);
		
		latest = new int[nlb];
		for (int n=0;n<nlb;n++) latest[n] = n;
		
		int[] segment = ObjectLabeling.listOrderedLabels(segmentation, nx, ny, nz);
		// make more labels than needed (faster)
		int nseg = segment[segment.length-1]+1;
		if (debug) System.out.println("max label of structures: "+nseg);
		
		float[][][] weight = new float[nx][ny][nz];
		
		if (slicedir.equals("X")) {
			int offset = 0;
			for (int x=0;x<nx;x++) {
				//if (debug) System.out.print("=");
				if (slicecount[x]>0) {
					float[] vol = new float[slicecount[x]];
					float[][] assoc = new float[slicecount[x]][nseg];
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							vol[l]++;
							assoc[l][(int)segmentation[x][y][z]]++;
						}
					}
					// find the best, assign overlap
					float[] bestproba = new float[slicecount[x]];
					float[] bestlabel = new float[slicecount[x]];
					for (int l=0;l<slicecount[x];l++) {
						bestproba[l] = 0.0f;
						bestlabel[l] = -1;
						// build the distance array
						for (int n=0;n<nseg;n++) if (assoc[l][n]>0) {
							
							float proba = assoc[l][n]/vol[l];
							if (proba > bestproba[l]) {
								bestlabel[l] = n;
								bestproba[l] = proba;
							}
						}
					}
					// write it back onto the image
					for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							image[x][y][z] = bestlabel[l];
							weight[x][y][z] = bestproba[l];
						}
					}
					
				}
				offset += slicecount[x];
			}
		} else
		if (slicedir.equals("Y")) {
			int offset = 1+slicecount[0];
			for (int y=1;y<ny-1;y++) {
				//if (debug) System.out.print("=");
				if (slicecount[y]>0) {
					float[] vol = new float[slicecount[y]];
					float[][] assoc = new float[slicecount[y]][nseg];
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							vol[l]++;
							assoc[l][(int)segmentation[x][y][z]]++;
						}
					}
					// find the best, assign overlap
					float[] bestproba = new float[slicecount[y]];
					float[] bestlabel = new float[slicecount[y]];
					for (int l=0;l<slicecount[y];l++) {
						bestproba[l] = 0.0f;
						bestlabel[l] = -1;
						// build the distance array
						for (int n=0;n<nseg;n++) if (assoc[l][n]>0) {
							
							float proba = assoc[l][n]/vol[l];
							if (proba > bestproba[l]) {
								bestlabel[l] = n;
								bestproba[l] = proba;
							}
						}
					}
					// write it back onto the image
					for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							image[x][y][z] = bestlabel[l];
							weight[x][y][z] = bestproba[l];
						}
					}
				}
				offset += slicecount[y];
			}
		} else
		if (slicedir.equals("Z")) {
			int offset = 1+slicecount[0];
			for (int z=1;z<nz-1;z++) {
				//if (debug) System.out.print("=");
				if (slicecount[z]>0) {
					float[] vol = new float[slicecount[z]];
					float[][] assoc = new float[slicecount[z]][nseg];
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							vol[l]++;
							assoc[l][(int)segmentation[x][y][z]]++;
						}
					}
					float[] bestproba = new float[slicecount[z]];
					float[] bestlabel = new float[slicecount[z]];
					for (int l=0;l<slicecount[z];l++) {
						bestproba[l] = 0.0f;
						bestlabel[l] = -1;
						// build the distance array
						for (int n=0;n<nseg;n++) if (assoc[l][n]>0) {
							
							float proba = assoc[l][n]/vol[l];
							if (proba > bestproba[l]) {
								bestlabel[l] = n;
								bestproba[l] = proba;
							}
						}
					}
					// write it back onto the image
					for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) {
						int l=-1;
						if (image[x][y][z]!=0) {
							l = (int)(image[x][y][z]-offset);
							image[x][y][z] = bestlabel[l];
							weight[x][y][z] = bestproba[l];
						}
					}
				}
				offset += slicecount[z];
			}
		}
		if (debug) System.out.println("number of changed labels: "+nchanged);
		
		return weight;		
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringSequential(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = latest[(int)image[x][y][z]];
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringConsistency(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = assoc[latest[(int)image[x][y][z]]][0].weight/assoc[latest[(int)image[x][y][z]]][0].size;
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportClusteringDifferences(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			double sum = 0.0;
			double den = 0.0;
			Triple[] node = assoc[latest[(int)image[x][y][z]]];
			for (int n=1;n<node.length;n++) {
				sum += node[n].weight;
				den += node[n].size;
			}
			tmp[x][y][z] = (float)(sum/den);
		}
		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] filteredClusters(int maxlb, float minratio, float maxratio) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// create volume measurements
		int[] volume = new int[2*nlb];
		int totalvol = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			volume[latest[(int)image[x][y][z]]]++;
			totalvol++;
		}
		
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (volume[latest[(int)image[x][y][z]]]>minratio*totalvol && volume[latest[(int)image[x][y][z]]]<maxratio*totalvol) {
				tmp[x][y][z] = latest[(int)image[x][y][z]]+1;
			} else {
				tmp[x][y][z] = 1;
			}
		}
		
		// second pass to get smaller label values?
		float[] lblist = ObjectLabeling.listLabels(tmp, nx, ny, nz);
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			boolean found = false;
			for (int n=0;n<lblist.length && !found;n++) {
				if (lblist[n]==tmp[x][y][z]) {
					tmp[x][y][z] = n;
					found = true;
				}
			}
		}

		return tmp;
	}
	
	// given a clustering result, update the labeling
	public final float[][][] exportClusterBoundaryScore() {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = boundaryScore[(int)image[x][y][z]];
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportIntensity() {
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = image[x][y][z];
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final float[][][] exportMergedClusters() {
		for (int n=nlb-1; n>0; n--) {
			int lb = latest[n];
			if (lb!=n) {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		float[][][] tmp = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			tmp[x][y][z] = latest[(int)image[x][y][z]];
		}
		return tmp;
	}
		
	// given a clustering result, update the labeling
	public final String clusteringStatisticsSequential(int maxlb) {
		// update the labels sequentially first from top to bottom (linear)
		for (int n=0;n<2*nlb;n++) latest[n] = 0;
		
		for (int n=2*nlb-maxlb-1; n>0; n--) {
			int lb = assoc[n][0].id;
			if (lb==n) latest[n] = lb;
			else {
				// recurse up to the top one
				while (latest[lb]>0 && latest[lb]!=lb && latest[lb]<2*nlb-maxlb) {
					lb = latest[lb];
				}
				latest[n] = lb;
			}
		}
		// create volume measurements
		int[] volume = new int[2*nlb];
		int totalvol = 0;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			volume[latest[(int)image[x][y][z]]]++;
			totalvol++;
		}
		int maxvol = Numerics.max(volume);
		int nclusters = 0;
		for (int n=0;n<2*nlb;n++) if (volume[n]>0) nclusters++; 
		String stats = "cluster number : "+nclusters+", ratio (%) = "+(100.0f*nclusters/(float)totalvol)+
		", \n volume : "+totalvol+", max cluster size (%) = "+(100.0f*maxvol/(float)totalvol);

		return stats;
	}

	private final float associationWeight(float overlap, float vol1, float vol2) {
		if (mode==DIRECT) return 1.0f-dataDistance(overlap, vol1, vol2);
		else if (mode==GAUSS) return gaussDistanceWeight(overlap, vol1, vol2);
		else if (mode==NP) return nonParamDistanceWeight(overlap, vol1, vol2);
		else return 0.0f;
	}

	private final float dataDistance(float overlap, float vol1, float vol2) {
		if (metric==DICE) return 1.0f-2.0f*overlap/(vol1+vol2);
		else if (metric==JACCARD) return 1.0f-overlap/(vol1+vol2-overlap);
		else if (metric==OVERLAP_SIZE) return 1.0f-overlap/(overlap+minsize);
		else if (metric==MIX) return 1.0f-overlap/(vol1+vol2-overlap)*overlap/(overlap+minsize);
		else if (metric==RATIO) return minsize/overlap;
		else return 0.5f;
	}
	
	private final float gaussDistanceWeight(float overlap, float vol1, float vol2) {
		return (float)FastMath.exp( -0.5*Numerics.square(dataDistance(overlap, vol1, vol2))/(scale*scale) );
	}
	private final float nonParamDistanceWeight(float overlap, float vol1, float vol2) {
		return distribution.getHistogramCount(dataDistance(overlap, vol1, vol2));
	}

}
