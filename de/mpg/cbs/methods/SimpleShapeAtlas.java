package de.mpg.cbs.methods;

import java.io.*;
import java.util.*;

import de.mpg.cbs.libraries.*;
import de.mpg.cbs.utilities.*;
import de.mpg.cbs.structures.*;

/**
 *
 *  This class handles basic shape atlas information
 *
 *	@version    Feb 2011
 *	@author     Pierre-Louis Bazin
 *		
 *
 */
 
public class SimpleShapeAtlas {
	
	// structures: basic information	
	private		int					nobj;			// number of structures in the atlas
	private		String[]			objName;				// their objNames
	private		byte[]				objLabel;				// their objLabels
	private		String[]			objType;			// their objType type
	
	private		String				atlasFile;			// the atlas file
	
	// atlas quantities
	private 	int 				nix,niy,niz; 			// image dimensions
	private 	float 				rix,riy,riz; 			// image resolutions
	private		int					orient,orix,oriy,oriz;		// image orientations
	private		float				x0i,y0i,z0i;		// the center of the image
	
	// spatial transformations
	private		float[]				transform;		// the transform parameters to get into image space
	private		float[][]			rotation;		// the associated rotation matrix
	private		float[][]			shapeTransform; // the global transform matrix (XI = A XP) for each shape
	private		float[][]			inverseShapeTransform; // the global transform matrix (XI = A XP) for each shape
	private		int					Nd;				// transform dimension
	private		float				maxscale = 4.0f;
	private		float				scalingFactor = 0.1f; // maximum variation due to scaling
	
	// shape maps
	private		float[][]			shape;				// the shape images
	private		String[]			shapeFile;			// the shape fileobjNames
	private		int					nax,nay,naz;		// the shape dimensions 
	private 	float				rax,ray,raz; 		// the shape resolutions
	private		float				x0a,y0a,z0a;		// the center of the shape image
	
	private		float				shapeScale;			//the slope of the sigmoid prior based on distance functions
	private		int					objLabelSamples;		// number of samples in the shape/distasnce/contact/direction model
	private		int[]				minx,miny,minz;		// the lowest image coordinate with non zero prior for each shape
	private		int[]				maxx,maxy,maxz;		// the highest image coordinate with non zero prior for each shape
	
	// topology template
	private		byte[]				template;		// the topology template image
	private		String				templateFile;	// the topology template file
	private		int					ntx,nty,ntz;	// the topology template dimensions
	private 	float 				rtx,rty,rtz; 	// the topology template resolutions
	private		float				x0t,y0t,z0t;	// the center of the topology template
	
	// intensity models
	private		float[][]			intensity;		// the intensity models, normalised between 0 and 1
	public		static final int	MPRAGE3T = 0;
	public		static final int	MP2RAGE7T = 1;
	public		static final int	T1MAP7T = 2;
	public		static final int	T2SW7T = 3;
	public		static final int	MP2RAGE3T = 4;
	public		static final int	T1MAP3T = 5;
	public		static final int	PVDURA  = 6;
	public		static final int	FILTER  = 7;
	public		static final int	HCPT1W  = 8;
	public		static final int	HCPT2W  = 9;
	public		static final int	QSM7T = 10;
	public		static final int	NORMMPRAGE = 11;
	
	private		static final int	INTENSITY = 12;	// the number of possible intensity models
	
	// new intensity models
	private		float[][][]		intensityMap = null;
	
	// other parameters
	private		float[]			regularizationFactor = null;
	
	// for registration
	private 	float[]			famap;
	private 	float[][]		mems;
	private 	short[][]		lbmems;
	private		int				nbest;
	private		int				dataSource;
	private	static final int	FA=1;
	private	static final int	MEMS=2;
	
	private		ParametricTransform		transformModel;	// the type of transform (from possible ones below)
	private		int						transformMode;	// the type of transform (from possible ones below)
	private static final	int   		NONE = 0;
	private static final	int   		PARAMETRIC = 1;
	private		boolean[]				registeredShape;
		
	// preset computation arrays for speed up
	private		float[]		XP;
	
	// constants
	private static final	float	PI2 = (float)(Math.PI/2.0);
	private static final	float	ISQRT2 =  (float)(1.0/Math.sqrt(2.0f));
	
	// for debug and display
	private static final boolean		debug=true;
	private static final boolean		verbose=false;
	
	
	// numerics
	private static final	float   INF=1e30f;
	private static final	float   ZERO=1e-30f;

	// convenience tags
	private static final	int   X=0;
	private static final	int   Y=1;
	private static final	int   Z=2;
	private static final	int   T=3;

	/**
	 *	constructor: load the atlas information from a file.
	 *	<p>
	 *	The atlas files follow a certain template; 
	 *  the separator between numbers is a tab, not a space.
	 */
	public SimpleShapeAtlas(String fileobjName) {
		
		transformMode = NONE;
		Nd = 0;
		transform = null;
		
		objLabelSamples = 0;
		loadAtlas(fileobjName);
		
		XP = new float[3];
	}
	
	/**
	 *	constructor: create an empty atlas.
	 */
	public SimpleShapeAtlas(int Nc_) {
		
		nobj = Nc_;
		
		transformMode = NONE;
		Nd = 0;
		transform = null;
		
		objLabelSamples = 0;
		
		// allocate everiything
		objName = new String[nobj];
		objLabel = new byte[nobj];
		objType = new String[nobj];
			
		shape = new float[nobj][];
		minx = new int[nobj];
		miny = new int[nobj];
		minz = new int[nobj];
		maxx = new int[nobj];
		maxy = new int[nobj];
		maxz = new int[nobj];
		
		intensity = new float[INTENSITY][nobj];
	}
	
	/** clean-up: destroy membership and centroid arrays */
	public final void finalize() {
		shape = null;
		template = null;
		System.gc();
	}
	
	/** link the variables */
	final public int 		getNumber() { return nobj; }
	final public String[] 	getNames() { return objName; }
	final public String[] 	getObjectType() { return objType; }
	final public byte[]		getLabels() { return objLabel; }
	
	final public int		getObject(String name) {
		for (int n=0;n<nobj;n++) if (objName[n].equals(name)) return n;
		return -1;
	}
	
	final public byte[] 	getTemplate() { return template; }
	final public void 	setTemplate(byte[] tpl) { template = tpl; }
	
	final public float[] 	getShape(int n) { return shape[n]; }
	final public float[][] 	getShapes() { return shape; }
	
	final public float[] 		getTransform() { 		
		return transform; 
	}
	final public void 		setTransform(float[] trans) { transform = trans; }
	
	final public boolean isRegistered(int k) {
		return registeredShape[k];
	}
		
	final public float[]	exportLabels() {
		float[] lb = new float[nobj];
		for (int k=0;k<nobj;k++) lb[k] = (float)objLabel[k];
		return lb;
	}

	final public int[] getTemplateDim() {
		int[] dim = new int[3];
		dim[0] = ntx;
		dim[1] = nty;
		dim[2] = ntz;
		return dim;
	}
	
   final public int[] getShapeDim() {
		int[] dim = new int[3];
		dim[0] = nax;
		dim[1] = nay;
		dim[2] = naz;
		return dim;
	}
	
   final public int getShapeSize() {
		return nax*nay*naz;
	}
	
   final public float[] getShapeRes() {
		float[] res = new float[3];
		res[0] = rax;
		res[1] = ray;
		res[2] = raz;
		return res;
	}
	
    final public int[] getImageDim() {
		int[] dim = new int[3];
		dim[0] = nix;
		dim[1] = niy;
		dim[2] = niz;
		return dim;
	}
	
    final public float[] getImageRes() {
		float[] res = new float[3];
		res[0] = rix;
		res[1] = riy;
		res[2] = riz;
		return res;
	}
	
	final public int[] getImageOrient() {
		int[] ori = new int[4];
		ori[0] = orient;
		ori[1] = orix;
		ori[2] = oriy;
		ori[3] = oriz;
		return ori;
	}
	final public boolean hasTopology() { return template!=null; }
	final public boolean hasShape(int id) { return shape[id]!=null; }
	
	final public float[][] 	getIntensityPriors(String[] modality, int nc) {
		float[][]	prior = new float[nc][nobj];
		for (int n=0;n<nc;n++) {
			if (contrastId(modality[n])>-1) {
				for (int k=0;k<nobj;k++) prior[n][k] = intensity[contrastId(modality[n])][k];
			}	
		}
		return prior;
	}

	final public float[][] 	getMap(int modal) {
		return intensityMap[modal];
	}

	final public float[] 	getRegularizationFactor() {
		return regularizationFactor;
	}

	final public void 	setMap(int modal, int n, int t, float val) {
		intensityMap[modal][n][t] = val;
	}

	final public void addLabelSample() { objLabelSamples++; }
	
    /** 
	 *  set image-related information for segmentation
	 */
	final public void setImageInfo(int nix_, int niy_, int niz_, float rix_, float riy_, float riz_, int orient_, int orix_, int oriy_, int oriz_) {
		nix = nix_; niy = niy_; niz = niz_;
		rix = rix_; riy = riy_; riz = riz_;
		orient = orient_;
		orix = orix_; oriy = oriy_; oriz = oriz_;
		
		x0i = nix/2.0f;
		y0i = niy/2.0f;
		z0i = niz/2.0f;
		
		if (debug) {
			System.out.print("dimensions: "+nix+", "+niy+", "+niz+"\n");
			System.out.print("resolutions: "+rix+", "+riy+", "+riz+"\n");
			System.out.print("orientation: "+orient+" | "+orix+", "+oriy+", "+oriz+"\n");
			System.out.print("Atlas\n");
			System.out.print("dimensions: "+nax+", "+nay+", "+naz+"\n");
			System.out.print("resolutions: "+rax+", "+ray+", "+raz+"\n");
		}
	}
	
   /** 
	 *  set image-related information for segmentation
	 */
	final public void adjustAtlasScale(float[][] image, int nimg) {
		// compute the brain size on the image and the atlas; update atlas scale to match
		double[] imgmax = new double[nimg];
		for (int xyz=0;xyz<nix*niy*niz;xyz++) {
			for (int i=0;i<nimg;i++) {
				if (Numerics.abs(image[i][xyz])>imgmax[i]) imgmax[i] = Numerics.abs(image[i][xyz]);
			}
		}
		double imgvol=0.0;
		for (int xyz=0;xyz<nix*niy*niz;xyz++) {
			boolean isMasked = true;
			for (int i=0;i<nimg;i++) {
				if (Numerics.abs(image[i][xyz])>0.1f*imgmax[i]) isMasked = false;
			}
			if (!isMasked) imgvol++;
		}
		System.out.print("Image-based brain volume: "+(imgvol*rix*riy*riz)+" mm^3\n");
		
		float atlasvol=0.0f;
		for (int xyz=0;xyz<nax*nay*naz;xyz++) {
			boolean isMasked = true;
			for (int n=1;n<nobj;n++) {
				if (shape[n][xyz]>shape[0][xyz]) isMasked = false;
			}
			if (!isMasked) atlasvol++;
		}
		System.out.print("Atlas-based brain volume: "+(atlasvol*rax*ray*raz)+" mm^3\n");
		
		float scaling = (float)Math.cbrt( (imgvol*rix*riy*riz)/(atlasvol*rax*ray*raz) );
		System.out.print("Scaling factor: "+scaling+" mm^3\n");
		
		rax *= scaling;
		ray *= scaling;
		raz *= scaling;
		
		rtx *= scaling;
		rty *= scaling;
		rtz *= scaling;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public byte[] generateTransformedClassification() {
		float dist,max,count,val;
		byte best;
		byte[] img = new byte[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			// compute each class probability : attribute the highest
			max = 0; best = -1;
			for (byte k=0;k<nobj;k++) {
				val = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
				if (val>max) {
					best = k;
					max = val;
				}
			}
			img[x+y*nix+z*nix*niy] = best;
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public int[] generateTransformedClassificationLabeling() {
		float dist,max,count,val;
		int best=0;
		int[] img = new int[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			// compute each class probability : attribute the highest
			max = 0; best = -1;
			for (int k=0;k<nobj;k++) {
				val = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
				if (val>max) {
					best = k;
					max = val;
				}
			}
			if (best>-1) img[x+y*nix+z*nix*niy] = objLabel[best];
			else img[x+y*nix+z*nix*niy] = 0;
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedClassificationFloat() {
		float dist,max,count,val;
		int best=0;
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			// compute each class probability : attribute the highest
			max = 0; best = -1;
			for (int k=0;k<nobj;k++) {
				val = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
				if (val>max) {
					best = k;
					max = val;
				}
			}
			if (best>-1) img[x+y*nix+z*nix*niy] = objLabel[best];
			else img[x+y*nix+z*nix*niy] = 0;
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[][] generateTransformedShapes() {
		float[][] img = new float[nobj][nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			for (int k=0;k<nobj;k++) {	
				img[k][x+y*nix+z*nix*niy] = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
			}
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedShape(int k) {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			img[x+y*nix+z*nix*niy] = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedTemplate() {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			img[x+y*nix+z*nix*niy] = ImageInterpolation.nearestNeighborInterpolation(template,(byte)1,XP[0],XP[1],XP[2],nax,nay,naz);	
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedObject(String type) {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		int xyz;
		
    	boolean[] iswm = new boolean[nobj];
    	for (int n=0;n<nobj;n++) 
    		if (objType[n].startsWith(type)) iswm[n] = true;
    		else iswm[n] = false;

    	for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			xyz = x+y*nix+z*nix*niy;
			img[xyz] = 0.0f;
			for (int k=0;k<nobj;k++) if (iswm[k]) {	
				img[xyz] = Numerics.max(img[xyz],ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz));	
			}
		}
		return img;
	}
	
    /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedMap(float[] map) {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP,x,y,z);
			
			img[x+y*nix+z*nix*niy] = ImageInterpolation.linearInterpolation(map,0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
		}
		return img;
	}
	
   /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedMap(byte[] map) {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP,x,y,z);
			
			img[x+y*nix+z*nix*niy] = ImageInterpolation.nearestNeighborInterpolation(map,(byte)0,XP[0],XP[1],XP[2],nax,nay,naz);	
		}
		return img;
	}
	
  /** 
	 *  generate atlas image from information
	 */
    final public float[] generateTransformedMap(int[] map) {
		float[] img = new float[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP,x,y,z);
			
			img[x+y*nix+z*nix*niy] = ImageInterpolation.nearestNeighborInterpolation(map,0,XP[0],XP[1],XP[2],nax,nay,naz);	
		}
		return img;
	}
	
	/** display the atlas data */
	final public String displayNames() {
		String output = "Structures \n";
		
		for (int k=0;k<nobj;k++) {
			output += objName[k]+" ("+objType[k]+")	"+objLabel[k]+"\n";
		}
		
		return output;	
	}
		/** display the atlas data */
	final public String displayRegisteredShapes() {
		String output = "Registered Shapes (0/1: true/false) \n";
		
		for (int k=0;k<nobj;k++) output += objLabel[k]+"	";
		output += "\n";
		for (int k=0;k<nobj;k++) {
			if (registeredShape[k]) output += "1	";
			else output += "0	";
		}
		output += "\n";
		
		return output;	
	}
		/** display the atlas data */
	final public String displayRegularizationFactor() {
		String output = "Regularization Factor \n";
		
		for (int k=0;k<nobj;k++) output += objLabel[k]+"	";
		output += "\n";
		for (int k=0;k<nobj;k++) {
			output += regularizationFactor[k]+"	";
		}
		output += "\n";
		
		return output;	
	}
	/**
	 *	read template image (the image must be in bytes)
	 */
	private byte[] loadTemplateImage(String fileobjName, int Nx, int Ny, int Nz) {
		// read the raw data
		byte[] buffer = null;
		try {
          File f = new File( fileobjName );
		  //System.out.println("exists ? "+f.exists());
          //System.out.println("can read ? "+f.canRead());
          FileInputStream fis = new FileInputStream( f );
            
		   buffer = new byte[Nx*Ny*Nz];
		   fis.read(buffer);
           fis.close();
		} catch (IOException io) {
           System.out.println("i/o pb: "+io.getMessage());
		}
		return buffer;
	}
	
	/**
	 *	read shape image (the image must be in float, little endian)
	 */
	private final float[] loadShapeImage(String fileobjName, int Nx, int Ny, int Nz) {
		// read the raw data
		byte[] buffer = null;
		try {
           File f = new File( fileobjName );
           FileInputStream fis = new FileInputStream( f );
            
		   buffer = new byte[4*Nx*Ny*Nz];
		   fis.read(buffer);
           fis.close();
		} catch (IOException io) {
           System.out.println("i/o pb: "+io.getMessage());
		}
		// convert to the image format
		float [] img  = new float[Nx*Ny*Nz];

		for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) {
			int b1 = buffer[4*(xyz)+0] & 0xff;
			int b2 = buffer[4*(xyz)+1] & 0xff;
			int b3 = buffer[4*(xyz)+2] & 0xff;
			int b4 = buffer[4*(xyz)+3] & 0xff;
			// big endian
			//int tmpInt = ((b1 << 24) | (b2 << 16) | (b3 << 8) | b4);
			// little endian
			int tmpInt = ((b4 << 24) | (b3 << 16) | (b2 << 8) | b1);

			img[xyz] = Float.intBitsToFloat(tmpInt);
		}
        buffer = null;
		
		return img;
	}
	
	/**
	 *	read shape image (the image must be in float, little endian)
	 */
	private final float[][] loadDirectionImage(String fileobjName, int Nx, int Ny, int Nz) {
		// read the raw data
		byte[] buffer = null;
		try {
           File f = new File( fileobjName );
           FileInputStream fis = new FileInputStream( f );
            
		   buffer = new byte[12*Nx*Ny*Nz];
		   fis.read(buffer);
           fis.close();
		} catch (IOException io) {
           System.out.println("i/o pb: "+io.getMessage());
		}
		// convert to the image format
		float[][] img  = new float[3][Nx*Ny*Nz];

		for (int xyz=0;xyz<Nx*Ny*Nz;xyz++) for (int d=0;d<3;d++) {
			int b1 = buffer[4*(xyz+Nx*Ny*Nz*d)+0] & 0xff;
			int b2 = buffer[4*(xyz+Nx*Ny*Nz*d)+1] & 0xff;
			int b3 = buffer[4*(xyz+Nx*Ny*Nz*d)+2] & 0xff;
			int b4 = buffer[4*(xyz+Nx*Ny*Nz*d)+3] & 0xff;
			// big endian
			//int tmpInt = ((b1 << 24) | (b2 << 16) | (b3 << 8) | b4);
			// little endian
			int tmpInt = ((b4 << 24) | (b3 << 16) | (b2 << 8) | b1);

			img[d][xyz] = Float.intBitsToFloat(tmpInt);
		}
        buffer = null;
		
		return img;
	}
	
	/** 
	 *	load the atlas data from a file. 
	 *  All associated images are loaded at this time
	 */
	final public void loadAtlas(String fileobjName) {
		if (verbose) System.out.println("loading atlas file: "+fileobjName);
		try {
            File f = new File(fileobjName);
			String dir = f.getParent();
            FileReader fr = new FileReader(f);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
			StringTokenizer st;
			String imageFile, objLabelFile;
            // Exact corresponding template
            if (!line.equals("Structure Atlas File (edit at your own risks)")) {
                System.out.println("not a proper Structure Atlas file");
                br.close();
                fr.close();
                return;
            }
			line = br.readLine();
			while (line!=null) {
				if (line.startsWith("Structures")) {
					//System.out.println(line);
					// Structures:	nobj	objLabel	objType
					st = new StringTokenizer(line, "	");
					st.nextToken();
					nobj = BasicInfo.getInt(st);
					objName = new String[nobj];
					objLabel = new byte[nobj];
					objType = new String[nobj];
					for (int n=0;n<nobj;n++) {
						// Name:objLabel:objType
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						objName[n] = st.nextToken();
						objLabel[n] = (byte)BasicInfo.getInt(st);
						objType[n] = st.nextToken();
					}
					// allocate other quantities
					shape = new float[nobj][];
					minx = new int[nobj];
					miny = new int[nobj];
					minz = new int[nobj];
					maxx = new int[nobj];
					maxy = new int[nobj];
					maxz = new int[nobj];
					shapeFile = new String[nobj];
					for (int n=0;n<nobj;n++) shapeFile[n] = null;
					templateFile = null;
					intensity = new float[INTENSITY][nobj];
					intensityMap = new float[INTENSITY][][];
					registeredShape = new boolean[nobj];
					regularizationFactor = new float[nobj];
					for (int n=0;n<nobj;n++) registeredShape[n] = true;
					registeredShape[0] = false;
					for (int n=0;n<nobj;n++) regularizationFactor[n] = 1.0f;
					if (debug) System.out.println(displayNames());
				} else
				if (line.startsWith("Topology Atlas")) {
					//System.out.println(line);
					// File: objName
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					imageFile = dir+File.separator+st.nextToken();
					if (debug) System.out.print("file: "+imageFile+"\n");
					// Dimensions: ntx nty ntz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					ntx = BasicInfo.getInt(st);
					nty = BasicInfo.getInt(st);
					ntz = BasicInfo.getInt(st);
					x0t = ntx/2.0f;
					y0t = nty/2.0f;
					z0t = ntz/2.0f;
					if (debug) System.out.print("dims: "+ntx+" "+nty+" "+ntz+"\n");
					// Resolutions: rtx rty rtz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					rtx = BasicInfo.getFloat(st);
					rty = BasicInfo.getFloat(st);
					rtz = BasicInfo.getFloat(st);
					if (debug) System.out.print("res: "+rtx+"x"+rty+"x"+rtz+"\n");
					template = loadTemplateImage(imageFile, ntx, nty, ntz);
					templateFile = imageFile;
				} else
				if (line.startsWith("Shape Atlas")) {
					//if (debug) System.out.println(line);
					// Shape:	objLabelSamples
					st = new StringTokenizer(line, "	");
					st.nextToken();
					objLabelSamples = BasicInfo.getInt(st);
					// Dimensions: nax nay naz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					nax = BasicInfo.getInt(st);
					nay = BasicInfo.getInt(st);
					naz = BasicInfo.getInt(st);
					x0a = nax/2.0f;
					y0a = nay/2.0f;
					z0a = naz/2.0f;
					if (debug) System.out.print("shape dim: "+nax+"x"+nay+"x"+naz+"\n");
					// Resolutions: rax ray raz
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					st.nextToken();
					rax = BasicInfo.getFloat(st);
					ray = BasicInfo.getFloat(st);
					raz = BasicInfo.getFloat(st);
					if (debug) System.out.print("shape res: "+rax+"x"+ray+"x"+raz+"\n");
					line = br.readLine();
					while (line.startsWith("Structure:")) {
						// find structure id
						st = new StringTokenizer(line, "	");
						st.nextToken();
						String title = st.nextToken();
						int id=-1;
						for (int n=0;n<nobj;n++) {
							if (title.equals(objName[n])) { id = n; break; }
						}
						if (debug) System.out.print("Shape: "+objName[id]+"\n");
						if (id==-1) {
							line = br.readLine();
						} else {
							// Proba: objName
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							st.nextToken();
							imageFile = dir+File.separator+st.nextToken();
							// min, max : initial values
							minx[id] = 0; miny[id] = 0; minz[id] = 0;
							maxx[id] = nax; maxy[id] = nay; maxz[id] = naz;
			
							shape[id] = loadShapeImage(imageFile, nax, nay, naz);
							shapeFile[id] = imageFile;
						}
						line = br.readLine();
					}
				} else
				if (line.startsWith("Intensity Atlas")) {
					line = br.readLine();
					// Type value value value...
					line = br.readLine();
					st = new StringTokenizer(line, "	");
					String type = st.nextToken();
					//if (debug) System.out.println(type);
					int id = contrastId(type);
					
					while (id !=-1 ) {
						for (int n=0;n<nobj;n++) {
							intensity[id][n] = BasicInfo.getFloat(st);
						}
						// search for next intensity profile
						line = br.readLine();
						st = new StringTokenizer(line, "	");
						type = st.nextToken();
						//if (debug) System.out.println(type);
						id =contrastId(type);
					}
					if (debug) System.out.println(displayIntensity());
				} else
				if (line.startsWith("Registered Shapes")) {
					if (debug) System.out.println(line);
					// Type value value
					line = br.readLine();
					if (!line.startsWith("0") && !line.startsWith("1")) line = br.readLine();
					//if (debug) System.out.println(line);
					st = new StringTokenizer(line, "	");
					for (int n=0;n<nobj;n++) {
						registeredShape[n] = (BasicInfo.getInt(st)==1);
					}
					if (debug) System.out.println(displayRegisteredShapes());
				} else
				if (line.startsWith("Regularization Factor")) {
					if (debug) System.out.println(line);
					// structure structure ... (not read)
					line = br.readLine();
					// value value...
					line = br.readLine();
					//if (debug) System.out.println(line);
					st = new StringTokenizer(line, "	");
					for (int n=0;n<nobj;n++) {
						regularizationFactor[n] = BasicInfo.getFloat(st);
					}
					if (debug) System.out.println(displayRegularizationFactor());
				} else
				if (line.startsWith("Intensity Prior:")) {
					st = new StringTokenizer(line, "	");
					st.nextToken();
					String title = st.nextToken();
					if (title.equals("T1map") || title.equals("T1map7T")) {
						if (debug) System.out.println("T1map 7T prior");
						// fill in t1 map values
						intensityMap[T1MAP7T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[T1MAP7T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[T1MAP7T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[T1MAP7T][id] = intensityMap[T1MAP7T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(T1MAP7T));
					} else
					if (title.equals("Iso") || title.equals("FlatImg") || title.equals("Mp2rage7T")) {
						if (debug) System.out.println("Mp2rage 7T prior");
						// fill in flat img values
						intensityMap[MP2RAGE7T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[MP2RAGE7T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[MP2RAGE7T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[MP2RAGE7T][id] = intensityMap[MP2RAGE7T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(MP2RAGE7T));
					} else
					if (title.equals("Flash7T") || title.equals("T2Sw7T") || title.equals("T2SW7T") || title.equals("T2sw7T")) {
						if (debug) System.out.println("T2*-weighted 7T prior");
						// fill in flat img values
						intensityMap[T2SW7T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[T2SW7T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[T2SW7T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[T2SW7T][id] = intensityMap[T2SW7T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(T2SW7T));
					} else
					if (title.equals("QSM7T") || title.equals("qsm7T")) {
						if (debug) System.out.println("QSM 7T prior");
						// fill in flat img values
						intensityMap[QSM7T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[QSM7T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[QSM7T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[QSM7T][id] = intensityMap[QSM7T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(QSM7T));
					} else
					if (title.equals("T1map3T")) {
						if (debug) System.out.println("T1map 3T prior");
						// fill in t1 map values
						intensityMap[T1MAP3T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[T1MAP3T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[T1MAP3T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[T1MAP3T][id] = intensityMap[T1MAP3T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(T1MAP3T));
					} else
					if (title.equals("Mp2rage3T")) {
						if (debug) System.out.println("Mp2rage 3T prior");
						// fill in flat img values
						intensityMap[MP2RAGE3T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[MP2RAGE3T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[MP2RAGE3T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[MP2RAGE3T][id] = intensityMap[MP2RAGE3T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(MP2RAGE3T));
					} else
					if (title.equals("T1w") || title.equals("Mprage3T")) {
						if (debug) System.out.println("MPRAGE 3T prior");
						// fill in flat img values
						intensityMap[MPRAGE3T] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[MPRAGE3T][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[MPRAGE3T][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[MPRAGE3T][id] = intensityMap[MPRAGE3T][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(MPRAGE3T));
					} else
					if (title.equals("NormMPRAGE")) {
						if (debug) System.out.println("B1-normalized MPRAGE 3T prior");
						// fill in flat img values
						intensityMap[NORMMPRAGE] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[NORMMPRAGE][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[NORMMPRAGE][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[NORMMPRAGE][id] = intensityMap[NORMMPRAGE][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(NORMMPRAGE));
					} else
					if (title.equals("PV") || title.equals("PVdura") || title.equals("PVDURA")) {
						if (debug) System.out.println("signed PV priors");
						// fill in values
						intensityMap[PVDURA] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3], [value1,2,3] ...
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[PVDURA][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[PVDURA][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[PVDURA][id] = intensityMap[PVDURA][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(PVDURA));
					} else
					if (title.equals("Filters")) {
						if (debug) System.out.println("filter priors");
						// fill in values
						intensityMap[FILTER] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3,...]
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[FILTER][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[FILTER][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[FILTER][id] = intensityMap[FILTER][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(FILTER));
					} else
					if (title.equals("HCPT1w")) {
						if (debug) System.out.println("HCP T1w priors");
						// fill in values
						intensityMap[HCPT1W] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3,...]
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[HCPT1W][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[HCPT1W][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[HCPT1W][id] = intensityMap[HCPT1W][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(HCPT1W));
					} else
					if (title.equals("HCPT2w")) {
						if (debug) System.out.println("HCP T2w priors");
						// fill in values
						intensityMap[HCPT2W] = new float[nobj][];
						for (int n=0;n<nobj;n++) {
							// label [value1,2,3,...]
							line = br.readLine();
							st = new StringTokenizer(line, "	");
							String name = st.nextToken();
							//if (debug) System.out.println(name);
							int id=-1;
							for (int m=0;m<nobj;m++) {
								if (name.equals(objName[n])) { id = n; break; }
							}
							//if (debug) System.out.print("Structure: "+objName[id]+"\n");
							if (id==-1) {
								continue;
							} else {
								int count = st.countTokens();
								intensityMap[HCPT2W][id] = new float[count];
								for (int t=0;t<count;t++)
									intensityMap[HCPT2W][id][t] = BasicInfo.getFloat(st);
								// set intensity maps to the first value
								intensity[HCPT2W][id] = intensityMap[HCPT2W][id][0];
							}
						}
						if (debug) System.out.println(displayMapIntensity(HCPT2W));
					}
				}

				line = br.readLine();
				if (debug) System.out.println(line);
			}		
			br.close();
            fr.close();
			atlasFile = fileobjName;
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

		if (debug) BasicInfo.displayMessage("initialisation\n");
	}

	public static final int contrastId(String type) {
		int id = -1;
		if (type.equals("T1w") || type.equals("T1_SPGR") || type.equals("T1_MPRAGE") || type.equals("MPRAGE3T")) id = MPRAGE3T;
		else if (type.equals("T2Sw") || type.equals("T2*") || type.equals("T2SW7T") || type.equals("T2sw7T") || type.equals("FLASH7T")) id = T2SW7T;
		else if (type.equals("ISO") || type.equals("Iso") || type.equals("FlatImg") || type.equals("MP2RAGE7T")) id = MP2RAGE7T;
		else if (type.equals("T1M") || type.equals("T1map") || type.equals("T1map7T") || type.equals("T1MAP7T")) id = T1MAP7T;
		else if (type.equals("MP2RAGE3T") || type.equals("Mp2rage3T")) id = MP2RAGE3T;
		else if (type.equals("T1map3T") || type.equals("T1MAP3T")) id = T1MAP3T;
		else if (type.equals("sPV") || type.equals("PV") || type.equals("PV_dura") || type.equals("PVDURA")) id = PVDURA;
		else if (type.equals("Filters") || type.equals("FILTER")) id = FILTER;
		else if (type.equals("HCPT1w") || type.equals("HCPT1W")) id = HCPT1W;
		else if (type.equals("HCPT2w") || type.equals("HCPT2W")) id = HCPT2W;
		else if (type.equals("QSM7T") || type.equals("qsm7T")) id = QSM7T;
		else if (type.equals("NORMMPRAGE") || type.equals("NormMPRAGE")) id = NORMMPRAGE;
		
		return id;
	}					
		
	public static final String displayContrastName(int id) {
		if (id==MPRAGE3T) return "MPRAGE3T";
		else if (id==T2SW7T) return "T2SW7T";
		else if (id==MP2RAGE7T) return "MP2RAGE7T";
		else if (id==T1MAP7T) return "T1MAP7T";
		else if (id==MP2RAGE3T) return "MP2RAGE3T";
		else if (id==T1MAP3T) return "T1MAP3T";
		else if (id==PVDURA) return "PVDURA";
		else if (id==FILTER) return "FILTER";
		else if (id==HCPT1W) return "HCPT1W";
		else if (id==HCPT2W) return "HCPT2W";
		else if (id==QSM7T) return "QSM7T";
		else if (id==NORMMPRAGE) return "NORMMPRAGE";
		else return "none";
	}					
		
	public static final boolean isIntensityContrast(int id) {
		if (id==MPRAGE3T || id==T2SW7T || id==MP2RAGE7T || id==T1MAP7T || id==MP2RAGE3T || id==T1MAP3T || id==HCPT1W || id==HCPT2W || id==QSM7T || id==NORMMPRAGE) return true;
		else return false;
	}					
		
	public static final boolean isQuantitativeContrast(int id) {
		if (id==T1MAP7T || id==T1MAP3T || id==QSM7T) return true;
		else return false;
	}					
		
	public static final boolean isProbabilityContrast(int id) {
		if (id==PVDURA || id==FILTER) return true;
		else return false;
	}					
		
	/** transformations: re-compute the template using the transform
	 */
	public final void computeTransformedTemplate() {
		float[] XP=new float[3];
		byte[] tmp = new byte[nix*niy*niz];
		
		boolean unknownLabel = false;
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			transformModel.imageToTemplate(XP,x,y,z,transform,rotation,1.0f);
			
			tmp[x+y*nix+z*nix*niy] = ImageInterpolation.nearestNeighborInterpolation(template,objLabel[0],XP[0],XP[1],XP[2],ntx,nty,ntz);
			// debug
			boolean wrong=true;
			for (int k=0;k<nobj;k++) {
				if (tmp[x+y*nix+z*nix*niy]==objLabel[k]) { wrong=false; }
			}
			if (wrong) unknownLabel = true;
		}
		if (unknownLabel) System.out.println("warning: incorrect objLabel image \n");

		template = tmp;
		
		return;
	}
					
    public final void alignObjectCenter(float[] img, String type) {
    	float xi = 0.0f, yi = 0.0f, zi = 0.0f, wi = 0.0f;
    	float xs = 0.0f, ys = 0.0f, zs = 0.0f, ws = 0.0f;
    	int xyz;
    	float pobj;
    	
    	boolean[] isobj = new boolean[nobj];
    	for (int n=0;n<nobj;n++) 
    		if (objType[n].startsWith(type)) isobj[n] = true;
    		else isobj[n] = false;
    	
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		xi += img[xyz]*x;
    		yi += img[xyz]*y;
    		zi += img[xyz]*z;
    		wi += img[xyz];
    		
    		pobj = 0.0f;
    		for (int n=0;n<nobj;n++) 
    			if (isobj[n]) 
    				pobj = Numerics.max(pobj,shape[n][xyz]);
    			
    		xs += pobj*x;
    		ys += pobj*y;
    		zs += pobj*z;
    		ws += pobj;
    	}
    	
    	xi /= wi;
    	yi /= wi;
    	zi /= wi;
    	
    	xs /= ws;
    	ys /= ws;
    	zs /= ws;
    	
    	// translation
    	transform[3] += (xs-xi)*rax;
    	transform[4] += (ys-yi)*ray;
    	transform[5] += (zs-zi)*raz;
    	
    	return;
    }
    public final float[] generateObjectImage(String type) {
    	float[] img = new float[nax*nay*naz];
    	boolean[] isobj = new boolean[nobj];
    	int xyz;
    	for (int n=0;n<nobj;n++) 
    		if (objType[n].startsWith(type)) isobj[n] = true;
    		else isobj[n] = false;
    	
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		img[xyz] = 0.0f;
    		for (int n=0;n<nobj;n++) if (isobj[n])
    			img[xyz] = Numerics.max(img[xyz],shape[n][xyz]);    			
    	}
    	return img;
    }
    public final float[] generateObjectSegmentation(String type) {
    	float[] img = new float[nax*nay*naz];
    	boolean[] isobj = new boolean[nobj];
    	int xyz;
    	float in,out;
    	for (int n=0;n<nobj;n++) 
    		if (objType[n].startsWith(type)) isobj[n] = true;
    		else isobj[n] = false;
    	
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		in = 0.0f;
    		out = 0.0f;
    		for (int n=0;n<nobj;n++) {
    			if (isobj[n]) in = Numerics.max(in,shape[n][xyz]);  
				else out = Numerics.max(out,shape[n][xyz]);  
			}
			if (in>out) img[xyz] = 1.0f;
			else img[xyz] = 0.0f;
    	}
    	return img;
    }
    public final float[] generateObjectSegmentation(String ptype, String ntype) {
    	float[] img = new float[nax*nay*naz];
    	boolean[] ispobj = new boolean[nobj];
    	boolean[] isnobj = new boolean[nobj];
    	int xyz;
    	float pos,neg,out;
    	for (int n=0;n<nobj;n++) {
    		if (objType[n].startsWith(ptype)) ispobj[n] = true;
    		else ispobj[n] = false;
    		if (objType[n].startsWith(ntype)) isnobj[n] = true;
    		else isnobj[n] = false;
    	}
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		pos = 0.0f;
    		neg = 0.0f;
    		out = 0.0f;
    		for (int n=0;n<nobj;n++) {
    			if (ispobj[n]) pos = Numerics.max(pos,shape[n][xyz]);  
				else if (isnobj[n]) neg = Numerics.max(neg,shape[n][xyz]);  
				else out = Numerics.max(out,shape[n][xyz]);  
			}
			if (pos>out && pos>neg) img[xyz] = 1.0f;
			else if (neg>out && neg>pos) img[xyz] = 2.0f;
			else img[xyz] = 0.0f;
    	}
    	return img;
    }
    public final float[] generateDifferentialObjectSegmentation(String ptype, String ntype) {
    	float[] img = new float[nax*nay*naz];
    	boolean[] ispobj = new boolean[nobj];
    	boolean[] isnobj = new boolean[nobj];
    	int xyz;
    	float pos,neg,out;
    	for (int n=0;n<nobj;n++) {
    		if (objType[n].startsWith(ptype)) ispobj[n] = true;
    		else ispobj[n] = false;
    		if (objType[n].startsWith(ntype)) isnobj[n] = true;
    		else isnobj[n] = false;
    	}
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		pos = 0.0f;
    		neg = 0.0f;
    		out = 0.0f;
    		for (int n=0;n<nobj;n++) {
    			if (ispobj[n]) pos = Numerics.max(pos,shape[n][xyz]);  
				else if (isnobj[n]) neg = Numerics.max(neg,shape[n][xyz]);  
				else out = Numerics.max(out,shape[n][xyz]);  
			}
			if (pos>out && pos>neg) img[xyz] = 2.0f;
			else if (neg>out && neg>pos) img[xyz] = 0.0f;
			else img[xyz] = 1.0f;
    	}
    	return img;
    }
    /** 
	 *  generate atlas image from information
	 */
    final public byte[] generateTransformedRestrictedClassification(String ptype, String ntype) {
		float dist,max,count,val;
		byte best;
		byte[] img = new byte[nix*niy*niz];
		float[] XP=new float[3];
		
		for (int x=0;x<nix;x++) for (int y=0;y<niy;y++) for (int z=0;z<niz;z++) {
			imageToShapeCoordinates(XP, x,y,z);
			
			// compute each class probability : attribute the highest
			max = 0; best = -1;
			for (byte k=0;k<nobj;k++) {
				val = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
				if (val>max) {
					best = k;
					max = val;
				}
			}
			if (best>-1 && (objType[best].startsWith(ptype) || objType[best].startsWith(ntype) ) )
				img[x+y*nix+z*nix*niy] = best;
			else img[x+y*nix+z*nix*niy] = -1;
		}
		return img;
	}
    public final float[] generateObjectSegmentation(String type1, String type2, String type3) {
    	float[] img = new float[nax*nay*naz];
    	boolean[] isobj1 = new boolean[nobj];
    	boolean[] isobj2 = new boolean[nobj];
    	boolean[] isobj3 = new boolean[nobj];
    	int xyz;
    	float val1, val2, val3, out;
    	for (int n=0;n<nobj;n++) {
    		if (objType[n].startsWith(type1)) isobj1[n] = true;
    		else isobj1[n] = false;
    		if (objType[n].startsWith(type2)) isobj2[n] = true;
    		else isobj2[n] = false;
    		if (objType[n].startsWith(type3)) isobj3[n] = true;
    		else isobj3[n] = false;
    	}
    	for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
    		xyz = x+nax*y+nax*nay*z;
    		val1 = 0.0f;
    		val2 = 0.0f;
    		val3 = 0.0f;
    		out = 0.0f;
    		for (int n=0;n<nobj;n++) {
    				 if (isobj1[n]) val1 = Numerics.max(val1,shape[n][xyz]);  
				else if (isobj2[n]) val2 = Numerics.max(val2,shape[n][xyz]);  
				else if (isobj3[n]) val3 = Numerics.max(val3,shape[n][xyz]);  
				else out = Numerics.max(out,shape[n][xyz]);  
			}
				 if (val1>val2 && val1>val3 && val1>out) img[xyz] = 1.0f;
			else if (val2>val1 && val2>val3 && val2>out) img[xyz] = 2.0f;
			else if (val3>val2 && val3>val1 && val3>out) img[xyz] = 3.0f;
			else img[xyz] = 0.0f;
    	}
    	return img;
    }
	/**
	 *	initialize registration parameters
	 */
	public final void initShapeMapping() {
		// transform mapping tools
		transformMode = PARAMETRIC;
		transformModel = new ParametricTransform("rigid", x0i,y0i,z0i, rix,riy,riz, nix,niy,niz, x0t,y0t,z0t, rtx,rty,rtz, ntx,nty,ntz);

		// init parameters		
		Nd = transformModel.getDimension();
		initRigidTransform();
		rotation = transformModel.computeRotation(transform);
			
		shapeTransform = new float[3][4];
		precomputeTransformMatrix(1.0f);	
		inverseShapeTransform = new float[3][4];
		precomputeInverseTransformMatrix(1.0f);	
	}		
	private final void initRigidTransform() {
		// orientation: initialize a rotation
		transform = new float[6];
		// note : assumes a rotation around the image center
		if (orient==BasicInfo.AXIAL) {
			transform[0] = 0.0f;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
			
			if (orix==BasicInfo.L2R) rix *= -1.0f;
			if (oriy==BasicInfo.P2A) riy *= -1.0f;
			if (oriz==BasicInfo.S2I) riz *= -1.0f;
		} else if (orient==BasicInfo.CORONAL) {
			transform[0] = -ISQRT2;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
			
			if (orix==BasicInfo.L2R) rix *= -1.0f;
			if (oriy==BasicInfo.I2S) riy *= -1.0f;
			if (oriz==BasicInfo.P2A) riz *= -1.0f;
		} else if (orient==BasicInfo.SAGITTAL) {
			transform[0] = -0.5f;
			transform[1] = -0.5f;
			transform[2] = 0.5f;
			
			if (orix==BasicInfo.P2A) rix *= -1.0f;
			if (oriy==BasicInfo.I2S) riy *= -1.0f;
			if (oriz==BasicInfo.R2L) riz *= -1.0f;
		} else {
			// default is axial
			transform[0] = 0.0f;
			transform[1] = 0.0f;
			transform[2] = 0.0f;
		}
		transform[3] = 0.0f;
		transform[4] = 0.0f;
		transform[5] = 0.0f;
	}
	
	public final void updateRigidTransform(float[] trans) {
		// 1. compose rotations
		float[][] rot0 = transformModel.computeRotation(transform);
		float[][] rot1 = transformModel.computeRotation(trans);
		
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			rotation[i][j] = 0.0f;	
			for (int k=0;k<3;k++) rotation[i][j] += rot1[i][k]*rot0[k][j];
			trans[3+i] += rot1[i][j]*transform[3+j];
		}
		RotationMatrix rmat = new RotationMatrix();
		rmat.setMatrix(rotation);
		
		for (int i=0;i<3;i++) {
			transform[i] = rmat.getParameter(i);
			transform[3+i] = trans[3+i];
		}
	}
	
	public final void updateNonRigidTransform(BasicDemonsWarping warp) {
		float[][] deformed = new float[nobj][nax*nay*naz];
		
		int xyz;
		float[] XP = new float[3];
		// deform the atlas to map tt the transform
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
			xyz = x+nax*y+nax*nay*z;
    		warp.getCurrentMapping(XP, x,y,z);
			for (int k=0;k<nobj;k++) {
				deformed[k][xyz] = ImageInterpolation.linearInterpolation(shape[k],0.0f,XP[0],XP[1],XP[2],nax,nay,naz);	
			}
		}
		for (xyz=0;xyz<nax*nay*naz;xyz++) for (int k=0;k<nobj;k++) {
			shape[k][xyz] = deformed[k][xyz];
		}
		deformed = null;
	}
	
	/**
	 *	recompute registration parameters
	 */
	public final void refreshShapeMapping() {
		
		rotation = transformModel.computeRotation(transform);
			
		precomputeTransformMatrix(1.0f);	
		precomputeInverseTransformMatrix(1.0f);	
	}		
	
	/** 
	 *	normalizes the priors into memberships
	 */
	public final void normalizeShapePriors() {
		for (int x=0;x<nax;x++) for (int y=0;y<nay;y++) for (int z=0;z<naz;z++) {
            float sum=0.0f;
			for (int k=0;k<nobj;k++) {
				sum += shape[k][x+nax*y+nax*nay*z];
			}
			if (sum>ZERO) for (int k=0;k<nobj;k++) {
				shape[k][x+nax*y+nax*nay*z] /= sum;
			}
		}
		return;
	}
	
	public final void imageToShapeCoordinates(float[] X, float x,float y,float z) {
		X[0] = shapeTransform[0][0]*x + shapeTransform[0][1]*y + shapeTransform[0][2]*z + shapeTransform[0][3];
		X[1] = shapeTransform[1][0]*x + shapeTransform[1][1]*y + shapeTransform[1][2]*z + shapeTransform[1][3];
		X[2] = shapeTransform[2][0]*x + shapeTransform[2][1]*y + shapeTransform[2][2]*z + shapeTransform[2][3];
		
		return;
	}
	public final void precomputeTransformMatrix(float scale) {
		float[][] rot = null;
		if (transformModel.useRotation())
			rot = transformModel.computeRotation(transform);
		
		transformModel.precomputeImageToTemplateMatrix(shapeTransform, transform, rot, scale);
	}

	public final void shapeToImageCoordinates(float[] X, int x,int y,int z) {
		X[0] = inverseShapeTransform[0][0]*x + inverseShapeTransform[0][1]*y + inverseShapeTransform[0][2]*z + inverseShapeTransform[0][3];
		X[1] = inverseShapeTransform[1][0]*x + inverseShapeTransform[1][1]*y + inverseShapeTransform[1][2]*z + inverseShapeTransform[1][3];
		X[2] = inverseShapeTransform[2][0]*x + inverseShapeTransform[2][1]*y + inverseShapeTransform[2][2]*z + inverseShapeTransform[2][3];
		
		return;
	}
	public final void precomputeInverseTransformMatrix(float scale) {
		float[][] rot = null;
		if (transformModel.useRotation())
			rot = transformModel.computeRotation(transform);
		
		transformModel.precomputeTemplateToImageMatrix(inverseShapeTransform, transform, rot, scale);
	}

	public final String displayTransform(float[] trans) {
		String info = "transform: (";
		for (int n=0;n<Nd-1;n++) info += trans[n]+", ";
		info += trans[Nd-1]+")\n";
		
		return info;
	}
	
	/** display the atlas data */
	public final String displayVector(float[] vect) {
		String info = "vector: (";
		for (int n=0;n<vect.length-1;n++) info += vect[n]+", ";
		info += vect[vect.length-1]+")\n";
		
		return info;
	}
	
		/** display the atlas data */
	final public String displayIntensity() {
		String output = "Intensity \n";
		
		output += "MPRAGE 3T : ";
		for (int k=0;k<nobj;k++) output += intensity[MPRAGE3T][k]+" ";
		output += "\n";
		output += "MP2RAGE 7T : ";
		for (int k=0;k<nobj;k++) output += intensity[MP2RAGE7T][k]+" ";
		output += "\n";
		output += "T1 map 7T : ";
		for (int k=0;k<nobj;k++) output += intensity[T1MAP7T][k]+" ";
		output += "\n";
		output += "MP2RAGE 3T : ";
		for (int k=0;k<nobj;k++) output += intensity[MP2RAGE3T][k]+" ";
		output += "\n";
		output += "T1 map 3T : ";
		for (int k=0;k<nobj;k++) output += intensity[T1MAP3T][k]+" ";
		output += "\n";
		output += "T2*-w 7T : ";
		for (int k=0;k<nobj;k++) output += intensity[T2SW7T][k]+" ";
		output += "\n";
		output += "PV/dura : ";
		for (int k=0;k<nobj;k++) output += intensity[PVDURA][k]+" ";
		output += "\n";
		output += "Filters : ";
		for (int k=0;k<nobj;k++) output += intensity[FILTER][k]+" ";
		output += "\n";
		output += "HCP T1w : ";
		for (int k=0;k<nobj;k++) output += intensity[HCPT1W][k]+" ";
		output += "\n";
		output += "HCP T2w : ";
		for (int k=0;k<nobj;k++) output += intensity[HCPT2W][k]+" ";
		output += "\n";
		
		return output;	
	}
	
	/** display the atlas data */
	final public String displayMapIntensity(int modal) {
		String output = displayContrastName(modal)+" Intensity Prior \n";
		
		for (int n=0;n<nobj;n++) {
			output += objName[n]+" : ";
			for (int k=0;k<intensityMap[modal][n].length;k++) output += intensityMap[modal][n][k]+" ";
			output += "\n";
		}
		
		return output;	
	}


}
