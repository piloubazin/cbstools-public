package de.mpg.cbs.libraries;


/**
 *
 *  This class computes the numbers of connected components of higher or lower values
 *  in the neighborhood of a point, with various connectivities.
 *
 *	@version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class NeighborConnectivity {
	
	// no data: used as a library of functions
	
	public final static int SUPERIOR = 0;
	public final static int SUPEQUAL = 1;
	public final static int INFERIOR = 2;
	public final static int INFEQUAL = 3;
	public final static int NONE = 4;
	
	/** 
	 *  builds a boolean object from the case in hand
	 */
	public static int countNeighborhoods(byte[][][] image, int x, int y, int z, int type, int connectivity, int dz) {
		
		boolean obj[][][] = objectFromNeighborhood(image,x,y,z,type);
		int N = connectedNeighborhood(obj,connectivity,dz);
		obj = null;
		
		return N;
	}
	/** 
	 *  builds a boolean object from the case in hand
	 */
	public static boolean[][][] objectFromNeighborhood(byte[][][] image, int x, int y, int z, int type) {
		boolean obj[][][] = new boolean[3][3][3];
		int i,j,l;
	
		for (i=0;i<3;i++) for (j=0;j<3;j++) for (l=0;l<3;l++) {
			obj[i][j][l] = true;
            if ( (type==INFERIOR) && (image[x+i-1][y+j-1][z+l-1] >=image[x][y][z]) ) obj[i][j][l] = false;
            if ( (type==INFEQUAL) && (image[x+i-1][y+j-1][z+l-1] > image[x][y][z]) ) obj[i][j][l] = false;
            if ( (type==SUPERIOR) && (image[x+i-1][y+j-1][z+l-1] <=image[x][y][z]) ) obj[i][j][l] = false;
            if ( (type==SUPEQUAL) && (image[x+i-1][y+j-1][z+l-1] < image[x][y][z]) ) obj[i][j][l] = false;
		}
		obj[1][1][1] = false;
		//System.out.println("o");
			
		return obj;
	}

   /**
     *  connected neighborhood: handles all cases
     */
    public static int connectedNeighborhood(boolean img[][][], int connectivity, int dz) {
		int Nlabel = 0;
		int[][][]   label = new int[3][3][3];
		int[]       lb = new int[26];
		int lbMin;
		int x,y,z,i,j,k,l,c,n;
		int Nlb;
		int[]   connect = new int[26];
		int Nconnect;
		int AddLabel;
		
		for (x=0;x<3;x++) for (y=0;y<3;y++) for (z=0;z<3;z++)
			label[x][y][z] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (x=0;x<3;x++) for (y=0;y<3;y++) for (z=0;z<3;z++) {
			if (img[x][y][z]) {
				// object point: neighbors ?
				Nconnect = 0;
				for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
					if ( (i*i+j*j+l*l < connectivity) && (i*i+j*j+l*l > 0) ){
						if ( (x+i>=0) && (x+i<3) && (y+j>=0) && (y+j<3) && (z+l>=0) && (z+l<3) ) {
							if (label[x+i][y+j][z+l] > 0) {
								connect[Nconnect] = lb[ label[x+i][y+j][z+l] ];
								Nconnect++;
							}
						}
					}
				}
				// if connected values, find the smallest lb label and attribute it
				// to all others (-> join labels)
				if (Nconnect>0) {
					lbMin = lb[connect[0]];
					for (l=0;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
					for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
					label[x][y][z] = lbMin;
				} else {
					// new, unconnected region
					label[x][y][z] = Nlabel;
					lb[Nlabel] = Nlabel;
					Nlabel++;
				}
			}
        }
        // only one level of labels
        for (k=1;k<Nlabel;k++) {
            c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
		// count the valid labels
		Nlb = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				// check the labels for neighbors only 
				AddLabel = 0;
				for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (l=-dz;l<=dz;l++) {
					if ( (i*i+j*j+l*l < connectivity) && (i*i+j*j+l*l > 0) ) {
						if (lb[label[i+1][j+1][l+1]] == k) AddLabel=1;
					}
				}
				if (AddLabel>0) Nlb++;
			}
		}
        // clean up
        lb = null;
        connect = null;
		label = null;
	   
		return Nlb;
    }

	/** 
	 *  2D images: 4-connectivity
	 */
	public static int connected4Neighborhood2D(boolean img[][]) {
		int Nlabel;
		int[][] label = new int[3][3];
		int[]   lb = new int[9];
		int lbMin;
		int i,j,k,c;
		int Nlb;
		int[]   connect = new int[4];
		int Nconnect;
		int AddLabel=0;
	
		// the input is a 3x3 binary image (0 out, 1 in)
		for (i=0;i<3;i++)
			for (j=0;j<3;j++)
				label[i][j] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (i=0;i<3;i++) {
			for (j=0;j<3;j++) {
				if (img[i][j]) {
					// object point: neighbors ?
					Nconnect = 0;
					if (i>0)
						if (label[i-1][j] > 0) {
							connect[Nconnect] = lb [label[i-1][j] ];
							Nconnect++;
						}
					if (j>0) 
						if (label[i][j-1] > 0) {
							connect[Nconnect] = lb[ label[i][j-1] ];
							Nconnect++;
						}
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						lbMin = lb[connect[0]];
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb[connect[k]]);
						for (k=0;k<Nconnect;k++) lb[connect[k]] = lbMin;
						label[i][j] = lbMin;
					} else {
						// new, unconnected region
						label[i][j] = Nlabel;
						lb[Nlabel] = Nlabel;
						Nlabel++;
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels
		Nlb = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				// check the labels for 4-neighbors only
				AddLabel = 0;
				if (lb[label[0][1]] == k) AddLabel=1;
				if (lb[label[1][0]] == k) AddLabel=1;
				if (lb[label[2][1]] == k) AddLabel=1;
				if (lb[label[1][2]] == k) AddLabel=1;
				if (AddLabel>0) Nlb++;
			}
		}
        // clean up
        label = null;
        lb = null;
        connect = null;
	   
		return Nlb;
	}
	
	/** 
	 *  2D images: 4-connectivity, positive  neighbors
	 */
	public static int positive4Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
		int i,j;
	
		for (i=0;i<3;i++) for (j=0;j<3;j++) {
			if (image[x+i-1][y+j-1][z] > image[x][y][z])
				img[i][j] = true;
			else
				img[i][j] = false;
		}
			
		return connected4Neighborhood2D(img);
	}
	/** 
	 *  2D images: 4-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual4Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if ( ( (i!=1) || (j!=1) ) && (image[x+i-1][y+j-1][z] >= image[x][y][z]) ) 
				img[i][j] = true;
			else
				img[i][j] = false;
		}
		return connected4Neighborhood2D(img);
	}
	/** 
	 *  2D images: 4-connectivity, negative neighbors
	 */
	public static int negative4Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];

		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if (image[x+i-1][y+j-1][z] < image[x][y][z])
				img[i][j] = true;
			else
				img[i][j] = false;
        }
		return connected4Neighborhood2D(img);
	}
	/** 
	 *  2D images: 4-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual4Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
		
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if ( ( (i!=1) || (j!=1) ) && (image[x+i-1][y+j-1][z] <= image[x][y][z]) )
				img[i][j] = true;
			else
				img[i][j] = false;
		}
			
		return connected4Neighborhood2D(img);
	}


	/** 
	 *  2D images: 8-neighborhood 
	 */
	public static int connected8Neighborhood2D(boolean img[][]) {
		int Nlabel;
		int[][] label = new int[3][3];
		int[]   lb = new int[9];
		int lbMin;
		int i,j,k,c;
		int Nlb;
		int[]   connect = new int[8];
		int Nconnect;
		int AddLabel=0;
	
		for (i=0;i<3;i++)
			for (j=0;j<3;j++)
				label[i][j] = 0;
		
		lb[0] = 0;
		Nlabel = 1;
		for (i=0;i<3;i++) {
			for (j=0;j<3;j++) {
				if (img[i][j]) {
					// object point: neighbors ?
					Nconnect = 0;
					if (i>0) 
						if (label[i-1][j] > 0) {
							connect[Nconnect] = lb[ label[i-1][j] ];
							Nconnect++;
						}
					if ( (i>0) && (j>0) )
						if (label[i-1][j-1] > 0) {
							connect[Nconnect] = lb[ label[i-1][j-1] ];
							Nconnect++;
						}
					if ( (i>0) && (j<2) )
						if (label[i-1][j+1] > 0) {
							connect[Nconnect] = lb[ label[i-1][j+1] ];
							Nconnect++;
						}
					if (j>0)
						if (label[i][j-1] > 0) {
							connect[Nconnect] = lb[ label[i][j-1] ];
							Nconnect++;
						}
					// if connected values, find the smallest lb label and attribute it
					// to all others (-> join labels)
					if (Nconnect>0) {
						//printf("c:%d",Nconnect);
						lbMin = lb[connect[0]];
						for (k=1;k<Nconnect;k++) lbMin = Math.min(lbMin,lb[connect[k]]);
						for (k=0;k<Nconnect;k++) lb[connect[k]] = lbMin;
						label[i][j] = lbMin;
					} else {
						// new, unconnected region
						label[i][j] = Nlabel;
						lb[Nlabel] = Nlabel;
						//printf("l:%d", Nlabel);
						Nlabel++;
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels : all neighbors
		Nlb = 0;
		for (k=1;k<Nlabel;k++)
			if (lb[k]==k) Nlb++;
        // clean up
        label = null;
        lb = null;
        connect = null;
		
		return Nlb;
	}
	/** 
	 *  2D images: 8-connectivity, positive  neighbors
	 */
	public static int positive8Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
		int i,j;
	
		for (i=0;i<3;i++) for (j=0;j<3;j++) {
			if (image[x+i-1][y+j-1][z] > image[x][y][z])
				img[i][j] = true;
			else
				img[i][j] = false;
		}
			
		return connected8Neighborhood2D(img);
	}
	/** 
	 *  2D images: 8-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual8Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if ( ( (i!=1) || (j!=1) ) && (image[x+i-1][y+j-1][z] >= image[x][y][z]) ) 
				img[i][j] = true;
			else
				img[i][j] = false;
		}
		return connected8Neighborhood2D(img);
	}
	/** 
	 *  2D images: 8-connectivity, negative neighbors
	 */
	public static int negative8Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];

		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if (image[x+i-1][y+j-1][z] < image[x][y][z])
				img[i][j] = true;
			else
				img[i][j] = false;
        }
		return connected8Neighborhood2D(img);
	}
	/** 
	 *  2D images: 8-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual8Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][] = new boolean[3][3];
		
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
			if ( ( (i!=1) || (j!=1) ) && (image[x+i-1][y+j-1][z] <= image[x][y][z]) )
				img[i][j] = true;
			else
				img[i][j] = false;
		}
			
		return connected8Neighborhood2D(img);
	}


	/** 
	 *  3D images: 6-neighborhood (implies 18-background) 
	 */
	public static int connected618Neighborhood3D(boolean img[][][]) {
		int Nlabel = 0;
		int[][][]   label = new int[3][3][3];
		int[]       lb = new int[27];
		int lbMin;
		int i,j,k,x,y,z,l,c,n;
		int Nlb;
		int[]   connect = new int[6];
		int Nconnect;
		int AddLabel;
		
		for (i=0;i<3;i++)
			for (j=0;j<3;j++)
				for (k=0;k<3;k++)
					label[i][j][k] = 0;
		
		// remove the 26-C corners ?
		img[0][0][0] = false;
		img[0][0][2] = false;
		img[0][2][0] = false;
		img[2][0][0] = false;
		img[0][2][2] = false;
		img[2][0][2] = false;
		img[2][2][0] = false;
		img[2][2][2] = false;
		
		lb[0] = 0;
		Nlabel = 1;
		for (i=0;i<3;i++) {
			for (j=0;j<3;j++) {
				for (k=0;k<3;k++) {
					if (img[i][j][k]) {
						// object point: neighbors ?
						Nconnect = 0;
						if (i>0) 
							if (label[i-1][j][k] > 0) {
								connect[Nconnect] = lb[ label[i-1][j][k] ];
								Nconnect++;
							}
						if (j>0)
							if (label[i][j-1][k] > 0) {
								connect[Nconnect] = lb[ label[i][j-1][k] ];
								Nconnect++;
							}
						if (k>0)
							if (label[i][j][k-1] > 0) {
								connect[Nconnect] = lb[ label[i][j][k-1] ];
								Nconnect++;
							}
						// if connected values, find the smallest lb label and attribute it
						// to all others (-> join labels)
						if (Nconnect>0) {
							//printf("c:%d",Nconnect);
							lbMin = lb[connect[0]];
							for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
							for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
							label[i][j][k] = lbMin;
						} else {
							// new, unconnected region
							label[i][j][k] = Nlabel;
							lb[Nlabel] = Nlabel;
							//printf("l:%d", Nlabel);
							Nlabel++;
						}
					}
				}
			}
		}
		// only one level of labels
		for (k=1;k<Nlabel;k++) {
			c = k;
			while (lb[c]!=c) c = lb[c];
			lb[k] = c;
		}
		// count the valid labels
		Nlb = 0;
		for (k=1;k<Nlabel;k++) {
			if (lb[k]==k) {
				// check the labels for 6-neighbors only 
				AddLabel = 0;
				if (lb[label[0][1][1]] == k) AddLabel=1;
				if (lb[label[1][0][1]] == k) AddLabel=1;
				if (lb[label[1][1][0]] == k) AddLabel=1;
				if (lb[label[2][1][1]] == k) AddLabel=1;
				if (lb[label[1][2][1]] == k) AddLabel=1;
				if (lb[label[1][1][2]] == k) AddLabel=1;
				if (AddLabel>0) Nlb++;
			}
		}
        // clean up
        label = null;
        lb = null;
        connect = null;

        return Nlb;
	}
	
	/** 
	 *  3D images: 6(18)-connectivity, positive  neighbors
	 */
	public static int positive618Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if (image[x+i-1][y+j-1][z+k-1] > image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
			
		return connected618Neighborhood3D(img);
	}

	/** 
	 *  3D images: 6(18)-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual618Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] >= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected618Neighborhood3D(img);
    }


	/** 
	 *  3D images: 6(18)-connectivity, negative neighbors
	 */
	public static int negative618Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if (image[x+i-1][y+j-1][z+k-1] < image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
    	return connected618Neighborhood3D(img);
    }

	/** 
	 *  3D images: 6(18)-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual618Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] <= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected618Neighborhood3D(img);
    }

    /** 
     *  3D images: 6-neighborhood  (implies 26-background) 
     */
    public static int connected626Neighborhood3D(boolean img[][][]) {
        int Nlabel = 0;
        int[][][]   label = new int[3][3][3];
        int[]       lb = new int[27];
        int lbMin;
        int i,j,k,l,c,n;
        int Nlb;
        int[]       connect = new int[6];
        int Nconnect;
        int AddLabel;
        
        for (i=0;i<3;i++)
            for (j=0;j<3;j++)
                for (k=0;k<3;k++)
                    label[i][j][k] = 0;
        
        lb[0] = 0;
        Nlabel = 1;
        // 6-neighbors only
        
        // -x
        if (img[0][1][1]) {
            // new, unconnected region
            label[0][1][1] = Nlabel;
            lb[Nlabel] = Nlabel;
            Nlabel++;
        }
        // -y
        if (img[1][0][1]) {
            // object point: neighbors ?
            Nconnect = 0;
            if ( (label[0][1][1] > 0) && (img[0][0][1]) ) {
                connect[Nconnect] = lb[ label[0][1][1] ];
                Nconnect++;
            }
            // if connected values, find the smallest lb label and attribute it
            // to all others (-> join labels)
            if (Nconnect>0) {
                lbMin = lb[connect[0]];
                for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                label[1][0][1] = lbMin;
            } else {
                // new, unconnected region
                label[1][0][1] = Nlabel;
                lb[Nlabel] = Nlabel;
                Nlabel++;
            }
        }
        // -z
        if (img[1][1][0]) {
            // object point: neighbors ?
            Nconnect = 0;
            if ( (label[0][1][1] > 0) && (img[0][1][0]) ) {
                connect[Nconnect] = lb[ label[0][1][1] ];
                Nconnect++;
            }
            if ( (label[1][0][1] > 0) && (img[1][0][0]) ) {
                connect[Nconnect] = lb[ label[1][0][1] ];
                Nconnect++;
            }
            // if connected values, find the smallest lb label and attribute it
            // to all others (-> join labels)
            if (Nconnect>0) {
                lbMin = lb[connect[0]];
                for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                label[1][1][0] = lbMin;
            } else {
                // new, unconnected region
                label[1][1][0] = Nlabel;
                lb[Nlabel] = Nlabel;
                Nlabel++;
            }
        }
        // +x
        if (img[2][1][1]) {
            // object point: neighbors ?
            Nconnect = 0;
            if ( (label[1][0][1] > 0) && (img[2][0][1]) ) {
                connect[Nconnect] = lb[ label[1][0][1] ];
                Nconnect++;
            }
            if ( (label[1][1][0] > 0) && (img[2][1][0]) ) {
                connect[Nconnect] = lb[ label[1][1][0] ];
                Nconnect++;
            }
            // if connected values, find the smallest lb label and attribute it
            // to all others (-> join labels)
            if (Nconnect>0) {
                lbMin = lb[connect[0]];
                for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                label[2][1][1] = lbMin;
            } else {
                // new, unconnected region
                label[2][1][1] = Nlabel;
                lb[Nlabel] = Nlabel;
                Nlabel++;
            }
        }
        // +y
        if (img[1][2][1]) {
            // object point: neighbors ?
            Nconnect = 0;
            if ( (label[0][1][1] > 0) && (img[0][2][1]) ) {
                connect[Nconnect] = lb[ label[0][1][1] ];
                Nconnect++;
            }
            if ( (label[1][1][0] > 0) && (img[1][2][0]) ) {
                connect[Nconnect] = lb[ label[1][1][0] ];
                Nconnect++;
            }
            if ( (label[2][1][1] > 0) && (img[2][2][1]) ) {
                connect[Nconnect] = lb[ label[2][1][1] ];
                Nconnect++;
            }
            // if connected values, find the smallest lb label and attribute it
            // to all others (-> join labels)
            if (Nconnect>0) {
                lbMin = lb[connect[0]];
                for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                label[1][2][1] = lbMin;
            } else {
                // new, unconnected region
                label[1][2][1] = Nlabel;
                lb[Nlabel] = Nlabel;
                Nlabel++;
            }
        }
        // +z
        if (img[1][1][2]) {
            // object point: neighbors ?
            Nconnect = 0;
            if ( (label[0][1][1] > 0) && (img[0][1][2]) ) {
                connect[Nconnect] = lb[ label[0][1][1] ];
                Nconnect++;
            }
            if ( (label[1][0][1] > 0) && (img[1][0][2]) ) {
                connect[Nconnect] = lb[ label[1][0][1] ];
                Nconnect++;
            }
            if ( (label[2][1][1] > 0) && (img[2][1][2]) ) {
                connect[Nconnect] = lb[ label[2][1][1] ];
                Nconnect++;
            }
            if ( (label[1][2][1] > 0) && (img[1][2][2]) ) {
                connect[Nconnect] = lb[ label[1][2][1] ];
                Nconnect++;
            }
            // if connected values, find the smallest lb label and attribute it
            // to all others (-> join labels)
            if (Nconnect>0) {
                lbMin = lb[connect[0]];
                for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                label[1][1][2] = lbMin;
            } else {
                // new, unconnected region
                label[1][1][2] = Nlabel;
                lb[Nlabel] = Nlabel;
                Nlabel++;
            }
        }
        // only one level of labels
        for (k=1;k<Nlabel;k++) {
            c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
        // count the valid labels
        Nlb = 0;
        for (k=1;k<Nlabel;k++) {
            if (lb[k]==k) {
                // check the labels for 6-neighbors only 
                AddLabel = 0;
                if (lb[label[0][1][1]] == k) AddLabel=1;
                if (lb[label[1][0][1]] == k) AddLabel=1;
                if (lb[label[1][1][0]] == k) AddLabel=1;
                if (lb[label[2][1][1]] == k) AddLabel=1;
                if (lb[label[1][2][1]] == k) AddLabel=1;
                if (lb[label[1][1][2]] == k) AddLabel=1;
                if (AddLabel>0) Nlb++;
            }
        }
        // clean up
        label = null;
        lb = null;
        connect = null;

        return Nlb;
    }

    /** 
	 *  3D images: 6(26)-connectivity, positive neighbors
	 */
	public static int positive626Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if (image[x+i-1][y+j-1][z+k-1] > image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
			
		return connected626Neighborhood3D(img);
	}

	/** 
	 *  3D images: 6(26)-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual626Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] >= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected626Neighborhood3D(img);
    }

	/** 
	 *  3D images: 6(26)-connectivity, negative neighbors
	 */
	public static int negative626Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if (image[x+i-1][y+j-1][z+k-1] < image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
    	return connected626Neighborhood3D(img);
    }

	/** 
	 *  3D images: 6(26)-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual626Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] <= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected626Neighborhood3D(img);
    }

    /**
     *  3D images: 18-neighborhood
     */
    public static int connected18Neighborhood3D(boolean img[][][]) {
        int Nlabel = 0;
        int[][][]   label = new int[3][3][3];
        int[]       lb = new int[27];
        int lbMin;
        int i,j,k,l,c;
        int Nlb;
        int[]   connect = new int[26];
        int Nconnect;
        
        for (i=0;i<3;i++)
            for (j=0;j<3;j++)
                for (k=0;k<3;k++)
                    label[i][j][k] = 0;
        
        lb[0] = 0;
        Nlabel = 1;
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                for (k=0;k<3;k++) {
                    if (img[i][j][k]) {
                        // object point: neighbors ?
                        Nconnect = 0;
						// 1
                        if (i>0) {
                            if (label[i-1][j][k] > 0) {
                                connect[Nconnect] = label[i-1][j][k];
                                Nconnect++;
                            }
							// 2
                            if (j>0) {
                                if (label[i-1][j-1][k] > 0) {
                                    connect[Nconnect] = lb[ label[i-1][j-1][k] ];
                                    Nconnect++;
                                }
                            }
							// 3
                            if (j<2) {
                                if (label[i-1][j+1][k] > 0) {
                                    connect[Nconnect] = lb[ label[i-1][j+1][k] ];
                                    Nconnect++;
                                }
                            }
							// 4
                            if (k>0) {
                                if (label[i-1][j][k-1] > 0) {
                                    connect[Nconnect] = lb[ label[i-1][j][k-1] ];
                                    Nconnect++;
                                }
                            }
							// 5
                            if (k<2) {
                                if (label[i-1][j][k+1] > 0) {
                                    connect[Nconnect] = lb[ label[i-1][j][k+1] ];
                                    Nconnect++;
                                }
                            }
                        }
						// 6
                        if (j>0) {
                            if (label[i][j-1][k] > 0) {
                                connect[Nconnect] = lb[ label[i][j-1][k] ];
                                Nconnect++;
                            }
							// 7
                            if (k>0) {
                                if (label[i][j-1][k-1] > 0) {
                                    connect[Nconnect] = lb[ label[i][j-1][k-1] ];
                                    Nconnect++;
                                }
                            }
							// 8
                            if (k<2) {
                                if (label[i][j-1][k+1] > 0) {
                                    connect[Nconnect] = lb[ label[i][j-1][k+1] ];
                                    Nconnect++;
                                }
                            }
                        }
						// 9
                        if (k>0) {
                            if (label[i][j][k-1] > 0) {
                                connect[Nconnect] = lb[ label[i][j][k-1] ];
                                Nconnect++;
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[i][j][k] = lbMin;
                        } else {
                            // new, unconnected region
                            label[i][j][k] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (k=1;k<Nlabel;k++) {
            c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
        // count the valid labels
        Nlb = 0;
        for (k=1;k<Nlabel;k++)
            if (lb[k]==k) Nlb++;
        //if (Nlb>0) printf("(%d,%d):%d\n",x,y,Nlb);
        // clean up
        label = null;
        lb = null;
        connect = null;
        
        return Nlb;
    }

    /** 
	 *  3D images: 18-connectivity, positive neighbors
	 */
	public static int positive18Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if (image[x+i-1][y+j-1][z+k-1] > image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
			
		return connected18Neighborhood3D(img);
	}

	/** 
	 *  3D images: 18-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual18Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] >= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected18Neighborhood3D(img);
    }

	/** 
	 *  3D images: 18-connectivity, negative neighbors
	 */
	public static int negative18Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if (image[x+i-1][y+j-1][z+k-1] < image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
    	return connected18Neighborhood3D(img);
    }

	/** 
	 *  3D images: 18-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual18Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] <= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected18Neighborhood3D(img);
    }
	
    /**
     *  3D images: 26-neighborhood
     */
    public static int connected26Neighborhood3D(boolean img[][][]) {
        int Nlabel = 0;
        int[][][]   label = new int[3][3][3];
        int[]       lb = new int[27];
        int lbMin;
        int i,j,k,x,y,z,l,c;
        int Nlb;
        int[]   connect = new int[26];
        int Nconnect;
        
        for (i=0;i<3;i++)
            for (j=0;j<3;j++)
                for (k=0;k<3;k++)
                    label[i][j][k] = 0;
        
        lb[0] = 0;
        Nlabel = 1;
		for (x=0;x<3;x++) {
			for (y=0;y<3;y++) {
				for (z=0;z<3;z++) {
					if (img[x][y][z]) {
                        // object point: neighbors ?
                        Nconnect = 0;
                        for (i=-1;i<=1;i++) for (j=-1;j<=1;j++) for (k=-1;k<=1;k++) {
                            if ( (x+i>=0) && (x+i<3) && (y+j>=0) && (y+j<3) && (z+k>=0) && (z+k<3) ) {
                                if (label[x+i][y+j][z+k] > 0) {
                                    connect[Nconnect] = lb[ label[x+i][y+j][z+k] ];
                                    Nconnect++;
                                }
                            }
                        }
                        // if connected values, find the smallest lb label and attribute it
                        // to all others (-> join labels)
                        if (Nconnect>0) {
                            //printf("c:%d",Nconnect);
                            lbMin = lb[connect[0]];
                            for (l=1;l<Nconnect;l++) lbMin = Math.min(lbMin,lb[connect[l]]);
                            for (l=0;l<Nconnect;l++) lb[connect[l]] = lbMin;
                            label[x][y][z] = lbMin;
                        } else {
                            // new, unconnected region
                            label[x][y][z] = Nlabel;
                            lb[Nlabel] = Nlabel;
                            //printf("l:%d", Nlabel);
                            Nlabel++;
                        }
                    }
                }
            }
        }
        // only one level of labels
        for (k=1;k<Nlabel;k++) {
            c = k;
            while (lb[c]!=c) c = lb[c];
            lb[k] = c;
        }
        // count the valid labels
        Nlb = 0;
        for (k=1;k<Nlabel;k++)
            if (lb[k]==k) Nlb++;
        //if (Nlb>0) printf("(%d,%d):%d\n",x,y,Nlb);
        // clean up
        label = null;
        lb = null;
        connect = null;
        
        return Nlb;
    }

    /** 
	 *  3D images: 26-connectivity, positive neighbors
	 */
	public static int positive26Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if (image[x+i-1][y+j-1][z+k-1] > image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
			
		return connected26Neighborhood3D(img);
	}

	/** 
	 *  3D images: 26-connectivity, positive or equal neighbors
	 */
	public static int positiveEqual26Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] >= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected26Neighborhood3D(img);
    }

	/** 
	 *  3D images: 26-connectivity, negative neighbors
	 */
	public static int negative26Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
            if (image[x+i-1][y+j-1][z+k-1] < image[x][y][z])
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
    	return connected26Neighborhood3D(img);
    }

	/** 
	 *  3D images: 26-connectivity, negative or equal neighbors
	 */
	public static int negativeEqual26Neighborhood3D(byte[][][] image, int x, int y, int z) {
		boolean img[][][] = new boolean[3][3][3];
	
		for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) {
			if ( ( (i!=1) || (j!=1) || (k!=1) ) && (image[x+i-1][y+j-1][z+k-1] <= image[x][y][z]) )
				img[i][j][k] = true;
			else
				img[i][j][k] = false;
		}
	
	    return connected26Neighborhood3D(img);
    }


}