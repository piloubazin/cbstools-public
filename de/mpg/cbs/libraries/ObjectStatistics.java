package de.mpg.cbs.libraries;

import java.io.*;
import java.util.*;

import de.mpg.cbs.utilities.*;

/**
 *
 *  This class computes labels and properties for sets of binary objects
 *	
 *  @version    July 2004
 *	@author     Pierre-Louis Bazin
 *		
 *
 */

public class ObjectStatistics {
	
	// no data: used as a library of functions
	private final static int MaxObject = 1000000;
	
	public final static int SUPERIOR =  10;
	public final static int SUPEQUAL =  11;
	public final static int INFERIOR =  12;
	public final static int INFEQUAL =  13;
	public final static int EQUAL =  	14;
	public final static int UNEQUAL =  	15;
	public final static int NONE =      25;
	public final static int AND =       26;
	public final static int OR =        27;
	public final static int XOR =       28;
	
	public final static float SQR2 = (float)Math.sqrt(2.0f);
	public final static float SQR3 = (float)Math.sqrt(3.0f);
   // object manipulation
    
	/**
     *  Counts the different boundaries
	 *	and returns the Euler characteristic
     */
    public static final int eulerCharacteristic(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int Nv = 0;
		int Ne = 0;
		int Nf = 0;
		
		// hyp: the object is away from the image boundary
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// count all boundaries : 6-C
			if (img[x][y][z]!=img[x-1][y][z]) Nf++;
			if (img[x][y][z]!=img[x][y-1][z]) Nf++;
			if (img[x][y][z]!=img[x][y][z-1]) Nf++;
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// edges : XY
			if (img[x][y][z]!=img[x-1][y-1][z]) Ne++;
			else if (img[x-1][y][z]!=img[x][y-1][z]) Ne++;
			else if ( (img[x][y][z]==img[x-1][y-1][z]) && (img[x][y][z]!=img[x][y-1][z]) && (img[x][y][z]!=img[x-1][y][z]) ) {
				Ne+=2;
			} else if ( (img[x-1][y][z]==img[x][y-1][z]) && (img[x-1][y][z]!=img[x][y][z]) && (img[x-1][y][z]!=img[x-1][y-1][z]) ) {
				Ne+=2;
			}
			// edges : YZ
			if (img[x][y][z]!=img[x][y-1][z-1]) Ne++;
			else if (img[x][y-1][z]!=img[x][y][z-1]) Ne++;
			else if ( (img[x][y][z]==img[x][y-1][z-1]) && (img[x][y][z]!=img[x][y-1][z]) && (img[x][y][z]!=img[x][y][z-1]) ) {
				Ne+=2;
			} else if ( (img[x][y-1][z]==img[x][y][z-1]) && (img[x][y-1][z]!=img[x][y][z]) && (img[x][y-1][z]!=img[x][y-1][z-1]) ) {
				Ne+=2;
			}
			// edges : ZX
			if (img[x][y][z]!=img[x-1][y][z-1]) Ne++;
			else if (img[x][y][z-1]!=img[x-1][y][z]) Ne++;
			else if ( (img[x][y][z]==img[x-1][y][z-1]) && (img[x][y][z]!=img[x][y][z-1]) && (img[x][y][z]!=img[x-1][y][z]) ) {
				Ne+=2;
			} else if ( (img[x][y][z-1]==img[x-1][y][z]) && (img[x][y][z-1]!=img[x][y][z]) && (img[x][y][z-1]!=img[x-1][y][z-1]) ) {
				Ne+=2;
			}
		}
		int Nobj,Nbg,N6,N18,N26;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// classify by cases : look into 2x2 neighborhood, count elements and connections
			Nobj=0;Nbg=0;N6=0;N18=0;N26=0;
			if (img[x][y][z]) Nobj++; else Nbg++;
			if (img[x-1][y][z]) Nobj++; else Nbg++;
			if (img[x][y-1][z]) Nobj++; else Nbg++;
			if (img[x][y][z-1]) Nobj++; else Nbg++;
			if (img[x][y-1][z-1]) Nobj++; else Nbg++;
			if (img[x-1][y][z-1]) Nobj++; else Nbg++;
			if (img[x-1][y-1][z]) Nobj++; else Nbg++;
			if (img[x-1][y-1][z-1]) Nobj++; else Nbg++;
			
			// use the smallest one
			if (Nobj<=Nbg) {
				if (Nobj==0) { 
					// do nothing 
				} else if (Nobj==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (img[x][y][z]) && (img[x-1][y][z]) ) N6++;
					if ( (img[x][y][z]) && (img[x][y-1][z]) ) N6++;
					if ( (img[x][y][z]) && (img[x][y][z-1]) ) N6++;
					if ( (img[x-1][y][z]) && (img[x-1][y-1][z]) ) N6++;
					if ( (img[x-1][y][z]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x][y-1][z]) && (img[x-1][y-1][z]) ) N6++;
					if ( (img[x][y-1][z]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x][y][z-1]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x][y][z-1]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x][y-1][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x-1][y][z-1]) ) N6++;
					if ( (img[x-1][y-1][z-1]) && (img[x-1][y-1][z]) ) N6++;
					//System.out.print("case: obj "+Nobj+", ("+N6+" |");
						
					if ( (Nobj==2) && (N6==1) ) Nv++;
					else if ( (Nobj==3) && (N6==2) ) Nv++;
					else if ( (Nobj==4) && (N6==4) ) Nv++;
					else if ( (Nobj==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (img[x][y][z]) && (img[x][y-1][z-1]) ) N18++;
						if ( (img[x][y][z]) && (img[x-1][y][z-1]) ) N18++;
						if ( (img[x][y][z]) && (img[x-1][y-1][z]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x][y-1][z]) ) N18++;
						if ( (img[x-1][y][z]) && (img[x][y][z-1]) ) N18++;
						if ( (img[x][y-1][z]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x][y-1][z]) && (img[x][y][z-1]) ) N18++;
						if ( (img[x][y][z-1]) && (img[x-1][y-1][z-1]) ) N18++;
						if ( (img[x][y-1][z-1]) && (img[x-1][y][z-1]) ) N18++;
						if ( (img[x-1][y][z-1]) && (img[x-1][y-1][z]) ) N18++;
						if ( (img[x-1][y-1][z]) && (img[x][y-1][z-1]) ) N18++;
						//System.out.println(" "+N18+")");
						
						if ( (Nobj==2) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==3) ) {
							if (cObj==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nobj==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nobj==2) {
								if (cObj==26) Nv+=0;
								else Nv+=2;
							} else {
								//System.out.println("!:"+Nobj+", "+N6+", "+N18);
							}
						}
					}
				}
			} else {
				if (Nbg==0) { 
					// do nothing 
				} else if (Nbg==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (!img[x][y][z]) && (!img[x-1][y][z]) ) N6++;
					if ( (!img[x][y][z]) && (!img[x][y-1][z]) ) N6++;
					if ( (!img[x][y][z]) && (!img[x][y][z-1]) ) N6++;
					if ( (!img[x-1][y][z]) && (!img[x-1][y-1][z]) ) N6++;
					if ( (!img[x-1][y][z]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x][y-1][z]) && (!img[x-1][y-1][z]) ) N6++;
					if ( (!img[x][y-1][z]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x][y][z-1]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x][y][z-1]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x][y-1][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x-1][y][z-1]) ) N6++;
					if ( (!img[x-1][y-1][z-1]) && (!img[x-1][y-1][z]) ) N6++;
					
					if ( (Nbg==2) && (N6==1) ) Nv++;
					else if ( (Nbg==3) && (N6==2) ) Nv++;
					else if ( (Nbg==4) && (N6==4) ) Nv++;
					else if ( (Nbg==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (!img[x][y][z]) && (!img[x][y-1][z-1]) ) N18++;
						if ( (!img[x][y][z]) && (!img[x-1][y][z-1]) ) N18++;
						if ( (!img[x][y][z]) && (!img[x-1][y-1][z]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x][y-1][z]) ) N18++;
						if ( (!img[x-1][y][z]) && (!img[x][y][z-1]) ) N18++;
						if ( (!img[x][y-1][z]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x][y-1][z]) && (!img[x][y][z-1]) ) N18++;
						if ( (!img[x][y][z-1]) && (!img[x-1][y-1][z-1]) ) N18++;
						if ( (!img[x][y-1][z-1]) && (!img[x-1][y][z-1]) ) N18++;
						if ( (!img[x-1][y][z-1]) && (!img[x-1][y-1][z]) ) N18++;
						if ( (!img[x-1][y-1][z]) && (!img[x][y-1][z-1]) ) N18++;
						
						if ( (Nbg==2) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==3) ) {
							if (cBg==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nbg==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nbg==2) {
								if (cBg==26) Nv+=0;
								else Nv+=2;
							} else {
								//System.out.println("!:"+Nbg+", "+N6+", "+N18);
							}
						}
					}
				}
			}

		}
		
		//System.out.println("Surface: "+Nv+" vertices, "+Ne+" edges, "+Nf+" faces");
			
		return Nv-Ne+Nf;
	}

	/**
     *  Counts the different boundaries
	 *	and returns the Euler characteristic
     */
    public static final int eulerCharacteristic(boolean img[], int nx, int ny, int nz, int cObj, int cBg) {
		int Nv = 0;
		int Ne = 0;
		int Nf = 0;
		
		// hyp: the object is away from the image boundary
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// count all boundaries : 6-C
			if (img[x+nx*y+nx*ny*z]!=img[(x-1)+nx*y+nx*ny*z]) Nf++;
			if (img[x+nx*y+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*z]) Nf++;
			if (img[x+nx*y+nx*ny*z]!=img[x+nx*y+nx*ny*(z-1)]) Nf++;
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// edges : XY
			if (img[x+nx*y+nx*ny*z]!=img[(x-1)+nx*(y-1)+nx*ny*z]) Ne++;
			else if (img[(x-1)+nx*y+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*z]) Ne++;
			else if ( (img[x+nx*y+nx*ny*z]==img[(x-1)+nx*(y-1)+nx*ny*z]) && (img[x+nx*y+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*z]) && (img[x+nx*y+nx*ny*z]!=img[(x-1)+nx*y+nx*ny*z]) ) {
				Ne+=2;
			} else if ( (img[(x-1)+nx*y+nx*ny*z]==img[x+nx*(y-1)+nx*ny*z]) && (img[(x-1)+nx*y+nx*ny*z]!=img[x+nx*y+nx*ny*z]) && (img[(x-1)+nx*y+nx*ny*z]!=img[(x-1)+nx*(y-1)+nx*ny*z]) ) {
				Ne+=2;
			}
			// edges : YZ
			if (img[x+nx*y+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*(z-1)]) Ne++;
			else if (img[x+nx*(y-1)+nx*ny*z]!=img[x+nx*y+nx*ny*(z-1)]) Ne++;
			else if ( (img[x+nx*y+nx*ny*z]==img[x+nx*(y-1)+nx*ny*(z-1)]) && (img[x+nx*y+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*z]) && (img[x+nx*y+nx*ny*z]!=img[x+nx*y+nx*ny*(z-1)]) ) {
				Ne+=2;
			} else if ( (img[x+nx*(y-1)+nx*ny*z]==img[x+nx*y+nx*ny*(z-1)]) && (img[x+nx*(y-1)+nx*ny*z]!=img[x+nx*y+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*z]!=img[x+nx*(y-1)+nx*ny*(z-1)]) ) {
				Ne+=2;
			}
			// edges : ZX
			if (img[x+nx*y+nx*ny*z]!=img[(x-1)+nx*y+nx*ny*(z-1)]) Ne++;
			else if (img[x+nx*y+nx*ny*(z-1)]!=img[(x-1)+nx*y+nx*ny*z]) Ne++;
			else if ( (img[x+nx*y+nx*ny*z]==img[(x-1)+nx*y+nx*ny*(z-1)]) && (img[x+nx*y+nx*ny*z]!=img[x+nx*y+nx*ny*(z-1)]) && (img[x+nx*y+nx*ny*z]!=img[(x-1)+nx*y+nx*ny*z]) ) {
				Ne+=2;
			} else if ( (img[x+nx*y+nx*ny*(z-1)]==img[(x-1)+nx*y+nx*ny*z]) && (img[x+nx*y+nx*ny*(z-1)]!=img[x+nx*y+nx*ny*z]) && (img[x+nx*y+nx*ny*(z-1)]!=img[(x-1)+nx*y+nx*ny*(z-1)]) ) {
				Ne+=2;
			}
		}
		int Nobj,Nbg,N6,N18,N26;
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			// classify by cases : look into 2x2 neighborhood, count elements and connections
			Nobj=0;Nbg=0;N6=0;N18=0;N26=0;
			if (img[x+nx*y+nx*ny*z]) Nobj++; else Nbg++;
			if (img[(x-1)+nx*y+nx*ny*z]) Nobj++; else Nbg++;
			if (img[x+nx*(y-1)+nx*ny*z]) Nobj++; else Nbg++;
			if (img[x+nx*y+nx*ny*(z-1)]) Nobj++; else Nbg++;
			if (img[x+nx*(y-1)+nx*ny*(z-1)]) Nobj++; else Nbg++;
			if (img[(x-1)+nx*y+nx*ny*(z-1)]) Nobj++; else Nbg++;
			if (img[(x-1)+nx*(y-1)+nx*ny*z]) Nobj++; else Nbg++;
			if (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) Nobj++; else Nbg++;
			
			// use the smallest one
			if (Nobj<=Nbg) {
				if (Nobj==0) { 
					// do nothing 
				} else if (Nobj==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (img[x+nx*y+nx*ny*z]) && (img[(x-1)+nx*y+nx*ny*z]) ) N6++;
					if ( (img[x+nx*y+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (img[x+nx*y+nx*ny*z]) && (img[x+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (img[(x-1)+nx*y+nx*ny*z]) && (img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (img[(x-1)+nx*y+nx*ny*z]) && (img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (img[x+nx*(y-1)+nx*ny*z]) && (img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (img[x+nx*(y-1)+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (img[x+nx*y+nx*ny*(z-1)]) && (img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (img[x+nx*y+nx*ny*(z-1)]) && (img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					//System.out.print("case: obj "+Nobj+", ("+N6+" |");
						
					if ( (Nobj==2) && (N6==1) ) Nv++;
					else if ( (Nobj==3) && (N6==2) ) Nv++;
					else if ( (Nobj==4) && (N6==4) ) Nv++;
					else if ( (Nobj==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (img[x+nx*y+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*y+nx*ny*z]) && (img[(x-1)+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*y+nx*ny*z]) && (img[(x-1)+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (img[(x-1)+nx*y+nx*ny*z]) && (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (img[(x-1)+nx*y+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (img[(x-1)+nx*y+nx*ny*z]) && (img[x+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*(y-1)+nx*ny*z]) && (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*(y-1)+nx*ny*z]) && (img[x+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*y+nx*ny*(z-1)]) && (img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (img[x+nx*(y-1)+nx*ny*(z-1)]) && (img[(x-1)+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (img[(x-1)+nx*y+nx*ny*(z-1)]) && (img[(x-1)+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (img[(x-1)+nx*(y-1)+nx*ny*z]) && (img[x+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						//System.out.println(" "+N18+")");
						
						if ( (Nobj==2) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==1) ) {
							if (cObj==6) Nv+=2;
							else Nv++;
						} else if ( (Nobj==3) && (N18==3) ) {
							if (cObj==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nobj==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nobj==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nobj==2) {
								if (cObj==26) Nv+=0;
								else Nv+=2;
							} else {
								//System.out.println("!:"+Nobj+", "+N6+", "+N18);
							}
						}
					}
				}
			} else {
				if (Nbg==0) { 
					// do nothing 
				} else if (Nbg==1) {
					Nv++;
				} else {
					//count 6-C connections
					if ( (!img[x+nx*y+nx*ny*z]) && (!img[(x-1)+nx*y+nx*ny*z]) ) N6++;
					if ( (!img[x+nx*y+nx*ny*z]) && (!img[x+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (!img[x+nx*y+nx*ny*z]) && (!img[x+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (!img[(x-1)+nx*y+nx*ny*z]) && (!img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (!img[(x-1)+nx*y+nx*ny*z]) && (!img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (!img[x+nx*(y-1)+nx*ny*z]) && (!img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					if ( (!img[x+nx*(y-1)+nx*ny*z]) && (!img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (!img[x+nx*y+nx*ny*(z-1)]) && (!img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (!img[x+nx*y+nx*ny*(z-1)]) && (!img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (!img[x+nx*(y-1)+nx*ny*(z-1)]) ) N6++;
					if ( (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (!img[(x-1)+nx*y+nx*ny*(z-1)]) ) N6++;
					if ( (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) && (!img[(x-1)+nx*(y-1)+nx*ny*z]) ) N6++;
					
					if ( (Nbg==2) && (N6==1) ) Nv++;
					else if ( (Nbg==3) && (N6==2) ) Nv++;
					else if ( (Nbg==4) && (N6==4) ) Nv++;
					else if ( (Nbg==4) && (N6==3) ) Nv++;
					else {
						// count 18-C connections
						if ( (!img[x+nx*y+nx*ny*z]) && (!img[x+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*y+nx*ny*z]) && (!img[(x-1)+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*y+nx*ny*z]) && (!img[(x-1)+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (!img[(x-1)+nx*y+nx*ny*z]) && (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (!img[(x-1)+nx*y+nx*ny*z]) && (!img[x+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (!img[(x-1)+nx*y+nx*ny*z]) && (!img[x+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*(y-1)+nx*ny*z]) && (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*(y-1)+nx*ny*z]) && (!img[x+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*y+nx*ny*(z-1)]) && (!img[(x-1)+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						if ( (!img[x+nx*(y-1)+nx*ny*(z-1)]) && (!img[(x-1)+nx*y+nx*ny*(z-1)]) ) N18++;
						if ( (!img[(x-1)+nx*y+nx*ny*(z-1)]) && (!img[(x-1)+nx*(y-1)+nx*ny*z]) ) N18++;
						if ( (!img[(x-1)+nx*(y-1)+nx*ny*z]) && (!img[x+nx*(y-1)+nx*ny*(z-1)]) ) N18++;
						
						if ( (Nbg==2) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==1) ) {
							if (cBg==6) Nv+=2;
							else Nv++;
						} else if ( (Nbg==3) && (N18==3) ) {
							if (cBg==6) Nv+=3;
							else Nv+=2;
						} else if ( (Nbg==4) && (N18==2) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==3) ) {
							Nv+=2;
						} else if ( (Nbg==4) && (N18==6) ) {
							Nv+=4;
						} else {
							// no need to count the 26-connections
							if (Nbg==2) {
								if (cBg==26) Nv+=0;
								else Nv+=2;
							} else {
								//System.out.println("!:"+Nbg+", "+N6+", "+N18);
							}
						}
					}
				}
			}

		}
		
		//System.out.println("Surface: "+Nv+" vertices, "+Ne+" edges, "+Nf+" faces");
			
		return Nv-Ne+Nf;
	}

	/**
     *  checks that the object has all points being 6 connected to each neighbors
	 *	(not properly functional yet)
     */
    public static final float[][][] fullConnection(boolean img[][][], int nx, int ny, int nz) {
		float [][][] res = new float[nx][ny][nz];
		int[] count = new int[8];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			res[x][y][z] = 0;
			if (img[x][y][z]) {
				res[x][y][z]++;
				for (int n=0;n<8;n++) count[n] = 0;
				// count the number of 8-voxel cubes
				for (int i=0;i<=1;i++) for (int j=0;j<=1;j++) for (int l=0;l<=1;l++) {
					if (img[x+i][y+j][z+l]) count[0]++;
					if (img[x-i][y+j][z+l]) count[1]++;
					if (img[x+i][y-j][z+l]) count[2]++;
					if (img[x+i][y+j][z-l]) count[3]++;
					if (img[x+i][y-j][z-l]) count[4]++;
					if (img[x-i][y+j][z-l]) count[5]++;
					if (img[x-i][y-j][z+l]) count[6]++;
					if (img[x-i][y-j][z-l]) count[7]++;
				}
				for (int n=0;n<8;n++) if (count[n]==8) res[x][y][z]++;
			}
		}
		return res;
	}
	
	/**
     *  compute the highest connectivity junction for each pixel.
	 *	<p>
	 *	The result is 6 (only 6-C neighbors), 18 (only 6 and 18-C neighbors)
	 *	or 26 (6, 18 and 26-C neighbors).
     */
    public static final float[][][] connectivity(boolean img[][][], int nx, int ny, int nz) {
		float [][][] obj = new float[nx][ny][nz];
		int[] Nb = new int[8];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]) {
				obj[x][y][z] = 6;
				if ( img[x+1][y+1][z] && !(img[x+1][y][z] || img[x][y+1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y+1][z+1] && !(img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 18;
				if ( img[x+1][y][z+1] && !(img[x][y][z+1] || img[x+1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x-1][y+1][z] && !(img[x-1][y][z] || img[x][y+1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y-1][z+1] && !(img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 18;
				if ( img[x+1][y][z-1] && !(img[x][y][z-1] || img[x+1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x+1][y-1][z] && !(img[x+1][y][z] || img[x][y-1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y+1][z-1] && !(img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 18;
				if ( img[x-1][y][z+1] && !(img[x][y][z+1] || img[x-1][y][z]) ) obj[x][y][z] = 18;
				
				if ( img[x-1][y-1][z] && !(img[x-1][y][z] || img[x][y-1][z]) ) obj[x][y][z] = 18;
				if ( img[x][y-1][z-1] && !(img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 18;
				if ( img[x-1][y][z-1] && !(img[x][y][z-1] || img[x-1][y][z]) ) obj[x][y][z] = 18;
				
				// note: for 26-C, test only the case with no neighbors 
				// (otherwise, one of the neighbors show up as 18-C)
				if ( img[x+1][y+1][z+1] && !(img[x+1][y][z] || img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y+1][z+1] && !(img[x-1][y][z] || img[x][y+1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y-1][z+1] && !(img[x+1][y][z] || img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y+1][z-1] && !(img[x+1][y][z] || img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x+1][y-1][z-1] && !(img[x+1][y][z] || img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y+1][z-1] && !(img[x-1][y][z] || img[x][y+1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y-1][z+1] && !(img[x-1][y][z] || img[x][y-1][z] || img[x][y][z+1]) ) obj[x][y][z] = 26;
				if ( img[x-1][y-1][z-1] && !(img[x-1][y][z] || img[x][y-1][z] || img[x][y][z-1]) ) obj[x][y][z] = 26;

			} else {
				obj[x][y][z] = 0;
			}
		}
		return obj;
	}
	
	
	/**
	 *  find well-composed and not well-composed regions: groups imgects with relations
	 */
    public static final float[][][] wellComposed(boolean[][][] img, int nx, int ny, int nz) {
		float[][][] obj = new float[nx][ny][nz];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) obj[x][y][z] = 0.0f;
			else if (isWellComposed(img,x,y,z)) obj[x][y][z] = 1.0f;
			else obj[x][y][z] = 2.0f;
		}
		return obj;
	}
				
				
	public static final boolean isWellComposed(boolean[][][] img, int x, int y, int z) {
		if (img[x][y][z]) {
			// 18-C
			if (img[x-1][y-1][z] && !(img[x-1][y][z] || img[x][y-1][z]) ) return false;
			if (img[x][y-1][z-1] && !(img[x][y-1][z] || img[x][y][z-1]) ) return false;
			if (img[x-1][y][z-1] && !(img[x][y][z-1] || img[x-1][y][z]) ) return false;
		
			if (img[x-1][y+1][z] && !(img[x-1][y][z] || img[x][y+1][z]) ) return false;
			if (img[x][y-1][z+1] && !(img[x][y-1][z] || img[x][y][z+1]) ) return false;
			if (img[x+1][y][z-1] && !(img[x][y][z-1] || img[x+1][y][z]) ) return false;
		
			if (img[x+1][y-1][z] && !(img[x+1][y][z] || img[x][y-1][z]) ) return false;
			if (img[x][y+1][z-1] && !(img[x][y+1][z] || img[x][y][z-1]) ) return false;
			if (img[x-1][y][z+1] && !(img[x][y][z+1] || img[x-1][y][z]) ) return false;
				
			if (img[x+1][y+1][z] && !(img[x+1][y][z] || img[x][y+1][z]) ) return false;
			if (img[x][y+1][z+1] && !(img[x][y+1][z] || img[x][y][z+1]) ) return false;
			if (img[x+1][y][z+1] && !(img[x][y][z+1] || img[x+1][y][z]) ) return false;
		
			// 26-C
			if (img[x-1][y-1][z-1] && !( (img[x-1][y][z] && img[x-1][y-1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z-1]) 
									  || (img[x][y-1][z] && img[x][y-1][z-1]) 
									  || (img[x][y-1][z] && img[x-1][y-1][z]) 
									  || (img[x][y][z-1] && img[x-1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y-1][z-1]) ) ) return false;
		
			if (img[x-1][y-1][z+1] && !( (img[x-1][y][z] && img[x-1][y-1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z+1]) 
									  || (img[x][y-1][z] && img[x][y-1][z+1]) 
									  || (img[x][y-1][z] && img[x-1][y-1][z]) 
									  || (img[x][y][z+1] && img[x-1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y-1][z+1]) ) ) return false;
		
			if (img[x-1][y+1][z-1] && !( (img[x-1][y][z] && img[x-1][y+1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z-1]) 
									  || (img[x][y+1][z] && img[x][y+1][z-1]) 
									  || (img[x][y+1][z] && img[x-1][y+1][z]) 
									  || (img[x][y][z-1] && img[x-1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y+1][z-1]) ) ) return false;
		
			if (img[x+1][y-1][z-1] && !( (img[x+1][y][z] && img[x+1][y-1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z-1]) 
									  || (img[x][y-1][z] && img[x][y-1][z-1]) 
									  || (img[x][y-1][z] && img[x+1][y-1][z]) 
									  || (img[x][y][z-1] && img[x+1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y-1][z-1]) ) ) return false;
		
			if (img[x-1][y+1][z+1] && !( (img[x-1][y][z] && img[x-1][y+1][z]) 
									  || (img[x-1][y][z] && img[x-1][y][z+1]) 
									  || (img[x][y+1][z] && img[x][y+1][z+1]) 
									  || (img[x][y+1][z] && img[x-1][y+1][z]) 
									  || (img[x][y][z+1] && img[x-1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y+1][z+1]) ) ) return false;
		
			if (img[x+1][y-1][z+1] && !( (img[x+1][y][z] && img[x+1][y-1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z+1]) 
									  || (img[x][y-1][z] && img[x][y-1][z+1]) 
									  || (img[x][y-1][z] && img[x+1][y-1][z]) 
									  || (img[x][y][z+1] && img[x+1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y-1][z+1]) ) ) return false;
		
			if (img[x+1][y+1][z-1] && !( (img[x+1][y][z] && img[x+1][y+1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z-1]) 
									  || (img[x][y+1][z] && img[x][y+1][z-1]) 
									  || (img[x][y+1][z] && img[x+1][y+1][z]) 
									  || (img[x][y][z-1] && img[x+1][y][z-1]) 
									  || (img[x][y][z-1] && img[x][y+1][z-1]) ) ) return false;
		
			if (img[x+1][y+1][z+1] && !( (img[x+1][y][z] && img[x+1][y+1][z]) 
									  || (img[x+1][y][z] && img[x+1][y][z+1]) 
									  || (img[x][y+1][z] && img[x][y+1][z+1]) 
									  || (img[x][y+1][z] && img[x+1][y+1][z]) 
									  || (img[x][y][z+1] && img[x+1][y][z+1]) 
									  || (img[x][y][z+1] && img[x][y+1][z+1]) ) ) return false;
		}
		
		return true;
	}
	
	public static final boolean isWellComposed(boolean[] img, int xyz, int dx, int dy, int dz) {
		if (img[xyz]) {
			// 18-C
			if (img[xyz-dx-dy] 		&& !(img[xyz-dx] 	|| img[xyz-dy]) ) 	return false;
			if (img[xyz-dy-dz] 		&& !(img[xyz-dy] 	|| img[xyz-dz]) ) 	return false;
			if (img[xyz-dz-dx] 		&& !(img[xyz-dz] 	|| img[xyz-dx]) ) 	return false;
		
			if (img[xyz-dx+dy] 		&& !(img[xyz-dx] 	|| img[xyz+dy]) ) 	return false;
			if (img[xyz-dy+dz] 		&& !(img[xyz-dy] 	|| img[xyz+dz]) ) 	return false;
			if (img[xyz-dz+dx] 		&& !(img[xyz-dz] 	|| img[xyz+dx]) ) 	return false;
		
			if (img[xyz+dx-dy] 		&& !(img[xyz+dx] 	|| img[xyz-dy]) ) 	return false;
			if (img[xyz+dy-dz] 		&& !(img[xyz+dy] 	|| img[xyz-dz]) ) 	return false;
			if (img[xyz+dz-dx] 		&& !(img[xyz+dz] 	|| img[xyz-dx]) ) 	return false;
				
			if (img[xyz+dx+dy] 		&& !(img[xyz+dx] 	|| img[xyz+dy]) ) 	return false;
			if (img[xyz+dy+dz] 		&& !(img[xyz+dy] 	|| img[xyz+dz]) ) 	return false;
			if (img[xyz+dz+dx] 		&& !(img[xyz+dz] 	|| img[xyz+dx]) ) 	return false;
		
			// 26-C
			if (img[xyz-dx-dy-dz] && !( (img[xyz-dx] 	&& img[xyz-dx-dy]) 
									  || (img[xyz-dx] 	&& img[xyz-dx-dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy-dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy-dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz-dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz-dy]) ) ) return false;
		
			if (img[xyz-dx-dy+dz] && !( (img[xyz-dx] 	&& img[xyz-dx-dy]) 
									  || (img[xyz-dx] 	&& img[xyz-dx+dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy+dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy-dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz-dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz-dy]) ) ) return false;
		
			if (img[xyz-dx+dy-dz] && !( (img[xyz-dx] 	&& img[xyz-dx+dy]) 
									  || (img[xyz-dx] 	&& img[xyz-dx-dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy-dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy-dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz-dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz+dy]) ) ) return false;
		
			if (img[xyz+dx-dy-dz] && !( (img[xyz+dx] 	&& img[xyz+dx-dy]) 
									  || (img[xyz+dx] 	&& img[xyz+dx-dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy-dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy+dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz+dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz-dy]) ) ) return false;
		
			if (img[xyz-dx+dy+dz] && !( (img[xyz-dx] 	&& img[xyz-dx+dy]) 
									  || (img[xyz-dx] 	&& img[xyz-dx+dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy+dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy-dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz-dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz+dy]) ) ) return false;
		
			if (img[xyz+dx-dy+dz] && !( (img[xyz+dx] 	&& img[xyz+dx-dy]) 
									  || (img[xyz+dx] 	&& img[xyz+dx+dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy+dz]) 
									  || (img[xyz-dy] 	&& img[xyz-dy+dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz+dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz-dy]) ) ) return false;
		
			if (img[xyz+dx+dy-dz] && !( (img[xyz+dx] 	&& img[xyz+dx+dy]) 
									  || (img[xyz+dx] 	&& img[xyz+dx-dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy-dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy+dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz+dx]) 
									  || (img[xyz-dz] 	&& img[xyz-dz+dy]) ) ) return false;

			if (img[xyz+dx+dy+dz] && !( (img[xyz+dx] 	&& img[xyz+dx+dy]) 
									  || (img[xyz+dx] 	&& img[xyz+dx+dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy+dz]) 
									  || (img[xyz+dy] 	&& img[xyz+dy+dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz+dx]) 
									  || (img[xyz+dz] 	&& img[xyz+dz+dy]) ) ) return false;
		
		}
		
		return true;
	}
	
	/**
	 *  compute a simplistic mean curvature function on the object boundary
	 */
	public static final float[][][] objectCurvature(boolean[][][] obj, int nx, int ny, int nz) {
		float   	num, den;
		float		u,xp,xm,yp,ym,zp,zm;
		float		xpyp,xpym,xpzp,xpzm,xmyp,xmym,xmzp,xmzm,ypzp,ypzm,ymzp,ymzm;
		float		xpypzp,xpypzm,xpymzp,xpymzm,xmypzp,xmypzm,xmymzp,xmymzm;
		float   	ux,uy,uz,uxx,uyy,uzz,uxy,uyz,uzx;
		float		val;
		float[][][]		curv = new float[nx][ny][nz];
		float[][][]		level = new float[nx][ny][nz];
		int 		boundary;
		
		// create a 'levelset' map
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (obj[x][y][z]) {
				// check for boundary
				boundary = 0;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if ( (i*i+j*j+l*l==1) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>6) ) boundary = 6;
					else if ( (i*i+j*j+l*l==2) && (!obj[x+i][y+j][z+l]) && (boundary==0 || boundary>18) ) boundary = 18;
					else if ( (i*i+j*j+l*l==3) && (!obj[x+i][y+j][z+l]) && (boundary==0) ) boundary = 26;
				}
				// attribute a levelset for the different boundaries
				if (boundary==0) level[x][y][z] = -1.0f;
				else level[x][y][z] = 0.0f;
				/*
				else if (boundary==6) level[x][y][z] = 0.0f;
				else if (boundary==18) level[x][y][z] = -SQR2+1.0f;
				else if (boundary==26) level[x][y][z] = -SQR3+1.0f;
				*/
			} else {
				level[x][y][z] = 1.0f;
			}
		}
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (level[x][y][z]==0) {
				// set the level set values
				u = level[x][y][z];
				
				// 6 - neighbors		
				xp = level[x+1][y][z];
				xm = level[x-1][y][z];
				yp = level[x][y+1][z];
				ym = level[x][y+1][z];
				zp = level[x][y][z+1];
				zm = level[x][y][z-1];
				
				// 18-neighbors		
				xpyp = level[x+1][y+1][z];
				xmyp = level[x-1][y+1][z];
				xpym = level[x+1][y-1][z];
				xmym = level[x-1][y-1][z];
				xpzp = level[x+1][y][z+1];
				xmzp = level[x-1][y][z+1];
				xpzm = level[x+1][y][z-1];
				xmzm = level[x-1][y][z-1];
				ypzp = level[x][y+1][z+1];
				ymzp = level[x][y-1][z+1];
				ypzm = level[x][y+1][z-1];
				ymzm = level[x][y-1][z-1];
				
				// 26-neighbors
				xpypzm = level[x+1][y+1][z-1];
				xpymzm = level[x+1][y-1][z-1];
				xmypzm = level[x-1][y+1][z-1];
				xmymzm = level[x-1][y-1][z-1];
				xpypzp = level[x+1][y+1][z+1];
				xpymzp = level[x+1][y-1][z+1];
				xmypzp = level[x-1][y+1][z+1];
				xmymzp = level[x-1][y-1][z+1];
			
				// central differences ?
				ux = 0.25f*( xp + xpyp + xpzp + xpypzp 
							-xm - xmyp - xmzp - xmypzp );
				
				uy = 0.25f*( yp + xpyp + ypzp + xpypzp 
							-ym - xpym - ymzp - xpymzp );
				
				uz = 0.25f*( zp + xpzp + ypzp + xpypzp 
							-zm - xpzm - ypzm - xpypzm );
		
				uxx = 0.0625f*( 4*xp + 2*xpyp + 2*xpzp + 2*xpym + 2*xpzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*xm + 2*xmyp + 2*xmzp + 2*xmym + 2*xmzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*yp + 2*zp + 2*ym + 2*zm + ypzp + ymzp + ypzm + ymzm);
							 
				uyy = 0.0625f*( 4*yp + 2*xpyp + 2*ypzp + 2*xmyp + 2*ypzm + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*ym + 2*xpym + 2*ymzp + 2*xmym + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*zp + 2*xm + 2*zm + xpzp + xmzp + xpzm + xmzm);
							 
				uzz = 0.0625f*( 4*zp + 2*xpzp + 2*ypzp + 2*xmzp + 2*ymzp + xpypzp + xpymzp + xpypzm + xpymzm
							   +4*zm + 2*xpzm + 2*ypzm + 2*xmzm + 2*ymzm + xmypzp + xmymzp + xmypzm + xmymzm )
					 -0.125f*( 4*u + 2*xp + 2*yp + 2*xm + 2*ym + xpyp + xmyp + xpym + xmym);
							 
		
				uxy = 0.0625f*( 2*xpyp + 2*xmym + xpypzp + xpypzm + xmymzp + xmymzm
							  - 2*xpym - 2*xmyp - xpymzp - xpymzm - xmypzp - xmypzm );
							 
				uyz = 0.0625f*( 2*ypzp + 2*ymzm + xpypzp + xmypzp + xpymzm + xmymzm
							  - 2*ymzp - 2*ypzm - xpymzp - xmymzp - xpypzm - xmypzm );
							 
				uzx = 0.0625f*( 2*xpzp + 2*xmzm + xpypzp + xpymzp + xmypzm + xmymzm
							  - 2*xmzp - 2*xpzm - xpypzm - xpymzm - xmypzp - xmymzp );
							 
				// 3D mean curvature
				num = ux*ux*(uyy+uzz) + uy*uy*(uzz+uxx) + uz*uz*(uxx+uyy) -2*ux*uy*uxy - 2*uy*uz*uyz -2*uz*ux*uzx;
				den = (float)Math.sqrt(ux*ux + uy*uy + uz*uz);
				den = 2.0f*den*den*den;
				
				if (den>0) {
					curv[x][y][z] = num/den;
				} else {
					curv[x][y][z] = 0.0f;
				}
			} else {
				curv[x][y][z] = 0.0f;
			}
		}
		return curv;
	}

	/**
     *  compute the object volume
     */
    public static final int volume(boolean img[][][], int nx, int ny, int nz) {
		int count=0;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) count++;
		}

		return count;
	}
	
	/**
     *  compute the object volume
     */
    public static final int volume(boolean img[], int nx, int ny, int nz) {
		int count=0;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]) count++;
		}

		return count;
	}
	
	/**
     *  compute the object volume
     */
    public static final int volume(int img[][][], int lb, int nx, int ny, int nz) {
		int count=0;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]==lb) count++;
		}

		return count;
	}
	
    public static final int volume(int img[], int lb, int nx, int ny, int nz) {
		int count=0;
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]==lb) count++;
		}

		return count;
	}
	
	
	/**
     *  compute object volumes
     */
    public static final int[] volumes(int img[][][], int nlb, int nx, int ny, int nz) {
		int[] count=new int[nlb];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]>0) count[img[x][y][z]-1]++;
		}

		return count;
	}
	
	/**
     *  compute object volumes
     */
    public static final int[] volumes(int img[], int nlb, int nx, int ny, int nz) {
		int[] count=new int[nlb];
		
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (img[xyz]>0) count[img[xyz]-1]++;
		}

		return count;
	}
	
	
	/**
     *  compute the number of separate object parts
     */
    public static final int countParts(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int[][][] lb = null;
		
			 if (cObj== 6) lb = ObjectLabeling.connected6Object3D(img,nx,ny,nz);
		else if (cObj==18) lb = ObjectLabeling.connected18Object3D(img,nx,ny,nz);
		else if (cObj==26) lb = ObjectLabeling.connected26Object3D(img,nx,ny,nz);
		
		return (ObjectLabeling.countLabels(lb,nx,ny,nz)-1);
	}
	
	/**
     *  compute the number of holes in the object
     */
    public static final int countHoles(boolean img[][][], int nx, int ny, int nz, int cObj, int cBg) {
		int[][][] lb = null;
		boolean[][][] obj = new boolean[nx][ny][nz];
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			obj[x][y][z] = !img[x][y][z];
		}
			 if (cBg== 6) lb = ObjectLabeling.connected6Object3D(obj,nx,ny,nz);
		else if (cBg==18) lb = ObjectLabeling.connected18Object3D(obj,nx,ny,nz);
		else if (cBg==26) lb = ObjectLabeling.connected26Object3D(obj,nx,ny,nz);

		return (ObjectLabeling.countLabels(lb,nx,ny,nz)-2);
	}
	
	/**
	 *	compute the Hausdorff distance
	 *  between two objects
	 */
	public static final float hausdorffDistance(boolean[][][] obj1, boolean[][][] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float maxdist = 0.0f;
		boolean[][][] bound1, bound2;
		float[][][] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (bound1[x][y][z]) {
				if (dist2[x][y][z]>maxdist) maxdist = dist2[x][y][z];
			}
			if (bound2[x][y][z]) {
				if (dist1[x][y][z]>maxdist) maxdist = dist1[x][y][z];
			}
		}
		return (float)Math.sqrt(maxdist);	
	}
		
	/**
	 *	compute the Hausdorff distance
	 *  between two objects
	 */
	public static final float hausdorffDistance(boolean[] obj1, boolean[] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float maxdist = 0.0f;
		boolean[] bound1, bound2;
		float[] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (bound1[xyz]) {
				if (dist2[xyz]>maxdist) maxdist = dist2[xyz];
			}
			if (bound2[xyz]) {
				if (dist1[xyz]>maxdist) maxdist = dist1[xyz];
			}
		}
		return (float)Math.sqrt(maxdist);	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSurfaceDistance(boolean[][][] obj1, boolean[][][] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[][][] bound1, bound2;
		float[][][] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (bound1[x][y][z]) {
				// same boundary
				avgdist += (float)Math.sqrt(dist2[x][y][z]);
				avgnb++;
			}
			if (bound2[x][y][z]) {
				// same boundary
				avgdist += (float)Math.sqrt(dist1[x][y][z]);
				avgnb++;
			}
		}
		return avgdist/avgnb;	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSurfaceDistance(boolean[] obj1, boolean[] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[] bound1, bound2;
		float[] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (bound1[xyz]) {
				// same boundary
				avgdist += (float)Math.sqrt(dist2[xyz]);
				avgnb++;
			}
			if (bound2[xyz]) {
				// same boundary
				avgdist += (float)Math.sqrt(dist1[xyz]);
				avgnb++;
			}
		}
		return avgdist/avgnb;	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSurfaceDifference(boolean[][][] obj1, boolean[][][] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[][][] bound1, bound2;
		float[][][] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (bound1[x][y][z]) {
				if (obj2[x][y][z]) {
					// inside: positive
					avgdist += (float)Math.sqrt(dist2[x][y][z]);
					avgnb++;
				} else {
					//outside: negative
					avgdist +=-(float)Math.sqrt(dist2[x][y][z]);
					avgnb++;
				}
			}
			if (bound2[x][y][z]) {
				if (obj1[x][y][z]) {
					// inside: negative
					avgdist +=-(float)Math.sqrt(dist1[x][y][z]);
					avgnb++;
				} else {
					// outside: positive
					avgdist += (float)Math.sqrt(dist1[x][y][z]);
					avgnb++;
				}					
			}
		}
		return avgdist/avgnb;	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSurfaceDifference(boolean[] obj1, boolean[] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[] bound1, bound2;
		float[] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// find all correspondences: compute both distances
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (bound1[xyz]) {
				if (obj2[xyz]) {
					// inside: positive
					avgdist += (float)Math.sqrt(dist2[xyz]);
					avgnb++;
				} else {
					//outside: negative
					avgdist +=-(float)Math.sqrt(dist2[xyz]);
					avgnb++;
				}
			}
			if (bound2[xyz]) {
				if (obj1[xyz]) {
					// inside: negative
					avgdist +=-(float)Math.sqrt(dist1[xyz]);
					avgnb++;
				} else {
					// outside: positive
					avgdist += (float)Math.sqrt(dist1[xyz]);
					avgnb++;
				}					
			}
		}
		return avgdist/avgnb;	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSquaredSurfaceDistance(boolean[][][] obj1, boolean[][][] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[][][] bound1, bound2;
		float[][][] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// no used maps: average
		
		// find all correspondences
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (bound1[x][y][z]) {
				// same boundary
				avgdist += dist2[x][y][z];
				avgnb++;
			}
			if (bound2[x][y][z]) {
				// same boundary
				avgdist += dist1[x][y][z];
				avgnb++;
			}
		}
		return (float)Math.sqrt(avgdist/avgnb);	
	}
		
	/**
	 *	compute the average surface distance
	 *  between two objects
	 */
	public static final float averageSquaredSurfaceDistance(boolean[] obj1, boolean[] obj2, float rx, float ry, float rz, int nx, int ny, int nz) {
		float avgdist = 0.0f;
		float avgnb = 0.0f;
		float minres = Numerics.min(Numerics.min(rx,ry),rz);
		boolean[] bound1, bound2;
		float[] dist1,dist2;
		
		// crop the objects for speed
		ImageCropping crop = new ImageCropping(nx,ny,nz);
		crop.setBorderSize(1);
		crop.findCroppingBoundaries(obj1,obj2);
		obj1 = crop.cropImage(obj1);
		obj2 = crop.cropImage(obj2);
		nx = crop.mx();
		ny = crop.my();
		nz = crop.mz();
		
		// check for empty objects
		if ( volume(obj1,nx,ny,nz)==0 || volume(obj2,nx,ny,nz)==0 ) return -1.0f;
		
		// convert to boundaries
		bound1 = ObjectGeometry.objectBoundary(obj1,nx,ny,nz);
		bound2 = ObjectGeometry.objectBoundary(obj2,nx,ny,nz);
		
		// use distance maps
		dist1 = ObjectTransforms.voronoiFeatureSquaredDistance(bound1,nx,ny,nz);
		dist2 = ObjectTransforms.voronoiFeatureSquaredDistance(bound2,nx,ny,nz);
		
		// no used maps: average
		
		// find all correspondences
		for (int xyz=0;xyz<nx*ny*nz;xyz++) {
			if (bound1[xyz]) {
				// same boundary
				avgdist += dist2[xyz];
				avgnb++;
			}
			if (bound2[xyz]) {
				// same boundary
				avgdist += dist1[xyz];
				avgnb++;
			}
		}
		return (float)Math.sqrt(avgdist/avgnb);	
	}
		
	/**
	 *	check for convex points / concave points:
	 *  count the number of added boundaries (26-C).
	 */
	private static final float convexityScore(boolean[][][] obj, int x, int y, int z) {
		float boundaries=0;
		
		// 6-C
		if (obj[x-1][y][z]) boundaries++;
		if (obj[x+1][y][z]) boundaries++;
		if (obj[x][y-1][z]) boundaries++;
		if (obj[x][y+1][z]) boundaries++;
		if (obj[x][y][z-1]) boundaries++;
		if (obj[x][y][z+1]) boundaries++;
		
		// 18-C + implied 18-C
		if (obj[x-1][y-1][z]) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) ) boundaries++;
		if (obj[x-1][y+1][z]) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) ) boundaries++;
		if (obj[x+1][y-1][z]) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) ) boundaries++;
		if (obj[x+1][y+1][z]) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) ) boundaries++;
		
		if (obj[x][y-1][z-1]) boundaries++;
		else if ( (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;
		if (obj[x][y-1][z+1]) boundaries++;
		else if ( (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;
		if (obj[x][y+1][z-1]) boundaries++;
		else if ( (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;
		if (obj[x][y+1][z+1]) boundaries++;
		else if ( (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;
		
		if (obj[x-1][y][z-1]) boundaries++;
		else if ( (obj[x][y][z-1]) && (obj[x-1][y][z]) ) boundaries++;
		if (obj[x+1][y][z-1]) boundaries++;
		else if ( (obj[x][y][z-1]) && (obj[x+1][y][z]) ) boundaries++;
		if (obj[x-1][y][z+1]) boundaries++;
		else if ( (obj[x][y][z+1]) && (obj[x-1][y][z]) ) boundaries++;
		if (obj[x+1][y][z+1]) boundaries++;
		else if ( (obj[x][y][z+1]) && (obj[x+1][y][z]) ) boundaries++;
		
		// 26-C and implied 26-C
		if (obj[x-1][y-1][z-1]) boundaries++;
		else if ( (obj[x-1][y-1][z]) && (obj[x][y-1][z-1]) ) boundaries++;
		else if ( (obj[x][y-1][z-1]) && (obj[x-1][y][z-1]) ) boundaries++;
		else if ( (obj[x-1][y][z-1]) && (obj[x-1][y-1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x+1][y-1][z-1]) boundaries++;
		else if ( (obj[x+1][y-1][z]) && (obj[x][y-1][z-1]) ) boundaries++;
		else if ( (obj[x][y-1][z-1]) && (obj[x+1][y][z-1]) ) boundaries++;
		else if ( (obj[x+1][y][z-1]) && (obj[x+1][y-1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x-1][y+1][z-1]) boundaries++;
		else if ( (obj[x-1][y+1][z]) && (obj[x][y+1][z-1]) ) boundaries++;
		else if ( (obj[x][y+1][z-1]) && (obj[x-1][y][z-1]) ) boundaries++;
		else if ( (obj[x-1][y][z-1]) && (obj[x-1][y+1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x-1][y-1][z+1]) boundaries++;
		else if ( (obj[x-1][y-1][z]) && (obj[x][y-1][z+1]) ) boundaries++;
		else if ( (obj[x][y-1][z+1]) && (obj[x-1][y][z+1]) ) boundaries++;
		else if ( (obj[x-1][y][z+1]) && (obj[x-1][y-1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x-1][y+1][z+1]) boundaries++;
		else if ( (obj[x-1][y+1][z]) && (obj[x][y+1][z+1]) ) boundaries++;
		else if ( (obj[x][y+1][z+1]) && (obj[x-1][y][z+1]) ) boundaries++;
		else if ( (obj[x-1][y][z+1]) && (obj[x-1][y+1][z]) ) boundaries++;
		else if ( (obj[x-1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x+1][y-1][z+1]) boundaries++;
		else if ( (obj[x+1][y-1][z]) && (obj[x][y-1][z+1]) ) boundaries++;
		else if ( (obj[x][y-1][z+1]) && (obj[x+1][y][z+1]) ) boundaries++;
		else if ( (obj[x+1][y][z+1]) && (obj[x+1][y-1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y-1][z]) && (obj[x][y][z+1]) ) boundaries++;

		if (obj[x+1][y+1][z-1]) boundaries++;
		else if ( (obj[x+1][y+1][z]) && (obj[x][y+1][z-1]) ) boundaries++;
		else if ( (obj[x][y+1][z-1]) && (obj[x+1][y][z-1]) ) boundaries++;
		else if ( (obj[x+1][y][z-1]) && (obj[x+1][y+1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z-1]) ) boundaries++;

		if (obj[x+1][y+1][z+1]) boundaries++;
		else if ( (obj[x+1][y+1][z]) && (obj[x][y+1][z+1]) ) boundaries++;
		else if ( (obj[x][y+1][z+1]) && (obj[x+1][y][z+1]) ) boundaries++;
		else if ( (obj[x+1][y][z+1]) && (obj[x+1][y+1][z]) ) boundaries++;
		else if ( (obj[x+1][y][z]) && (obj[x][y+1][z]) && (obj[x][y][z+1]) ) boundaries++;

		return boundaries/26.0f;
	}

	/**
     *  compute the object convexity map.
     */
    public static final float[][][] convexityMap(boolean img[][][], int nx, int ny, int nz) {
		float[][][] score = new float[nx][ny][nz];
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			score[x][y][z] = convexityScore(img,x,y,z);
		}

		return score;
	}
	
	/**
     *  compute the object mean convexity
     */
    public static final float meanConvexity(boolean img[][][], int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += convexityScore(img,x,y,z);
					count++;
				}
			}
		}

		return score/(float)count;
	}
	
	/**
     *  compute the object mean absolute convexity
     */
    public static final float absConvexity(boolean img[][][], int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += Numerics.abs(convexityScore(img,x,y,z)-0.5f);
					count++;
				}
			}
		}

		return score/(float)count;
	}
	
	/**
     *  compute the object convexity variance
     */
    public static final float stdConvexity(boolean img[][][], float mean, int nx, int ny, int nz) {
		float score = 0.0f;
		int count=0;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (!img[x][y][z]) {
				boolean boundary=false;
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) {
					if (img[x+i][y+j][z+l]) { boundary = true; break; }
				}
				if (boundary) {
					score += (convexityScore(img,x,y,z)-mean)*(convexityScore(img,x,y,z)-mean);
					count++;
				}
			}
		}

		return (float)Math.sqrt(score/(float)(count-1));
	}
	
	/**
     *  compute the object center
     */
    public static final float[] center(boolean img[][][], int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(int img[][][], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(int img[], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		int ind = -1;
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			ind = x + y*nx + z*nx*ny;
			if (img[ind]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
    public static final float[] center(byte img[], int lb, int nx, int ny, int nz) {
		float[] center = new float[3];
		float count=0.0f;
		int ind = -1;
		for (int n=0;n<3;n++) center[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			ind = x + y*nx + z*nx*ny;
			if (img[ind]==lb) {
				center[0] += x;
				center[1] += y;
				center[2] += z;
				count++;
			}
		}
		if (count>0) {
			center[0] = center[0]/count;
			center[1] = center[1]/count;
			center[2] = center[2]/count;
		}
		return center;
	}
	
	/**
     *  compute the object deviation from the center
     */
    public static final float[] deviation(boolean img[][][], float[] center, int nx, int ny, int nz) {
		float[] std = new float[3];
		float count=0.0f;
		
		for (int n=0;n<3;n++) std[n] = 0.0f;
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				std[0] += (x-center[0])*(x-center[0]);
				std[1] += (y-center[1])*(y-center[1]);
				std[2] += (z-center[2])*(z-center[2]);
				count++;
			}
		}
		if (count>1) {
			std[0] = (float)Math.sqrt( std[0]/(count-1) );
			std[1] = (float)Math.sqrt( std[1]/(count-1) );
			std[2] = (float)Math.sqrt( std[2]/(count-1) );
		}
		return std;
	}

	/**
     *  compute the mean object distance to a point
     */
    public static final float meanDistance(boolean img[][][], float x0, float y0, float z0, int nx, int ny, int nz, float rx, float ry, float rz) {
		float dist = 0.0f;
		float count = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist += Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				count++;
			}
		}
		if (count>0) dist = dist/count;
		return dist;
	}

	/**
     *  compute the variance of object distance to a point
     */
    public static final float stdDistance(boolean img[][][], float x0, float y0, float z0, float mean, int nx, int ny, int nz, float rx, float ry, float rz) {
		float var = 0.0f;
		float count = 0.0f;
		float dist;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				var += (dist-mean)*(dist-mean);
				count++;
			}
		}
		if (count>1) var = var/(count-1);
		return (float)Math.sqrt(var);
	}

	/**
     *  compute the mean object distance to a point
     */
    public static final float[] meanDirection(boolean img[][][], float x0, float y0, float z0, int nx, int ny, int nz, float rx, float ry, float rz) {
		float[] dir = new float[3];
		for (int i=0;i<3;i++) dir[i] = 0.0f;
		float count = 0.0f;
		float dist = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				dir[0] += (x-x0)*rx/dist;
				dir[1] += (y-y0)*ry/dist;
				dir[2] += (z-z0)*rz/dist;
				count++;
			}
		}
		if (count>0) {
			dist = (float)Math.sqrt( dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2] );
			for (int i=0;i<3;i++) dir[i] = dir[i]/dist;
		}
		return dir;
	}

	/**
     *  compute the variance of object distance to a point
     */
    public static final float[] stdDirection(boolean img[][][], float x0, float y0, float z0, float[] mean, int nx, int ny, int nz, float rx, float ry, float rz) {
		float[] var = new float[3];
		for (int i=0;i<3;i++) var[i] = 0.0f;
		float count = 0.0f, dist;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				dist = (float)Math.sqrt( (x-x0)*rx*(x-x0)*rx+(y-y0)*ry*(y-y0)*ry+(z-z0)*rz*(z-z0)*rz );
				var[0] += ( (x-x0)*rx/dist - mean[0] )*( (x-x0)*rx/dist - mean[0] );
				var[1] += ( (y-y0)*ry/dist - mean[1] )*( (y-y0)*ry/dist - mean[1] );
				var[2] += ( (z-z0)*rz/dist - mean[2] )*( (z-z0)*rz/dist - mean[2] );
				count++;
			}
		}
		if (count>1) {
			for (int i=0;i<3;i++) var[i] = (float)Math.sqrt(var[i]/(count-1));
		}
		return var;
	}

	/**
     *  compute the area on the boundary between two objects.
	 */
    public static final float sharedBoundaryArea(float img[][][], float id1, float id2, int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]==id1) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if (img[x+i][y+j][z+l]==id2) {
						// shared boundary point: add the corresponding surface patch
						if (i*i==1) {
							area += ry*rz;
						} else if (j*j==1) {
							area += rz*rx;
						} else if (l*l==1) {
							area += rx*ry;
						}
					}
				}
			}
		}
		return area;
	}

	/**
     *  compute the area on the boundary between two objects.
	 */
    public static final float sharedBoundaryArea(byte img[][][], byte id1, byte id2, int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++) {
			if (img[x][y][z]==id1) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if (img[x+i][y+j][z+l]==id2) {
						// shared boundary point: add the corresponding surface patch
						if (i*i==1) {
							area += ry*rz;
						} else if (j*j==1) {
							area += rz*rx;
						} else if (l*l==1) {
							area += rx*ry;
						}
					}
				}
			}
		}
		return area;
	}

	/**
     *  compute the area on the object boundary
	 */
    public static final float boundaryArea(boolean img[][][], int nx, int ny, int nz, float rx, float ry, float rz) {
		float area = 0.0f;
		
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			if (img[x][y][z]) {
				for (int i=-1;i<=1;i++) for (int j=-1;j<=1;j++) for (int l=-1;l<=1;l++) if (i*i+j*j+l*l==1) {
					if ( x+i>=0 && x+i<nx && y+j>=0 && y+j<ny && z+l>=0 && z+l<nz ) {
						if (!img[x+i][y+j][z+l]) {
							// boundary point: add the corresponding surface patch
							if (i*i==1) {
								area += ry*rz;
							} else if (j*j==1) {
								area += rz*rx;
							} else if (l*l==1) {
								area += rx*ry;
							}
						}
					}
				}
			}
		}
		return area;
	}

}
