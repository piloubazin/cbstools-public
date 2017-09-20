package de.mpg.cbs.utilities;


import java.awt.*;
import java.awt.event.*;
import java.net.URL;

import javax.swing.*;

import java.io.*;
import java.net.*;
import java.util.*;
import java.util.jar.*;
import java.util.zip.*;


/**
 * A collection of utilities for I/O to bypass MIPAV and JIST settings
 * 
 * @author Pierre-Louis Bazin
 */
public class BasicInfo {
	/**
	 * Display message.
	 * 
	 * @param message the message
	 */
	public static final void displayMessage(String message){
		// console output
		System.out.print(message);
		System.out.flush();
	}
	
	/**
	 * Display error.
	 * 
	 * @param message the message
	 */
	public static final void displayError(String message){
		// console output
		System.err.print(message);
		System.err.flush();
	}
    
    /** return the maximum common base name for two strings (assumes the same beginning) */
    public static final String commonNameBase(String name1, String name2) {
    	int id=1;
    	while (name1.regionMatches(0, name2, 0, id)) id++;
    	
    	return name1.substring(0,id-1);
    }
    /** return the part of the base name not present in the second string (assumes the same beginning) */
    public static final String differentNameEnding(String name1, String name2) {
    	int id=1;
    	while (name1.regionMatches(0, name2, 0, id)) id++;
    	
    	return name1.substring(id-1);
    }
    
    // internal copy of MIPAV's internal orientation, etc. labels
    // hard-coded...
    public static final int AXIAL = 	0;	// FileInfoBase.AXIAL;
    public static final int CORONAL = 	1;	// FileInfoBase.CORONAL;
    public static final int SAGITTAL = 	2;	// FileInfoBase.SAGITTAL;
    
    public static final int R2L = 1; // FileInfoBase.ORI_R2L_TYPE;
    public static final int L2R = 2; // FileInfoBase.ORI_L2R_TYPE;
    public static final int P2A = 3; // FileInfoBase.ORI_P2A_TYPE;
    public static final int A2P = 4; // FileInfoBase.ORI_A2P_TYPE;
    public static final int I2S = 5; // FileInfoBase.ORI_I2S_TYPE;
    public static final int S2I = 6; // FileInfoBase.ORI_S2I_TYPE;
    
    // imported / adapted convenience classes from MIPAV, copied for self-containedness
    public static final boolean getBoolean(final StringTokenizer st) {
        final String str = st.nextToken();

        // returns true if str.equalsIgnoreCase( "true" )
        return new Boolean(str).booleanValue();
    }

    public static final float getFloat(final StringTokenizer st) {
        final String str = st.nextToken();

        return Float.parseFloat(str);
    }

    public static final double getDouble(final StringTokenizer st) {
        final String str = st.nextToken();

        return Double.parseDouble(str);
    }

    public static final int getInt(final StringTokenizer st) {
        final String str = st.nextToken();

        return Integer.parseInt(str);
    }

}
