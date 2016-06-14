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

import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.MipavUtil;
import gov.nih.mipav.view.ViewUserInterface;
import gov.nih.mipav.model.structures.*;
import gov.nih.mipav.model.file.*;

import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.PrintWriter;

import java.net.*;
import java.io.*;
import java.util.*;
import java.lang.*;
import java.text.*;


// TODO: Auto-generated Javadoc
/**
 * A class of utilities for csv I/O with MIPAV, etc.
 * 
 * @author Juliane Dinse
 */
public class CSVInOut {

	public static void main(String[] args){
		
		int lengthA = args.length;
		String 
		
	}
	
	private String delim = ",";
	
	private final void addStatisticsToFile(String name, ArrayList<String> output, String lbline) {

		// open the file
		ArrayList<String> 	previous = loadStatisticsFile(name);

		// merge the output
		appendStatistics(previous, output, lbline);

		// save the result
		writeStatisticsFile(previous, name);
	}

	private	final void appendStatistics(ArrayList<String> main, ArrayList<String> added, String lbline) {
		for (int n=0;n<added.size();n++) {
			// extract statistics type
			String type = added.get(n).substring(0,added.get(n).indexOf(delim));
			System.out.println(added.get(n));
			System.out.println(type);

			// find the last line with this type
			int last=-1;
			for (int m=0;m<main.size();m++) {
				if (main.get(m).indexOf(delim)>-1)
					if (main.get(m).substring(0,main.get(m).indexOf(delim)).equals(type)) last = m;
			}
			if (last>-1) {
				main.add(last+1, added.get(n));
			} else {
				// add a space after each different statistic as well as labels before
				main.add(lbline);
				main.add(added.get(n));
				main.add(" \n");
			}
		}
	}

	private final ArrayList<String> loadStatisticsFile(String name) {
		ArrayList<String> list = new ArrayList();
		try {
			System.out.println("reading previous statistic file: "+name);
			File f = new File(name);
			FileReader fr = new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();

			// Exact corresponding template for first line ?
			if (!line.startsWith("MIPAV Volumetric Statistics File")) {
				System.out.println("not a proper MIPAV statistics file");
				br.close();
				fr.close();
				return null;
			}
			line = br.readLine();
			while (line!=null) {
				list.add(line+"\n");
				line = br.readLine();
				System.out.println(line);
			}
			br.close();
			fr.close();
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
		return list;	
	}

	private final void writeStatisticsFile(ArrayList<String> list, String name) {
		try {
			File f = new File(name);
			FileWriter fw = new FileWriter(f);
			PrintWriter pw = new PrintWriter( fw );
			pw.write("MIPAV Volumetric Statistics File\n");

			for (int n=0;n<list.size();n++) {
				pw.write(list.get(n));
				System.out.print(list.get(n));
			}
			pw.close();
			fw.close();
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
	}	





}
