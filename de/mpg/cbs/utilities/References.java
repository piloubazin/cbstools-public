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

import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataUByte;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataFloat;
import edu.jhu.ece.iacl.jist.structures.image.VoxelType;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamVolume;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;

// TODO: Auto-generated Javadoc
/**
 * A collection of utilities and database for including references to articles, authors, etc.
 * 
 * @author Pierre-Louis Bazin
 */
public class References {

	private static final String[] authors = new String[]{
		"Pierre-Louis Bazin", "bazin@cbs.mpg.de","http://www.cbs.mpg.de/",
		"Christine Tardif", "tardif@cbs.mpg.de","http://www.cbs.mpg.de/",
		"Christopher Steele", "steele@cbs.mpg.de","http://www.cbs.mpg.de/",
		"unknown"
	};
	
	private static final String[] articles = new String[]{
		"bazin-nimg14", "Bazin et al., Computational Framework etc.",
		"unknown"
	};
	
	public static final AlgorithmAuthor getAuthor(String name) {
		for (int n=0;n<authors.length;n+=3) {
			if (name.equals(authors[n])) {
				return new AlgorithmAuthor(authors[n],authors[n+1],authors[n+2]);
			}
		}
		return new AlgorithmAuthor(name, "", "");
	}
}
