/**
 * Java Image Science Toolkit (JIST)
 *
 * Image Analysis and Communications Laboratory &
 * Laboratory for Medical Image Computing &
 * The Johns Hopkins University
 * 
 * http://www.nitrc.org/projects/jist/
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.  The license is available for reading at:
 * http://www.gnu.org/copyleft/lgpl.html
 *
 */
package edu.jhu.ece.iacl.jist.pipeline;

import java.awt.Component;
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.Vector;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.ProgressMonitor;
import javax.swing.SwingWorker;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;

import jist.modules.JistModulesLocation;

import org.apache.commons.io.FileUtils;
import org.json.JSONException;

import edu.jhu.ece.iacl.jist.io.FileReaderWriter;
import edu.jhu.ece.iacl.jist.io.MipavController;
import edu.jhu.ece.iacl.jist.pipeline.graph.PipeAlgorithmFactory;
import edu.jhu.ece.iacl.jist.pipeline.graph.PipeModuleFactory;
import edu.jhu.ece.iacl.jist.pipeline.parameter.ParamCollection;
import edu.jhu.ece.iacl.jist.pipeline.parser.JSONPluginParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.LoniPipeParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.MipavPluginParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.MipavScriptParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.ScriptParser;
import edu.jhu.ece.iacl.jist.pipeline.parser.ShellScriptXmlParser;
import edu.jhu.ece.iacl.jist.pipeline.tree.PackageNode;
import edu.jhu.ece.iacl.jist.pipeline.view.input.Refresher;
import edu.jhu.ece.iacl.jist.plugins.MedicAlgorithmExecutableAdapter;
import edu.jhu.ece.iacl.jist.plugins.MedicAlgorithmGeneralScriptAdapter;
import edu.jhu.ece.iacl.jist.plugins.MedicAlgorithmMipavAdapter;
import edu.jhu.ece.iacl.jist.utility.FileUtil;
import edu.jhu.ece.iacl.jist.utility.JistLogger;
import gov.nih.mipav.view.MipavUtil;
import gov.nih.mipav.view.Preferences;
import gov.nih.mipav.view.ViewUserInterface;
import gov.nih.mipav.view.dialogs.ActionDiscovery;

// TODO: Auto-generated Javadoc
/**
 * Interface to open, edit, and save module definitions.
 * 
 * @author Blake Lucas (bclucas@jhu.edu)
 */
public class PipeLibrary {

	/** The Constant RootPackageName. */
	//	private static final String RootPackageName="MAPS";

	/**
	 * Compare two classes according to their simple names.
	 * 
	 * @author Blake Lucas
	 */
	private static class ClassComparator implements Comparator<Class> {

		/**
		 * Compare classes alphabetically.
		 * 
		 * @param c0 the c0
		 * @param c1 the c1
		 * 
		 * @return the int
		 */
		public int compare(Class c0, Class c1) {
			return c0.getSimpleName().compareTo(c1.getSimpleName());
		}
	}

	/**
	 * Gets the host name.
	 * 
	 * @return the host name
	 */
	public static String getHostName(){
		String host=System.getenv("HOSTNAME");
		if(host!=null){
			return host;
		} else {
			return "localhost";
		}
	}

	/**
	 * Worker to find and initialize module definitions in library directory.
	 * 
	 * @author Blake Lucas (bclucas@jhu.edu)
	 */
	public static class InitializeTask extends SwingWorker<Void, Void> {

		/** The parent. */
		protected Component parent;

		/**
		 * Instantiates a new initialize task.
		 * 
		 * @param parent the parent
		 */
		public InitializeTask(Component parent) {
			this.parent = parent;
		}

		/**
		 * Load library in background.
		 * 
		 * @return always returns null
		 */
		protected Void doInBackground() {
			try {
				ProgressMonitor monitor = new ProgressMonitor(parent, "Initializing Plug-in Library",
						"Initializing...", 0, 100);
				monitor.setMillisToDecideToPopup(0);
				Refresher.getInstance().disable();
				PipeLibrary.getInstance().initModules(monitor);
				Refresher.getInstance().enable();
				monitor.close();
			} catch (Exception e) {
				MipavUtil.displayError(e.getMessage());
				e.printStackTrace();
			}
			return null;
		}
	}

	/**
	 * Worker to rebuild library from module definitions in the library
	 * directory.
	 * 
	 * @author Blake Lucas (bclucas@jhu.edu)
	 */
	public static class RebuildTask extends SwingWorker<Void, Void> {

		/** The parent. */
		protected Component parent;

		/** The remove. */
		protected boolean remove;

		/**
		 * Instantiates a new rebuild task.
		 * 
		 * @param parent the parent
		 * @param remove the remove
		 */
		public RebuildTask(Component parent, boolean remove) {
			this.parent = parent;
			this.remove = remove;
		}

		/**
		 * Rebuild library in background.
		 * 
		 * @return always returns null
		 */
		protected Void doInBackground() {
			try {
				ProgressMonitor monitor = new ProgressMonitor(parent, "Rebuilding Plug-in Library", "Rebuilding...", 0,
						100);
				monitor.setMillisToDecideToPopup(0);
				JistLogger.logOutput(JistLogger.INFO, getClass().getCanonicalName()+"\t"+"Remove Existing Modules");
				if (remove) {
					PipeLibrary.getInstance().removeAtomicModules(monitor);
				}
				PipeLibrary.getInstance().discoverAllModules(monitor);
				PipeLibrary.getInstance().initModules(monitor);
				monitor.close();
			} catch (Exception e) {
				e.printStackTrace();
				MipavUtil.displayError("Error Occurred while Rebuilding Libarary:"+e.getMessage());

			}
			return null;
		}
	}

	/** text string to identify library location in MIPAV parameters. */
	private static final String LIBRARY_KEY = "JIST_Library_Directory";
	private static final String JSON_LIBRARY_KEY="JIST_External_Library_Directory";

	/** Singleton reference to library. */
	protected static PipeLibrary library = null;

	/**
	 * Open file chooser to select directory.
	 * 
	 * @param oldDir the old directory
	 * 
	 * @return absolute path of the file
	 */
	public static File askForLibraryDirectory(File oldDir) {
//				RuntimeException e = new RuntimeException("askForLibraryDirectory: ARG!");
//				StackTraceElement []el = e.getStackTrace();
//				for(StackTraceElement m:el)
		
//				JistLogger.logError(JistLogger.SEVERE, m.toString());
//				e.printStackTrace(System.out);/
//				System.out.flush();

		
		JFileChooser openDialog = new JFileChooser();
		openDialog.setSelectedFile(MipavController.getDefaultWorkingDirectory());
		openDialog.setDialogTitle("Select Library Directory");
		openDialog.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		
		if (oldDir != null) {
			openDialog.setCurrentDirectory(oldDir);
		}
		openDialog.setDialogType(JFileChooser.OPEN_DIALOG);
		int returnVal = openDialog.showOpenDialog(null);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			MipavController.setDefaultWorkingDirectory(openDialog.getSelectedFile().getParentFile());
			return openDialog.getSelectedFile();
		} else {
			return null;
		}
	}

	/**
	 * Gets the all classes.
	 * 
	 * @param rootPckg the root pckg
	 * 
	 * @return the all classes
	 */
	public static Vector<Class> getAllClasses(String rootPckg) {
		return getAllClasses(rootPckg, null);
	}
	
	

	/**
	 * Get all classes starting from root package name.
	 * 
	 * @param rootPckg root package name
	 * @param simpleName the simple name
	 * 
	 * @return the all classes
	 */
	public static Vector<Class> getAllClasses(String rootPckg,String simpleName) {
		String[] strs = System.getProperty("java.class.path").split(File.pathSeparator);
		LinkedList<File> dirs = new LinkedList<File>();
		LinkedList<String> packs = new LinkedList<String>();
		Vector<Class> classes = new Vector<Class>();
		String pckg;
		for (String str : strs) {
			File f = new File(str, rootPckg.replace(".", File.separator));
			if (f.exists() && f.isDirectory()) {
				dirs.add(f);
				packs.add(rootPckg);
			}
		}
		while (dirs.size() > 0) {
			File dir = dirs.removeFirst();
			String pack = packs.removeFirst();
			File[] files = dir.listFiles();
			// Test all files in directory
			for (File f : files) {
				if (f.isDirectory()) {
					dirs.add(f);
					if (pack.length() > 0) {
						packs.add(pack + "." + FileReaderWriter.getFileName(f));
					} else {
						packs.add(FileReaderWriter.getFileName(f));
					}
				} else if (FileReaderWriter.getFileExtension(f).equals("class")) {
					if (pack.length() > 0) {
						pckg = pack + "." + FileReaderWriter.getFileName(f);
					} else {
						pckg = FileReaderWriter.getFileName(f);
					}
					try {
						Class c = Class.forName(pckg);
						if (c != null) {
							if(simpleName==null||c.getSimpleName().equals(simpleName)){
								classes.add(c);
							}
						}
					} catch (ClassNotFoundException e) {
					} catch (UnsatisfiedLinkError e) {
					} catch (VerifyError e) {
					} catch (Error e) {
					}
				}
			}
		}
		return classes;
	}

	/**
	 * Get singleton reference to library.
	 * 
	 * @return library
	 */
	public static PipeLibrary getInstance() {
		if (library == null) {
			library = new PipeLibrary();
		}
		return library;
	}

	/** Location of library directory. */
	private File libraryPath = null;
	private File jsonLibraryPath = null;


	/** List of algorithms discovered from library. */
	private Vector<PipeAlgorithmFactory> algos;

	/** Tree node to represent directory. */
	private PackageNode algoRoot = null;

	/**
	 * Default constructor.
	 */
	protected PipeLibrary() {
	}

	/**
	 * Determine of library already contains an instance of this algorithm.
	 * 
	 * @param algo algorithm
	 * 
	 * @return true if algorithm exists in library
	 */
	private boolean containsAlgorithm(PipeAlgorithm algo) {
		Class algoClass = algo.getAlgorithmClass();
		if (MedicAlgorithmExecutableAdapter.class.equals(algoClass)
				|| MedicAlgorithmMipavAdapter.class.equals(algoClass)
				|| MedicAlgorithmGeneralScriptAdapter.class.equals(algoClass)) {
			return false;
		}
		for (PipeAlgorithmFactory child : algos) {
			Class childClass = child.getModuleClass();
			if ((childClass != null) && (algoClass != null) && childClass.equals(algoClass)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * The Class CompareFiles.
	 */
	public static class CompareFiles implements Comparator<File>{

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(File f1, File f2) {
			if(f1.isDirectory()){
				if(f2.isDirectory()){
					return f1.getName().compareTo(f2.getName());
				} else {
					return 1;
				}
			} else {
				if(f2.isDirectory()){
					return -1;
				} else {
					return f1.getName().compareTo(f2.getName());
				}
			}
		}

	}

	/**
	 * The Class CompareFactories.
	 */
	public static class CompareFactories implements Comparator<PipeModuleFactory>{

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(PipeModuleFactory fact1, PipeModuleFactory fact2) {
			return fact1.getModuleName().compareTo(fact2.getModuleName());
		}

	}

	/**
	 * Discover all MAPS and LONI files.
	 * 
	 * @param monitor progress monitor
	 */

	public void discoverAllMapsAndLoniModules(ProgressMonitor monitor) {
		LoniPipeParser loniParser = new LoniPipeParser();
		ScriptParser mipavParser = new MipavScriptParser();
		ShellScriptXmlParser scriptParser = new ShellScriptXmlParser();
		LinkedList<File> dirs = new LinkedList<File>();
		dirs.add(libraryPath);
		File d;
		monitor.setProgress(0);
		monitor.setMaximum(1);
		while (dirs.size() > 0) {
			if ((monitor != null) && monitor.isCanceled()) {
				break;
			}
			d = dirs.remove();
			if (monitor != null) {
				if (dirs.size() > monitor.getMaximum()) {
					monitor.setMaximum(dirs.size());
				}
				monitor.setNote(d.getName());
			}
			if (d == null) {
				continue;
			}
			File[] files = d.listFiles();
			if (files == null) {
				continue;
			}
			if (monitor != null) {
				monitor.setMaximum(files.length);
			}
			for (File f : files) {
				if ((monitor != null) && monitor.isCanceled()) {
					break;
				}
				if (!f.isHidden()) {
					if (f.isDirectory()) {
						// Add directory to search path
						dirs.add(f);
					} else {
						String ext = FileReaderWriter.getFileExtension(f);
						if (ext.equals("pipe")) {
							PipeAlgorithm algo = loniParser.parsePipeAlgorithm(f);
							if (algo != null) {
								f = new File(f.getParentFile(), FileReaderWriter.getFileName(f) + "."+JistPreferences.getDefaultModuleExtension());
								if (!f.exists()) {
									algo.write(f);
								}
							}
						} else if (ext.equals("sct")) {
							PipeAlgorithm algo = mipavParser.parsePipeAlgorithm(f);
							if (algo != null) {
								f = new File(f.getParentFile(), FileReaderWriter.getFileName(f) + "."+JistPreferences.getDefaultModuleExtension());
								if (!f.exists()) {
									algo.write(f);
								}
							}
						}else{
							try{
								PipeAlgorithm algo = scriptParser.parsePipeAlgorithm(f);
								if (algo != null) {
									f = new File(f.getParentFile(), FileReaderWriter.getFileName(f) + "."+JistPreferences.getDefaultModuleExtension());
									if (!f.exists()) {
										algo.write(f);
									}
								}
							}catch(Exception e){
								JistLogger.logError(JistLogger.SEVERE, "Unable to load: " + f +" as a script");
							}
						}
					}
				}
			}
			if (monitor != null) {
				monitor.setProgress(monitor.getMaximum() - dirs.size());
			}
		}
	}

	/**
	 * Discover all modules that can be loaded by MAPS.
	 * 
	 * @param monitor progress monitor
	 */
	public void discoverAllModules(ProgressMonitor monitor) {
		loadPreferences();
		algos = new Vector<PipeAlgorithmFactory>();
		algoRoot = new PackageNode(getLibraryPath(), "Algorithms");
		algoRoot.setEditable(false);

		try {
			JistLogger.logOutput(JistLogger.INFO, "Discover: Maps&Loni");
			discoverAllMapsAndLoniModules(monitor);
		} catch(Throwable e) {
			JistLogger.logError(JistLogger.SEVERE, e.toString());
		}
		try {
			JistLogger.logOutput(JistLogger.INFO, "Discover: MIPAV");
			discoverMIPAVModules(monitor);
		} catch(Throwable e) {
			JistLogger.logError(JistLogger.SEVERE, e.toString());
		}
		try{
			JistLogger.logOutput(JistLogger.INFO, "Discover: INTERNAL");
			discoverInternalMapsModules(monitor);
		} catch(Throwable e) {			
			JistLogger.logError(JistLogger.SEVERE, e.toString());
		}
		try{
			JistLogger.logOutput(JistLogger.INFO, "Discover: EXTERNAL");
			discoverExternalWrappedModules(monitor);
		} catch(Throwable e) {			
			JistLogger.logError(JistLogger.SEVERE, e.toString());
		}
		try{
			JistLogger.logOutput(JistLogger.INFO, "Discover: EXTERNAL JSON Modules");
			discoverExternalJSONWrappedModules(monitor);
		} catch(Throwable e) {			
			JistLogger.logError(JistLogger.SEVERE, e.toString());
		}
		
		JistLogger.logOutput(JistLogger.INFO, "Discover: DONE");
	}

	private void discoverMIPAVModules(ProgressMonitor monitor) {
		JistLogger.logOutput(JistLogger.FINE,"discoverMIPAVModules:");
		final Vector<Class<ActionDiscovery>> classes = ViewUserInterface.getDiscoverableActionList();
		File f;	
		System.out.println(libraryPath.getAbsolutePath());
		File discoveredPath = libraryPath;//new File(libraryPath, File.separatorChar+RootPackageName*/);
		if (!discoveredPath.exists()) {
			discoveredPath.mkdir();
		} 
		monitor.setMaximum(classes.size());
		monitor.setProgress(0);        

		for (final Class<ActionDiscovery> c : classes) {
			
			if (!c.getCanonicalName().equals("gov.nih.mipav.view.dialogs.JDialogVOILogicalOperations") && !c.getCanonicalName().equals("gov.nih.mipav.view.dialogs.JDialogSM2")){
				JistLogger.logOutput(JistLogger.FINE,"discoverMIPAVModules:"+c.getCanonicalName());
	
				if ((monitor != null) && monitor.isCanceled()) {
					break;
				}
				monitor.setNote(c.toString());
				JistLogger.logOutput(JistLogger.INFO, getClass().getCanonicalName()+"\t"+c.getName());
				PipeAlgorithm algo;
				try {
					algo = MipavPluginParser.parsePipeAlgorithm((ActionDiscovery)c.newInstance());
	
					File path=createPath(discoveredPath, algo.getPackage()+(algo.getCategory()!=null?"."+algo.getCategory():""));
					ProcessingAlgorithm method=algo.getAlgorithm();
					if(method!=null){
	//					f = new File(path, FileUtil.forceSafeFilename(method.getClass().getSimpleName()) + "."+JistPreferences.getDefaultModuleExtension());
						f = new File(path, FileUtil.forceSafeFilename(algo.getName()) + "."+JistPreferences.getDefaultModuleExtension());
						algo.write(f); 
					}
				} catch (InstantiationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}//new PipeAlgorithm((ActionDiscovery) obj);
			}

		}
	}



	/**
	 * The Class CompareSimpleNames.
	 */
	private static class CompareSimpleNames implements Comparator<Class>{

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Class cls1, Class cls2) {
			return (cls1.getSimpleName().compareTo(cls2.getSimpleName()));
		}

	}

	/**
	 * Discover new algorithms along source path.
	 * 
	 * @param monitor the monitor
	 */
	public void discoverInternalMapsModules(ProgressMonitor monitor) {
		JistLogger.logOutput(JistLogger.FINE,"discoverInternalMapsModules");
		Vector<Class> classes = PipeLibrary.getAllClasses("jist.modules");
		for(Class c : classes) {
			JistLogger.logOutput(JistLogger.FINE,"Found: "+c.getName());
			Object o;
			try {
				o = c.newInstance();
				if(o instanceof JistModulesLocation) {
					JistModulesLocation jml = (JistModulesLocation)o;
					String []classPath = jml.getValidModuleJavaPaths();
					for(String path : classPath) {
						JistLogger.logOutput(JistLogger.FINE,"Found path: "+path);
						discoverInternalMapsModulesPrivate(monitor,path);
					}
				}
			} catch (InstantiationException e) {
			} catch (IllegalAccessException e) {
			}
		}
	}
	

	private void discoverInternalMapsModulesPrivate(ProgressMonitor monitor, String searchJavaClassPath) {
		Vector<Class> classes = PipeLibrary.getAllClasses(searchJavaClassPath);
		File f;	
		File discoveredPath = libraryPath;//new File(libraryPath, File.separatorChar+RootPackageName*/);
		if (!discoveredPath.exists()) {
			discoveredPath.mkdir();
		} else {
			for (int i = 0; i < algoRoot.getChildCount(); i++) {
				TreeNode n = algoRoot.getChildAt(i);
				if (n instanceof PackageNode) {
					Object obj = ((PackageNode) n).getUserObject();
				}
			}
		}
		monitor.setMaximum(classes.size());
		monitor.setProgress(0);

		for (Class c : classes) {
			if ((monitor != null) && monitor.isCanceled()) {
				break;
			}
			monitor.setNote(c.toString());
			
			Class c2 = c.getSuperclass();
			boolean cont_searching = true;
			
			// we can skip abstract classes
			if (Modifier.isAbstract(c.getModifiers()))
				cont_searching = false;
			
			while (cont_searching) {
				if (ProcessingAlgorithm.class.equals(c2)) {

					JistLogger.logOutput(JistLogger.INFO, getClass().getCanonicalName()+"\t"+c.getName());
					Object obj=null;
					try {
						obj = c.newInstance();
					} catch (InstantiationException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IllegalAccessException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (Exception e) {
						e.printStackTrace();
						obj = null;
						
					}
					if(obj!=null){
						PipeAlgorithm algo = new PipeAlgorithm((ProcessingAlgorithm) obj);
						if (true) {//!containsAlgorithm(algo)) {
							File path=createPath(discoveredPath, algo.getPackage()+(algo.getCategory()!=null?"."+algo.getCategory():""));
							ProcessingAlgorithm method=algo.getAlgorithm();
							if(method!=null){
								f = new File(path, FileUtil.forceSafeFilename(method.getClass().getSimpleName()) + "."+JistPreferences.getDefaultModuleExtension());
								algo.write(f); 
							}
						}
					}
				}
				
				// continue down the rabbit hole
				if (c2 != null)
					c2 = c2.getSuperclass();
				else
					cont_searching = false;
			}
		}
	}
public void discoverExternalJSONWrappedModules(ProgressMonitor monitor) throws FileNotFoundException, JSONException {
		
		
		File extModuleDir = new File(Preferences.getProperty("JIST_External_Library_Directory"));
		if (!extModuleDir.exists() || !extModuleDir.isDirectory())
			return;
		
		
		Collection<File> moduleJSONList = FileUtils.listFiles(extModuleDir, new String[]{"json"}, true);
		JistLogger.logOutput(JistLogger.INFO, "found " + moduleJSONList.size() + " moduleJSON files");
		for(File f:moduleJSONList){
			JSONPluginParser parser = new JSONPluginParser(f);
			JistLogger.logOutput(JistLogger.INFO, "moduleJSON file: " + f);
			parser.pluginInitTest(parser.pluginInitTest);
			parser.writePlugin();
		
		}
	}
	public void discoverExternalWrappedModules(ProgressMonitor monitor) {
		
		File discoveredPath = libraryPath;//new File(libraryPath, File.separatorChar+RootPackageName*/);
		if (!discoveredPath.exists()) {
			discoveredPath.mkdir();
		} else {
			for (int i = 0; i < algoRoot.getChildCount(); i++) {
				TreeNode n = algoRoot.getChildAt(i);
				if (n instanceof PackageNode) {
					Object obj = ((PackageNode) n).getUserObject();
				}
			}
		}
		
		File extModuleDir = JistPreferences.getPreferences().getExternalModuleDirectory();
		
		if (!extModuleDir.exists() || !extModuleDir.isDirectory())
			return;
		
		ShellScriptXmlParser parser = new ShellScriptXmlParser();
		String ext = parser.getFileExtensionFilter().getPreferredExtension();
		Collection<File> moduleXmlList = FileUtils.listFiles(extModuleDir, new String[]{ext}, true);
		JistLogger.logOutput(JistLogger.INFO, "found " + moduleXmlList.size() + " moduleXML files");
		for(File f:moduleXmlList){
			JistLogger.logOutput(JistLogger.INFO, "moduleXML file: " + f);
			PipeAlgorithm algo = parser.parsePipeAlgorithm(f);
			String relativePath = f.getParent().replaceFirst(extModuleDir.getPath()+File.separator, "");
			File path=createPath(discoveredPath, "External"+File.separator+relativePath);
			if(!path.exists()){
				path.mkdir();
			}
			ProcessingAlgorithm method=algo.getAlgorithm();
			if(method!=null){
				File fout = new File(path,FileUtil.forceSafeFilename(algo.getName()) + "."+JistPreferences.getDefaultModuleExtension());
				algo.write(fout); 
			}
		}
	}

	/**
	 * Creates the path.
	 * 
	 * @param discoveredPath the discovered path
	 * @param category the category
	 * 
	 * @return the file
	 */
	public File createPath(File discoveredPath,String category){
		if(category==null){
			return discoveredPath;
		}
		File f=new File(discoveredPath,category.replace('.', File.separatorChar).replace('_',' '));
		if(!f.exists()){
			f.mkdirs();
		}
		return f;
	}

	/**
	 * Get location of library.
	 * 
	 * @return library directory
	 */
	public File getLibraryPath() {
		return libraryPath;
	}


	/**
	 * Get tree representing algorithm factory hierarchy.
	 * 
	 * @return root tree node
	 */
	public DefaultMutableTreeNode getTreeRoot() {
		return algoRoot;
	}

	/**
	 * Find and initialize modules from library directory.
	 * 
	 * @param monitor progress monitor
	 */
	public void initModules(ProgressMonitor monitor) {
		
		if (loadPreferences()) {
			discoverAllModules(monitor);
		}
		algos = new Vector<PipeAlgorithmFactory>();
		algoRoot = new PackageNode(getLibraryPath(), "Algorithms");
		algoRoot.setEditable(false);
		LinkedList<File> dirs = new LinkedList<File>();
		LinkedList<DefaultMutableTreeNode> nodes = new LinkedList<DefaultMutableTreeNode>();
		dirs.add(libraryPath);
		nodes.add(algoRoot);
		File d;
		DefaultMutableTreeNode n;
		PipeAlgorithmFactory factory;
		monitor.setProgress(0);
		monitor.setMaximum(1);
		LinkedList<File> brokenList=new LinkedList<File>();
		while (dirs.size() > 0) {
			if ((monitor != null) && monitor.isCanceled()) {
				break;
			}
			
			d = dirs.remove();
			
			if (d == null || !d.isDirectory())
				continue;
			
			if (monitor != null) {
				if (dirs.size() > monitor.getMaximum()) {
					monitor.setMaximum(dirs.size());
				}
				monitor.setNote(d.getName());
			}
			n = nodes.remove();
			
			File[] files = d.listFiles();
			if (files == null || files.length == 0)
				continue;

			Arrays.sort(files,new CompareFiles());

			Vector<PipeAlgorithmFactory> addList=new Vector<PipeAlgorithmFactory>();
			for (File f : files) {
				if ((monitor != null) && monitor.isCanceled()) {
					break;
				}
				if (!f.isHidden()) {
					if (f.isDirectory()) {
						// Add directory to search path
						dirs.add(f);
						PackageNode child = new PackageNode(f);
						n.add(child);
						nodes.add(child);
					} else {
						String ext = FileReaderWriter.getFileExtension(f);
						if (ext.equals(JistPreferences.getDefaultModuleExtension())) {
							try {
								factory = new PipeAlgorithmFactory(f);
								addList.add(factory);
							} catch (ClassNotFoundException e) {
								brokenList.add(f);
								JistLogger.logError(JistLogger.SEVERE, getClass().getCanonicalName()+e.getMessage());
							} 
						}
					}
				}
			}

			Collections.sort(addList,new CompareFactories());
			for(PipeAlgorithmFactory fact:addList){
				n.add(fact.createTreeNode());
				algos.add(fact);
			}
			if (monitor != null) {
				monitor.setProgress(monitor.getMaximum() - dirs.size());
			}

		}
		if(brokenList.size()>0){
			int opt = JOptionPane.showOptionDialog(null, "Warning: Found "+brokenList.size()+" broken module(s). Do you want to delete them now?",null,JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
			if(opt==0){
				for(File f:brokenList){
					if(f.delete()){
						JistLogger.logOutput(JistLogger.INFO, getClass().getCanonicalName()+"\t"+"Deleted "+f.getAbsolutePath());
					}
				}
			}
		}
	}

	/**
	 * Load library path from MIPAV preferences. If entry does not exist, ask
	 * user to specify library.
	 * 
	 * Display a user popup dialog if the directory cannot be found
	 * 
	 * @return true if library had to be specified by user
	 */
	public boolean loadLibraryPath() {
		return loadLibraryPath(true);
	}
	
	public boolean loadJSONLibraryPath() {
		return loadJSONLibraryPath(true);
	}
	/**
	 * 
	 * 
	 * Load library path from MIPAV preferences. If entry does not exist, ask
	 * user to specify library.
	 * 
	 * @return true if library had to be specified by user
	 */
	public boolean loadLibraryPath(boolean silentFailure) {
		if (libraryPath != null) {
			return false;
		}
		String storedPath = Preferences.getProperty(LIBRARY_KEY);
		File defaultLibPath = new File(System.getProperty("user.home"));		
		if (((storedPath == null) || !(new File(storedPath)).exists())) {
			if (!silentFailure && !MipavController.isQuiet()) {
				File dir = MipavController.getDefaultWorkingDirectory();
				libraryPath = askForLibraryDirectory(dir.getParentFile());
				if (libraryPath == null) {
					libraryPath = new File(defaultLibPath,"jist-lib/");
					if (!libraryPath.exists())
						libraryPath.mkdirs();
					JistLogger.logError(JistLogger.WARNING, "Could not find library, using working directory "
									+ libraryPath.getAbsolutePath());
				} else {
					Preferences.setProperty(LIBRARY_KEY,
							libraryPath.getAbsolutePath());
					JistLogger.logError(JistLogger.WARNING, "Could not find library, using user specified directory "
									+ libraryPath.getAbsolutePath());
				}
			} else {
				libraryPath = new File(defaultLibPath,"jist-lib/");
				if (!libraryPath.exists())
					libraryPath.mkdirs();
				JistLogger.logError(JistLogger.WARNING, "Could not find library, using working directory "
								+ libraryPath.getAbsolutePath());
			}

			return true;
		} else {
			libraryPath = new File(storedPath);
			return false;
		}
	}
	
	public boolean loadJSONLibraryPath(boolean silentFailure) {
		if (jsonLibraryPath != null) {
			return false;
		}
		String storedPath = Preferences.getProperty(JSON_LIBRARY_KEY);
		File defaultLibPath = new File(System.getProperty("user.home"));		
		if (((storedPath == null) || !(new File(storedPath)).exists())) {
			if (!silentFailure && !MipavController.isQuiet()) {
				File dir = MipavController.getDefaultWorkingDirectory();
				libraryPath = askForLibraryDirectory(dir.getParentFile());
				if (libraryPath == null) {
					libraryPath = new File(defaultLibPath,"jist-json-lib/");
					if (!libraryPath.exists())
						libraryPath.mkdirs();
					JistLogger.logError(JistLogger.WARNING, "Could not find library, using working directory "
									+ libraryPath.getAbsolutePath());
				} else {
					Preferences.setProperty(JSON_LIBRARY_KEY,
							libraryPath.getAbsolutePath());
					JistLogger.logError(JistLogger.WARNING, "Could not find library, using user specified directory "
									+ libraryPath.getAbsolutePath());
				}
			} else {
				libraryPath = new File(defaultLibPath,"jist-json-lib/");
				if (!libraryPath.exists())
					libraryPath.mkdirs();
				JistLogger.logError(JistLogger.WARNING, "Could not find library, using working directory "
								+ libraryPath.getAbsolutePath());
			}

			return true;
		} else {
			libraryPath = new File(storedPath);
			return false;
		}
	}
//		if(libraryPath!=null){
//			return false;
//		}
//		//FIXME Should not use Preferences.NON-STATIC-METHODS
//		String storedPath = Preferences.getProperty(LIBRARY_KEY);
//		if ((storedPath == null) || !(new File(storedPath)).exists()) {
//			// File dir = new
//			// File(getClass().getClassLoader().getResource("PlugInMapsPlugInSelector.class").getFile());
//			if(!silentFailure) {
//				File dir = MipavController.getDefaultWorkingDirectory();
//				libraryPath = askForLibraryDirectory(dir.getParentFile());
//				if (libraryPath == null) {
//					libraryPath = dir.getParentFile();
//					if (!libraryPath.exists()) {
//						libraryPath.mkdir();
//					}
//				}
//				//FIXME Should not use Preferences.NON-STATIC-METHODS
//				Preferences.setProperty(LIBRARY_KEY, libraryPath.getAbsolutePath());
//				return true;
//			} else {
//				throw new RuntimeException("PipeLibrary: Cannot locate JIST library directory.");
//			}
//		} else {
//			libraryPath = new File(storedPath);
//			//FIXME Should not use Preferences.NON-STATIC-METHODS
//			Preferences.setProperty(LIBRARY_KEY, libraryPath.getAbsolutePath());
//			return false;
//		}
//	}

	/**
	 * Load existing modules and add them to the tree.
	 * 
	 * @return true, if load preferences
	 */
	public boolean loadPreferences()  {
		return loadPreferences(true);
	}

	/**
	 * Load existing modules and add them to the tree.
	 * 
	 * @return true, if load preferences
	 */
	public boolean loadPreferences(boolean silentFailure) {
		boolean ret = loadLibraryPath(silentFailure);
		JistPreferences.loadPreferences(silentFailure);
		return ret;
	}

	/**
	 * Remove non-group modules that can be rebuilt from existing script/pipe
	 * files or from the source tree.
	 * 
	 * @param monitor the monitor
	 */
	public void removeAtomicModules(ProgressMonitor monitor) {
		if(algos==null)return;
		monitor.setProgress(0);
		monitor.setMaximum(1);
		for (PipeAlgorithmFactory factory : algos) {
			Class cls = factory.getModuleClass();
			if ((cls != null) && !cls.equals(PipeAlgorithmGroup.class)) {
				if (cls.equals(MedicAlgorithmExecutableAdapter.class) || cls.equals(MedicAlgorithmMipavAdapter.class)
						 || cls.equals(MedicAlgorithmGeneralScriptAdapter.class)) {
					// These are Loni files that can always be rebuilt from pipe
					// files
					monitor.setNote("Deleting "+factory.getMapFile().getName());
					factory.getMapFile().delete();
				} else {
					PipeAlgorithm algo;
					try {
						// Delete all atomic modules. Don't worry! We can always
						// rebuild these automatically.
						monitor.setNote("Deleting "+factory.getMapFile().getName());
						factory.getMapFile().delete();
						algo = new PipeAlgorithm((ProcessingAlgorithm) cls.newInstance());
						// Write map back to directory
						algo.write(factory.getMapFile());

					} catch (InstantiationException e) {
						e.printStackTrace();
					} catch (IllegalAccessException e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	/**
	 * Set library directory.
	 * 
	 * @param libraryPath library path
	 */
	public void setLibraryPath(File libraryPath) {
		this.libraryPath = libraryPath;
		//FIXME Should not use Preferences.NON-STATIC-METHODS
		Preferences.setProperty(LIBRARY_KEY, libraryPath.getAbsolutePath());
	}
	
	
}
