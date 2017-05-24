package jist.modules;

public class JistDefaultBase implements JistModulesLocation {

	final static String[] default_paths = new String[]{
		"de.mpg.cbs.jist",
		"edu.jhu.ece.iacl.jist.plugins",
		"edu.jhu.ece.iacl.plugins"
		};
	
	@Override
	public String[] getValidModuleJavaPaths() {

		return default_paths;
	}
 
}
