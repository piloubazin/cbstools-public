package de.mpg.cbs.jist.cortex;

import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.AlgorithmAuthor;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation.Citation;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.image.*;
import edu.jhu.ece.iacl.algorithms.PrinceGroupAuthors;
import edu.jhu.ece.iacl.algorithms.ReferencedPapers;
import edu.jhu.ece.iacl.algorithms.ace.AnatomicallyConsistentEnhancement;

/**
 * ACE : Anatomically Consistent Enhancement Algorithm (with convenient option wrappers)
 * 
 * @author Blake Lucas
 * @author Pierre-Louis Bazin
 * 
 */
public class JistCortexACE extends ProcessingAlgorithm {
	
	private ParamVolume gmVol;
	private ParamVolume wmVol;
	private ParamVolume csfVol;

	/* use defaults of 126.9 
	private ParamFloat thresh;
	*/
	
	private ParamVolume thinSkelVol;
	private ParamVolume aceGmVol;
	private ParamVolume aceCsfVol;

	private static final String revnum = AnatomicallyConsistentEnhancement.getVersion();

	/**
	 * Create input parameters for ACE: gray matter,white matter and intensity
	 * threshold
	 */
	protected void createInputParameters(ParamCollection inputParams) {
		wmVol = new ParamVolume(VoxelType.FLOAT);
		wmVol.setName("White Matter Probability Image");
		gmVol = new ParamVolume(VoxelType.FLOAT);
		gmVol.setName("Gray Matter Probability Image");
		csfVol = new ParamVolume(VoxelType.FLOAT);
		csfVol.setName("Sulcal CSF Probability Image");
		
		inputParams.add(wmVol);
		inputParams.add(gmVol);
		inputParams.add(csfVol);
		
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Cortex Processing.prev");
		inputParams.setLabel("ACE CSF Enhancement");
		inputParams.setName("ACE_CSF_Enhancement");
	
		AlgorithmInformation info = getAlgorithmInformation();
		info.add(new Citation("C. Xu, X. Han, J.L. Prince, "
								+"Improving Cortical Surface Reconstruction Accuracy Using an Anatomically Consistent Gray Matter Representation, "
								+"NeuroImage Human Brain Mapping 2000 Meeting, Poster N. 581, NeuroImage Vol. 11, No. 5, May 2000."));
		info.add(PrinceGroupAuthors.chenyangXu);
		info.add(PrinceGroupAuthors.xiaoHan);
		info.add(PrinceGroupAuthors.blakeLucas);
		info.setDescription("ACE : Anatomically Consistent Enhancement (enhances deep sulcal CSF based on WM, GM and CSF memberships).");
		info.setLongDescription("ACE takes WM, GM and CSF memberships. Sometimes sulcal CSF is blurred by partial volume effect, because the thickness of CSF is often quite small. ACE tries to undo the partial volume effect and creates thin CSF to reveal sulcal banks.");
		info.setWebsite("http://www.iacl.ece.jhu.edu/");
		info.setVersion("2.0f");
		info.setEditable(false);
		info.setStatus(DevelopmentStatus.RC);
		//info.setHelpID("10002");
	}

	/**
	 * Create output parameter for ACE: enhanced gray matter, thinned csf,
	 * thinned gray matter
	 */
	protected void createOutputParameters(ParamCollection outputParams) {
		aceGmVol = new ParamVolume(VoxelType.FLOAT,-1,-1,-1,1);
		aceGmVol.setName("Enhanced GM Probability");
		
		aceCsfVol = new ParamVolume(VoxelType.FLOAT,-1,-1,-1,1);
		aceCsfVol.setName("Enhanced CSF Probability");

		thinSkelVol = new ParamVolume(VoxelType.INT,-1,-1,-1,1);
		thinSkelVol.setName("Thinned Skeleton");

		outputParams.setName("ace");
		outputParams.setLabel("ACE");
		outputParams.add(aceGmVol);
		outputParams.add(aceCsfVol);
		outputParams.add(thinSkelVol);
	}
	
	/**
	 * Execute ACE Algorithm
	 */
	protected void execute(CalculationMonitor monitor) {
		
		// pre-process: compute memberships and scale into [0,254]
		ImageDataFloat wm = new ImageDataFloat(wmVol.getImageData());
		float[][][] wmimg = wm.toArray3d();
		ImageDataFloat gm = new ImageDataFloat(gmVol.getImageData());
		float[][][] gmimg = gm.toArray3d();
		ImageDataFloat csf = new ImageDataFloat(csfVol.getImageData());
		float[][][] csfimg = csf.toArray3d();
		
		int nx = wm.getRows();
		int ny= wm.getCols();
		int nz = wm.getSlices();

		float[][][] sum = new float[nx][ny][nz];
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			// memberships
			sum[x][y][z] = 0.0f;
			sum[x][y][z] += wmimg[x][y][z];
			sum[x][y][z] += gmimg[x][y][z];
			sum[x][y][z] += csfimg[x][y][z];
			
			if (sum[x][y][z]>0.0001f) {
				wmimg[x][y][z] /= sum[x][y][z] ;
				gmimg[x][y][z] /= sum[x][y][z] ;
				csfimg[x][y][z] /= sum[x][y][z] ;
			}
			
			// scaling
			wmimg[x][y][z] *= 254;
			gmimg[x][y][z] *= 254;
			csfimg[x][y][z] *= 254;
		}
		
		// main algorithm
		AnatomicallyConsistentEnhancement ace = new AnatomicallyConsistentEnhancement();
		monitor.observe(ace);
		ace.solve(wm, gm, 126.9f);

		ImageDataFloat aceGm = (ImageDataFloat)ace.getEnhancedGM();

		// change and rescale
		float[][][] acegmimg = aceGm.toArray3d();
		for (int x=0;x<nx;x++) for (int y=0;y<ny;y++) for (int z=0;z<nz;z++) {
			csfimg[x][y][z] += gmimg[x][y][z]-acegmimg[x][y][z];
			
			acegmimg[x][y][z] *= sum[x][y][z]/254;
			csfimg[x][y][z] *= sum[x][y][z]/254;
		}

		aceGm.setHeader(gm.getHeader());
		aceGm.setName(gm.getName()+"_ace");
		aceGmVol.setValue(aceGm);

		ImageDataFloat aceCsf = new ImageDataFloat(csfimg);		
		aceCsf.setHeader(csf.getHeader());
		aceCsf.setName(csf.getName()+"_ace");
		aceCsfVol.setValue(aceCsf);
		
		ImageData thinSkel = ace.getThinnedSkeleton();
		thinSkel.setHeader(gm.getHeader());
		thinSkel.setName(gm.getName()+"_ace-skel");
		thinSkelVol.setValue(thinSkel);

	}

}
