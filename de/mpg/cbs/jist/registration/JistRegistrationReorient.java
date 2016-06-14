package de.mpg.cbs.jist.registration;

import edu.jhu.ece.iacl.jist.pipeline.AbstractCalculation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmInformation;
import edu.jhu.ece.iacl.jist.pipeline.AlgorithmRuntimeException;
import edu.jhu.ece.iacl.jist.pipeline.CalculationMonitor;
import edu.jhu.ece.iacl.jist.pipeline.DevelopmentStatus;
import edu.jhu.ece.iacl.jist.pipeline.ProcessingAlgorithm;
import edu.jhu.ece.iacl.jist.pipeline.parameter.*;
import edu.jhu.ece.iacl.jist.structures.image.ImageData;
import edu.jhu.ece.iacl.jist.structures.image.ImageDataMipavWrapper;
import edu.jhu.ece.iacl.jist.structures.image.ImageHeader;
import gov.nih.mipav.model.algorithms.AlgorithmTransform;
import gov.nih.mipav.model.file.FileInfoBase;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.structures.TransMatrix;
import Jama.Matrix;

public class JistRegistrationReorient extends ProcessingAlgorithm{
	protected ParamVolume template;
	protected ParamVolumeCollection source,result;
	protected ParamOption orientDicom;
	protected ParamOption orientX;
	protected ParamOption orientY;
	protected ParamOption orientZ;
	protected ParamOption userdefImageOrient;
	protected ParamOption interpolation;
	protected ParamOption resolution;
	protected ParamMatrix  transmatrix;
	
	private final String[] axisOrientStrs = new String[]{"UNKNOWN","R2L_TYPE","L2R_TYPE","P2A_TYPE","A2P_TYPE","I2S_TYPE","S2I_TYPE"};
	private final String[] imageOrientStrs = new String[]{"AXIAL","CORONAL","SAGITTAL","UNKNOWN"};
	
	protected void createInputParameters(ParamCollection inputParams) {
		String[] dicomStrs= {"Dicom axial", "Dicom coronal", "Dicom sagittal", "User defined"};
		String[] imageStrs= {"Axial", "Coronal", "Sagittal", "Unknown"};
		String[] orientStrs= {"Unknown", "Patient Right to Left", "Patient Left to Right", 
								"Patient Posterior to Anterior", "Patient Anterior to Posterior",
								"Patient Inferior to Superior", "Patient Superior to Inferior"};
		String[] resolStrs = {"Unchanged", "Finest cubic", "Coarsest cubic", "Same as template"};
		String[] interpStrs= {"Nearest Neighbor", "Trilinear", "Bspline 3rd order", "Bspline 4th order",
								"Cubic Lagrangian", "Quintic Lagrangian", "Heptic Lagrangian", "Windowed Sinc"};

		inputParams.add(source=new ParamVolumeCollection("Source"));
		inputParams.add(template=new ParamVolume("Template"));		
		template.setMandatory(false);
		inputParams.add(orientDicom=new ParamOption("New image orientation",dicomStrs));
		inputParams.add(orientX=new ParamOption("User defined X-axis orientation (image left to right)",orientStrs));
		inputParams.add(orientY=new ParamOption("User defined Y-axis orientation (image top to bottom)",orientStrs));
		inputParams.add(orientZ=new ParamOption("User defined Z-axis orientation (into the screen)",orientStrs));
		inputParams.add(userdefImageOrient=new ParamOption("User defined Image Orientation",imageStrs));
		inputParams.add(interpolation=new ParamOption("Interpolation",interpStrs));
		inputParams.add(resolution=new ParamOption("Resolution",resolStrs));
		inputParams.setName("reorient");
		inputParams.setLabel("Reorient Volume");
		inputParams.setPackage("CBS Tools");
		inputParams.setCategory("Registration");
		AlgorithmInformation info=getAlgorithmInformation();
		info.setVersion("3.0.8");
		info.setStatus(DevelopmentStatus.BETA);
		info.setDescription("Reorient a volume to a particular anatomical orientation.");
		info.setLongDescription("Images are placed in a cannonical orientation, with or without resampling (wrapper for the MIPAV utility).");
	}
	
	protected void createOutputParameters(ParamCollection outputParams) {
		outputParams.add(result=new ParamVolumeCollection("Reoriented Volume"));
		outputParams.add(transmatrix=new ParamMatrix("Reorientation Matrix",4,4));
	}
	
	protected void execute(CalculationMonitor monitor) 			throws AlgorithmRuntimeException {
		reorientWrapper algo_run = new reorientWrapper();
		monitor.observe(algo_run);
		algo_run.execute();
	}


	protected class reorientWrapper extends AbstractCalculation {
		public reorientWrapper() {
			setLabel("Reorient (MIPAV)");
		}


		public void execute() {

			
			for(ImageData sourceImg : source.getImageDataList()){
	
			ModelImage image = sourceImg.getModelImageCopy();
			ImageHeader hdr = sourceImg.getHeader();
			ImageHeader.AxisOrientation[] axor = hdr.getAxisOrientation();
			ImageHeader.ImageOrientation imgor = hdr.getImageOrientation();
			ModelImage templ = null;
			if(template.getImageData()!=null)
				templ = template.getImageData().getModelImageCopy();
			
			int i, j;
			boolean found;
			int newOrient;
			float ri[] = new float[3];
			int   ni[] = new int[3];
			float r0[] = new float[3];
			int   n0[] = new int[3];
			
			FileInfoBase fileInfo = (FileInfoBase)(image.getFileInfo()[0].clone());
			
			// set resampled resolutions, dimensions
			ri[0] = image.getFileInfo()[0].getResolutions()[0];
			ri[1] = image.getFileInfo()[0].getResolutions()[1];
			ri[2] = image.getFileInfo()[0].getResolutions()[2];
			
			ni[0] = image.getExtents()[0];
			ni[1] = image.getExtents()[1];
			ni[2] = image.getExtents()[2];
			
			float r[] = new float[3];
			int   n[] = new int[3];
			for (i = 0; i <= 2; i++) {
				r[i] = ri[i];
				n[i] = ni[i];
			}
			
			int or[] = new int[3];
			or[0] = image.getFileInfo()[0].getAxisOrientation()[0];
			or[1] = image.getFileInfo()[0].getAxisOrientation()[1];
			or[2] = image.getFileInfo()[0].getAxisOrientation()[2];
        
			
			int newOr[] = new int[3];
			int newOrientIndex = orientDicom.getIndex();
            if (newOrientIndex == 0) {  // DICOM AXIAL
                newOr[0] = FileInfoBase.ORI_R2L_TYPE;  
                newOr[1] = FileInfoBase.ORI_A2P_TYPE;
                newOr[2] = FileInfoBase.ORI_I2S_TYPE;
                
                axor[0]=ImageHeader.AxisOrientation.R2L_TYPE;
                axor[1]=ImageHeader.AxisOrientation.A2P_TYPE;
                axor[2]=ImageHeader.AxisOrientation.I2S_TYPE;
                imgor = ImageHeader.ImageOrientation.AXIAL;
            }
            else if (newOrientIndex == 1) {// DICOM CORONAL
                newOr[0] = FileInfoBase.ORI_R2L_TYPE;
                newOr[1] = FileInfoBase.ORI_S2I_TYPE;
                newOr[2] = FileInfoBase.ORI_A2P_TYPE;
                
                axor[0]=ImageHeader.AxisOrientation.R2L_TYPE;
                axor[1]=ImageHeader.AxisOrientation.S2I_TYPE;
                axor[2]=ImageHeader.AxisOrientation.A2P_TYPE;
                imgor = ImageHeader.ImageOrientation.CORONAL;
            }
            else if (newOrientIndex == 2) { // DICOM SAGITTAL
                newOr[0] = FileInfoBase.ORI_A2P_TYPE;
                newOr[1] = FileInfoBase.ORI_S2I_TYPE;
                newOr[2] = FileInfoBase.ORI_L2R_TYPE;
                
                axor[0]=ImageHeader.AxisOrientation.A2P_TYPE;
                axor[1]=ImageHeader.AxisOrientation.S2I_TYPE;
                axor[2]=ImageHeader.AxisOrientation.L2R_TYPE;
                imgor = ImageHeader.ImageOrientation.SAGITTAL;
            }
            else { // USER DEFINED
				newOr[0] = orientX.getIndex();
				newOr[1] = orientY.getIndex();
				newOr[2] = orientZ.getIndex();
				
				axor[0]=ImageHeader.AxisOrientation.valueOf(axisOrientStrs[orientX.getIndex()]);
				axor[1]=ImageHeader.AxisOrientation.valueOf(axisOrientStrs[orientY.getIndex()]);
				axor[2]=ImageHeader.AxisOrientation.valueOf(axisOrientStrs[orientZ.getIndex()]);
				imgor = ImageHeader.ImageOrientation.valueOf(imageOrientStrs[userdefImageOrient.getIndex()]);
			}
			int resolutionIndex = resolution.getIndex();
			String interpType = "";
			if (interpolation.getIndex()>=0)
				interpType = (String)interpolation.getValue();
		
			if (resolutionIndex == 1) {
				// Finest cubic
				float rn = Math.min(r[0],Math.min(r[1],r[2]));
				n[0] = (int)Math.ceil(n[0]*r[0]/rn);
				r[0] = rn;
				n[1] = (int)Math.ceil(n[1]*r[1]/rn);
				r[1] = rn;
				n[2] = (int)Math.ceil(n[2]*r[2]/rn);
				r[2] = rn;
			} else if (resolutionIndex == 2) {
				// Coarsest cubic
				float rn = Math.max(r[0],Math.max(r[1],r[2]));
				n[0] = (int)Math.ceil(n[0]*r[0]/rn);
				r[0] = rn;
				n[1] = (int)Math.ceil(n[1]*r[1]/rn);
				r[1] = rn;
				n[2] = (int)Math.ceil(n[2]*r[2]/rn);
				r[2] = rn;
			} else if (resolutionIndex == 3) {
				// Same as template
				r[0] = templ.getFileInfo()[0].getResolutions()[0];
				r[1] = templ.getFileInfo()[0].getResolutions()[1];
				r[2] = templ.getFileInfo()[0].getResolutions()[2];
				n[0] = templ.getExtents()[0];
				n[1] = templ.getExtents()[1];
				n[2] = templ.getExtents()[2];
			}
			
			double X[][] = new double[4][4];
			for (j = 0; j <= 2; j++) {
				switch (or[j]) {
					case FileInfoBase.ORI_R2L_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_R2L_TYPE) {
								
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_L2R_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
					case FileInfoBase.ORI_L2R_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_L2R_TYPE) {
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_R2L_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
					case FileInfoBase.ORI_A2P_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_A2P_TYPE) {
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_P2A_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
					case FileInfoBase.ORI_P2A_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_P2A_TYPE) {
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_A2P_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
					case FileInfoBase.ORI_I2S_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_I2S_TYPE) {
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_S2I_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
					case FileInfoBase.ORI_S2I_TYPE:
						found = false;
						for (i = 0; (i <= 2) && (!found); i++) {
							if (newOr[i] == FileInfoBase.ORI_S2I_TYPE) {
								found = true;
								X[i][j] = 1.0;
								r0[i] = r[j];
								n0[i] = n[j];
							}
							else if (newOr[i] == FileInfoBase.ORI_I2S_TYPE) {
								found = true;
								X[i][j] = -1.0;
								X[i][3] = ri[j]*(ni[j] - 1);
								r0[i] = r[j];
								n0[i] = n[j];
							}
						}
						break;
				}
			} // for (j = 0; j <= 2; j++)
			
			for (i = 0; i <= 2; i++) {
				fileInfo.setResolutions(r0[i], i);   
				fileInfo.setExtents(n0[i], i);
				fileInfo.setAxisOrientation(newOr[i], i);
			}
			
			if ((newOr[2] == FileInfoBase.ORI_I2S_TYPE) || (newOr[2] == FileInfoBase.ORI_S2I_TYPE)) {
				newOrient = FileInfoBase.AXIAL;
			}
			else if ((newOr[2] == FileInfoBase.ORI_A2P_TYPE) || (newOr[2] == FileInfoBase.ORI_P2A_TYPE)) {
				newOrient = FileInfoBase.CORONAL;
			}
			else if ((newOr[2] == FileInfoBase.ORI_L2R_TYPE) || (newOr[2] == FileInfoBase.ORI_R2L_TYPE)) {
				newOrient = FileInfoBase.SAGITTAL;
			}
			else {
				newOrient = FileInfoBase.UNKNOWN_ORIENT;
			} 
			fileInfo.setImageOrientation(newOrient);
			
			TransMatrix transform = new TransMatrix(4);
			transform.setMatrix(0, 2, 0, 3, X);
			
			System.out.println(transform.toString());
			
			int interp = AlgorithmTransform.TRILINEAR;
			if (interpType.equals("Nearest Neighbor")) {
				interp = AlgorithmTransform.NEAREST_NEIGHBOR;
			} else if (interpType.equals("Trilinear")) {
				interp = AlgorithmTransform.TRILINEAR;
			} else if (interpType.equals("Bspline 3rd order")) {
				interp = AlgorithmTransform.BSPLINE3;
			} else if (interpType.equals("Bspline 4th order")) {
				interp = AlgorithmTransform.BSPLINE4;
			} else if (interpType.equals("Cubic Lagrangian")) {
				interp = AlgorithmTransform.CUBIC_LAGRANGIAN;
			} else if (interpType.equals("Quintic Lagrangian")) {
				interp = AlgorithmTransform.QUINTIC_LAGRANGIAN;
			} else if (interpType.equals("Heptic Lagrangian")) {
				interp = AlgorithmTransform.HEPTIC_LAGRANGIAN;
			} else if  (interpType.equals("Windowed Sinc")) {
				interp = AlgorithmTransform.WSINC;
			}
				
			AlgorithmTransform algoTrans = new AlgorithmTransform(image, transform, interp, r0[0], r0[1], r0[2], n0[0], n0[1], n0[2], 
												true, false, false);
			algoTrans.setUpdateOriginFlag(true);
			
			algoTrans.run();
			
			ImageDataMipavWrapper img = new ImageDataMipavWrapper(algoTrans.getTransformedImage());
			hdr = new ImageHeader(fileInfo);
			//hdr.setAxisOrientation(axor);
			//hdr.setImageOrientation(imgor);
			
			img.setHeader(hdr);
			
			result.add(img);
			
			// renormalize the matrix X
			/*
			for (i=0;i<3;i++) for (j=0;j<3;j++) {
				if (X[i][j]<0) X[i][3] /= ri[j];
			}
			*/
			X[3][3] = 1.0;
			
			transmatrix.setValue(new Matrix(X));
			
			image.disposeLocal(); image = null;
			if (templ!=null) templ.disposeLocal(); templ = null;
			}
		}
	}

}
