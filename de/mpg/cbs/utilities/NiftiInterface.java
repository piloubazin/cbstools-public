// 
// Decompiled by Procyon v0.5.30
// 

package de.mpg.cbs.utilities;

import gov.nih.mipav.model.structures.MatrixHolder;
import gov.nih.mipav.model.file.FileInfoDicom;
import Jama.Matrix;
import java.awt.Frame;
import gov.nih.mipav.view.dialogs.JDialogNIFTIChoice;
import gov.nih.mipav.view.ViewUserInterface;
import gov.nih.mipav.model.file.FileUtility;
import gov.nih.mipav.model.file.FileWriteOptions;
import gov.nih.mipav.model.algorithms.utilities.AlgorithmChangeType;
import gov.nih.mipav.model.file.FileRaw;
import java.text.DecimalFormat;
import gov.nih.mipav.model.file.FileInfoBase;
import java.io.BufferedInputStream;
import gov.nih.mipav.view.Preferences;
import java.io.InputStream;
import java.io.IOException;
import java.io.FileNotFoundException;
import gov.nih.mipav.view.MipavUtil;
import java.io.RandomAccessFile;
import gov.nih.mipav.model.file.CBZip2InputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;
import java.io.FileInputStream;
import gov.nih.mipav.model.structures.TransMatrix;
import gov.nih.mipav.model.structures.ModelImage;
import gov.nih.mipav.model.file.FileInfoNIFTI;
import java.io.File;
import gov.nih.mipav.model.file.FileBase;

public class NiftiInterface extends FileBase
{
    private int[] axisOrientation;
    private int[] axisOrientation2;
    private byte[] bufferByte;
    private short coord_code;
    private String fileDir;
    private File fileHeader;
    private FileInfoNIFTI fileInfo;
    private String fileName;
    private int freq_dim;
    private int headerSize;
    private ModelImage image;
    private double imageMax;
    private double imageMin;
    private String intentName;
    private float[] LPSOrigin;
    private float[] LPSOrigin2;
    private TransMatrix matrix;
    private TransMatrix matrixTwoDim;
    private TransMatrix matrix2;
    private double newMax;
    private double newMin;
    private boolean oneFile;
    private float[] origin;
    private int phase_dim;
    private float[] pixdim;
    private float qfac;
    private short qform_code;
    public float qoffset_x;
    public float qoffset_y;
    public float qoffset_z;
    private float quatern_a;
    private float quatern_b;
    private float quatern_c;
    private float quatern_d;
    private double r00;
    private double r01;
    private double r02;
    private double r10;
    private double r11;
    private double r12;
    private double r20;
    private double r21;
    private double r22;
    private String patientOrientationString;
    private float[] resolutions;
    private float scl_slope;
    private float scl_inter;
    private short sform_code;
    private int slice_dim;
    private byte sliceCode;
    private float sliceDuration;
    private short sliceEnd;
    private short sliceStart;
    private short sourceBitPix;
    private short sourceType;
    private int spaceUnits;
    private float[] srow_x;
    private float[] srow_y;
    private float[] srow_z;
    private int timeUnits;
    private float tOffset;
    private float vox_offset;
    private int esize;
    private int ecode;
    private int[] esizeArray;
    private int[] ecodeArray;
    private String[] mindIdentArray;
    private float[] bValueArray;
    private float[] azimuthArray;
    private float[] zenithArray;
    private int[][] dtComponentArray;
    private int[] degreeArray;
    private int[] orderArray;
    private String[] afniGroupArray;
    private String[] asciiTextArray;
    private String[] caretArray;
    private File file;
    private FileInputStream fis;
    private ZipInputStream zin;
    private GZIPInputStream gzin;
    private CBZip2InputStream bz2in;
    
    public NiftiInterface(String fName, final String fDir) {
        this.bufferByte = null;
        this.fileInfo = null;
        this.freq_dim = 0;
        this.headerSize = 348;
        this.matrix = new TransMatrix(4);
        this.matrixTwoDim = new TransMatrix(3);
        this.matrix2 = null;
        this.phase_dim = 0;
        this.patientOrientationString = null;
        this.slice_dim = 0;
        this.spaceUnits = 0;
        this.timeUnits = 0;
        this.vox_offset = 0.0f;
        this.esize = 0;
        this.esizeArray = null;
        this.ecodeArray = null;
        this.mindIdentArray = null;
        this.bValueArray = null;
        this.azimuthArray = null;
        this.zenithArray = null;
        this.dtComponentArray = null;
        this.degreeArray = null;
        this.orderArray = null;
        this.afniGroupArray = null;
        this.asciiTextArray = null;
        this.caretArray = null;
        int index = fName.length();
        for (int i = fName.length() - 1; i >= 0; --i) {
            if (fName.charAt(i) == '.') {
                index = i;
                break;
            }
        }
        if (fName.substring(index).equalsIgnoreCase(".HDR")) {
            String fileDataName = fName.substring(0, index) + ".img";
            File fileData = new File(fDir + fileDataName);
            if (fileData.exists()) {
                fName = fileDataName;
            }
            else {
                fileDataName = fName.substring(0, index) + ".IMG";
                fileData = new File(fDir + fileDataName);
                if (fileData.exists()) {
                    fName = fileDataName;
                }
            }
        }
        this.fileName = fName;
        this.fileDir = fDir;
    }
    
    public static boolean isNIFTI(final String fName, final String fDir) {
        try {
            int index = fName.length();
            for (int i = fName.length() - 1; i >= 0; --i) {
                if (fName.charAt(i) == '.') {
                    index = i;
                    break;
                }
            }
            String fileHeaderName = fName.substring(0, index) + ".hdr";
            File fileHeader = new File(fDir + fileHeaderName);
            if (!fileHeader.exists()) {
                fileHeaderName = fName.substring(0, index) + ".HDR";
                fileHeader = new File(fDir + fileHeaderName);
                if (!fileHeader.exists()) {
                    fileHeaderName = fName.substring(0, index) + ".nii";
                    fileHeader = new File(fDir + fileHeaderName);
                    if (!fileHeader.exists()) {
                        return false;
                    }
                }
            }
            final RandomAccessFile raFile = new RandomAccessFile(fileHeader, "r");
            raFile.seek(344L);
            final byte[] b = new byte[4];
            raFile.read(b);
            final String tmpString = new String(b);
            raFile.close();
            return tmpString.equals("ni1\u0000") || tmpString.equals("n+1\u0000");
        }
        catch (FileNotFoundException e) {
            MipavUtil.displayError("FileNIFTI: Error reading file.");
        }
        catch (IOException e2) {
            MipavUtil.displayError("FileNIFTI: Error reading file.");
        }
        return false;
    }
    
    public void absoluteValue(final ModelImage image) throws IOException {
        try {
            float[] buffer = null;
            int bufferSize;
            if (image.getNDims() > 1) {
                bufferSize = image.getSliceSize();
            }
            else {
                bufferSize = image.getExtents()[0];
            }
            int nBuffers;
            if (image.getNDims() == 5) {
                nBuffers = image.getExtents()[4] * image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 4) {
                nBuffers = image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 3) {
                nBuffers = image.getExtents()[2];
            }
            else {
                nBuffers = 1;
            }
            if (!image.isColorImage()) {
                buffer = new float[bufferSize];
                final int xDim = image.getExtents()[0];
                final int yDim = image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    image.exportData(k * bufferSize, bufferSize, buffer);
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; ++i) {
                            buffer[j * xDim + i] = Math.abs(buffer[j * xDim + i]);
                        }
                    }
                    image.importData(k * bufferSize, buffer, false);
                }
            }
        }
        catch (IOException error) {
            throw new IOException("FileNIFTI.absoluteValue: " + error);
        }
        catch (OutOfMemoryError error2) {
            throw error2;
        }
    }
    
    public void flipTopBottom(final ModelImage image) throws IOException {
        try {
            float[] buffer = null;
            float[] resultBuffer = null;
            int bufferSize;
            if (image.getNDims() > 1) {
                bufferSize = image.getSliceSize();
            }
            else {
                bufferSize = image.getExtents()[0];
            }
            int nBuffers;
            if (image.getNDims() == 5) {
                nBuffers = image.getExtents()[4] * image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 4) {
                nBuffers = image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 3) {
                nBuffers = image.getExtents()[2];
            }
            else {
                nBuffers = 1;
            }
            if (image.isColorImage()) {
                buffer = new float[bufferSize * 4];
                resultBuffer = new float[bufferSize * 4];
                bufferSize *= 4;
                final int xDim = image.getExtents()[0] * 4;
                final int yDim = image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    image.exportData(k * bufferSize, bufferSize, buffer);
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; i += 4) {
                            resultBuffer[j * xDim + i] = 255.0f;
                            resultBuffer[j * xDim + i + 1] = buffer[(yDim - 1 - j) * xDim + i + 1];
                            resultBuffer[j * xDim + i + 2] = buffer[(yDim - 1 - j) * xDim + i + 2];
                            resultBuffer[j * xDim + i + 3] = buffer[(yDim - 1 - j) * xDim + i + 3];
                        }
                    }
                    image.importData(k * bufferSize, resultBuffer, false);
                }
            }
            else {
                buffer = new float[bufferSize];
                resultBuffer = new float[bufferSize];
                final int xDim = image.getExtents()[0];
                final int yDim = image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    image.exportData(k * bufferSize, bufferSize, buffer);
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; ++i) {
                            resultBuffer[j * xDim + i] = buffer[(yDim - 1 - j) * xDim + i];
                        }
                    }
                    image.importData(k * bufferSize, resultBuffer, false);
                }
            }
        }
        catch (IOException error) {
            throw new IOException("FileNIFTI.flipTopBottom: " + error);
        }
        catch (OutOfMemoryError error2) {
            throw error2;
        }
    }
    
    public void flipTopBottom(final float[] buffer, final FileInfoNIFTI fileInfo) throws IOException {
        float[] resultBuffer = null;
        try {
            int bufferSize;
            if (fileInfo.getExtents().length - 1 > 1) {
                bufferSize = fileInfo.getExtents()[0] * fileInfo.getExtents()[1];
            }
            else {
                bufferSize = fileInfo.getExtents()[0];
            }
            int nBuffers;
            if (fileInfo.getExtents().length - 1 == 5) {
                nBuffers = fileInfo.getExtents()[4] * fileInfo.getExtents()[3] * fileInfo.getExtents()[2];
            }
            else if (fileInfo.getExtents().length - 1 == 4) {
                nBuffers = fileInfo.getExtents()[3] * fileInfo.getExtents()[2];
            }
            else if (fileInfo.getExtents().length - 1 == 3) {
                nBuffers = fileInfo.getExtents()[2];
            }
            else {
                nBuffers = 1;
            }
            if (fileInfo.getDataType() == 9) {
                resultBuffer = new float[bufferSize * 4];
                bufferSize *= 4;
                final int xDim = this.image.getExtents()[0] * 4;
                final int yDim = this.image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; i += 4) {
                            resultBuffer[j * xDim + i] = 255.0f;
                            resultBuffer[k * bufferSize + j * xDim + i + 1] = buffer[k * bufferSize + (yDim - 1 - j) * xDim + i + 1];
                            resultBuffer[k * bufferSize + j * xDim + i + 2] = buffer[k * bufferSize + (yDim - 1 - j) * xDim + i + 2];
                            resultBuffer[k * bufferSize + j * xDim + i + 3] = buffer[k * bufferSize + (yDim - 1 - j) * xDim + i + 3];
                        }
                    }
                }
            }
            else {
                resultBuffer = new float[buffer.length];
                final int xDim = fileInfo.getExtents()[0];
                final int yDim = fileInfo.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; ++i) {
                            resultBuffer[k * bufferSize + j * xDim + i] = buffer[k * bufferSize + (yDim - 1 - j) * xDim + i];
                        }
                    }
                }
            }
        }
        catch (OutOfMemoryError error) {
            throw error;
        }
        System.arraycopy(resultBuffer, 0, buffer, 0, buffer.length);
    }
    
    public void finalize() {
        this.axisOrientation = null;
        this.axisOrientation2 = null;
        this.bufferByte = null;
        this.LPSOrigin = null;
        this.LPSOrigin2 = null;
        this.matrix = null;
        this.matrixTwoDim = null;
        this.matrix2 = null;
        this.origin = null;
        this.pixdim = null;
        this.resolutions = null;
        this.srow_x = null;
        this.srow_y = null;
        this.srow_z = null;
        this.intentName = null;
        this.fileDir = null;
        this.fileName = null;
        this.fileInfo = null;
        this.image = null;
        try {
            super.finalize();
        }
        catch (Throwable t) {}
    }
    
    public FileInfoNIFTI getFileInfo() {
        return this.fileInfo;
    }
    
    public TransMatrix getMatrix() {
        return this.matrix;
    }
    
    public TransMatrix getMatrix2() {
        return this.matrix2;
    }
    
    public boolean readHeader(final String imageFileName, final String fileDir, final boolean niftiCompressed) throws IOException {
        long fileLength = 348L;
        boolean endianess = NiftiInterface.BIG_ENDIAN;
        final int[] niftiExtents = new int[5];
        int numDims = 0;
        boolean isQform = true;
        byte extension0 = 0;
        int ecodeNumber = 0;
        int mindIdentNumber = 0;
        int mindIdentIndex = 0;
        int bValueNumber = 0;
        int bValueIndex = 0;
        int sphericalDirectionNumber = 0;
        int sphericalDirectionIndex = 0;
        int dtComponentNumber = 0;
        int dtComponentIndex = 0;
        int sphericalHarmonicNumber = 0;
        int sphericalHarmonicIndex = 0;
        int afniGroupNumber = 0;
        int afniGroupIndex = 0;
        int asciiTextNumber = 0;
        int asciiTextIndex = 0;
        int caretNumber = 0;
        int caretIndex = 0;
        this.bufferByte = new byte[this.headerSize];
        this.fileInfo = new FileInfoNIFTI(this.fileName, fileDir, 40);
        if (niftiCompressed) {
            this.oneFile = true;
            byte[] buffer = new byte[this.headerSize];
            byte[] buffer2 = null;
            int bytesRead2 = 0;
            this.file = new File(fileDir + File.separator + this.fileName);
            if (this.fileName.endsWith("zip")) {
                this.zin.getNextEntry();
                final int bytesRead3 = this.zin.read(buffer);
                if (bytesRead3 != this.headerSize) {
                    buffer = this.getFullBuffer(this.zin, buffer, bytesRead3, this.headerSize);
                }
                endianess = NiftiInterface.BIG_ENDIAN;
                if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.BIG_ENDIAN)) {
                    endianess = NiftiInterface.LITTLE_ENDIAN;
                    if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.LITTLE_ENDIAN)) {
                        Preferences.debug("FileNIFTI:readHeader NIFTI header length != 348.\n", 2);
                        return false;
                    }
                }
                else {
                    Preferences.debug("FileNIFTI:readHeader Endianess = Big endian.\n", 2);
                }
                this.fileInfo.setEndianess(endianess);
                this.fileInfo.setSizeOfHeader(this.headerSize);
                this.vox_offset = this.getBufferFloat(buffer, 108, endianess);
                this.fileInfo.setVoxOffset(this.vox_offset);
                Preferences.debug("vox_offset = " + this.vox_offset + "\n", 2);
                int dataStart = Math.round(this.vox_offset);
                if (dataStart > this.headerSize) {
                    buffer2 = new byte[dataStart - this.headerSize];
                    bytesRead2 = this.zin.read(buffer2);
                    if (bytesRead2 != dataStart - this.headerSize) {
                        buffer2 = this.getFullBuffer(this.zin, buffer2, bytesRead2, dataStart - this.headerSize);
                    }
                }
                else {
                    dataStart = this.headerSize;
                }
                this.bufferByte = new byte[dataStart];
                for (int i = 0; i < this.headerSize; ++i) {
                    this.bufferByte[i] = buffer[i];
                }
                for (int i = 0; i < dataStart - this.headerSize; ++i) {
                    this.bufferByte[i + this.headerSize] = buffer2[i];
                }
            }
            else if (this.fileName.endsWith("gz")) {
                try {
                    this.fis = new FileInputStream(this.file);
                }
                catch (FileNotFoundException e) {
                    MipavUtil.displayError("File not found exception on fis = new FileInputStream(file) for " + fileDir + File.separator + this.fileName);
                    return false;
                }
                try {
                    this.gzin = new GZIPInputStream(new BufferedInputStream(this.fis));
                }
                catch (IOException e2) {
                    MipavUtil.displayError("IOException on GZIPInputStream for " + this.fileName);
                    return false;
                }
                final int bytesRead3 = this.gzin.read(buffer);
                if (bytesRead3 != this.headerSize) {
                    buffer = this.getFullBuffer(this.gzin, buffer, bytesRead3, this.headerSize);
                }
                endianess = NiftiInterface.BIG_ENDIAN;
                if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.BIG_ENDIAN)) {
                    endianess = NiftiInterface.LITTLE_ENDIAN;
                    if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.LITTLE_ENDIAN)) {
                        Preferences.debug("FileNIFTI:readHeader NIFTI header length != 348.\n", 2);
                        return false;
                    }
                }
                else {
                    Preferences.debug("FileNIFTI:readHeader Endianess = Big endian.\n", 2);
                }
                this.fileInfo.setEndianess(endianess);
                this.fileInfo.setSizeOfHeader(this.headerSize);
                this.vox_offset = this.getBufferFloat(buffer, 108, endianess);
                this.fileInfo.setVoxOffset(this.vox_offset);
                Preferences.debug("vox_offset = " + this.vox_offset + "\n", 2);
                int dataStart = Math.round(this.vox_offset);
                if (dataStart > this.headerSize) {
                    buffer2 = new byte[dataStart - this.headerSize];
                    bytesRead2 = this.gzin.read(buffer2);
                    if (bytesRead2 != dataStart - this.headerSize) {
                        buffer2 = this.getFullBuffer(this.gzin, buffer2, bytesRead2, dataStart - this.headerSize);
                    }
                }
                else {
                    dataStart = this.headerSize;
                }
                this.bufferByte = new byte[dataStart];
                for (int i = 0; i < this.headerSize; ++i) {
                    this.bufferByte[i] = buffer[i];
                }
                for (int i = 0; i < dataStart - this.headerSize; ++i) {
                    this.bufferByte[i + this.headerSize] = buffer2[i];
                }
            }
            else if (this.fileName.endsWith("bz2")) {
                final int bytesRead3 = this.bz2in.read(buffer);
                if (bytesRead3 != this.headerSize) {
                    buffer = this.getFullBuffer((InputStream)this.bz2in, buffer, bytesRead3, this.headerSize);
                }
                endianess = NiftiInterface.BIG_ENDIAN;
                if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.BIG_ENDIAN)) {
                    endianess = NiftiInterface.LITTLE_ENDIAN;
                    if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.LITTLE_ENDIAN)) {
                        Preferences.debug("FileNIFTI:readHeader NIFTI header length != 348.\n", 2);
                        return false;
                    }
                }
                else {
                    Preferences.debug("FileNIFTI:readHeader Endianess = Big endian.\n", 2);
                }
                this.fileInfo.setEndianess(endianess);
                this.fileInfo.setSizeOfHeader(this.headerSize);
                this.vox_offset = this.getBufferFloat(buffer, 108, endianess);
                this.fileInfo.setVoxOffset(this.vox_offset);
                Preferences.debug("vox_offset = " + this.vox_offset + "\n", 2);
                int dataStart = Math.round(this.vox_offset);
                if (dataStart > this.headerSize) {
                    buffer2 = new byte[dataStart - this.headerSize];
                    bytesRead2 = this.bz2in.read(buffer2);
                    if (bytesRead2 != dataStart - this.headerSize) {
                        buffer2 = this.getFullBuffer((InputStream)this.bz2in, buffer2, bytesRead2, dataStart - this.headerSize);
                    }
                }
                else {
                    dataStart = this.headerSize;
                }
                this.bufferByte = new byte[dataStart];
                for (int i = 0; i < this.headerSize; ++i) {
                    this.bufferByte[i] = buffer[i];
                }
                for (int i = 0; i < dataStart - this.headerSize; ++i) {
                    this.bufferByte[i + this.headerSize] = buffer2[i];
                }
            }
        }
        else {
            final int index = this.fileName.lastIndexOf(".");
            if (index == -1) {
                this.oneFile = false;
                String fileHeaderName = this.fileName + ".hdr";
                this.fileHeader = new File(fileDir + fileHeaderName);
                if (!this.fileHeader.exists()) {
                    fileHeaderName = this.fileName + ".HDR";
                    this.fileHeader = new File(fileDir + fileHeaderName);
                }
                if (!this.fileHeader.exists()) {
                    return false;
                }
            }
            else {
                String fileHeaderName;
                if (this.fileName.substring(index + 1).equalsIgnoreCase("nii")) {
                    this.oneFile = true;
                    fileHeaderName = this.fileName;
                }
                else {
                    this.oneFile = false;
                    fileHeaderName = this.fileName.substring(0, index) + ".hdr";
                }
                this.fileHeader = new File(fileDir + fileHeaderName);
                if (!this.fileHeader.exists()) {
                    fileHeaderName = this.fileName.substring(0, index) + ".HDR";
                    this.fileHeader = new File(fileDir + fileHeaderName);
                    if (!this.fileHeader.exists()) {
                        return false;
                    }
                }
            }
            if (this.fileInfo == null) {
                this.fileInfo = new FileInfoNIFTI(imageFileName, fileDir, 40);
            }
            this.raFile = new RandomAccessFile(this.fileHeader, "r");
            final byte[] buffer = new byte[this.headerSize];
            this.raFile.read(buffer);
            endianess = NiftiInterface.BIG_ENDIAN;
            if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.BIG_ENDIAN)) {
                endianess = NiftiInterface.LITTLE_ENDIAN;
                if (this.headerSize != this.getBufferInt(buffer, 0, NiftiInterface.LITTLE_ENDIAN)) {
                    Preferences.debug("FileNIFTI:readHeader NIFTI header length != 348.\n", 2);
                    return false;
                }
            }
            else {
                Preferences.debug("FileNIFTI:readHeader Endianess = Big endian.\n", 2);
            }
            this.fileInfo.setEndianess(endianess);
            this.fileInfo.setSizeOfHeader(this.headerSize);
            this.vox_offset = this.getBufferFloat(buffer, 108, endianess);
            this.fileInfo.setVoxOffset(this.vox_offset);
            Preferences.debug("vox_offset = " + this.vox_offset + "\n", 2);
            int dataStart2 = Math.round(this.vox_offset);
            if (dataStart2 < this.headerSize) {
                dataStart2 = this.headerSize;
            }
            this.bufferByte = new byte[dataStart2];
            this.raFile.seek(0L);
            this.raFile.read(this.bufferByte);
            fileLength = this.raFile.length();
            Preferences.debug("\nThe size of the file with the header information = " + fileLength + "\n", 2);
        }
        switch (this.freq_dim = (this.bufferByte[39] & 0x3)) {
            case 0: {
                Preferences.debug("No frequency encoding direction is present\n", 2);
                break;
            }
            case 1: {
                Preferences.debug("Frequency encoding in the x direction\n", 2);
                break;
            }
            case 2: {
                Preferences.debug("Frequency encoding in the y direction\n", 2);
                break;
            }
            case 3: {
                Preferences.debug("Frequency encoding in the z direction\n", 2);
                break;
            }
        }
        this.fileInfo.setFreqDim(this.freq_dim);
        switch (this.phase_dim = (this.bufferByte[39] >> 2 & 0x3)) {
            case 0: {
                Preferences.debug("No phase encoding direction is present\n", 2);
                break;
            }
            case 1: {
                Preferences.debug("Phase encoding in the x direction\n", 2);
                break;
            }
            case 2: {
                Preferences.debug("Phase encoding in the y direction\n", 2);
                break;
            }
            case 3: {
                Preferences.debug("Phase encoding in the z direction\n", 2);
                break;
            }
        }
        this.fileInfo.setPhaseDim(this.phase_dim);
        switch (this.slice_dim = (this.bufferByte[39] >> 4 & 0x3)) {
            case 0: {
                Preferences.debug("No slice acquisition direction is present\n", 2);
                break;
            }
            case 1: {
                Preferences.debug("Slice acquisition in the x direction\n", 2);
                break;
            }
            case 2: {
                Preferences.debug("Slice acquisition in the y direction\n", 2);
                break;
            }
            case 3: {
                Preferences.debug("Slice acquisition in the z direction\n", 2);
                break;
            }
        }
        this.fileInfo.setSliceDim(this.slice_dim);
        final int dims = this.getBufferShort(this.bufferByte, 40, endianess);
        Preferences.debug("FileNIFTI:readHeader. Number of dimensions = " + dims + "\n", 2);
        for (int i = 0; i < dims; ++i) {
            niftiExtents[i] = this.getBufferShort(this.bufferByte, 42 + 2 * i, endianess);
            Preferences.debug("FileNIFTI:readHeader. Dimension " + (i + 1) + " = " + niftiExtents[i] + "\n", 2);
            if (niftiExtents[i] > 1) {
                ++numDims;
            }
        }
        int spatialDims = 0;
        for (int i = 0; i < Math.min(dims, 3); ++i) {
            if (niftiExtents[i] > 1) {
                ++spatialDims;
            }
        }
        final int[] extents = new int[numDims];
        int i = 0;
        int j = 0;
        while (i < dims) {
            if (niftiExtents[i] > 1) {
                extents[j++] = niftiExtents[i];
            }
            ++i;
        }
        this.fileInfo.setExtents(extents);
        final float intentP1 = this.getBufferFloat(this.bufferByte, 56, endianess);
        this.fileInfo.setIntentP1(intentP1);
        Preferences.debug("FileNIFTI:readHeader. intentP1 = " + this.fileInfo.getIntentP1() + "\n", 2);
        final float intentP2 = this.getBufferFloat(this.bufferByte, 60, endianess);
        this.fileInfo.setIntentP2(intentP2);
        Preferences.debug("FileNIFTI:readHeader. statPar2 = " + this.fileInfo.getIntentP2() + "\n", 2);
        final float intentP3 = this.getBufferFloat(this.bufferByte, 64, endianess);
        this.fileInfo.setIntentP3(intentP3);
        Preferences.debug("FileNIFTI:readHeader. intentP3 = " + this.fileInfo.getIntentP3() + "\n", 2);
        final short intentCode = this.getBufferShort(this.bufferByte, 68, endianess);
        this.fileInfo.setIntentCode(intentCode);
        Preferences.debug("FileNIFTI:readHeader. intentCode = " + intentCode + "\n", 2);
        switch (intentCode) {
            case 0: {
                Preferences.debug("No intention\n");
                break;
            }
            case 2: {
                Preferences.debug("Correlation coefficient R\n");
                Preferences.debug("Degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 2) {
                    Preferences.debug("Dimension " + numDims + " has the Correlation Coefficient R\n", 2);
                    Preferences.debug("in the first data plane and degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane\n", 2);
                    break;
                }
                break;
            }
            case 3: {
                Preferences.debug("Student t statistic\n", 2);
                Preferences.debug("Degress of freedom = " + Math.round(intentP1) + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 2) {
                    Preferences.debug("Dimension " + numDims + " has the Student t statistic\n", 2);
                    Preferences.debug("in the first data plane and degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane\n", 2);
                    break;
                }
                break;
            }
            case 4: {
                Preferences.debug("Fisher F statistic\n");
                Preferences.debug("Numerator degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                Preferences.debug("Denominator degrees of freedom = " + Math.round(intentP2) + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Fisher F statistic\n", 2);
                    Preferences.debug("in the first data plane, numerator degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane, and denominator degrees of freedom in the\n", 2);
                    Preferences.debug("third data plane\n", 2);
                    break;
                }
                break;
            }
            case 5: {
                Preferences.debug("Standard normal - N(0,1) distributed\n", 2);
                break;
            }
            case 6: {
                Preferences.debug("Chi - squared\n");
                Preferences.debug("Degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 2) {
                    Preferences.debug("Dimension " + numDims + " has Chi-squared\n", 2);
                    Preferences.debug("in the first data plane and degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane\n", 2);
                    break;
                }
                break;
            }
            case 7: {
                Preferences.debug("Beta distribution\n", 2);
                Preferences.debug("a parameter = " + intentP1 + "\n", 2);
                Preferences.debug("b parameter = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Beta distribution\n", 2);
                    Preferences.debug("in the first data plane, the a parameter in the\n", 2);
                    Preferences.debug("second data plane, and the b parameter in the third\n", 2);
                    Preferences.debug("third data plane\n", 2);
                    break;
                }
                break;
            }
            case 8: {
                Preferences.debug("Binomial distribution\n", 2);
                Preferences.debug("Number of trials = " + Math.round(intentP1) + "\n", 2);
                Preferences.debug("Probability per trial = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Binomial distribution\n", 2);
                    Preferences.debug("in the first data plane, the number of trials in the\n", 2);
                    Preferences.debug("second data plane, and the probability per trial in the\n", 2);
                    Preferences.debug("third data plane\n", 2);
                    break;
                }
                break;
            }
            case 9: {
                Preferences.debug("Gamma with PDF = x^(shape-1) * exp(-Scale*x)\n", 2);
                Preferences.debug("for x >= 0\n", 2);
                Preferences.debug("Shape = " + intentP1 + "\n", 2);
                Preferences.debug("Scale = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has Gamma\n", 2);
                    Preferences.debug("in the first data plane, shape in the\n", 2);
                    Preferences.debug("second data plane, and scale in the third\n", 2);
                    Preferences.debug("data plane\n", 2);
                    break;
                }
                break;
            }
            case 10: {
                Preferences.debug("Poisson distribution\n", 2);
                Preferences.debug("Mean = " + intentP1 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 2) {
                    Preferences.debug("Dimension " + numDims + " has the Poisson distribution\n", 2);
                    Preferences.debug("in the first data plane and the mean in the\n", 2);
                    Preferences.debug("second data plane\n", 2);
                    break;
                }
                break;
            }
            case 11: {
                Preferences.debug("Normal distribution\n", 2);
                Preferences.debug("Mean = " + intentP1 + "\n", 2);
                Preferences.debug("Standard deviation = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Normal distribution\n", 2);
                    Preferences.debug("in the first data plane, the mean in the\n", 2);
                    Preferences.debug("second data plane, and the standard deviation\n", 2);
                    Preferences.debug("in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 12: {
                Preferences.debug("Noncentral F statistic\n", 2);
                Preferences.debug("Numerator degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                Preferences.debug("Denominator degrees of freedom = " + Math.round(intentP2) + "\n", 2);
                Preferences.debug("Numerator noncentrality parameter= " + intentP3 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 4) {
                    Preferences.debug("Dimension " + numDims + " has the Noncentral F statistic\n", 2);
                    Preferences.debug("in the first data plane, numerator degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane, denominator degrees of freedom in the\n", 2);
                    Preferences.debug("third data plane, and the numerator noncentrality parameter\n", 2);
                    Preferences.debug("in the fourth data plane\n", 2);
                    break;
                }
                break;
            }
            case 13: {
                Preferences.debug("Noncentral chi-squared statistic\n", 2);
                Preferences.debug("Degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                Preferences.debug("Noncentrality parameter = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Noncentral chi-squared\n", 2);
                    Preferences.debug("statistic in the first data plane, degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane, and the noncentrality parameter in the\n", 2);
                    Preferences.debug("third data plane\n", 2);
                    break;
                }
                break;
            }
            case 14: {
                Preferences.debug("Logistic distribution\n", 2);
                Preferences.debug("Location = " + intentP1 + "\n", 2);
                Preferences.debug("Scale = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Logistic distribution\n", 2);
                    Preferences.debug("in the first data plane, location in the second\n", 2);
                    Preferences.debug("data plane, and scale in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 15: {
                Preferences.debug("Laplace distribution\n", 2);
                Preferences.debug("Location = " + intentP1 + "\n", 2);
                Preferences.debug("Scale = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Laplace distribution\n", 2);
                    Preferences.debug("in the first data plane, location in the second\n", 2);
                    Preferences.debug("data plane, and scale in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 16: {
                Preferences.debug("Uniform distribution\n", 2);
                Preferences.debug("Start = " + intentP1 + "\n", 2);
                Preferences.debug("End = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Uniform distribution\n", 2);
                    Preferences.debug("in the first data plane, start in the second data\n", 2);
                    Preferences.debug("plane, and end in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 17: {
                Preferences.debug("Noncentral t statistic\n", 2);
                Preferences.debug("Degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                Preferences.debug("Noncentrality parameter = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Noncentral t statistic\n", 2);
                    Preferences.debug("in the first data plane, degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane, and the noncentrality parameter in\n", 2);
                    Preferences.debug("third data plane\n", 2);
                    break;
                }
                break;
            }
            case 18: {
                Preferences.debug("Weibull distribution\n", 2);
                Preferences.debug("Location = " + intentP1 + "\n", 2);
                Preferences.debug("Scale = " + intentP2 + "\n", 2);
                Preferences.debug("Power = " + intentP3 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 4) {
                    Preferences.debug("Dimension " + numDims + " has the Weibull distribution\n", 2);
                    Preferences.debug("in the first data plane, location in the second\n", 2);
                    Preferences.debug("data plane, scale in the third data plane, and power\n", 2);
                    Preferences.debug("in the fourth data plane\n", 2);
                    break;
                }
                break;
            }
            case 19: {
                Preferences.debug("Chi distribution\n", 2);
                Preferences.debug("Degrees of freedom = " + Math.round(intentP1) + "\n", 2);
                final int p1 = Math.round(intentP1);
                if (p1 == 1) {
                    Preferences.debug("dof = 1 = half normal distribution\n", 2);
                }
                else if (p1 == 2) {
                    Preferences.debug("dof = 2 = Rayleigh distribution\n", 2);
                }
                else if (p1 == 3) {
                    Preferences.debug("dof = 3 = Maxwell-Boltzmann distribution\n", 2);
                }
                if (dims == 5 && extents[numDims - 1] == 2) {
                    Preferences.debug("Dimension " + numDims + " has the Chi distribution\n", 2);
                    Preferences.debug("in the first data plane and degrees of freedom in the\n", 2);
                    Preferences.debug("second data plane\n", 2);
                    break;
                }
                break;
            }
            case 20: {
                Preferences.debug("Inverse Gaussian\n", 2);
                Preferences.debug("Mu = " + intentP1 + "\n", 2);
                Preferences.debug("Lambda = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has the Inverse Gaussian\n", 2);
                    Preferences.debug("in the first data plane, mu in the second data\n", 2);
                    Preferences.debug("plane, and lambda in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 21: {
                Preferences.debug("Extreme value type 1\n", 2);
                Preferences.debug("Location = " + intentP1 + "\n", 2);
                Preferences.debug("Scale = " + intentP2 + "\n", 2);
                if (dims == 5 && extents[numDims - 1] == 3) {
                    Preferences.debug("Dimension " + numDims + " has Extreme value type 1\n", 2);
                    Preferences.debug("in the first data plane, location in the second\n", 2);
                    Preferences.debug("data plane, and scale in the third data plane\n", 2);
                    break;
                }
                break;
            }
            case 22: {
                Preferences.debug("Data is a p-value\n", 2);
                break;
            }
            case 23: {
                Preferences.debug("Data is ln(p-value)\n", 2);
                break;
            }
            case 24: {
                Preferences.debug("Data is log10(p-value)\n", 2);
                break;
            }
            case 1001: {
                Preferences.debug("Each voxel is an estimate of some parameter\n", 2);
                break;
            }
            case 1002: {
                Preferences.debug("Each voxel is an index into some set of labels\n", 2);
                break;
            }
            case 1003: {
                Preferences.debug("Each voxel is an index into the NeuroNames label set\n", 2);
                break;
            }
            case 1004: {
                Preferences.debug("Each voxel has a M x N matrix\n", 2);
                break;
            }
            case 1005: {
                Preferences.debug("Each voxel has a NxN symmetric matrix\n", 2);
                break;
            }
            case 1006: {
                Preferences.debug("Each voxel has a displacement vector\n", 2);
                break;
            }
            case 1007: {
                Preferences.debug("Each voxel has a vector\n", 2);
                break;
            }
            case 1008: {
                Preferences.debug("Each voxel has a spatial coordinate\n", 2);
                break;
            }
            case 1009: {
                Preferences.debug("Each voxel has a triple of indexes\n", 2);
                break;
            }
            case 1010: {
                Preferences.debug("Each voxel has a quarternion\n", 2);
                break;
            }
            case 1011: {
                Preferences.debug("Each voxel is a dimensionless value\n", 2);
                break;
            }
            default: {
                Preferences.debug("intentCode = " + intentCode + " is not a recognized value\n", 2);
                break;
            }
        }
        this.sourceType = this.getBufferShort(this.bufferByte, 70, endianess);
        this.fileInfo.setSourceType(this.sourceType);
        Preferences.debug("Original unscaled source data type:\n", 2);
        switch (this.sourceType) {
            case 0: {
                Preferences.debug("Unknown data type\n", 2);
                MipavUtil.displayError("Mipav cannot handle data type DT_UNKNOWN");
                return false;
            }
            case 1: {
                Preferences.debug("Binary data\n", 2);
                break;
            }
            case 256: {
                Preferences.debug("Signed byte data\n", 2);
                break;
            }
            case 2: {
                Preferences.debug("Unsigned byte data\n", 2);
                break;
            }
            case 4: {
                Preferences.debug("Signed short data\n", 2);
                break;
            }
            case 512: {
                Preferences.debug("Unsigned short data\n", 2);
                break;
            }
            case 8: {
                Preferences.debug("Signed integer data\n", 2);
                break;
            }
            case 768: {
                Preferences.debug("Unsigned integer data\n", 2);
                break;
            }
            case 1024: {
                Preferences.debug("Signed long data\n", 2);
                break;
            }
            case 1280: {
                Preferences.debug("Unsigned long data\n", 2);
                break;
            }
            case 16: {
                Preferences.debug("32 bit float data\n", 2);
                break;
            }
            case 64: {
                Preferences.debug("64 bit double data\n", 2);
                break;
            }
            case 1536: {
                Preferences.debug("128 bit float data\n");
                MipavUtil.displayError("MIPAV cannot handle 128 bit floating point data");
                return false;
            }
            case 128: {
                Preferences.debug("RGB 24 bit data\n", 2);
                break;
            }
            case 32: {
                Preferences.debug("64 bit complex data\n", 2);
                break;
            }
            case 1792: {
                Preferences.debug("128 bit DCOMPLEX data\n", 2);
                break;
            }
            case 2048: {
                Preferences.debug("256 bit complex data\n", 2);
                MipavUtil.displayError("MIPAV cannot handle 256 bit complex data");
                return false;
            }
            default: {
                Preferences.debug("Unknown datatype code = " + this.sourceType + "\n", 2);
                MipavUtil.displayError("Unknown datatype code = " + this.sourceType);
                return false;
            }
        }
        this.sourceBitPix = this.getBufferShort(this.bufferByte, 72, endianess);
        this.fileInfo.setSourceBitPix(this.sourceBitPix);
        Preferences.debug("FileNIFTI:readHeader. source bits per pixel = " + this.sourceBitPix + "\n", 2);
        this.sliceStart = this.getBufferShort(this.bufferByte, 74, endianess);
        this.pixdim = new float[dims + 1];
        this.resolutions = new float[Math.max(3, numDims)];
        i = 0;
        j = 0;
        while (i < dims + 1) {
            this.pixdim[i] = this.getBufferFloat(this.bufferByte, 76 + 4 * i, endianess);
            if (i >= 1 && niftiExtents[i - 1] > 1) {
                this.resolutions[j] = Math.abs(this.pixdim[i]);
                Preferences.debug("FileNIFTI:readHeader. Resolutions " + (j + 1) + " = " + this.resolutions[j] + "\n", 2);
                ++j;
            }
            ++i;
        }
        this.fileInfo.setResolutions(this.resolutions);
        this.scl_slope = this.getBufferFloat(this.bufferByte, 112, endianess);
        this.fileInfo.setSclSlope(this.scl_slope);
        Preferences.debug("Data scaling slope = " + this.scl_slope + "\n", 2);
        this.scl_inter = this.getBufferFloat(this.bufferByte, 116, endianess);
        this.fileInfo.setSclInter(this.scl_inter);
        Preferences.debug("Data offset = " + this.scl_inter + "\n", 2);
        this.sliceEnd = this.getBufferShort(this.bufferByte, 120, endianess);
        this.sliceCode = this.bufferByte[122];
        if (this.sliceCode > 0 && this.sliceStart > 0) {
            this.fileInfo.setSliceStart(this.sliceStart);
            Preferences.debug("Slice timing pattern starts with slice = " + (this.sliceStart + 1) + "\n", 2);
        }
        if (this.sliceCode > 0 && this.sliceEnd > this.sliceStart) {
            this.fileInfo.setSliceEnd(this.sliceEnd);
            Preferences.debug("Slice timing pattern ends with slice = " + (this.sliceEnd + 1) + "\n", 2);
        }
        int unitMeasure = 0;
        if (spatialDims == 0) {
            Preferences.debug("No x, y, or z dimensions are present\n", 2);
        }
        else {
            switch (this.spaceUnits = (this.bufferByte[123] & 0x7)) {
                case 0: {
                    Preferences.debug("Spatial units are unknown\n", 2);
                    unitMeasure = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
                    break;
                }
                case 1: {
                    Preferences.debug("Spatial units are meters\n", 2);
                    unitMeasure = FileInfoBase.Unit.METERS.getLegacyNum();
                    break;
                }
                case 2: {
                    Preferences.debug("Spatial units are millimeters\n", 2);
                    unitMeasure = FileInfoBase.Unit.MILLIMETERS.getLegacyNum();
                    break;
                }
                case 3: {
                    Preferences.debug("Spatial units are micrometers\n", 2);
                    unitMeasure = FileInfoBase.Unit.MICROMETERS.getLegacyNum();
                    break;
                }
                default: {
                    Preferences.debug("Spatial units are an illegal " + this.spaceUnits + "\n", 2);
                    unitMeasure = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
                    break;
                }
            }
            for (i = 0; i < spatialDims; ++i) {
                this.fileInfo.setUnitsOfMeasure(unitMeasure, i);
            }
        }
        if (dims >= 4 && niftiExtents[3] > 1) {
            switch (this.timeUnits = (this.bufferByte[123] & 0x38)) {
                case 0: {
                    Preferences.debug("Time units are unknown\n", 2);
                    unitMeasure = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
                    break;
                }
                case 8: {
                    Preferences.debug("Time units are seconds\n", 2);
                    unitMeasure = FileInfoBase.Unit.SECONDS.getLegacyNum();
                    break;
                }
                case 16: {
                    Preferences.debug("Time units are milliseconds\n", 2);
                    unitMeasure = FileInfoBase.Unit.MILLISEC.getLegacyNum();
                    break;
                }
                case 24: {
                    Preferences.debug("Time units are microseconds\n", 2);
                    unitMeasure = FileInfoBase.Unit.MICROSEC.getLegacyNum();
                    break;
                }
                case 32: {
                    Preferences.debug("Time units are hertz\n", 2);
                    unitMeasure = FileInfoBase.Unit.HZ.getLegacyNum();
                    break;
                }
                case 40: {
                    Preferences.debug("Time units are parts per million\n", 2);
                    unitMeasure = FileInfoBase.Unit.PPM.getLegacyNum();
                    break;
                }
                case 48: {
                    Preferences.debug("Time units are radians per second\n", 2);
                    unitMeasure = FileInfoBase.Unit.RADS.getLegacyNum();
                    break;
                }
                default: {
                    Preferences.debug("Time units are an illegal = " + this.timeUnits + "\n", 2);
                    unitMeasure = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
                    break;
                }
            }
            this.fileInfo.setUnitsOfMeasure(unitMeasure, spatialDims);
        }
        this.fileInfo.setCalMax(this.getBufferFloat(this.bufferByte, 124, endianess));
        this.fileInfo.setCalMin(this.getBufferFloat(this.bufferByte, 128, endianess));
        this.sliceDuration = this.getBufferFloat(this.bufferByte, 132, endianess);
        if (this.sliceDuration > 0.0f && this.slice_dim > 0) {
            this.fileInfo.setSliceDuration(this.sliceDuration);
            Preferences.debug("Time used to acquire 1 slice = " + this.sliceDuration + "\n", 2);
        }
        if (this.sliceCode > 0 && this.slice_dim > 0 && this.sliceDuration > 0.0f) {
            if (this.sliceCode == 1) {
                Preferences.debug("Slice timing order is sequentially increasing\n", 2);
            }
            else if (this.sliceCode == 2) {
                Preferences.debug("Slice timing order is sequentially decreasing\n", 2);
            }
            else if (this.sliceCode == 3) {
                Preferences.debug("Slice timing order is alternately increasing\n", 2);
            }
            else if (this.sliceCode == 4) {
                Preferences.debug("Slice timing order is alternately decreasing\n", 2);
            }
            else if (this.sliceCode == 5) {
                Preferences.debug("Slice timing order is alternately increasing #2\n", 2);
            }
            else if (this.sliceCode == 6) {
                Preferences.debug("Slice timing order is alternately decreasing #2\n", 2);
            }
            else {
                Preferences.debug("slice code has an illegal value = " + this.sliceCode + "\n", 2);
            }
        }
        else {
            Preferences.debug("Slice timing order is not specified\n", 2);
        }
        this.tOffset = this.getBufferFloat(this.bufferByte, 136, endianess);
        this.fileInfo.setOrigin(this.tOffset, 3);
        Preferences.debug("tOffset = " + this.tOffset + "\n", 2);
        switch (this.sourceType) {
            case 0: {
                return false;
            }
            case 1: {
                this.fileInfo.setDataType(0);
                this.fileInfo.setBitPix((short)1);
                break;
            }
            case 2: {
                this.fileInfo.setDataType(2);
                this.fileInfo.setBitPix((short)8);
                break;
            }
            case 4: {
                this.fileInfo.setDataType(3);
                this.fileInfo.setBitPix((short)16);
                break;
            }
            case 8: {
                this.fileInfo.setDataType(5);
                this.fileInfo.setBitPix((short)32);
                break;
            }
            case 16: {
                this.fileInfo.setDataType(7);
                this.fileInfo.setBitPix((short)32);
                break;
            }
            case 32: {
                this.fileInfo.setDataType(12);
                this.fileInfo.setBitPix((short)64);
                break;
            }
            case 64: {
                this.fileInfo.setDataType(8);
                this.fileInfo.setBitPix((short)64);
                break;
            }
            case 128: {
                this.fileInfo.setDataType(9);
                this.fileInfo.setBitPix((short)24);
                break;
            }
            case 256: {
                this.fileInfo.setDataType(1);
                this.fileInfo.setBitPix((short)8);
                break;
            }
            case 512: {
                this.fileInfo.setDataType(4);
                this.fileInfo.setBitPix((short)16);
                break;
            }
            case 768: {
                this.fileInfo.setDataType(14);
                this.fileInfo.setBitPix((short)32);
                break;
            }
            case 1024: {
                this.fileInfo.setDataType(6);
                this.fileInfo.setBitPix((short)64);
                break;
            }
            case 1280: {
                this.fileInfo.setDataType(6);
                this.fileInfo.setBitPix((short)64);
                break;
            }
            case 1792: {
                this.fileInfo.setDataType(13);
                this.fileInfo.setBitPix((short)128);
                break;
            }
            default: {
                return false;
            }
        }
        this.fileInfo.setDescription(new String(this.bufferByte, 148, 80));
        this.fileInfo.setModality(FileInfoBase.getModalityFromStr(this.fileInfo.getDescription()));
        this.fileInfo.setAuxFile(new String(this.bufferByte, 228, 24));
        this.qform_code = this.getBufferShort(this.bufferByte, 252, endianess);
        this.sform_code = this.getBufferShort(this.bufferByte, 254, endianess);
        if (this.pixdim[0] >= 0.0f) {
            this.qfac = 1.0f;
        }
        else {
            this.qfac = -1.0f;
        }
        if (this.qform_code == 0 && this.sform_code == 0) {
            this.fileInfo.setImageOrientation(3);
            (this.axisOrientation = new int[3])[0] = 0;
            this.axisOrientation[1] = 0;
            this.axisOrientation[2] = 0;
            this.fileInfo.setAxisOrientation(this.axisOrientation);
            this.matrix.set(0, 0, this.resolutions[0]);
            this.matrix.set(1, 1, this.resolutions[1]);
            this.matrix.set(2, 2, this.resolutions[2]);
        }
        Preferences.debug("qform_code = " + this.qform_code + "\n", 2);
        Preferences.debug("sform_code = " + this.sform_code + "\n", 2);
        if (this.qform_code > 0) {
            this.coord_code = this.qform_code;
            isQform = true;
        }
        else if (this.sform_code > 0) {
            this.coord_code = this.sform_code;
            isQform = false;
        }
        this.fileInfo.setCoordCode(this.coord_code);
        if (this.qform_code > 0 && this.sform_code > 0) {
            this.fileInfo.setCoordCode2(this.sform_code);
        }
        if (this.coord_code > 0) {
            this.matrix.setIsNIFTI(true);
            this.matrix.setIsQform(isQform);
            switch (this.coord_code) {
                case 0: {
                    Preferences.debug("Arbitrary X,Y,Z coordinate system\n", 2);
                    this.matrix.setTransformID(0);
                    break;
                }
                case 1: {
                    Preferences.debug("Scanner based anatomical coordinates\n", 2);
                    this.matrix.setTransformID(6);
                    break;
                }
                case 2: {
                    Preferences.debug("Coordinates aligned to another file's or to anatomical truth\n", 2);
                    this.matrix.setTransformID(2);
                    break;
                }
                case 3: {
                    Preferences.debug("Talairach X,Y,Z coordinate system\n", 2);
                    this.matrix.setTransformID(3);
                    break;
                }
                case 4: {
                    this.matrix.setTransformID(4);
                    Preferences.debug("MNI 152 normalized X,Y,Z coordinates\n", 2);
                    break;
                }
                default: {
                    this.matrix.setTransformID(0);
                    Preferences.debug("Unknown coord_code = " + this.coord_code + "\n", 2);
                    break;
                }
            }
        }
        if (this.qform_code > 0 && this.sform_code > 0) {
            (this.matrix2 = new TransMatrix(4)).setIsNIFTI(true);
            this.matrix2.setIsQform(false);
            switch (this.sform_code) {
                case 0: {
                    Preferences.debug("Matrix 2 arbitrary X,Y,Z coordinate system\n", 2);
                    this.matrix2.setTransformID(0);
                    break;
                }
                case 1: {
                    Preferences.debug("Matrix 2 scanner based anatomical coordinates\n", 2);
                    this.matrix2.setTransformID(6);
                    break;
                }
                case 2: {
                    Preferences.debug("Matrix 2 coordinates aligned to another file's or to anatomical truth\n", 2);
                    this.matrix2.setTransformID(2);
                    break;
                }
                case 3: {
                    Preferences.debug("Matrix 2 Talairach X,Y,Z coordinate system\n", 2);
                    this.matrix2.setTransformID(3);
                    break;
                }
                case 4: {
                    this.matrix2.setTransformID(4);
                    Preferences.debug("Matrix 2 MNI 152 normalized X,Y,Z coordinates\n", 2);
                    break;
                }
                default: {
                    this.matrix2.setTransformID(0);
                    Preferences.debug("Unknown sform_code = " + this.sform_code + "\n", 2);
                    break;
                }
            }
        }
        if (this.qform_code > 0) {
            this.quatern_b = this.getBufferFloat(this.bufferByte, 256, endianess);
            double b = this.quatern_b;
            this.quatern_c = this.getBufferFloat(this.bufferByte, 260, endianess);
            double c = this.quatern_c;
            this.quatern_d = this.getBufferFloat(this.bufferByte, 264, endianess);
            double d = this.quatern_d;
            double a = 1.0 - b * b - c * c - d * d;
            if (a < 1.0E-7) {
                a = 1.0 / Math.sqrt(b * b + c * c + d * d);
                b *= a;
                c *= a;
                d *= a;
                a = 0.0;
            }
            else {
                a = Math.sqrt(a);
            }
            this.r00 = a * a + b * b - c * c - d * d;
            this.matrix.set(0, 0, -this.r00 * this.resolutions[0]);
            this.r01 = 2.0 * (b * c - a * d);
            this.matrix.set(0, 1, -this.r01 * this.resolutions[1]);
            this.r02 = 2.0 * (b * d + a * c);
            this.matrix.set(0, 2, -this.r02 * this.qfac * this.resolutions[2]);
            this.r10 = 2.0 * (b * c + a * d);
            this.matrix.set(1, 0, -this.r10 * this.resolutions[0]);
            this.r11 = a * a + c * c - b * b - d * d;
            this.matrix.set(1, 1, -this.r11 * this.resolutions[1]);
            this.r12 = 2.0 * (c * d - a * b);
            this.r20 = 2.0 * (b * d - a * c);
            this.matrix.set(2, 0, this.r20 * this.resolutions[0]);
            this.r21 = 2.0 * (c * d + a * b);
            this.matrix.set(2, 1, this.r21 * this.resolutions[1]);
            this.r22 = a * a + d * d - c * c - b * b;
            this.matrix.set(2, 2, this.r22 * this.qfac * this.resolutions[2]);
            this.patientOrientationString = new String();
            final DecimalFormat nf = new DecimalFormat("##0.0000000");
            this.patientOrientationString = nf.format(-this.r00) + "\\" + nf.format(-this.r10) + "\\" + nf.format(this.r20) + "\\" + nf.format(-this.r01) + "\\" + nf.format(-this.r11) + "\\" + nf.format(this.r21);
            this.fileInfo.setPatientOrientationString(this.patientOrientationString);
            this.matrix.set(1, 2, -this.r12 * this.qfac * this.resolutions[2]);
            this.qoffset_x = this.getBufferFloat(this.bufferByte, 268, endianess);
            this.qoffset_y = this.getBufferFloat(this.bufferByte, 272, endianess);
            this.qoffset_z = this.getBufferFloat(this.bufferByte, 276, endianess);
            this.LPSOrigin = new float[3];
            this.axisOrientation = this.getAxisOrientation(this.matrix);
            Preferences.debug("axisOrientation = " + this.axisOrientation[0] + "  " + this.axisOrientation[1] + "  " + this.axisOrientation[2] + "\n", 2);
            this.fileInfo.setAxisOrientation(this.axisOrientation);
            for (j = 0; j < 3; ++j) {
                if (this.axisOrientation[j] == 1) {
                    this.LPSOrigin[j] = -Math.abs(this.qoffset_x);
                }
                else if (this.axisOrientation[j] == 2) {
                    this.LPSOrigin[j] = Math.abs(this.qoffset_x);
                }
                else if (this.axisOrientation[j] == 4) {
                    this.LPSOrigin[j] = -Math.abs(this.qoffset_y);
                }
                else if (this.axisOrientation[j] == 3) {
                    this.LPSOrigin[j] = Math.abs(this.qoffset_y);
                }
                else if (this.axisOrientation[j] == 5) {
                    this.LPSOrigin[j] = -Math.abs(this.qoffset_z);
                }
                else if (this.axisOrientation[j] == 6) {
                    this.LPSOrigin[j] = Math.abs(this.qoffset_z);
                }
            }
            this.fileInfo.setOrigin(this.LPSOrigin[0], 0);
            this.fileInfo.setOrigin(this.LPSOrigin[1], 1);
            this.fileInfo.setOrigin(this.LPSOrigin[2], 2);
            this.matrix.set(0, 3, (double)this.LPSOrigin[0]);
            this.matrix.set(1, 3, (double)this.LPSOrigin[1]);
            this.matrix.set(2, 3, (double)this.LPSOrigin[2]);
            if (this.axisOrientation[2] == 1 || this.axisOrientation[2] == 2) {
                this.fileInfo.setImageOrientation(2);
            }
            else if (this.axisOrientation[2] == 4 || this.axisOrientation[2] == 3) {
                this.fileInfo.setImageOrientation(1);
            }
            else {
                this.fileInfo.setImageOrientation(0);
            }
            Preferences.debug("matrix = \n" + this.matrix + "\n", 2);
            Preferences.debug("quatern_a = " + this.quatern_a + "\n", 2);
            Preferences.debug("quatern_b = " + this.quatern_b + "\n", 2);
            Preferences.debug("quatern_c = " + this.quatern_c + "\n", 2);
            Preferences.debug("quatern_d = " + this.quatern_d + "\n", 2);
            Preferences.debug("qoffset_x = " + this.qoffset_x + "\n", 2);
            Preferences.debug("qoffset_y = " + this.qoffset_y + "\n", 2);
            Preferences.debug("qoffset_z = " + this.qoffset_z + "\n", 2);
            Preferences.debug("qfac = " + this.qfac + "\n", 2);
            Preferences.debug("r00 = " + this.r00 + "\n", 2);
            Preferences.debug("r01 = " + this.r01 + "\n", 2);
            Preferences.debug("r02 = " + this.r02 + "\n", 2);
            Preferences.debug("r10 = " + this.r10 + "\n", 2);
            Preferences.debug("r11 = " + this.r11 + "\n", 2);
            Preferences.debug("r12 = " + this.r12 + "\n", 2);
            Preferences.debug("r20 = " + this.r20 + "\n", 2);
            Preferences.debug("r21 = " + this.r21 + "\n", 2);
            Preferences.debug("r22 = " + this.r22 + "\n", 2);
        }
        else if (this.sform_code > 0) {
            this.srow_x = new float[4];
            this.srow_y = new float[4];
            this.srow_z = new float[4];
            this.srow_x[0] = this.getBufferFloat(this.bufferByte, 280, endianess);
            this.srow_x[1] = this.getBufferFloat(this.bufferByte, 284, endianess);
            this.srow_x[2] = this.getBufferFloat(this.bufferByte, 288, endianess);
            this.srow_x[3] = this.getBufferFloat(this.bufferByte, 292, endianess);
            this.srow_y[0] = this.getBufferFloat(this.bufferByte, 296, endianess);
            this.srow_y[1] = this.getBufferFloat(this.bufferByte, 300, endianess);
            this.srow_y[2] = this.getBufferFloat(this.bufferByte, 304, endianess);
            this.srow_y[3] = this.getBufferFloat(this.bufferByte, 308, endianess);
            this.srow_z[0] = this.getBufferFloat(this.bufferByte, 312, endianess);
            this.srow_z[1] = this.getBufferFloat(this.bufferByte, 316, endianess);
            this.srow_z[2] = this.getBufferFloat(this.bufferByte, 320, endianess);
            this.srow_z[3] = this.getBufferFloat(this.bufferByte, 324, endianess);
            this.matrix.set(0, 0, (double)(-this.srow_x[0]));
            this.matrix.set(0, 1, (double)(-this.srow_x[1]));
            this.matrix.set(0, 2, (double)(-this.srow_x[2]));
            this.matrix.set(1, 0, (double)(-this.srow_y[0]));
            this.matrix.set(1, 1, (double)(-this.srow_y[1]));
            this.matrix.set(1, 2, (double)(-this.srow_y[2]));
            this.matrix.set(2, 0, (double)this.srow_z[0]);
            this.matrix.set(2, 1, (double)this.srow_z[1]);
            this.matrix.set(2, 2, (double)this.srow_z[2]);
            this.r00 = -this.matrix.get(0, 0) / this.resolutions[0];
            this.r10 = -this.matrix.get(1, 0) / this.resolutions[0];
            this.r20 = this.matrix.get(2, 0) / this.resolutions[0];
            this.r01 = -this.matrix.get(0, 1) / this.resolutions[1];
            this.r11 = -this.matrix.get(1, 1) / this.resolutions[1];
            this.r21 = this.matrix.get(2, 1) / this.resolutions[1];
            this.patientOrientationString = new String();
            final DecimalFormat nf = new DecimalFormat("##0.0000000");
            this.patientOrientationString = nf.format(-this.r00) + "\\" + nf.format(-this.r10) + "\\" + nf.format(this.r20) + "\\" + nf.format(-this.r01) + "\\" + nf.format(-this.r11) + "\\" + nf.format(this.r21);
            this.fileInfo.setPatientOrientationString(this.patientOrientationString);
            this.LPSOrigin = new float[3];
            this.axisOrientation = this.getAxisOrientation(this.matrix);
            Preferences.debug("axisOrientation = " + this.axisOrientation[0] + "  " + this.axisOrientation[1] + "  " + this.axisOrientation[2] + "\n", 2);
            this.fileInfo.setAxisOrientation(this.axisOrientation);
            for (j = 0; j < 3; ++j) {
                if (this.axisOrientation[j] == 1) {
                    this.LPSOrigin[j] = -Math.abs(this.srow_x[3]);
                }
                else if (this.axisOrientation[j] == 2) {
                    this.LPSOrigin[j] = Math.abs(this.srow_x[3]);
                }
                else if (this.axisOrientation[j] == 4) {
                    this.LPSOrigin[j] = -Math.abs(this.srow_y[3]);
                }
                else if (this.axisOrientation[j] == 3) {
                    this.LPSOrigin[j] = Math.abs(this.srow_y[3]);
                }
                else if (this.axisOrientation[j] == 5) {
                    this.LPSOrigin[j] = -Math.abs(this.srow_z[3]);
                }
                else if (this.axisOrientation[j] == 6) {
                    this.LPSOrigin[j] = Math.abs(this.srow_z[3]);
                }
            }
            this.fileInfo.setOrigin(this.LPSOrigin[0], 0);
            this.fileInfo.setOrigin(this.LPSOrigin[1], 1);
            this.fileInfo.setOrigin(this.LPSOrigin[2], 2);
            this.matrix.set(0, 3, (double)this.LPSOrigin[0]);
            this.matrix.set(1, 3, (double)this.LPSOrigin[1]);
            this.matrix.set(2, 3, (double)this.LPSOrigin[2]);
            if (this.axisOrientation[2] == 1 || this.axisOrientation[2] == 2) {
                this.fileInfo.setImageOrientation(2);
            }
            else if (this.axisOrientation[2] == 4 || this.axisOrientation[2] == 3) {
                this.fileInfo.setImageOrientation(1);
            }
            else {
                this.fileInfo.setImageOrientation(0);
            }
            Preferences.debug("matrix = \n" + this.matrix + "\n", 2);
            Preferences.debug("srow_x = " + this.srow_x[0] + "  " + this.srow_x[1] + "  " + this.srow_x[2] + "  " + this.srow_x[3] + "\n", 2);
            Preferences.debug("srow_y = " + this.srow_y[0] + "  " + this.srow_y[1] + "  " + this.srow_y[2] + "  " + this.srow_y[3] + "\n", 2);
            Preferences.debug("srow_z = " + this.srow_z[0] + "  " + this.srow_z[1] + "  " + this.srow_z[2] + "  " + this.srow_z[3] + "\n", 2);
        }
        if (this.matrix2 != null) {
            this.srow_x = new float[4];
            this.srow_y = new float[4];
            this.srow_z = new float[4];
            this.srow_x[0] = this.getBufferFloat(this.bufferByte, 280, endianess);
            this.srow_x[1] = this.getBufferFloat(this.bufferByte, 284, endianess);
            this.srow_x[2] = this.getBufferFloat(this.bufferByte, 288, endianess);
            this.srow_x[3] = this.getBufferFloat(this.bufferByte, 292, endianess);
            this.srow_y[0] = this.getBufferFloat(this.bufferByte, 296, endianess);
            this.srow_y[1] = this.getBufferFloat(this.bufferByte, 300, endianess);
            this.srow_y[2] = this.getBufferFloat(this.bufferByte, 304, endianess);
            this.srow_y[3] = this.getBufferFloat(this.bufferByte, 308, endianess);
            this.srow_z[0] = this.getBufferFloat(this.bufferByte, 312, endianess);
            this.srow_z[1] = this.getBufferFloat(this.bufferByte, 316, endianess);
            this.srow_z[2] = this.getBufferFloat(this.bufferByte, 320, endianess);
            this.srow_z[3] = this.getBufferFloat(this.bufferByte, 324, endianess);
            this.matrix2.set(0, 0, (double)(-this.srow_x[0]));
            this.matrix2.set(0, 1, (double)(-this.srow_x[1]));
            this.matrix2.set(0, 2, (double)(-this.srow_x[2]));
            this.matrix2.set(1, 0, (double)(-this.srow_y[0]));
            this.matrix2.set(1, 1, (double)(-this.srow_y[1]));
            this.matrix2.set(1, 2, (double)(-this.srow_y[2]));
            this.matrix2.set(2, 0, (double)this.srow_z[0]);
            this.matrix2.set(2, 1, (double)this.srow_z[1]);
            this.matrix2.set(2, 2, (double)this.srow_z[2]);
            this.LPSOrigin2 = new float[3];
            this.axisOrientation2 = this.getAxisOrientation(this.matrix2);
            Preferences.debug("axisOrientation2 = " + this.axisOrientation2[0] + "  " + this.axisOrientation2[1] + "  " + this.axisOrientation2[2] + "\n", 2);
            for (j = 0; j < 3; ++j) {
                if (this.axisOrientation2[j] == 1) {
                    this.LPSOrigin2[j] = -Math.abs(this.srow_x[3]);
                }
                else if (this.axisOrientation2[j] == 2) {
                    this.LPSOrigin2[j] = Math.abs(this.srow_x[3]);
                }
                else if (this.axisOrientation2[j] == 4) {
                    this.LPSOrigin2[j] = -Math.abs(this.srow_y[3]);
                }
                else if (this.axisOrientation2[j] == 3) {
                    this.LPSOrigin2[j] = Math.abs(this.srow_y[3]);
                }
                else if (this.axisOrientation2[j] == 5) {
                    this.LPSOrigin2[j] = -Math.abs(this.srow_z[3]);
                }
                else if (this.axisOrientation2[j] == 6) {
                    this.LPSOrigin2[j] = Math.abs(this.srow_z[3]);
                }
            }
            this.matrix2.set(0, 3, (double)this.LPSOrigin2[0]);
            this.matrix2.set(1, 3, (double)this.LPSOrigin2[1]);
            this.matrix2.set(2, 3, (double)this.LPSOrigin2[2]);
            Preferences.debug("matrix2 = \n" + this.matrix2 + "\n", 2);
            Preferences.debug("srow_x = " + this.srow_x[0] + "  " + this.srow_x[1] + "  " + this.srow_x[2] + "  " + this.srow_x[3] + "\n", 2);
            Preferences.debug("srow_y = " + this.srow_y[0] + "  " + this.srow_y[1] + "  " + this.srow_y[2] + "  " + this.srow_y[3] + "\n", 2);
            Preferences.debug("srow_z = " + this.srow_z[0] + "  " + this.srow_z[1] + "  " + this.srow_z[2] + "  " + this.srow_z[3] + "\n", 2);
        }
        if (numDims == 2) {
            if (this.axisOrientation[0] == 2 || this.axisOrientation[0] == 3 || this.axisOrientation[0] == 6) {
                this.matrixTwoDim.set(0, 0, (double)(-this.resolutions[0]));
            }
            else {
                this.matrixTwoDim.set(0, 0, (double)this.resolutions[0]);
            }
            if (this.axisOrientation[1] == 2 || this.axisOrientation[1] == 3 || this.axisOrientation[1] == 6) {
                this.matrixTwoDim.set(1, 1, (double)(-this.resolutions[1]));
            }
            else {
                this.matrixTwoDim.set(1, 1, (double)this.resolutions[1]);
            }
            if (this.LPSOrigin != null) {
                this.matrixTwoDim.set(0, 2, (double)this.LPSOrigin[0]);
                this.matrixTwoDim.set(1, 2, (double)this.LPSOrigin[1]);
            }
        }
        if (this.qform_code > 0 && this.sform_code > 0) {
            this.fileInfo.setMatrixQ(this.matrix);
            this.fileInfo.setMatrixS(this.matrix2);
        }
        else if (this.qform_code > 0) {
            this.fileInfo.setMatrixQ(this.matrix);
        }
        else if (this.sform_code > 0) {
            this.fileInfo.setMatrixS(this.matrix);
        }
        this.intentName = new String(this.bufferByte, 328, 16);
        Preferences.debug("Name or meaning of data = " + this.intentName + "\n", 2);
        this.fileInfo.setIntentName(this.intentName.trim());
        if (this.bufferByte.length > 348 && this.vox_offset > 352.0f) {
            extension0 = this.bufferByte[348];
            Preferences.debug("First byte in extension array = " + extension0 + "\n", 2);
        }
        if (extension0 == 0) {
            Preferences.debug("No extended header information is present\n", 2);
        }
        else {
            Preferences.debug("This indicates a header extension follows the extension array\n", 2);
            int currentAddress;
            final int extendedHeaderStart = currentAddress = 352;
            this.esize = 8;
            while (this.bufferByte.length >= currentAddress + this.esize && (!this.oneFile || this.vox_offset >= currentAddress + this.esize)) {
                this.esize = this.getBufferInt(this.bufferByte, currentAddress, endianess);
                this.ecode = this.getBufferInt(this.bufferByte, currentAddress + 4, endianess);
                if (this.ecode == 4) {
                    ++afniGroupNumber;
                }
                else if (this.ecode == 6) {
                    ++asciiTextNumber;
                }
                else if (this.ecode == 18) {
                    ++mindIdentNumber;
                }
                else if (this.ecode == 20) {
                    ++bValueNumber;
                }
                else if (this.ecode == 22) {
                    ++sphericalDirectionNumber;
                }
                else if (this.ecode == 24) {
                    ++dtComponentNumber;
                }
                else if (this.ecode == 26) {
                    ++sphericalHarmonicNumber;
                }
                else if (this.ecode == 30) {
                    ++caretNumber;
                }
                ++ecodeNumber;
                currentAddress += this.esize;
            }
            Preferences.debug("The number of header fields = " + ecodeNumber + "\n", 2);
            if (ecodeNumber >= 1) {
                (this.esizeArray = new int[ecodeNumber])[0] = 8;
                this.ecodeArray = new int[ecodeNumber];
                this.mindIdentArray = new String[mindIdentNumber];
                this.bValueArray = new float[bValueNumber];
                this.azimuthArray = new float[sphericalDirectionNumber];
                this.zenithArray = new float[sphericalDirectionNumber];
                this.dtComponentArray = new int[dtComponentNumber][];
                this.degreeArray = new int[sphericalHarmonicNumber];
                this.orderArray = new int[sphericalHarmonicNumber];
                this.afniGroupArray = new String[afniGroupNumber];
                this.asciiTextArray = new String[asciiTextNumber];
                this.caretArray = new String[caretNumber];
                for (currentAddress = extendedHeaderStart, ecodeNumber = 0; this.bufferByte.length >= currentAddress + this.esizeArray[Math.max(0, ecodeNumber - 1)] && (!this.oneFile || this.vox_offset >= currentAddress + this.esizeArray[Math.max(0, ecodeNumber - 1)]); currentAddress += this.esizeArray[ecodeNumber], ++ecodeNumber) {
                    this.esizeArray[ecodeNumber] = this.getBufferInt(this.bufferByte, currentAddress, endianess);
                    Preferences.debug("Header field number " + (ecodeNumber + 1) + " size in bytes = " + this.esizeArray[ecodeNumber] + "\n", 2);
                    this.ecodeArray[ecodeNumber] = this.getBufferInt(this.bufferByte, currentAddress + 4, endianess);
                    Preferences.debug("Header field number " + (ecodeNumber + 1) + " has ", 2);
                    switch (this.ecodeArray[ecodeNumber]) {
                        case 0: {
                            Preferences.debug("ecode = 0 for an unknown private format\n", 2);
                            break;
                        }
                        case 2: {
                            Preferences.debug("ecode = 2 for DICOM format (i.e., attribute tags and values)\n", 2);
                            break;
                        }
                        case 4: {
                            Preferences.debug("ecode = 4 for AFNI group (i.e., ASCII XML-ish elements)\n", 2);
                            this.afniGroupArray[afniGroupIndex] = new String(this.bufferByte, currentAddress + 8, this.esizeArray[ecodeNumber] - 8);
                            ++afniGroupIndex;
                            break;
                        }
                        case 6: {
                            Preferences.debug("ecode = 6 for comment: arbitrary non-NUL ASCII text\n", 2);
                            this.asciiTextArray[asciiTextIndex] = new String(this.bufferByte, currentAddress + 8, this.esizeArray[ecodeNumber] - 8);
                            ++asciiTextIndex;
                            break;
                        }
                        case 8: {
                            Preferences.debug("ecode = 8 for XCEDE metadata\n", 2);
                            break;
                        }
                        case 10: {
                            Preferences.debug("ecode = 10 for dimensional information for JIM software(XML format)\n", 2);
                            break;
                        }
                        case 12: {
                            Preferences.debug("ecode = 12 for Fiswidget XML pipeline descriptions\n", 2);
                            break;
                        }
                        case 18: {
                            Preferences.debug("ecode = 18 for MIND_IDENT field with character data\n", 2);
                            this.mindIdentArray[mindIdentIndex] = new String(this.bufferByte, currentAddress + 8, this.esizeArray[ecodeNumber] - 8);
                            Preferences.debug("Mind Ident field = " + this.mindIdentArray[mindIdentIndex] + "\n", 2);
                            ++mindIdentIndex;
                            break;
                        }
                        case 20: {
                            Preferences.debug("ecode = 20 for B_VALUE for b-value in units of s/mm-squared\n", 2);
                            this.bValueArray[bValueIndex] = this.getBufferFloat(this.bufferByte, currentAddress + 8, endianess);
                            Preferences.debug("b-value = " + this.bValueArray[bValueIndex] + " s/(mm*mm)\n", 2);
                            ++bValueIndex;
                            break;
                        }
                        case 22: {
                            Preferences.debug("ecode = 22 for SPHERICAL_DIRECTION with spherical coordinates\n", 2);
                            this.azimuthArray[sphericalDirectionIndex] = this.getBufferFloat(this.bufferByte, currentAddress + 8, endianess);
                            Preferences.debug("Azimuthal angle = " + this.azimuthArray[sphericalDirectionIndex] + " radians\n", 2);
                            this.zenithArray[sphericalDirectionIndex] = this.getBufferFloat(this.bufferByte, currentAddress + 12, endianess);
                            Preferences.debug("Zenith angle = " + this.zenithArray[sphericalDirectionIndex] + " radians\n", 2);
                            ++sphericalDirectionIndex;
                            break;
                        }
                        case 24: {
                            Preferences.debug("ecode = 24 for DT_COMPONENT specifying the indicies of a single diffusion tensor component\n", 2);
                            final int dtComponents = (this.esizeArray[ecodeNumber] - 8) / 4;
                            this.dtComponentArray[dtComponentIndex] = new int[dtComponents];
                            for (i = 0; i < dtComponents; ++i) {
                                this.dtComponentArray[dtComponentIndex][i] = this.getBufferInt(this.bufferByte, currentAddress + 8 + 4 * i, endianess);
                                Preferences.debug("DT component index " + (i + 1) + " = " + this.dtComponentArray[dtComponentIndex][i] + "\n", 2);
                            }
                            ++dtComponentIndex;
                            break;
                        }
                        case 26: {
                            Preferences.debug("ecode = 26 for SHC_DEGREEORDER specifying degree and order of a spherical harmonic basis function\n", 2);
                            this.degreeArray[sphericalHarmonicIndex] = this.getBufferInt(this.bufferByte, currentAddress + 8, endianess);
                            Preferences.debug("Degree = " + this.degreeArray[sphericalHarmonicIndex] + "\n", 2);
                            this.orderArray[sphericalHarmonicIndex] = this.getBufferInt(this.bufferByte, currentAddress + 12, endianess);
                            Preferences.debug("Order = " + this.orderArray[sphericalHarmonicIndex] + "\n", 2);
                            ++sphericalHarmonicIndex;
                            break;
                        }
                        case 30: {
                            Preferences.debug("ecode = 30 for CARET an XML extension\n", 2);
                            this.caretArray[caretIndex] = new String(this.bufferByte, currentAddress + 8, this.esizeArray[ecodeNumber] - 8);
                            Preferences.debug("Caret field = " + this.caretArray[caretIndex] + "\n", 2);
                            ++caretIndex;
                            break;
                        }
                        default: {
                            Preferences.debug("ecode = " + this.ecodeArray[ecodeNumber] + " an unspecified ecode value\n", 2);
                            break;
                        }
                    }
                }
                this.fileInfo.setEsize(this.esizeArray);
                this.fileInfo.setEcode(this.ecodeArray);
                this.fileInfo.setMindIdent(this.mindIdentArray);
                this.fileInfo.setBValue(this.bValueArray);
                this.fileInfo.setAzimuth(this.azimuthArray);
                this.fileInfo.setZenith(this.zenithArray);
                this.fileInfo.setDTComponent(this.dtComponentArray);
                this.fileInfo.setDegree(this.degreeArray);
                this.fileInfo.setOrder(this.orderArray);
                this.fileInfo.setAfniGroup(this.afniGroupArray);
                this.fileInfo.setAsciiText(this.asciiTextArray);
                this.fileInfo.setCaret(this.caretArray);
            }
        }
        if (this.raFile != null) {
            this.raFile.close();
        }
        return true;
    }
    
    public ModelImage readImage(final boolean one, final boolean niftiCompressed) throws IOException, OutOfMemoryError {
        this.fileInfo = new FileInfoNIFTI(this.fileName, this.fileDir, 40);
        boolean flip = false;
        if (niftiCompressed) {
            this.file = new File(this.fileDir + File.separator + this.fileName);
            final String ext = this.fileName.substring(this.fileName.lastIndexOf(".") + 1, this.fileName.length());
            if (ext.equalsIgnoreCase("zip")) {
                try {
                    this.fis = new FileInputStream(this.file);
                }
                catch (FileNotFoundException e3) {
                    MipavUtil.displayError("File not found exception on fis = new FileInputStream(file) for " + this.fileDir + File.separator + this.fileName);
                    return null;
                }
                try {
                    this.zin = new ZipInputStream(new BufferedInputStream(this.fis));
                }
                catch (Exception e4) {
                    MipavUtil.displayError("Exception on ZipInputStream for " + this.fileName);
                    return null;
                }
                if (!this.readHeader(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), true)) {
                    throw new IOException(" NIFTI header file error");
                }
            }
            else if (ext.equalsIgnoreCase("gz")) {
                try {
                    this.fis = new FileInputStream(this.file);
                }
                catch (FileNotFoundException e3) {
                    MipavUtil.displayError("File not found exception on fis = new FileInputStream(file) for " + this.fileDir + File.separator + this.fileName);
                    return null;
                }
                try {
                    this.gzin = new GZIPInputStream(new BufferedInputStream(this.fis));
                }
                catch (IOException e5) {
                    MipavUtil.displayError("IOException on GZIPInputStream for " + this.fileName);
                    return null;
                }
                if (!this.readHeader(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), true)) {
                    throw new IOException(" NIFTI header file error");
                }
            }
            else if (ext.equalsIgnoreCase("bz2")) {
                try {
                    this.fis = new FileInputStream(this.file);
                }
                catch (FileNotFoundException e3) {
                    MipavUtil.displayError("File not found exception on fis = new FileInputStream(file) for " + this.fileDir + File.separator + this.fileName);
                    return null;
                }
                try {
                    this.fis.read();
                }
                catch (IOException e5) {
                    MipavUtil.displayError("IOException on fis.read() trying to read byte B");
                    return null;
                }
                try {
                    this.fis.read();
                }
                catch (IOException e5) {
                    MipavUtil.displayError("IOException on fis.read() trying to read byte Z");
                    return null;
                }
                try {
                    this.bz2in = new CBZip2InputStream((InputStream)new BufferedInputStream(this.fis));
                }
                catch (Exception e4) {
                    MipavUtil.displayError("Exception on CBZip2InputStream for " + this.fileName);
                    return null;
                }
                if (!this.readHeader(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), true)) {
                    throw new IOException(" NIFTI header file error");
                }
            }
        }
        else if (!this.readHeader(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), false)) {
            throw new IOException(" NIFTI header file error");
        }
        int[] extents = null;
        try {
            if (one) {
                extents = new int[this.fileInfo.getExtents().length];
                for (int i = 0; i < extents.length; ++i) {
                    extents[i] = this.fileInfo.getExtents()[i];
                }
                this.image = new ModelImage(this.fileInfo.getDataType(), new int[] { extents[0], extents[1] }, this.fileInfo.getFileName());
            }
            else {
                this.image = new ModelImage(this.fileInfo.getDataType(), this.fileInfo.getExtents(), this.fileInfo.getFileName());
            }
        }
        catch (OutOfMemoryError error) {
            throw error;
        }
        this.axisOrientation = this.fileInfo.getAxisOrientation();
        if (Preferences.is("FlipNIFTIRead") && (this.axisOrientation[1] == 3 || this.axisOrientation[1] == 5)) {
            final int orient = this.axisOrientation[1];
            flip = true;
            if (this.axisOrientation[1] == 3) {
                this.axisOrientation[1] = 4;
            }
            else {
                this.axisOrientation[1] = 6;
            }
            this.fileInfo.setAxisOrientation(this.axisOrientation);
            this.LPSOrigin = this.fileInfo.getOrigin();
            if (orient == 1 || orient == 4 || orient == 5) {
                this.LPSOrigin[1] += (this.fileInfo.getExtents()[1] - 1) * this.fileInfo.getResolutions()[1];
            }
            else {
                this.LPSOrigin[1] -= (this.fileInfo.getExtents()[1] - 1) * this.fileInfo.getResolutions()[1];
            }
            this.fileInfo.setOrigin(this.LPSOrigin);
            this.matrix.set(0, 1, -this.matrix.get(0, 1));
            this.matrix.set(1, 1, -this.matrix.get(1, 1));
            this.matrix.set(2, 1, -this.matrix.get(2, 1));
            this.matrix.set(0, 3, this.LPSOrigin[0]);
            this.matrix.set(1, 3, this.LPSOrigin[1]);
            this.matrix.set(2, 3, this.LPSOrigin[2]);
            if (this.matrix2 != null) {
                this.matrix2.set(0, 1, -this.matrix2.get(0, 1));
                this.matrix2.set(1, 1, -this.matrix2.get(1, 1));
                this.matrix2.set(2, 1, -this.matrix2.get(2, 1));
                this.matrix2.set(0, 3, this.LPSOrigin[0]);
                this.matrix2.set(1, 3, this.LPSOrigin[1]);
                this.matrix2.set(2, 3, this.LPSOrigin[2]);
            }
        }
        extents = this.fileInfo.getExtents();
        if (this.image.getNDims() == 2) {
            this.image.setFileInfo((FileInfoBase)this.fileInfo, 0);
        }
        else if (this.image.getNDims() == 3) {
            for (int i = 0; i < extents[2]; ++i) {
                final FileInfoNIFTI newFileInfo = (FileInfoNIFTI)this.fileInfo.clone();
                newFileInfo.setOrigin(this.fileInfo.getOriginAtSlice(i));
                this.image.setFileInfo((FileInfoBase)newFileInfo, i);
            }
        }
        else if (this.image.getNDims() == 4) {
            for (int i = 0; i < extents[2] * extents[3]; ++i) {
                final FileInfoNIFTI newFileInfo = (FileInfoNIFTI)this.fileInfo.clone();
                newFileInfo.setOrigin(this.fileInfo.getOriginAtSlice(i));
                this.image.setFileInfo((FileInfoBase)newFileInfo, i);
            }
        }
        this.updateorigins(this.image.getFileInfo());
        if (this.image.getNDims() >= 3) {
            this.image.setMatrix(this.matrix);
        }
        else {
            this.image.setMatrix(this.matrixTwoDim);
        }
        if (this.matrix2 != null) {
            this.image.getMatrixHolder().addMatrix(this.matrix2);
        }
        try {
            long offset;
            if (this.oneFile) {
                offset = (long)Math.abs(this.vox_offset);
                if (offset < this.headerSize) {
                    offset = this.headerSize;
                }
            }
            else {
                offset = 0L;
            }
            if (one && this.fileInfo.getExtents().length > 2) {
                offset += this.getOffset(this.fileInfo);
            }
            if (niftiCompressed) {
                byte[] buffer;
                if (this.image.getType() == 9) {
                    buffer = new byte[255];
                }
                else if (this.image.getType() == 10) {
                    buffer = new byte[4];
                }
                else if (this.image.getType() == 11) {
                    buffer = new byte[8];
                }
                else {
                    buffer = new byte[256];
                }
                final String ext2 = this.fileName.substring(this.fileName.lastIndexOf(".") + 1, this.fileName.length());
                if (ext2.equalsIgnoreCase("zip")) {
                    int start = 0;
                    final boolean endianness = this.fileInfo.getEndianess();
                    final int type = this.image.getType();
                    try {
                        while (true) {
                            final int bytesRead = this.zin.read(buffer);
                            if (bytesRead == -1) {
                                break;
                            }
                            if (this.image.getType() == 9) {
                                if (bytesRead != 255) {
                                    buffer = this.getFullBuffer(this.zin, buffer, bytesRead, 255);
                                }
                            }
                            else if (this.image.getType() == 10) {
                                if (bytesRead != 4) {
                                    buffer = this.getFullBuffer(this.zin, buffer, bytesRead, 4);
                                }
                            }
                            else if (this.image.getType() == 11) {
                                if (bytesRead != 8) {
                                    buffer = this.getFullBuffer(this.zin, buffer, bytesRead, 8);
                                }
                            }
                            else if (bytesRead != 256) {
                                buffer = this.getFullBuffer(this.zin, buffer, bytesRead, 256);
                            }
                            if (type == 1 || type == 2 || type == 0) {
                                this.image.importData(start, buffer, false);
                                if (start >= this.image.getDataSize()) {
                                    continue;
                                }
                                start += buffer.length;
                            }
                            else if (type == 3 || type == 4) {
                                final short[] shortBuff = new short[buffer.length / 2];
                                for (int m = 0, k = 0; m < buffer.length; m += 2, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1] };
                                    shortBuff[k] = FileBase.bytesToShort(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, shortBuff, false);
                                }
                                start += shortBuff.length;
                            }
                            else if (type == 5 || type == 14) {
                                final int[] intBuff = new int[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    intBuff[k] = FileBase.bytesToInt(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, intBuff, false);
                                }
                                start += intBuff.length;
                            }
                            else if (type == 7) {
                                final float[] floatBuff = new float[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    floatBuff[k] = FileBase.bytesToFloat(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, floatBuff, false);
                                }
                                start += floatBuff.length;
                            }
                            else if (type == 8) {
                                final double[] doubleBuff = new double[buffer.length / 8];
                                for (int m = 0, k = 0; m < buffer.length; m += 8, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3], buffer[m + 4], buffer[m + 5], buffer[m + 6], buffer[m + 7] };
                                    doubleBuff[k] = FileBase.bytesToDouble(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, doubleBuff, false);
                                }
                                start += doubleBuff.length;
                            }
                            else if (type == 9) {
                                final byte[] buff2 = new byte[buffer.length + buffer.length / 3];
                                if (start < this.image.getDataSize()) {
                                    int counter = 0;
                                    for (int j = 0; j < buffer.length; j += 3) {
                                        buff2[counter] = 1;
                                        buff2[counter + 1] = buffer[j];
                                        buff2[counter + 2] = buffer[j + 1];
                                        buff2[counter + 3] = buffer[j + 2];
                                        counter += 4;
                                    }
                                    this.image.importData(start, buff2, false);
                                }
                                start += buff2.length;
                            }
                            else if (type == 10) {
                                final short[] shortBuff2 = { 1, 0, 0 };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 2, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1] };
                                        shortBuff2[k] = FileBase.bytesToShort(endianness, 0, b);
                                    }
                                    this.image.importData(start, shortBuff2, false);
                                }
                                start += shortBuff2.length;
                            }
                            else {
                                if (type != 11) {
                                    continue;
                                }
                                final float[] floatBuff2 = { 1.0f, 0.0f, 0.0f };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 4, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                        floatBuff2[k] = FileBase.bytesToFloat(endianness, 0, b);
                                    }
                                    this.image.importData(start, floatBuff2, false);
                                }
                                start += floatBuff2.length;
                            }
                        }
                    }
                    catch (IOException e) {
                        e.printStackTrace();
                        MipavUtil.displayError("IOException on gzin.read(buffer) for " + this.fileName);
                        return null;
                    }
                }
                else if (ext2.equalsIgnoreCase("gz")) {
                    int start = 0;
                    final boolean endianness = this.fileInfo.getEndianess();
                    final int type = this.image.getType();
                    try {
                        while (true) {
                            final int bytesRead = this.gzin.read(buffer);
                            if (bytesRead == -1) {
                                break;
                            }
                            if (this.image.getType() == 9) {
                                if (bytesRead != 255) {
                                    buffer = this.getFullBuffer(this.gzin, buffer, bytesRead, 255);
                                }
                            }
                            else if (this.image.getType() == 10) {
                                if (bytesRead != 4) {
                                    buffer = this.getFullBuffer(this.gzin, buffer, bytesRead, 4);
                                }
                            }
                            else if (this.image.getType() == 11) {
                                if (bytesRead != 8) {
                                    buffer = this.getFullBuffer(this.gzin, buffer, bytesRead, 8);
                                }
                            }
                            else if (bytesRead != 256) {
                                buffer = this.getFullBuffer(this.gzin, buffer, bytesRead, 256);
                            }
                            if (type == 1 || type == 2 || type == 0) {
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, buffer, false);
                                }
                                start += buffer.length;
                            }
                            else if (type == 3 || type == 4) {
                                final short[] shortBuff = new short[buffer.length / 2];
                                for (int m = 0, k = 0; m < buffer.length; m += 2, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1] };
                                    shortBuff[k] = FileBase.bytesToShort(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, shortBuff, false);
                                }
                                start += shortBuff.length;
                            }
                            else if (type == 5 || type == 14) {
                                final int[] intBuff = new int[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    intBuff[k] = FileBase.bytesToInt(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, intBuff, false);
                                }
                                start += intBuff.length;
                            }
                            else if (type == 7) {
                                final float[] floatBuff = new float[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    floatBuff[k] = FileBase.bytesToFloat(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, floatBuff, false);
                                }
                                start += floatBuff.length;
                            }
                            else if (type == 8) {
                                final double[] doubleBuff = new double[buffer.length / 8];
                                for (int m = 0, k = 0; m < buffer.length; m += 8, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3], buffer[m + 4], buffer[m + 5], buffer[m + 6], buffer[m + 7] };
                                    doubleBuff[k] = FileBase.bytesToDouble(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, doubleBuff, false);
                                }
                                start += doubleBuff.length;
                            }
                            else if (type == 9) {
                                final byte[] buff2 = new byte[buffer.length + buffer.length / 3];
                                if (start < this.image.getDataSize()) {
                                    int counter = 0;
                                    for (int j = 0; j < buffer.length; j += 3) {
                                        buff2[counter] = 1;
                                        buff2[counter + 1] = buffer[j];
                                        buff2[counter + 2] = buffer[j + 1];
                                        buff2[counter + 3] = buffer[j + 2];
                                        counter += 4;
                                    }
                                    this.image.importData(start, buff2, false);
                                }
                                start += buff2.length;
                            }
                            else if (type == 10) {
                                final short[] shortBuff2 = { 1, 0, 0 };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 2, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1] };
                                        shortBuff2[k] = FileBase.bytesToShort(endianness, 0, b);
                                    }
                                    this.image.importData(start, shortBuff2, false);
                                }
                                start += shortBuff2.length;
                            }
                            else {
                                if (type != 11) {
                                    continue;
                                }
                                final float[] floatBuff2 = { 1.0f, 0.0f, 0.0f };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 4, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                        floatBuff2[k] = FileBase.bytesToFloat(endianness, 0, b);
                                    }
                                    this.image.importData(start, floatBuff2, false);
                                }
                                start += floatBuff2.length;
                            }
                        }
                    }
                    catch (IOException e) {
                        Preferences.debug("IOException on gzin.read(buffer) for " + this.fileName, 2);
                        throw e;
                    }
                }
                else if (ext2.equalsIgnoreCase("bz2")) {
                    int start = 0;
                    final boolean endianness = this.fileInfo.getEndianess();
                    final int type = this.image.getType();
                    try {
                        while (true) {
                            final int bytesRead = this.bz2in.read(buffer);
                            if (bytesRead == -1) {
                                break;
                            }
                            if (this.image.getType() == 9) {
                                if (bytesRead != 255) {
                                    buffer = this.getFullBuffer((InputStream)this.bz2in, buffer, bytesRead, 255);
                                }
                            }
                            else if (this.image.getType() == 10) {
                                if (bytesRead != 4) {
                                    buffer = this.getFullBuffer((InputStream)this.bz2in, buffer, bytesRead, 4);
                                }
                            }
                            else if (this.image.getType() == 11) {
                                if (bytesRead != 8) {
                                    buffer = this.getFullBuffer((InputStream)this.bz2in, buffer, bytesRead, 8);
                                }
                            }
                            else if (bytesRead != 256) {
                                buffer = this.getFullBuffer((InputStream)this.bz2in, buffer, bytesRead, 256);
                            }
                            if (type == 1 || type == 2 || type == 0) {
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, buffer, false);
                                }
                                start += buffer.length;
                            }
                            else if (type == 3 || type == 4) {
                                final short[] shortBuff = new short[buffer.length / 2];
                                for (int m = 0, k = 0; m < buffer.length; m += 2, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1] };
                                    shortBuff[k] = FileBase.bytesToShort(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, shortBuff, false);
                                }
                                start += shortBuff.length;
                            }
                            else if (type == 5 || type == 14) {
                                final int[] intBuff = new int[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    intBuff[k] = FileBase.bytesToInt(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, intBuff, false);
                                }
                                start += intBuff.length;
                            }
                            else if (type == 7) {
                                final float[] floatBuff = new float[buffer.length / 4];
                                for (int m = 0, k = 0; m < buffer.length; m += 4, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                    floatBuff[k] = FileBase.bytesToFloat(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, floatBuff, false);
                                }
                                start += floatBuff.length;
                            }
                            else if (type == 8) {
                                final double[] doubleBuff = new double[buffer.length / 8];
                                for (int m = 0, k = 0; m < buffer.length; m += 8, ++k) {
                                    final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3], buffer[m + 4], buffer[m + 5], buffer[m + 6], buffer[m + 7] };
                                    doubleBuff[k] = FileBase.bytesToDouble(endianness, 0, b);
                                }
                                if (start < this.image.getDataSize()) {
                                    this.image.importData(start, doubleBuff, false);
                                }
                                start += doubleBuff.length;
                            }
                            else if (type == 9) {
                                final byte[] buff2 = new byte[buffer.length + buffer.length / 3];
                                if (start < this.image.getDataSize()) {
                                    int counter = 0;
                                    for (int j = 0; j < buffer.length; j += 3) {
                                        buff2[counter] = 1;
                                        buff2[counter + 1] = buffer[j];
                                        buff2[counter + 2] = buffer[j + 1];
                                        buff2[counter + 3] = buffer[j + 2];
                                        counter += 4;
                                    }
                                    this.image.importData(start, buff2, false);
                                }
                                start += buff2.length;
                            }
                            else if (type == 10) {
                                final short[] shortBuff2 = { 1, 0, 0 };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 2, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1] };
                                        shortBuff2[k] = FileBase.bytesToShort(endianness, 0, b);
                                    }
                                    this.image.importData(start, shortBuff2, false);
                                }
                                start += shortBuff2.length;
                            }
                            else {
                                if (type != 11) {
                                    continue;
                                }
                                final float[] floatBuff2 = { 1.0f, 0.0f, 0.0f };
                                if (start < this.image.getDataSize()) {
                                    for (int m = 0, k = 1; m < buffer.length; m += 4, ++k) {
                                        final byte[] b = { buffer[m], buffer[m + 1], buffer[m + 2], buffer[m + 3] };
                                        floatBuff2[k] = FileBase.bytesToFloat(endianness, 0, b);
                                    }
                                    this.image.importData(start, floatBuff2, false);
                                }
                                start += floatBuff2.length;
                            }
                        }
                    }
                    catch (IOException e) {
                        Preferences.debug("IOException on gzin.read(buffer) for " + this.fileName, 2);
                        throw e;
                    }
                }
            }
            else {
                final FileRaw rawFile = new FileRaw(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), (FileInfoBase)this.fileInfo, 0);
                this.linkProgress((FileBase)rawFile);
                rawFile.readImage(this.image, offset);
            }
            if (this.zin != null) {
                this.zin.close();
            }
            if (this.gzin != null) {
                this.gzin.close();
            }
            if (this.bz2in != null) {
                this.bz2in.close();
            }
            if (this.vox_offset < 0.0f) {
                this.absoluteValue(this.image);
            }
            if (this.scl_slope != 0.0 && (this.scl_slope != 1.0f || this.scl_inter != 0.0f) && this.sourceType != 32 && this.sourceType != 64 && this.sourceType != 128) {
                if (this.image.getType() == 12 || this.image.getType() == 13) {
                    this.image.calcMinMaxMag(Preferences.is("LogMagDisplay"));
                }
                else {
                    this.image.calcMinMax();
                }
                this.imageMin = this.image.getMin();
                this.imageMax = this.image.getMax();
                int newType = this.image.getType();
                final double m2 = this.scl_slope * this.imageMin + this.scl_inter;
                final double m3 = this.scl_slope * this.imageMax + this.scl_inter;
                this.newMin = Math.min(m2, m3);
                this.newMax = Math.max(m2, m3);
                boolean needFloat = false;
                boolean doChangeType = false;
                if (this.scl_slope != Math.round(this.scl_slope) || this.scl_inter != Math.round(this.scl_inter)) {
                    needFloat = true;
                }
                if (needFloat && this.newMax <= 3.4028234663852886E38 && this.newMin >= -3.4028234663852886E38 && (this.sourceType == 1 || this.sourceType == 256 || this.sourceType == 4 || this.sourceType == 8 || this.sourceType == 2 || this.sourceType == 512 || this.sourceType == 768)) {
                    newType = 7;
                    doChangeType = true;
                    this.fileInfo.setBitPix((short)32);
                }
                else if (!needFloat || this.newMax > 3.4028234663852886E38 || this.newMin < -3.4028234663852886E38 || this.sourceType != 16) {
                    if (needFloat) {
                        newType = 8;
                        doChangeType = true;
                        this.fileInfo.setBitPix((short)64);
                    }
                    else if (this.newMax > this.imageMax || this.newMin < this.imageMin) {
                        if (this.newMin >= -128.0 && this.newMax <= 127.0) {
                            newType = 1;
                            doChangeType = true;
                            this.fileInfo.setBitPix((short)8);
                        }
                        else if (this.newMin >= -32768.0 && this.newMax <= 32767.0) {
                            newType = 3;
                            doChangeType = true;
                            this.fileInfo.setBitPix((short)16);
                        }
                        else if (this.newMin >= -2.147483648E9 && this.newMax <= 2.147483647E9) {
                            newType = 5;
                            doChangeType = true;
                            this.fileInfo.setBitPix((short)32);
                        }
                        else if (this.newMin >= -9.223372036854776E18 && this.newMax <= 9.223372036854776E18) {
                            newType = 6;
                            doChangeType = true;
                            this.fileInfo.setBitPix((short)64);
                        }
                        else {
                            newType = 8;
                            doChangeType = true;
                            this.fileInfo.setBitPix((short)64);
                        }
                    }
                }
                if (doChangeType) {
                    AlgorithmChangeType changeTypeAlgo = new AlgorithmChangeType(this.image, newType, this.imageMin, this.imageMax, this.imageMin, this.imageMax, false);
                    changeTypeAlgo.run();
                    changeTypeAlgo.finalize();
                    changeTypeAlgo = null;
                }
                this.scale(this.image);
            }
            if (flip) {
                this.flipTopBottom(this.image);
            }
            if (one) {
                this.fileInfo.setExtents(extents);
            }
            this.fireProgressStateChanged(100);
        }
        catch (IOException error2) {
            throw error2;
        }
        catch (OutOfMemoryError e2) {
            throw e2;
        }
        return this.image;
    }
    
    private byte[] getFullBuffer(final InputStream in, final byte[] buff, final int off, final int fullBufferSize) throws IOException {
        int bytesRead = 0;
        int offset = off;
        while (true) {
            if (offset == fullBufferSize) {
                if (bytesRead == -1) {
                    break;
                }
            }
            try {
                bytesRead = in.read(buff, offset, fullBufferSize - offset);
                if (bytesRead != -1) {
                    if (bytesRead != 0) {
                        offset += bytesRead;
                        if (offset != fullBufferSize) {
                            continue;
                        }
                    }
                }
            }
            catch (IOException e) {
                throw e;
            }
            break;
        }
        if (offset != fullBufferSize) {
            final byte[] buffShortened = new byte[offset];
            for (int i = 0; i < offset; ++i) {
                buffShortened[i] = buff[i];
            }
            return buffShortened;
        }
        return buff;
    }
    
    public void readImage(final float[] buffer) throws IOException, OutOfMemoryError {
        this.readImage(buffer, 0L);
    }
    
    public void readImage(final float[] buffer, final long userOffset) throws IOException, OutOfMemoryError {
        if (this.fileInfo == null) {
            this.fileInfo = new FileInfoNIFTI(this.fileName, this.fileDir, 40);
            if (!this.readHeader(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), false)) {
                throw new IOException("Cannot read image because of NIFTI header file error");
            }
        }
        try {
            final FileRaw rawFile = new FileRaw(this.fileInfo.getFileName(), this.fileInfo.getFileDirectory(), (FileInfoBase)this.fileInfo, 0);
            this.linkProgress((FileBase)rawFile);
            long offset;
            if (this.oneFile) {
                offset = userOffset + (long)Math.abs(this.vox_offset);
                if (offset < this.headerSize) {
                    offset = userOffset + this.headerSize;
                }
            }
            else {
                offset = userOffset;
            }
            rawFile.readImage(buffer, offset, this.fileInfo.getDataType());
            if (this.vox_offset < 0.0f) {
                for (int i = 0; i < buffer.length; ++i) {
                    buffer[i] = Math.abs(buffer[i]);
                }
            }
            if (this.scl_slope != 0.0 && (this.scl_slope != 1.0f || this.scl_inter != 0.0f)) {
                for (int i = 0; i < buffer.length; ++i) {
                    buffer[i] = buffer[i] * this.scl_slope + this.scl_inter;
                }
            }
            this.axisOrientation = this.fileInfo.getAxisOrientation();
            if (Preferences.is("FlipNIFTIRead") && (this.axisOrientation[1] == 3 || this.axisOrientation[1] == 5)) {
                final int orient = this.axisOrientation[1];
                if (this.axisOrientation[1] == 3) {
                    this.axisOrientation[1] = 4;
                }
                else {
                    this.axisOrientation[1] = 6;
                }
                this.fileInfo.setAxisOrientation(this.axisOrientation);
                this.LPSOrigin = this.fileInfo.getOrigin();
                if (orient == 1 || orient == 4 || orient == 5) {
                    this.LPSOrigin[1] += (this.fileInfo.getExtents()[1] - 1) * this.fileInfo.getResolutions()[1];
                }
                else {
                    this.LPSOrigin[1] -= (this.fileInfo.getExtents()[1] - 1) * this.fileInfo.getResolutions()[1];
                }
                this.fileInfo.setOrigin(this.LPSOrigin);
                this.matrix.set(0, 1, -this.matrix.get(0, 1));
                this.matrix.set(1, 1, -this.matrix.get(1, 1));
                this.matrix.set(2, 1, -this.matrix.get(2, 1));
                this.matrix.set(0, 3, this.LPSOrigin[0]);
                this.matrix.set(1, 3, this.LPSOrigin[1]);
                this.matrix.set(2, 3, this.LPSOrigin[2]);
                if (this.matrix2 != null) {
                    this.matrix2.set(0, 1, -this.matrix2.get(0, 1));
                    this.matrix2.set(1, 1, -this.matrix2.get(1, 1));
                    this.matrix2.set(2, 1, -this.matrix2.get(2, 1));
                    this.matrix2.set(0, 3, this.LPSOrigin[0]);
                    this.matrix2.set(1, 3, this.LPSOrigin[1]);
                    this.matrix2.set(2, 3, this.LPSOrigin[2]);
                }
                this.flipTopBottom(buffer, this.fileInfo);
            }
            rawFile.raFile.close();
            this.fireProgressStateChanged(100);
        }
        catch (IOException error) {
            error.printStackTrace();
            throw new IOException("FileNIFTI: " + error);
        }
        catch (OutOfMemoryError e) {
            throw e;
        }
    }
    
    public void scale(final ModelImage image) throws IOException {
        try {
            float[] buffer = null;
            int bufferSize;
            if (image.getNDims() > 1) {
                bufferSize = image.getSliceSize();
            }
            else {
                bufferSize = image.getExtents()[0];
            }
            int nBuffers;
            if (image.getNDims() == 5) {
                nBuffers = image.getExtents()[4] * image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 4) {
                nBuffers = image.getExtents()[3] * image.getExtents()[2];
            }
            else if (image.getNDims() == 3) {
                nBuffers = image.getExtents()[2];
            }
            else {
                nBuffers = 1;
            }
            if (image.isColorImage()) {
                buffer = new float[bufferSize * 4];
                bufferSize *= 4;
                final int xDim = image.getExtents()[0] * 4;
                final int yDim = image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    image.exportData(k * bufferSize, bufferSize, buffer);
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; i += 4) {
                            buffer[j * xDim + i] = 255.0f;
                            buffer[j * xDim + i + 1] = this.scl_slope * buffer[j * xDim + i + 1] + this.scl_inter;
                            buffer[j * xDim + i + 2] = this.scl_slope * buffer[j * xDim + i + 2] + this.scl_inter;
                            buffer[j * xDim + i + 3] = this.scl_slope * buffer[j * xDim + i + 3] + this.scl_inter;
                        }
                    }
                    image.importData(k * bufferSize, buffer, false);
                }
            }
            else {
                buffer = new float[bufferSize];
                final int xDim = image.getExtents()[0];
                final int yDim = image.getExtents()[1];
                for (int k = 0; k < nBuffers; ++k) {
                    image.exportData(k * bufferSize, bufferSize, buffer);
                    for (int j = 0; j < yDim; ++j) {
                        for (int i = 0; i < xDim; ++i) {
                            buffer[j * xDim + i] = this.scl_slope * buffer[j * xDim + i] + this.scl_inter;
                        }
                    }
                    image.importData(k * bufferSize, buffer, false);
                }
            }
        }
        catch (IOException error) {
            throw new IOException("FileNIFTI.scale: " + error);
        }
        catch (OutOfMemoryError error2) {
            throw error2;
        }
    }
    
    public void writeImage(final ModelImage image, final FileWriteOptions options) throws IOException {
        final String suffix = FileUtility.getExtension(this.fileName);
        if (suffix.equalsIgnoreCase(".nii")) {
            this.oneFile = true;
        }
        else if (suffix.equalsIgnoreCase(".img")) {
            this.oneFile = false;
        }
        else if (suffix.equalsIgnoreCase(".hdr")) {
            this.oneFile = false;
        }
        else {
            final JDialogNIFTIChoice choice = new JDialogNIFTIChoice((Frame)ViewUserInterface.getReference().getMainFrame());
            if (!choice.okayPressed()) {
                throw new IOException("FileNIFTIWrite dialog error");
            }
            this.oneFile = choice.getOneFile();
        }
        final int index = this.fileName.lastIndexOf(".");
        String fhName;
        if (index != -1) {
            fhName = this.fileName.substring(0, index);
            if (suffix.equalsIgnoreCase(".hdr")) {
                this.fileName = fhName + ".img";
            }
        }
        else {
            fhName = this.fileName.substring(0);
        }
        if (options.isMultiFile()) {
            final FileRaw rawFile = new FileRaw(image.getFileInfo(0));
            rawFile.setZeroLengthFlag(true);
            this.linkProgress((FileBase)rawFile);
            if (this.oneFile) {
                rawFile.setStartPosition(352L);
                if (image.getNDims() == 3) {
                    rawFile.writeImage3DTo2D(image, options, ".nii");
                    this.writeHeader3DTo2D(image, fhName, this.fileDir, options, this.oneFile);
                }
                else if (image.getNDims() == 4) {
                    rawFile.writeImage4DTo3D(image, options, ".nii");
                    this.writeHeader4DTo3D(image, fhName, this.fileDir, options, this.oneFile);
                }
            }
            else {
                rawFile.setStartPosition(0L);
                if (image.getNDims() == 3) {
                    rawFile.writeImage3DTo2D(image, options, ".img");
                    this.writeHeader3DTo2D(image, fhName, this.fileDir, options, this.oneFile);
                }
                else if (image.getNDims() == 4) {
                    rawFile.writeImage4DTo3D(image, options, ".img");
                    this.writeHeader4DTo3D(image, fhName, this.fileDir, options, this.oneFile);
                }
            }
        }
        else {
            try {
                final FileRaw rawFile = new FileRaw(this.fileName, this.fileDir, image.getFileInfo(0), 1);
                rawFile.setZeroLengthFlag(true);
                this.linkProgress((FileBase)rawFile);
                if (this.oneFile) {
                    rawFile.setStartPosition(352L);
                }
                else {
                    rawFile.setStartPosition(0L);
                }
                rawFile.writeImage(image, options);
                final int nImagesSaved = rawFile.getNImages();
                final int nTimePeriodsSaved = rawFile.getNTimePeriods();
                if (nImagesSaved != 0) {
                    this.writeHeader(image, nImagesSaved, nTimePeriodsSaved, fhName, this.fileDir, false, this.oneFile);
                }
            }
            catch (IOException error) {
                throw new IOException("FileNIFTIWrite: " + error);
            }
            catch (OutOfMemoryError error2) {
                throw error2;
            }
        }
        this.fireProgressStateChanged(100);
    }
    
    private int[] getAxisOrientation(final TransMatrix mat) {
        final int[] axisOrientation = new int[3];
        int k = 0;
        double xi = mat.get(0, 0);
        double xj = mat.get(0, 1);
        double xk = mat.get(0, 2);
        double yi = mat.get(1, 0);
        double yj = mat.get(1, 1);
        double yk = mat.get(1, 2);
        double zi = mat.get(2, 0);
        double zj = mat.get(2, 1);
        double zk = mat.get(2, 2);
        axisOrientation[0] = 0;
        double val = Math.sqrt(xi * xi + yi * yi + zi * zi);
        if (val == 0.0) {
            return axisOrientation;
        }
        xi /= val;
        yi /= val;
        zi /= val;
        val = Math.sqrt(xj * xj + yj * yj + zj * zj);
        if (val == 0.0) {
            return axisOrientation;
        }
        xj /= val;
        yj /= val;
        zj /= val;
        val = xi * xj + yi * yj + zi * zj;
        if (Math.abs(val) > 1.0E-4) {
            xj -= val * xi;
            yj -= val * yi;
            zj -= val * zi;
            val = Math.sqrt(xj * xj + yj * yj + zj * zj);
            if (val == 0.0) {
                return axisOrientation;
            }
            xj /= val;
            yj /= val;
            zj /= val;
        }
        val = Math.sqrt(xk * xk + yk * yk + zk * zk);
        if (val == 0.0) {
            xk = yi * zj - zi * yj;
            yk = zi * xj - zj * xi;
            zk = xi * yj - yi * xj;
        }
        else {
            xk /= val;
            yk /= val;
            zk /= val;
        }
        val = xi * xk + yi * yk + zi * zk;
        if (Math.abs(val) > 1.0E-4) {
            xk -= val * xi;
            yk -= val * yi;
            zk -= val * zi;
            val = Math.sqrt(xk * xk + yk * yk + zk * zk);
            if (val == 0.0) {
                return axisOrientation;
            }
            xk /= val;
            yk /= val;
            zk /= val;
        }
        val = xj * xk + yj * yk + zj * zk;
        if (Math.abs(val) > 1.0E-4) {
            xk -= val * xj;
            yk -= val * yj;
            zk -= val * zj;
            val = Math.sqrt(xk * xk + yk * yk + zk * zk);
            if (val == 0.0) {
                return axisOrientation;
            }
            xk /= val;
            yk /= val;
            zk /= val;
        }
        final Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, xi);
        Q.set(0, 1, xj);
        Q.set(0, 2, xk);
        Q.set(1, 0, yi);
        Q.set(1, 1, yj);
        Q.set(1, 2, yk);
        Q.set(2, 0, zi);
        Q.set(2, 1, zj);
        Q.set(2, 2, zk);
        final double detQ = Q.det();
        if (detQ == 0.0) {
            MipavUtil.displayError("detQ == 0.0 in getAxisOrientation");
            return axisOrientation;
        }
        final Matrix P = new Matrix(3, 3);
        double vbest = -666.0;
        int rbest;
        int qbest;
        int ibest;
        int pbest = ibest = (qbest = (rbest = 1));
        int jbest = 2;
        int kbest = 3;
        int i;
        int j = 0;
        for (i = 1; i <= 3; ++i) {
            for (j = 1; j <= 3; ++j) {
                if (i != j) {
                    for (k = 1; k <= 3; ++k) {
                        if (i != k) {
                            if (j != k) {
                                P.set(0, 0, 0.0);
                                P.set(0, 1, 0.0);
                                P.set(0, 2, 0.0);
                                P.set(1, 0, 0.0);
                                P.set(1, 1, 0.0);
                                P.set(1, 2, 0.0);
                                P.set(2, 0, 0.0);
                                P.set(2, 1, 0.0);
                                P.set(2, 2, 0.0);
                                for (int p = -1; p <= 1; p += 2) {
                                    for (int q = -1; q <= 1; q += 2) {
                                        for (int r = -1; r <= 1; r += 2) {
                                            P.set(0, i - 1, (double)p);
                                            P.set(1, j - 1, (double)q);
                                            P.set(2, k - 1, (double)r);
                                            final double detP = P.det();
                                            if (detP * detQ > 0.0) {
                                                final Matrix M = P.times(Q);
                                                val = M.get(0, 0) + M.get(1, 1) + M.get(2, 2);
                                                if (val > vbest) {
                                                    vbest = val;
                                                    ibest = i;
                                                    jbest = j;
                                                    kbest = k;
                                                    pbest = p;
                                                    qbest = q;
                                                    rbest = r;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        switch (ibest * pbest) {
            case 1: {
                i = 1;
                break;
            }
            case -1: {
                i = 2;
                break;
            }
            case 2: {
                i = 4;
                break;
            }
            case -2: {
                i = 3;
                break;
            }
            case 3: {
                i = 5;
                break;
            }
            case -3: {
                i = 6;
                break;
            }
        }
        switch (jbest * qbest) {
            case 1: {
                j = 1;
                break;
            }
            case -1: {
                j = 2;
                break;
            }
            case 2: {
                j = 4;
                break;
            }
            case -2: {
                j = 3;
                break;
            }
            case 3: {
                j = 5;
                break;
            }
            case -3: {
                j = 6;
                break;
            }
            default: {
                j = 1;
                break;
            }
        }
        switch (kbest * rbest) {
            case 1: {
                k = 1;
                break;
            }
            case -1: {
                k = 2;
                break;
            }
            case 2: {
                k = 4;
                break;
            }
            case -2: {
                k = 3;
                break;
            }
            case 3: {
                k = 5;
                break;
            }
            case -3: {
                k = 6;
                break;
            }
        }
        axisOrientation[0] = i;
        axisOrientation[1] = j;
        axisOrientation[2] = k;
        return axisOrientation;
    }
    
    private int getOffset(final FileInfoNIFTI fileInfo) {
        int offset = fileInfo.getExtents()[0] * fileInfo.getExtents()[1] * (fileInfo.getExtents()[2] / 2);
        switch (fileInfo.getSourceType()) {
            case 4:
            case 512: {
                offset *= 2;
                break;
            }
            case 8:
            case 16: {
                offset *= 4;
                break;
            }
            case 32:
            case 64:
            case 1024: {
                offset *= 8;
                break;
            }
            case 1792: {
                offset *= 16;
                break;
            }
            case 128: {
                offset *= 3;
                break;
            }
        }
        return offset;
    }
    
    private double mat33_colnorm(final Matrix A) {
        double r1 = Math.abs(A.get(0, 0)) + Math.abs(A.get(1, 0)) + Math.abs(A.get(2, 0));
        final double r2 = Math.abs(A.get(0, 1)) + Math.abs(A.get(1, 1)) + Math.abs(A.get(2, 1));
        final double r3 = Math.abs(A.get(0, 2)) + Math.abs(A.get(1, 2)) + Math.abs(A.get(2, 2));
        if (r1 < r2) {
            r1 = r2;
        }
        if (r1 < r3) {
            r1 = r3;
        }
        return r1;
    }
    
    private Matrix mat33_polar(final Matrix A) {
        double dif = 1.0;
        int k = 0;
        Matrix X = A.copy();
        final Matrix Z = new Matrix(3, 3);
        for (double gam = X.det(); gam == 0.0; gam = X.det()) {
            gam = 1.0E-5 * (0.001 + this.mat33_rownorm(X));
            double val = X.get(0, 0);
            X.set(0, 0, val + gam);
            val = X.get(1, 1);
            X.set(1, 1, val + gam);
            val = X.get(2, 2);
            X.set(2, 2, val + gam);
        }
        while (true) {
            final Matrix Y = X.inverse();
            double gam;
            double gmi;
            if (dif > 0.3) {
                final double alp = Math.sqrt(this.mat33_rownorm(X) * this.mat33_colnorm(X));
                final double bet = Math.sqrt(this.mat33_rownorm(Y) * this.mat33_colnorm(Y));
                gam = Math.sqrt(bet / alp);
                gmi = 1.0 / gam;
            }
            else {
                gmi = (gam = 1.0);
            }
            Z.set(0, 0, 0.5 * (gam * X.get(0, 0) + gmi * Y.get(0, 0)));
            Z.set(0, 1, 0.5 * (gam * X.get(0, 1) + gmi * Y.get(1, 0)));
            Z.set(0, 2, 0.5 * (gam * X.get(0, 2) + gmi * Y.get(2, 0)));
            Z.set(1, 0, 0.5 * (gam * X.get(1, 0) + gmi * Y.get(0, 1)));
            Z.set(1, 1, 0.5 * (gam * X.get(1, 1) + gmi * Y.get(1, 1)));
            Z.set(1, 2, 0.5 * (gam * X.get(1, 2) + gmi * Y.get(2, 1)));
            Z.set(2, 0, 0.5 * (gam * X.get(2, 0) + gmi * Y.get(0, 2)));
            Z.set(2, 1, 0.5 * (gam * X.get(2, 1) + gmi * Y.get(1, 2)));
            Z.set(2, 2, 0.5 * (gam * X.get(2, 2) + gmi * Y.get(2, 2)));
            dif = Math.abs(Z.get(0, 0) - X.get(0, 0)) + Math.abs(Z.get(0, 1) - X.get(0, 1)) + Math.abs(Z.get(0, 2) - X.get(0, 2)) + Math.abs(Z.get(1, 0) - X.get(1, 0)) + Math.abs(Z.get(1, 1) - X.get(1, 1)) + Math.abs(Z.get(1, 2) - X.get(1, 2)) + Math.abs(Z.get(2, 0) - X.get(2, 0)) + Math.abs(Z.get(2, 1) - X.get(2, 1)) + Math.abs(Z.get(2, 2) - X.get(2, 2));
            ++k;
            if (k > 100 || dif < 3.0E-6) {
                break;
            }
            X = Z.copy();
        }
        return Z;
    }
    
    private double mat33_rownorm(final Matrix A) {
        double r1 = Math.abs(A.get(0, 0)) + Math.abs(A.get(0, 1)) + Math.abs(A.get(0, 2));
        final double r2 = Math.abs(A.get(1, 0)) + Math.abs(A.get(1, 1)) + Math.abs(A.get(1, 2));
        final double r3 = Math.abs(A.get(2, 0)) + Math.abs(A.get(2, 1)) + Math.abs(A.get(2, 2));
        if (r1 < r2) {
            r1 = r2;
        }
        if (r1 < r3) {
            r1 = r3;
        }
        return r1;
    }
    
    private void updateorigins(final FileInfoBase[] fileInfo) {
        final float[] origin = fileInfo[0].getOrigin().clone();
        final float[] resolutions = fileInfo[0].getResolutions();
        if (this.image.getNDims() == 3) {
            for (int i = 0; i < this.image.getExtents()[2]; ++i) {
                fileInfo[i].setOrigin(origin[0] + this.matrix.get(0, 2) * i, 0);
                fileInfo[i].setOrigin(origin[1] + this.matrix.get(1, 2) * i, 1);
                fileInfo[i].setOrigin(origin[2] + this.matrix.get(2, 2) * i, 2);
            }
        }
        else if (this.image.getNDims() == 4) {
            for (int i = 0; i < this.image.getExtents()[3]; ++i) {
                for (int j = 0; j < this.image.getExtents()[2]; ++j) {
                    fileInfo[i * this.image.getExtents()[2] + j].setOrigin(origin[0] + this.matrix.get(0, 2) * j, 0);
                    fileInfo[i * this.image.getExtents()[2] + j].setOrigin(origin[1] + this.matrix.get(1, 2) * j, 1);
                    fileInfo[i * this.image.getExtents()[2] + j].setOrigin(origin[2] + this.matrix.get(2, 2) * j, 2);
                    fileInfo[i * this.image.getExtents()[2] + j].setOrigin(origin[3] + i * resolutions[3], 3);
                }
            }
        }
    }
    
    public boolean writeHeader(final ModelImage image, final int nImagesSaved, final int nTimeSaved, final String fileName, final String fileDir, final boolean doGzip, final boolean oneFile) throws IOException {
        boolean isNIFTI = true;
        int firstSpatialUnits = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
        int firstTimeUnits = FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum();
        int niftiSpatialUnits = 0;
        int niftiTimeUnits = 0;
        final float[] niftiOrigin = new float[3];
        float[] niftiOriginS = null;
        MatrixHolder matHolder = null;
        TransMatrix[] matrixArray = null;
        TransMatrix matrixQ = null;
        TransMatrix matrixS = null;
        int transformIDQ = 0;
        int transformIDS = 0;
        int qform_code = 0;
        int sform_code = 0;
        double xDel = 0.0;
        double yDel = 0.0;
        double zDel = 0.0;
        TransMatrix dicomMatrix = null;
        final FileInfoBase myFileInfo = image.getFileInfo(0);
        final boolean endianess = myFileInfo.getEndianess();
        try {
            this.fileInfo = (FileInfoNIFTI)image.getFileInfo(0);
        }
        catch (ClassCastException e) {
            this.fileInfo = new FileInfoNIFTI(fileName, fileDir, 40);
            isNIFTI = false;
            try {
                final FileInfoDicom fileInfoDicom = (FileInfoDicom)image.getFileInfo(0);
                dicomMatrix = fileInfoDicom.getPatientOrientation();
                final String orientation = (String)fileInfoDicom.getTagTable().getValue("0020,0032");
                if (orientation != null) {
                    int index1 = -1;
                    int index2 = -1;
                    for (int i = 0; i < orientation.length(); ++i) {
                        if (orientation.charAt(i) == '\\') {
                            if (index1 == -1) {
                                index1 = i;
                            }
                            else {
                                index2 = i;
                            }
                        }
                    }
                    final double[] coord = { Double.valueOf(orientation.substring(0, index1)), Double.valueOf(orientation.substring(index1 + 1, index2)), Double.valueOf(orientation.substring(index2 + 1)) };
                    final FileInfoDicom fileInfoDicom2 = (FileInfoDicom)image.getFileInfo(1);
                    final String orientation2 = (String)fileInfoDicom2.getTagTable().getValue("0020,0032");
                    if (orientation2 != null) {
                        index1 = -1;
                        index2 = -1;
                        for (int i = 0; i < orientation2.length(); ++i) {
                            if (orientation2.charAt(i) == '\\') {
                                if (index1 == -1) {
                                    index1 = i;
                                }
                                else {
                                    index2 = i;
                                }
                            }
                        }
                        final double[] coord2 = { Double.valueOf(orientation2.substring(0, index1)), Double.valueOf(orientation2.substring(index1 + 1, index2)), Double.valueOf(orientation2.substring(index2 + 1)) };
                        xDel = coord2[0] - coord[0];
                        yDel = coord2[1] - coord[1];
                        zDel = coord2[2] - coord[2];
                    }
                }
            }
            catch (ClassCastException ex) {}
        }
        String fileHeaderName;
        if (oneFile) {
            fileHeaderName = fileName + ".nii";
        }
        else {
            fileHeaderName = fileName + ".hdr";
        }
        this.fileHeader = new File(fileDir + fileHeaderName);
        if (!doGzip) {
            this.raFile = new RandomAccessFile(this.fileHeader, "rw");
        }
        if (!oneFile && !doGzip) {
            this.raFile.setLength(0L);
        }
        this.bufferByte = new byte[this.headerSize + 4];
        this.fileInfo.setSizeOfHeader(this.headerSize);
        final int[] extents = myFileInfo.getExtents();
        final int[] unitsOfMeasure = myFileInfo.getUnitsOfMeasure();
        final float[] resols = myFileInfo.getResolutions();
        this.origin = myFileInfo.getOrigin();
        this.axisOrientation = myFileInfo.getAxisOrientation();
        boolean found = false;
        int firstSpatialDim = -1;
        for (int i = 0; i < unitsOfMeasure.length && !found; ++i) {
            if (unitsOfMeasure[i] == FileInfoBase.Unit.UNKNOWN_MEASURE.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.INCHES.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MILS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.CENTIMETERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.ANGSTROMS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.NANOMETERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MICROMETERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MILLIMETERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.METERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.KILOMETERS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MILES.getLegacyNum()) {
                found = true;
                firstSpatialDim = i;
                firstSpatialUnits = unitsOfMeasure[i];
            }
        }
        found = false;
        int firstTimeDim = -1;
        for (int i = 0; i < unitsOfMeasure.length && !found; ++i) {
            if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOSEC.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MICROSEC.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MILLISEC.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.SECONDS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.MINUTES.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.HOURS.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.HZ.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.PPM.getLegacyNum() || unitsOfMeasure[i] == FileInfoBase.Unit.RADS.getLegacyNum()) {
                found = true;
                firstTimeDim = i;
                firstTimeUnits = unitsOfMeasure[i];
            }
        }
        if (firstSpatialDim >= 0) {
            switch (FileInfoBase.Unit.getUnitFromLegacyNum(firstSpatialUnits)) {
                case UNKNOWN_MEASURE: {
                    niftiSpatialUnits = 0;
                    break;
                }
                case INCHES:
                case MILS:
                case CENTIMETERS:
                case MILLIMETERS: {
                    niftiSpatialUnits = 2;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.INCHES.getLegacyNum()) {
                            resols[i] *= 25.4f;
                            this.origin[i] *= 25.4f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILS.getLegacyNum()) {
                            resols[i] *= 0.0254f;
                            this.origin[i] *= 0.0254f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.CENTIMETERS.getLegacyNum()) {
                            resols[i] *= 10.0f;
                            this.origin[i] *= 10.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.ANGSTROMS.getLegacyNum()) {
                            resols[i] *= 1.0E-7f;
                            this.origin[i] *= 1.0E-7f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOMETERS.getLegacyNum()) {
                            resols[i] *= 1.0E-6f;
                            this.origin[i] *= 1.0E-6f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MICROMETERS.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.METERS.getLegacyNum()) {
                            resols[i] *= 1000.0f;
                            this.origin[i] *= 1000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.KILOMETERS.getLegacyNum()) {
                            resols[i] *= 1000000.0f;
                            this.origin[i] *= 1000000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILES.getLegacyNum()) {
                            resols[i] *= 1609300.0f;
                            this.origin[i] *= 1609300.0f;
                        }
                    }
                    break;
                }
                case ANGSTROMS:
                case NANOMETERS:
                case MICROMETERS: {
                    niftiSpatialUnits = 3;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.INCHES.getLegacyNum()) {
                            resols[i] *= 254000.0f;
                            this.origin[i] *= 254000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILS.getLegacyNum()) {
                            resols[i] *= 254.0f;
                            this.origin[i] *= 254.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.CENTIMETERS.getLegacyNum()) {
                            resols[i] *= 100000.0f;
                            this.origin[i] *= 100000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.ANGSTROMS.getLegacyNum()) {
                            resols[i] *= 1.0E-4f;
                            this.origin[i] *= 1.0E-4f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOMETERS.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILLIMETERS.getLegacyNum()) {
                            resols[i] *= 1000.0f;
                            this.origin[i] *= 1000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.METERS.getLegacyNum()) {
                            resols[i] *= 1000000.0f;
                            this.origin[i] *= 1000000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.KILOMETERS.getLegacyNum()) {
                            resols[i] *= 1.0E9f;
                            this.origin[i] *= 1.0E9f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILES.getLegacyNum()) {
                            resols[i] *= 1.60929997E9f;
                            this.origin[i] *= 1.60929997E9f;
                        }
                    }
                    break;
                }
                case METERS:
                case KILOMETERS:
                case MILES: {
                    niftiSpatialUnits = 1;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.INCHES.getLegacyNum()) {
                            resols[i] *= 0.0254f;
                            this.origin[i] *= 0.0254f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILS.getLegacyNum()) {
                            resols[i] *= 2.54E-5f;
                            this.origin[i] *= 2.54E-5f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.CENTIMETERS.getLegacyNum()) {
                            resols[i] *= 0.01f;
                            this.origin[i] *= 0.01f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.ANGSTROMS.getLegacyNum()) {
                            resols[i] *= 1.0E-10f;
                            this.origin[i] *= 1.0E-10f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOMETERS.getLegacyNum()) {
                            resols[i] *= 1.0E-9f;
                            this.origin[i] *= 1.0E-9f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MICROMETERS.getLegacyNum()) {
                            resols[i] *= 1.0E-6f;
                            this.origin[i] *= 1.0E-6f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILLIMETERS.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.KILOMETERS.getLegacyNum()) {
                            resols[i] *= 1000.0f;
                            this.origin[i] *= 1000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILES.getLegacyNum()) {
                            resols[i] *= 1609.3f;
                            this.origin[i] *= 1609.3f;
                        }
                    }
                    break;
                }
            }
        }
        if (firstTimeDim >= 0) {
            switch (FileInfoBase.Unit.getUnitFromLegacyNum(firstTimeUnits)) {
                case NANOSEC:
                case MICROSEC: {
                    niftiTimeUnits = 24;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOSEC.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILLISEC.getLegacyNum()) {
                            resols[i] *= 1000.0f;
                            this.origin[i] *= 1000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.SECONDS.getLegacyNum()) {
                            resols[i] *= 1000000.0f;
                            this.origin[i] *= 1000000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MINUTES.getLegacyNum()) {
                            resols[i] *= 6.0E7f;
                            this.origin[i] *= 6.0E7f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.HOURS.getLegacyNum()) {
                            resols[i] *= 3.6E9f;
                            this.origin[i] *= 3.6E9f;
                        }
                    }
                    break;
                }
                case MILLISEC: {
                    niftiTimeUnits = 16;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOSEC.getLegacyNum()) {
                            resols[i] *= 1.0E-6f;
                            this.origin[i] *= 1.0E-6f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MICROSEC.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.SECONDS.getLegacyNum()) {
                            resols[i] *= 1000.0f;
                            this.origin[i] *= 1000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MINUTES.getLegacyNum()) {
                            resols[i] *= 60000.0f;
                            this.origin[i] *= 60000.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.HOURS.getLegacyNum()) {
                            resols[i] *= 3600000.0f;
                            this.origin[i] *= 3600000.0f;
                        }
                    }
                    break;
                }
                case SECONDS:
                case MINUTES:
                case HOURS: {
                    niftiTimeUnits = 8;
                    for (int i = 0; i < extents.length; ++i) {
                        if (unitsOfMeasure[i] == FileInfoBase.Unit.NANOSEC.getLegacyNum()) {
                            resols[i] *= 1.0E-9f;
                            this.origin[i] *= 1.0E-9f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MICROSEC.getLegacyNum()) {
                            resols[i] *= 1.0E-6f;
                            this.origin[i] *= 1.0E-6f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MILLISEC.getLegacyNum()) {
                            resols[i] *= 0.001f;
                            this.origin[i] *= 0.001f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.MINUTES.getLegacyNum()) {
                            resols[i] *= 60.0f;
                            this.origin[i] *= 60.0f;
                        }
                        else if (unitsOfMeasure[i] == FileInfoBase.Unit.HOURS.getLegacyNum()) {
                            resols[i] *= 3600.0f;
                            this.origin[i] *= 3600.0f;
                        }
                    }
                    break;
                }
                case HZ: {
                    niftiTimeUnits = 32;
                    break;
                }
                case PPM: {
                    niftiTimeUnits = 40;
                    break;
                }
                case RADS: {
                    niftiTimeUnits = 48;
                    break;
                }
            }
        }
        int nDims;
        int nDimsLength1;
        if (firstTimeDim < 0) {
            nDims = extents.length;
            nDimsLength1 = 0;
        }
        else {
            nDims = extents.length + 3 - firstTimeDim;
            nDimsLength1 = 3 - firstTimeDim;
        }
        Preferences.debug("FileNIFTI:writeHeader - nImagesSaved = " + nImagesSaved + "\n", 2);
        Preferences.debug("FileNIFTI:writeHeader - nDims = " + nDims + "\n", 2);
        final int[] niftiExtents = new int[nDims + 1];
        niftiExtents[0] = nDims;
        for (int i = 1; i <= nDimsLength1; ++i) {
            niftiExtents[i] = 1;
        }
        for (int i = 1; i <= nDims - nDimsLength1; ++i) {
            if (i == 3) {
                niftiExtents[i] = nImagesSaved;
            }
            else if (i == 4) {
                niftiExtents[i] = nTimeSaved;
            }
            else {
                niftiExtents[i] = extents[i - 1];
            }
        }
        for (int i = 0; i < niftiExtents.length; ++i) {
            Preferences.debug("FileNIFTI:writeHeader - i = " + i + " dim = " + niftiExtents[i] + "\n", 2);
        }
        switch (image.getType()) {
            case 0: {
                this.sourceType = 1;
                break;
            }
            case 1: {
                this.sourceType = 256;
                break;
            }
            case 2: {
                this.sourceType = 2;
                break;
            }
            case 3: {
                this.sourceType = 4;
                break;
            }
            case 4: {
                this.sourceType = 512;
                break;
            }
            case 5: {
                this.sourceType = 8;
                break;
            }
            case 14: {
                this.sourceType = 768;
                break;
            }
            case 6: {
                this.sourceType = 1024;
                break;
            }
            case 7: {
                this.sourceType = 16;
                break;
            }
            case 8: {
                this.sourceType = 64;
                break;
            }
            case 9: {
                this.sourceType = 128;
                break;
            }
            case 12: {
                this.sourceType = 32;
                break;
            }
            case 13: {
                this.sourceType = 1792;
                break;
            }
            default: {
                return false;
            }
        }
        this.fileInfo.setSourceType(this.sourceType);
        matHolder = image.getMatrixHolder();
        if (matHolder != null) {
            matrixArray = matHolder.getNIFTICompositeMatrices();
            if (matrixArray != null) {
                if (matrixArray.length >= 1 && matrixArray[0] != null) {
                    if (matrixArray[0].isQform()) {
                        matrixQ = matrixArray[0];
                        transformIDQ = matrixArray[0].getTransformID();
                    }
                    else {
                        matrixS = matrixArray[0];
                        transformIDS = matrixArray[0].getTransformID();
                    }
                }
                if (matrixArray.length >= 2 && matrixArray[1] != null) {
                    if (matrixArray[1].isQform()) {
                        matrixQ = matrixArray[1];
                        transformIDQ = matrixArray[1].getTransformID();
                    }
                    else {
                        matrixS = matrixArray[1];
                        transformIDS = matrixArray[1].getTransformID();
                    }
                }
            }
        }
        if (matrixQ == null && matrixS == null && isNIFTI) {
            matrixQ = image.getMatrix();
            transformIDQ = matrixQ.getTransformID();
        }
        if (matrixQ != null || (matrixQ == null && matrixS == null)) {
            double r00;
            double r2;
            double r3;
            double r4;
            double r5;
            double r6;
            double r7;
            double r8;
            double r9;
            if (matrixQ != null) {
                Preferences.debug("matrixQ on write entry = " + matrixQ + "\n", 2);
                switch (transformIDQ) {
                    case 6: {
                        qform_code = 1;
                        break;
                    }
                    case 2: {
                        qform_code = 2;
                        break;
                    }
                    case 3: {
                        qform_code = 3;
                        break;
                    }
                    case 4: {
                        qform_code = 4;
                        break;
                    }
                    default: {
                        qform_code = 1;
                        break;
                    }
                }
                if (image.getNDims() >= 3) {
                    this.axisOrientation = this.getAxisOrientation(matrixQ);
                    r00 = -matrixQ.get(0, 0) / resols[0];
                    r2 = -matrixQ.get(0, 1) / resols[1];
                    r3 = -matrixQ.get(0, 2) / resols[2];
                    r4 = -matrixQ.get(1, 0) / resols[0];
                    r5 = -matrixQ.get(1, 1) / resols[1];
                    r6 = -matrixQ.get(1, 2) / resols[2];
                    r7 = matrixQ.get(2, 0) / resols[0];
                    r8 = matrixQ.get(2, 1) / resols[1];
                    r9 = matrixQ.get(2, 2) / resols[2];
                }
                else {
                    if (this.origin != null && this.origin.length > 2) {
                        niftiOrigin[2] = this.origin[2];
                    }
                    if (this.axisOrientation[0] != 0 && this.axisOrientation[1] != 0 && this.axisOrientation[2] == 0) {
                        if (this.axisOrientation[0] != 1 && this.axisOrientation[0] != 2 && this.axisOrientation[1] != 1 && this.axisOrientation[1] != 2) {
                            this.axisOrientation[2] = 1;
                        }
                        else if (this.axisOrientation[0] != 4 && this.axisOrientation[0] != 3 && this.axisOrientation[1] != 4 && this.axisOrientation[1] != 3) {
                            this.axisOrientation[2] = 4;
                        }
                        else {
                            this.axisOrientation[2] = 5;
                        }
                    }
                    if (this.axisOrientation[0] == 0 || this.axisOrientation[1] == 0) {
                        r00 = -1.0;
                        r2 = 0.0;
                        r3 = 0.0;
                        r4 = 0.0;
                        r5 = -1.0;
                        r6 = 0.0;
                        r7 = 0.0;
                        r8 = 0.0;
                        r9 = 1.0;
                    }
                    else {
                        r00 = 0.0;
                        r2 = 0.0;
                        r3 = 0.0;
                        r4 = 0.0;
                        r5 = 0.0;
                        r6 = 0.0;
                        r7 = 0.0;
                        r8 = 0.0;
                        r9 = 0.0;
                        switch (this.axisOrientation[0]) {
                            case 1: {
                                r00 = -1.0;
                                break;
                            }
                            case 2: {
                                r00 = 1.0;
                                break;
                            }
                            case 4: {
                                r4 = -1.0;
                                break;
                            }
                            case 3: {
                                r4 = 1.0;
                                break;
                            }
                            case 5: {
                                r7 = 1.0;
                                break;
                            }
                            case 6: {
                                r7 = -1.0;
                                break;
                            }
                        }
                        switch (this.axisOrientation[1]) {
                            case 1: {
                                r2 = -1.0;
                                break;
                            }
                            case 2: {
                                r2 = 1.0;
                                break;
                            }
                            case 4: {
                                r5 = -1.0;
                                break;
                            }
                            case 3: {
                                r5 = 1.0;
                                break;
                            }
                            case 5: {
                                r8 = 1.0;
                                break;
                            }
                            case 6: {
                                r8 = -1.0;
                                break;
                            }
                        }
                        switch (this.axisOrientation[2]) {
                            case 1: {
                                r3 = -1.0;
                                break;
                            }
                            case 2: {
                                r3 = 1.0;
                                break;
                            }
                            case 4: {
                                r6 = -1.0;
                                break;
                            }
                            case 3: {
                                r6 = 1.0;
                                break;
                            }
                            case 5: {
                                r9 = 1.0;
                                break;
                            }
                            case 6: {
                                r9 = -1.0;
                                break;
                            }
                        }
                    }
                }
                for (int j = 0; j < Math.min(3, image.getNDims()); ++j) {
                    if (this.axisOrientation[j] == 2) {
                        niftiOrigin[0] = -Math.abs(matrixQ.get(j, 3));
                    }
                    else if (this.axisOrientation[j] == 1) {
                        niftiOrigin[0] = Math.abs(matrixQ.get(j, 3));
                    }
                    else if (this.axisOrientation[j] == 3) {
                        niftiOrigin[1] = -Math.abs(matrixQ.get(j, 3));
                    }
                    else if (this.axisOrientation[j] == 4) {
                        niftiOrigin[1] = Math.abs(matrixQ.get(j, 3));
                    }
                    else if (this.axisOrientation[j] == 5) {
                        niftiOrigin[2] = -Math.abs(matrixQ.get(j, 3));
                    }
                    else if (this.axisOrientation[j] == 6) {
                        niftiOrigin[2] = Math.abs(matrixQ.get(j, 3));
                    }
                }
            }
            else if (dicomMatrix != null) {
                final TransMatrix transposeMatrix = new TransMatrix(4);
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        transposeMatrix.set(i, j, dicomMatrix.get(j, i));
                    }
                }
                if ((transposeMatrix.get(0, 2) <= 0.0f && xDel > 0.0) || (transposeMatrix.get(0, 2) >= 0.0f && xDel < 0.0)) {
                    transposeMatrix.set(0, 2, -transposeMatrix.get(0, 2));
                }
                if ((transposeMatrix.get(1, 2) <= 0.0f && yDel > 0.0) || (transposeMatrix.get(1, 2) >= 0.0f && yDel < 0.0)) {
                    transposeMatrix.set(1, 2, -transposeMatrix.get(1, 2));
                }
                if ((transposeMatrix.get(2, 2) <= 0.0f && zDel > 0.0) || (transposeMatrix.get(2, 2) >= 0.0f && zDel < 0.0)) {
                    transposeMatrix.set(2, 2, -transposeMatrix.get(2, 2));
                }
                this.axisOrientation = this.getAxisOrientation(transposeMatrix);
                qform_code = 1;
                r00 = -transposeMatrix.get(0, 0);
                r2 = -transposeMatrix.get(0, 1);
                r3 = -transposeMatrix.get(0, 2);
                r4 = -transposeMatrix.get(1, 0);
                r5 = -transposeMatrix.get(1, 1);
                r6 = -transposeMatrix.get(1, 2);
                r7 = transposeMatrix.get(2, 0);
                r8 = transposeMatrix.get(2, 1);
                r9 = transposeMatrix.get(2, 2);
                for (int j = 0; j < 3; ++j) {
                    if (this.axisOrientation[j] == 2) {
                        niftiOrigin[0] = -Math.abs(this.origin[j]);
                    }
                    else if (this.axisOrientation[j] == 1) {
                        niftiOrigin[0] = Math.abs(this.origin[j]);
                    }
                    else if (this.axisOrientation[j] == 3) {
                        niftiOrigin[1] = -Math.abs(this.origin[j]);
                    }
                    else if (this.axisOrientation[j] == 4) {
                        niftiOrigin[1] = Math.abs(this.origin[j]);
                    }
                    else if (this.axisOrientation[j] == 5) {
                        niftiOrigin[2] = -Math.abs(this.origin[j]);
                    }
                    else if (this.axisOrientation[j] == 6) {
                        niftiOrigin[2] = Math.abs(this.origin[j]);
                    }
                }
            }
            else {
                qform_code = 1;
                if (this.axisOrientation != null && this.axisOrientation[0] != 0 && this.axisOrientation[1] != 0) {
                    if (this.axisOrientation[2] == 0) {
                        if (this.axisOrientation[0] != 1 && this.axisOrientation[0] != 2 && this.axisOrientation[1] != 1 && this.axisOrientation[1] != 2) {
                            this.axisOrientation[2] = 1;
                        }
                        else if (this.axisOrientation[0] != 4 && this.axisOrientation[0] != 3 && this.axisOrientation[1] != 4 && this.axisOrientation[1] != 3) {
                            this.axisOrientation[2] = 4;
                        }
                        else {
                            this.axisOrientation[2] = 5;
                        }
                    }
                    r00 = 0.0;
                    r2 = 0.0;
                    r3 = 0.0;
                    r4 = 0.0;
                    r5 = 0.0;
                    r6 = 0.0;
                    r7 = 0.0;
                    r8 = 0.0;
                    r9 = 0.0;
                    switch (this.axisOrientation[0]) {
                        case 1: {
                            r00 = -1.0;
                            break;
                        }
                        case 2: {
                            r00 = 1.0;
                            break;
                        }
                        case 4: {
                            r4 = -1.0;
                            break;
                        }
                        case 3: {
                            r4 = 1.0;
                            break;
                        }
                        case 5: {
                            r7 = 1.0;
                            break;
                        }
                        case 6: {
                            r7 = -1.0;
                            break;
                        }
                    }
                    switch (this.axisOrientation[1]) {
                        case 1: {
                            r2 = -1.0;
                            break;
                        }
                        case 2: {
                            r2 = 1.0;
                            break;
                        }
                        case 4: {
                            r5 = -1.0;
                            break;
                        }
                        case 3: {
                            r5 = 1.0;
                            break;
                        }
                        case 5: {
                            r8 = 1.0;
                            break;
                        }
                        case 6: {
                            r8 = -1.0;
                            break;
                        }
                    }
                    switch (this.axisOrientation[2]) {
                        case 1: {
                            r3 = -1.0;
                            break;
                        }
                        case 2: {
                            r3 = 1.0;
                            break;
                        }
                        case 4: {
                            r6 = -1.0;
                            break;
                        }
                        case 3: {
                            r6 = 1.0;
                            break;
                        }
                        case 5: {
                            r9 = 1.0;
                            break;
                        }
                        case 6: {
                            r9 = -1.0;
                            break;
                        }
                    }
                    for (int j = 0; j < 3; ++j) {
                        if (this.axisOrientation[j] == 2) {
                            niftiOrigin[0] = -Math.abs(this.origin[j]);
                        }
                        else if (this.axisOrientation[j] == 1) {
                            niftiOrigin[0] = Math.abs(this.origin[j]);
                        }
                        else if (this.axisOrientation[j] == 3) {
                            niftiOrigin[1] = -Math.abs(this.origin[j]);
                        }
                        else if (this.axisOrientation[j] == 4) {
                            niftiOrigin[1] = Math.abs(this.origin[j]);
                        }
                        else if (this.axisOrientation[j] == 5) {
                            niftiOrigin[2] = -Math.abs(this.origin[j]);
                        }
                        else if (this.axisOrientation[j] == 6) {
                            niftiOrigin[2] = Math.abs(this.origin[j]);
                        }
                    }
                }
                else {
                    r00 = -1.0;
                    r2 = 0.0;
                    r3 = 0.0;
                    r4 = 0.0;
                    r5 = -1.0;
                    r6 = 0.0;
                    r7 = 0.0;
                    r8 = 0.0;
                    r9 = 1.0;
                }
            }
            double xd = Math.sqrt(r00 * r00 + r4 * r4 + r7 * r7);
            double yd = Math.sqrt(r2 * r2 + r5 * r5 + r8 * r8);
            double zd = Math.sqrt(r3 * r3 + r6 * r6 + r9 * r9);
            if (xd == 0.0) {
                r00 = 1.0;
                r4 = 0.0;
                r7 = 0.0;
                xd = 1.0;
            }
            if (yd == 0.0) {
                r2 = 0.0;
                r5 = 1.0;
                r8 = 0.0;
                yd = 1.0;
            }
            if (zd == 0.0) {
                r3 = 0.0;
                r6 = 0.0;
                r9 = 1.0;
                zd = 1.0;
            }
            r00 /= xd;
            r4 /= xd;
            r7 /= xd;
            r2 /= yd;
            r5 /= yd;
            r8 /= yd;
            r3 /= zd;
            r6 /= zd;
            r9 /= zd;
            final Matrix Q = new Matrix(3, 3);
            Q.set(0, 0, r00);
            Q.set(0, 1, r2);
            Q.set(0, 2, r3);
            Q.set(1, 0, r4);
            Q.set(1, 1, r5);
            Q.set(1, 2, r6);
            Q.set(2, 0, r7);
            Q.set(2, 1, r8);
            Q.set(2, 2, r9);
            final Matrix P = this.mat33_polar(Q);
            r00 = P.get(0, 0);
            r2 = P.get(0, 1);
            r3 = P.get(0, 2);
            r4 = P.get(1, 0);
            r5 = P.get(1, 1);
            r6 = P.get(1, 2);
            r7 = P.get(2, 0);
            r8 = P.get(2, 1);
            r9 = P.get(2, 2);
            zd = P.det();
            if (zd > 0.0) {
                this.qfac = 1.0f;
            }
            else {
                this.qfac = -1.0f;
                r3 = -r3;
                r6 = -r6;
                r9 = -r9;
            }
            double a = r00 + r5 + r9 + 1.0;
            double b;
            double c;
            double d;
            if (a > 0.5) {
                a = 0.5 * Math.sqrt(a);
                b = 0.25 * (r8 - r6) / a;
                c = 0.25 * (r3 - r7) / a;
                d = 0.25 * (r4 - r2) / a;
            }
            else {
                xd = 1.0 + r00 - (r5 + r9);
                yd = 1.0 + r5 - (r00 + r9);
                zd = 1.0 + r9 - (r00 + r5);
                if (xd > 1.0) {
                    b = 0.5 * Math.sqrt(xd);
                    c = 0.25 * (r2 + r4) / b;
                    d = 0.25 * (r3 + r7) / b;
                    a = 0.25 * (r8 - r6) / b;
                }
                else if (yd > 1.0) {
                    c = 0.5 * Math.sqrt(yd);
                    b = 0.25 * (r2 + r4) / c;
                    d = 0.25 * (r6 + r8) / c;
                    a = 0.25 * (r3 - r7) / c;
                }
                else {
                    d = 0.5 * Math.sqrt(zd);
                    b = 0.25 * (r3 + r7) / d;
                    c = 0.25 * (r6 + r8) / d;
                    a = 0.25 * (r4 - r2) / d;
                }
                if (a < 0.0) {
                    a = -a;
                    b = -b;
                    c = -c;
                    d = -d;
                }
            }
            this.quatern_a = (float)a;
            this.quatern_b = (float)b;
            this.quatern_c = (float)c;
            this.quatern_d = (float)d;
        }
        if (matrixS != null) {
            Preferences.debug("matrixS on write entry = " + matrixS + "\n", 2);
            switch (transformIDS) {
                case 6: {
                    sform_code = 1;
                    break;
                }
                case 2: {
                    sform_code = 2;
                    break;
                }
                case 3: {
                    sform_code = 3;
                    break;
                }
                case 4: {
                    sform_code = 4;
                    break;
                }
                default: {
                    sform_code = 1;
                    break;
                }
            }
            niftiOriginS = new float[3];
            this.axisOrientation = this.getAxisOrientation(matrixS);
            for (int j = 0; j < 3; ++j) {
                if (this.axisOrientation[j] == 2) {
                    niftiOriginS[0] = -Math.abs(matrixS.get(j, 3));
                }
                else if (this.axisOrientation[j] == 1) {
                    niftiOriginS[0] = Math.abs(matrixS.get(j, 3));
                }
                else if (this.axisOrientation[j] == 3) {
                    niftiOriginS[1] = -Math.abs(matrixS.get(j, 3));
                }
                else if (this.axisOrientation[j] == 4) {
                    niftiOriginS[1] = Math.abs(matrixS.get(j, 3));
                }
                else if (this.axisOrientation[j] == 5) {
                    niftiOriginS[2] = -Math.abs(matrixS.get(j, 3));
                }
                else if (this.axisOrientation[j] == 6) {
                    niftiOriginS[2] = Math.abs(matrixS.get(j, 3));
                }
            }
        }
        if (isNIFTI) {
            this.setBufferInt(this.bufferByte, this.fileInfo.getSizeOfHeader(), 0, endianess);
            this.setBufferString(this.bufferByte, fileName, 14);
            this.setBufferInt(this.bufferByte, 0, 32, endianess);
            this.bufferByte[38] = 114;
            this.freq_dim = this.fileInfo.getFreqDim();
            this.phase_dim = this.fileInfo.getPhaseDim();
            this.slice_dim = this.fileInfo.getSliceDim();
            final byte dim_info = (byte)(this.freq_dim | this.phase_dim << 2 | this.slice_dim << 4);
            this.bufferByte[39] = dim_info;
            for (int i = 0; i < niftiExtents.length; ++i) {
                this.setBufferShort(this.bufferByte, (short)niftiExtents[i], 40 + i * 2, endianess);
            }
            for (int i = niftiExtents.length; i < 8; ++i) {
                this.setBufferShort(this.bufferByte, (short)1, 40 + i * 2, endianess);
            }
            this.setBufferFloat(this.bufferByte, this.fileInfo.getIntentP1(), 56, endianess);
            this.setBufferFloat(this.bufferByte, this.fileInfo.getIntentP2(), 60, endianess);
            this.setBufferFloat(this.bufferByte, this.fileInfo.getIntentP3(), 64, endianess);
            this.setBufferShort(this.bufferByte, this.fileInfo.getIntentCode(), 68, endianess);
            this.setBufferShort(this.bufferByte, this.sourceType, 70, endianess);
            this.setBufferShort(this.bufferByte, this.fileInfo.getBitPix(), 72, endianess);
            this.setBufferShort(this.bufferByte, this.fileInfo.getSliceStart(), 74, endianess);
            this.setBufferFloat(this.bufferByte, this.qfac, 76, endianess);
            final float[] niftiResols = new float[nDims];
            for (int i = 0; i < nDimsLength1; ++i) {
                niftiResols[i] = 1.0f;
            }
            for (int i = 0; i < nDims - nDimsLength1; ++i) {
                niftiResols[i + nDimsLength1] = resols[i];
            }
            for (int i = 0; i < nDims; ++i) {
                this.setBufferFloat(this.bufferByte, niftiResols[i], 80 + i * 4, endianess);
            }
            for (int i = nDims; i < 7; ++i) {
                this.setBufferFloat(this.bufferByte, 0.0f, 80 + i * 4, endianess);
            }
            if (oneFile) {
                this.setBufferFloat(this.bufferByte, 352.0f, 108, endianess);
            }
            else {
                this.setBufferFloat(this.bufferByte, 0.0f, 108, endianess);
            }
            this.setBufferFloat(this.bufferByte, 1.0f, 112, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 116, endianess);
            this.setBufferShort(this.bufferByte, this.fileInfo.getSliceEnd(), 120, endianess);
            this.bufferByte[122] = this.fileInfo.getSliceCode();
            this.bufferByte[123] = (byte)(niftiSpatialUnits | niftiTimeUnits);
            this.setBufferFloat(this.bufferByte, this.fileInfo.getCalMax(), 124, endianess);
            this.setBufferFloat(this.bufferByte, this.fileInfo.getCalMin(), 128, endianess);
            this.setBufferFloat(this.bufferByte, this.fileInfo.getSliceDuration(), 132, endianess);
            if (this.origin.length >= 4) {
                this.setBufferFloat(this.bufferByte, this.origin[3], 136, endianess);
            }
            else {
                this.setBufferFloat(this.bufferByte, 0.0f, 136, endianess);
            }
            if (image.getType() == 12 || image.getType() == 13) {
                image.calcMinMaxMag(Preferences.is("LogMagDisplay"));
            }
            else {
                image.calcMinMax();
            }
            int imageMax;
            int imageMin;
            if (image.isColorImage()) {
                imageMax = (int)Math.max(image.getMaxR(), Math.max(image.getMaxG(), image.getMaxB()));
                imageMin = (int)Math.min(image.getMinR(), Math.min(image.getMinG(), image.getMinB()));
            }
            else {
                imageMax = (int)image.getMax();
                imageMin = (int)image.getMin();
            }
            this.setBufferInt(this.bufferByte, imageMax, 140, endianess);
            this.setBufferInt(this.bufferByte, imageMin, 144, endianess);
            final int modality = this.fileInfo.getModality();
            if (modality != 0) {
                this.fileInfo.setDescription(FileInfoBase.getModalityStr(modality));
            }
            if (this.fileInfo.getDescription() != null) {
                this.setBufferString(this.bufferByte, this.fileInfo.getDescription(), 148);
            }
            if (this.fileInfo.getAuxFile() != null) {
                this.setBufferString(this.bufferByte, this.fileInfo.getAuxFile(), 228);
            }
            this.setBufferShort(this.bufferByte, (short)qform_code, 252, endianess);
            this.setBufferShort(this.bufferByte, (short)sform_code, 254, endianess);
            if (qform_code > 0) {
                Preferences.debug("Writing quatern_b = " + this.quatern_b + "\n", 2);
                this.setBufferFloat(this.bufferByte, this.quatern_b, 256, endianess);
                Preferences.debug("Writing quatern_c = " + this.quatern_c + "\n", 2);
                this.setBufferFloat(this.bufferByte, this.quatern_c, 260, endianess);
                Preferences.debug("Writing quatern_d = " + this.quatern_d + "\n");
                this.setBufferFloat(this.bufferByte, this.quatern_d, 264, endianess);
                this.setBufferFloat(this.bufferByte, niftiOrigin[0], 268, endianess);
                this.setBufferFloat(this.bufferByte, niftiOrigin[1], 272, endianess);
                this.setBufferFloat(this.bufferByte, niftiOrigin[2], 276, endianess);
            }
            if (matrixS != null) {
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 0), 280, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 1), 284, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 2), 288, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[0], 292, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 0), 296, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 1), 300, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 2), 304, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[1], 308, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 0), 312, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 1), 316, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 2), 320, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[2], 324, endianess);
            }
            if (this.fileInfo.getIntentName() != null) {
                this.setBufferString(this.bufferByte, this.fileInfo.getIntentName(), 328);
            }
        }
        else {
            this.setBufferInt(this.bufferByte, this.fileInfo.getSizeOfHeader(), 0, endianess);
            this.setBufferString(this.bufferByte, "         \n", 4);
            this.setBufferString(this.bufferByte, fileName + "\n", 14);
            this.setBufferInt(this.bufferByte, 0, 32, endianess);
            this.setBufferShort(this.bufferByte, (short)0, 36, endianess);
            this.bufferByte[38] = 114;
            this.bufferByte[39] = 0;
            for (int i = 0; i < niftiExtents.length; ++i) {
                this.setBufferShort(this.bufferByte, (short)niftiExtents[i], 40 + i * 2, endianess);
            }
            for (int i = niftiExtents.length; i < 8; ++i) {
                this.setBufferShort(this.bufferByte, (short)1, 40 + i * 2, endianess);
            }
            final int[] units = myFileInfo.getUnitsOfMeasure();
            this.fileInfo.setUnitsOfMeasure(units);
            this.setBufferFloat(this.bufferByte, 0.0f, 56, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 60, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 64, endianess);
            this.setBufferShort(this.bufferByte, (short)0, 68, endianess);
            Preferences.debug("FileNIFTI:writeHeader(simple): data type = " + this.sourceType + "\n", 2);
            this.setBufferShort(this.bufferByte, this.sourceType, 70, endianess);
            switch (image.getType()) {
                case 0: {
                    this.fileInfo.setBitPix((short)1);
                    break;
                }
                case 1: {
                    this.fileInfo.setBitPix((short)8);
                    break;
                }
                case 2: {
                    this.fileInfo.setBitPix((short)8);
                    break;
                }
                case 3: {
                    this.fileInfo.setBitPix((short)16);
                    break;
                }
                case 4: {
                    this.fileInfo.setBitPix((short)16);
                    break;
                }
                case 5: {
                    this.fileInfo.setBitPix((short)32);
                    break;
                }
                case 6: {
                    this.fileInfo.setBitPix((short)64);
                    break;
                }
                case 7: {
                    this.fileInfo.setBitPix((short)32);
                    break;
                }
                case 8: {
                    this.fileInfo.setBitPix((short)64);
                    break;
                }
                case 9: {
                    this.fileInfo.setBitPix((short)24);
                    break;
                }
                case 12: {
                    this.fileInfo.setBitPix((short)64);
                    break;
                }
                case 13: {
                    this.fileInfo.setBitPix((short)128);
                    break;
                }
                default: {
                    return false;
                }
            }
            Preferences.debug("FileNIFTI:writeHeader(simple): bits per pixel = " + this.fileInfo.getBitPix() + "\n", 2);
            this.setBufferShort(this.bufferByte, this.fileInfo.getBitPix(), 72, endianess);
            this.setBufferShort(this.bufferByte, (short)0, 74, endianess);
            this.setBufferFloat(this.bufferByte, this.qfac, 76, endianess);
            final float[] niftiResols = new float[nDims];
            for (int i = 0; i < nDimsLength1; ++i) {
                niftiResols[i] = 1.0f;
            }
            for (int i = 0; i < nDims - nDimsLength1; ++i) {
                niftiResols[i + nDimsLength1] = resols[i];
            }
            for (int i = 0; i < nDims; ++i) {
                this.setBufferFloat(this.bufferByte, niftiResols[i], 80 + i * 4, endianess);
            }
            for (int i = nDims; i < 7; ++i) {
                this.setBufferFloat(this.bufferByte, 0.0f, 80 + i * 4, endianess);
            }
            if (oneFile) {
                this.setBufferFloat(this.bufferByte, 352.0f, 108, endianess);
            }
            else {
                this.setBufferFloat(this.bufferByte, 0.0f, 108, endianess);
            }
            this.setBufferFloat(this.bufferByte, 1.0f, 112, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 116, endianess);
            this.setBufferShort(this.bufferByte, (short)0, 120, endianess);
            this.bufferByte[122] = 0;
            this.bufferByte[123] = (byte)(niftiSpatialUnits | niftiTimeUnits);
            this.setBufferFloat(this.bufferByte, 0.0f, 124, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 128, endianess);
            this.setBufferFloat(this.bufferByte, 0.0f, 132, endianess);
            if (this.origin.length >= 4) {
                this.setBufferFloat(this.bufferByte, this.origin[3], 136, endianess);
            }
            else {
                this.setBufferFloat(this.bufferByte, 0.0f, 136, endianess);
            }
            if (image.getType() == 12 || image.getType() == 13) {
                image.calcMinMaxMag(Preferences.is("LogMagDisplay"));
            }
            else {
                image.calcMinMax();
            }
            int imageMax;
            int imageMin;
            if (image.isColorImage()) {
                imageMax = (int)Math.max(image.getMaxR(), Math.max(image.getMaxG(), image.getMaxB()));
                imageMin = (int)Math.min(image.getMinR(), Math.min(image.getMinG(), image.getMinB()));
            }
            else {
                imageMax = (int)image.getMax();
                imageMin = (int)image.getMin();
            }
            this.setBufferInt(this.bufferByte, imageMax, 140, endianess);
            this.setBufferInt(this.bufferByte, imageMin, 144, endianess);
            final int modality2 = myFileInfo.getModality();
            this.fileInfo.setModality(modality2);
            this.setBufferString(this.bufferByte, FileInfoBase.getModalityStr(modality2), 148);
            this.setBufferString(this.bufferByte, " ", 228);
            this.setBufferShort(this.bufferByte, (short)qform_code, 252, endianess);
            this.setBufferShort(this.bufferByte, (short)sform_code, 254, endianess);
            Preferences.debug("Writing quatern_b = " + this.quatern_b + "\n", 2);
            this.setBufferFloat(this.bufferByte, this.quatern_b, 256, endianess);
            Preferences.debug("Writing quatern_c = " + this.quatern_c + "\n", 2);
            this.setBufferFloat(this.bufferByte, this.quatern_c, 260, endianess);
            Preferences.debug("Writing quatern_d = " + this.quatern_d + "\n", 2);
            this.setBufferFloat(this.bufferByte, this.quatern_d, 264, endianess);
            this.setBufferFloat(this.bufferByte, niftiOrigin[0], 268, endianess);
            this.setBufferFloat(this.bufferByte, niftiOrigin[1], 272, endianess);
            this.setBufferFloat(this.bufferByte, niftiOrigin[2], 276, endianess);
            if (matrixS != null) {
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 0), 280, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 1), 284, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(0, 2), 288, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[0], 292, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 0), 296, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 1), 300, endianess);
                this.setBufferFloat(this.bufferByte, -matrixS.get(1, 2), 304, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[1], 308, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 0), 312, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 1), 316, endianess);
                this.setBufferFloat(this.bufferByte, matrixS.get(2, 2), 320, endianess);
                this.setBufferFloat(this.bufferByte, niftiOriginS[2], 324, endianess);
            }
            this.setBufferString(this.bufferByte, " ", 328);
        }
        if (oneFile) {
            this.setBufferString(this.bufferByte, "n+1\u0000", 344);
        }
        else {
            this.setBufferString(this.bufferByte, "ni1\u0000", 344);
        }
        this.bufferByte[348] = 0;
        this.bufferByte[349] = 0;
        this.bufferByte[350] = 0;
        this.bufferByte[351] = 0;
        if (!doGzip) {
            this.raFile.write(this.bufferByte);
            this.raFile.close();
        }
        return true;
    }
    
    private void writeHeader3DTo2D(final ModelImage image, String fileName, final String fileDir, final FileWriteOptions options, final boolean oneFile) throws IOException {
        final int beginSlice = options.getBeginSlice();
        final int endSlice = options.getEndSlice();
        final String origName = new String(fileName);
        for (int k = beginSlice, seq = options.getStartNumber(); k <= endSlice; ++k, ++seq) {
            fileName = origName;
            if (options.getDigitNumber() == 1) {
                fileName += Integer.toString(seq);
            }
            else if (options.getDigitNumber() == 2) {
                if (seq < 10) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            else if (options.getDigitNumber() == 3) {
                if (seq < 10) {
                    fileName = fileName + "00" + Integer.toString(seq);
                }
                else if (seq < 100) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            else if (options.getDigitNumber() == 4) {
                if (seq < 10) {
                    fileName = fileName + "000" + Integer.toString(seq);
                }
                else if (seq < 100) {
                    fileName = fileName + "00" + Integer.toString(seq);
                }
                else if (seq < 1000) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            this.writeHeader(image, 1, 1, fileName, fileDir, false, oneFile);
        }
    }
    
    private void writeHeader4DTo3D(final ModelImage image, String fileName, final String fileDir, final FileWriteOptions options, final boolean oneFile) throws IOException {
        final int beginTime = options.getBeginTime();
        final int endTime = options.getEndTime();
        final String origName = new String(fileName);
        for (int k = beginTime, seq = options.getStartNumber(); k <= endTime; ++k, ++seq) {
            fileName = origName;
            if (options.getDigitNumber() == 1) {
                fileName += Integer.toString(seq);
            }
            else if (options.getDigitNumber() == 2) {
                if (seq < 10) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            else if (options.getDigitNumber() == 3) {
                if (seq < 10) {
                    fileName = fileName + "00" + Integer.toString(seq);
                }
                else if (seq < 100) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            else if (options.getDigitNumber() == 4) {
                if (seq < 10) {
                    fileName = fileName + "000" + Integer.toString(seq);
                }
                else if (seq < 100) {
                    fileName = fileName + "00" + Integer.toString(seq);
                }
                else if (seq < 1000) {
                    fileName = fileName + "0" + Integer.toString(seq);
                }
                else {
                    fileName += Integer.toString(seq);
                }
            }
            this.writeHeader(image, image.getExtents()[2], 1, fileName, fileDir, false, oneFile);
        }
    }
    
    public byte[] getBufferByte() {
        return this.bufferByte;
    }
}
