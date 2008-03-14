/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

package affymetrix.calvin.data;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.exception.DataGroupNotFoundException;
import affymetrix.calvin.exception.DataSetNotFoundException;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.DataGroupHeaderReader;
import affymetrix.calvin.parsers.DataSetHeaderReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides interfaces to store analysis results and data. */
public class GenericData {

	/** Creates a new instance of GenericData */
	public GenericData() {
		header = new FileHeader();
	}

	/**
	 * Returns a reference to the file identifier.
	 */
	public AffymetrixGuidType getFileIdentifier() {
		return header.getGenericDataHdr().getFileId();
	}

	/**
	 * Returns a reference to the parent array identifier.
	 */
	public AffymetrixGuidType getArrayIdentifier() {
		// Find the parent array file generic header
		GenericDataHeader hdr = header.getGenericDataHdr();
		int n = hdr.getParentCnt();
		for (int i = 0; i < n; i++) {
			GenericDataHeader dh = hdr.getParent(i);
			return dh.getFileId();
		}
		return null;
	}

	/**
	 * Returns a reference to the file header object
	 * 
	 * @return File header object
	 */
	public FileHeader getHeader() {
		return header;
	}

	/**
	 * Return the number of DataGroups in the GenericData object.
	 * 
	 * @return Number of DataGroups.
	 */
	public int getDataGroupCnt() {
		return header.getDataGroupCnt();
	}

	/**
	 * Return the names of the DataGroups in the generic data object.
	 * 
	 * @return vector that will receive the names of all DataGroups.
	 */
	public List<String> getDataGroupNames() {
		List<String> names = new ArrayList<String>();
		int n = header.getDataGroupCnt();
		for (int i = 0; i < n; i++) {
			DataGroupHeader dgh = header.getDataGroup(i);
			names.set(i, dgh.getName());
		}
		return names;
	}

	/**
	 * Return the number of DataSets in the DataGroup referenced by index.
	 * 
	 * @param dataGroupIdx
	 *          Index of the DataGroup.
	 * @return Number of DataSets associated with the DataGroup.
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 */
	public int getDataSetCnt(int dataGroupIdx) throws DataGroupNotFoundException {
		DataGroupHeader dch = findDataGroupHeader(dataGroupIdx);
		if (dch != null) {
			return dch.getDataSetCnt();
		}
		else {
			throw new DataGroupNotFoundException();
		}
	}

	/**
	 * Return the number of DataSets in the DataGroup referenced by name.
	 * 
	 * @param dataGroupName
	 *          Name of the DataGroup.
	 * @return Number of DataSets associated with the DataGroup.
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 */
	public int getDataSetCnt(String dataGroupName) throws DataGroupNotFoundException {
		DataGroupHeader dch = findDataGroupHeader(dataGroupName);
		if (dch != null) {
			return dch.getDataSetCnt();
		}
		else {
			throw new DataGroupNotFoundException();
		}
	}

	/**
	 * Return the DataSet names associated with a DataGroup.
	 * 
	 * @param dataGroupIdx
	 *          Index of the DataGroup from which to retrieve the DataSet names.
	 * @return Vector that will receive the names of all DataSets.
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 */
	public List<String> getDataSetNames(int dataGroupIdx) throws DataGroupNotFoundException {
		DataGroupHeader dch = findDataGroupHeader(dataGroupIdx);
		if (dch == null) {
			throw new DataGroupNotFoundException();
		}
		List<String> names = new ArrayList<String>();
		int n = dch.getDataSetCnt();
		for (int i = 0; i < n; i++) {
			DataSetHeader dsh = dch.getDataSet(i);
			names.set(i, dsh.getName());
		}
		return names;
	}

	/**
	 * Return the DataSet names associated with a DataGroup.
	 * 
	 * @param dataGroupName
	 *          Name of the DataGroup from which to retrieve the DataSet names.
	 * @return vector that will receive the names of all DataSets.
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 */
	public List<String> getDataSetNames(String dataGroupName) throws DataGroupNotFoundException {
		DataGroupHeader dch = findDataGroupHeader(dataGroupName);
		if (dch == null) {
			throw new DataGroupNotFoundException();
		}
		List<String> names = new ArrayList<String>();
		int n = dch.getDataSetCnt();
		for (int i = 0; i < n; i++) {
			DataSetHeader dsh = dch.getDataSet(i);
			names.set(i, dsh.getName());
		}
		return names;
	}

	/**
	 * Returns a pointer to the DataSet object by DataGroup and DataSet index. Each call will return a new DataSet object.
	 * The caller should call Delete when finished with the DataSet.
	 * 
	 * @param dataGroupIdx
	 *          The index of the DataGroup from which to find the DataSet.
	 * @param dataSetIdx
	 *          The index of the DataSet to return.
	 * @return DataSet
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 * @exception affymetrix_calvin_exceptions::DataSetNotFoundException
	 *              DataSet not found.
	 */
	public DataSet getDataSet(int dataGroupIdx, int dataSetIdx) throws DataSetNotFoundException,
			DataGroupNotFoundException, FileNotFoundException {
		open();

		DataGroupHeader dch = findDataGroupHeader(dataGroupIdx);
		if (dch != null) {
			DataSetHeader dph = findDataSetHeader(dch, dataSetIdx);
			if (dph != null) {
				return createDataSet(dph);
			}
			else {
				throw new DataSetNotFoundException();
			}
		}
		else {
			throw new DataGroupNotFoundException();
		}
	}

	/**
	 * Returns a pointer to the DataSet object by DataGroup and DataSet name. Each call will return a new DataSet object.
	 * The caller should call Delete when finished with the DataSet.
	 * 
	 * @param dataGroupName
	 *          The name of the DataGroup from which to find the DataSet.
	 * @param dataSetName
	 *          The name of the DataSet to return.
	 * @return DataSet
	 * @exception affymetrix_calvin_exceptions::DataGroupNotFoundException
	 *              DataGroup not found.
	 * @exception affymetrix_calvin_exceptions::DataSetNotFoundException
	 *              DataSet not found.
	 */
	public DataSet getDataSet(String dataGroupName, String dataSetName) throws FileNotFoundException,
			DataSetNotFoundException, DataGroupNotFoundException {
		open();

		DataGroupHeader dch = findDataGroupHeader(dataGroupName);
		if (dch != null) {
			DataSetHeader dph = findDataSetHeader(dch, dataSetName);
			if (dph != null) {
				return createDataSet(dph);
			}
			else {
				throw new DataSetNotFoundException();
			}
		}
		else {
			throw new DataGroupNotFoundException();
		}
	}

	/**
	 * Creates a new DataSet
	 * 
	 * @param dsh
	 *          The DataSetHeader of the DataSet to create.
	 * @return The new DataSet
	 */
	public DataSet createDataSet(DataSetHeader dsh) {
		readFullDataSetHeader(dsh);
		return new DataSet(header.getFilename(), dsh, fileStream);
	}

	/**
	 * Returns a DataGroup object based on a DataGroup file position. This is useful when there are many DataGroups and
	 * the file position of each DataGroup is known (Calvin CDF). In this case the GenericFileReader::ReadHeader() method
	 * should be called with the ReadNoDataGroupHeader flag.
	 * 
	 * @param dataGroupFilePos
	 *          File position of the DataGroup in the current file
	 * @return DataGroup object.
	 */
	public DataGroup getDataGroup(int dataGroupFilePos) throws IOException, FileNotFoundException,
			UnsignedOutOfLimitsException {

		// Open a file stream
		open();

		// Read the DataGroupHeader and all DataSetHeaders
		FileChannel fc = fileStream.getChannel();
		fc.position(dataGroupFilePos);
		DataGroupHeader dch = new DataGroupHeader();
		DataGroupHeaderReader reader = new DataGroupHeaderReader();
		reader.read(fileStream, dch);
		return new DataGroup(header.getFilename(), dch, fileStream);
	}

	/**
	 * Read the full DataSetHeader if it has only been parially read.
	 * 
	 * @param dph
	 *          Pointer to the DataSetHeader to read
	 */
	public void readFullDataSetHeader(DataSetHeader dph) {

		// Check if the DataSet has been read fully.
		if (isDSHPartiallyRead(dph)) {
			// Open a file stream
			try {
				open();

				// Read the DataGroupHeader and all DataSetHeaders
				FileChannel fc = fileStream.getChannel();
				fc.position(dph.getHeaderStartFilePos().toLong());

				// Read the header
				DataSetHeaderReader reader = new DataSetHeaderReader();
				reader.read(fileStream, dph);
			}
			catch (Throwable t) {
			}
		}
	}

	/**
	 * Determine if the DataSetHeader has been partially read.
	 * 
	 * @param dph
	 *          Pointer to the DataSetHeader to check
	 * @return true if the dph has only been partially read or is 0, otherwise false.
	 */
	public boolean isDSHPartiallyRead(DataSetHeader dph) {
		if (dph == null) {
			return false;
		}
		if ((dph.getRowCnt() == 0) && (dph.getColumnCnt() == 0)) {
			return true;
		}
		return false;
	}

	/**
	 * Finds a DataSetHeader by name.
	 * 
	 * @param name
	 *          The name of the DataGroup
	 * @return A pointer to the DataGroupHeader. If not found, the return is 0.
	 */
	public DataGroupHeader findDataGroupHeader(String name) {
		return header.findDataGroupHeader(name);
	}

	/**
	 * Finds a DataGroupHeader by index.
	 * 
	 * @param index
	 *          The index of the DataGroup.
	 * @return A pointer to the DataGroupHeader. If not found, the return is 0.
	 */
	public DataGroupHeader findDataGroupHeader(int index) {
		DataGroupHeader dch = null;
		if ((index >= 0) && (index < header.getDataGroupCnt())) {
			dch = header.getDataGroup(index);
		}
		return dch;
	}

	/**
	 * Finds a DataSetHeader by index.
	 * 
	 * @param dch
	 *          The DataGroupHeader of the DataGroup to which the DataSet belongs.
	 * @param dataSetIdx
	 *          The DataSet index of the DataSetHeader to find.
	 * @return A pointer to the DataSetHeader if it is found, otherwise 0.
	 */
	public static DataSetHeader findDataSetHeader(DataGroupHeader dch, int dataSetIdx) {
		DataSetHeader dph = null;
		if (dch != null) {
			if ((dataSetIdx >= 0) && (dataSetIdx < dch.getDataSetCnt())) {
				dph = dch.getDataSet(dataSetIdx);
			}
		}
		return dph;
	}

	/**
	 * Finds a DataSetHeader by name.
	 * 
	 * @param dch
	 *          The DataGroupHeader of the DataGroup to which the DataSet belongs.
	 * @param dataSetName
	 *          The DataSet name of the DataSetHeader to find.
	 * @return A pointer to the DataSetHeader if it is found, otherwise 0.
	 */
	public static DataSetHeader findDataSetHeader(DataGroupHeader dch, String dataSetName) {
		DataSetHeader dph = null;
		if (dch != null) {
			dph = dch.findDataSetHeader(dataSetName);
		}
		return dph;
	}

	/**
	 * Opens the file for access.
	 * 
	 * @return True if opened successfully.
	 */
	public void open() throws FileNotFoundException {
		// Open a file stream
		if (fileStream == null) {
			fileStream = new FileInputStream(header.getFilename());
		}
	}

	/**
	 * Closes the file.
	 */
	public void close() {
		fileStream = null;
	}

	/** The header and generic header objects */
	protected FileHeader header;

	/** The map of the file. */
	private FileInputStream fileStream;
}
