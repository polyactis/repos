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

import affymetrix.calvin.exception.DataSetNotFoundException;

/** This class provides methods to get a DataSet in a DataGroup. */
public class DataGroup {

	/**
	 * Constructor
	 * 
	 * @param filename_
	 *          The name of the generic file to access.
	 * @param dch
	 *          The DataGroupHeader of the DataGroup to access.
	 * @param handle_
	 *          A handle to the file mapping object
	 */
	public DataGroup(String filename_, DataGroupHeader dch, FileInputStream handle_) {
		filename = filename_;
		dataGroupHeader = dch;
		handle = handle_;
	}

	/**
	 * Method to get a reference to the DataGroupHeader
	 * 
	 * @return A reference to the DataGroupHeader.
	 */
	public DataGroupHeader getHeader() {
		return dataGroupHeader;
	}

	/**
	 * Returns a pointer to the DataSet object by DataSet index. Each call will return a new DataSet object. The caller
	 * should call Delete when finished with the DataSet.
	 * 
	 * @param dataSetIdx
	 *          The index of the DataSet to return.
	 * @return DataSet
	 * @exception affymetrix_calvin_exceptions::DataSetNotFoundException
	 *              DataSet not found.
	 */
	public DataSet getDataSet(int dataSetIdx) throws DataSetNotFoundException {
		DataSetHeader dph = GenericData.findDataSetHeader(dataGroupHeader, dataSetIdx);
		if (dph != null) {
			return new DataSet(filename, dph, handle);
		}
		else {
			throw new DataSetNotFoundException();
		}
	}

	/**
	 * Returns a pointer to the DataSet object by DataSet name. Each call will return a new DataSet object. The caller
	 * should call Delete when finished with the DataSet.
	 * 
	 * @param dataSetName
	 *          The name of the DataSet to return.
	 * @return DataSet
	 * @exception affymetrix_calvin_exceptions::DataSetNotFoundException
	 *              DataSet not found.
	 */
	public DataSet getDataSet(String dataSetName) throws DataSetNotFoundException {
		DataSetHeader dph = GenericData.findDataSetHeader(dataGroupHeader, dataSetName);
		if (dph != null) {
			return new DataSet(filename, dph, handle);
		}
		else {
			throw new DataSetNotFoundException();
		}
	}

	/** Name of the generic file to access */
	protected String filename;

	/** DataGroupHeader of the DataGroup from which to get DataSets */
	protected DataGroupHeader dataGroupHeader;

	/** File input handle */
	protected FileInputStream handle;

}
