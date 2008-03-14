////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////

package affymetrix.gcos.cel;

import java.util.ArrayList;
import java.util.List;

import affymetrix.gcos.CoordinatePoint;
import affymetrix.gcos.GridCoordinates;
import affymetrix.gcos.TagValuePair;

/** Stores information in the header of a GCOS CEL file. */
public class CELFileHeaderData {

	/** Magic number for identifying XDA format */
	private int magic;

	public int getMagic() {
		return magic;
	}

	public void setMagic(int value) {
		magic = value;
	}

	/** CEL file format version number */
	private int version;

	public int getVersion() {
		return version;
	}

	public void setVersion(int value) {
		version = value;
	}

	/** Number of columns in array */
	private int nCols;

	public int getCols() {
		return nCols;
	}

	public void setCols(int value) {
		nCols = value;
	}

	/** Number of rows in array */
	private int nRows;

	public int getRows() {
		return nRows;
	}

	public void setRows(int value) {
		nRows = value;
	}

	/** Number of cells in array */
	private int nCells;

	public int getCells() {
		return nCells;
	}

	public void setCells(int value) {
		nCells = value;
	}

	/** Header information concatenated in a string */
	private String header;

	public String getHeader() {
		return header;
	}

	public void setHeader(String value) {
		header = value;
	}

	/** Algorithm name */
	private String alg;

	public String getAlg() {
		return alg;
	}

	public void setAlg(String value) {
		alg = value;
	}

	/** Chip type of array */
	private String chipType;

	public String getChipType() {
		return chipType;
	}

	public void setChipType(String value) {
		chipType = value;
	}

	/** DAT header string */
	private String datHeader;

	public String getDatHeader() {
		return datHeader;
	}

	public void setDatHeader(String value) {
		datHeader = value;
	}

	/** Cell margin */
	private int margin;

	public int getMargin() {
		return margin;
	}

	public void setMargin(int value) {
		margin = value;
	}

	/** Number of outliers */
	private int nOutliers;

	public int getOutliers() {
		return nOutliers;
	}

	public void setOutliers(int value) {
		nOutliers = value;
	}

	/** Number of masked cells */
	private int nMasked;

	public int getMasked() {
		return nMasked;
	}

	public void setMasked(int value) {
		nMasked = value;
	}

	/** Grid coordinates of array */
	private GridCoordinates cellGrid;

	public GridCoordinates getGrid() {
		return cellGrid;
	}

	public void setGrid(GridCoordinates value) {
		cellGrid = new GridCoordinates(value);
	}

	/** Algorithm parameters */
	private List<TagValuePair> parameters = new ArrayList<TagValuePair>();

	public List<TagValuePair> getParameters() {
		return parameters;
	}

	public void setParameters(List<TagValuePair> value) {
		parameters.clear();
		for (int i = 0; i < value.size(); i++) {
			TagValuePair param = new TagValuePair(value.get(i));
			parameters.add(param);
		}
	}

	/**
	 * Parses the parameters from the input string.
	 * 
	 * @param params
	 *          The parameters in the format of tag=value or tag:value
	 */
	public void setParameters(String params) {
		parameters.clear();
		String p1 = params.replace(';', ' ');
		p1 = p1.replace(':', '=');
		String s[] = p1.split(" ");
		for (int i = 0; i < s.length; i++) {
			String[] tagvalue = s[i].split("=");
			TagValuePair p = new TagValuePair();
			p.setTag(tagvalue[0]);
			p.setValue(tagvalue[1]);
			parameters.add(p);
		}
	}

	/** Parses the chip type from the header. */
	public void parseChipType() {
		int index = header.indexOf(".1sq");
		if (index == -1) {
			return;
		}
		int first = header.lastIndexOf(' ', index);
		chipType = header.substring(first + 1, index);
	}

	/** Parses the cell margin from the header. */
	public void parseMargin() {
		int index = header.indexOf("CellMargin");
		if (index == -1) {
			return;
		}
		int len = "CellMargin".length();
		margin = Integer.parseInt(header.substring(index + len + 1, index + len + 2).trim());
	}

	private CoordinatePoint parseCoordinate(String startString, String endString) {
		int start = header.indexOf(startString) + startString.length();
		int end = header.indexOf(endString, start);
		String[] scoord = header.substring(start, end).trim().split(" ");
		CoordinatePoint coord = new CoordinatePoint();
		coord.setX(Integer.parseInt(scoord[0].trim()));
		coord.setY(Integer.parseInt(scoord[1].trim()));
		return coord;
	}

	/** Parses the grid from the header. */
	public void parseGrid() {
		cellGrid = new GridCoordinates();
		cellGrid.setUpperLeft(parseCoordinate("GridCornerUL=", "GridCornerUR="));
		cellGrid.setUpperRight(parseCoordinate("GridCornerUR=", "GridCornerLR="));
		cellGrid.setLowerRight(parseCoordinate("GridCornerLR=", "GridCornerLL="));
		cellGrid.setLowerLeft(parseCoordinate("GridCornerLL=", "Axis"));
	}

	/** Parses the DatHeader part from the header. */
	public void parseDatHeader() {
		int index = header.indexOf("DatHeader=");
		if (index == -1) {
			return;
		}
		int last = header.indexOf("Algorithm=", index);
		String str = header.substring(index + 10, last - 1);
		datHeader = new String(str);

	}

	/** Creates a new instance of CELFileHeaderData */
	public CELFileHeaderData() {
		clear();
	}

	/** Clears the members. */
	public void clear() {
		magic = 0;
		version = 0;
		nCols = 0;
		nRows = 0;
		nCells = 0;
		header = "";
		alg = "";
		chipType = "";
		datHeader = "";
		margin = 0;
		nOutliers = 0;
		nMasked = 0;
		cellGrid = null;
		parameters.clear();
	}

}
