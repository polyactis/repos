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

package affymetrix.gcos.psi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/** This class is used to read and access data stored in a PSI file. */
public class PSIFileData {

	/** The name of the PSI file. */
	private String fileName;

	/** The array of probe set names. */
	private List<ProbeSetInfo> probeSets = new ArrayList<ProbeSetInfo>();

	/**
	 * Gets the name of the PSI file.
	 * 
	 * @return The file name.
	 */
	public String getFileName() {
		return fileName;
	}

	/**
	 * Sets the name of the PSI file.
	 * 
	 * @param name
	 *          The name of the file.
	 */
	public void setFileName(String name) {
		fileName = name;
	}

	/**
	 * The number of probe sets in the file.
	 * 
	 * @return The number of probe sets.
	 */
	public int getProbeSetCount() {
		return probeSets.size();
	}

	/**
	 * Gets the probe set name.
	 * 
	 * @param index
	 *          The 0 based index to the probe set array.
	 * @return The probe set name.
	 */
	public String getProbeSetName(int index) {
		ProbeSetInfo ps = probeSets.get(index);
		return ps.getProbeSetName();
	}

	/**
	 * Gets the number of probe pairs in a set.
	 * 
	 * @param index
	 *          The 0 based index to the probe set array.
	 * @return The number of pairs.
	 */
	public int getProbePairs(int index) {
		ProbeSetInfo ps = probeSets.get(index);
		return ps.getNumberPairs();
	}

	/**
	 * Reads the contents of the PSI file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		BufferedReader b = null;
		clear();
		try {
			b = new BufferedReader(new FileReader(fileName));

			// Check that the first line is what we expect.
			String line;
			line = b.readLine();
			String header = "#Probe Sets: ";
			if (!line.startsWith(header)) {
				return false;
			}

			// Get the number of probe sets
			// int nsets = Integer.parseInt(line.substring(header.length()).trim());
			probeSets.clear();
			// probeSets.setSize(nsets);

			// Get the probe set names and #pairs
			while ((line = b.readLine()) != null) {
				String[] s = line.split("\t", 3);
				ProbeSetInfo info = new ProbeSetInfo();
				info.setProbeSetName(s[1]);
				info.setNumberPairs(Integer.parseInt(s[2].trim()));
				probeSets.add(info);
			}
		}
		catch (Throwable t) {
			clear();
			return false;
		}
		finally {
			try {
				if (b != null) {
					b.close();
				}
			}
			catch (Exception e) {
			}
		}
		return true;
	}

	/**
	 * Checks for the existance of a file.
	 * 
	 * @return True if exists.
	 */
	public boolean exists() {
		return new File(fileName).exists();
	}

	/** Clears the members. */
	public void clear() {
		probeSets.clear();
	}

	/** Creates a new instance of PSIFileData */
	public PSIFileData() {
		fileName = "";
		clear();
	}

}
