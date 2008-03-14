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

package affymetrix.calvin.parameter;

/** Defines constant names for static attributes. */
public class AffymetrixParameterConsts {

	/** Defines the static attribute name for the probe array type of the physical array. */
	public static final String ARRAY_TYPE_PARAM_NAME = "affymetrix-array-type";

	/** Defines the static attribute name for the master file. */
	public static final String MASTER_FILE_PARAM_NAME = "affymetrix-master-file";

	/** Defines the static attribute name for the library package. */
	public static final String LIBRARY_PACKAGE_PARAM_NAME = "affymetrix-library-package";

	// Defines the number of characters to reserve in the parameter list for the array type name
	public static final int ARRAY_TYPE_MAX_LEN = 100;

	/** Defines the static attribute name for the barcode of the physical array. */
	public static final String ARRAY_BARCODE_PARAM_NAME = "affymetrix-array-barcode";

	/** Defines the static attribute name for the array lot number of the physical array. */
	public static final String ARRAY_LOT_PARAM_NAME = "affymetrix-array-lot";

	/** Defines the static attribute name for the expiration date of the physical array. */
	public static final String ARRAY_EXPIRE_DATE_PARAM_NAME = "affymetrix-array-expiration-date";

	/** Defines the static attribute name for the barcode of the array plate. */
	public static final String PLATE_BARCODE_PARAM_NAME = "affymetrix-plate-barcode";

	/** Defines the static attribute name for the plate type of the array plate. */
	public static final String PLATE_TYPE_PARAM_NAME = "affymetrix-plate-type";

	/** Defines the static attribute name for the row location of the array plate well. */
	public static final String PLATE_WELL_ROW_PARAM_NAME = "affymetrix-plate-well-row";

	/** Defines the static attribute name for the column location of the array plate well. */
	public static final String PLATE_WELL_COL_PARAM_NAME = "affymetrix-plate-well-col";

	/** Defines the static attribute name for the algorithm name. */
	public static final String ALGORITHM_NAME_PARAM_NAME = "affymetrix-algorithm-name";

	/** Defines the static attribute for the algorithm version. */
	public static final String ALG_VERSION_PARAM_NAME = "affymetrix-algorithm-version";

	/** Defines the static attribute prefix for algorithm parameter names */
	public static final String ALGORITHM_PARAM_NAME_PREFIX = "affymetrix-algorithm-param-";

	/** A prefix for chip summary parameter ids. */
	public static final String CHIP_SUMMARY_PARAM_NAME_PREFIX = "affymetrix-chipsummary-";

	/** Defines the static attribute name for the DATHeader */
	public static final String DAT_HEADER_PARAM_NAME = "affymetrix-dat-header";

	/** Defines the static attribute name for the partial DATHeader */
	public static final String PARTIAL_DAT_HEADER_PARAM_NAME = "affymetrix-partial-dat-header";

	/** Defines the static attribute name for the max pixel intensity */
	public static final String MAX_PIXEL_INTENSITY_PARAM_NAME = "affymetrix-max-pixel-intensity";

	/** Defines the static attribute name for the min pixel intensity */
	public static final String MIN_PIXEL_INTENSITY_PARAM_NAME = "affymetrix-min-pixel-intensity";

	/** Defines the static attribute name for the orientation */
	public static final String ORIENTATION_PARAM_NAME = "affymetrix-image-orientation";

	/** Defines the static attribute name for the file version. This is not the file format version. */
	public static final String FILE_VERSION_PARAM_NAME = "affymetrix-file-version";

	/** Defines the static attribute name for the flip-flag which indicates if an image is flipped about the y-axis. */
	public static final String FLIP_FLAG_PARAM_NAME = "affymetrix-image-flip-flag";

	/** CDF Data Type Expression */
	public static final String AFFY_EXPR_PS = "affymetrix-expression-probesets";

	/** CDF Data Type Genotyping */
	public static final String AFFY_GENO_PS = "affymetrix-genotyping-probesets";

	/** CDF Data Type Tag */
	public static final String AFFY_TAG_PS = "affymetrix-tag-probesets";

	/** CDF Data Type Resequencing */
	public static final String AFFY_RESEQ_PS = "affymetrix-resequencing-probesets";

	/** CDF Data Type Control */
	public static final String AFFY_CNTRL_PS = "affymetrix-control-probesets";

	/** Defines US English locale. */
	public static final String US_ENGLISH_LOCALE = "en-US";

	/** Defines an identifier for the scan acquisition data file. */
	public static final String SCAN_ACQUISITION_DATA_TYPE = "affymetrix-calvin-scan-acquisition";

}
