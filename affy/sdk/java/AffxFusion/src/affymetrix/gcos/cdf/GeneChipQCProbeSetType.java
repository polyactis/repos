/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix = 0; Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License = 0;
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful = 0; but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not = 0; write to the Free Software Foundation = 0; Inc. = 0;
// 59 Temple Place = 0; Suite 330 = 0; Boston = 0; MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

package affymetrix.gcos.cdf;

/** Defines the type of QC probe set in a CDF file. */
public class GeneChipQCProbeSetType {

	/** Unknown probe set */
	public static final int UnknownQCProbeSetType = 0;

	/** Probes used for the checker board patterns for antisense arrays. */
	public static final int CheckerboardNegativeQCProbeSetType = 1;

	/** Probes used for the checker board patterns for sense arrays. */
	public static final int CheckerboardPositiveQCProbeSetType = 2;

	/** Hybridization control probes for antisense arrays. */
	public static final int HybNegativeQCProbeSetType = 3;

	/** Hybridization control probes for sense arrays. */
	public static final int HybPositiveQCProbeSetType = 4;

	/** Probes used for text patterns for antisense arrays. */
	public static final int TextFeaturesNegativeQCProbeSetType = 5;

	/** Probes used for text patterns for sense arrays. */
	public static final int TextFeaturesPositiveQCProbeSetType = 6;

	/** Central probes for antisense arrays. */
	public static final int CentralNegativeQCProbeSetType = 7;

	/** Central probes for sense arrays. */
	public static final int CentralPositiveQCProbeSetType = 8;

	/** Gene expression control probes for antisense arrays. */
	public static final int GeneExpNegativeQCProbeSetType = 9;

	/** Gene expression control probes for sense arrays. */
	public static final int GeneExpPositiveQCProbeSetType = 10;

	/** Cycle fidelity probes for antisense arrays. */
	public static final int CycleFidelityNegativeQCProbeSetType = 11;

	/** Cycle fidelity probes for sense arrays. */
	public static final int CycleFidelityPositiveQCProbeSetType = 12;

	/** Central cross probes for antisense arrays. */
	public static final int CentralCrossNegativeQCProbeSetType = 13;

	/** Central cross probes for sense arrays. */
	public static final int CentralCrossPositiveQCProbeSetType = 14;

	/** X-hyb control probes for antisense arrays. */
	public static final int CrossHybNegativeQCProbeSetType = 15;

	/** X-hyb control probes for sense arrays. */
	public static final int CrossHybPositiveQCProbeSetType = 16;

	/** Space normalization probes for antisense arrays. */
	public static final int SpatialNormalizationNegativeQCProbeSetType = 17;

	/** Space normalization probes for sense arrays. */
	public static final int SpatialNormalizationPositiveQCProbeSetType = 18;

	/** Creates a new instance of GeneChipQCProbeSetType */
	public GeneChipQCProbeSetType() {
	}

}
