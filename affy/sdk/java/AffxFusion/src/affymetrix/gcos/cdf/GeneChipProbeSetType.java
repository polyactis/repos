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

package affymetrix.gcos.cdf;

/** The types of probe sets in a CHP file. */
public class GeneChipProbeSetType {

	/** Unknown probe set */
	public static final int UnknownProbeSetType = 0;

	/** Expression probe set */
	public static final int ExpressionProbeSetType = 1;

	/** Genotyping probe set */
	public static final int GenotypingProbeSetType = 2;

	/** Resequencing probe set */
	public static final int ResequencingProbeSetType = 3;

	/** Tag probe set */
	public static final int TagProbeSetType = 4;

	/** Tag probe set */
	public static final int CopyNumberProbeSetType = 5;

	/** Tag probe set */
	public static final int GenotypeControlProbeSetType = 6;

	/** Tag probe set */
	public static final int ExpressionControlProbeSetType = 7;

	/** Creates a new instance of GeneChipProbeSetType */
	public GeneChipProbeSetType() {
	}

}
