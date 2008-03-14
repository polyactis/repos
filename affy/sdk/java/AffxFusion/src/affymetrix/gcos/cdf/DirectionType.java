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

/** The direction of the target the probes are designed to interrogate. */
public class DirectionType {

	/** No direction specified */
	public static final int NoDirection = 0;

	/** Sense */
	public static final int SenseDirection = 1;

	/** Anti sense */
	public static final int AntiSenseDirection = 2;

	/** Either */
	public static final int EitherDirection = 3;

	/** Creates a new instance of DirectionType */
	public DirectionType() {
	}

}
