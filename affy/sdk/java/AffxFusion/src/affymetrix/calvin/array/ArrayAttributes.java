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

package affymetrix.calvin.array;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.array.ArrayMedia.ArrayMediaType;
import affymetrix.calvin.array.CreateStep.CreateStepType;
import affymetrix.calvin.array.PATAssignmentMethod.PATAssignmentMethodType;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides interfaces to store physical array attributes. */
public class ArrayAttributes {

	/** A unique idendifier for the array object */
	private AffymetrixGuidType id;

	/** The array attributes */
	private List<ParameterNameValue> attributes = new ArrayList<ParameterNameValue>();

	/** The array name. */
	private String arrayName;

	/** The barcode on the array cartridge. */
	private String arrayBarcode;

	/** The type of assembly. */
	private ArrayMediaType media;

	/** The row number of the plate or strip. */
	private int mediaRow;

	/** The column number of the plate or strip. */
	private int mediaCol;
	
	/** The name of the media file. */
	private String mediaFileName;
	
	/** The media file guid. */
	private AffymetrixGuidType mediaFileGUID;

	/** A customer barcode. */
	private String customerBarcode;

	/** The associated master file. */
	private String masterFile;

	/** A unique idendifier for the master file */
	private AffymetrixGuidType masterFileId;

	/** The method the probe array type was assigned. */
	private PATAssignmentMethodType patAssignment;

	/** The date the array object was created. */
	private String creationDateTime;

	/** The user who created the data object. */
	private String createdBy;

	/** A user comment. */
	private String comment = null;

	/** The step in Calvin that created the array set data. */
	private CreateStepType createdStep = CreateStepType.NoStep;
	
	/*! The name of the library file package. */
	private String libraryPackageName = null;

	/**
	 * The unique idendifier for the object.
	 * 
	 * @return The unique idendifier for the object.
	 */
	public AffymetrixGuidType getIdentifier() {
		return id;
	}

	/**
	 * The unique idendifier for the object.
	 * 
	 * @param value
	 *          The unique idendifier for the object.
	 */
	public void setIdentifier(AffymetrixGuidType value) {
		id = value;
	}

	/**
	 * The array name.
	 * 
	 * @return The array name.
	 */
	public String getArrayName() {
		return arrayName;
	}

	/**
	 * The array name.
	 * 
	 * @param value
	 *          The array name.
	 */
	public void setArrayName(String value) {
		arrayName = value;
	}

	/**
	 * The barcode on the array cartridge.
	 * 
	 * @return The barcode.
	 */
	public String getArrayBarcode() {
		return arrayBarcode;
	}

	/**
	 * The barcode on the array cartridge.
	 * 
	 * @param value
	 *          The barcode.
	 */
	public void setArrayBarcode(String value) {
		arrayBarcode = value;
	}

	/**
	 * The type of assembly.
	 * 
	 * @return The assembly type.
	 */
	public ArrayMediaType getMedia() {
		return media;
	}

	/**
	 * The type of assembly.
	 * 
	 * @param value
	 *          The assembly type.
	 */
	public void setMedia(ArrayMediaType value) {
		media = value;
	}

	/**
	 * The row number of the media or strip.
	 * 
	 * @return The row.
	 */
	public int getMediaRow() {
		return mediaRow;
	}

	/**
	 * The row number of the media or strip.
	 * 
	 * @param value
	 *          The row.
	 */
	public void setMediaRow(int value) {
		mediaRow = value;
	}

	/**
	 * The column number of the media or strip.
	 * 
	 * @return The column.
	 */
	public int getMediaCol() {
		return mediaCol;
	}

	/**
	 * The column number of the media or strip.
	 * 
	 * @param value
	 *          The column.
	 */
	public void setMediaCol(int value) {
		mediaCol = value;
	}

	/**
	 * A customer barcode.
	 * 
	 * @return The barcode.
	 */
	public String getCustomerBarcode() {
		return customerBarcode;
	}

	/**
	 * A customer barcode.
	 * 
	 * @param value
	 *          The barcode.
	 */
	public void setCustomerBarcode(String value) {
		customerBarcode = value;
	}

	/**
	 * The associated master file.
	 * 
	 * @return The master file name.
	 */
	public String getMasterFile() {
		return masterFile;
	}

	/**
	 * The associated master file.
	 * 
	 * @param value
	 *          The master file name.
	 */
	public void setMasterFile(String value) {
		masterFile = value;
	}

	/**
	 * The unique idendifier for the master file.
	 * 
	 * @return The master file guid.
	 */
	public AffymetrixGuidType getMasterFileId() {
		return masterFileId;
	}

	/**
	 * The unique idendifier for the master file.
	 * 
	 * @param value
	 *          The master file guid.
	 */
	public void setMasterFileId(AffymetrixGuidType value) {
		masterFileId = value;
	}

	/**
	 * The method the probe array type was assigned.
	 * 
	 * @return The assignment method.
	 */
	public PATAssignmentMethodType getPATAssignment() {
		return patAssignment;
	}

	/**
	 * The method the probe array type was assigned.
	 * 
	 * @param value
	 *          The assignment method.
	 */
	public void setPATAssignment(PATAssignmentMethodType value) {
		patAssignment = value;
	}

	/**
	 * The date and time of initial creation.
	 * 
	 * @return The creation date and time.
	 */
	public String getCreationDateTime() {
		return creationDateTime;
	}

	/**
	 * The date and time of initial creation.
	 * 
	 * @param value
	 *          The creation date and time.
	 */
	public void setCreationDateTime(String value) {
		creationDateTime = value;
	}

	/**
	 * The user who created the data object.
	 * 
	 * @return The user who created the data object.
	 */
	public String getCreatedBy() {
		return createdBy;
	}

	/**
	 * The user who created the data object.
	 * 
	 * @param value
	 *          The user who created the data object.
	 */
	public void setCreatedBy(String value) {
		createdBy = value;
	}

	/**
	 * A user comment.
	 * 
	 * @return A user comment.
	 */
	public String getComment() {
		return comment;
	}

	/**
	 * A user comment.
	 * 
	 * @param value
	 *          A user comment.
	 */
	public void setComment(String value) {
		comment = value;
	}

	/**
	 * The step in Calvin that created the array set data.
	 * 
	 * @return The step in calvin that create the array set data.
	 */
	public CreateStepType getCreatedStep() {
		return createdStep;
	}

	/**
	 * The step in Calvin that created the array set data.
	 * 
	 * @param value
	 *          The step in calvin that create the array set data.
	 */
	public void setCreatedStep(CreateStepType value) {
		createdStep = value;
	}

	/**
	 * The array attributes.
	 * 
	 * @return The vector of array attributes.
	 */
	public List<ParameterNameValue> getAttributes() {
		return attributes;
	}

	/** Clears the member objects. */
	public void clear() {
		attributes.clear();
		id = null;
		arrayName = null;
		arrayBarcode = null;
		media = ArrayMediaType.CartridgeMedia;
		mediaRow = 0;
		mediaCol = 0;
		customerBarcode = null;
		masterFile = null;
		masterFileId = null;
		patAssignment = PATAssignmentMethodType.NoAssignment;
		creationDateTime = null;
		createdBy = null;
		createdStep = CreateStepType.NoStep;
		comment = null;
	}

	/** Constructs a new ArrayAttributes object. */
	public ArrayAttributes() {
		clear();
	}

	public String getMediaFileName() {
		return mediaFileName;
	}

	public void setMediaFileName(String mediaFileName) {
		this.mediaFileName = mediaFileName;
	}

	public AffymetrixGuidType getMediaFileGUID() {
		return mediaFileGUID;
	}

	public void setMediaFileGUID(AffymetrixGuidType mediaFileGUID) {
		this.mediaFileGUID = mediaFileGUID;
	}

	public String getLibraryPackageName() {
		return libraryPackageName;
	}

	public void setLibraryPackageName(String libraryPackageName) {
		this.libraryPackageName = libraryPackageName;
	}
}
