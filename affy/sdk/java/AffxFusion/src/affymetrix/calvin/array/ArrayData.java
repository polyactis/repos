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

import affymetrix.calvin.array.CreateStep.CreateStepType;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.calvin.utils.IOUtils;

/** This class provides interfaces to store array information. */
public class ArrayData {

	/** A unique idendifier for the array set object */
	private AffymetrixGuidType fileId = null;

	/** An identifier to the type of data stored in the file */
	private AffymetrixGuidType dataTypeId = null;

	/** The step in Calvin that created the array set data. */
	private CreateStepType createdStep = CreateStepType.NoStep;

	/** The name of the project that initially created the array set data. */
	private String initialProject = IOUtils.EMPTY;

	/** The date and time of initial creation. */
	private String creationDateTime = IOUtils.EMPTY;

	/** The user who created the data object. */
	private String createdBy = IOUtils.EMPTY;

	/** The arrays attributes for the arrays in the set */
	private List<ArrayAttributes> physicalArraysAttributes = new ArrayList<ArrayAttributes>();

	/** The user attributes */
	private List<ParameterNameValueDefaultRequired> userAttributes = new ArrayList<ParameterNameValueDefaultRequired>();

	/** Creates a new instance of ArrayData. */
	public ArrayData() {
	}

	/**
	 * The unique idendifier for the array set.
	 * 
	 * @return The unique idendifier for the array set.
	 */
	public AffymetrixGuidType getArraySetFileIdentifier() {
		return fileId;
	}

	/**
	 * The unique idendifier for the array set.
	 * 
	 * @param value
	 *          The unique idendifier for the array set.
	 */
	public void setArraySetFileIdentifier(AffymetrixGuidType value) {
		fileId = value;
	}

	/**
	 * The identifier of the type of data stored in the file.
	 * 
	 * @return The identifier of the type of data.
	 */
	public AffymetrixGuidType getDataTypeIdentifier() {
		return dataTypeId;
	}

	/**
	 * The identifier of the type of data stored in the file.
	 * 
	 * @param value
	 *          The identifier of the type of data.
	 */
	public void setDataTypeIdentifier(AffymetrixGuidType value) {
		dataTypeId = value;
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
	 * The name of the project that initially created the array set data.
	 * 
	 * @return The project name.
	 */
	public String getInitialProject() {
		return initialProject;
	}

	/**
	 * The name of the project that initially created the array set data.
	 * 
	 * @param value
	 *          The project name.
	 */
	public void setInitialProject(String value) {
		initialProject = value;
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
	 * @return The user name.
	 */
	public String getCreatedBy() {
		return createdBy;
	}

	/**
	 * The user who created the data object.
	 * 
	 * @param value
	 *          The user name.
	 */
	public void setCreatedBy(String value) {
		createdBy = value;
	}

	/**
	 * The arrays attributes. Each array in a set will have its own attributes.
	 * 
	 * @return The vector of arrays attributes.
	 */
	public List<ArrayAttributes> getPhysicalArraysAttributes() {
		return physicalArraysAttributes;
	}

	/**
	 * The arrays attributes. Each array in a set will have its own attributes.
	 * 
	 * @param value
	 *          The vector of arrays attributes.
	 */
	public void addPhysicalArraysAttributes(List<ArrayAttributes> value) {
		physicalArraysAttributes.addAll(value);
	}

	/**
	 * The user attributes.
	 * 
	 * @return The vector of user attributes.
	 */
	public List<ParameterNameValueDefaultRequired> getUserAttributes() {
		return userAttributes;
	}

	/**
	 * The user attributes.
	 * 
	 * @param value
	 *          The vector of user attributes.
	 */
	public void addUserAttributes(List<ParameterNameValueDefaultRequired> value) {
		userAttributes.addAll(value);
	}

	/** Clears the member objects. */
	public void clear() {
		physicalArraysAttributes.clear();
		userAttributes.clear();
		fileId = null;
		dataTypeId = null;
		createdStep = CreateStepType.NoStep;
		initialProject = "";
		creationDateTime = "";
		createdBy = "";
	}
}
