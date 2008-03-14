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

/** Defines the steps that can create array data. */
public class CreateStep {

	/** The none step. */
	public static final String CREATED_STEP_NONE = "None";

	/** The array registration step. */
	public static final String CREATED_STEP_ARRAY_REG = "ArrayRegistration";

	/** The scan step. */
	public static final String CREATED_STEP_SCAN = "Scanning";

	/** The grid alignment step. */
	public static final String CREATED_STEP_GRID = "Gridding";

	/** The cel analysis step. */
	public static final String CREATED_STEP_CEL = "CELAnalysis";

	/** Other undefined step. */
	public static final String CREATED_STEP_OTHER = "Other";

	/** From step. */
	public static final String CREATED_STEP_FROM = "From";

	/** Job order step. */
	public static final String CREATED_STEP_JOB_ORDER = "JobOrderServer";

	/** File indexer step. */
	public static final String CREATED_STEP_FILE_INDEXER = "FileIndexer";

	/** Defines the steps that can create array data. */
	public enum CreateStepType {
		NoStep, // No step.
		ArrayRegistrationStep, // Array registration.
		ScanningStep, // Scanning.
		GriddingStep, // Grid analysis.
		CELAnalysisStep, // CEL file analysis.
		FromStep, // From.
		JobOrderServerStep, // Job order.
		FileIndexerStep, // File indexer.
		OtherStep
	}

	/** Ctor */
	private CreateStep() {
	}

	
	/**
	 * Converts the step type to a string.
	 * 
	 * @return The string representation.
	 */
	public static String toString(CreateStepType t) {
		switch (t) {
		case NoStep:
			return CREATED_STEP_NONE;

		case ArrayRegistrationStep:
			return CREATED_STEP_ARRAY_REG;

		case ScanningStep:
			return CREATED_STEP_SCAN;

		case GriddingStep:
			return CREATED_STEP_GRID;

		case CELAnalysisStep:
			return CREATED_STEP_CEL;

		case FromStep:
			return CREATED_STEP_FROM;

		case JobOrderServerStep:
			return CREATED_STEP_JOB_ORDER;

		case FileIndexerStep:
			return CREATED_STEP_FILE_INDEXER;

		case OtherStep:
			return CREATED_STEP_OTHER;
		}
		return null;
	}

	/**
	 * Converts a string to step type.
	 * 
	 * @param s
	 *          The string representation.
	 */
	public static CreateStepType toCreateStepType(String s) {
		if (s.equals(CREATED_STEP_NONE)) {
			return CreateStepType.NoStep;
		}
		else if (s.equals(CREATED_STEP_ARRAY_REG)) {
			return CreateStepType.ArrayRegistrationStep;
		}
		else if (s.equals(CREATED_STEP_SCAN)) {
			return CreateStepType.ScanningStep;
		}
		else if (s.equals(CREATED_STEP_GRID)) {
			return CreateStepType.GriddingStep;
		}
		else if (s.equals(CREATED_STEP_CEL)) {
			return CreateStepType.CELAnalysisStep;
		}
		else if (s.equals(CREATED_STEP_FROM)) {
			return CreateStepType.FromStep;
		}
		else if (s.equals(CREATED_STEP_JOB_ORDER)) {
			return CreateStepType.JobOrderServerStep;
		}
		else if (s.equals(CREATED_STEP_FILE_INDEXER)) {
			return CreateStepType.FileIndexerStep;
		}
		return CreateStepType.OtherStep;
	}
}
