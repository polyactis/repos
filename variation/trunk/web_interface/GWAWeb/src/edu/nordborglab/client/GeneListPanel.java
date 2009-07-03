/**
 * 2009-7-2
 * panel served in GeneListView
 */

package edu.nordborglab.client;

import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.visualization.client.AbstractDataTable;

public class GeneListPanel  extends CustomVerticalPanel {
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private HTML statusReport;
	private AbstractDataTable geneInfoDataTable;
	private WrapVisualizationTable geneInfoTable;
	
	private String geneListTypeID;	//not used
	private String fetchGeneListURL;
	
	GeneListPanel(AccessionConstants constants,
			DisplayJSONObject jsonErrorDialog, String geneListTypeID, String fetchGeneListURL) {
		super(constants, jsonErrorDialog, constants.GeneListPanelHelpID());
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		this.geneListTypeID = geneListTypeID;
		this.fetchGeneListURL = fetchGeneListURL;
		
		statusReport = new HTML(constants.LoadingText());
		statusReport.setVisible(false);
		add(statusReport);
		
		geneInfoTable = new WrapVisualizationTable();

		add(geneInfoTable);
		Common.fillInTable(fetchGeneListURL, geneInfoTable, jsonErrorDialog, statusReport, geneInfoDataTable);
		this.setSpacing(5);
	}
}
