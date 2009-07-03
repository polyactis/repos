/**
 * 2009-7-2
 * 	class to view information of genes coming from a specified gene list
 */
package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.ui.RootPanel;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class GeneListView implements EntryPoint {

	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
		
	
	private String geneListTypeID;
	private String fetchGeneListURL;
	private GeneListPanel geneListPanel;
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		geneListTypeID = Common.getGeneListTypeID();
		fetchGeneListURL = Common.getFetchGeneListURL();		
		geneListPanel = new GeneListPanel(constants, jsonErrorDialog, geneListTypeID, fetchGeneListURL);
		RootPanel.get("geneListView").add(geneListPanel);

	}
}
