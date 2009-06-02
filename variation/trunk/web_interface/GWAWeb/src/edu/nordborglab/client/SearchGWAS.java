package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class SearchGWAS implements EntryPoint {
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	private SearchGWASByGeneName searchGWASByGeneName;
	private TabPanel tPanel;
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		tPanel = new TabPanel();
		
		//tPanel.selectTab(1);
		//tPanel.setWidth("1000px");
		tPanel.setAnimationEnabled(true);
		searchGWASByGeneName = new SearchGWASByGeneName(constants, jsonErrorDialog);
		tPanel.add(searchGWASByGeneName, "By Gene Name");
		tPanel.selectTab(0);
		
		// Add it to the root panel.
		RootPanel.get("gwt").add(tPanel);
	}
	
}
