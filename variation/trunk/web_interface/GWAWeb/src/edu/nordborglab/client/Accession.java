package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Accession implements EntryPoint {

	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
	    // Create the constants
	    AccessionConstants constants = (AccessionConstants) GWT.create(AccessionConstants.class);
	    
		TabPanel tp = new TabPanel();
		tp.add(new HTML("Accessions"), "Accessions");
		tp.add(new AccessionByName(constants), "Query By Name");
		tp.add(new HTML("Bar"), "bar");
		tp.add(new HTML("Baz"), "baz");
		tp.add(new HTML("Under construction"), "Raw SQL Query");
		
		// Show the 'bar' tab initially.
		tp.selectTab(1);
		tp.setWidth("800px");
		tp.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("accession").add(tp);
	}
}
