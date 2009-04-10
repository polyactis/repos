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
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.AbstractVisualization.VisualizationFactory;

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
		DisplayJSONObject jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		/*
		AbstractVisualization.registerVisualization("MapWithPhenotype",
				new VisualizationFactory() {

			public AbstractVisualization<?> create() {
				return new MapWithPhenotype(constants, jsonErrorDialog);
			}
		});
		*/

		TabPanel tp = new TabPanel();
		tp.add(new HTML("Accessions"), "Accessions");
		tp.add(new AccessionByName(constants, jsonErrorDialog), "By Name");
		tp.add(new HTML("Under construction"), "By EcotypeID");
		tp.add(new HTML("Under construction"), "By Genetic Distance");
		tp.add(new HTML("Under construction"), "By Geographic Distance");
		tp.add(new HTML("Under construction"), "By Country");
		// Show the 'bar' tab initially.
		tp.selectTab(1);
		tp.setWidth("1000px");
		tp.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("accession").add(tp);
	}
}
