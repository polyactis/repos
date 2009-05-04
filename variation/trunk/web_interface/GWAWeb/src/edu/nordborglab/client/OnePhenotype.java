package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;

import com.google.gwt.user.client.ui.RootPanel;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class OnePhenotype implements EntryPoint {
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private OnePhenotypePanel onePhenotypePanel;
	private RootPanel rootPanel;
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		rootPanel = RootPanel.get("onePhenotype");
		
		onePhenotypePanel = new OnePhenotypePanel(constants, jsonErrorDialog);
		rootPanel.add(onePhenotypePanel);
		
	}


};
