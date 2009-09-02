/*
 * 2009-8-31 a utilities framework to hold several applications.
 * 
 */
package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Utils implements EntryPoint {
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	private MotionChartAppMCPanel motionChartFileUploadPanel;
	private String motionChartFormActionURL;
	private TabPanel tPanel;
	
	public native String get_motionChartFormActionURL()/*-{ return $wnd.motionChartFormActionURL; }-*/;
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		tPanel = new TabPanel();
		
		//tPanel.selectTab(1);
		//tPanel.setWidth("1100px");
		tPanel.setAnimationEnabled(true);
		
		motionChartFormActionURL = get_motionChartFormActionURL();
		motionChartFileUploadPanel = new MotionChartAppMCPanel(constants, jsonErrorDialog, motionChartFormActionURL);
		tPanel.add(motionChartFileUploadPanel, "MotionChart App");
		tPanel.selectTab(0);
		
		// Add it to the root panel.
		RootPanel.get("gwt").add(tPanel);
		
	}
}
