package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.event.dom.client.ChangeEvent;
import com.google.gwt.event.dom.client.ChangeHandler;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.visualization.client.visualizations.Table;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class AssociationOverlap implements EntryPoint {
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	
	/**
	 * Create a remote service proxy to talk to the server-side Greeting
	 * service.
	 */
	//private final GreetingServiceAsync greetingService = GWT
	//		.create(GreetingService.class);
	
	private TabPanel tPanel;
	private String callMethodID;
	private String callMethodLsJsonURL;
	private String overlappingDataAcrossPhenotypesURL;
	private String noOfTopSNPsLsJsonURL;
	private String noOfTopSNPs;
	
	private CustomVerticalPanel summaryVPanel;
	private ListBox noOfTopSNPsListBox;
	private MotionChart motionChart;
	private DisclosurePanel motionChartDPanel;
	private WrapVisualizationTable summaryTable;
	private DataTable summaryDataTable;
	
	private HTML statusReport;
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		tPanel = new TabPanel();
		summaryVPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.AssociationOverlapSummaryHelpID());
		
		statusReport = new HTML("Waiting ...");
		statusReport.setVisible(false);
		summaryVPanel.add(statusReport);
		noOfTopSNPsListBox = new ListBox();
		summaryVPanel.add(noOfTopSNPsListBox);
		
		motionChart = new MotionChart();
		motionChart.setSize("1000px", "800px");	 	//doesn't work. have to be set in the options before draw()
		
		motionChartDPanel = new DisclosurePanel("MotionChart");
		motionChartDPanel.setAnimationEnabled(true);
		motionChartDPanel.setContent(motionChart);
		motionChartDPanel.setOpen(true);
		summaryVPanel.add(motionChartDPanel);
		
		summaryTable = new WrapVisualizationTable();
		summaryVPanel.add(summaryTable);
		
		tPanel.add(summaryVPanel, "summary");
		
		tPanel.selectTab(0);
		//tPanel.setWidth("1000px");
		tPanel.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("gwt").add(tPanel);
		
		
		callMethodLsJsonURL = getCallMethodLsJsonURL();
		overlappingDataAcrossPhenotypesURL = getOverlappingDataAcrossPhenotypesURL();
		noOfTopSNPsLsJsonURL = getNoOfTopSNPsLsJsonURL();
		
		// Create the popup dialog box
		final DialogBox dialogBox = new DialogBox();
		dialogBox.setText("Choose a Dataset");
		dialogBox.setAnimationEnabled(true);
		final ListBox callMethodListBox = new ListBox();
		
		VerticalPanel dialogVPanel = new VerticalPanel();
		dialogVPanel.setHorizontalAlignment(VerticalPanel.ALIGN_RIGHT);
		dialogVPanel.add(callMethodListBox);
		dialogBox.setWidget(dialogVPanel);
		dialogBox.setModal(true);
		dialogBox.show();
		dialogBox.center();
		// Add a handler to close the DialogBox
		callMethodListBox.addChangeHandler(new ChangeHandler() {
			public void onChange(ChangeEvent event) {
				dialogBox.setModal(false);
				dialogBox.hide();
				callMethodID = callMethodListBox.getValue(callMethodListBox.getSelectedIndex());
				Common.fillSelectBox(noOfTopSNPsLsJsonURL+"?call_method_id="+callMethodID, noOfTopSNPsListBox, jsonErrorDialog);
				//populateData(callMethodID);
			}
		});
		
		
		Common.fillSelectBox(callMethodLsJsonURL, callMethodListBox, jsonErrorDialog);
		
		noOfTopSNPsListBox.addChangeHandler(new ChangeHandler() {
			public void onChange(ChangeEvent event) {
				noOfTopSNPs = noOfTopSNPsListBox.getValue(noOfTopSNPsListBox.getSelectedIndex());
				populateData();
			}
		});
	}
	
	public native String getCallMethodLsJsonURL()/*-{ return $wnd.getCallMethodLsJsonURL; }-*/;
	public native String getOverlappingDataAcrossPhenotypesURL()/*-{ return $wnd.getOverlappingDataAcrossPhenotypesURL; }-*/;
	public native String getNoOfTopSNPsLsJsonURL()/*-{ return $wnd.getNoOfTopSNPsLsJsonURL; }-*/;
	
	public void populateData()
	{
		statusReport.setVisible(true);
		String url = URL.encode(overlappingDataAcrossPhenotypesURL+"?no_of_top_snps="+noOfTopSNPs+"&call_method_id="+callMethodID);
		Common.fillInTable(url, summaryTable, jsonErrorDialog, statusReport, summaryDataTable);
		Common.fillInMotionChart(url, motionChart, jsonErrorDialog, statusReport, summaryDataTable);
	}
	
	public void resetTitle()
	{
		statusReport.setVisible(false);
		//statusReport.setText("");
		//DOM.getElementById("title").setInnerText(TITLE_DEFAULT_TEXT);
	}
}
