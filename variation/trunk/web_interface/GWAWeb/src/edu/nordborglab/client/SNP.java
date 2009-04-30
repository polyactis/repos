package edu.nordborglab.client;


import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.core.client.JsArray;
import com.google.gwt.core.client.JsArrayString;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONArray;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONValue;

import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Frame;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;


/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class SNP implements EntryPoint {
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	
	private TabPanel tPanel;
	
	private int callMethodID;
	private int phenotypeMethodID;
	private int analysisMethodID;
	private String pageTitle;
	private JsArrayString  snpSpace;
	
	private String snpSummaryQueryURL;
	private String snpSignificantHitsQueryURL;
	private String ecotypeAllelePhenotypeURL;
	private String GBrowseURLJS;
	private String GBrowseURL;
	
	private CustomVerticalPanel vPanel;
	private Frame GBrowseFrame;
	
	private String TITLE_DEFAULT_TEXT = "SNP";
	private String TITLE_WAITING_TEXT = "Waiting...";
	
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
		
		
		// Add it to the root panel.
		RootPanel.get("snp").add(tPanel);
		
		snpSpace = getSNPSpace();
		snpSummaryQueryURL = snpSpace.get(0);
		snpSignificantHitsQueryURL =  snpSpace.get(1);
		ecotypeAllelePhenotypeURL = snpSpace.get(2);
		GBrowseURL = snpSpace.get(3);
		
		//callMethodID = getCallMethodID();
		//phenotypeMethodID = getPhenotypeMethodID();
		//analysisMethodID = getAnalysisMethodID();
		pageTitle = getPageTitle();
		TITLE_DEFAULT_TEXT = pageTitle;
		TITLE_WAITING_TEXT = TITLE_WAITING_TEXT + pageTitle;
		
		//tPanel.setSize("100%", "100%");
		tPanel.setWidth("100%");
		vPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.GBrowseHelpID());
		vPanel.setSize("100%", "100%");
		GBrowseFrame = new Frame(GBrowseURL);
		
		GBrowseFrame.setSize("100%", "800px");
		//GBrowseHTML.getElement().getId();
		vPanel.add(GBrowseFrame);
		tPanel.add(vPanel, "GBrowse");
		tPanel.selectTab(0);
		
		RootPanel SNPSummaryDiv = RootPanel.get("SNPSummary");
		RootPanel.detachNow(SNPSummaryDiv);
		CustomVerticalPanel SNPSummaryPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.SNPSummaryPanelHelpID());
		SNPSummaryPanel.add(SNPSummaryDiv);
		//SNPSummaryDiv.setVisible(true);
		tPanel.add(SNPSummaryPanel, "SNP Summary Info");
		
		RootPanel SignificantHitsInAllPhenotypeDiv = RootPanel.get("SignificantHitsInAllPhenotype");
		RootPanel.detachNow(SignificantHitsInAllPhenotypeDiv);
		CustomVerticalPanel SignificantHitsInAllPhenotypePanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.SignificantHitsInAllPhenotypeHelpID());
		SignificantHitsInAllPhenotypePanel.add(SignificantHitsInAllPhenotypeDiv);
		tPanel.add(SignificantHitsInAllPhenotypePanel, "All Phenotypes in which SNP is significant");
		
		RootPanel EcotypeAlleleMotionChartDiv = RootPanel.get("EcotypeAlleleMotionChart");
		RootPanel.detachNow(EcotypeAlleleMotionChartDiv);
		CustomVerticalPanel EcotypeAlleleMotionChartPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.EcotypeAlleleMotionChartHelpID());
		EcotypeAlleleMotionChartPanel.add(EcotypeAlleleMotionChartDiv);
		tPanel.add(EcotypeAlleleMotionChartPanel, "Ecotype Allele Phenotype BarChart");
		
		RootPanel EcotypeAlleleTableDiv = RootPanel.get("EcotypeAlleleTable");
		RootPanel.detachNow(EcotypeAlleleTableDiv);
		CustomVerticalPanel EcotypeAlleleTablePanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.EcotypeAlleleTableHelpID());
		EcotypeAlleleTablePanel.add(EcotypeAlleleTableDiv);
		tPanel.add(EcotypeAlleleTablePanel, "Ecotype Allele Phenotype Table");
		
		resetTitle();
		//fetchGBrowseHTML();
	}
	
	public native int getCallMethodID()/*-{	return $wnd.call_method_id; }-*/;

	public native int getPhenotypeMethodID()/*-{ return $wnd.phenotype_method_id; }-*/;
	public native int getAnalysisMethodID()/*-{ return $wnd.analysis_method_id; }-*/;
	public native String getPageTitle()/*-{ return $wnd.pageTitle; }-*/;
	public native JsArrayString getSNPSpace()/*-{
	return [$wnd.snpSummaryQueryURL, $wnd.snpSignificantHitsQueryURL, $wnd.ecotypeAllelePhenotypeURL, $wnd.GBrowseURL];
	}-*/;
	
	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}
		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				//GBrowseHTML.setHTML(responseText);
				tPanel.selectTab(0);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	public void fetchGBrowseHTML()
	{
		/*
		 * each category is gonna be a tab.
		 * fetch all categories from server (which further fetches from db)
		 * 		add a tab for each category
		 */
		setIntoWaitState();
		GBrowseFrame.setUrl(GBrowseURL);
		/*
		String url = URL.encode(GBrowseURL);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
		*/
		resetTitle();
	}
	
	public void setIntoWaitState()
	{
		Window.setTitle(TITLE_WAITING_TEXT);
		//editSaveButton.setEnabled(false);
	}
	
	public void resetTitle()
	{
		Window.setTitle(TITLE_DEFAULT_TEXT);
		//DOM.getElementById("title").setInnerText(TITLE_DEFAULT_TEXT);
	}
	
}
