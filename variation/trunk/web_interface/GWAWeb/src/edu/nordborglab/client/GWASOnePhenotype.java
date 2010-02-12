package edu.nordborglab.client;

import java.util.ArrayList;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONArray;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONValue;
import com.google.gwt.json.client.JSONObject;

import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DecoratedPopupPanel;
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
public class GWASOnePhenotype implements EntryPoint {

	
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	private ArrayList<ArrayList<String>> phenotypeCategoryLs;
	
	private OnePhenotypePanel onePhenotypePanel;
	
	private TabPanel tPanel;
	
	
	private static final String TITLE_DEFAULT_TEXT = "Phenotype";
	private static final String TITLE_WAITING_TEXT = "Waiting...";
	
	private int callMethodID;
	private int phenotypeMethodID;
	private String GWABaseURL;
	private String SNPBaseURL;
	private String getAnalysisMethodLsURL;
	
	
	private DecoratedPopupPanel popupLink = new DecoratedPopupPanel();
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
		RootPanel.get("gwt").add(tPanel);
		
		
		onePhenotypePanel = new OnePhenotypePanel(constants, jsonErrorDialog);
		
		// 2009-5-2 not necessary, use OnePhenotypePanel
		//move the onePhenotype division into tabPanel
		//RootPanel onePhenotypeDiv = RootPanel.get("onePhenotype");
		//RootPanel.detachNow(onePhenotypeDiv);	//2009-4-23 detach the element now in order to be attached in another place.
		
		tPanel.add(onePhenotypePanel, "Phenotype");
		
		tPanel.selectTab(0);
		
		
		callMethodID = getCallMethodID();
		phenotypeMethodID = getPhenotypeMethodID();
		GWABaseURL = getGWABaseURL();
		SNPBaseURL = getSNPBaseURL();
		getAnalysisMethodLsURL = getAnalysisMethodLsURL();
		
		addAssociationResults();
	}
	public native int getCallMethodID()/*-{
		return $wnd.call_method_id;
	}-*/;
	
	public native int getPhenotypeMethodID()/*-{
		return $wnd.phenotype_method_id;
	}-*/;
	
	public native String getGWABaseURL()/*-{ return $wnd.GWABaseURL}-*/;
	public native String getSNPBaseURL()/*-{ return $wnd.SNPBaseURL}-*/;
	public native String getAnalysisMethodLsURL()/*-{ return $wnd.getAnalysisMethodLsURL}-*/;
	
	
	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}
		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONValue jsonValue = JSONParser.parse(response.getText());				
				JSONObject jsonObject = jsonValue.isObject();
				JSONArray analysisMethodLs = jsonObject.get("options").isArray();
				for (int i = 0; i < analysisMethodLs.size(); i++) {
					JSONObject object = analysisMethodLs.get(i).isObject();
					int analysisMethodID = (int) object.get("value").isNumber().doubleValue();
					if (analysisMethodID != 0 )
					{
						
						String analysisMethodName = object.get("id").isString().stringValue();
						String analysisMethodDescription = object.get("description").isString().stringValue();
						//int analysis_method_id = Integer.parseInt(analysisMethodName);
						GWASOneResult gwaOneResult = new GWASOneResult(constants, jsonErrorDialog, popupLink, 
								analysisMethodID, analysisMethodDescription, GWABaseURL, SNPBaseURL+"&analysis_method_id="+analysisMethodID);
						tPanel.add(gwaOneResult, analysisMethodName);
					}
				}
				tPanel.selectTab(0);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	public void addAssociationResults()
	{
		/*
		 * each category is gonna be a tab.
		 * fetch all categories from server (which further fetches from db)
		 * 		add a tab for each category
		 */
		setIntoWaitState();
		String url = URL.encode(getAnalysisMethodLsURL);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
	}
	
	public void setIntoWaitState()
	{
		//Window.setTitle(TITLE_WAITING_TEXT);
		//editSaveButton.setEnabled(false);
	}
	
	public void resetTitle()
	{
		//Window.setTitle(TITLE_DEFAULT_TEXT);
		//DOM.getElementById("title").setInnerText(TITLE_DEFAULT_TEXT);
	}
}

