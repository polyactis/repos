package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.Response;
import com.google.gwt.json.client.JSONException;

import com.google.gwt.user.client.DOM;
import com.google.gwt.user.client.Window;

import com.google.gwt.user.client.ui.DecoratedPopupPanel;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.TabPanel;


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
import com.google.gwt.json.client.JSONString;
import com.google.gwt.json.client.JSONValue;



import java.util.ArrayList;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class GWASPhenotypes implements EntryPoint {

	/**
	 * This is the entry point method.
	 */
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	private ArrayList<ArrayList<String>> phenotypeCategoryLs;
	
	private TabPanel tPanel;
	
	
	private static final String TITLE_DEFAULT_TEXT = "All Phenotypes with GWAS ...";
	private static final String TITLE_WAITING_TEXT = "Waiting all phenotypes ...";
	
	private int callMethodID;
	private String fetchPhenotypeCategoryLsURL;
	private String fetchPhenotypeTableDataURL;
	private String fetchGWAURL;
	private DecoratedPopupPanel popupLink = new DecoratedPopupPanel();
	
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		tPanel = new TabPanel();
		
		//tPanel.selectTab(1);
		//tPanel.setWidth("1000px");
		tPanel.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("gwt").add(tPanel);
		
		callMethodID = getCallMethodID();
		fetchPhenotypeCategoryLsURL = getPhenotypeCategoryLsURL();
		fetchPhenotypeTableDataURL = getPhenotypeTableDataURL();
		fetchGWAURL = getGWAURL();
		fillInPhenotypeCategories();
		
	}
	public native int getCallMethodID()/*-{
		return $wnd.call_method_id;
	}-*/;
	
	public native String getPhenotypeCategoryLsURL()/*-{
		return $wnd.getPhenotypeCategoryLsURL;
	}-*/;
	
	public native String getPhenotypeTableDataURL()/*-{
		return $wnd.getPhenotypeTableDataURL;
	}-*/;
	
	public native String getGWAURL()/*-{
		return $wnd.getGWAURL;
	}-*/;
	
	
	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONValue jsonValue = JSONParser.parse(response.getText());				
				JSONArray array = jsonValue.isArray();	//jsonValue.isArray();
				//Window.alert("json gets: "+ array);
				for (int i = 0; i < array.size(); i++) {
					//JSONObject object = array.get(i).isObject();
					//String movieName = object.isString().stringValue();
					JSONArray oneCategoryTuple = array.get(i).isArray();
					String phenotypeCategoryID = oneCategoryTuple.get(0).isString().stringValue();
					String phenotypeCategoryName = oneCategoryTuple.get(1).isString().stringValue();
					addOnePhenotypeCategory(phenotypeCategoryID, phenotypeCategoryName);
				}
				tPanel.selectTab(0);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	public void fillInPhenotypeCategories()
	{
		Window.setTitle(TITLE_WAITING_TEXT);
		//DOM.getElementById("title").setInnerText(TITLE_WAITING_TEXT);
		String url = URL.encode(fetchPhenotypeCategoryLsURL);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
	}
	
	public void addOnePhenotypeCategory(String phenotypeCategoryID, String phenotypeCategoryName)
	{
		PhenotypeTable pTable = new PhenotypeTable(constants, jsonErrorDialog, popupLink, phenotypeCategoryID, phenotypeCategoryName, 
				fetchPhenotypeTableDataURL, callMethodID, fetchGWAURL);
		tPanel.add(pTable, phenotypeCategoryName);
		
	}
	
	public void resetTitle()
	{
		Window.setTitle(TITLE_DEFAULT_TEXT);
		//DOM.getElementById("title").setInnerText(TITLE_DEFAULT_TEXT);
	}
}