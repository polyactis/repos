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
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONValue;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;


/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class OnePhenotypeGWASGene implements EntryPoint {

	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	private ArrayList<ArrayList<String>> phenotypeCategoryLs;
	
	private OnePhenotypePanel onePhenotypePanel;
	
	private TabPanel tPanel;
	
	
	private static final String TITLE_DEFAULT_TEXT = "Phenotype";
	private static final String TITLE_WAITING_TEXT = "Waiting...";
	
	private String snpGeneAssociationTypeID;
	private String callMethodID;
	private String phenotypeMethodID;
	private String getAnalysisMethodLsURL;
	
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
		RootPanel.get("onePhenotypeGWASGene").add(tPanel);
		
		
		onePhenotypePanel = new OnePhenotypePanel(constants, jsonErrorDialog);
		
		// 2009-5-2 not necessary, use OnePhenotypePanel
		//move the onePhenotype division into tabPanel
		//RootPanel onePhenotypeDiv = RootPanel.get("onePhenotype");
		//RootPanel.detachNow(onePhenotypeDiv);	//2009-4-23 detach the element now in order to be attached in another place.
		
		tPanel.add(onePhenotypePanel, "Phenotype");
		
		tPanel.selectTab(0);
		
		snpGeneAssociationTypeID = Common.getSNPGeneAssociationTypeID();
		callMethodID = Common.getCallMethodID();
		phenotypeMethodID = Common.getPhenotypeMethodID();
		getAnalysisMethodLsURL = Common.getAnalysisMethodLsURL();
		
		addGeneAssociationTabs();
	}
	
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
						Integer analysisMethodIDInt = (Integer) analysisMethodID; 
						OneGWASGene oneGWASGene = new OneGWASGene(constants, jsonErrorDialog, callMethodID, phenotypeMethodID, 
								analysisMethodIDInt.toString(), snpGeneAssociationTypeID);
						tPanel.add(oneGWASGene, analysisMethodName);
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
	
	public void addGeneAssociationTabs()
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
