package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

import com.google.gwt.visualization.client.DataTable;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class HaploGroup implements EntryPoint {

	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	private MapTableTree contentTree;
	
	private static final String SUGGEST_BUTTON_DEFAULT_TEXT = "Search";
	private static final String SUGGEST_BUTTON_WAITING_TEXT = "Waiting...";
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		
		VerticalPanel vPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.HaploGroupHelpID());
		vPanel.setWidth("100%");
		//vPanel.setHorizontalAlignment(VerticalPanel.ALIGN_CENTER);
		
		contentTree = new MapTableTree(constants, jsonErrorDialog);
		//contentTree.populateData(getDataTable());
		vPanel.add(contentTree);
		
		// Add image and button to the RootPanel
		RootPanel.get("gwt").add(vPanel);
		doFetchURL();
	}
	public native int getHaploGroupID()/*-{
		return $wnd.haplo_group_id;
	}-*/;
	public native String getHaploGroupURL()/*-{
		return $wnd.getHaploGroupURL;
	}-*/;
	
	// 2009-4-19 getDataTable() is useless because dataTable is now fetched thru asynchronous request.
	public native DataTable getDataTable()/*-{
		dataTable = new $wnd.google.visualization.DataTable($wnd.dataTable, 0.5);
		return dataTable; 
	}-*/;
	
	
	private final native DataTable asDataTable(String json) /*-{
		dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
		return dataTable;
	}-*/;

	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			contentTree.resetTablePanelHeaderText();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				//JSONValue jsonValue = JSONParser.parse(responseText);
				//displayJSONObject(jsonValue);
				DataTable data = asDataTable(responseText);	//DataTable.create();//new google.visualization.DataTable(eval("("+response+")"), 0.5)
				contentTree.populateData(data);
				
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			contentTree.resetTablePanelHeaderText();
		}
	}

	private void doFetchURL() {
		contentTree.setTablePanelHeaderText(SUGGEST_BUTTON_WAITING_TEXT);
		contentTree.setMapPanelHeaderText(SUGGEST_BUTTON_WAITING_TEXT);	//2009-4-19 no need to reset it cuz it's reset in contentTree.populateData().
		String url = URL.encode(getHaploGroupURL());
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			contentTree.resetTablePanelHeaderText();
		}
		
	}
}
