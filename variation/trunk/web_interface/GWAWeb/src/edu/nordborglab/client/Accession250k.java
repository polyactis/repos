package edu.nordborglab.client;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HTML;


import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;


import com.google.gwt.visualization.client.DataTable;
/*
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.core.client.GWT;
import com.google.gwt.core.client.JavaScriptException;
import com.google.gwt.user.client.Window;
import com.google.gwt.visualization.client.AbstractDataTable.ColumnType;
import java.util.Date;

import com.google.gwt.visualization.client.Query;
import com.google.gwt.visualization.client.QueryResponse;
import com.google.gwt.visualization.client.Query.Callback;
*/


public class Accession250k extends Sink{
	private VerticalPanel vpanel;
	
	private DisplayJSONObject jsonErrorDialog;
	private MapTableTree contentTree;	
	private String find250kAccessionsURL;
	/**
	 * An instance of the constants.
	 */
	private AccessionConstants constants;
	
	private static final String SUGGEST_BUTTON_DEFAULT_TEXT = "Search";
	private static final String SUGGEST_BUTTON_WAITING_TEXT = "Waiting...";
	/**
	 * Constructor.
	 * 
	 * @param constants the constants
	 */
	public Accession250k(AccessionConstants constants, DisplayJSONObject jsonErrorDialog) {
		//super(constants);
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		vpanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.Accession250kHelpID());
		
		find250kAccessionsURL = this.getFetchURL();
		
		contentTree = new MapTableTree(constants, jsonErrorDialog);
		vpanel.add(contentTree);
		
		// All composites must call initWidget() in their constructors.
		initWidget(vpanel);

		// Give the overall composite a style name.
		setStyleName("Accession250k");
		doFetchURL();
	}
	
	//@Override
	public String getDescription() {
		return constants.cwAccession250kDescription();
	}
	
	@Override
	public HTML getDescriptionHTML() {
		return new HTML("<p>" + getDescription() + "</p>");
	}
	
	@Override
	public String getName() {
		return constants.cwAccession250kName();
	}
	
	private final native DataTable asDataTable(String json) /*-{
		dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
		return dataTable;
	}-*/;
	
	private final native String getFetchURL() /*-{
		return $wnd.find250kAccessionsURL; 
	}-*/;
	
	/**
	 * Class for handling the response text associated with a request for a JSON
	 * object.
	 * 
	 */
	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetSearchButtonCaption();
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
			resetSearchButtonCaption();
		}
	}

	private void doFetchURL() {
		contentTree.setTablePanelHeaderText(SUGGEST_BUTTON_WAITING_TEXT);
		String url = URL.encode(find250kAccessionsURL);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
		
	}

	private void resetSearchButtonCaption() {
		contentTree.resetTablePanelHeaderText();
	}
	public void resetSize()
	{
		contentTree.resetSize();	//2009-4-17 this would avoid map partial coverup.		
	}
}
