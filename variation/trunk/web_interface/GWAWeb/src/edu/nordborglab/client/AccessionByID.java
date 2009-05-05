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


public class AccessionByID extends Sink implements ClickListener{

	private TextBox idBox = new TextBox();
	private Button submitButton = new Button();
	private HorizontalPanel submitPanel = new HorizontalPanel();
	private VerticalPanel vpanel;
	
	private DisplayJSONObject jsonErrorDialog;
	public MapTableTree contentTree;	
	
	/**
	 * An instance of the constants.
	 */
	private AccessionConstants constants;
	private static final String SUGGEST_BUTTON_DEFAULT_TEXT = "Search";
	private static final String SUGGEST_BUTTON_WAITING_TEXT = "Waiting...";

	/*
	private class SuggestBoxChangeListener implements ChangeListener {
		public void onChange(Widget sender) {
			oracle.requestSuggestions(SuggestOracle.Request request,
					new AccessionSuggestOracleCallback());
		}
	}

	private class OracleRequestSuggestionsCallBack implements SuggestOracle.Callback {
		public void onSuggestionsReady(SuggestOracle.Request request, SuggestOracle.Response response) {
			String responseText = response.toString();	//getText();
			try {
				JSONValue jsonValue = JSONParser.parse(responseText);
				displayJSONObject(jsonValue);
			} catch (JSONException e) {
				displayParseError(responseText);
			}
		}
	}
	*/
	
	
	/**
	 * Constructor.
	 * 
	 * @param constants the constants
	 */
	public AccessionByID(AccessionConstants constants, DisplayJSONObject jsonErrorDialog) {
		//super(constants);
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		vpanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.AccessionByIDHelpID());
		
		//Label lbl = new Label(constants.cwAccessionByNameLabel());
		//suggestBox.addChangeListener(new SuggestBoxChangeListener());
		submitButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
		submitButton.addClickListener(this);

		submitPanel.add(new HTML(constants.cwAccessionByIDLabel()));
		idBox.setTitle(constants.cwAccessionByIDDescription());
		submitPanel.add(idBox);
		
		submitPanel.add(submitButton);
		submitPanel.setSpacing(5);
		
		vpanel.add(submitPanel);
		//panel.add(textBox);
		
		contentTree = new MapTableTree(constants, jsonErrorDialog);
		vpanel.add(contentTree);
		
		
		// All composites must call initWidget() in their constructors.
		initWidget(vpanel);

		// Give the overall composite a style name.
		setStyleName("AccessionByID");
	}

	public void onClick(Widget sender) {
		doFetchURL();
	}
	
	
	//@Override
	public String getDescription() {
		return constants.cwAccessionByIDDescription();
	}
	
	@Override
	public HTML getDescriptionHTML() {
		return new HTML("<p>" + getDescription() + "</p>");
	}
	
	@Override
	public String getName() {
		return constants.cwAccessionByIDName();
	}
	
	private final native DataTable asDataTable(String json) /*-{
		dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
		return dataTable;
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
		submitButton.setText(SUGGEST_BUTTON_WAITING_TEXT);
		String url = URL.encode(constants.AccessionByIDURL() + idBox.getText());
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
		
	}

	private void resetSearchButtonCaption() {
		submitButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
	}
	public void resetSize()
	{
		contentTree.resetSize();
	}
	
}
