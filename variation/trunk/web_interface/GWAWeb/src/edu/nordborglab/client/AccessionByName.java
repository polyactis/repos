package edu.nordborglab.client;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HTML;

import com.google.gwt.user.client.ui.SuggestBox;

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

import com.google.gwt.user.client.ui.SuggestionHandler;
import com.google.gwt.user.client.ui.SuggestionEvent;

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


public class AccessionByName extends Sink implements ClickListener{

	private CustomSuggestOracle oracle;
	//String[] words = constants.cwSuggestBoxWords();
	//for (int i = 0; i < words.length; ++i) {
	//	oracle.add(words[i]);
	//}

	// Create the suggest box
	private SuggestBox suggestBox;
	private Button suggestButton;
	//suggestBox.ensureDebugId("cwSuggestBox");
	private HorizontalPanel suggestPanel;
	private VerticalPanel panel;
	
	//private String dataUrl = "http://spreadsheets.google.com/tq?key=prll1aQH05yQqp_DKPP9TNg&pub=1";
	//private Query query = Query.create(dataUrl);
	private DisplayJSONObject jsonErrorDialog;
	private MapTableTree contentTree;
	

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
	public AccessionByName(AccessionConstants constants, DisplayJSONObject jsonErrorDialog) {
		//super(constants);
		this.constants = constants;
		oracle = new CustomSuggestOracle(this.constants.AccessionSuggestOracleURL() + "?namelike=");
		
		this.jsonErrorDialog = jsonErrorDialog;
		
		panel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.AccessionByNameHelpID());
		
		suggestBox = new SuggestBox(oracle);
		//Label lbl = new Label(constants.cwAccessionByNameLabel());
		//suggestBox.addChangeListener(new SuggestBoxChangeListener());
		suggestBox.addEventHandler(new SuggestBoxSuggestionHandler());
		suggestBox.setTitle(constants.cwAccessionByNameDescription());
		
		suggestButton = new Button();
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
		suggestButton.addClickListener(this);
		
		suggestPanel = new HorizontalPanel();
		suggestPanel.add(new HTML(constants.cwAccessionByNameLabel()));
		suggestPanel.add(suggestBox);
		suggestPanel.add(suggestButton);
		suggestPanel.setSpacing(5);
		
		panel.add(suggestPanel);
		//panel.add(textBox);
		
		contentTree = new MapTableTree(constants, jsonErrorDialog);
		panel.add(contentTree);
		
		// All composites must call initWidget() in their constructors.
		initWidget(panel);

		// Give the overall composite a style name.
		setStyleName("AccessionByName");
	}

	public void onClick(Widget sender) {
		doFetchURL();
	}
	
	private class SuggestBoxSuggestionHandler implements SuggestionHandler {
		public void onSuggestionSelected(SuggestionEvent event){
			doFetchURL();
		}
	}
	
	//@Override
	public String getDescription() {
		return constants.cwAccessionByNameDescription();
	}
	
	@Override
	public HTML getDescriptionHTML() {
		return new HTML("<p>" + getDescription() + "</p>");
	}
	
	@Override
	public String getName() {
		return constants.cwAccessionByNameName();
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
		suggestButton.setText(SUGGEST_BUTTON_WAITING_TEXT);
		String url = URL.encode(constants.AccessionByNameURL() + suggestBox.getText());
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
		
	}

	private void resetSearchButtonCaption() {
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
	}
	public void resetSize()
	{
		contentTree.resetSize();
	}
}
