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
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.Query;
import com.google.gwt.visualization.client.QueryResponse;
import com.google.gwt.visualization.client.Query.Callback;

import java.util.Set;

public class AccessionByName extends Composite implements ClickListener{

	private AccessionSuggestOracle oracle = new AccessionSuggestOracle();
	//String[] words = constants.cwSuggestBoxWords();
	//for (int i = 0; i < words.length; ++i) {
	//	oracle.add(words[i]);
	//}

	// Create the suggest box
	private SuggestBox suggestBox = new SuggestBox(oracle);
	private Button suggestButton = new Button();
	//suggestBox.ensureDebugId("cwSuggestBox");
	private HorizontalPanel suggestPanel = new HorizontalPanel();
	private VerticalPanel panel = new VerticalPanel();
	private Tree jsonTree = new Tree();
	
	private Table accessionTable = new Table();
	//private String dataUrl = "http://spreadsheets.google.com/tq?key=prll1aQH05yQqp_DKPP9TNg&pub=1";
	//private Query query = Query.create(dataUrl);
	
	
	private DialogBox dialogBox = new DialogBox();

	public void onClick(Widget sender) {
		//jsonTree.setVisible(false);
		doFetchURL();
	}

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


	private class SuggestBoxSuggestionHandler implements SuggestionHandler {
		public void onSuggestionSelected(SuggestionEvent event){
			jsonTree.setVisible(false);
			doFetchURL();
		}
	}
	/**
	 * Constructor.
	 * 
	 * @param constants the constants
	 */
	public AccessionByName(AccessionConstants constants) {
		//super(constants);
		this.constants = constants;
		oracle.setConstants(constants);
		
		//Label lbl = new Label(constants.cwAccessionByNameLabel());
		//suggestBox.addChangeListener(new SuggestBoxChangeListener());
		suggestBox.addEventHandler(new SuggestBoxSuggestionHandler());

		suggestPanel.add(new HTML(constants.cwAccessionByNameLabel()));
		suggestPanel.add(suggestBox);
		suggestPanel.add(suggestButton);
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
		suggestButton.addClickListener(this);

		panel.add(suggestPanel);
		//panel.add(textBox);

		panel.add(jsonTree);
		//Avoids showing an "empty" cell
		jsonTree.setVisible(false);

		panel.add(accessionTable);
		
		/*
		 * 2009-4-3 construction and more has to be put into constructors, not under the class.
		 */
		dialogBox.setText("Welcome to GWT!");
		dialogBox.setAnimationEnabled(true);
		Button closeButton = new Button("close");
		VerticalPanel dialogVPanel = new VerticalPanel();
		dialogVPanel.setWidth("100%");
		dialogVPanel.setHorizontalAlignment(VerticalPanel.ALIGN_CENTER);
		dialogVPanel.add(closeButton);

		closeButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				dialogBox.hide();
			}
		});
		
		// All composites must call initWidget() in their constructors.
		initWidget(panel);

		// Give the overall composite a style name.
		setStyleName("AccessionByName");
	}

	//@Override
	public String getDescription() {
		return constants.cwAccessionByNameDescription();
	}

	//@Override
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
			displayRequestError(exception.toString());
			resetSearchButtonCaption();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				//JSONValue jsonValue = JSONParser.parse(responseText);
				//displayJSONObject(jsonValue);
				DataTable data = asDataTable(responseText);	//DataTable.create();//new google.visualization.DataTable(eval("("+response+")"), 0.5)
				Table.Options options = Table.Options.create();
				options.setShowRowNumber(true);
				accessionTable.draw(data, options);
			} catch (JSONException e) {
				displayParseError(responseText);
			}
			resetSearchButtonCaption();
		}
	}

	/*
	 * Add the object presented by the JSONValue as a children to the requested
	 * TreeItem.
	 */
	private void addChildren(TreeItem treeItem, JSONValue jsonValue) {
		JSONArray jsonArray;
		JSONObject jsonObject;
		JSONString jsonString;

		if ((jsonArray = jsonValue.isArray()) != null) {
			for (int i = 0; i < jsonArray.size(); ++i) {
				TreeItem child = treeItem.addItem(getChildText("["
						+ Integer.toString(i) + "]"));
				addChildren(child, jsonArray.get(i));
			}
		} else if ((jsonObject = jsonValue.isObject()) != null) {
			Set<String> keys = jsonObject.keySet();
			for (String key : keys) {
				TreeItem child = treeItem.addItem(getChildText(key));
				addChildren(child, jsonObject.get(key));
			}
		} else if ((jsonString = jsonValue.isString()) != null) {
			// Use stringValue instead of toString() because we don't want escaping
			treeItem.addItem(jsonString.stringValue());
		} else {
			// JSONBoolean, JSONNumber, and JSONNull work well with toString().
			treeItem.addItem(getChildText(jsonValue.toString()));
		}
	}

	/*
	 * Causes the text of child elements to wrap.
	 */
	private String getChildText(String text) {
		return "<span style='white-space:normal'>" + text + "</span>";
	}

	private void displayJSONObject(JSONValue jsonValue) {
		jsonTree.removeItems();
		jsonTree.setVisible(true);
		TreeItem treeItem = jsonTree.addItem("JSON Response");
		addChildren(treeItem, jsonValue);
		treeItem.setStyleName("JSON-JSONResponseObject");
		treeItem.setState(true);
	}

	private void displayParseError(String responseText) {
		displayError("Failed to parse JSON response", responseText);
	}

	private void displayRequestError(String message) {
		displayError("Request failed.", message);
	}

	private void displaySendError(String message) {
		displayError("Failed to send the request.", message);
	}

	private void displayError(String errorType, String errorMessage) {
		jsonTree.removeItems();
		jsonTree.setVisible(true);
		TreeItem treeItem = jsonTree.addItem(errorType);
		treeItem.addItem(errorMessage);
		treeItem.setStyleName("JSON-JSONResponseObject");
		treeItem.setState(true);
	}


	private void doFetchURL() {
		suggestButton.setText(SUGGEST_BUTTON_WAITING_TEXT);
		String url = URL.encode(constants.AccessionByNameURL() + "?name=" + suggestBox.getText());
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
		/*
		query.send(new Callback() {
			public void onResponse(QueryResponse response) {
				if (response.isError()) {
					dialogBox.setText("Error in query: " + response.getMessage() + ' '
							+ response.getDetailedMessage());
					dialogBox.center();
					dialogBox.show();
					return;
				}
				Table.Options options = Table.Options.create();
				options.setShowRowNumber(true);
				accessionTable.draw(response.getDataTable(), options);
			}
		});
		*/
		
	}

	private void resetSearchButtonCaption() {
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
	}

}
