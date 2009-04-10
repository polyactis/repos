package edu.nordborglab.client;

import java.util.Set;

import com.google.gwt.json.client.JSONArray;
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONString;
import com.google.gwt.json.client.JSONValue;

import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.VerticalPanel;

import com.google.gwt.user.client.ui.ClickListener;

public class DisplayJSONObject extends DialogBox{
	private Tree jsonTree = new Tree();
	private VerticalPanel dialogVPanel = new VerticalPanel();
	
	public DisplayJSONObject(String dialogTitle)
	{
		// Set the dialog box's caption.
		setText(dialogTitle);
		setSize("600px", "600px");
		
		dialogVPanel.add(jsonTree);
		//Avoids showing an "empty" cell
		jsonTree.setVisible(false);
		
		Button closeButton = new Button("Close");
		closeButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				DisplayJSONObject.this.hide();
			}
		});
		dialogVPanel.add(closeButton);
		setWidget(dialogVPanel);
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

	public void displayJSONObject(JSONValue jsonValue) {
		this.show();
		jsonTree.removeItems();
		jsonTree.setVisible(true);
		TreeItem treeItem = jsonTree.addItem("JSON Response");
		addChildren(treeItem, jsonValue);
		treeItem.setStyleName("JSON-JSONResponseObject");
		treeItem.setState(true);
	}

	public void displayParseError(String responseText) {
		displayError("Failed to parse JSON response", responseText);
	}

	public void displayRequestError(String message) {
		displayError("Request failed.", message);
	}

	public void displaySendError(String message) {
		displayError("Failed to send the request.", message);
	}

	public void displayError(String errorType, String errorMessage) {
		this.show();
		jsonTree.removeItems();
		jsonTree.setVisible(true);
		TreeItem treeItem = jsonTree.addItem(errorType);
		treeItem.addItem(errorMessage);
		treeItem.setStyleName("JSON-JSONResponseObject");
		treeItem.setState(true);
	}

}
