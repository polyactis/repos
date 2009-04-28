package edu.nordborglab.client;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DecoratorPanel;
import com.google.gwt.user.client.ui.FlexTable;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.RichTextArea;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.PopupPanel;
import com.google.gwt.user.client.ui.VerticalPanel;


import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.PopupListener;

import com.google.gwt.user.client.ui.KeyboardListener;
import com.google.gwt.user.client.ui.MouseListener;


import edu.nordborglab.module.text.RichTextToolbar;

public class HelpDialog extends DialogBox implements ClickListener{
	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	private String helpID;
	
	private Button closeButton;
	private Button editSaveButton;

	private HorizontalPanel hPanel;
	private VerticalPanel vPanel;
	
	private DecoratorPanel wrapHelpPanel;
	private FlexTable helpLayout;
	private RichTextArea helpTextArea;
	private RichTextToolbar toolbar;
	
	private String fetchURL;
	private String saveURL;

	private String helpContent;
	private int contentSaved=0;	//flag indicating whether helpContent is correctly saved in the db.

	private final static String editCaption = "edit";
	private static final String saveCaption = "save";
	private static final String viewCaption = "view";

	private String WAITING_TEXT = "";
	private String DEFAULT_TEXT = "";
	
	private int maxWidth = 400;
	
	HelpDialog(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, String helpID)
	{
		super(false, false);	//false for both arguments in DialogBox(boolean autoHide, boolean modal)
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		this.helpID = helpID;

		this.fetchURL = this.constants.HelpDialogFetchURL();
		this.saveURL = this.constants.HelpDialogSaveURL();

		WAITING_TEXT = "Saving the help of "+helpID + "...";
		// Set the dialog box's caption.
		DEFAULT_TEXT = "Help for " + this.helpID; 
		setText(DEFAULT_TEXT);
		
		
		helpLayout = new FlexTable();
		
		helpTextArea = new RichTextArea();
		helpTextArea.addKeyboardListener(new HelpTextAreaKeyboardListener());
		helpTextArea.addMouseListener(new HelpTextAreaMouseListener());
		helpTextArea.setEnabled(false);	//disable the text area in the beginning, just let user to view. maybe should use setFocus().
		helpTextArea.setSize("100%", "14em");
		helpTextArea.setVisible(false);
		toolbar = new RichTextToolbar(helpTextArea);
		toolbar.setVisible(false);	//toolbar not visible in the beginning.
		
		helpLayout.setWidget(0, 0, toolbar);
		helpLayout.setWidget(1, 0, helpTextArea);
		
		wrapHelpPanel = new DecoratorPanel();
		wrapHelpPanel.setWidget(helpLayout);
		
		vPanel = new VerticalPanel();		
		vPanel.add(wrapHelpPanel);
		this.addPopupListener(new DialogPopupListener());

		closeButton = new Button("close");
		closeButton.addClickListener(this);

		editSaveButton = new Button(editCaption);
		editSaveButton.addClickListener(new EditSaveClickListener());

		hPanel = new HorizontalPanel();
		hPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);
		hPanel.add(editSaveButton);
		hPanel.add(closeButton);
		vPanel.add(hPanel);
		setWidget(vPanel);
		fetchHelpContent();
		//if (this.getOffsetWidth()>maxWidth)
		this.setWidth(maxWidth+"px");
	}
	public void onClick(Widget sender)
	{
		/*
		 * close button is hit. hide the window.
		 */
		this.hide();
	}
	private void transition2ViewState(){
		toolbar.setVisible(false);
		helpTextArea.setEnabled(false);		
		helpTextArea.setVisible(false);
		editSaveButton.setText(editCaption);
		helpLayout.setHTML(2, 0, helpContent);
		//if (this.getOffsetWidth()>maxWidth)
		this.setWidth(maxWidth+"px");
	}
	
	private void transition2EditState(){
		helpTextArea.setEnabled(true);
		helpTextArea.setVisible(true);
		toolbar.setVisible(true);
		editSaveButton.setText(viewCaption);
		helpLayout.removeRow(2);
		//helpLayout.setHTML(2, 0, helpTextArea.getHTML());
	}
	
	private class DialogPopupListener implements PopupListener
	{
		public void onPopupClosed(PopupPanel sender, boolean autoClosed)
		{
			//save the help content if it's modified or ask the user whether he/she wants to save it?
			//
			if (editSaveButton.getText().equals(saveCaption))
			{
				saveHelpContent();
				transition2ViewState();
			}
			else if (editSaveButton.getText().equals(viewCaption))
			{
				transition2ViewState();
			}
		}
	}

	private class EditSaveClickListener implements ClickListener
	{
		public void onClick(Widget sender)
		{
			/*
			 * if its caption is 'edit' which means it's in viewing mode, ready to switch to edit mode
			 * 		1. add toolbar, 2. make text area editable, 3. change the caption to 'view', 4. caption would be changed to 'save' upon content change
			 * elif its caption is 'view' (in edit mode nothing is changed or change is saved, switch to view mode)
			 * 		1. remove toolbar, 2. make text area uneditable, 3. change the caption to 'edit',
			 * elif its caption is 'save' (in edit mode, something is changed, to save stuff)
			 * 		1. disable the button, 2. save the content, 3. enable the button, 4. change its caption to view
			 */
			
			if (editSaveButton.getText().equals(editCaption))
			{
				transition2EditState();				
			}
			else if(editSaveButton.getText().equals(saveCaption))
			{
				editSaveButton.setEnabled(false);
				saveHelpContent();
				editSaveButton.setEnabled(true);
				if (contentSaved==1)
					editSaveButton.setText(viewCaption);
			}
			else if(editSaveButton.getText().equals(viewCaption))
			{
				transition2ViewState();
			}
		}
	}

	private class HelpTextAreaKeyboardListener implements KeyboardListener
	{
		public void onClick(Widget sender)
		{

		}
		public void onKeyDown(Widget sender, char keyCode, int modifiers)
		{

		}
		public void onKeyPress(Widget sender, char keyCode, int modifiers) {

		}

		public void onKeyUp(Widget sender, char keyCode, int modifiers) {
			if (sender == helpTextArea) {
				/*
				 * something is entered in the text area. change the button caption to save.
				 */
				if (!helpTextArea.getHTML().equals(helpContent))	//set the caption into saveCaption only if the content is changed
					editSaveButton.setText(saveCaption);
				else
					editSaveButton.setText(viewCaption);	//set the button caption to viewCaption if no content is changed
			}
		}
	}
	
	private class HelpTextAreaMouseListener implements MouseListener
	{
		public void onMouseMove(Widget sender, int x, int y) 
		{
		}
		public void onMouseDown(Widget sender, int x, int y)
		{

		}
		public void onMouseEnter(Widget sender)
		{

		}
		public void onMouseLeave(Widget sender) {
		
		}

		public void onMouseUp(Widget sender, int x, int y) {
			if (sender == helpTextArea) {
				/*
				 * something is entered in the text area. change the button caption to save.
				 */
				if (!helpTextArea.getHTML().equals(helpContent))	//set the caption into saveCaption only if the content is changed
					editSaveButton.setText(saveCaption);
				else
					editSaveButton.setText(viewCaption);	//set the button caption to viewCaption if no content is changed
			}
		}
	}
	
	private class FetchHelpContentHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				helpContent = responseText;
				helpTextArea.setHTML(helpContent);
				helpLayout.setHTML(2, 0, helpContent);
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}

	private void fetchHelpContent() {
		setIntoWaitState();
		String url = URL.encode(this.fetchURL + "?helpID="+this.helpID);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new FetchHelpContentHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}

	}
	
	private class SaveHelpContentHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				if (responseText.equals("1"))
				{
					helpContent = helpTextArea.getHTML();	//update the helpContent which holds the text
					contentSaved = 1;
				}
				else
				{
					jsonErrorDialog.displayError("Help Content Save Error", responseText);
					contentSaved = 0;
				}

			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	private void saveHelpContent()
	{
		setIntoWaitState();
		if (!helpTextArea.getHTML().equals(helpContent))	//only save when the content is different from original
		{
			String url = URL.encode(this.saveURL + "?helpID="+this.helpID);
			RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.POST, url);
			//requestBuilder.setHeader("helpContent", helpTextArea.getHTML());
			try {
				requestBuilder.sendRequest(helpTextArea.getHTML(), new SaveHelpContentHandler());	//1st argument is requestData, which appears as request.body in pylons.
			} catch (RequestException ex) {
				jsonErrorDialog.displaySendError(ex.toString());
				resetTitle();
			}
		}
		resetTitle();
	}
	
	public void setIntoWaitState()
	{
		setText(WAITING_TEXT);
		editSaveButton.setEnabled(false);
	}
	public void resetTitle()
	{
		setText(DEFAULT_TEXT);
		editSaveButton.setEnabled(true);
	}

}
