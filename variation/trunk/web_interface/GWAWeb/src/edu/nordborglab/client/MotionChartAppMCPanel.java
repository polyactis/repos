/*
 * 2009-8-31 app to view a user-uploaded data matrix in motion chart
 * 
 */
package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.FileUpload;
import com.google.gwt.user.client.ui.FormHandler;
import com.google.gwt.user.client.ui.FormPanel;
import com.google.gwt.user.client.ui.FormSubmitEvent;
import com.google.gwt.user.client.ui.FormSubmitCompleteEvent;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.Window;

import com.google.gwt.json.client.JSONException;


import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.DataView;
import com.google.gwt.visualization.client.visualizations.ColumnChart;
import com.google.gwt.visualization.client.visualizations.MotionChart;

public class MotionChartAppMCPanel extends Composite {
	
	final private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private FileUpload upload;
	private Button submitButton;
	private DisclosurePanel motionChartDisclosurePanel;
	private MotionChart motionChart;
	
	private CustomVerticalPanel vPanel;
	
	private String geneNameAutoCompleteURL;
	
	private static final String BUTTON_DEFAULT_TEXT = "Upload";
	private static final String BUTTON_WAITING_TEXT = "Waiting...";
	
	public MotionChartAppMCPanel(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, String formActionURL)
	{
		//super(constants);
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		vPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.MotionChartAppMCPanelHelpID());
		
		final FormPanel form = new FormPanel();
		form.setAction(formActionURL);
		
		// Because we're going to add a FileUpload widget, we'll need to set the
		// form to use the POST method, and multipart MIME encoding.
		form.setEncoding(FormPanel.ENCODING_MULTIPART);
		form.setMethod(FormPanel.METHOD_POST);

		// Create a panel to hold all of the form widgets.
		VerticalPanel panel = new VerticalPanel();
		form.setWidget(panel);
		
		// Create a FileUpload widget.
		upload = new FileUpload();
		upload.setName("dataMatrixFile");
		
		// create a 'submit' button
		submitButton = new Button(BUTTON_DEFAULT_TEXT);
		submitButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				form.submit();
			}
		});
		
		// a horizontal panel to hold FileUpload and Button 
		HorizontalPanel hPanel = new HorizontalPanel();
		panel.add(hPanel);
		
		hPanel.setSpacing(20);
		hPanel.add(upload);
		hPanel.add(submitButton);

		motionChartDisclosurePanel = new DisclosurePanel("MotionChart (Click to refresh)");
		motionChart = new MotionChart();
		//motionChart.setSize("800px", "1000px");	 	//doesn't work. have to be set in the options before draw()
		motionChartDisclosurePanel.setContent(motionChart);
		motionChartDisclosurePanel.setAnimationEnabled(true);
		motionChartDisclosurePanel.setOpen(true);
		panel.add(motionChartDisclosurePanel);		
		
		// Add an event handler to the form.
		form.addFormHandler(new FormHandler() {
			public void onSubmit(FormSubmitEvent event) {
				// This event is fired just before the form is submitted. We can
				// take
				// this opportunity to perform validation.
				if (upload.getFilename().length() == 0) {
					Window.alert("No File is selected.");
					event.setCancelled(true);
				}
			}

			public void onSubmitComplete(FormSubmitCompleteEvent event) {
				// When the form submission is successfully completed, this
				// event is
				// fired. Assuming the service returned a response of type
				// text/html,
				// we can get the result text here (see the FormPanel
				// documentation for
				// further explanation).
				// Window.alert(event.getResults());
				String responseText = event.getResults();
				if (responseText.startsWith("Error"))	// 2009-9-1 server returns an error message.
				{
					DisplayJSONObject jsonErrorDialog;
					jsonErrorDialog = new DisplayJSONObject("Error Dialog");
					jsonErrorDialog.displayParseError(responseText);
				}
				else	// 2009-9-1 presumably a good json data structure
				{
					try {
						motionChart.setVisible(false);
						DataTable motionData = Common.asDataTable(responseText);
						MotionChart.Options options = MotionChart.Options.create();
						//com.google.gwt.core.client.JavaScriptObject value;
						//2009-4-22 stateStr written same as in javascript won't work because it's encoded into the following crap in client's browser and GWT doesnt' do the encoding for you.
						options.setHeight(600);
						options.setWidth(1000);
						motionChart.draw(motionData, options);
						motionChart.setVisible(true);
					}
					catch (JSONException e) {
						DisplayJSONObject jsonErrorDialog;
						jsonErrorDialog = new DisplayJSONObject("Error Dialog");
						jsonErrorDialog.displayParseError(responseText);
					}					
				}
			}
		});

		vPanel.add(form);
		// All composites must call initWidget() in their constructors.
		initWidget(vPanel);
	}
}
