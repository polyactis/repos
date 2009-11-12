package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.FlexTable;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.DecoratedPopupPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.Hyperlink;
import com.google.gwt.user.client.ui.HasHorizontalAlignment.HorizontalAlignmentConstant;

import com.google.gwt.user.client.Event;
import com.google.gwt.user.client.Window;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;



import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Table.Options.Policy;

import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.Query;
import com.google.gwt.visualization.client.QueryResponse;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.AbstractVisualization.VisualizationFactory;
import com.google.gwt.visualization.client.Query.Callback;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;

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


public class PhenotypeTable extends CustomVerticalPanel implements CustomClickListener{

	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	private String phenotypeCategoryID;
	private String phenotypeCategoryName;
	private String fetchPhenotypeTableDataURL;
	private int callMethodID;
	private String fetchGWAURL;
	private String displayResultsGeneURL;
	private String fetchPhenotypeURL;
	
	private HTML statusReport;
	private CustomVisualizationTable accessionTable;

	private DataTable dataTable;
	private DecoratedPopupPanel popupLink;

	private int phenotype_id_idx = -1;	//will find it according to the column id automatically later, which doesn't work in GWT shell.
	private int phenotype_name_idx = -1;
	private int association_results_idx = -1;
	
	public PhenotypeTable(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, DecoratedPopupPanel popupLink, String phenotypeCategoryID, 
			String phenotypeCategoryName, String fetchPhenotypeTableDataURL, int callMethodID, String fetchGWAURL, String displayResultsGeneURL,
			String fetchPhenotypeURL)
	{
		super(constants, jsonErrorDialog, constants.PhenotypeTableHelpID());
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		this.phenotypeCategoryID = phenotypeCategoryID;
		this.phenotypeCategoryName = phenotypeCategoryName;
		this.fetchPhenotypeTableDataURL = fetchPhenotypeTableDataURL;
		this.callMethodID = callMethodID;
		this.fetchGWAURL = fetchGWAURL;
		this.displayResultsGeneURL = displayResultsGeneURL;
		this.fetchPhenotypeURL = fetchPhenotypeURL;
		
		this.popupLink = popupLink;

		statusReport = new HTML("Waiting ...");
		this.add(statusReport);

		accessionTable = new CustomVisualizationTable();
		accessionTable.addClickListener(this);
		this.add(accessionTable);
		populateData();
	}
	public void findColIndex(AbstractDataTable dataTable)
	{
		int no_of_cols = dataTable.getNumberOfColumns();

		for (int i =0; i<no_of_cols; i++)
		{
			String col_id = dataTable.getColumnId(i);
			if (col_id.equals("id"))
				phenotype_id_idx = i;
			else if (col_id.equals("short_name"))
				phenotype_name_idx = i;
			else if (col_id.equals("association_results"))
				association_results_idx=i;
		}
	}

	private class JSONResponseTextHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {

				//JSONValue jsonValue = JSONParser.parse(responseText);
				//displayJSONObject(jsonValue);
				dataTable = Common.asDataTable(responseText);	//DataTable.create();//new google.visualization.DataTable(eval("("+response+")"), 0.5)
				findColIndex(dataTable);
				Table.Options options = Table.Options.create();
				options.setShowRowNumber(true);
				options.setAllowHtml(true);
				accessionTable.draw(dataTable, options);
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}


	public void populateData()
	{
		statusReport.setVisible(true);
		String url = URL.encode(fetchPhenotypeTableDataURL+"?biology_category_id="+phenotypeCategoryID + "&call_method_id="+callMethodID);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new JSONResponseTextHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
	}
	public void resetTitle()
	{
		statusReport.setVisible(false);
		statusReport.setText("");
		//DOM.getElementById("title").setInnerText(TITLE_DEFAULT_TEXT);
	}

	public void onClick(Widget sender, Event evt) {

		JsArray<Selection> selectionLs = accessionTable.getSelections();
		for (int i=0; i< selectionLs.length(); i++ )
		{
			popupLink.clear();
			popupLink.hide();
			Selection s = selectionLs.get(i);
			String association_results = dataTable.getFormattedValue(s.getRow(), association_results_idx);

			VerticalPanel vpanel = new VerticalPanel();
			FlexTable layout = new FlexTable();
			final String phenotype_name = dataTable.getFormattedValue(s.getRow(), phenotype_name_idx);
			String phenotype_id = dataTable.getFormattedValue(s.getRow(), phenotype_id_idx);
			HTML report = new HTML("Click below to open in new window.");
			//vpanel.add(report);
			layout.setWidget(0, 0, report);
			
			final String _fetchPhenotypeURL = fetchPhenotypeURL + "?phenotype_method_id=" + phenotype_id;
			HTML phenotypeLink = new HTML("<a href=" + _fetchPhenotypeURL + " target='_blank'>Information for phenotype: " + 
						phenotype_id + " " + phenotype_name + "</a>");
			layout.setWidget(1, 0, phenotypeLink);
			
			if (!association_results.isEmpty())
			{
				final String _fetchGWAURL = fetchGWAURL + "?call_method_id="+callMethodID+"&phenotype_method_id=" + phenotype_id;
				HTML gwaLink = new HTML("<a href=" + _fetchGWAURL + " target='_blank'>GWA plot for phenotype: " + phenotype_id + " " +
						phenotype_name + "</a>");
				/*
				 * gwaLink.addClickListener(new ClickListener() {
					public void onClick(Widget sender) {
						Window.open(_fetchGWAURL, phenotype_name, "");
					}
				});
				*/
				//vpanel.add(new HTML("<p></p>"));
				layout.setWidget(2, 0, gwaLink);
				
				final String _displayResultsGeneURL = displayResultsGeneURL + "?call_method_id="+callMethodID+"&phenotype_method_id=" + phenotype_id;
				HTML resultsGeneLink = new HTML("<a href=" + _displayResultsGeneURL + " target='_blank'>Genes within 20kb of Top SNPs for phenotype: " +
							phenotype_id + " " + phenotype_name + "</a>");
				/*
				resultsGeneLink.addClickListener(new ClickListener() {
					public void onClick(Widget sender) {
						Window.open(_displayResultsGeneURL, phenotype_name, "");
					}
				});
				*/
				//vpanel.add(new HTML("<p></p>"));
				layout.setWidget(3, 0, resultsGeneLink);
				
				/*
				if (association_results.contains("KW"))
				{
					final String kw_displayResultsGeneURL = displayResultsGeneURL + "?call_method_id="+callMethodID+"&phenotype_method_id=" + phenotype_id+"&analysis_method_id=1";
					HTML kw_resultsGeneLink = new HTML("<a>Genes associated with top 100 SNPs by KW for " + phenotype_id + " " + phenotype_name + "</a>");
					kw_resultsGeneLink.addClickListener(new ClickListener() {
						public void onClick(Widget sender) {
							Window.open(kw_displayResultsGeneURL, phenotype_name, "");
						}
					});
					layout.setWidget(2, 0, kw_resultsGeneLink);
				}
				if (association_results.contains("Emma"))
				{
					final String emma_displayResultsGeneURL = displayResultsGeneURL + "?call_method_id="+callMethodID+"&phenotype_method_id=" + phenotype_id+"&analysis_method_id=7";
					HTML emma_resultsGeneLink = new HTML("<a>Genes associated with top 100 SNPs by Emma for " + phenotype_id + " " + phenotype_name + "</a>");
					emma_resultsGeneLink.addClickListener(new ClickListener() {
						public void onClick(Widget sender) {
							Window.open(emma_displayResultsGeneURL, phenotype_name, "");
						}
					});
					layout.setWidget(3, 0, emma_resultsGeneLink);
				}
				*/
			}				
			vpanel.add(layout);
			Button closeButton = new Button("close",
					new ClickListener() {
				public void onClick(Widget sender) {
					popupLink.hide();
				}
			});
			vpanel.add(closeButton);
			vpanel.setCellHorizontalAlignment(closeButton, HasHorizontalAlignment.ALIGN_RIGHT);
			popupLink.add(vpanel);
			int left = evt.getClientX() + Window.getScrollLeft();
			int top = evt.getClientY() + Window.getScrollTop();
			popupLink.setPopupPosition(left, top);
			// Show the popup
			popupLink.show();

		}
	}

	public void onClick(Widget sender)
	{
	}
}