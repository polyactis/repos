package edu.nordborglab.client;

import java.util.Set;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONValue;
import com.google.gwt.json.client.JSONString;

import com.google.gwt.user.client.ui.DecoratedPopupPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.visualization.client.DataTable;


public class GWASOneResult extends CustomVerticalPanel{
	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	
	private int analysisMethodID;
	private String analysisMethodDescription;
	private String GWABaseURL;
	private String SNPBaseURL;
	
	private HTML statusReport;
	private String[] colors = {"blue", "green", "red", "cyan", "purple"};
	
	private DataTable dataTable;
	private DecoratedPopupPanel popupLink;
	
	public GWASOneResult(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, DecoratedPopupPanel popupLink, 
			int analysisMethodID, String analysisMethodDescription, String GWABaseURL, String SNPBaseURL)
	{
		super(constants, jsonErrorDialog, constants.GWASOneResultHelpID());
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		this.analysisMethodID = analysisMethodID;
		this.analysisMethodDescription = analysisMethodDescription;
		this.GWABaseURL = GWABaseURL;
		this.SNPBaseURL = SNPBaseURL;
		
		this.popupLink = popupLink;
		
		statusReport = new HTML(this.constants.LoadingText());
		this.add(statusReport);
		
		this.add(new HTML(this.analysisMethodDescription));
		
		loadGWA();
	}
	
	private class LoadGWAResponseHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONValue jsonValue = JSONParser.parse(response.getText());
				JSONObject serverData = jsonValue.isObject();	//jsonValue.isArray();
				JSONObject chr2data = serverData.get("chr2data").isObject();
				JSONObject chr2length = serverData.get("chr2length").isObject();
				double max_value = serverData.get("max_value").isNumber().doubleValue();
				int max_length = (int) serverData.get("max_length").isNumber().doubleValue();
				Set<String> keys = chr2data.keySet();
				int i =0;
				for (String chromosome : keys) 
				//for (i=0; i<keys.size(); i++)
				{
					//String chromosome = String.valueOf(i);
					//JSONString data = chr2data.get(chromosome).isString();
					String data = chr2data.get(chromosome).toString();
					
					//jsonErrorDialog.displayJSONObject(chr2data.get(chromosome));
					int chrLength = (int) chr2length.get(chromosome).isNumber().doubleValue();
					//jsonErrorDialog.displayRequestError("chromosome "+ chromosome + "length: " + chrLength);
					String color = colors[i%colors.length];
					AssociationScatterChart associationChart = new AssociationScatterChart(constants, jsonErrorDialog, chromosome,
							color, chrLength, max_length, max_value, SNPBaseURL);
					add(associationChart);
					dataTable = Common.asDataTable(data.substring(1, data.length()-1));	//2009-4-25 data has extra " on both ends
					associationChart.draw_gwas(dataTable);
					i += 1;
				}
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	private void loadGWA()
	{
		setIntoWaitState();
		String url = URL.encode(this.GWABaseURL + "&analysis_method_id="+this.analysisMethodID);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new LoadGWAResponseHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
	}
	
	public void setIntoWaitState()
	{
		statusReport.setVisible(true);
	}
	public void resetTitle()
	{
		statusReport.setVisible(false);
	}
}
