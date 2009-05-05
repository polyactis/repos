package edu.nordborglab.client;

import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.DataView;
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.visualization.client.visualizations.Table;

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



public class Common {
	public static native DataTable asDataTable(String json) /*-{
		dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
		return dataTable;
	}-*/;
	
	public static native DataTable jsonString2DataTable(JSONString json) /*-{
	dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
	return dataTable;
	}-*/;
	
	/*
	 * 2009-4-25 ``double jsonify`` means the data structure from the server gets jsonified twice on the server end. 
	 * JSONParser.parse() would return it with double quote (") on both ends, need to remove.
	 */
	public static native DataTable doubleJsonify2DataTable(String json) /*-{
	dataTable = new $wnd.google.visualization.DataTable(eval("("+json.substring(1, json.length-1)+")"), 0.5); 
	return dataTable;
	}-*/;
	
	private static class FillSelectBoxJSONResponseHandler implements RequestCallback {
		private ListBox selectBox;
		private DisplayJSONObject jsonErrorDialog;
		FillSelectBoxJSONResponseHandler(ListBox selectBox, DisplayJSONObject jsonErrorDialog)
		{
			this.selectBox = selectBox;
			this.jsonErrorDialog = jsonErrorDialog;
		}
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONObject options = JSONParser.parse(responseText).isObject();
				JSONArray optionArray = options.get("options").isArray();
				for (int i = 0; i < optionArray.size(); i++) {
					JSONObject idValueDict= optionArray.get(i).isObject();
					String item = idValueDict.get("id").isString().stringValue();
					//Integer id_int = id.intValue();
					
					Double value = idValueDict.get("value").isNumber().doubleValue();
					Integer valueInt = value.intValue();
					this.selectBox.addItem(item, valueInt.toString());
				}

			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
		}
	}

	public static void fillSelectBox(String url, ListBox selectBox, DisplayJSONObject jsonErrorDialog)
	{
		//setText(DIALOG_WAITING_TEXT);
		String _url = URL.encode(url);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, _url);
		try {
			requestBuilder.sendRequest(null, new FillSelectBoxJSONResponseHandler(selectBox, jsonErrorDialog));
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
		}
	}
	
	public static void onSelectChange(String url, ListBox selectBox, ListBox targetListBox, DisplayJSONObject jsonErrorDialog)
	{
		//setText(DIALOG_WAITING_TEXT);
		String _url = URL.encode(url);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, _url);
		try {
			requestBuilder.sendRequest(null, new FillSelectBoxJSONResponseHandler(targetListBox, jsonErrorDialog));
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
		}
	}
	
	private static class MotionChartDataResponseHandler implements RequestCallback {
		private MotionChart motionChart;
		private DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		MotionChartDataResponseHandler(MotionChart motionChart, DisplayJSONObject jsonErrorDialog, Widget statusReport, AbstractDataTable dataTable)
		{
			this.motionChart = motionChart;
			this.jsonErrorDialog = jsonErrorDialog;
			this.statusReport = statusReport;
			this.dataTable  = dataTable;
		}
		
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			statusReport.setVisible(false);
		}
	
		public void onResponseReceived(Request request, Response response) {
			statusReport.setVisible(true);
			String responseText = response.getText();
			try {
				dataTable = Common.asDataTable(responseText);
				MotionChart.Options options = MotionChart.Options.create();
				//com.google.gwt.core.client.JavaScriptObject value;
				//options.setOption("state", "{'colorOption':4,'sizeOption':'_UNISIZE'};");
				//String stateStr = new "{'time':'notime','iconType':'BUBBLE','xZoomedDataMin':null,'yZoomedDataMax':null,'xZoomedIn':false,'iconKeySettings':[],'showTrails':true,'xAxisOption':'2','colorOption':'4','yAxisOption':'3','playDuration':15,'xZoomedDataMax':null,'orderedByX':false,'duration':{'multiplier':1,'timeUnit':'none'},'xLambda':1,'orderedByY':false,'sizeOption':'_UNISIZE','yZoomedDataMin':null,'nonSelectedAlpha':0.4,'stateVersion':3,'dimensions':{'iconDimensions':['dim0']},'yLambda':1,'yZoomedIn':false};";
				//2009-4-22 stateStr written same as in javascript won't work because it's encoded into the following crap in client's browser and GWT doesnt' do the encoding for you.
				//options.setOption("state", "%7B%22time%22%3A%22notime%22%2C%22iconType%22%3A%22BUBBLE%22%2C%22xZoomedDataMin%22%3Anull%2C%22yZoomedDataMax%22%3Anull%2C%22xZoomedIn%22%3Afalse%2C%22iconKeySettings%22%3A%5B%5D%2C%22showTrails%22%3Atrue%2C%22xAxisOption%22%3A%222%22%2C%22colorOption%22%3A%224%22%2C%22yAxisOption%22%3A%223%22%2C%22playDuration%22%3A15%2C%22xZoomedDataMax%22%3Anull%2C%22orderedByX%22%3Afalse%2C%22duration%22%3A%7B%22multiplier%22%3A1%2C%22timeUnit%22%3A%22none%22%7D%2C%22xLambda%22%3A1%2C%22orderedByY%22%3Afalse%2C%22sizeOption%22%3A%22_UNISIZE%22%2C%22yZoomedDataMin%22%3Anull%2C%22nonSelectedAlpha%22%3A0.4%2C%22stateVersion%22%3A3%2C%22dimensions%22%3A%7B%22iconDimensions%22%3A%5B%22dim0%22%5D%7D%2C%22yLambda%22%3A1%2C%22yZoomedIn%22%3Afalse%7D%3B");
				options.setHeight(600);
				options.setWidth(1000);
				motionChart.draw(dataTable, options);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			statusReport.setVisible(false);
		}
	}
	
	public static void fillInMotionChart(String url, MotionChart motionChart, DisplayJSONObject jsonErrorDialog, Widget statusReport,
			AbstractDataTable dataTable)
	{
		statusReport.setVisible(true);
		String _url = URL.encode(url);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, _url);
		try {
			requestBuilder.sendRequest(null, new MotionChartDataResponseHandler(motionChart, jsonErrorDialog, statusReport, dataTable));
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
		}
	}
	
	
	private static class TableResponseTextHandler implements RequestCallback {
		private Table table;
		private DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		TableResponseTextHandler(Table table, DisplayJSONObject jsonErrorDialog, Widget statusReport, AbstractDataTable dataTable)
		{
			this.table = table;
			this.jsonErrorDialog = jsonErrorDialog;
			this.statusReport = statusReport;
			this.dataTable  = dataTable;
		}
		
		
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			statusReport.setVisible(false);
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {

				//JSONValue jsonValue = JSONParser.parse(responseText);
				//displayJSONObject(jsonValue);
				dataTable = Common.asDataTable(responseText);	//DataTable.create();//new google.visualization.DataTable(eval("("+response+")"), 0.5)
				Table.Options options = Table.Options.create();
				options.setShowRowNumber(true);
				options.setAllowHtml(true);
				table.draw(dataTable, options);
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			statusReport.setVisible(false);
		}
	}


	public static void fillInTable(String url, Table table, DisplayJSONObject jsonErrorDialog, Widget statusReport,
			AbstractDataTable dataTable)
	{
		statusReport.setVisible(true);
		String _url = URL.encode(url);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, _url);
		try {
			requestBuilder.sendRequest(null, new TableResponseTextHandler(table, jsonErrorDialog, statusReport, dataTable));
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			statusReport.setVisible(false);
		}
	}
	
}
