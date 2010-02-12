package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ChangeListener;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Visualization;

import com.google.gwt.core.client.JsArray;
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
	
	/**
	 * 2009-6-29 used internally in fillSelectBox() & onSelectChange()
	 * @author crocea
	 *
	 */
	public static class FillSelectBoxJSONResponseHandler implements RequestCallback {
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
			this.selectBox.clear();
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
	
	/**
	 * 2009-6-30 API to fill a select list box with stuff fetched from one URL
	 * @param url
	 * @param selectBox
	 * @param jsonErrorDialog
	 */
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
	
	/**
	 * 2009-6-30 not sure why i wrote the function below or how to use it
	 * @param url
	 * @param selectBox
	 * @param targetListBox
	 * @param jsonErrorDialog
	 */
	public static void fillTargetListBoxOnSelectChange(String url, ListBox selectBox, ListBox targetListBox, DisplayJSONObject jsonErrorDialog)
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
	
	/**
	 * 2009-7-1 a generic change listener for one list box. upon its change, another list box is to be filled by content fetched from a baseURL
	 * 	the value of the 1st list box must be String
	 * @author crocea
	 *
	 */
	
	public static class FillTargetListOnListChangeListener implements ChangeListener
	{
		String _url;
		String argumentName;
		ListBox listBox;
		DisplayJSONObject jsonErrorDialog;
		FillTargetListOnListChangeListener(String baseURL, String argumentName, ListBox listBox, DisplayJSONObject jsonErrorDialog)
		{
			_url = baseURL;
			this.argumentName = argumentName;
			this.listBox = listBox;
			this.jsonErrorDialog = jsonErrorDialog;
		}
		public void onChange(Widget sender) {
			ListBox senderListBox = (ListBox) sender;
			String selectedValue = senderListBox.getValue(senderListBox.getSelectedIndex());
			String url = _url + "?" + argumentName + "=" + selectedValue;
			Common.fillSelectBox(url, listBox, jsonErrorDialog);
		}
	}
	
	/**
	 * 2009-7-1 a generic change listener for one list box. upon its change, 
	 * 	another list box is to be filled by content fetched from a baseURL + values from two list boxes.
	 * 	the value of the 1st list box must be String
	 * @author crocea
	 *
	 */
	
	public static class FillTargetListBasedOnTwoListsChange implements ChangeListener
	{
		String _url;
		ListBox listBox1;
		String argumentName1;
		ListBox listBox2;
		String argumentName2;
		ListBox targetListBox;
		DisplayJSONObject jsonErrorDialog;
		FillTargetListBasedOnTwoListsChange(String baseURL, ListBox listBox1, String argumentName1, ListBox listBox2, String argumentName2, 
				ListBox targetListBox, DisplayJSONObject jsonErrorDialog)
		{
			_url = baseURL;
			this.listBox1 = listBox1;
			this.argumentName1 = argumentName1;
			this.listBox2 = listBox2;
			this.argumentName2 = argumentName2;
			this.targetListBox = targetListBox;
			this.jsonErrorDialog = jsonErrorDialog;
		}
		public void onChange(Widget sender) {
			String urlArguments = "";
			String selectedValue1 = listBox1.getValue(listBox1.getSelectedIndex());
			urlArguments += argumentName1 +"=" + selectedValue1;
			String selectedValue2 = listBox2.getValue(listBox2.getSelectedIndex());
			urlArguments += "&" + argumentName2 +"=" + selectedValue2;
			String url = _url + "?" + urlArguments;
			Common.fillSelectBox(url, targetListBox, jsonErrorDialog);
		}
	}
	
	/**
	 * 2009-7-1 MotionChartDataResponseHandler used in fillInMotionChart()
	 * @author crocea
	 *
	 */
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
	
	/**
	 * 2009-7-1 a function to fill in motion chart with data fetched from a url
	 * @param url
	 * @param motionChart
	 * @param jsonErrorDialog
	 * @param statusReport
	 * @param dataTable
	 */
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
	
	/**
	 * 2009-7-1 class used in fillInTable()
	 * @author crocea
	 *
	 */
	private static class TableResponseTextHandler implements RequestCallback {
		private AbstractVisualization<Table.Options> abstractVTable;
		private Visualization<Table.Options> vTable;
		private DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		private int isVisualizationAbstract = 0;
		TableResponseTextHandler(AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, Widget statusReport, AbstractDataTable dataTable)
		{
			this.abstractVTable = table;
			this.isVisualizationAbstract = 1;
			this.jsonErrorDialog = jsonErrorDialog;
			this.statusReport = statusReport;
			this.dataTable  = dataTable;
		}
		
		TableResponseTextHandler(Visualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, Widget statusReport, AbstractDataTable dataTable)
		{
			this.vTable = table;
			this.isVisualizationAbstract = 0;
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
				if (this.isVisualizationAbstract==1)
					this.abstractVTable.draw(dataTable, options);
				else
					this.vTable.draw(dataTable, options);
			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			statusReport.setVisible(false);
		}
	}

	/**
	 * 2009-7-1 a function to fill in a google abstract visualization (custom class WrapVisualizationTable) with data fetched from a url
	 * @param url
	 * @param table
	 * @param jsonErrorDialog
	 * @param statusReport
	 * @param dataTable
	 */
	public static void fillInTable(String url, AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, Widget statusReport,
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
	
	/**
	 * 2009-7-1 a function to fill in google visualization table with data fetched from a url
	 * @param url
	 * @param table
	 * @param jsonErrorDialog
	 * @param statusReport
	 * @param dataTable
	 */
	public static void fillInTable(String url, Visualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, Widget statusReport,
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
	
	/**
	 * 2009-7-1 refactored out of HaplotypeSingleView
	 * @param _url
	 * @param argumentName
	 * @param selectBox
	 * @param defaultValue
	 * @return
	 */
	public static String getURLFromSelectBox(String _url, String argumentName, ListBox selectBox, DisplayJSONObject jsonErrorDialog, 
			String defaultValue)
	{
		String selectedValue;
		if (selectBox.getItemCount()==0)
			if (defaultValue.isEmpty() || defaultValue.equals(null))
			{
				jsonErrorDialog.displayError("Form Error", argumentName + " hasn't been selected.");
				return "";
			}
			else
				selectedValue = defaultValue;
		else
			selectedValue = selectBox.getValue(selectBox.getSelectedIndex());
		
		if ((!argumentName.equals("list_type_id") && selectedValue.equals("0")) || selectedValue.equals("-1"))	// list_type_id is allowed to be 0.
		{
			if (defaultValue.isEmpty() || defaultValue.equals(null))
			{
				jsonErrorDialog.displayError("Form Error", argumentName + ", "+ selectedValue + ", is invalid.");
				return "";
			}
			else
				selectedValue = defaultValue;
		}

		if (_url.isEmpty())
		{
			_url += argumentName+"="+selectedValue;
		}
		else
			_url += "&"+argumentName+"="+selectedValue;
		return _url;

	}
	
	/**
	 * 2009-7-1 a common change listener that fills the WrapVisualizationTable based on the selected value of a list Box.
	 * 	the value must be String. baseURL contains everything except the selected value.
	 * 
	 * - user shall customize constructURLArgument()
	 * 
	 * @author crocea
	 *
	 */
	public static class FillWrapVTableOnListChangeListener implements ChangeListener
	{
		String _url;
		private AbstractVisualization<Table.Options> abstractVTable;
		DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		FillWrapVTableOnListChangeListener(String baseURL, AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, 
					Widget statusReport, AbstractDataTable dataTable)
		{
			_url = baseURL;
			this.abstractVTable = table;
			this.jsonErrorDialog = jsonErrorDialog;
			this.statusReport = statusReport;
			this.dataTable  = dataTable;
		}
		/**
		 * to be customized
		 * @return
		 */
		public String constructURLArgument()
		{
			return "";
		}
		public void onChange(Widget sender) {
			String urlArguments = constructURLArgument();
			if (!urlArguments.isEmpty())
			{
				String url = _url + "?"+ urlArguments;
				Common.fillInTable(url, abstractVTable, jsonErrorDialog, statusReport, dataTable);
			}
		}
	}
	
	/**
	 * 2009-7-1 a common class to fill the WrapVisualizationTable when the button is clicked.
	 * 	the value must be String. baseURL contains everything.
	 * 
	 * - user shall customize constructURLArgument()
	 * 
	 * @author crocea
	 *
	 */
	public static class FillWrapVTableOnClickListener implements ClickListener
	{
		String _url;
		private AbstractVisualization<Table.Options> abstractVTable;
		DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		FillWrapVTableOnClickListener(String baseURL, AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, 
					Widget statusReport, AbstractDataTable dataTable)
		{
			_url = baseURL;
			this.abstractVTable = table;
			this.jsonErrorDialog = jsonErrorDialog;
			this.statusReport = statusReport;
			this.dataTable  = dataTable;
		}
		/**
		 * to be customized
		 * @return
		 */
		public String constructURLArgument()
		{
			return "";
		}
		public void onClick(Widget sender) {
			Button button = (Button) sender;
			button.setEnabled(false);
			String urlArguments = constructURLArgument();
			if (!urlArguments.isEmpty())
			{
				String url = _url + "?"+ urlArguments;
				Common.fillInTable(url, abstractVTable, jsonErrorDialog, statusReport, dataTable);
			}
			button.setEnabled(true);
		}
	}
	
	public static void setCertainValuedItemSelectedInListBox(ListBox selectBox, String value)
	{
		for (int i =0; i<selectBox.getItemCount(); i++)
		{
			String itemValue = selectBox.getValue(i);
			if (itemValue==value)
			{
				selectBox.setItemSelected(i, true);
				break;
			}
			//newListBox.addItem(oldListBox.getItemText(i), oldListBox.getValue(i));
		}
	}
	public static native String getCallMethodID()/*-{ return $wnd.call_method_id;}-*/;
	public static native String getPhenotypeMethodID()/*-{ return $wnd.phenotype_method_id; }-*/;
	public static native String getAnalysisMethodID()/*-{ return $wnd.analysis_method_id; }-*/;
	public static native String getSNPGeneAssociationTypeID()/*-{ return $wnd.snpGeneAssociationTypeID; }-*/;
	public static native String getPhenotypeMethodShortName()/*-{ return $wnd.phenotype_method_short_name;}-*/;
	public static native String getPhenotypeMethodDescription()/*-{ return $wnd.phenotype_method_description; }-*/;
	public static native String getCallInfoURL()/*-{ return $wnd.callInfoURL;}-*/;
	public static native String getPhenotypeHistImageURL()/*-{ return $wnd.phenotypeHistImageURL;}-*/;	
	public static native String getCallPhenotypeQQImageURL()/*-{ return $wnd.callPhenotypeQQImageURL;}-*/;
	public static native String getPhenotypeHistogramDataURL()/*-{ return $wnd.phenotypeHistogramDataURL;}-*/;
	public static native JsArray getStaticPlotNameArr() /*-{ return $wnd.static_plot_name_arr; }-*/;
	public static native String getGWABaseURL()/*-{ return $wnd.GWABaseURL}-*/;
	public static native String getSNPBaseURL()/*-{ return $wnd.SNPBaseURL}-*/;
	public static native String getAnalysisMethodLsURL()/*-{ return $wnd.getAnalysisMethodLsURL}-*/;
	public static native String getCallMethodLsURL()/*-{ return $wnd.callMethodLsURL; }-*/;
	public static native String getGeneListLsURL()/*-{ return $wnd.geneListLsURL; }-*/;
	public static native String getCallMethodOnChangeURL()/*-{ return $wnd.callMethodOnChangeURL; }-*/;
	public static native String getHaplotypeImgURL()/*-{ return $wnd.haplotypeImgURL; }-*/;
	public static native String getSNPGeneAssociationTypeLsURL() /*-{ return $wnd.snpGeneAssociationTypeLsURL; }-*/;
	public static native String getSNPGeneAssociationOnChangeURL() /*-{ return $wnd.snpGeneAssociationOnChangeURL; }-*/;
	public static native String getPhenotypeMethodOnChangeURL() /*-{ return $wnd.phenotypeMethodOnChangeURL; }-*/;
	public static native String getFetchResultsGeneURL() /*-{ return $wnd.fetchResultsGeneURL; }-*/;
	public static native String getCandidateGeneListURL() /*-{ return $wnd.candidateGeneListURL; }-*/;
	public static native String getGeneListTypeID() /*-{ return $wnd.geneListTypeID; }-*/;
	public static native String getFetchGeneListURL() /*-{ return $wnd.fetchGeneListURL; }-*/;
}