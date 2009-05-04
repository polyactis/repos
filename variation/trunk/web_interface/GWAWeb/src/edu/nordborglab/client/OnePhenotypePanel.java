package edu.nordborglab.client;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.core.client.JsArrayString;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.DataView;
import com.google.gwt.visualization.client.visualizations.ColumnChart;
import com.google.gwt.visualization.client.visualizations.MotionChart;


public class OnePhenotypePanel extends CustomVerticalPanel{
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private ColumnChart phenotypeHistChart;
	private DisclosurePanel phenotypeHistPanel;
	private HorizontalPanel hPanel;
	private DisclosurePanel hPanelDisclosurePanel;
	private MotionChart motionChart;
	private DisclosurePanel motionChartDisclosurePanel;

	private int callMethodID;
	private int phenotypeMethodID;
	private String phenotypeMethodShortName;
	private String phenotypeMethodDescription;

	private String callInfoURL;
	private String phenotypeHistImageURL;
	private String callPhenotypeQQImageURL;
	private String phenotypeHistogramDataURL;
	private JsArray staticPlotNameArr;
	
	private HTML report;
	private Tree infoTree;
	
	private DialogBox imageBox;
	
	OnePhenotypePanel(AccessionConstants constants, DisplayJSONObject jsonErrorDialog)
	{
		super(constants, jsonErrorDialog, constants.OnePhenotypeHelpID());
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		callMethodID = getCallMethodID();
		phenotypeMethodID = getPhenotypeMethodID();
		phenotypeMethodShortName = getphenotypeMethodShortName();
		phenotypeMethodDescription = getphenotypeMethodDescription();

		Window.setTitle(phenotypeMethodID + " " + phenotypeMethodShortName);
		infoTree = new Tree();
		infoTree.setVisible(true);
		TreeItem treeItem = infoTree.addItem(new HTML("<b>Phenotype ID:</b> "+phenotypeMethodID));
		treeItem.setState(true);
		TreeItem nameTreeItem = infoTree.addItem(new HTML("<b>Name:</b>" + phenotypeMethodShortName));
		//nameTreeItem.addItem(phenotypeMethodShortName);
		nameTreeItem.setState(true);
		TreeItem descTreeItem = infoTree.addItem(new HTML("<b>Description:</b>"));
		descTreeItem.addItem(phenotypeMethodDescription);
		descTreeItem.setState(true);
		
		add(infoTree);
		report = new HTML(constants.LoadingText());
		add(report);

		callInfoURL = getCallInfoURL();
		phenotypeHistImageURL = getPhenotypeHistImageURL();
		callPhenotypeQQImageURL = getCallPhenotypeQQImageURL();
		phenotypeHistogramDataURL = getPhenotypeHistogramDataURL();
		
		staticPlotNameArr = getstaticPlotNameArr();
		
		phenotypeHistChart = new ColumnChart();
		phenotypeHistChart.setSize("800px", "400px");
		phenotypeHistPanel = new DisclosurePanel(constants.phenotypeHistPanelLabel());
		phenotypeHistPanel.setAnimationEnabled(true);
		phenotypeHistPanel.setContent(phenotypeHistChart);
		phenotypeHistPanel.setOpen(true);
		add(phenotypeHistPanel);
		
		hPanel = new HorizontalPanel();
		hPanelDisclosurePanel = new DisclosurePanel(constants.staticHistQQPlotPanelLabel());
		hPanelDisclosurePanel.setAnimationEnabled(true);
		hPanelDisclosurePanel.setContent(hPanel);
		//hPanelDisclosurePanel.setOpen(true);
		add(hPanelDisclosurePanel);

		motionChart = new MotionChart();
		motionChart.setSize("1000px", "800px");	 	//doesn't work. have to be set in the options before draw()
		motionChartDisclosurePanel = new DisclosurePanel(constants.motionChartDisclosurePanelLabel());
		motionChartDisclosurePanel.setContent(motionChart);
		motionChartDisclosurePanel.setAnimationEnabled(true);
		//motionChartDisclosurePanel.setOpen(true);
		add(motionChartDisclosurePanel);
		
		imageBox = new DialogBox(true, false);
		imageBox.setText("click outside to close");
		fillIn();
	}
	
	public native int getCallMethodID()/*-{
	return $wnd.call_method_id;
	}-*/;
	
	public native int getPhenotypeMethodID()/*-{
		return $wnd.phenotype_method_id;
	}-*/;
	
	public native String getphenotypeMethodShortName()/*-{
	return $wnd.phenotype_method_short_name;
	}-*/;
	
	public native String getphenotypeMethodDescription()/*-{
		return $wnd.phenotype_method_description;
	}-*/;
	
	public native String getCallInfoURL()/*-{
		return $wnd.callInfoURL;
	}-*/;
	
	public native String getPhenotypeHistImageURL()/*-{
		return $wnd.phenotypeHistImageURL;
	}-*/;	
	
	public native String getCallPhenotypeQQImageURL()/*-{
	return $wnd.callPhenotypeQQImageURL;
	}-*/;
	
	public native String getPhenotypeHistogramDataURL()/*-{
		return $wnd.phenotypeHistogramDataURL;
	}-*/;
	
	public native JsArray getstaticPlotNameArr()/*-{
		return $wnd.static_plot_name_arr;
	}-*/;
	
	private class PhenotypeHistogramDataJSONResponseHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}
	
		public void onResponseReceived(Request request, Response response) {
			report.setVisible(true);
			String responseText = response.getText();
			try {
				DataTable dataTable = Common.asDataTable(responseText);
				ColumnChart.Options options = ColumnChart.Options.create();
				options.setTitle("Phenotype Histogram");
				//options.setWidth(600);
				options.set3D(true);
				phenotypeHistChart.draw(dataTable, options);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	public int[] findColIndex(AbstractDataTable dataTable)
	{
		int label_idx = 2;	//will find it according to the column id automatically later, which doesn't work in GWT shell.
		int date_idx = 0;
		int lon_idx = 4;
		int lat_idx = 3;
		int phenotype_idx = 8;
		int name_idx = 5;
		int pc1_idx = 6;
		int pc2_idx = 7;
		int no_of_cols = dataTable.getNumberOfColumns();
		for (int i =0; i<no_of_cols; i++)
		{
			String col_id = dataTable.getColumnId(i);
			if (col_id.equals("label"))
				label_idx = i;
			else if (col_id.equals("date"))
				date_idx = i;
			else if (col_id.equals("lon"))
				lon_idx=i;
			else if (col_id.equals("lat"))
				lat_idx=i;
			else if (col_id.equals("phenotype"))
				phenotype_idx=i;
			else if (col_id.equals("name"))
				name_idx=i;
			else if (col_id.equals("pc1"))
				pc1_idx=i;
			else if (col_id.equals("pc2"))
				pc2_idx=i;
		}
		int [] columnIndices = {label_idx, date_idx, lon_idx, lat_idx, phenotype_idx, name_idx, pc1_idx, pc2_idx};
		return columnIndices;
	}
	
	private class MotionChartDataResponseHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetTitle();
		}
	
		public void onResponseReceived(Request request, Response response) {
			report.setVisible(true);
			String responseText = response.getText();
			try {
				DataTable callInfoData = Common.asDataTable(responseText);
				DataView motionView = DataView.create(callInfoData);
				int[] columnIndices = findColIndex(callInfoData); // {2,0,4,3,8,5,6,7};
				motionView.setColumns(columnIndices);
				MotionChart.Options options = MotionChart.Options.create();
				//com.google.gwt.core.client.JavaScriptObject value;
				//options.setOption("state", "{'colorOption':4,'sizeOption':'_UNISIZE'};");
				//String stateStr = new "{'time':'notime','iconType':'BUBBLE','xZoomedDataMin':null,'yZoomedDataMax':null,'xZoomedIn':false,'iconKeySettings':[],'showTrails':true,'xAxisOption':'2','colorOption':'4','yAxisOption':'3','playDuration':15,'xZoomedDataMax':null,'orderedByX':false,'duration':{'multiplier':1,'timeUnit':'none'},'xLambda':1,'orderedByY':false,'sizeOption':'_UNISIZE','yZoomedDataMin':null,'nonSelectedAlpha':0.4,'stateVersion':3,'dimensions':{'iconDimensions':['dim0']},'yLambda':1,'yZoomedIn':false};";
				//2009-4-22 stateStr written same as in javascript won't work because it's encoded into the following crap in client's browser and GWT doesnt' do the encoding for you.
				options.setOption("state", "%7B%22time%22%3A%22notime%22%2C%22iconType%22%3A%22BUBBLE%22%2C%22xZoomedDataMin%22%3Anull%2C%22yZoomedDataMax%22%3Anull%2C%22xZoomedIn%22%3Afalse%2C%22iconKeySettings%22%3A%5B%5D%2C%22showTrails%22%3Atrue%2C%22xAxisOption%22%3A%222%22%2C%22colorOption%22%3A%224%22%2C%22yAxisOption%22%3A%223%22%2C%22playDuration%22%3A15%2C%22xZoomedDataMax%22%3Anull%2C%22orderedByX%22%3Afalse%2C%22duration%22%3A%7B%22multiplier%22%3A1%2C%22timeUnit%22%3A%22none%22%7D%2C%22xLambda%22%3A1%2C%22orderedByY%22%3Afalse%2C%22sizeOption%22%3A%22_UNISIZE%22%2C%22yZoomedDataMin%22%3Anull%2C%22nonSelectedAlpha%22%3A0.4%2C%22stateVersion%22%3A3%2C%22dimensions%22%3A%7B%22iconDimensions%22%3A%5B%22dim0%22%5D%7D%2C%22yLambda%22%3A1%2C%22yZoomedIn%22%3Afalse%7D%3B");
				options.setHeight(600);
				options.setWidth(1000);
				motionChart.draw(motionView, options);
			}
			catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetTitle();
		}
	}
	
	public void fillIn()
	{
		/*
		 * fill in the phenotype histogram (column chart), 4 static images, motion chart
		 */
		String url = URL.encode(phenotypeHistogramDataURL);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new PhenotypeHistogramDataJSONResponseHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
	
		url = URL.encode(callInfoURL);
		requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new MotionChartDataResponseHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetTitle();
		}
		
		for (int i=0; i<staticPlotNameArr.length(); i++)
		{
			JsArrayString onePlotNameArr = (JsArrayString) staticPlotNameArr.get(i);
			String img_type = onePlotNameArr.get(1);
			String fullImg = onePlotNameArr.get(2);
			String imageURL;
			final String fullImageURL;
			if (onePlotNameArr.get(0)=="getPhenotypeHistImage")
			{
				imageURL = phenotypeHistImageURL + "&img_type="+img_type;
				fullImageURL = phenotypeHistImageURL + "&img_type="+fullImg;
			}
			else
			{
				imageURL = callPhenotypeQQImageURL + "&img_type="+img_type;
				fullImageURL = callPhenotypeQQImageURL + "&img_type="+fullImg;
			}
			Image image = new Image();
			image.setUrl(imageURL);
			image.setTitle(constants.staticImageToolTip());
			image.addClickListener(new ClickListener(){
				public void onClick(Widget sender)
				{
					//imageBox.clear();	// not necessary because setWidth() is called.
					Image fullImage = new Image();
					fullImage.setUrl(fullImageURL);
					imageBox.setWidget(fullImage);
					imageBox.center();
					imageBox.show();
				}
			});
			hPanel.add(image);
		}
	}
	
	public void resetTitle()
	{
		report.setVisible(false);
		//Window.setTitle(TITLE_DEFAULT_TEXT);
	}
}
