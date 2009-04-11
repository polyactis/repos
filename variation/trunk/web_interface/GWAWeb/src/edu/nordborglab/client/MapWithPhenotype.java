/*
 * 2009-4-1 yh: map widget, part of MapTableTree widget.
 *  1. wiggle the latitude, longitudes of strains with same lat,lon so they occupy different place
 *  2. map and the table in MapTableTree would respond to each other's select event
 *  3. color the accessions according to phenotypes selected or haplotype structure 
 */
package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.RadioButton;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.VerticalPanel;

import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.ChangeListener;

import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.ColumnChart;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.core.client.JsArray;
import com.google.gwt.visualization.client.AbstractDrawOptions;
import com.google.gwt.visualization.client.ArrayHelper;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;
import com.google.gwt.visualization.client.visualizations.Table;

import com.google.gwt.maps.client.InfoWindowContent;
import com.google.gwt.maps.client.MapWidget;
import com.google.gwt.maps.client.MapType;
import com.google.gwt.maps.client.control.LargeMapControl;
import com.google.gwt.maps.client.control.MapTypeControl;
import com.google.gwt.maps.client.geom.LatLng;
import com.google.gwt.maps.client.overlay.Marker;
import com.google.gwt.maps.client.overlay.MarkerOptions;

import com.google.gwt.maps.client.event.MarkerClickHandler;
import com.google.gwt.maps.client.event.MarkerMouseOverHandler;
import com.google.gwt.maps.client.InfoWindow;
import com.google.gwt.maps.client.InfoWindowContent;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONArray;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONString;
import com.google.gwt.json.client.JSONValue;

import java.util.HashMap;


public class MapWithPhenotype extends AbstractVisualization<MapWithPhenotype.CustomVisualizationDrawOptions> implements Selectable{
	/**
	 * Drawing options supported by this visualization.
	 */
	public static class CustomVisualizationDrawOptions extends
	AbstractDrawOptions {
		protected CustomVisualizationDrawOptions() {
		}
	}
	
	
	private AccessionConstants constants;
	private AbstractDataTable dataTable;

	private ListBox phenotypeSelectBox = new ListBox(false);
	//private RadioButton displayOptionRadioButton;
	private ListBox displayOptionSelectBox = new ListBox(false);
	private HorizontalPanel topHPanel = new HorizontalPanel();

	private MapWidget map;
	private VerticalPanel dialogVPanel = new VerticalPanel();
	private DisplayJSONObject jsonErrorDialog;

	private Button closeButton = new Button("Close");
	private static final String DIALOG_DEFAULT_TEXT = "Map With Phenotype";
	private static final String DIALOG_WAITING_TEXT = "Waiting...";

	private int latitude_idx = 6;
	private int longitude_idx = 7;
	private int nativename_idx = 2;
	private int ecotypeid_idx = 0;

	private final HashMap<Integer, Marker> rowIndex2Marker = new HashMap<Integer, Marker>();
	private int selectedRow=-1;
	private SelectHandler selectHandler;
	
	@Override
	public void draw(AbstractDataTable dataTable, CustomVisualizationDrawOptions options)
	{
		/*
		 * 2009-4-9 required for AbstractVisualization
		 */
		addMarkers(dataTable);
	}
	
	public MapWithPhenotype(AccessionConstants constants, DisplayJSONObject jsonErrorDialog) {
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		// Set the dialog box's caption.
		//setText(DIALOG_DEFAULT_TEXT);

		// setSize("800px", "600px");

		//fillPhenotypeSelectBox(phenotypeSelectBox);
		phenotypeSelectBox.addChangeListener(new ChangeListener() {
			public void onChange(Widget sender) {
				refreshMap(map, phenotypeSelectBox.getSelectedIndex(), displayOptionSelectBox.getSelectedIndex());
				//multiBox.ensureDebugId("cwListBox-multiBox");
			}
		});

		String[] displayOptions = {constants.MapWithPhenotypeDisplayOption1(), constants.MapWithPhenotypeDisplayOption2(), constants.MapWithPhenotypeDisplayOption3()};
		for (int i = 0; i < displayOptions.length; i++) {
			displayOptionSelectBox.addItem(displayOptions[i]);
		}
		displayOptionSelectBox.addChangeListener(new ChangeListener() {
			public void onChange(Widget sender) {
				refreshMap(map, phenotypeSelectBox.getSelectedIndex(), displayOptionSelectBox.getSelectedIndex());
				//multiBox.ensureDebugId("cwListBox-multiBox");
			}
		});
		topHPanel.add(phenotypeSelectBox);
		topHPanel.add(displayOptionSelectBox);

		map = new MapWidget();
		map.setSize("900px", "500px");

		// Add some controls for the zoom level
		map.addControl(new LargeMapControl());
		map.addControl(new MapTypeControl());
		map.addMapType(MapType.getPhysicalMap());	//2009-4-9 add the terrain map
		map.setCurrentMapType(MapType.getPhysicalMap());	//2009-4-9 set the terrain map as default
		map.setScrollWheelZoomEnabled(true);
		
		dialogVPanel.add(topHPanel);
		dialogVPanel.add(map);
		initWidget(dialogVPanel);

	}

	public void findLatLongCol(AbstractDataTable dataTable)
	{
		int no_of_cols = dataTable.getNumberOfColumns();

		for (int i =0; i<no_of_cols; i++)
		{
			String col_id = dataTable.getColumnId(i);
			if (col_id=="latitude")
				latitude_idx = i;
			else if (col_id=="longitude")
				longitude_idx = i;
			else if (col_id=="nativename")
				nativename_idx = i;
			else if (col_id=="tg_ecotypeid")
				ecotypeid_idx = i;
		}
	}

	public void addMarkers(AbstractDataTable dataTable)
	{
		findLatLongCol(dataTable);
		JSONObject ecotype_id2phenotype_value = null;
		double min_value = 0;
		double max_value = 0;
		int phenotype_method_id = -1;
		int displayOption = 0;
		addMarkers(dataTable, ecotype_id2phenotype_value, min_value, max_value, phenotype_method_id, displayOption);
	}

	public void addMarkers(AbstractDataTable dataTable, JSONObject ecotype_id2phenotype_value, double min_value, double max_value, int phenotype_method_id, int displayOption)
	{
		/*
		 * 1. create 3 data tables
		 * 2. for loop to add markers (with url for the icon) 
		 * 3. add click callback on each marker
		 */
		map.clearOverlays();

		double latitude;
		double longitude;
		if (ecotype_id2phenotype_value== null)
		{
			for (int i=0; i<dataTable.getNumberOfRows(); i++)
			{
				latitude = dataTable.getValueDouble(i, latitude_idx);
				longitude = dataTable.getValueDouble(i, longitude_idx);
				LatLng point = LatLng.newInstance(latitude, longitude);
				final String markerLabel = dataTable.getValueString(i, nativename_idx)+" ID: " + dataTable.getValueInt(i, ecotypeid_idx);
				final MarkerOptions markerOption = MarkerOptions.newInstance();
				markerOption.setTitle(markerLabel);	// title shows up as tooltip
				final Marker marker = new Marker(point, markerOption);
				final int rowIndex = i;
				rowIndex2Marker.put(rowIndex, marker);

				map.addOverlay(marker);
				marker.addMarkerClickHandler(new MarkerClickHandler() {
					public void onClick(MarkerClickEvent event) {
						InfoWindow info = map.getInfoWindow();
						info.open(marker, new InfoWindowContent(markerLabel));
						selectedRow = rowIndex;
						//Selection.triggerSelection(MapWithPhenotype.this, getSelections());	//2009-4-9  doesn't work
						if (selectHandler!=null)	//2009-4-9 check if selectHandler is initialized.
						{
							SelectEvent e = new SelectEvent();
							selectHandler.onSelect(e);
						}
						//MapWithPhenotype.this.fireSelectionEvent();	//2009-4-9 fireSelectionEvent() below doesn't work 
					}
				});
				/*
				marker.addMarkerMouseOverHandler(new MarkerMouseOverHandler(){
					public void onMouseOver(MarkerMouseOverEvent event) {
						InfoWindow info = map.getInfoWindow();
						info.open(marker, new InfoWindowContent(markerLabel));
					}
				});
				 */
			}
		}

	}

	private class PhenotypeChangeJSONResponseHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetSearchButtonCaption();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONValue jsonValue = JSONParser.parse(responseText);
				JSONObject jsonObject = jsonValue.isObject();
				double min_value = jsonObject.get("min_value").isNumber().doubleValue();
				double max_value = jsonObject.get("max_value").isNumber().doubleValue();
				JSONObject ecotype_id2phenotype_value = jsonObject.get("ecotype_id2phenotype_value").isObject();
				addMarkers(dataTable, ecotype_id2phenotype_value, min_value, max_value, phenotypeSelectBox.getSelectedIndex(), displayOptionSelectBox.getSelectedIndex());

				//displayJSONObject(jsonValue);

			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetSearchButtonCaption();
		}
	}

	public void refreshMap(MapWidget map, int phenotype_method_id, int displayOption)
	{
		/*
		 * 1. fetch data (phenotype-value + icon urls)
		 * 3. addMarkers
		 */
		//setText(DIALOG_WAITING_TEXT);
		String url = URL.encode(constants.GetPhenotypeValueURL() + "/" + phenotype_method_id);
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new PhenotypeChangeJSONResponseHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
	}

	private class FillPhenotypeJSONResponseHandler implements RequestCallback {
		public void onError(Request request, Throwable exception) {
			jsonErrorDialog.displayRequestError(exception.toString());
			resetSearchButtonCaption();
		}

		public void onResponseReceived(Request request, Response response) {
			String responseText = response.getText();
			try {
				JSONArray phenotypeMethodArray = JSONParser.parse(responseText).isArray();
				for (int i = 0; i < phenotypeMethodArray.size(); i++) {
					JSONArray phenotypeTuple = phenotypeMethodArray.get(i).isArray();
					String phenotype_id = phenotypeTuple.get(0).isString().stringValue();
					String phenotype_name = phenotypeTuple.get(1).isString().stringValue();
					phenotypeSelectBox.addItem(phenotype_name, phenotype_id);

				}

			} catch (JSONException e) {
				jsonErrorDialog.displayParseError(responseText);
			}
			resetSearchButtonCaption();
		}
	}

	public void fillPhenotypeSelectBox(ListBox phenotypeSelectBox)
	{
		//setText(DIALOG_WAITING_TEXT);
		String url = URL.encode(constants.MapWithPhenotypeGetPhenotypeMethodLsURL());
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			requestBuilder.sendRequest(null, new FillPhenotypeJSONResponseHandler());
		} catch (RequestException ex) {
			jsonErrorDialog.displaySendError(ex.toString());
			resetSearchButtonCaption();
		}
	}	
	private void resetSearchButtonCaption() {
		//setText(DIALOG_DEFAULT_TEXT);
	}
	public void setSelections(JsArray<Selection> s)
	{
		for (int i = 0; i < s.length(); ++i) {
			Integer rowIndex=null;
			if (s.get(i).isCell()) {
				rowIndex = s.get(i).getRow();
			} else if (s.get(i).isRow()) {
				rowIndex = s.get(i).getRow();
			}
			if (rowIndex!=null)
			{
				Marker marker = rowIndex2Marker.get(rowIndex);
				map.setCenter(marker.getLatLng());
				InfoWindow info = map.getInfoWindow();
				info.open(marker, new InfoWindowContent(marker.getTitle()));
			}
		}
	}
	public final void addSelectHandler(SelectHandler handler)
	{
		this.selectHandler = handler;	//2009-4-9 since fireSelectionEvent() doesn't work, stuff below
		//SelectEvent event = new SelectEvent();
		//handler.onSelect(event);
		//Selection.addSelectHandler(this, handler);	//2009-4-9 doesn't work
	}
	
	public JsArray<Selection> getSelections()
	{
		//JsArray<Selection> sel = new JsArray<Selection>();
		return ArrayHelper.toJsArray(Selection.createRowSelection(selectedRow));
	}
}