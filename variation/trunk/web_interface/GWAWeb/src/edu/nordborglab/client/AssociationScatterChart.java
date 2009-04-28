package edu.nordborglab.client;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.http.client.URL;

import com.google.gwt.user.client.Window;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.ArrayHelper;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;
import com.google.gwt.visualization.client.visualizations.ScatterChart;
import com.google.gwt.visualization.client.visualizations.ScatterChart.Options;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.LegendPosition;


public class AssociationScatterChart extends AbstractVisualization<Options> implements Selectable {
	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	
	private ScatterChart scatterChart;
	
	private String chromosome;
	private String color;
	private int chrLength;
	private int maxChrLength;
	private double maxY;
	private String SNPBaseURL;
	
	private int pointSize =3;
	
	private AbstractDataTable dataTable;
	
	AssociationScatterChart(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, String chromosome,
			String color, int chrLength, int maxChrLength, double maxY, String SNPBaseURL)
	{
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		this.chromosome = chromosome;
		this.color = color;
		this.chrLength = chrLength;
		this.maxChrLength = maxChrLength;
		this.maxY = maxY;
		this.SNPBaseURL = SNPBaseURL;
		
		scatterChart = new ScatterChart();
		scatterChart.addSelectHandler(new ScatterSelectionHandler(scatterChart));
		initWidget(scatterChart);
	}
	
	public void draw(AbstractDataTable data){
		this.dataTable = data;
		Options options = Options.create();
		options = setOptions(options);
		this.scatterChart.draw(this.dataTable, options);
	}
	public Options setOptions(Options options){
		options.setHeight(200);
		//options.setWidth((int) 1000*this.chrLength/this.maxChrLength);	// adjust chromosome length according to maximum chromosome length
		options.setWidth(1000);
		options.setColors(this.color);
		options.setOption("max", this.maxY);
		options.setOption("titleX", "Chr "+this.chromosome);
		options.setOption("titleY", "-log Pvalue");
		options.setPointSize(this.pointSize);
		options.setLegend(LegendPosition.NONE);
		return options;
	}
	public void draw(AbstractDataTable data, Options options){
		this.dataTable = data;
		options = setOptions(options);
		this.scatterChart.draw(data, options);
	}
	
	class ScatterSelectionHandler extends SelectHandler {
		private final Selectable viz;

		ScatterSelectionHandler(Selectable viz) {
			this.viz = viz;
		}

		@Override
		public void onSelect(SelectEvent event) {
			//open a new window pointing to the SNP page
			JsArray<Selection> selectionLs = getSelections();
			for (int i=0; i< selectionLs.length(); i++ )
			{
				Selection s = selectionLs.get(i);
				int position = dataTable.getValueInt(s.getRow(), 0);
				double score = dataTable.getValueDouble(s.getRow(), 1);
				final String _SNPURL = URL.encode(SNPBaseURL + "&chromosome="+chromosome+"&position=" + position +"&score="+score);
				Window.open(_SNPURL, "", "");
			}
		}

		private JsArray<Selection> getSelections() {
			return viz.getSelections();
		}
	}
	
	public void setSelections(JsArray<Selection> s)
	{
		this.scatterChart.setSelections(s);
	}
	public final void addSelectHandler(SelectHandler handler)
	{
		this.scatterChart.addSelectHandler(handler);
		//this.selectHandler = handler;	//2009-4-9 since fireSelectionEvent() doesn't work, stuff below
		//SelectEvent event = new SelectEvent();
		//handler.onSelect(event);
		//Selection.addSelectHandler(this, handler);	//2009-4-9 doesn't work
	}
	
	public JsArray<Selection> getSelections()
	{
		//return ArrayHelper.toJsArray(Selection.createRowSelection(selectedRow));
		return this.scatterChart.getSelections();
	}
	
	
}
