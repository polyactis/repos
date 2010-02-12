/**
 * a customized scatterchart used to display association mapping results.

 *   add a select event handler which leads to the opening of a new window showing SNP information
 *   add mouseover/mouseout handler to display popup
 *   sink the MOUSEEVENTS in order to pass X,Y to visualization's mouseover/mouseout event handler
 */
package edu.nordborglab.client;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.http.client.URL;

import com.google.gwt.user.client.DOM;
import com.google.gwt.user.client.Event;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.DecoratedPopupPanel;
import com.google.gwt.user.client.ui.FlexTable;
import com.google.gwt.user.client.ui.MouseListener;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;
import com.google.gwt.visualization.client.visualizations.ScatterChart;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.LegendPosition;
import com.google.gwt.visualization.client.events.OnMouseOverHandler;
import com.google.gwt.visualization.client.events.OnMouseOutHandler;

public class AssociationScatterChart extends ScatterChart {
	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	
	//private ScatterChart scatterChart;
	
	private String chromosome;
	private String color;
	private int chrLength;
	private int maxChrLength;
	private double maxY;
	private String SNPBaseURL;
	
	private int pointSize =3;
	
	private AbstractDataTable dataTable;
	
	private DecoratedPopupPanel popup = new DecoratedPopupPanel(true);
	private int mouseX;
	private int mouseY;
	private MouseListener mListener;
	
	AssociationScatterChart(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, String chromosome,
			String color, int chrLength, int maxChrLength, double maxY, String SNPBaseURL)
	{
		super();
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		this.chromosome = chromosome;
		this.color = color;
		this.chrLength = chrLength;
		this.maxChrLength = maxChrLength;
		this.maxY = maxY;
		this.SNPBaseURL = SNPBaseURL;
		
		//scatterChart = new ScatterChart();
		addSelectHandler(new ScatterSelectionHandler(this));
		addOnMouseOverHandler(new ScatterMouseOverHandler());
		addOnMouseOutHandler(new ScatterMouseOutHandler());
		sinkEvents(Event.MOUSEEVENTS);	//2009-5-21 sink the mouse event because visualization's MouseOverHandler() doesn't have x,y coordinates passed
		
	}
	
	public void onBrowserEvent(Event event) {
		super.onBrowserEvent(event);
		switch (DOM.eventGetType(event)) {
		case Event.ONMOUSEOVER:
			// mListener.onMouseMove(this.getParent(), evt.getClientX(), evt.getClientY());
			//Window.alert("mouseX ");
			mouseX = event.getClientX();
			mouseY = event.getClientY();
			//Window.alert("mouseX " + mouseX + " mouseY " + mouseY);
			break;
		case Event.ONMOUSEMOVE:
			// mListener.onMouseMove(this.getParent(), evt.getClientX(), evt.getClientY());
			//Window.alert("mouseX ");
			mouseX = event.getClientX();
			mouseY = event.getClientY();
			//Window.alert("mouseX " + mouseX + " mouseY " + mouseY);
			break;
		}
	}
	
	public void addMouseListener(MouseListener mListener)
	{
		this.mListener = mListener;
	}
	
	public void draw_gwas(AbstractDataTable data){
		this.dataTable = data;
		Options options = Options.create();
		options = setOptions(options);
		this.draw(this.dataTable, options);
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
	public void draw_gwas(AbstractDataTable data, Options options){
		this.dataTable = data;
		options = setOptions(options);
		this.draw(data, options);
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
	
	/**
	 * 2009-5-21 display the popup reporting position and score upon mouse over the point
	 * @author crocea
	 *
	 */
	class ScatterMouseOverHandler extends OnMouseOverHandler {
		ScatterMouseOverHandler() {
		}

		@Override
		public void onMouseOverEvent(OnMouseOverHandler.OnMouseOverEvent event)
		{
			//int left = mouseX + Window.getScrollLeft();
			//int top = mouseY + Window.getScrollTop();
			popup.setPopupPosition(10, Window.getScrollTop()+10);
			//popup.center();
			popup.clear();
			int position = dataTable.getValueInt(event.getRow(), 0);
			double score = dataTable.getValueDouble(event.getRow(), 1);
			FlexTable layout = new FlexTable();
			layout.setWidget(0, 0, new HTML("Position: "+ position));
			layout.setWidget(1, 0, new HTML("-log Pvalue: "+ score));
			/*
			layout.setWidget(2, 0, new HTML("Row: "+ event.getRow()));
			layout.setWidget(3, 0, new HTML("Column: "+ event.getColumn()));
			layout.setWidget(4, 0, new HTML("mouseX: "+ mouseX));
			layout.setWidget(5, 0, new HTML("mouseY: "+ mouseY));
			*/
			popup.add(layout);
			// Show the popup
			popup.show();
		}
	}
	
	/**
	 * 2009-5-21 hide the popup when the mouse leaves the scatter point
	 * @author crocea
	 *
	 */
	class ScatterMouseOutHandler extends OnMouseOutHandler {
		ScatterMouseOutHandler() {
		}

		@Override
		public void onMouseOutEvent(OnMouseOutHandler.OnMouseOutEvent event)
		{
			// Hide the popup
			popup.hide();
		}
	}
	
	/**
	 * 2009-5-21 a mouse listener just to record the x,y positions of the mouse,
	 * 		which is done onBrowserEvent(). this handler is not used.
	 * @author crocea
	 *
	 */
	private class ScatterMouseListener implements MouseListener
	{
		public void onMouseMove(Widget sender, int x, int y)
		{
			mouseX = x;
			mouseY = y;
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
		}
	}
	
	
}
