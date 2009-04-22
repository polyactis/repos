package edu.nordborglab.client;

import com.google.gwt.user.client.Event;

import com.google.gwt.user.client.DOM;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Widget;

import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Table.Options.Policy;

public class CustomVisualizationTable extends Table{
	private CustomClickListener cListener;
	
	public CustomVisualizationTable() {
		super();
		//setElement(DOM.createDiv());
		sinkEvents(Event.ONCLICK);
	}

	public void onBrowserEvent(Event evt) {
		switch (DOM.eventGetType(evt)) {
		case Event.ONCLICK:
			//evt.getClientX();
			cListener.onClick(this.getParent(), evt);
			break;
		}
	}
	
	public void addClickListener(CustomClickListener cListener)
	{
		this.cListener = cListener;
	}
}