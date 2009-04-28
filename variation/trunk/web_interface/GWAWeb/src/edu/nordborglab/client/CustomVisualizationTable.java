/*
 * 2009-4-26 This custom google visualization table sinks the click event (besides the inherent select event).
 * 		One unique feature about calling the click handler here is the Event is passed along with the Widget. 
 */

package edu.nordborglab.client;

import com.google.gwt.user.client.Event;

import com.google.gwt.user.client.DOM;

import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;

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
			cListener.onClick(this.getParent(), evt);
			break;
		}
	}
	
	public void addClickListener(CustomClickListener cListener)
	{
		this.cListener = cListener;
	}
}