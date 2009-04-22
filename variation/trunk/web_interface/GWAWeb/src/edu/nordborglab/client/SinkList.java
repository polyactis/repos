/*
 * 2009-4-14
 * modelled after DemoList.java from the HelloMaps sample in gwt-maps
 */

package edu.nordborglab.client;

import com.google.gwt.user.client.History;
import com.google.gwt.user.client.ui.TabListener;
import com.google.gwt.user.client.ui.SourcesTabEvents;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.Widget;

import java.util.ArrayList;
import java.util.HashMap;
/**
 * The left panel that contains all of the sinks, along with a short description
 * of each.
 */
public class SinkList extends Composite {

	//	private HorizontalPanel list = new HorizontalPanel();
	private HashMap<String, Integer> sinkName2index = new HashMap<String, Integer>();
	private ArrayList<Sink> list = new ArrayList<Sink>();
	
	private int selectedSinkIndex = -1;
	public DisplayJSONObject jsonErrorDialog;
	
	public SinkList(DisplayJSONObject jsonErrorDialog) {
		this.jsonErrorDialog = jsonErrorDialog;
		//initWidget(tp);
		//this.addTabListener((new tabSelectListener()));
	}

	public int addSink(final Sink sink) {
		String name = sink.getName();
		if (sinkName2index.containsKey(name))
		{
			return sinkName2index.get(name).intValue();
		}
		else
		{
			list.add(sink);
			sinkName2index.put(name, list.size()-1);
			return sinkName2index.get(name).intValue();
		}
	}

	public int find(String sinkName) {
		//jsonErrorDialog.displayRequestError(sinkName);
		//Integer no_of_widgets = (Integer)this.getWidgetCount();
		//jsonErrorDialog.displayRequestError(no_of_widgets.toString());
		
		if (sinkName2index.containsKey(sinkName))
		{
			return sinkName2index.get(sinkName);			
		}
		else
			return -1;
	}
	
	public Sink getSelectedSink(int i) {
		selectedSinkIndex = i;
		return list.get(selectedSinkIndex);
	}
	
	public int getNumberOfSinks()
	{
		return list.size();
	}
	public String getSinkName(int i)
	{
		return list.get(i).getName();
	}
	public Sink getLastSink()
	{
		return list.get(getNumberOfSinks()-1);
	}
}