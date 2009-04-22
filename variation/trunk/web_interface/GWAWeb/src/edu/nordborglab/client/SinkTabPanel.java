/*
 * 2009-4-14
 * copied from DemoList.java from the HelloMaps sample in gwt-maps
 */

package edu.nordborglab.client;

import com.google.gwt.user.client.History;
import com.google.gwt.user.client.ui.TabListener;
import com.google.gwt.user.client.ui.SourcesTabEvents;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.Widget;

import java.util.ArrayList;

/**
 * The left panel that contains all of the sinks, along with a short description
 * of each.
 */
public class SinkTabPanel extends TabPanel {

	//	private HorizontalPanel list = new HorizontalPanel();
	//private TabPanel tp = new TabPanel();
	private ArrayList<String> sinkNames = new ArrayList<String>();
	
	private int selectedSink = -1;
	public DisplayJSONObject jsonErrorDialog;
	public class tabSelectListener implements TabListener
	{
		public boolean onBeforeTabSelected(SourcesTabEvents sender, int tabIndex)
		{
			//History.newItem(sinkNames.get(tabIndex));
			return true;
		}
		
		public void onTabSelected(SourcesTabEvents sender, int tabIndex)
		{
			History.newItem(sinkNames.get(tabIndex));
		}
	}
	
	public SinkTabPanel(DisplayJSONObject jsonErrorDialog) {
		this.jsonErrorDialog = jsonErrorDialog;
		//initWidget(tp);
		this.addTabListener((new tabSelectListener()));
	}

	public void addSink(final Sink sink) {
		String name = sink.getName();
		this.add(sink, name);
		//list.addItem(name);
		sinkNames.add(name);
	}

	public int find(String sinkName) {
		//jsonErrorDialog.displayRequestError(sinkName);
		//Integer no_of_widgets = (Integer)this.getWidgetCount();
		//jsonErrorDialog.displayRequestError(no_of_widgets.toString());
		for (int i = 0; i < sinkNames.size(); ++i) {
			//Sink sink = (Sink) this.getWidget(i);
			//jsonErrorDialog.displayRequestError(sink.getName());
			if (sinkNames.get(i).equals(sinkName)) {
				return i;
			}
		}
		return -1;
	}

	public String getSinkName(int i)
	{
		return sinkNames.get(i);
	}
	
	public void setSinkSelection(int i) {
		this.selectTab(i);
		Sink sink = (Sink)this.getWidget(i);
		sink.resetSize();	//2009-4-17 in most sinks, this resizes the map size and would avoid map partial blank.
		/*
		for (int i = 0; i < this.getWidgetCount(); ++i) {
			Sink sink = (Sink) this.getWidget(i);
			if (sink.getName().equals(name)) {
				selectedSink = i;
				this.selectTab(selectedSink);
				return;
			}
		}
		*/
		
		/*		
		for (int i = 0; i < sinks.size(); ++i) {
			Sink sink = (Sink) sinks.get(i);
			if (sink.getName().equals(name)) {
				selectedSink = i;
				this.selectTab(selectedSink);
				return;
			}
		}

		for (int i = 0; i < list.getItemCount(); i++) {
			if (name.equals(list.getItemText(i))) {
				list.setSelectedIndex(i);
				break;
			}
		}
		*/
	}
}