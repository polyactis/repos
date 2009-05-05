package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.AbstractVisualization.VisualizationFactory;

import com.google.gwt.user.client.History;
import com.google.gwt.user.client.HistoryListener;

import com.google.gwt.user.client.Window;
/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Accession implements EntryPoint, HistoryListener{

	/**
	 * This is the entry point method.
	 */
	protected SinkTabPanel tp;
	private int curSinkIndex=-1;
	
	public DisplayJSONObject jsonErrorDialog;
	public AccessionConstants constants;
	public class ToySink extends Sink
	{
		private String name;
		private HTML widget;
		public ToySink(String name, String content)
		{
			this.name = name;
			widget  = new HTML(content);
			initWidget(widget);
		}
		public String getName()
		{
			return name;
		}

	}
	public void onModuleLoad() {

		// Create the constants
		
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		/*
		AbstractVisualization.registerVisualization("MapWithPhenotype",
				new VisualizationFactory() {

			public AbstractVisualization<?> create() {
				return new MapWithPhenotype(constants, jsonErrorDialog);
			}
		});
		*/
		
		tp = new SinkTabPanel(jsonErrorDialog);
		tp.addSink(new Accession250k(constants, jsonErrorDialog));
		tp.addSink(new AccessionByName(constants, jsonErrorDialog));
		tp.addSink(new AccessionByID(constants, jsonErrorDialog));
		//tp.addSink(new ToySink("By Genetic Distance", "Under construction"));
		//tp.addSink(new ToySink("By Geographic Distance", "Under construction"));
		// tp.addSink(new ToySink("By Country", "Under construction"));
		// Show the 'bar' tab initially.
		tp.selectTab(1);
		tp.setWidth("1000px");
		tp.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("accession").add(tp);
		
		//		Add history listener
		
		History.addHistoryListener(this);
		

		String initToken = History.getToken();
		if (initToken.length() > 0) {
			onHistoryChanged(initToken);
		} else {
			showInfo();
		}
		
		
		// Now that we've setup our listener, fire the initial history state.
		//History.fireCurrentHistoryState();
		
	}
	public void onHistoryChanged(String token) {
		// This method is called whenever the application's history changes. Set
		// the label to reflect the current history token.
		//lbl.setText("The current history token is: " + historyToken);
		//	 Find the MapsDemoInfo associated with the history context. If one is
	    // found, show it (It may not be found, for example, when the user mis-
	    // types a URL, or on startup, when the first context will be "").
	    
		//jsonErrorDialog.displayRequestError("The current history token is: " + token);
		
		int i = tp.find(token);
	    if (i == -1) {
	      showInfo();
	      Window.alert("Couldn't find " + token);
	      return;
	    }
	    show(i, false);
	}

	public void show(int selectedSink, boolean affectHistory) {
		// Don't bother re-displaying the existing MapsDemo. This can be an issue
		// in practice, because when the history context is set, our
		// onHistoryChanged() handler will attempt to show the currently-visible
		// MapsDemo.
		if (selectedSink == curSinkIndex) {
			return;
		}
		curSinkIndex = selectedSink;
		
		tp.setSinkSelection(curSinkIndex);
		
		// If affectHistory is set, create a new item on the history stack. This
		// will ultimately result in onHistoryChanged() being called. It will call
		// show() again, but nothing will happen because it will request the exact
		// same MapsDemo we're already showing.
		if (affectHistory) {
			History.newItem(tp.getSinkName(curSinkIndex));
		}

		//curSink.onShow();
	}
	
	private void showInfo() {
		int i = tp.find("By Name");
		show(i, false);
	}
}
