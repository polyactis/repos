package edu.nordborglab.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.visualization.client.AbstractDataTable;
/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class PhenotypeTrend implements EntryPoint {
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private HTML status1;
	private HTML status2;
	
	private String phenotypeBaseURL;
	private String trendDataURL;
	private String trendAnnotationDataURL;
	private String multiPhenotypeDataURL;
	
	private CustomVerticalPanel trendVPanel;
	private CustomVerticalPanel multiPhenotypeVPanel;
	
	private MotionChart trendMotionChart;
	private AbstractDataTable trendData;
	private CustomVisualizationTable trendAnnotationTable;
	private AbstractDataTable trendAnnotationData;
	private DisclosurePanel trendAnnotationTablePanel;
	private MotionChart multiPhenotypeMotionChart;
	private AbstractDataTable multiPhenotypeData;
	
	private TabPanel tPanel;
	private RootPanel rootPanel;
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		
		constants = (AccessionConstants) GWT.create(AccessionConstants.class);
		jsonErrorDialog = new DisplayJSONObject("Error Dialog");
		
		rootPanel = RootPanel.get("phenotypeTrend");
		
		tPanel = new TabPanel();
		tPanel.setAnimationEnabled(true);
		
		rootPanel.add(tPanel);
		
		phenotypeBaseURL = getPhenotypeBaseURL();
		trendDataURL = getTrendDataURL();
		trendAnnotationDataURL = getTrendAnnotationDataURL();
		multiPhenotypeDataURL = getMultiPhenotypeDataURL();
		
		trendVPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.AccessionByIDHelpID());
		tPanel.add(trendVPanel, "trend");
		multiPhenotypeVPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.AccessionByIDHelpID());
		tPanel.add(multiPhenotypeVPanel, "Match Two Phenotypes");
		tPanel.selectTab(1);
		
		status1 = new HTML(constants.LoadingText());
		trendVPanel.add(status1);
		
		status2 = new HTML(constants.LoadingText());
		multiPhenotypeVPanel.add(status2);
		
		trendMotionChart = new MotionChart();
		trendVPanel.add(trendMotionChart);
		trendAnnotationTable = new CustomVisualizationTable();

		trendAnnotationTablePanel = new DisclosurePanel("Phenotype Annotation Table");
		trendAnnotationTablePanel.setAnimationEnabled(true);
		trendAnnotationTablePanel.setContent(trendAnnotationTable);
		trendAnnotationTablePanel.setOpen(true);
		trendVPanel.add(trendAnnotationTablePanel);

		multiPhenotypeMotionChart = new MotionChart();
		multiPhenotypeVPanel.add(multiPhenotypeMotionChart);
		
		//Common.fillInMotionChart(trendDataURL, trendMotionChart, jsonErrorDialog, status1, trendData);
		Common.fillInMotionChart(multiPhenotypeDataURL, multiPhenotypeMotionChart, jsonErrorDialog, status2, multiPhenotypeData);
		Common.fillInTable(trendAnnotationDataURL, trendAnnotationTable, jsonErrorDialog, status1, trendAnnotationData);
	}
	
	public native String getPhenotypeBaseURL()/*-{ return $wnd.phenotypeBaseURL; }-*/;
	public native String getTrendDataURL()/*-{ return $wnd.trendDataURL; }-*/;
	
	public native String getTrendAnnotationDataURL()/*-{ return $wnd.trendAnnotationDataURL; }-*/;	
	
	public native String getMultiPhenotypeDataURL()/*-{ return $wnd.multiPhenotypeDataURL; }-*/;

}
