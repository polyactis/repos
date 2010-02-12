/**
 * 2009-7-2
 * class serves as a tab in OnePhenotypeGWASGene
 */
package edu.nordborglab.client;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.ChangeListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Hyperlink;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.visualizations.ColumnChart;
import com.google.gwt.visualization.client.visualizations.MotionChart;
import com.google.gwt.visualization.client.visualizations.Table;

public class OneGWASGene extends CustomVerticalPanel {
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;

	public ListBox callMethodListBox;
	public ListBox phenotypeMethodListBox;
	public ListBox analysisMethodListBox;
	public ListBox snpGeneAssociationTypeListBox;

	public Button submitButton;
	public TextBox topNumberTxtBox;
	public Button refreshButton;
	public ListBox geneListTypeListBox;
	public HTML linkToCandidateGeneList;
	private HTML statusReport;
	
	private AbstractDataTable resultsGeneDataTable;
	private WrapVisualizationTable resultsGeneTable;
	
	private VerticalPanel formVPanel;
	private DisclosurePanel formPanel;
	private DisclosurePanel optionPanel;
	private DisclosurePanel tablePanel;

	private String callMethodID;	//not used
	private String phenotypeMethodID;
	private String analysisMethodID;
	private String snpGeneAssociationTypeID;

	private String snpGeneAssociationTypeLsURL;
	private String snpGeneAssociationOnChangeURL;
	private String phenotypeMethodOnChangeURL;
	private String callMethodLsURL;
	private String geneListLsURL;
	private String callMethodOnChangeURL;
	private String fetchResultsGeneURL;
	private String candidateGeneListURL;

	private static final String SUBMIT_BUTTON_DEFAULT_TEXT = "Refresh";
	private static final String SUBMIT_BUTTON_WAITING_TEXT = "Waiting...";

	OneGWASGene(AccessionConstants constants,
			DisplayJSONObject jsonErrorDialog, String callMethodID,
			String phenotypeMethodID, String analysisMethodID, String typeID) {
		super(constants, jsonErrorDialog, constants.OneGWASGeneHelpID());
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;

		this.callMethodID = callMethodID;	//not used
		this.phenotypeMethodID = phenotypeMethodID;
		this.analysisMethodID = analysisMethodID;
		this.snpGeneAssociationTypeID = typeID;

		formVPanel = new VerticalPanel();
		/*
		 * callMethodListBox = new ListBox(); HorizontalPanel hpanel1 = new
		 * HorizontalPanel(); hpanel1.add(new
		 * HTML(constants.callMethodLabel())); hpanel1.add(callMethodListBox);
		 * hpanel1.setSpacing(5); formVPanel.add(hpanel1);
		 */
		
		snpGeneAssociationTypeListBox = new ListBox();
		HorizontalPanel hpanelSNPGeneAssociationType = new HorizontalPanel();
		hpanelSNPGeneAssociationType .add(new HTML(constants.snpGeneAssociationTypeLabel()));
		hpanelSNPGeneAssociationType.add(snpGeneAssociationTypeListBox);
		hpanelSNPGeneAssociationType.setSpacing(5);
		formVPanel.add(hpanelSNPGeneAssociationType);
		
		phenotypeMethodListBox = new ListBox();
		HorizontalPanel hpanel2 = new HorizontalPanel();
		hpanel2.add(new HTML(constants.phenotypeMethodLabel()));
		hpanel2.add(phenotypeMethodListBox);
		hpanel2.setSpacing(5);
		formVPanel.add(hpanel2);

		analysisMethodListBox = new ListBox();
		HorizontalPanel hpanel3 = new HorizontalPanel();
		hpanel3.add(new HTML(constants.analysisMethodLabel()));
		hpanel3.add(analysisMethodListBox);
		hpanel3.setSpacing(5);
		formVPanel.add(hpanel3);

		formPanel = new DisclosurePanel("Form");
		formPanel.setAnimationEnabled(false);
		formPanel.setContent(formVPanel);
		formPanel.setOpen(false);
		
		add(formPanel);
		/*
		optionPanel = new DisclosurePanel();
		optionPanel.setAnimationEnabled(true);
		optionPanel.setOpen(true);
		 */
		topNumberTxtBox = new TextBox();
		topNumberTxtBox.setText("50");
		topNumberTxtBox.setWidth("50px");
		submitButton = new Button();
		submitButton.setText(SUBMIT_BUTTON_DEFAULT_TEXT);
		statusReport = new HTML(constants.LoadingText());
		statusReport.setVisible(false);

		HorizontalPanel hpanelTopNumber = new HorizontalPanel();
		hpanelTopNumber.add(new HTML("Top Number: "));
		hpanelTopNumber.add(topNumberTxtBox);
		hpanelTopNumber.add(submitButton);
		hpanelTopNumber.add(statusReport);
		hpanelTopNumber.setSpacing(5);
		add(hpanelTopNumber);

		geneListTypeListBox = new ListBox();
		HorizontalPanel hpanelGeneList = new HorizontalPanel();
		hpanelGeneList.add(new HTML(constants.geneListTypeLabel()));
		hpanelGeneList.add(geneListTypeListBox);
		linkToCandidateGeneList = new HTML();
		linkToCandidateGeneList.setVisible(false);
		hpanelGeneList.add(linkToCandidateGeneList);
		hpanelGeneList.setSpacing(5);

		add(hpanelGeneList);

		resultsGeneTable = new WrapVisualizationTable();

		//add(optionPanel);
		add(resultsGeneTable);

		snpGeneAssociationTypeLsURL = Common.getSNPGeneAssociationTypeLsURL();
		snpGeneAssociationOnChangeURL = Common.getSNPGeneAssociationOnChangeURL();
		phenotypeMethodOnChangeURL = Common.getPhenotypeMethodOnChangeURL();
		fetchResultsGeneURL = Common.getFetchResultsGeneURL();
		geneListLsURL = Common.getGeneListLsURL();
		candidateGeneListURL = Common.getCandidateGeneListURL();
		
		Common.fillSelectBox(snpGeneAssociationTypeLsURL,
				snpGeneAssociationTypeListBox, jsonErrorDialog);
		Common.fillSelectBox(geneListLsURL, geneListTypeListBox,
				jsonErrorDialog);

		snpGeneAssociationTypeListBox.addChangeListener(new Common.FillTargetListOnListChangeListener(
						snpGeneAssociationOnChangeURL, "type_id", phenotypeMethodListBox,
						jsonErrorDialog));
		// callMethodListBox.addChangeListener(new
		// Common.ListBoxChangeListener(callMethodOnChangeURL,
		// phenotypeMethodListBox, jsonErrorDialog));
		phenotypeMethodListBox.addChangeListener(new Common.FillTargetListBasedOnTwoListsChange(
						phenotypeMethodOnChangeURL, snpGeneAssociationTypeListBox, "type_id", phenotypeMethodListBox, 
						"phenotype_method_id", analysisMethodListBox, jsonErrorDialog));
		
		analysisMethodListBox.addChangeListener(new FillWrapVTableOnListChangeListener(fetchResultsGeneURL, resultsGeneTable, jsonErrorDialog,
				statusReport, resultsGeneDataTable));
		submitButton.addClickListener(new FillWrapVTableOnClickListener(fetchResultsGeneURL, resultsGeneTable, jsonErrorDialog,
				statusReport, resultsGeneDataTable));
		geneListTypeListBox.addChangeListener(new FillWrapVTableOnListChangeListener(fetchResultsGeneURL, resultsGeneTable, jsonErrorDialog,
				statusReport, resultsGeneDataTable));
		
		// 2009-7-1 one more listener for geneListTypeListBox to update the link to the candidate gene list
		geneListTypeListBox.addChangeListener(new ChangeListener() {
			public void onChange(Widget sender) {
				ListBox senderListBox = (ListBox) sender;
				String selectedValue = senderListBox.getValue(senderListBox.getSelectedIndex());
				if (selectedValue.equals("0"))
				{
					linkToCandidateGeneList.setVisible(false);
				}
				else
				{
					String url = candidateGeneListURL + "?list_type_id=" + selectedValue;
					linkToCandidateGeneList.setHTML("<a href="+ url +"  target='_blank'>Link To Gene List "+selectedValue+"</a>");
					linkToCandidateGeneList.setVisible(true);
				}
			}
		});
		
		
		if ((snpGeneAssociationTypeID!="0"||snpGeneAssociationTypeID!=null) && (phenotypeMethodID !="0" ||phenotypeMethodID!=null)
				&& (analysisMethodID != "0" || analysisMethodID!=null)) {
			// set various values
			// fetch the data and fill the table
			// 2009-7-2 doesn't work cuz at this moment all list boxes are empty (http requests haven't been issued yet. 
			Common.setCertainValuedItemSelectedInListBox(snpGeneAssociationTypeListBox, snpGeneAssociationTypeID);
			Common.setCertainValuedItemSelectedInListBox(phenotypeMethodListBox, phenotypeMethodID);
			Common.setCertainValuedItemSelectedInListBox(analysisMethodListBox, analysisMethodID);
			String urlArguments = _constructURLArgument();
			if (!urlArguments.isEmpty())
			{
				String url = fetchResultsGeneURL + "?"+ urlArguments;
				Common.fillInTable(url, resultsGeneTable, jsonErrorDialog, statusReport, resultsGeneDataTable);
			}
		}
		
	}
	
	/**
	 * 2009-7-1 fill the WrapVisualizationTable based on the selected value of a list Box.
	 * 	the value must be String. baseURL contains no arguments at all.
	 * 
	 * @author crocea
	 *
	 */
	public class FillWrapVTableOnListChangeListener extends Common.FillWrapVTableOnListChangeListener
	{
		String _url;
		private AbstractVisualization<Table.Options> abstractVTable;
		DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		FillWrapVTableOnListChangeListener(String baseURL, AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, 
					Widget statusReport, AbstractDataTable dataTable)
		{
			super(baseURL, table, jsonErrorDialog, statusReport, dataTable);
		}
		public String constructURLArgument()
		{
			return _constructURLArgument();
		}
	}
	
	/**
	 * 2009-7-1 fill the WrapVisualizationTable when the button is clicked.
	 * 	the value must be String. baseURL contains no arguments at all.
	 * @author crocea
	 *
	 */
	public class FillWrapVTableOnClickListener extends Common.FillWrapVTableOnClickListener
	{
		String _url;
		private AbstractVisualization<Table.Options> abstractVTable;
		DisplayJSONObject jsonErrorDialog;
		private Widget statusReport;
		private AbstractDataTable dataTable;
		FillWrapVTableOnClickListener(String baseURL, AbstractVisualization<Table.Options> table, DisplayJSONObject jsonErrorDialog, 
					Widget statusReport, AbstractDataTable dataTable)
		{
			super(baseURL, table, jsonErrorDialog, statusReport, dataTable);
		}
		public String constructURLArgument()
		{
			return _constructURLArgument();
		}
	}
	
	public String _constructURLArgument()
	{
		String _url = "";
		
		String newArgument;
		newArgument = Common.getURLFromSelectBox(_url, "type_id", snpGeneAssociationTypeListBox, jsonErrorDialog, snpGeneAssociationTypeID);
		if (newArgument.isEmpty())
			return "";
		else
			_url = newArgument;
		
		newArgument = Common.getURLFromSelectBox(_url, "phenotype_method_id", phenotypeMethodListBox, jsonErrorDialog, phenotypeMethodID);
		if (newArgument.isEmpty())
			return "";
		else
			_url = newArgument;
		
		newArgument = Common.getURLFromSelectBox(_url, "analysis_method_id", analysisMethodListBox, jsonErrorDialog, analysisMethodID);
		if (newArgument.isEmpty())
			return "";
		else
			_url = newArgument;
		
		newArgument = Common.getURLFromSelectBox(_url, "list_type_id", geneListTypeListBox, jsonErrorDialog, "0");
		if (newArgument.isEmpty())
			return "";
		else
			_url = newArgument;
		
		if (topNumberTxtBox.getText().isEmpty())
		{
			jsonErrorDialog.displayError("Form Error", "Top number is not given.");
			return "";
		}
		else
			_url += "&" + "max_rank=" + topNumberTxtBox.getText();
		
		return _url;
	}

}
