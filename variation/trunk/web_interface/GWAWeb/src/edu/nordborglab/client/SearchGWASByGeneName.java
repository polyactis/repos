/**
 * 2009-5-26
 * 		class to search GWAS results by gene name
 */
package edu.nordborglab.client;

import com.google.gwt.http.client.URL;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.SuggestBox;
import com.google.gwt.user.client.ui.SuggestionEvent;
import com.google.gwt.user.client.ui.SuggestionHandler;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;

public class SearchGWASByGeneName extends Sink implements ClickListener{
	
	private CustomSuggestOracle oracle;
	//String[] words = constants.cwSuggestBoxWords();
	//for (int i = 0; i < words.length; ++i) {
	//	oracle.add(words[i]);
	//}

	// Create the suggest box
	private SuggestBox suggestBox;
	private Button suggestButton;
	//suggestBox.ensureDebugId("cwSuggestBox");
	private HorizontalPanel suggestPanel;
	public ListBox typeListBox;
	public TextBox minScoreTxtBox;
	public TextBox maxRankTxtBox;
	private HTML statusReport;
	
	private VerticalPanel formVPanel;
	private DisclosurePanel formPanel;
	
	private AbstractDataTable geneSummaryDataTable;
	private CustomVisualizationTable geneSummaryTable;
	private AbstractDataTable associationDataTable;
	private WrapVisualizationTable associationTable;
	private DisclosurePanel geneSummaryPanel;
	
	private CustomVerticalPanel vPanel;
	
	//private String dataUrl = "http://spreadsheets.google.com/tq?key=prll1aQH05yQqp_DKPP9TNg&pub=1";
	//private Query query = Query.create(dataUrl);
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	
	private int site_public;
	private int snp_gene_association_type_id;
	private String snpGeneAssociationTypeURL;
	private String geneNameAutoCompleteURL;
	private String searchURL;
	
	private static final String SUGGEST_BUTTON_DEFAULT_TEXT = "Search";
	private static final String SUGGEST_BUTTON_WAITING_TEXT = "Waiting...";

	
	public SearchGWASByGeneName(AccessionConstants constants, DisplayJSONObject jsonErrorDialog)
	{
		//super(constants);
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		site_public = get_site_public();
		snp_gene_association_type_id = get_snp_gene_association_type_id();
		snpGeneAssociationTypeURL = getsnpGeneAssociationTypeURL();
		geneNameAutoCompleteURL = getgeneNameAutoCompleteURL();
		searchURL = getsearchURL();
		
		oracle = new CustomSuggestOracle(geneNameAutoCompleteURL + "?namelike=");
		
		vPanel = new CustomVerticalPanel(constants, jsonErrorDialog, constants.SearchGWASByGeneNameHelpID());
		
		formVPanel = new VerticalPanel();
		typeListBox = new ListBox();
		HorizontalPanel hpanel1 = new HorizontalPanel(); 
		hpanel1.add(new HTML("Type: "));
		hpanel1.add(typeListBox);
		hpanel1.setSpacing(5);
		if (site_public==1)	// hide it if the site is public.
			hpanel1.setVisible(false);
		else
			Common.fillSelectBox(snpGeneAssociationTypeURL, typeListBox, jsonErrorDialog);
		formVPanel.add(hpanel1);
		
		minScoreTxtBox = new TextBox();
		minScoreTxtBox.setWidth("50px");
		maxRankTxtBox = new TextBox();
		maxRankTxtBox.setWidth("50px");
		maxRankTxtBox.setText("50");
		HorizontalPanel hpanel4 = new HorizontalPanel();
		hpanel4.add(new HTML("min score: "));
		hpanel4.add(minScoreTxtBox);
		hpanel4.add(new HTML("max rank: "));
		hpanel4.add(maxRankTxtBox);
		hpanel4.setSpacing(5);
		formVPanel.add(hpanel4);
		
		suggestBox = new SuggestBox(oracle);
		//Label lbl = new Label(constants.cwAccessionByNameLabel());
		//suggestBox.addChangeListener(new SuggestBoxChangeListener());
		suggestBox.addEventHandler(new SuggestBoxSuggestionHandler());
		suggestBox.setTitle(constants.SearchGWASByGeneNameSuggestBoxTitle());
		
		suggestButton = new Button();
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
		suggestButton.addClickListener(this);
		
		statusReport = new HTML(this.constants.LoadingText());
		statusReport.setVisible(false);
		
		suggestPanel = new HorizontalPanel();
		suggestPanel.add(new HTML("Gene ID/Name: "));
		suggestPanel.add(suggestBox);
		suggestPanel.add(suggestButton);
		suggestPanel.add(statusReport);
		suggestPanel.setSpacing(5);
		
		formVPanel.add(suggestPanel);
		
		formPanel = new DisclosurePanel("Form");
		formPanel.setAnimationEnabled(true);
		formPanel.setContent(formVPanel);
		formPanel.setOpen(true);
		
		
		geneSummaryTable = new CustomVisualizationTable();
		geneSummaryPanel = new DisclosurePanel("Gene Summary Info");
		geneSummaryPanel.setAnimationEnabled(true);
		geneSummaryPanel.add(geneSummaryTable);
		geneSummaryPanel.setOpen(true);
		

		associationTable = new WrapVisualizationTable();
		
		vPanel.add(formPanel);
		//vPanel.add(geneSummaryTable);
		vPanel.add(geneSummaryPanel);
		vPanel.add(associationTable);
		// All composites must call initWidget() in their constructors.
		initWidget(vPanel);
		
	}
	
	public native int get_site_public()/*-{ return $wnd.site_public; }-*/;
	
	public native int get_snp_gene_association_type_id()/*-{ return $wnd.snp_gene_association_type_id; }-*/;
	
	public native String getsnpGeneAssociationTypeURL()/*-{ return $wnd.snpGeneAssociationTypeURL; }-*/;
	
	public native String getgeneNameAutoCompleteURL()/*-{ return $wnd.geneNameAutoCompleteURL; }-*/;
	
	public native String getsearchURL()/*-{ return $wnd.searchURL; }-*/;
	
	public void onClick(Widget sender) {
		doFetchURL();
	}
	
	private class SuggestBoxSuggestionHandler implements SuggestionHandler {
		public void onSuggestionSelected(SuggestionEvent event){
			doFetchURL();
		}
	}
	
	@Override
	public String getName() {
		return "SearchGWASByGeneName";
	}
	
	private String getTypeID()
	{
		if (site_public==1)
		{
			Integer typeID = (Integer) snp_gene_association_type_id;
			return typeID.toString();
		}
		else
		{
			return typeListBox.getValue(typeListBox.getSelectedIndex());
		}
	}

	private void doFetchURL() {
		setIntoWaitState();
		suggestButton.setText(SUGGEST_BUTTON_WAITING_TEXT);
		String url;
		String _url = "";
		_url += "typeID=" + getTypeID();
		_url += "&minScore=" + minScoreTxtBox.getText();
		_url += "&maxRank=" + maxRankTxtBox.getText();
		_url += "&geneName=" + suggestBox.getText();
		
		url = URL.encode(searchURL + "?"+ _url + "&mode=1");
		Common.fillInTable(url, geneSummaryTable, jsonErrorDialog, statusReport, geneSummaryDataTable);
		
		url = URL.encode(searchURL + "?"+ _url + "&mode=2");
		Common.fillInTable(url, associationTable, jsonErrorDialog, statusReport, associationDataTable);
		resetTitle();
	}
	
	public void setIntoWaitState()
	{
		suggestButton.setText(SUGGEST_BUTTON_WAITING_TEXT);
		suggestButton.setEnabled(false);
	}
	public void resetTitle()
	{
		suggestButton.setText(SUGGEST_BUTTON_DEFAULT_TEXT);
		suggestButton.setEnabled(true);
	}
}
