package edu.nordborglab.client;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.user.client.Event;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.DisclosurePanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Table.Options;


public class WrapVisualizationTable extends Composite implements Selectable{
	private DisclosurePanel tablePanel;
	private Button pageRefreshButton;
	private TextBox entriesPerPageBox;
	
	private CustomVisualizationTable visualizationTable;
	private VerticalPanel vpanel;
	
	private AbstractDataTable dataTable;
	
	private String DEFAULT_TABLE_PANEL_TITLE = "Table";
	
	public WrapVisualizationTable() {
		
		entriesPerPageBox = new TextBox();
		entriesPerPageBox.setText("20");
		
		pageRefreshButton = new Button("refresh");
		pageRefreshButton.setVisible(false);	//invisible upon initilization cuz accessionTable is still empty.
		pageRefreshButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				fillInTable(dataTable);
			}
		});
		
		HorizontalPanel hpanel = new HorizontalPanel();
		hpanel.add(entriesPerPageBox);
		hpanel.add(new HTML("entries per page."));
		hpanel.add(pageRefreshButton);
		hpanel.setSpacing(5);
		//hpanel.setVisible(false);
		
		visualizationTable = new CustomVisualizationTable();
		vpanel = new VerticalPanel();
		vpanel.add(visualizationTable);
		vpanel.add(hpanel);
		
		tablePanel = new DisclosurePanel(DEFAULT_TABLE_PANEL_TITLE);
		tablePanel.setAnimationEnabled(true);
		tablePanel.setContent(vpanel);
		tablePanel.setOpen(true);
		initWidget(tablePanel);
	}
	
	public void fillInTable(AbstractDataTable dataTable)
	{
		if (dataTable!=null)
		{
			Table.Options options = Table.Options.create();
			options.setShowRowNumber(true);
			options.setAllowHtml(true);
			options.setPage(Options.Policy.ENABLE);
			int page_size = Integer.parseInt(entriesPerPageBox.getText());
			options.setPageSize(page_size);
			visualizationTable.draw(dataTable, options);
			pageRefreshButton.setVisible(true);
			DEFAULT_TABLE_PANEL_TITLE = dataTable.getNumberOfRows() + " Row(s) in Total";
			resetTablePanelHeaderText();
		}
	}
	public void draw(AbstractDataTable data){
		this.dataTable = data;
		this.visualizationTable.draw(this.dataTable);
	}
	
	public void draw(AbstractDataTable data, Options options){
		this.dataTable = data;
		this.visualizationTable.draw(this.dataTable, options);
	}
	
	public final void addSelectHandler(SelectHandler handler)
	{
		this.visualizationTable.addSelectHandler(handler);
	}
	
	public void addClickListener(CustomClickListener cListener)
	{
		this.visualizationTable.addClickListener(cListener);
	}
	
	public void setSelections(JsArray<Selection> s)
	{
		this.visualizationTable.setSelections(s);
	}
	public JsArray<Selection> getSelections()
	{
		//return ArrayHelper.toJsArray(Selection.createRowSelection(selectedRow));
		return this.visualizationTable.getSelections();
	}
	
	public void setTablePanelHeaderText(String caption)
	{
		tablePanel.getHeaderTextAccessor().setText(caption);
	}
	
	public void resetTablePanelHeaderText()
	{
		tablePanel.getHeaderTextAccessor().setText(DEFAULT_TABLE_PANEL_TITLE);
	}
}
