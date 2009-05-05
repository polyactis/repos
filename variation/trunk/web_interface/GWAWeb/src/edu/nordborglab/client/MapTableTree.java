package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.DisclosurePanel;

import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Table.Options.Policy;

import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.AbstractVisualization.VisualizationFactory;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;

import com.google.gwt.core.client.JsArray;


public class MapTableTree extends VerticalPanel{
	private MapWithPhenotype mapWidget;

	private HTML accessionReport = new HTML();

	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;

	private WrapVisualizationTable accessionTable;
	
	/*
	private TextBox entriesPerPageBox;
	private Table accessionTable = new Table();
	private VerticalPanel vpanel = new VerticalPanel();
	private Button pageRefreshButton;
	private DisclosurePanel tablePanel;
	*/
	private AbstractDataTable dataTable;
	private DisclosurePanel mapPanel;
	
	class TableSelectionHandler extends SelectHandler {
		private final Selectable viz;
		private final Selectable targetWidget;

		TableSelectionHandler(Selectable viz, Selectable targetWidget) {
			this.viz = viz;
			this.targetWidget = targetWidget;
		}

		@Override
		public void onSelect(SelectEvent event) {
			JsArray<Selection> s = getSelections();
			targetWidget.setSelections(s);
		}

		private JsArray<Selection> getSelections() {
			return viz.getSelections();
		}
	}

	public MapTableTree(AccessionConstants constants, DisplayJSONObject jsonErrorDialog)
	{
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		accessionTable = new WrapVisualizationTable();
		
		mapWidget = new MapWithPhenotype(constants, jsonErrorDialog);
		mapWidget.addSelectHandler(new TableSelectionHandler(mapWidget, accessionTable));
		mapPanel = new DisclosurePanel("Map");
		mapPanel.setAnimationEnabled(true);
		mapPanel.setContent(mapWidget);
		//mapPanel.setHeader(accessionReport);
		mapPanel.setOpen(true);
		
		//TreeItem treeItem = this.addItem(accessionReport);
		//treeItem.addItem(mapWidget);
		//treeItem.setStyleName("JSON-JSONResponseObject");
		//treeItem.setState(true);

		
		
		accessionTable.addSelectHandler(new TableSelectionHandler(accessionTable, mapWidget));
		
		/*
		TreeItem treeItemTable = this.addItem("Table");
		treeItemTable.addItem(vpanel);
		treeItemTable.setState(true);
		*/
		// make the whole tree visible
		//this.setVisible(true);
		//this.setFocus(false);
		this.add(mapPanel);
		this.add(accessionTable);
	}
	
	public void populateData(AbstractDataTable dataTable)
	{
		this.dataTable = dataTable;
		//accessionReport.setHTML("Found <b>" + dataTable.getNumberOfRows() + "</b> Accessions.");
		mapPanel.getHeaderTextAccessor().setText("Found " + dataTable.getNumberOfRows() + " Accessions.");
		accessionTable.fillInTable(dataTable);
		mapWidget.addMarkers(dataTable);
	}
	
	public void setMapPanelHeaderText(String caption)
	{
		mapPanel.getHeaderTextAccessor().setText(caption);
	}
	
	public void setTablePanelHeaderText(String caption)
	{
		accessionTable.setTablePanelHeaderText(caption);
	}
	
	public void resetTablePanelHeaderText()
	{
		accessionTable.resetTablePanelHeaderText();
	}
	public void resetSize()
	{
		mapWidget.resetMapSize();
	}
}
