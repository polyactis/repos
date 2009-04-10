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

import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.visualizations.Table;
import com.google.gwt.visualization.client.visualizations.Table.Options.Policy;

import com.google.gwt.core.client.JsArray;
import com.google.gwt.visualization.client.AbstractVisualization;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.AbstractVisualization.VisualizationFactory;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.events.SelectHandler.SelectEvent;

public class MapTableTree extends Tree{
	private MapWithPhenotype mapWidget;

	private HTML accessionReport = new HTML();

	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;

	private TextBox entriesPerPageBox;
	private Table accessionTable = new Table();
	private VerticalPanel vpanel = new VerticalPanel();
	
	private AbstractDataTable dataTable;
	private Button pageRefreshButton;
	
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
		
		mapWidget = new MapWithPhenotype(constants, jsonErrorDialog);
		mapWidget.addSelectHandler(new TableSelectionHandler(mapWidget, accessionTable));
		TreeItem treeItem = this.addItem(accessionReport);
		treeItem.addItem(mapWidget);
		//treeItem.setStyleName("JSON-JSONResponseObject");
		treeItem.setState(true);


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
		
		vpanel.add(accessionTable);
		accessionTable.addSelectHandler(new TableSelectionHandler(accessionTable, mapWidget));
		vpanel.add(hpanel);
		
		TreeItem treeItemTable = this.addItem("Table");
		treeItemTable.addItem(vpanel);
		treeItemTable.setState(true);
		// make the whole tree visible
		this.setVisible(true);
		//this.setFocus(false);

	}
	
	public void fillInTable(AbstractDataTable dataTable)
	{
		if (dataTable!=null)
		{
			Table.Options options = Table.Options.create();
			options.setShowRowNumber(true);
			options.setAllowHtml(true);
			options.setPage(Policy.ENABLE);
			int page_size = Integer.parseInt(entriesPerPageBox.getText());
			options.setPageSize(page_size);
			accessionTable.draw(dataTable, options);
		}
	}
	
	public void populateData(AbstractDataTable dataTable)
	{
		this.dataTable = dataTable;
		accessionReport.setHTML("Found <b>" + dataTable.getNumberOfRows() + "</b> Accessions.");
		fillInTable(dataTable);
		pageRefreshButton.setVisible(true);
		mapWidget.addMarkers(dataTable);		
	}
}
