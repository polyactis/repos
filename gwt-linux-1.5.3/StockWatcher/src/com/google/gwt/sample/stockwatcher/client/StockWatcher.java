package com.google.gwt.sample.stockwatcher.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

import com.google.gwt.user.client.ui.FlexTable;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.TextBox;

import com.google.gwt.user.client.ui.KeyboardListenerAdapter;
import com.google.gwt.user.client.Window;

import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.HTML;

import java.util.ArrayList;

import com.google.gwt.user.client.Timer;
import com.google.gwt.user.client.Random;
import com.google.gwt.i18n.client.NumberFormat;

import com.google.gwt.i18n.client.DateTimeFormat;

import java.util.Date;


/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class StockWatcher implements EntryPoint {

	/**
	 * This is the entry point method.
	 */
	private static final int REFRESH_INTERVAL = 5000; // ms
	
	private VerticalPanel mainPanel = new VerticalPanel();

	private FlexTable stocksFlexTable = new FlexTable();

	private HorizontalPanel addPanel = new HorizontalPanel();

	private TextBox newSymbolTextBox = new TextBox();

	private Button addStockButton = new Button("Add");

	private Label lastUpdatedLabel = new Label();
	private ArrayList<String> stocks = new ArrayList<String>();
	
	public void onModuleLoad() {

		// Create table for stock data.
		stocksFlexTable.setText(0, 0, "Symbol");
		stocksFlexTable.setText(0, 1, "Price");
		stocksFlexTable.setText(0, 2, "Change");
		stocksFlexTable.setText(0, 3, "Remove");

		 // Add styles to elements in the stock list table.
		stocksFlexTable.getRowFormatter().addStyleName(0, "watchListHeader");
		stocksFlexTable.addStyleName("watchList");
		stocksFlexTable.getCellFormatter().addStyleName(0, 1,
				"watchListNumericColumn");
		stocksFlexTable.getCellFormatter().addStyleName(0, 2,
				"watchListNumericColumn");
		stocksFlexTable.getCellFormatter().addStyleName(0, 3,
				"watchListRemoveColumn");

		// Assemble Add Stock panel.
		addPanel.add(newSymbolTextBox);
		addPanel.add(addStockButton);
		addPanel.addStyleName("addPanel");

		// Assemble Main panel.
		mainPanel.add(stocksFlexTable);
		mainPanel.add(addPanel);
		mainPanel.add(lastUpdatedLabel);

		// Associate the Main panel with the HTML host page.
		RootPanel.get("stockList").add(mainPanel);
		// Move cursor focus to the input box.
		newSymbolTextBox.setFocus(true);
		// Setup timer to refresh list automatically.
		Timer refreshTimer = new Timer() {
			public void run() {
				refreshWatchList();
			}
		};
		refreshTimer.scheduleRepeating(REFRESH_INTERVAL);

		// Listen for mouse events on the Add button.
		addStockButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				addStock();
				}
			});
		
		newSymbolTextBox.addKeyboardListener(new KeyboardListenerAdapter() {
			@Override
			public void onKeyDown(Widget sender, char keyCode, int modifiers) {
				if (keyCode == KEY_ENTER) {
					addStock();
					}
				}
			});

		TabPanel tp = new TabPanel();
		tp.add(new HTML("Foo"), "foo");
		tp.add(new HTML("Bar"), "bar");
		tp.add(new HTML("Baz"), "baz");

		// Show the 'bar' tab initially.
		tp.selectTab(1);
		tp.setWidth("400px");
		tp.setAnimationEnabled(true);

		// Add it to the root panel.
		RootPanel.get("stockList").add(tp);
		
		Image img = new Image(
				"http://code.google.com/webtoolkit/logo-185x175.png");
		Button button = new Button("Click me");

		// We can add style names
		button.addStyleName("pc-template-btn");
		// or we can set an id on a specific element for styling
		img.getElement().setId("pc-template-img");

		VerticalPanel vPanel = new VerticalPanel();
		vPanel.setWidth("100%");
		vPanel.setHorizontalAlignment(VerticalPanel.ALIGN_CENTER);
		vPanel.add(img);
		vPanel.add(button);

		// Add image and button to the RootPanel
		RootPanel.get("stockList").add(vPanel);

		// Create the dialog box
		final DialogBox dialogBox = new DialogBox();
		dialogBox.setText("Welcome to GWT!");
		dialogBox.setAnimationEnabled(true);
		Button closeButton = new Button("close");
		VerticalPanel dialogVPanel = new VerticalPanel();
		dialogVPanel.setWidth("100%");
		dialogVPanel.setHorizontalAlignment(VerticalPanel.ALIGN_CENTER);
		dialogVPanel.add(closeButton);

		closeButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				dialogBox.hide();
			}
		});

		// Set the contents of the Widget
		dialogBox.setWidget(dialogVPanel);

		button.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				dialogBox.center();
				dialogBox.show();
			}
		});
	}
	
	/**
	* Add stock to FlexTable. Executed when the user clicks the
	* addStockButton or presses enter in the newSymbolTextBox.
	*/
	

	private void addStock() {
		// TODO Auto-generated method stub
		final String symbol = newSymbolTextBox.getText().toUpperCase().trim();
		newSymbolTextBox.setFocus(true);
		
		// Stock code must be between 1 and 10 chars that are numbers, letters, or dots.
		if (!symbol.matches("^[0-9a-zA-Z\\.]{1,10}$")) {
			Window.alert("'" + symbol + "' is not a valid symbol.");
			newSymbolTextBox.selectAll();
			return;
		}
		
		if (stocks.contains(symbol))
			return;
		// Add the stock to the table.
		int row = stocksFlexTable.getRowCount();
		stocks.add(symbol);
		stocksFlexTable.setText(row, 0, symbol);
		stocksFlexTable.setWidget(row, 2, new Label());
		
		stocksFlexTable.getCellFormatter().addStyleName(row, 1, "watchListNumericColumn");
		stocksFlexTable.getCellFormatter().addStyleName(row, 2, "watchListNumericColumn");
		stocksFlexTable.getCellFormatter().addStyleName(row, 3, "watchListRemoveColumn");
		
		// Add the stock to the table.
		Button removeStockButton = new Button("x");
		removeStockButton.addStyleDependentName("remove");
		
		removeStockButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				int removedIndex = stocks.indexOf(symbol);
				stocks.remove(removedIndex);
				stocksFlexTable.removeRow(removedIndex + 1);
			}
		});
		stocksFlexTable.setWidget(row, 3, removeStockButton);
		
		// Get the stock price.
		refreshWatchList();

		newSymbolTextBox.setText("");

		// TODO Don't add the stock if it's already in the table.

		// TODO Add the stock to the table.

		// TODO Add a button to remove this stock from the table.

		// TODO Get the stock price.

	}
	
	private void refreshWatchList() {
		// TODO Auto-generated method stub
		final double MAX_PRICE = 100.0; // $100.00
		final double MAX_PRICE_CHANGE = 0.02; // +/- 2%

		StockPrice[] prices = new StockPrice[stocks.size()];
		for (int i = 0; i < stocks.size(); i++) {
			double price = Random.nextDouble() * MAX_PRICE;
			double change = price * MAX_PRICE_CHANGE
					* (Random.nextDouble() * 2.0 - 1.0);

			prices[i] = new StockPrice((String) stocks.get(i), price, change);
		}

		updateTable(prices);

	}
	
	private void updateTable(StockPrice[] prices) {
		// TODO Auto-generated method stub
		for (int i = 0; i < prices.length; i++) {
			updateTable(prices[i]);
		}
		// Display timestamp showing last refresh.
		lastUpdatedLabel.setText("Last update : "
				+ DateTimeFormat.getMediumDateTimeFormat().format(new Date()));

	}
	
	private void updateTable(StockPrice price) {
		// Make sure the stock is still in the stock table.
		if (!stocks.contains(price.getSymbol())) {
			return;
		}

		int row = stocks.indexOf(price.getSymbol()) + 1;

		// Format the data in the Price and Change fields.
		String priceText = NumberFormat.getFormat("#,##0.00").format(
				price.getPrice());
		NumberFormat changeFormat = NumberFormat
				.getFormat("+#,##0.00;-#,##0.00");
		String changeText = changeFormat.format(price.getChange());
		String changePercentText = changeFormat.format(price.getChangePercent());
		
		// Populate the Price and Change fields with new data.
		stocksFlexTable.setText(row, 1, priceText);
		// stocksFlexTable.setText(row, 2, changeText + " (" + changePercentText + "%)");
		Label changeWidget = (Label)stocksFlexTable.getWidget(row, 2);
		changeWidget.setText(changeText + " (" + changePercentText + "%)");
		// Change the color of text in the Change field based on its value.
		String changeStyleName = "noChange";
		if (price.getChangePercent() < -0.1f) {
			changeStyleName = "negativeChange";
		} else if (price.getChangePercent() > 0.1f) {
			changeStyleName = "positiveChange";
		}

		changeWidget.setStyleName(changeStyleName);
		
	} 
	
}
