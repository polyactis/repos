package edu.nordborglab.client;

import com.google.gwt.visualization.client.DataTable;

public class Common {
	public static native DataTable asDataTable(String json) /*-{
		dataTable = new $wnd.google.visualization.DataTable(eval("("+json+")"), 0.5); 
		return dataTable;
	}-*/;
}
