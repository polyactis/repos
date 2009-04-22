/*
 * 2009-4-14
 * copied from MapsDemo.java from the HelloMaps sample in gwt-maps
 */
package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.HTML;

/**
 * All HelloMaps demos extend this class.
 */
public abstract class Sink extends Composite {

	public HTML getDescriptionHTML() {
		return new HTML("<p><i>Description not provided.</i></p>\n"
				+ "<p>(Add an implementation of <code>getDescriptionHTML()</code> "
				+ "for this demo)</p>");
	}
	public abstract String getName();


	/**
	 * Method that gets called by the main demo when this demo is now active on
	 * the screen.
	 */
	public void onShow() {
	}
	
	public void resetSize()
	{
		
	}
}
