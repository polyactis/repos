/*
 * 2009-4-26 This custom click listener adds Event as an argument to the onClick() besides the Widget which sends the event. 
 */
package edu.nordborglab.client;
import com.google.gwt.user.client.Event;

import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Widget;

public interface CustomClickListener extends ClickListener{
	abstract void onClick(Widget sender, Event evt);

}
