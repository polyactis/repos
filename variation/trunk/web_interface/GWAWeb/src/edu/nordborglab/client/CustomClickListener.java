package edu.nordborglab.client;
import com.google.gwt.user.client.Event;

import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Widget;

public interface CustomClickListener extends ClickListener{
	abstract void onClick(Widget sender, Event evt);

}
