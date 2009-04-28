package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.FlexTable;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.Widget;

public class CustomVerticalPanel extends VerticalPanel{
	private AccessionConstants constants;
	private DisplayJSONObject jsonErrorDialog;
	
	
	private HelpDialog helpDialog;
	private String helpID;
	
	private FlexTable layout;
	private Button helpButton;
	
	CustomVerticalPanel(AccessionConstants constants, DisplayJSONObject jsonErrorDialog, String helpID)
	{
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		this.helpID = helpID;
		
		
		helpDialog = new HelpDialog(constants, jsonErrorDialog, helpID);
		helpDialog.hide();
		
		helpButton = new Button("help");
		helpButton.addClickListener(new ClickListener(){
			public void onClick(Widget sender)
			{
				helpDialog.center();
				helpDialog.show();
			}
		});
		
		layout = new FlexTable();
		layout.setWidget(0, 0, helpButton);
		this.add(layout);
		//this.setCellHorizontalAlignment(layout, HasHorizontalAlignment.ALIGN_LEFT);
	}
}
