package edu.nordborglab.client;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.user.client.ui.DisclosurePanel;

import com.google.gwt.user.client.ui.ClickListener;
import com.google.gwt.user.client.ui.LoadListener;
import com.google.gwt.user.client.ui.ChangeListener;

public class HaplotypeSingleView extends Sink {
	public ListBox callMethodListBox;
	public ListBox phenotypeMethodListBox;
	public ListBox geneListTypeListBox;
	
	public TextBox chrTxtBox;
	public TextBox startTxtBox;
	public TextBox stopTxtBox;
	public ListBox whichPCBox;
	public Button submitButton;
	
	private VerticalPanel formVPanel;
	public Image image;
	private DisclosurePanel formPanel;
	private DisclosurePanel imagePanel;
	private VerticalPanel vPanel;
	
	private DisplayJSONObject jsonErrorDialog;
	private AccessionConstants constants;
	private static final String SUBMIT_BUTTON_DEFAULT_TEXT = "Submit";
	private static final String SUBMIT_BUTTON_WAITING_TEXT = "Waiting...";
	
	public HaplotypeSingleView(AccessionConstants constants, DisplayJSONObject jsonErrorDialog){
		this.constants = constants;
		this.jsonErrorDialog = jsonErrorDialog;
		
		vPanel = new VerticalPanel();
		formVPanel = new VerticalPanel();
		
		callMethodListBox = new ListBox();
		callMethodListBox.getElement().getId();
		/*
		callMethodListBox.addChangeListener(new ChangeListener() {
			public void onChange(Widget sender) {
				//refreshMap(map, phenotypeSelectBox.getSelectedIndex(), displayOptionSelectBox.getSelectedIndex());
				//multiBox.ensureDebugId("cwListBox-multiBox");
			}
		});
		*/
		HorizontalPanel hpanel1 = new HorizontalPanel(); 
		hpanel1.add(new HTML(constants.callMethodLabel()));
		hpanel1.add(callMethodListBox);
		hpanel1.setSpacing(5);
		formVPanel.add(hpanel1);
		
		phenotypeMethodListBox = new ListBox();
		HorizontalPanel hpanel2 = new HorizontalPanel();
		hpanel2.add(new HTML(constants.phenotypeMethodLabel()));
		hpanel2.add(phenotypeMethodListBox);
		hpanel2.setSpacing(5);
		formVPanel.add(hpanel2);
		
		geneListTypeListBox = new ListBox();
		HorizontalPanel hpanel3 = new HorizontalPanel();
		hpanel3.add(new HTML(constants.geneListTypeLabel()));
		hpanel3.add(geneListTypeListBox);
		hpanel3.setSpacing(5);
		formVPanel.add(hpanel3);
		
		chrTxtBox = new TextBox();
		chrTxtBox.setWidth("50px");
		startTxtBox = new TextBox();
		stopTxtBox = new TextBox();
		HorizontalPanel hpanel4 = new HorizontalPanel();
		hpanel4.add(new HTML("chromosome: "));
		hpanel4.add(chrTxtBox);
		hpanel4.add(new HTML("start: "));
		hpanel4.add(startTxtBox);
		hpanel4.add(new HTML("stop: "));
		hpanel4.add(stopTxtBox);
		hpanel4.setSpacing(5);
		formVPanel.add(hpanel4);
		
		whichPCBox = new ListBox();
		HorizontalPanel hpanel5 = new HorizontalPanel();
		hpanel5.add(new HTML(constants.whichPCLabel()));
		hpanel5.add(whichPCBox);
		hpanel5.setSpacing(5);
		formVPanel.add(hpanel5);
		
		submitButton = new Button();
		submitButton.setText(SUBMIT_BUTTON_DEFAULT_TEXT);
		formVPanel.add(submitButton);
		
		formPanel = new DisclosurePanel("Form");
		formPanel.setAnimationEnabled(true);
		formPanel.setContent(formVPanel);
		formPanel.setOpen(true);
		
		image = new Image();
		image.addLoadListener(new LoadListener() {
			public void onError(Widget sender) {
				submitButton.setText("Error occurred while loading.");
			}

			public void onLoad(Widget sender) {
				// submitButton.setText(SUBMIT_BUTTON_WAITING_TEXT);
			}
		});
		image.setSize("1000px", "800px");
		image.setVisible(false);
		/* When the user clicks this button, we want to clip the image.
		submitButton.addClickListener(new ClickListener() {
			public void onClick(Widget sender) {
				image.setVisibleRect(70, 0, 47, 110);
			}
		});
		*/
		
		imagePanel = new DisclosurePanel();
		imagePanel.setAnimationEnabled(true);
		imagePanel.setContent(image);
		imagePanel.setOpen(true);
		
		vPanel.add(formPanel);
		vPanel.add(imagePanel);
		initWidget(vPanel);
	}
	
	@Override
	public String getName() {
		String name;
		if (callMethodListBox.getSelectedIndex()>0)
		{
			name = callMethodListBox.getValue(callMethodListBox.getSelectedIndex());
		}
		else
			name = "_";
		
		if (phenotypeMethodListBox.getSelectedIndex()>=0) 
		{
			String phenotypeMethodValue = phenotypeMethodListBox.getValue(phenotypeMethodListBox.getSelectedIndex());
			name += phenotypeMethodValue;
		}
		if (geneListTypeListBox.getSelectedIndex()>=0)
		{
			name += geneListTypeListBox.getValue(geneListTypeListBox.getSelectedIndex());
		}
		name += chrTxtBox.getText();
		name += startTxtBox.getText();
		name += stopTxtBox.getText();
		if (whichPCBox.getSelectedIndex()>=0)
		{
			name += whichPCBox.getValue(whichPCBox.getSelectedIndex());
		}
		return name;
	}
	
	public void resetSubmitButtonCaption() {
		submitButton.setText(SUBMIT_BUTTON_DEFAULT_TEXT);
	}
	
	public void copyListBox(ListBox oldListBox, ListBox newListBox)
	{
		for (int i =0; i<oldListBox.getItemCount(); i++)
		{
			this.callMethodListBox.addItem(oldListBox.getItemText(i), oldListBox.getValue(i));
		}
		if (oldListBox.getSelectedIndex()>0)
			newListBox.setSelectedIndex(oldListBox.getSelectedIndex());
	}
	
	public void copySetting(HaplotypeSingleView haplotypeSingleView)
	{
		copyListBox(haplotypeSingleView.callMethodListBox, this.callMethodListBox);
		copyListBox(haplotypeSingleView.phenotypeMethodListBox, this.phenotypeMethodListBox);
		copyListBox(haplotypeSingleView.geneListTypeListBox, this.geneListTypeListBox);
		copyListBox(haplotypeSingleView.whichPCBox, this.whichPCBox);
		this.chrTxtBox.setText(haplotypeSingleView.chrTxtBox.getText());
		this.startTxtBox.setText(haplotypeSingleView.startTxtBox.getText());
		this.stopTxtBox.setText(haplotypeSingleView.stopTxtBox.getText());
	}
	public void setSubmitButtonInProgress() {
		submitButton.setText(SUBMIT_BUTTON_WAITING_TEXT);
	}
}
