package edu.nordborglab.client;

import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.URL;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.SuggestOracle;

public class AccessionSuggestOracle extends SuggestOracle {
	private AccessionConstants constants;

	public AccessionSuggestOracle() {
		super();
	}
	
	public AccessionSuggestOracle(AccessionConstants constants) {
		super();
		this.constants = constants;
	}
	public void setConstants(AccessionConstants constants){
		this.constants = constants;
	}
	
	public boolean isDisplayStringHTML() {
		return false;
	}

	public void requestSuggestions(SuggestOracle.Request request, SuggestOracle.Callback callback) {
		String url = URL.encode(constants.AccessionSuggestOracleURL() + "?namelike=" + request.getQuery());
		RequestBuilder builder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			builder.sendRequest(null, new AccessionSuggestOracleCallback(request, callback));
		} catch (RequestException e) {
			Window.alert("Couldn't retrieve JSON: " + e);
		}
	}
}