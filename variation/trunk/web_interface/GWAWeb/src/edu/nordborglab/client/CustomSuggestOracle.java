package edu.nordborglab.client;

import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.URL;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.SuggestOracle;

public class CustomSuggestOracle extends SuggestOracle {
	private String suggestURL;
	
	public CustomSuggestOracle() {
		super();
		this.suggestURL = "";
	}
	
	public CustomSuggestOracle(String suggestURL) {
		super();
		this.suggestURL = suggestURL;
	}
	
	public boolean isDisplayStringHTML() {
		return false;
	}

	public void requestSuggestions(SuggestOracle.Request request, SuggestOracle.Callback callback) {
		String url = URL.encode(this.suggestURL+ request.getQuery());
		RequestBuilder builder = new RequestBuilder(RequestBuilder.GET, url);
		try {
			builder.sendRequest(null, new CustomSuggestOracleCallback(request, callback));
		} catch (RequestException e) {
			Window.alert("Couldn't retrieve JSON: " + e);
		}
	}
}