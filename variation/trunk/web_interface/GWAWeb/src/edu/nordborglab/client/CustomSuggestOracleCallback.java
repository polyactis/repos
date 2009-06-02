package edu.nordborglab.client;

import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.Response;
import com.google.gwt.json.client.JSONArray;
import com.google.gwt.json.client.JSONException;
import com.google.gwt.json.client.JSONObject;
import com.google.gwt.json.client.JSONParser;
import com.google.gwt.json.client.JSONValue;
import com.google.gwt.user.client.ui.SuggestOracle;
import com.google.gwt.user.client.Window;

import java.util.ArrayList;

public class CustomSuggestOracleCallback implements RequestCallback {
	private SuggestOracle.Callback callback;
	private SuggestOracle.Request suggestRequest;
	
	public CustomSuggestOracleCallback(SuggestOracle.Request suggestRequest, SuggestOracle.Callback callback) {
		this.callback = callback;
		this.suggestRequest = suggestRequest;
	}

	public void onError(Request request, Throwable exception) {
		Window.alert("Couldn't retrieve JSON: " + exception);
	}

	public void onResponseReceived(Request request, Response response) {
		if (response.getStatusCode() == 200) {
			try {
				JSONValue jsonValue = JSONParser.parse(response.getText());
				JSONObject jsonObject = jsonValue.isObject();
				
				JSONArray array = jsonObject.get("result").isArray();	//jsonValue.isArray();
				//Window.alert("json gets: "+ array);
				ArrayList<MovieSuggestion> suggestions = new ArrayList<MovieSuggestion>();
				for (int i = 0; i < array.size(); i++) {
					//JSONObject object = array.get(i).isObject();
					//String movieName = object.isString().stringValue();
					String movieName = array.get(i).isString().stringValue();
					suggestions.add(new MovieSuggestion(movieName));
				}
				callback.onSuggestionsReady(suggestRequest, new SuggestOracle.Response(suggestions));
			} catch (JSONException e) {
				Window.alert("Failed to parse JSON response: " + e);
			}
		} else {
			Window.alert("Couldn't retrieve JSON: " + response.getStatusText());
		}
	}

	private class MovieSuggestion implements SuggestOracle.Suggestion {
		private String movieName;

		public MovieSuggestion(String movieName) {
			this.movieName = movieName;
		}

		public String getDisplayString() {
			return movieName;
		}

		public String getReplacementString() {
			return movieName;
		}
	}
}

