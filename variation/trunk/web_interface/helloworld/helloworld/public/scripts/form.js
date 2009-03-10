//2009-1-30 form-related javascript functions

function SelectChange(url, field_ls, replace){
	var callback = {
		success: function(o) {
			// YAHOO.util.Dom.get(replace).innerHTML = o.responseText;
			var parsed_options = YAHOO.lang.JSON.parse(o.responseText);
			var replace_elem = document.getElementById(replace);
			// Remove current options
			while(replace_elem.hasChildNodes() === true)
			{
				replace_elem.removeChild(replace_elem.childNodes[0]);
			}
			// Add new options
			for (var i=0; i<parsed_options.options.length; i++) {
				var new_option = document.createElement('option');
				new_option.text = parsed_options.options[i].id;
				new_option.value =  parsed_options.options[i].value;
				replace_elem.appendChild(new_option);
			}
			if (!YAHOO.lang.isUndefined(parsed_options.results_id)){
				var elem = document.getElementById('results_id');
				elem.value = parsed_options.results_id;
			}
			if (!YAHOO.lang.isUndefined(parsed_options.call_method_id)){
				var elem = document.getElementById('call_method_id');
				elem.value = parsed_options.call_method_id;
			}
		},
		failure: function(o) {
			alert("No data fetched under current selection.");
			var replace_elem = document.getElementById(replace);
			// Remove current options
			while(replace_elem.hasChildNodes() === true)
			{
				replace_elem.removeChild(replace_elem.childNodes[0]);
			}
         }
     }
	var url_query = "";
	for (var i = 0; i < field_ls.length; i += 1) {
		field = field_ls[i];
		if (i!=0){
			url_query = url_query + '&';
		}
		url_query = url_query + field + '=' + YAHOO.util.Dom.get(field).value;
		//url = url +'?'+field+'='+YAHOO.util.Dom.get(field).value;
	}
	url = url +'?' + url_query;
	var transaction = YAHOO.util.Connect.asyncRequest('GET', url, callback, null);
}

function removeChildNodes(domElem){
	while(domElem.hasChildNodes() === true)
	{
		domElem.removeChild(domElem.childNodes[0]);
	}
}