google.load("visualization", "1", {packages:["scatterchart", 'map', 'motionchart', "table"]});


//var field_ls = ['chromosome', 'position', 'call_method_id', 'phenotype_method_id', 'analysis_method_id', 'score'];

function loadTable(baseURL, elemName) {
	var displayTableCallBack = {
			success: function(o) {
				var replace_elem = document.getElementById(elemName);
				var table = new google.visualization.Table(replace_elem);
				var response = o.responseText;
				var snpSummaryData = new google.visualization.DataTable(eval("("+response+")"), 0.5);
				table.draw(snpSummaryData, {allowHtml: true, showRowNumber: true});
				//replace_elem.innerHTML = "<a href="+displayResultsSpace.imgQueryURL+"><img src="+displayResultsSpace.imgQueryURL+"></a>";
				replace_elem.style.display = 'block';
			},
			failure: function(o) {
				alert("Fetching data for "+elemName+" failed.");
			}
	}
	var transaction = YAHOO.util.Connect.asyncRequest('GET', baseURL, displayTableCallBack, null);
}

function loadSignificantHits(baseURL, elemName) {
	var snpSignificantHitsCallBack = {
			success: function(o) {
				var replace_elem = document.getElementById(elemName);
				var table = new google.visualization.Table(replace_elem);
				var response = o.responseText;
				//var callInfoData = response.getDataTable();
				var snpSummaryData = new google.visualization.DataTable(eval("("+response+")"), 0.5);
				table.draw(snpSummaryData, {showRowNumber: true});
				//replace_elem.innerHTML = "<a href="+displayResultsSpace.imgQueryURL+"><img src="+displayResultsSpace.imgQueryURL+"></a>";
				replace_elem.style.display = 'block';
			},
			failure: function(o) {
				alert("Fetching SNP significant hits failed.");
			}
	}
	var transaction = YAHOO.util.Connect.asyncRequest('GET', baseURL, snpSignificantHitsCallBack, null);
	
}