google.load("visualization", "1", {packages:["scatterchart", 'map', 'motionchart', "table"]});


//var field_ls = ['chromosome', 'position', 'call_method_id', 'phenotype_method_id', 'analysis_method_id', 'score'];

function loadTable(baseURL, elemName, columnIndexArray, sndElem) {
	var displayTableCallBack = {
			success: function(o) {
				var replace_elem = document.getElementById(elemName);
				var table = new google.visualization.Table(replace_elem);
				var response = o.responseText;
				var jsData = new google.visualization.DataTable(eval("("+response+")"), 0.5);
				if (!YAHOO.lang.isUndefined(columnIndexArray)){
					var jsView = new google.visualization.DataView(jsData);		
					jsView.setColumns(columnIndexArray);
				}
				else{
					var jsView = jsData;
				}
				table.draw(jsView, {allowHtml: true, showRowNumber: true});
				//replace_elem.innerHTML = "<a href="+displayResultsSpace.imgQueryURL+"><img src="+displayResultsSpace.imgQueryURL+"></a>";
				replace_elem.style.display = 'block';
				
				if (!YAHOO.lang.isUndefined(sndElem)){
					var replace_elem = document.getElementById(sndElem);
					//var jsView = new google.visualization.DataView(jsData);		
					//jsView.setColumns(columnIndexArray);
					var motionChartDiv = document.getElementById(sndElem);
					var motionChart = new google.visualization.MotionChart(motionChartDiv);
					motionChart.draw(jsData, {width: 1000, height:700});
					motionChartDiv.style.display = 'block';
				}
			},
			failure: function(o) {
				alert("Fetching data for "+elemName+" failed.");
			}
	}
	var transaction = YAHOO.util.Connect.asyncRequest('GET', baseURL, displayTableCallBack, null);
}