
//var field_ls = ['chromosome', 'position', 'call_method_id', 'phenotype_method_id', 'analysis_method_id', 'score'];

function loadTable(baseURL, elemName, columnIndexArray, sndElem) {
	var displayTableCallBack = {
			success: function(o) {
				var replace_elem = document.getElementById(elemName);
				var status_elem = document.getElementById(elemName+"Status");
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
				status_elem.style.display = 'none';
				
				if (!YAHOO.lang.isUndefined(sndElem)){
					var replace_elem = document.getElementById(sndElem);
					var status_elem = document.getElementById(sndElem+"Status");
					//var jsView = new google.visualization.DataView(jsData);		
					//jsView.setColumns(columnIndexArray);
					var motionChartDiv = document.getElementById(sndElem);
					var motionChart = new google.visualization.MotionChart(motionChartDiv);
					var options = {};
					options['state'] = '{"iconType":"BAR","showTrails":false,"xAxisOption":"6","colorOption":"7","yAxisOption":"6","playDuration":15,"orderedByX":true,"orderedByY":false,"sizeOption":"_UNISIZE","nonSelectedAlpha":0.4,"stateVersion":3,"dimensions":{"iconDimensions":["dim0"]},"yLambda":1,"yZoomedIn":false};';
					//{"iconType":"BAR","xZoomedDataMin":null,"yZoomedDataMax":null,"xZoomedIn":false,"iconKeySettings":[],"showTrails":false,"xAxisOption":"6","colorOption":"7","yAxisOption":"6","playDuration":15,"xZoomedDataMax":null,"orderedByX":true,"duration":{"multiplier":1},"xLambda":1,"orderedByY":false,"sizeOption":"_UNISIZE","yZoomedDataMin":null,"nonSelectedAlpha":0.4,"stateVersion":3,"dimensions":{"iconDimensions":["dim0"]},"yLambda":1,"yZoomedIn":false};
					options['width'] = 1000;
					options['height'] = 700;
					motionChart.draw(jsData, options);
					motionChartDiv.style.display = 'block';
					status_elem.style.display = 'none';
				}
			},
			failure: function(o) {
				alert("Fetching data for "+elemName+" failed.");
			}
	}
	var transaction = YAHOO.util.Connect.asyncRequest('GET', baseURL, displayTableCallBack, null);
}