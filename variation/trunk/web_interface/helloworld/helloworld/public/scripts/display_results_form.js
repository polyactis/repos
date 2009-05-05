var displayResultsSpace = {};
//02/17/09 initial values
displayResultsSpace.call_method_id = null;
displayResultsSpace.phenotype_method_id = null;
displayResultsSpace.analysis_method_id = null;
displayResultsSpace.static_plot_name_arr = [['getPhenotypeHistImage', 'hist_thumb','hist_plot'],
		['getPhenotypeHistImage', 'hist_log_thumb','hist_log_plot'],
		['getCallPhenotypeQQImage', 'qq_thumb', 'qq_plot'],
		['getCallPhenotypeQQImage', 'qq_log_thumb', 'qq_log_plot']]

google.load("visualization", "1", {packages:["scatterchart", 'map', 'motionchart']});

// The selection handler.
// Loop through all items in the selection and concatenate
// a single message from all of them.

function generate_handler(chr){
	return function(e) {
		this.chr = chr;
		//var elTarget = YAHOO.util.Event.getTarget(e);		//infer dom element from regular event
		//var chr = elTarget.id[elTarget.id.length-1];
		//alert(arguments.length);
		//alert(arguments[0]);
		//alert(this.chr);
		var selection = chr2chart[this.chr].getSelection();
		if (selection.length>0){
			var item = selection[0];
			var position = chr2data[this.chr].getValue(item.row, 0);
			var score = chr2data[this.chr].getValue(item.row, 1);
			//alert('You selected chr='+this.chr +" pos="+position);
			var start_pos = position-50000;
			var stop_pos = position+50000;
			var track_id = displayResultsSpace.field_value_ls.join("_");
			//window.open(displayResultsSpace.GBrowseURLJS.format(start_pos, stop_pos, this.chr)+track_id+"-"+track_id+"_SNP");
			window.open("/SNP/?chromosome="+this.chr+"&position="+position+"&call_method_id="+YAHOO.util.Dom.get("call_method_id").value+"&phenotype_method_id="+YAHOO.util.Dom.get("phenotype_method_id").value+"&analysis_method_id="+YAHOO.util.Dom.get("analysis_method_id").value+"&score="+score);
			//chr2chart[this.chr].setSelection([{row:null, column:null}]);
		}
	};
}

//2009-1-30 construct a handler for each chromosome, each handler is aware of the chromosome it corresponds to.
//	used later to addListener()
var no_of_chrs = 5;
chr2handler = {};
colors = ["blue", "green", "red", "cyan", "purple"];
chr2color = {};
for (var i=0; i<no_of_chrs; i++) {
	var chr = i+1;
	chr2handler[chr] = generate_handler(chr);
	chr2color[chr] = colors[i];
	}

//2009-1-29 establish the click event handler to each chromosome region. The selectHandler for google visualization can't distinguish which visualization gets clicked.  	
//YAHOO.util.Event.addListener("display_general", "click", clickHandler);

function drawGWChart(o)
{
	var chr2div = {};
	chr2chart = {};
	chr2data = {};
	var no_of_chrs = 5;
	var server_data = YAHOO.lang.JSON.parse(o.responseText);
	var server_chr2data = server_data["chr2data"];
	var max_value = server_data["max_value"];
	var chr2length = server_data["chr2length"];
	var max_length = server_data["max_length"];
	//var server_data = o.responseText;
	var titleY = '-log Pvalue';
	for (var i=0; i<no_of_chrs; i++) {
		var chr = i+1;
		chr2div[chr] = YAHOO.util.Dom.get('display_chr'+chr);
		//if (remove_old_plots==="on")
		//{
		removeChildNodes(chr2div[chr]);
		chr2chart[chr] = new google.visualization.ScatterChart(chr2div[chr]);
		chr2data[chr] = new google.visualization.DataTable( eval("("+server_chr2data[chr]+")"), 0.5);
		//width = 1000*chr2length[chr]/max_length
		width = 1000
		chr2chart[chr].draw(chr2data[chr], {width: width, height: 200, titleX: 'Chr'+chr, titleY: titleY, 
			legend: 'none', pointSize: 3, colors:[chr2color[chr],], max: max_value});
		chr2div[chr].style.display = 'block';
		google.visualization.events.addListener(chr2chart[chr], 'select', chr2handler[chr]);	//2009-1-29 use the clickHandler to each chr division
	}
	//<!--
	//var server_data = YAHOO.lang.JSON.parse(o.responseText);
	//for (var i=0; i<server_data.gwr.length; i++) {
	//	var chr = server_data.gwr[i].chromosome;
	//	var data_index = chr2data[chr].addRow();
	//	chr2data[chr].setValue(data_index, 0, server_data.gwr[i].position);
	//	chr2data[chr].setValue(data_index, 1, server_data.gwr[i].value);
	//}
	//-->
}

function showCallInfoData(o)
{
	//if (response.isError()) {
	//	alert('Error in call info query')
	//}
	//var response = YAHOO.lang.JSON.parse(o.responseText);
	var response = o.responseText;
	//var callInfoData = response.getDataTable();
	var callInfoData = new google.visualization.DataTable(eval("("+response+")"), 0.5);
	var strain_pca_div = document.getElementById('strain_pca_div');
	removeChildNodes(strain_pca_div);
	var strain_pca_chart = new google.visualization.ScatterChart(strain_pca_div);
	var pcaView = new google.visualization.DataView(callInfoData);
	pcaView.setColumns([7,6]);
	//alert("pcaView with "+pcaView.getNumberOfRows()+" rows and "+pcaView.getNumberOfColumns()+" cols.");
	strain_pca_chart.draw(pcaView, {width: 500, height: 300, titleX: 'PC2', titleY: 'PC1', legend: 'none', pointSize: 2});
	strain_pca_div.style.display = 'block';
	
	var geoView = new google.visualization.DataView(callInfoData);		
	geoView.setColumns([3,4,2]);
	var strain_map_div = document.getElementById('strain_map_div');
	var map = new google.visualization.Map(strain_map_div);
	map.draw(geoView, {showTip:true, enableScrollWheel:true});
	strain_map_div.style.display = 'block';
	
	// Set a 'select' event listener for the strain_pca_chart.
	// When the strain_pca_chart is selected,
	// we set the selection on the map.
	google.visualization.events.addListener(strain_pca_chart, 'select',
		function() {
			var selection = strain_pca_chart.getSelection();	// 2/3/09 scatter chart selection has non-null row&col
			for (var i=0; i<selection.length; i++){
				selection[i].column=null;	// 2/3/09 set col to null to be compatible with Map selection
				var item = selection[i];
				var label = callInfoData.getValue(item.row, 2);
				var pca_select_div = YAHOO.util.Dom.get('pca_select_div');
				pca_select_div.innerHTML = label;
				pca_select_div.style.display = 'block';
			}
			map.setSelection(selection);
		});

	// Set a 'select' event listener for the map.
	// When the map is selected,
	// we set the selection on the strain_pca_chart.
	google.visualization.events.addListener(map, 'select',
		function() {
			var selection = map.getSelection();	// 2/3/09 map selection row only 
			for (var i=0; i<selection.length; i++){
				selection[i].column=1;	// 2/3/09 set the column index to be compatible with strain_pca_chart selection
			}
			strain_pca_chart.setSelection(selection);
	});
	
	var motionView = new google.visualization.DataView(callInfoData);		
	motionView.setColumns([2,0,4,3,8,5,6,7]);
	var strain_motion_chart_div = document.getElementById('strain_motion_chart_div');
	var strain_motion_chart = new google.visualization.MotionChart(strain_motion_chart_div);
	var options = {};
	options['state'] = '{"time":"notime","iconType":"BUBBLE","xZoomedDataMin":null,"yZoomedDataMax":null,"xZoomedIn":false,"iconKeySettings":[],"showTrails":true,"xAxisOption":"2","colorOption":"4","yAxisOption":"3","playDuration":15,"xZoomedDataMax":null,"orderedByX":false,"duration":{"multiplier":1,"timeUnit":"none"},"xLambda":1,"orderedByY":false,"sizeOption":"_UNISIZE","yZoomedDataMin":null,"nonSelectedAlpha":0.4,"stateVersion":3,"dimensions":{"iconDimensions":["dim0"]},"yLambda":1,"yZoomedIn":false};';
	//options['state'] = "{'colorOption':4,'sizeOption':'_UNISIZE'};";
	options['width'] = 1000;
	options['height'] = 700;
	strain_motion_chart.draw(motionView, options);
	strain_motion_chart_div.style.display = 'block';

}

function handleResponse(urlArray, field_ls) {
	GWABaseURL=urlArray[0];
	CallInfoBaseURL=urlArray[1];
	PhenotypeHistImageBaseURL=urlArray[2];
	CallPhenotypeQQImageBaseURL=urlArray[3];
	displayResultsSpace.GBrowseURLJS=urlArray[4];
	var callback = {
			success: function(o) {
				//google.setOnLoadCallback(drawGWChart(o));
				drawGWChart(o);
			},
			failure: function(o) {
				alert("Failed to retrieve GenomeWideAssociation information.");
			}
	}
	var url_query = "";
	displayResultsSpace.field_value_ls = new Array();
	for (var i = 0; i < field_ls.length; i += 1) {
		field = field_ls[i];
		if (i!=0){
			url_query = url_query + '&';
		}
		url_query = url_query + field + '=' + YAHOO.util.Dom.get(field).value;
		displayResultsSpace.field_value_ls.push(YAHOO.util.Dom.get(field).value);
		//url = url +'?'+field+'='+YAHOO.util.Dom.get(field).value;
	}
	GWABaseURL = GWABaseURL +'?' + url_query;
	
	var transaction = YAHOO.util.Connect.asyncRequest('GET', GWABaseURL, callback, null);
	
	var display_general = YAHOO.util.Dom.get('display_general');
	display_general.innerHTML = displayResultsSpace.field_value_ls.join(" ");
	display_general.style.display = 'block';
	
	displayResultsSpace.callInfoQueryURL = CallInfoBaseURL+'?'+url_query;
					//2009-2-3 same query, but analysis_method_id is not used.
	var callInfoCallBack = {
			success: function(o) {
				showCallInfoData(o);
			},
			failure: function(o) {
				alert("Failed to retrieve Strain/Phenotype/Geo information.");
			}
	}
	var transaction = YAHOO.util.Connect.asyncRequest('GET', displayResultsSpace.callInfoQueryURL, callInfoCallBack, null);
	
	//var call_info_query = new google.visualization.Query(displayResultsSpace.callInfoQueryURL);	//2009-2-3 google query is slower than yahoo asyncRequest
	//call_info_query.send(showCallInfoData);
	for (var i = 0; i < displayResultsSpace.static_plot_name_arr.length; i+=1){
		var action_name = displayResultsSpace.static_plot_name_arr[i][0];
		var img_type = displayResultsSpace.static_plot_name_arr[i][1];
		var fullImg = displayResultsSpace.static_plot_name_arr[i][2];
		var div_name = 'display_'+img_type;
		//02/17/09 test if necessary to update plots			
		if (action_name==="getPhenotypeHistImage"){
			if (displayResultsSpace.phenotype_method_id===YAHOO.util.Dom.get("phenotype_method_id").value){
				continue;
			}
			displayResultsSpace.imgQueryURL = PhenotypeHistImageBaseURL+'?'+url_query+'&img_type='+img_type;
			displayResultsSpace.fullImgURL = PhenotypeHistImageBaseURL+'?'+url_query+'&img_type='+fullImg;
		}
		else if (action_name==="getCallPhenotypeQQImage"){
			if (displayResultsSpace.phenotype_method_id===YAHOO.util.Dom.get("phenotype_method_id").value && displayResultsSpace.call_method_id===YAHOO.util.Dom.get("call_method_id").value ){
				continue;
			}
			displayResultsSpace.imgQueryURL = CallPhenotypeQQImageBaseURL+'?'+url_query+'&img_type='+img_type;
			displayResultsSpace.fullImgURL = CallPhenotypeQQImageBaseURL+'?'+url_query+'&img_type='+fullImg;
		}
		var replace_elem = document.getElementById(div_name);
		replace_elem.innerHTML = "<a href="+displayResultsSpace.fullImgURL+"><img src="+displayResultsSpace.imgQueryURL+"></a>";
		replace_elem.style.display = 'block';
		//displayStaticPlot(displayResultsSpace.imgQueryURL, div_name);
	}
	displayResultsSpace.phenotype_method_id=YAHOO.util.Dom.get("phenotype_method_id").value;
	displayResultsSpace.call_method_id=YAHOO.util.Dom.get("call_method_id").value;
	displayResultsSpace.analysis_method_id=YAHOO.util.Dom.get("analysis_method_id").value;
	return false;	//2008-12-30 return false to prevent form action 'server.html' from being invoked
}