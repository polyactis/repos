//2009-1-30 try to write a prototype to extend google.visualization.ScatterChart
var example = {'name':'me'};

//2009-1-30 legendary function talked in Douglas Crockford's talk, but not that useful
function object(o) {
	function F() {}
	F.prototype = o;
	return new F();
	}

// Class constructor. Parameter container is a DOM elementon the client that
// that will contain the visualization.
google.load("visualization", "1", {packages:["scatterchart"]});

example.yh_scatter = google.visualization.ScatterChart;


//example.yh_scatter = function(container) {
//  this.containerElement = container;
//}

// Main drawing logic.
// Parameters:
//   data is data to display, type google.visualization.DataTable.
//   options is a name/value map of options. Our example takes one option.
//example.yh_scatter.prototype.draw = function(data, options) {
//	google.visualization.ScatterChart(data, options);
//}