//2009-3-5 common javascript functions


//string format for javascript: var result = "Hello {0}! This is {1}.".format("world","foo bar");
//from http://snipplr.com/view/8984/sprintf-in-javascript-string-format/
String.prototype.format = function(){
	var pattern = /\{\d+\}/g;
	var args = arguments;
	return this.replace(pattern, function(capture){ return args[capture.match(/\d+/)]; });
}