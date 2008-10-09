<h2>
Server info for ${request.host}
</h2>

<p>
The URL you called: ${h.url_for()}
</p>



<p>
the name you set is: ${c.name}
</p>

<%
    context.write('<p>Here is an example:</p>')
%>
<%doc>
<p>
% for key in context.keys():
The key is <tt>${key}</tt>, the value is ${str(context.get(key))}. <br />
% endfor
</p>
</%doc>


<p>
The WSGI environ:<br />
<pre>${c.pretty_environ}</pre>
</p>

