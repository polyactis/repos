from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223257472.276562
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/serverinfo.mako'
_template_uri='/serverinfo.mako'
_template_cache=cache.Cache(__name__, _modified_time)
_source_encoding=None
_exports = []


def render_body(context,**pageargs):
    context.caller_stack.push_frame()
    try:
        __M_locals = dict(pageargs=pageargs)
        h = context.get('h', UNDEFINED)
        c = context.get('c', UNDEFINED)
        request = context.get('request', UNDEFINED)
        # SOURCE LINE 1
        context.write(u'<h2>\nServer info for ')
        # SOURCE LINE 2
        context.write(unicode(request.host))
        context.write(u'\n</h2>\n\n<p>\nThe URL you called: ')
        # SOURCE LINE 6
        context.write(unicode(h.url_for()))
        context.write(u'\n</p>\n\n\n\n<p>\nthe name you set is: ')
        # SOURCE LINE 12
        context.write(unicode(c.name))
        context.write(u'\n</p>\n\n')
        # SOURCE LINE 15

        context.write('<p>Here is an example:</p>')
        
        
        # SOURCE LINE 17
        context.write(u'\n<p>\n')
        # SOURCE LINE 19
        for key in context.keys():
            # SOURCE LINE 20
            context.write(u'The key is <tt>')
            context.write(unicode(key))
            context.write(u'</tt>, the value is ')
            context.write(unicode(str(context.get(key))))
            context.write(u'. <br />\n')
        # SOURCE LINE 22
        context.write(u'</p>\n\n\n<p>\nThe WSGI environ:<br />\n<pre>')
        # SOURCE LINE 27
        context.write(unicode(c.pretty_environ))
        context.write(u'</pre>\n</p>\n\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


