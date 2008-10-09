from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223260413.4872279
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/navigation.html'
_template_uri='/navigation.html'
_template_cache=cache.Cache(__name__, _modified_time)
_source_encoding=None
_exports = ['navigation_links']


def render_body(context,**pageargs):
    context.caller_stack.push_frame()
    try:
        __M_locals = dict(pageargs=pageargs)
        def navigation_links(selected,links):
            return render_navigation_links(context.locals_(__M_locals),selected,links)
        # SOURCE LINE 1

        items = [
            ('James', 'http://jimmyg.org'),
            ('Ben', 'http://groovie.org'),
            ('Philip', ''),
        ]
        
        
        __M_locals.update(dict([(__M_key, locals()[__M_key]) for __M_key in ['items'] if __M_key in locals()]))
        # SOURCE LINE 7
        context.write(u'\n\n')
        # SOURCE LINE 9
        context.write(unicode(navigation_links('James', items)))
        context.write(u'\n\n')
        # SOURCE LINE 31
        context.write(u'\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


def render_navigation_links(context,selected,links):
    context.caller_stack.push_frame()
    try:
        def link(label,url):
            context.caller_stack.push_frame()
            try:
                # SOURCE LINE 12
                context.write(u'\n')
                # SOURCE LINE 13
                if url:
                    # SOURCE LINE 14
                    context.write(u'            <a href="')
                    context.write(unicode(url))
                    context.write(u'">')
                    context.write(unicode(label))
                    context.write(u'</a>\n')
                    # SOURCE LINE 15
                else:
                    # SOURCE LINE 16
                    context.write(u'            ')
                    context.write(unicode(label))
                    context.write(u'\n')
                # SOURCE LINE 18
                context.write(u'    ')
                return ''
            finally:
                context.caller_stack.pop_frame()
        # SOURCE LINE 11
        context.write(u'\n    ')
        # SOURCE LINE 18
        context.write(u'\n\n    <ul>\n')
        # SOURCE LINE 21
        for item in links:
            # SOURCE LINE 22
            context.write(u'        <li>')
            # SOURCE LINE 23
            if item[0] == selected:
                # SOURCE LINE 24
                context.write(u'        <b>')
                context.write(unicode(link(item[0], item[1])))
                context.write(u'</b>')
                # SOURCE LINE 25
            else:
                # SOURCE LINE 26
                context.write(u'        ')
                context.write(unicode(link(item[0], item[1])))
                context.write(u'')
            # SOURCE LINE 28
            context.write(u'        </li>\n')
        # SOURCE LINE 30
        context.write(u'    </ul>\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


