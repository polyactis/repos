from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223276022.8891349
_template_filename=u'/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/base.html'
_template_uri=u'/base.html'
_template_cache=cache.Cache(__name__, _modified_time)
_source_encoding=None
_exports = []


# SOURCE LINE 1

import datetime


def render_body(context,**pageargs):
    context.caller_stack.push_frame()
    try:
        __M_locals = dict(pageargs=pageargs)
        self = context.get('self', UNDEFINED)
        # SOURCE LINE 3
        context.write(u'\n\n<html>\n    <head>\n        <title>')
        # SOURCE LINE 7
        context.write(unicode(self.title()))
        context.write(u'</title>\n    </head>\n    <body>\n        ')
        # SOURCE LINE 10
        context.write(unicode(self.body()))
        context.write(u'\n        <div class="footer">\n            <p>Data presented from Database Stock_250kDB in <a href=http://papaya.usc.edu/>Nordborg Lab</a>. Page generated at ')
        # SOURCE LINE 12
        context.write(unicode(str(datetime.datetime.now())))
        context.write(u'.</p>\n        </div>\n    </body>\n</html>\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


