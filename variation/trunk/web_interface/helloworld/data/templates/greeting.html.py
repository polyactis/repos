from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223262439.258507
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/greeting.html'
_template_uri='/greeting.html'
_template_cache=cache.Cache(__name__, _modified_time)
_source_encoding=None
_exports = ['title']


def _mako_get_namespace(context, name):
    try:
        return context.namespaces[(__name__, name)]
    except KeyError:
        _mako_generate_namespaces(context)
        return context.namespaces[(__name__, name)]
def _mako_generate_namespaces(context):
    pass
def _mako_inherit(template, context):
    _mako_generate_namespaces(context)
    return runtime._inherit_from(context, u'/base.html', _template_uri)
def render_body(context,**pageargs):
    context.caller_stack.push_frame()
    try:
        __M_locals = dict(pageargs=pageargs)
        c = context.get('c', UNDEFINED)
        # SOURCE LINE 1
        context.write(u'')
        # SOURCE LINE 2
        context.write(u'\n<h1>Greetings</h1>\n\n<p>')
        # SOURCE LINE 5
        context.write(unicode(c.greeting))
        context.write(u' ')
        context.write(unicode(c.name))
        context.write(u'!</p>\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


def render_title(context):
    context.caller_stack.push_frame()
    try:
        # SOURCE LINE 2
        context.write(u'Greetings')
        return ''
    finally:
        context.caller_stack.pop_frame()


