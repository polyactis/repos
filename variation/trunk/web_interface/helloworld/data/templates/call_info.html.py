from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223278736.6544781
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/call_info.html'
_template_uri='/call_info.html'
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
        context.write(u'\n\n')
        # SOURCE LINE 3
        context.write(u'\n\n<p> ')
        # SOURCE LINE 5
        context.write(unicode(len(c.call_info_ls)))
        context.write(u' Genotype Calls In Total.</p>\n\n<table border=1>\n\n<tr>\n<th>No.</th>\n<th>Call Info ID</th>\n<th>Ecotype ID</th>\n<th>Nativename</th>\n<th>Stockparent</th>\n<th>Call Method ID</th>\n<th>Filename</th>\n<th>Array ID</th>\n<th>Original Array Filename</th>\n\n</tr>\n')
        # SOURCE LINE 21
        for call_info in c.call_info_ls:
            # SOURCE LINE 22
            context.write(u'<tr>\n<td>')
            # SOURCE LINE 23
            context.write(unicode(call_info.no))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 24
            context.write(unicode(call_info.call_info_id))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 25
            context.write(unicode(call_info.ecotype_id))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 26
            context.write(unicode(call_info.nativename.decode('utf-8', 'ignore')))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 27
            context.write(unicode(call_info.stockparent))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 28
            context.write(unicode(call_info.call_method_id))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 29
            context.write(unicode(call_info.filename))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 30
            context.write(unicode(call_info.array_id))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 31
            context.write(unicode(call_info.original_filename))
            context.write(u'</td>\n\n</tr>\n')
        # SOURCE LINE 35
        context.write(u'\n\n</table>\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


def render_title(context):
    context.caller_stack.push_frame()
    try:
        # SOURCE LINE 3
        context.write(u'Call Info')
        return ''
    finally:
        context.caller_stack.pop_frame()


