from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223338392.2885959
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/snp_region_plot.html'
_template_uri='/snp_region_plot.html'
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
        h = context.get('h', UNDEFINED)
        c = context.get('c', UNDEFINED)
        # SOURCE LINE 1
        context.write(u'\n\n')
        # SOURCE LINE 3
        context.write(u'\n\n\n<p> ')
        # SOURCE LINE 6
        context.write(unicode(len(c.snp_region_plots)))
        context.write(u' SNP Region Plots In Total.</p>\n\n<table border=1>\n\n<tr>\n<th>ID</th>\n<th>Chromosome</th>\n<th>Start</th>\n<th>Stop</th>\n<th>Center SNP</th>\n<th>Phenotype</th>\n<th>No of Genes</th>\n<th>Plot Type</th>\n<th>created by</th>\n<th>date created</th>\n</tr>\n')
        # SOURCE LINE 22
        for snp_region_plot in c.snp_region_plots:
            # SOURCE LINE 23
            context.write(u'<tr>\n<td><a href=')
            # SOURCE LINE 24
            context.write(unicode(h.url_for(controller='SNPRegionPlot', action='show_plot', id=snp_region_plot.id)))
            context.write(u'>')
            context.write(unicode(snp_region_plot.id))
            context.write(u'</a></td>\n<td>')
            # SOURCE LINE 25
            context.write(unicode(snp_region_plot.chromosome))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 26
            context.write(unicode(snp_region_plot.start))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 27
            context.write(unicode(snp_region_plot.stop))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 28
            context.write(unicode(snp_region_plot.center_snp_position))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 29
            context.write(unicode(snp_region_plot.phenotype_method.short_name))
            context.write(u' (id=')
            context.write(unicode(snp_region_plot.phenotype_method_id))
            context.write(u')</td>\n<td>')
            # SOURCE LINE 30
            context.write(unicode(len(snp_region_plot.plot2gene_ls)))
            context.write(u'</td>\n<td><a href=')
            # SOURCE LINE 31
            context.write(unicode(h.url_for(controller='SNPRegionPlot', action='type', id=snp_region_plot.plot_type_id)))
            context.write(u'>')
            context.write(unicode(snp_region_plot.plot_type.short_name))
            context.write(u'</a></td>\n<td>')
            # SOURCE LINE 32
            context.write(unicode(snp_region_plot.created_by))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 33
            context.write(unicode(snp_region_plot.date_created))
            context.write(u'</td>\n</tr>\n')
        # SOURCE LINE 36
        context.write(u'\n\n</table>\n')
        return ''
    finally:
        context.caller_stack.pop_frame()


def render_title(context):
    context.caller_stack.push_frame()
    try:
        # SOURCE LINE 3
        context.write(u'SNP Region Plots')
        return ''
    finally:
        context.caller_stack.pop_frame()


