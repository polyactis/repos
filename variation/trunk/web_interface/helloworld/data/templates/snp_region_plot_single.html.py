from mako import runtime, filters, cache
UNDEFINED = runtime.UNDEFINED
_magic_number = 2
_modified_time = 1223329143.736851
_template_filename='/usr/local/home_ubuntu/crocea/script/variation/web_interface/helloworld/helloworld/templates/snp_region_plot_single.html'
_template_uri='/snp_region_plot_single.html'
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
        context.write(u'\n')
        # SOURCE LINE 2
        context.write(u'\n\n<a href=')
        # SOURCE LINE 4
        context.write(unicode(h.url_for(controller='SNPRegionPlot', action='getImage', id=c.snp_region_plot.id)))
        context.write(u'><img width="1280" src=\'/SNPRegionPlot/getImage/')
        context.write(unicode(c.snp_region_plot.id))
        context.write(u"'></a>\n\n\n<table border=1>\n\n<tr>\n<th>Gene ID</th>\n<th>Symbol</th>\n<th>Type of Gene</th>\n<th>Chr</th>\n<th>Start</th>\n<th>Stop</th>\n<th>Protein Label</th>\n<th>Protein Comment</th>\n<th>Protein Text</th>\n</tr>\n")
        # SOURCE LINE 20
        for gene_desc_ls in c.matrix_of_gene_descriptions:
            # SOURCE LINE 21
            context.write(u'<tr>\n<td><a href=http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=gene&term=')
            # SOURCE LINE 22
            context.write(unicode(gene_desc_ls[0]))
            context.write(u'[uid]>')
            context.write(unicode(gene_desc_ls[0]))
            context.write(u'</a></td>\n<td>')
            # SOURCE LINE 23
            context.write(unicode(gene_desc_ls[1]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 24
            context.write(unicode(gene_desc_ls[2]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 25
            context.write(unicode(gene_desc_ls[3]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 26
            context.write(unicode(gene_desc_ls[4]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 27
            context.write(unicode(gene_desc_ls[5]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 28
            context.write(unicode(gene_desc_ls[6]))
            context.write(u'</td>\n<td>')
            # SOURCE LINE 29
            context.write(unicode(gene_desc_ls[7]))
            context.write(u'</td>\n</tr>\n')
        # SOURCE LINE 32
        context.write(u'\n\n</table>')
        return ''
    finally:
        context.caller_stack.pop_frame()


def render_title(context):
    context.caller_stack.push_frame()
    try:
        c = context.get('c', UNDEFINED)
        # SOURCE LINE 2
        context.write(u'Plot for SNPs at chr ')
        context.write(unicode(c.snp_region_plot.chromosome))
        context.write(u'. ')
        context.write(unicode(c.snp_region_plot.start))
        context.write(u'-')
        context.write(unicode(c.snp_region_plot.stop))
        context.write(u'. Phenotype ')
        context.write(unicode(c.snp_region_plot.phenotype_method.short_name))
        context.write(u' (id=')
        context.write(unicode(c.snp_region_plot.phenotype_method_id))
        context.write(u')')
        return ''
    finally:
        context.caller_stack.pop_frame()


