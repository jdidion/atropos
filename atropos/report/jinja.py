# import jinja2

# def generate_mako(stats, outfile, template, **kwargs):
#     from mako import Template
#     if not os.path.isabs(template):
#         from atropos import get_package_data
#         template = get_package_data('templates', template)
#     template = Template(filename=template)
#     args = dict(stats)
#     args.update(kwargs)
#     with open(outfile, 'wt') as out:
#         out.write(template.render(**args))
#
#
#     # Load the report template
#     try:
#         env = jinja2.Environment(loader=jinja2.FileSystemLoader(tmp_dir))
#         env.globals['include_file'] = include_file
#         j_template = env.get_template(template_mod.base_fn)
#     except:
#         raise IOError ("Could not load {} template file '{}'".format(config.template, template_mod.base_fn))
#
#     # Use jinja2 to render the template and overwrite
#     config.analysis_dir = [os.path.realpath(d) for d in config.analysis_dir]
#     report_output = j_template.render(report=report, config=config)
#     if filename == 'stdout':
#         print(report_output.encode('utf-8'), file = sys.stdout)
#     else:
#         try:
#             with io.open (config.output_fn, "w", encoding='utf-8') as f:
#                 print(report_output, file=f)
#         except IOError as e:
#             raise IOError ("Could not print report to '{}' - {}".format(config.output_fn, IOError(e)))
#
#         # Copy over files if requested by the theme
#         try:
#             for f in template_mod.copy_files:
#                 fn = os.path.join(tmp_dir, f)
#                 dest_dir = os.path.join( os.path.dirname(config.output_fn), f)
#                 copy_tree(fn, dest_dir)
#         except AttributeError:
#             pass # No files to copy
