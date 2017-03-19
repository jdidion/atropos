import jinja2
from atropos import get_package_path

def generate_report(
        fmt, summary, outfile, template_name=None, template_paths=None,
        template_globals=None):
    if not template_name:
        template_name = 'template.{}'.format(fmt)
    if not template_paths:
        template_paths = []
    template_paths.append(get_package_path('reports', 'templates'))
    
    # Load the report template
    try:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_paths))
        if template_globals:
            env.globals.update(template_globals)
        template = env.get_template(template_name)
    except:
        raise IOError("Could not load template file '{}'".format(template_name))

    # Render the template
    report_output = template.render(summary=summary)
    
    # Write to output
    is_path = isinstance(outfile, str)
    if is_path:
        outfile = open(outfile, 'wt', encoding='utf-8')
    else:
        report_output = report_output.encode('utf-8')
    
    try:
        print(report_output, file=outfile)
    except IOError as e:
        raise IOError("Could not print report to '{}' - {}".format(outfile, e))
    finally:
        if is_path:
            outfile.close()
