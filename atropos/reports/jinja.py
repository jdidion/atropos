"""Generate a report from a Jinja2 template.
"""
import jinja2
from atropos import get_package_path

def generate_report(
        fmt, summary, outfile, template_name=None, template_paths=None,
        template_globals=None):
    """Generate a report.
    
    Args:
        fmt: Report format. Must match the file extension on a discoverable
            template.
        summary: The summary dict.
        outfile: The output file name/prefix.
        template_name: A template name to use, rather than auto-discover.
        template_paths: Sequence of paths to search for templates.
        template_globals: Dict of additional globals to add to the template
            environment.
    """
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
    except IOError as err:
        raise IOError(
            "Could not print report to '{}' - {}".format(outfile, err))
    finally:
        if is_path:
            outfile.close()
