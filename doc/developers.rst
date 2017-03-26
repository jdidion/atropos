Style-wise, we try to adhere to the [Google python style guidelines](https://google.github.io/styleguide/pyguide.html).

We use Google-style docstrings, which are formatted by the [Napoleon Sphinx Plugin](https://pypi.python.org/pypi/sphinxcontrib-napoleon).

We run pylint as part of each build and strive to maintain a 10/10 score. However, we disable some pylint checks:

* Function annotations: pylint does not properly handle whitespace around function annotations (https://github.com/PyCQA/pylint/issues/238).
* White space on empty lines: we use white space as a visual guide to the structure of the code. Each blank line should have whitespace matching the indent level of the next non-blank line.
* Checks that are arbitrary/overly restrictive (e.g. 'too-many-xxx'; see .pylintrc for full list)
