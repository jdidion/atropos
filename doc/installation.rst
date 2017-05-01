============
Installation
============

Quickstart
----------

The easiest way to install atropos is to use ``pip`` on the command line::

    pip install --user --upgrade atropos

This will download the software from `PyPI (the Python packaging
index) <https://pypi.python.org/pypi/atropos/>`_, and
install the atropos binary into ``$HOME/.local/bin``. If an old version of
atropos exists on your system, the ``--upgrade`` parameter is required in order
to install a newer version. You can then run the program like this::

    ~/.local/bin/atropos --help

If you want to avoid typing the full path, add the directory
``$HOME/.local/bin`` to your ``$PATH`` environment variable.


Conda
-----

Atropos can also be installed using conda::

    conda install atropos


Docker
------

We provide a Docker container that can be executed using a Docker or Singularity
engine:

    docker run jdidion/atropos <args>

or 

    singularity run docker://jdidion/atropos


Manual installation
-------------------

Dependencies
~~~~~~~~~~~~

Atropos requires this software to be installed:

* One of Python 3.3+.
* A C compiler.
* Cython 0.25.2+

Under Ubuntu, you may need to install the packages ``build-essential`` and
``python-dev``.


Installation
~~~~~~~~~~~~

If you have already downloaded and unpacked the ``.tar.gz`` file, then
installation is done like this::

    python3 setup.py install --user

If you get an error message::

    error: command 'gcc' failed with exit status 1

Then check the entire error message. If it says something about a missing ``Python.h``
file, then you need to install the Python development packages. The
appropriate package is called ``python3-dev`` in Ubuntu.sour

We also provide a Makefile that can optionally run unit tests for you. For this 
you will need to have make installed and also the Python pytest library (if you 
want to run the tests). To both test and install, run::

    make installargs='--user'

To only install, run::

    make install installargs='--user'

System-wide installation
~~~~~~~~~~~~~~~~~~~~~~~~

If you have root access, then you can install atropos system-wide by running::

    sudo pip install atropos

This installs atropos into `/usr/local/bin`.

If you want to upgrade from an older version, use this command instead::

    sudo pip install --upgrade atropos

To use the Makefile, simply run::

    make

or (to skip running tests)

    make install

Use without installation
~~~~~~~~~~~~~~~~~~~~~~~~

Build the C extension module (you can try to skip this step -- a
compiled version of the module for Linux x86\_64 is already included)::

    python setup.py build_ext -i

Then simply run the script from where it is, similar to this::

    bin/atropos --help

If you get any errors, first try to explicitly request a specific Python
version by running atropos like this::

    python3 bin/atropos --help
