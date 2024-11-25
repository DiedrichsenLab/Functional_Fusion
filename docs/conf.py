# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import inspect
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../Functional_Fusion'))
import Functional_Fusion

# -- Project information -----------------------------------------------------

project = 'Functional Fusion Framework'
copyright = '2023, Diedrichsenlab'
author = 'Diedrichsenlab'

# The full version, including alpha/beta/rc tags
release = 'v.1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.napoleon',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosectionlabel',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx',
              'sphinx.ext.doctest',
              'nbsphinx',
              'sphinx.ext.viewcode']

#autodoc_member_order = 'bysource'

napoleon_custom_sections = [('Returns', 'params_style')]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['build']

# Make sure that class constructors are documented
autoclass_content = 'both'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

#html_theme_options = {
#    'collapse_navigation': True,
#    'sticky_navigation': True,
#    'navigation_depth': 4
#}

# Source Code available via sphinx
extensions.append('sphinx.ext.linkcode')

def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None

    obj = info['module']
    submodule = info['module']
    fullname = info['fullname']
    try:
        # Import the module and get the object
        mod = __import__(submodule, fromlist=[''])
        obj = mod
        for part in fullname.split('.'):
            obj = getattr(obj, part)
        
        # Get the source file and line number
        file = inspect.getsourcefile(obj)
        source, lineno = inspect.getsourcelines(obj)
        file = os.path.relpath(file, start=os.path.dirname(Functional_Fusion.__file__))
        return f"https://github.com/DiedrichsenLab/Functional_Fusion/blob/main/Functional_Fusion/{file}#L{lineno}"
    except Exception as e:
        print(f"Error in linkcode_resolve: {e}")
        return None