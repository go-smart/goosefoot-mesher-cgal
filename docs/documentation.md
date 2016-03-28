# Documentation generation

The documentation may be compiled and sent to Github pages by the following commands.

~~~~~~~~~~{.sh}
cd docs/
doxygen Doxyfile
ghp-import -n html
git push origin gh-pages
~~~~~~~~~~

Note that this requires the Python [ghp-import script](https://pypi.python.org/pypi/ghp-import).
