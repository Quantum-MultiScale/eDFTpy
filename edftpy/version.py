
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+gitaI.00220da"
__version__ = version
full_version = version

git_revision = "00220da8b9211bda4f56a4c8c17bda99326bacbe"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
