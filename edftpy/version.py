
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+gitaI.548fdb1"
__version__ = version
full_version = version

git_revision = "548fdb111a41969508901835812096f25bbf49c7"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
