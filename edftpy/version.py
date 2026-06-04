
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+git20260603.e66a3ba"
__version__ = version
full_version = version

git_revision = "e66a3ba7f6b3c8c3847a9701a47c68e20a16152f"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
