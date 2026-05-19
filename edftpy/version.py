
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+gitaI.926fb13"
__version__ = version
full_version = version

git_revision = "926fb13e51a4488bc55adc37eb925bf7cd392579"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
