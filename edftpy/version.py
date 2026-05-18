
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+gitaI.f0b7245"
__version__ = version
full_version = version

git_revision = "f0b7245dd201a44fe1b4ef78bf54fc259db2666a"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
