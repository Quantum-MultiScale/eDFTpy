
"""
Module to expose more detailed version info for the installed `dftpy`
"""
version = "0.0.1dev0+gitaI.67b55a2"
__version__ = version
full_version = version

git_revision = "67b55a2d76fc94bed7277efec87f6e35511f763b"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
