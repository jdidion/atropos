import pkg_resources


UNKNOWN_VERSION = "Unknown"


try:
    __version__ = pkg_resources.get_distribution(__name__).version
except pkg_resources.DistributionNotFound:
    __version__ = UNKNOWN_VERSION
