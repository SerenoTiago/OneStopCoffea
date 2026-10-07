"""Analysis modules used only by explicitly selected diagnostic configs."""

import pkgutil

__all__ = []
for _loader, _module_name, _is_pkg in pkgutil.walk_packages(__path__):
    __all__.append(_module_name)

from . import *  # noqa: E402,F401,F403
