__version__ = "0.1"

import contextlib
from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from pyCGMapPol.pyCGMapPol import main_app

    #with contextlib.ExitStack() as ctx:
    return main_app(__version__)
