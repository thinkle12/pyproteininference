#! /usr/bin/env python
# macOS packaging support
from multiprocessing import freeze_support  # noqa

freeze_support()  # noqa

from nicegui import native, ui

import pyproteininference.gui.gui


def main():
    ui.run(native=True, title="pyProteinInference", reload=False, port=native.find_open_port())


if __name__ == "__main__":
    main()
