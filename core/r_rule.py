#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Title

This file(script) can also be imported as a module and contains the following
functions:

    * main - the main function of the script
    * function - returns the column headers of the file
"""

# Standard library

# Third party library

# Private module


def conversion_function_symbol(script: str) -> str:
    """Replace curly braces"""

    return script.replace("@<", "{").replace(">@", "}")


def list_to_vector(ls: list) -> str:
    return "c(" + ", ".join(["".join(["'", x, "'"]) for x in ls]) + ")"
