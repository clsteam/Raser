#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Exception

"""


class CustomException(Exception):
    """custom data file is corrupted or mismatched"""
    pass


class DependencyException(Exception):
    """Operating environment dependence"""
    pass


class SettingException(Exception):
    """There is a problem in the settings file"""
    pass


class ConfigException(Exception):
    """There is a problem in the configure file"""
    pass


class ShellCommandsException(Exception):
    """Shell command cannot run"""
    pass
