#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

import ConfigParser;
import os
import sys

def return_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

class configparser:
    def __init__(self,conf_file=return_script_path()+"/reference.conf"):
        self.conf_file=conf_file
        self.config = ConfigParser.ConfigParser()
        self.config.readfp(open(conf_file))

    def add_values(self,a,b,c):
        self.config.set(a,b,c)
        with open(self.conf_file,'wb') as configfile:
            self.config.write(configfile)

    def get_values(self,a,b):
        return self.config.get(a,b)

    def has_option(self,a,b):
        return self.config.has_option(a,b)

    def get_bool(self,a,b):
        return self.config.getboolean(a,b)