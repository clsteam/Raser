#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
import re
import commands
import random
import subprocess
from configparser import *

global conf_subject
conf_subject = configparser()

class Color:

    def __init__(self):
        self.conf_subject = configparser()
        self.Master_switch=conf_subject.get_bool("Log_output_switch", "color_display")
        self.warning=self.conf_subject.get_bool("Log_output_switch", "show_warning")
        self.debug=self.conf_subject.get_bool("Log_output_switch","show_debug")

    def print_warning(self,text):
        if self.warning:
            text="WARNING:"+text+"\n"
            if self.Master_switch:
                sys.stderr.write(self.convert(text,"yellow"))
                sys.stderr.flush()
            else:
                sys.stderr.write(text)
                sys.stderr.flush()

    def print_debug(self,text):
        if self.debug:
            text=text+"\n"
            if self.Master_switch:
                sys.stdout.write(self.convert(text,"purple"))
                sys.stdout.flush()
            else:
                sys.stdout.write(text)
                sys.stdout.flush()

    def flush_print(self,text):
        text=text+"\n"
        if self.Master_switch:
            sys.stdout.write(self.convert(text,"sky blue"))
            sys.stdout.flush()
        else:
            sys.stdout.write(text)
            sys.stdout.flush()

    def fatal_error(self,text):
        text="FATAL ERROR:"+text+"\n"
        if self.Master_switch:
            sys.stdout.write(self.convert(text,"red"))
            sys.stdout.flush()
        else:
            sys.stdout.write(text)
            sys.stdout.flush()

    def convert(self,text,mark="yellow"):
        if mark == "yellow":
            return "\033[33m" + text + "\033[0m"
        elif mark == "sky blue":
            return "\033[36m" + text + "\033[0m"
        elif mark == "green":
            return "\033[32m" + text + "\033[0m"
        elif mark == "purple":
            return "\033[35m" + text + "\033[0m"
        elif mark == "blue":
            return "\033[34m" + text + "\033[0m"
        elif mark == "red":
            return "\033[31m" + text + "\033[0m"
        else:  #Blue background Red word
            return "\033[34;31m" + text + "\033[0m"

def get_absolute_path(in_path):
    if in_path.startswith("~"):
        return os.path.expanduser(in_path)
    else:
        return os.path.abspath(in_path)

def free_nodes(n):
    # pang nodes(need to apply)
    # search
    # ['', '', '', 'offline', 'unknown', 'comput3-24/24','comput43-14/32','comput43-14/32']
    # available_nodes={1:[comput3,comput4],5:[comput40,comput25]}
    pang_nodes=["comput3","comput43","comput54"]
    total_nodes_num=52
    available_nodes={}
    pat=re.compile('3\dm\s(.*?)\s')
    node_message=re.findall(pat,commands.getoutput('pnodes'))
    for i in range(len(node_message)):
        if re.match(re.compile('^comput\d+-\d+/\d+$'),node_message[i]):
            node_message=node_message[i:]
            break
    for mess in node_message:
        node=mess.split("-")[0]
        cup_message=mess.split("-")[1].split("/")
        free_cpu_num=int(cup_message[1])-int(cup_message[0])
        if node not in pang_nodes and free_cpu_num>0:
            if available_nodes.has_key(free_cpu_num):
                available_nodes[free_cpu_num].append(node)
            else:
                available_nodes[free_cpu_num]=[node]
    if available_nodes=={}:
        return "comput"+str(random.randint(3,42))
    for i in range(n,25):
        if available_nodes.has_key(i):
            return available_nodes[i][0]
        else:
            continue
    return 0

def run(command):
    global conf_subject
    devnull = open(os.devnull, 'w')
    if conf_subject.get_bool("Log_output_switch","show_tools_original_log"):
        subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)
    else:
        subprocess.call(command, shell=True)


