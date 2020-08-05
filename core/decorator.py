#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author  : Yao

""" run shell commands

This file(script) can also be imported as a module and contains the following
class:
functions:
    * runshell - Run other tools on the command line and choose whether to keep its logs.

"""

# Standard library
import subprocess
import os
from functools import wraps

from core.exception import ShellCommandsException
from core.abstract import AbstractRun


devnull = open(os.devnull, "w")


def runshell(func):
    """
    func: return loggerfile, *args
        loggerfile: logger file handle
        args: command you want to run in shell (allow multiple command list/tuple which separated by ",")
    """

    @wraps(func)
    def _wrap(*args, **kwargs):
        shell_commands = func(*args, **kwargs)
        if shell_commands:
            if isinstance(shell_commands, str):
                shell_commands = (shell_commands, )
            shell_commands = list(shell_commands)
            while "" in shell_commands:
                shell_commands.remove("")
            # error_code == 0 means OK
            ok_code = [0]
            if isinstance(shell_commands[-1], int):
                ok_code.append(shell_commands[-1])
                shell_commands = shell_commands[:-1]
            if args and isinstance(args[0], AbstractRun):
                ini = args[0].ini
                accession = args[0].dm.id
                log_file = ini.pipe_log_file(accession)
                workflow_name = args[0].__class__.__base__.__name__
            else:
                workflow_name = kwargs.get("name")
                accession = kwargs.get("accession")
                ini = kwargs.get("ini")
                log_file = ini.pipe_log_file(accession)

            with open(os.path.join(ini.log_path, "RASERCMD"), "a+") as handle:
                handle.writelines("# {0}\n".format(workflow_name.upper()))
                handle.writelines("\n".join(shell_commands) + "\n\n")
            error_code = subprocess.call("&&".join(shell_commands), shell=True, stdout=open(log_file.o, "a+"), stderr=open(log_file.e, "a+"))
            if error_code not in ok_code:
                print("errorcode={0}".format(error_code), flush=True)
                print("\n".join(shell_commands), flush=True)
                raise ShellCommandsException("[{0}] ShellCommandsError: {1}".format(accession, workflow_name))
        return None
    return _wrap
