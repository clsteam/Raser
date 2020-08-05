# -*- coding:utf-8 -*-
import functools


def workflow_stateful(cls):

    old__init__ = cls.__init__
    if hasattr(cls, '__getattr__'):
        old__getattr__ = getattr(cls, '__getattr__')
    else:
        old__getattr__ = getattr(cls, '__getattribute__')

    def __init__(self, *args, **kwargs):
        self.__tool__ = None
        return old__init__(self, *args, **kwargs)

    def __getattr__(self, name):
        if self.__dict__["__tool__"] is None:
            raise Exception(
                "{0} doesn't has initialize.".format(self.__dict__["workflow_name"]))
        try:
            old__getattr__(self, name)
        except AttributeError:
            pass

        try:
            f = getattr(curr(self), name)
        except AttributeError as e:
            raise e

        if callable(f):
            return functools.partial(f, self)
        return f

    cls.__init__ = __init__
    cls.__getattr__ = __getattr__
    return cls


def curr(host):
    return host.__tool__


def switch_tools(host, new_tool):
    host.__tool__ = new_tool
