from __future__ import absolute_import
import logging

def getLogger():
    return logging.getLogger('wphase')

def debug(*a, **kw):
    getLogger().debug(*a, **kw)

def info(*a, **kw):
    getLogger().info(*a, **kw)

def warning(*a, **kw):
    getLogger().warning(*a, **kw)

def error(*a, **kw):
    getLogger().error(*a, **kw)
