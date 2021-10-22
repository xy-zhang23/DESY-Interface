# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 17:05:47 2020

@author: lixiangk
"""

from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name = 'interface',
   version = '1.0.0',
   description = 'An interface for particle tracking with Astra',
   license = "Free License",
   long_description = long_description,
   author = 'Xiangkun Li',
   author_email = 'xiangkun.li@desy.de',
   url = "https://bitbucket.org/XiangkunLi/interface/",
   packages = ['interface'],  # same as name
   install_requires = [],     # external packages as dependencies
   examples = [
            'examples/dipole_demo',
            'examples/injector_demo',
            'examples/injector_optimization_demo'
           ]
)

