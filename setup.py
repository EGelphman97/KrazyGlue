"""
Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  September 8, 2017

  setup.py file needed to build the _Fulmine module
"""

from distutils.core import setup, Extension

setup(name="fulmine", version="1.0", ext_modules=[Extension("fulmine", ["Fulmine.c"])])
