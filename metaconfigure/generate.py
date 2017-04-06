#! /usr/bin/env python3

import sys

import os.path
sys.path.append( os.path.abspath( os.getcwd() ) )
import metaconfigure.description as description
import metaconfigure.cmake as cmake

previous = description.deserialize()

if previous['update']:
    description.generate( **previous )

if sys.argv[1] == 'cmake':
    cmake.generate()
