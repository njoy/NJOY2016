#! /usr/bin/env python3

import sys

import os.path
sys.path.append( os.path.abspath( os.getcwd() ) )
import metaconfigure.description as description
import metaconfigure.cmake as cmake

previous = description.deserialize()
keyword_arguments = {}
keywords = [ "include_path", "language_revision" ]

for keyword in keywords:
    if keyword in previous:
        keyword_arguments[keyword] = previous[keyword]

description.generate( previous["name"],
                      previous["target"],
                      previous["language"],
                      previous["version"],
                      previous["is_external_project"],
                      **keyword_arguments )

if sys.argv[1] == 'cmake':
    cmake.generate()
