#! /usr/bin/env python
"""
A script to generate a flat subproject directory from a tree of dependency directories
"""
import os
import subprocess
import textwrap
import json
import shutil
import sys
import warnings

def project_name():
    return os.path.split( os.getcwd() )[1]

def dependency_directory():
    """
    A function isolating the path specification to the project's dependency directory
    """
    return os.path.join( os.getcwd(), "dependencies" )

def traverse_dependencies( destination, traversed ):
    """
    collect unique links to dependency projects in a destination folder
    """
    if not os.path.isdir( dependency_directory() ):
        return
    
    os.chdir( dependency_directory() )

    for dependency in os.listdir( os.getcwd() ) :
        if os.path.isdir( dependency ) and not dependency in traversed :
            traversed.add( dependency )
            os.chdir( dependency )
            if not os.path.isdir( os.path.join( destination, dependency ) ):
                try:
                    os.symlink( os.getcwd(),
                                os.path.join( destination, dependency ) )

                except OSError:
                    warnings.warn( "Could not create symbolic "
                                   "link from {} to subprojects directory."\
                                   .format( os.getcwd() ) )
                    warnings.warn( "Copying directory contents instead" )
                    shutil.copytree( os.getcwd(),
                                     os.path.join( destination, dependency ),
                                     ignore = shutil.ignore_patterns("dependencies") )

            traverse_dependencies( destination, traversed )
            os.chdir( ".." )
            
    os.chdir( os.path.join( ".." ) )

def collect_subprojects():
    destination = os.path.join( os.getcwd(), "subprojects" )
    if not os.path.isdir( destination ):
        os.makedirs( destination )
        
    traverse_dependencies( destination, set() )

collect_subprojects()
