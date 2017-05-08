#! /usr/bin/env python
"""
"""

import os
import json
import subprocess
import sys
import time

def project_signature( git ):
    project_name = os.path.split( os.getcwd() )[1]
    invocation = [ "git", "rev-parse", "HEAD" ]
    if os.name == "nt":
        invocation.insert( 0, "powershell" )

    process = subprocess.Popen( invocation, stdout=subprocess.PIPE )
    hash_value = process.communicate()
    git[ project_name ] = str(hash_value[0].strip())

def generate( name = None ):
    if name is None:
        project_name = os.path.split( os.getcwd() )[1]
        name = project_name + str( time.time() )

    git = {}
    project_signature( git )
    if os.path.isdir( os.path.join( os.getcwd(), 'subprojects' ) ):
        os.chdir( 'subprojects' )
        root = os.getcwd()
        for subproject in os.listdir( os.getcwd() ):
            if os.path.isdir( os.path.join( os.getcwd(), subproject ) ):
                os.chdir( subproject )
                project_signature( git )
                os.chdir( root )

        os.chdir( '..' )

    with open ( name + ".json", "w" ) as json_file:
        json_file.write( json.dumps( { 'git' : git }, indent=0 ) )

if len(sys.argv) > 1:
    generate( sys.argv[1] )
else:
    generate()
