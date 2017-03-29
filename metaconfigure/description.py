"""
A script to generate a serialized project description for build script generation.
"""
import os
import glob
import json
from metaconfigure import configuration

def serialize( state ):
    with open ("metaconfigure/description.json", "w") as json_file:
        project_path = state.pop('project path')
        subprojects = state.pop('subprojects')
        json_file.write( json.dumps( state, indent=2, sort_keys=True ) )
        state['subprojects'] = subprojects
        state['project path'] = project_path

def deserialize():
    with open ("metaconfigure/description.json", "r") as json_file:
        state = json.loads( json_file.read() )
        state['project path'] = os.getcwd()
        state['subprojects'] = {}
        return state

def collect_subprojects( state, root ):
    if os.path.isdir( os.path.join( os.getcwd(), 'dependencies' ) ):
        os.chdir('dependencies')
        anchor = os.getcwd()
        for name in os.listdir( os.getcwd() ) :
            if os.path.isdir( os.path.join( os.getcwd(), name ) ) :
                os.chdir( os.path.join( root, 'subprojects', name ) )
                subproject = deserialize()
                collect_subprojects( subproject, root )
                state['subprojects'][name] = subproject
                os.chdir( anchor )
            
        os.chdir('..')

def reconstruct_dependency_graph( state ):
    dependencies = { state['name'] : state['subprojects'].keys() }
    for name in state['subprojects'].keys():
        dependencies.update(
            reconstruct_dependency_graph( state['subprojects'][name] ) )
      
    return dependencies

def reconstruct_build_queue( state ):
    queue = []
    graph = reconstruct_dependency_graph( state )
    while len( graph ) > 0:
        found_next = False
        next_node = None
        for node, edges in graph.items():
            queueable = True
            for edge in edges:
                if not edge in queue :
                    queueable = False
                    break
              
            if queueable:
                next_node = node
                found_next = True
                break
          
        if not found_next:
            raise RuntimeError('Cyclic Dependencies?!')
        
        else:
            queue.append(next_node)
            graph.pop(next_node)
        
    return queue

def relative_path( state ):
    return os.path.relpath( os.getcwd(), state['project path'] )

def evaluate_test_directory( state ):
    path, name = os.path.split( os.path.split( os.getcwd() )[0] )
    while True:
        path, directory = os.path.split( path )
        if not directory == "src":
            name = '.'.join( [directory, name] )
        else:
            break

    files = []
    for extension in state['file extension']['implementation files'] :
        filenames = glob.glob( '*.' + extension )
        file_paths = [ os.path.join( directory, filename ) for filename in filenames ]
        files.extend( file_paths )
            
    state['tests'][name] = files
    os.chdir('..')

def evaluate_directory( state ):
    path = relative_path( state )
    for group in state['file extension']:
        files = []
        for extension in state['file extension'][group] :
            filenames = glob.glob( '*.' + extension )
            file_paths = [ os.path.join( path, filename ) for filename in filenames ]
            files.extend( file_paths )
            
        state[group].extend( files )
        
    for name in os.listdir( os.getcwd() ):
        if os.path.isdir( os.path.join( os.getcwd(), name ) ):
            os.chdir(name)
            if name == 'test':
                evaluate_test_directory( state)
            else:
                evaluate_directory( state )
          
    os.chdir('..')

def collect_driver( state ):
    if 'driver' not in state:
        driver = None
        for extension in state['file extension']['implementation files']:
            possible_driver = 'src/main.' + extension
            if possible_driver in state['implementation files']:
                driver = possible_driver
                break
            return
        
        state['driver'] = driver
        
    state['implementation files'].remove(driver)
    
def generate( name, language, is_external_project = False, **kwargs ):
    state = locals()
    state.update( state.pop( 'kwargs' ) )
    if 'initialized' in state:
        state.update( configuration.languages[ state['language'] ] )
        state['initialized'] = True
    
    state['project path'] = os.getcwd()
    state['tests'] = {}
    state['subprojects'] = {}
    if 'strict' not in state:
        strict = True
    
    for group in state['file extension']:    
        state[group] = []
    
    os.chdir( 'src' )
    evaluate_directory( state )
    collect_driver( state )
    serialize( state )
