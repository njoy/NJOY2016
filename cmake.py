import os
import textwrap
from . import description

language = {'c' : 'C', 'c++' : 'CXX', 'fortran' : 'Fortran'}
vendor = {'gcc' : 'GNU',
          'g++' : 'GNU',
          'gfortran' : 'GNU',
          'llvm clang' : 'Clang',
          'llvm clang++' : 'Clang',
          'apple clang' : 'AppleClang',
          'apple clang++' : 'AppleClang'}

platform = {'linux':'Linux',
            'osx':'Darwin',
            'windows':'Windows'}

def fetch_subprojects( state ):
  contents = ""
  if state['subprojects']:
    contents += textwrap.dedent(
    """
    if( NOT ROOT_DIRECTORY )
        set( ROOT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
        if ( NOT fetched_subprojects )
            if ( NOT PYTHON_EXECUTABLE )
                find_package( PythonInterp )
                if ( NOT PYTHONINTERP_FOUND )
                    message( FATAL_ERROR "Python interpeter installation was not found." )
                endif()
            endif()
            execute_process( COMMAND ${PYTHON_EXECUTABLE} "./metaconfigure/fetch_subprojects.py"
                             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
                             RESULT_VARIABLE fetch_failure )
            if ( NOT fetch_failure )
                set( fetched_subprojects TRUE CACHE BOOL "fetch script ran" )
            else()
                message( FATAL_ERROR "Failed to fetch dependencies" )
            endif()
        endif()
    endif()
    """)
    
    return contents

def subproject_languages( state ):
    languages = set()
    if 'subprojects' in state:
        for subproject in state['subprojects'].values():
            languages.add( subproject['language'] )
            languages = languages | subproject_languages( subproject )
            
    return languages

def has_library( state ):
    return bool(state['implementation files'])

def has_executable( state ):
    return 'driver' in state 

def has_tests( state ):
    return bool(state['tests']) or os.path.isdir("tests") 
    
def project_statement( state ):
    project_language = language[ state['language'] ]
    contents = "\nproject( {name} LANGUAGES {language} )\n".format(name = state['name'],
                                                                 language = project_language)
    subproject_languages_ = subproject_languages( state )
    if state['language'] in subproject_languages_:
        subproject_languages_.remove( state['language'] )
    
    for addition in subproject_languages_:
        contents += "enable_language( {0} )\n".format( language[addition] )
    
    return contents

def make_aux_directories( state ):
    contents = """
    if ( profile_generate )
        file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profiling" )
    endif()"""
    if state['language'] == 'fortran':
        contents += """

    file( MAKE_DIRECTORY ${fortran_module_directory} )
    set( fortran_module_directory "${CMAKE_BINARY_DIR}/modules" CACHE PATH "directory for fortran modules" ) """
        
    contents = textwrap.dedent( contents ) 
    return contents

def define_options( state ):
    if state['is_external_project']:
        return ''
    
    contents = """

    # general properties
    option( strict "Compile time warnings are converted to errors" {strict} )
    
    # binary instrumentation
    option( coverage "Enable binary instrumentation to collect test coverage information in the DEBUG configuration" )
    option( profile_generate "Enable binary instrumentation to generation execution profiles in the RELEASE configuration which may be used to guide later optimization" )
    
    # additional optimizations
    option( link_time_optimization "Enable link time optimization in the RELEASE configuration" )
    option( profile_use "In the RELEASE configuration, leverage previously generated exeution profile to inform optimization decisions" )
    option( nonportable_optimization "Enable optimizations which compromise portability of resulting binary in the RELEASE configuration" )
    
    # libraries and linking
    option( static "Statically link component and environment libraries")
    option( static_libaries "Statically link component libraries" ${{static}} )"""
    if has_library(state):
        contents += """
    option( static_{name} "Statically the {name} component library" ${{static_libaries}} ) """
        
    if has_tests(state):
        contents += """
        
    option( unit_tests "Compile the {name} unit tests and integrate with ctest" ON ) """
        
    contents += """
    
    if (profile_generate AND profile_use)
        message( FATAL_ERROR "Cannot both generate and use execution profile in the same configuration" )
    endif()""" 
    
    contents = textwrap.dedent(contents.format(name=state['name'], strict = 'ON' if state['strict'] else 'OFF'))
    return contents

def traverse_subprojects( state ):
    contents = ''
    if state['subprojects'] and not state['is_external_project']:
        contents += '\nget_directory_property( is_subproject PARENT_DIRECTORY )'
        build_queue = description.reconstruct_build_queue( state )
        build_queue.pop()
        for subproject in build_queue:
            contents += textwrap.dedent("""

            if( NOT TARGET {subproject} )
                add_subdirectory( ${{ROOT_DIRECTORY}}/subprojects/{subproject} )
            endif()""".format(subproject = subproject))
        contents += '\n'
    return contents

def define_compiler_flags( state ):
    contents = ""
    for compiler in state['compiler'].keys():        
        for operating_system in set(['linux','windows','osx']).intersection(state['compiler'][compiler].keys()):
            print(operating_system)
            environment = state['compiler'][compiler][operating_system]
            print(environment)
            flags = environment['flags']
            print(flags)
            args ={ 'name' : state['name'],
                    'vendor' : vendor[compiler],
                    'platform' : platform[operating_system],
                    'profile_path' : '${CMAKE_BINARY_DIR}/profiling' }
            
            args['common_flags'] = ' '.join( ['"{0}"'.format(flag) for flag in
                                              flags['common'] + flags['warning']  ] )
            if 'standard' in state:
                args['common_flags'] += ' "{0}"'.format( environment['standard'][ state['standard'] ] )
                                                
            args['debug_flags'] = ' '.join( ['"{0}"'.format(flag) for flag in flags['debug'] ] ) 
            args['release_flags'] = ' '.join( ['"{0}"'.format(flag) for flag in flags['optimization'] ] )
            
            keys = ['strict', 'coverage', 'profile generate',
                    'link time optimization', 'profile use',
                    'nonportable optimization', 'static' ]
            
            options = [ (key).replace(' ', '_') for key in keys ]
            
            for pair in zip(keys, options):
                args[ pair[1] + '_flags' ] = ' '.join( ['"{0}"'.format(flag) for flag in flags[pair[0]] ] )
                
            block = "\nset( {name}_{vendor}_{platform}_common_flags {common_flags} )"
            block += "\nset( {name}_{vendor}_{platform}_DEBUG_flags {debug_flags} )" 
            block += "\nset( {name}_{vendor}_{platform}_RELEASE_flags {release_flags} )"
            for option in options:
                block += "\nset( {{name}}_{{vendor}}_{{platform}}_{option}_flags {{{option}_flags}} )".format(option = option)
                
            contents += block.format(**args)
            
    return contents

def target_flags_expression( state ):
    contents=""
    for entry in state['compiler'].keys():
        template = " ${{{{{{name}}_{{vendor}}_${{{{CMAKE_SYSTEM_NAME}}}}_{0}_flags}}}}"
        common = template.format('common')
        debug = template.format('DEBUG')
        release = template.format('RELEASE')
        
        option_template = "\n         $< $<BOOL:${{{{{0}}}}}> : ${{{{{{name}}_{{vendor}}_${{{{CMAKE_SYSTEM_NAME}}}}_{0}_flags}}}} >"
        strict = option_template.format('strict')
        coverage = option_template.format('coverage')
        profile_generate = option_template.format('profile_generate')
        link_time_optimization = option_template.format('link_time_optimization')
        profile_use = option_template.format('profile_use')
        nonportable_optimization = option_template.format('nonportable_optimization')
        static = option_template.format('static') 
        
        addition = " $< $<{language}_COMPILER_ID:{vendor}> :" + common + strict + static \
                   + "\n     $< $<CONFIG:DEBUG> :" + debug + coverage + ' >' \
                   + "\n     $< $<CONFIG:RELEASE> :" + release + profile_generate + profile_use + link_time_optimization + nonportable_optimization + " > >\n"
        contents += addition.format(language = language[state['language']], name = state['name'], vendor = vendor[entry])
        
    contents += " ${{{language}_appended_flags}} ${{{name}_appended_flags}}".format( language=language[state['language']], name=state['name'] )
    return contents

def test_flags_expression( state ):
    contents = ''
    for entry in state['compiler'].keys():
        template = " ${{{{{{name}}_{{vendor}}_${{{{CMAKE_SYSTEM_NAME}}}}_{0}_flags}}}}"
        common = template.format('common')
        debug = template.format('DEBUG')
        release = template.format('RELEASE')
        
        option_template = "\n         $< $<BOOL:${{{{{0}}}}}> : ${{{{{{name}}_{{vendor}}_${{{{CMAKE_SYSTEM_NAME}}}}_{0}_flags}}}} >"
        strict = option_template.format('strict')
        link_time_optimization = option_template.format('link_time_optimization')
        nonportable_optimization = option_template.format('nonportable_optimization')
        
        addition = " $< $<{language}_COMPILER_ID:{vendor}>:" + common + strict \
                   + "\n     $< $<CONFIG:DEBUG> :" + debug + ' >' \
                   + "\n     $< $<CONFIG:RELEASE> :" + release + link_time_optimization + nonportable_optimization + " > >\n"
        contents += addition.format(language = language[state['language']], name = state['name'], vendor = vendor[entry])
        
    contents += " ${{{language}_appended_flags}} ${{{name}_appended_flags}}".format( language=language[state['language']], name=state['name'] )
    return contents

def set_library_type( state ):
    if not has_library( state ):
        return ''
    
    return textwrap.dedent( """

    if ( static_{name} )
        set( {name}_library_linkage STATIC )
    else ()
        set( {name}_library_linkage SHARED )
    endif () """ ).format(name = state['name'])
    
def collect_revision_info( state ):
  return textwrap.dedent("""

    if ( NOT GIT_EXECUTABLE )
        find_package( Git )
        if ( NOT GIT_FOUND )
            message( FATAL_ERROR "git installation was not found." )
        endif()
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    ) """)
     
def print_banner( state ):
    return textwrap.dedent("""

        message( STATUS "" )
        message( STATUS "-----------------------------------------------------------" )
        message( STATUS "" )
        message( STATUS "{name}" )
        message( STATUS "Git current branch: ${{GIT_BRANCH}}" )
        message( STATUS "Git commit hash: ${{GIT_HASH}}" )
        message( STATUS "" )
        message( STATUS "-----------------------------------------------------------" )
        """.format(**state))

def link_dependencies( state ):
    contents = ''
    name = state['name']
    if len( state['subprojects'] ) > 0 :
        contents += '\ntarget_link_libraries( {}'
        for name, subproject in state['subprojects'].items():
            contents += (' INTERFACE {}' if has_library( subproject ) else ' PUBLIC {}').format(name)
            
        contents += ' ) \n'
    return contents

def add_targets( state ):
    sources = []
    for group in state['file extension']:
        sources.extend(state[group])
        
    sources = '\n                    '.join( sources )
    target = "${{{name}_library_linkage}}" if has_library(state) else "INTERFACE"
    policy = "PUBLIC" if has_library(state) else "INTERFACE"
    expression = target_flags_expression(state)
    contents = """
    
    add_library( {name} {target} )
    target_sources( name {policy} {sources} )
    """
    if state['language'] == 'fortran':
        contents += "target_include_directory( {name} PUBLIC ${{fortran_module_directory}} )\n"
        
    if state['include path']:
        contents += "target_include_directory( {name} {policy} {include_path} )\n"
        
    if has_library( state ):
        contents += textwrap.dedent(
        """target_compile_options( {name} PRIVATE 
                                   {expression} )
        """)

    contents += link_dependencies( state )
        
    if has_executable( state ):
        contents += textwrap.dedent(
        """
        if ( NOT is_subproject )
            add_executable( {name}_executable {driver} )
            set_target_prooperties( {name}_executable PROPERTIES OUTPUT_NAME {name} )
            target_link_libraries( {name}_executable {policy} {name} )
            target_compile_options( {name}_executable PRIVATE 
                                    {expression} )
        endif()
        """)
        
    return textwrap.dedent( contents.format( name = state['name'],
                                             driver = ( state['driver'] if 'driver' in state else ''),
                                             target = target,
                                             policy = policy,
                                             sources = sources,
                                             expression = expression,
                                             include_path = state['include path'] ) )

def add_tests( state ):
    contents = ''
    if not state['is_external_project']:
        if has_tests( state ):
            name = state['name']
            contents = """ 
            if( NOT is_subproject )
                enable_testing() """
            if state['tests']:
                expression = test_flags_expression(state)
                contents += """
                if (unit_tests)"""
                for test_name, sources in state['tests'].items():
                    executable_name = test_name + '.test'
                    directory = sources[0]
                    contents += """
                    add_subdirectory( {} )""".format(directory)
                    test_contents = """
                    add_executable( {executable_name}"""
                    split = '\n                ' if len(sources) > 1 else ' '
                    test_contents += split + split.join( [ os.path.split( entry )[1] for entry in sources ] ) + ' )'
                    test_contents += """
                    target_compile_options( {executable_name} PRIVATE {expression} )
                    target_link_libraries( {executable_name} PUBLIC {name} ) """
                    if os.path.isdir( os.path.join( directory, 'resources' ) ):
                        test_contents += """
                    file( GLOB resources "resources/*" )"""
                        test_contents += """
                    foreach( resource ${{resources}} )"""
                        test_contents += """
                        file( COPY "${{resource}}" DESTINATION "${{CMAKE_CURRENT_BINARY_DIR}}" )"""
                        test_contents += """
                    endforeach()"""
                        
                    test_contents += """
                    add_test( NAME {test_name} COMMAND {executable_name} )"""
                    test_contents = textwrap.dedent(test_contents).format( executable_name = executable_name,
                                                                           test_name = test_name,
                                                                           expression = expression,
                                                                           name = name )
                    with open( os.path.join( directory, 'CMakeLists.txt' ), 'w') as TestCMakeFile:
                        TestCMakeFile.write( test_contents )
                    
                contents += """
                endif() """
            if os.path.isdir( os.path.join( os.getcwd(), 'test' ) ) :
                contents += """
                add_subdirectory( test )"""
          
            contents += """
            endif()"""
            contents = textwrap.dedent(contents)
            
    return contents

def is_subdirectory( child, parent ):
    true_child = os.realpath( child )
    true_parent = os.realpath( parent )
    while true_child:
        true_child = os.path.split(true_child)[0]
        if true_child == true_parent:
            return True

    return False

def install( state ):
    targets = []
    contents = ""
    if has_library( state ):
        targets.append( "{name}" )
        
    if has_executable( state ):
        targets.append( " {name}_executable " )
        
    if targets:
        targets = ' '.join(targets)
        contents += """

        install( TARGETS {targets} 
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib
                 PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
                             GROUP_EXECUTE GROUP_READ 
                             WORLD_EXECUTE WORLD_READ """
        if "group id" in state:
            contents += """
                             SETGID {gid}"""
            
        contents += """
        )
        
        """
        contents = contents.format( name = state['name'],
                                    targets = targets,
                                    gid = state['group id'] if 'group id' in state else '' )
        
    if 'include path' in state and is_subdirectory( state['include path'], os.cwd() ):
        if 'header files' in state['file extension']:
            regex = []
            for extension in state['file extension']['header files']:
                regex.append( ".*\.{0}".format(extension) )
                
            regex = '|'.join( regex )
                
            contents += """

        install( DIRECTORY {include_path} DESTINATION include
                 FILES_MATCHING REGEX "{regex}"
                 PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
                             GROUP_EXECUTE GROUP_READ 
                             WORLD_EXECUTE WORLD_READ"""
            if "group id" in state:
                contents += """
                SETGID {gid}"""
            
            contents += """ )
                
                """
                    
        contents = contents.format( include_path = state['include_path'],
                                    regex = regex,
                                    gid = state['group id'] if 'group id' in state else '' )
        
    if state['language'] == 'fortran':
        contents += """

        install( DIRECTORY ${{fortran_module_directory}} DESTINATION include
                 PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
                             GROUP_EXECUTE GROUP_READ 
                             WORLD_EXECUTE WORLD_READ"""
        if "group id" in state:
            contents += """
                             SETGID {gid}"""
        
        contents += """ )
        
        """
        contents = contents.format( gid = state['group id'] if 'group id' in state else '' )
        
    contents = textwrap.dedent(contents)
    return contents

def generate():
    state = description.deserialize()
    description.collect_subprojects( state, state['project path'] )
    contents = "cmake_minimum_required( VERSION 3.2 ) \n"
    contents += fetch_subprojects( state )
    contents += project_statement( state )
    contents += make_aux_directories( state )
    contents += define_options( state )    
    contents += traverse_subprojects( state )
    contents += define_compiler_flags( state )
    contents += set_library_type( state )
    contents += collect_revision_info( state )
    contents += print_banner( state )
    contents += add_targets( state )
    contents += add_tests( state )
    contents += install( state )    
    with open('CMakeLists.txt', 'w') as CMakeFile:
        CMakeFile.write(contents)
