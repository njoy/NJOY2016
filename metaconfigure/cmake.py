import os
import textwrap
from . import description

language = {'c' : 'C', 'c++' : 'CXX', 'fortran' : 'Fortran'}
platform = {'linux':'Linux', 'osx':'Darwin', 'windows':'Windows'}
vendor = {'gcc' : 'GNU',
          'g++' : 'GNU',
          'gfortran' : 'GNU',
          'llvm clang' : 'Clang',
          'llvm clang++' : 'Clang',
          'apple clang' : 'AppleClang',
          'apple clang++' : 'AppleClang'}


def fetch_subprojects(state):
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


def subproject_languages(state):
    languages=set()
    if 'subprojects' in state:
        for subproject in state['subprojects'].values():
            languages.add(subproject['language'])
            languages=languages | subproject_languages(subproject)

    return languages


def has_library(state):
    return bool(state['implementation files'])


def has_executable(state):
    return 'driver' in state 


def has_unit_tests(state):
    return 'tests' in state and state['tests']


def has_tests(state):
    return has_unit_tests(state) or os.path.isdir("tests") 


def project_statement(state):
    project_language = language[state['language']]
    contents = "\nproject( {name} LANGUAGES {language} )".format(name=state['name'],
                                                               language=project_language)
    contents += "\nget_directory_property( is_subproject PARENT_DIRECTORY )"
    contents += "\ninclude( CMakeDependentOption REQUIRED )"
    additional_languages = subproject_languages(state)
    if state['language'] in additional_languages:
        additional_languages.remove(state['language'])

    for addition in additional_languages:
        contents += "\nenable_language( {0} )".format(language[addition])

    return contents


def compiler_minimum(state):
    contents=''
    any_=False
    for name, compiler in state['compiler'].items():
        if 'minimum version' in compiler:
            any_=True
            contents += "\nset( {name}_{vendor}_minimum_version {version} )".format(name=state['name'],
                                                                                    vendor=vendor[name],
                                                                                    version=compiler['minimum version'])

    if any_:
        contents += """

if( {name}_${{CMAKE_{language}_COMPILER_ID}}_minimum_version )
    if( CMAKE_{language}_COMPILER_VERSION AND
        CMAKE_{language}_COMPILER_VERSION VERSION_LESS
        ${{{name}_${{CMAKE_{language}_COMPILER_ID}}_minimum_version}} )
        message( FATAL_ERROR "${{CMAKE_{language}_COMPILER_ID}} version must be greater than ${{{name}_${{CMAKE_{language}_COMPILER_ID}}_minimum_version}}" )
    endif()
endif()""".format(language=language[state['language']],
                  name=state['name'])

    return '\n' + contents


def make_aux_directories(state):
    contents="""

    if ( profile_generate )
        file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/profiling" )
    endif()"""
    if state['language'] == 'fortran':
        contents += """

    set( CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules" CACHE PATH "directory for fortran modules" )
    file( MAKE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}" ) """

    contents = textwrap.dedent(contents) 
    return contents


def define_options(state):
    if state['is external project']:
        return ''

    contents="""

    # general properties
    option( {name}_strict "Compile time warnings are converted to errors" {strict} )
    
    # binary instrumentation
    option( coverage "Enable binary instrumentation to collect test coverage information in the DEBUG configuration" )
    option( profile_generate "Enable binary instrumentation to generation execution profiles in the RELEASE configuration which may be used to guide later optimization" )
    
    # additional optimizations
    option( link_time_optimization "Enable link time optimization in the RELEASE configuration" )
    option( profile_use "In the RELEASE configuration, leverage previously generated exeution profile to inform optimization decisions" )
    option( nonportable_optimization "Enable optimizations which compromise portability of resulting binary in the RELEASE configuration" )
    
    # libraries and linking
    option( static "Statically link component and environment libraries" OFF )
    if ( static AND ( "${{CMAKE_SYSTEM_NAME}}" STREQUAL "Darwin" ) )
        message( FATAL_ERROR "Static binaries not supported on OSX" )
    endif()

    CMAKE_DEPENDENT_OPTION( static_libraries "Statically link component libraries" OFF "NOT static" ON )"""

    if has_library(state):
        contents += """
    CMAKE_DEPENDENT_OPTION( static_{name} "Statically link the {name} component library" OFF "NOT static;NOT static_libraries" ON ) """

    if has_unit_tests(state):
        contents += """
        
    option( unit_tests "Compile the {name} unit tests and integrate with ctest" ON ) """

    contents += """
    
    if ( profile_generate AND profile_use )
        message( FATAL_ERROR "Cannot both generate and use execution profile in the same configuration" )
    endif()""" 

    contents=textwrap.dedent(contents.format(name=state['name'], strict='ON' if state['strict'] else 'OFF'))
    return contents


def traverse_subprojects(state):
    contents=''
    if state['subprojects'] and not state['is external project']:
        contents += '\nget_directory_property( is_subproject PARENT_DIRECTORY )'
        build_queue = description.reconstruct_build_queue(state)
        build_queue.pop()
        for subproject in build_queue:
            contents += textwrap.dedent("""

            if( NOT TARGET {subproject} )
                add_subdirectory( ${{ROOT_DIRECTORY}}/subprojects/{subproject} )
            endif()""".format(subproject=subproject))

        contents += '\n'

    return contents


def define_compiler_flags(state):
    contents="\n"
    for compiler in state['compiler'].keys():        
        for operating_system in set(['linux','windows','osx']).intersection(state['compiler'][compiler].keys()):
            environment=state['compiler'][compiler][operating_system]
            flags=environment['flags']
            args ={ 'name' : state['name'],
                    'vendor' : vendor[compiler],
                    'platform' : platform[operating_system],
                    'profile_path' : '${CMAKE_BINARY_DIR}/profiling' }

            args['common_flags']=' '.join(['"{0}"'.format(flag) for flag in
                                              flags['common'] + flags['warning'] ])
            if 'standard' in state:
                args['common_flags'] += ' "{0}"'.format(environment['standard'][state['standard']])

            args['debug_flags']=' '.join(['"{0}"'.format(flag) for flag in flags['debug']]) 
            args['release_flags']=' '.join(['"{0}"'.format(flag) for flag in flags['optimization']])

            keys=['strict', 'coverage', 'subproject', 'base project',
                    'profile generate', 'link time optimization', 'profile use',
                    'nonportable optimization', 'static']

            options=[(key).replace(' ', '_') for key in keys]

            for pair in zip(keys, options):
                args[pair[1] + '_flags']=' '.join(['"{0}"'.format(flag) for flag in flags[pair[0]]])

            block="\nset( {name}_{vendor}_{platform}_common_flags {common_flags} )"
            block += "\nset( {name}_{vendor}_{platform}_DEBUG_flags {debug_flags} )" 
            block += "\nset( {name}_{vendor}_{platform}_RELEASE_flags {release_flags} )"

            for option in options:
                block += "\nset( {{name}}_{{vendor}}_{{platform}}_{option}_flags {{{option}_flags}} )".format(option=option)

            contents += block.format(**args).format(**args) 

    return contents


def lto_flags_expression(state):
    contents = ""
    release = "${{${{PREFIX}}_RELEASE_flags}};".format(language=language[state['language']], name=state['name'])
    link_time_optimization = "${{${{PREFIX}}_link_time_optimization_flags}}".format(language=language[state['language']], name=state['name'])
    option_template = "$<$<BOOL:${{{{{0}}}}}>:${{{{${{{{PREFIX}}}}_{0}_flags}}}};>"

    profile_generate = option_template.format('profile_generate').format(language=language[state['language']], name=state['name'])
    profile_use = option_template.format('profile_use').format(language=language[state['language']], name=state['name'])
    nonportable_optimization = option_template.format('nonportable_optimization').format(language=language[state['language']], name=state['name'])
    coverage = option_template.format('coverage').format(language=language[state['language']], name=state['name'])
    language_appended_flags = "$<$<BOOL:{0}_appended_flags>:${{{0}_appended_flags}};>".format(language[state['language']])
    project_appended_flags = "$<$<BOOL:{0}_appended_flags>:${{{0}_appended_flags}};>".format(state['name'])
    contents = "\"$<$<AND:$<CONFIG:RELEASE>,$<BOOL:${{link_time_optimization}}>>:{release}{link_time_optimization}{profile_generate}{profile_use}{nonportable_optimization}>$<$<CONFIG:DEBUG>:{coverage}>{language_appended_flags}{project_appended_flags}\""
    contents = contents.format(release=release,
                               link_time_optimization=link_time_optimization,
                               profile_generate=profile_generate,
                               profile_use=profile_use,
                               nonportable_optimization=nonportable_optimization,
                               coverage=coverage,
                               language_appended_flags=language_appended_flags,
                               project_appended_flags=project_appended_flags)

    return contents

    
def target_flags_expression(state):
    contents = ""
    template = "\n${{{{${{{{PREFIX}}}}_{0}_flags}}}}"
    common = template.format('common')
    debug = template.format('DEBUG')
    release = template.format('RELEASE')

    option_template = "\n$<$<BOOL:${{{{{0}}}}}>:${{{{${{{{PREFIX}}}}_{0}_flags}}}}>"
    strict = option_template.format('{name}_strict')
    coverage = option_template.format('coverage')
    profile_generate = option_template.format('profile_generate')
    link_time_optimization = option_template.format('link_time_optimization')
    profile_use = option_template.format('profile_use')
    nonportable_optimization = option_template.format('nonportable_optimization')
    static = option_template.format('static')

    subproject = "\n$<$<BOOL:${{is_subproject}}>:${{${{PREFIX}}_subproject_flags}}>"
    base_project = "\n$<$<NOT:$<BOOL:${{is_subproject}}>>:${{${{PREFIX}}_base_project_flags}}>"
    addition = common + strict + static + subproject + base_project \
                + "\n$<$<CONFIG:DEBUG>:" + debug + ' ' + coverage + '>' \
                + "\n$<$<CONFIG:RELEASE>:"\
                + release + ' '\
                + profile_generate + ' '\
                + profile_use + ' '\
                + link_time_optimization\
                + ' ' + nonportable_optimization + ">"
    contents += addition.format(language=language[state['language']], name=state['name'])        
    contents += "\n${{{language}_appended_flags}} ${{{name}_appended_flags}}".format(language=language[state['language']], name=state['name'])
    return contents


def test_flags_expression(state):
    contents = ''
    template = "${{{{${{{{PREFIX}}}}_{0}_flags}}}}"
    common = template.format('common')
    debug = template.format('DEBUG')
    release = template.format('RELEASE')
        
    option_template = "\n$<$<BOOL:${{{{{0}}}}}>:${{{{${{{{PREFIX}}}}_{0}_flags}}}}>"
    coverage = option_template.format('coverage')
    strict = option_template.format('strict')
    link_time_optimization = option_template.format('link_time_optimization')
    nonportable_optimization = option_template.format('nonportable_optimization')
        
    addition = common + strict \
                + "$<$<CONFIG:DEBUG>:\n" + debug + coverage + '>' \
                + "\n$<$<CONFIG:RELEASE>:\n" + release + link_time_optimization + nonportable_optimization + ">\n"
    contents += addition.format(language=language[state['language']], name=state['name'])
        
    contents += "\n${{{language}_appended_flags}} ${{{name}_appended_flags}}".format(language=language[state['language']], name=state['name'])
    return contents


def set_library_type(state):
    if not has_library(state):
        return ''
    
    return textwrap.dedent("""

    if ( static_{name} )
        set( {name}_library_linkage STATIC )
    else ()
        set( {name}_library_linkage SHARED )
    endif () 

    set( CMAKE_SKIP_BUILD_RPATH FALSE )
    set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
    if ( CMAKE_SYSTEM_NAME STREQUAL "Darwin" )
        set( rpath_prefix "@loader_path" )
    else()
        set( rpath_prefix "\\\\$ORIGIN" )
    endif()
    list( INSERT 0 CMAKE_INSTALL_RPATH "${{rpath_prefix}}/../lib" )
    set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )""").format(name=state['name'])


def collect_revision_info(state):
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

     
def print_banner(state):
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


def link_dependencies(state):
    contents=''
    name=state['name']
    if len(state['subprojects']) > 0 :
        contents += '\ntarget_link_libraries( {name}'
        for name, subproject in state['subprojects'].items():
            contents += (' PUBLIC {}' if has_library(state) else ' INTERFACE {}').format(name)

        contents += ' )\n'

    return contents


def add_targets(state):
    sources=[]
    for group in state['file extension']:
        sources.extend(state[group])
        
    sources='"\n             "${CMAKE_CURRENT_SOURCE_DIR}/'.join(sources)
    policy="PUBLIC" if has_library(state) else "INTERFACE"
    compile_flags=target_flags_expression(state)
    link_flags=lto_flags_expression(state)
    
    if has_library(state):
        contents="""
add_library( {name} ${{{name}_library_linkage}} 
             "${{CMAKE_CURRENT_SOURCE_DIR}}/{sources}" )
        """
    else:
        contents="""
add_library( {name} INTERFACE )
target_sources( {name} INTERFACE "${{CMAKE_CURRENT_SOURCE_DIR}}/{sources}" )
        """

    if state['language'] == 'fortran':
        contents += """
target_include_directories( {name} PUBLIC "${{CMAKE_Fortran_MODULE_DIRECTORY}}" )
        """
        
    if 'include path' in state:
        contents += """
target_include_directories( {name} {policy} {include_path} )
        """

    if has_library(state) or has_executable(state) or has_tests(state):
        contents += """
set( PREFIX {name}_${{CMAKE_{language}_COMPILER_ID}}_${{CMAKE_SYSTEM_NAME}} )
        """
        
    if has_library(state):
        contents += """
target_compile_options( {name} PRIVATE {compile_flags} )
        """

    contents += """
target_link_libraries( {name} {policy} {link_flags} )
    """
    contents += link_dependencies(state)
        
    if has_executable(state):
        contents += textwrap.dedent(
        """
if ( NOT is_subproject )
    add_executable( {name}_executable {driver} )
    set_target_properties( {name}_executable PROPERTIES OUTPUT_NAME {name} )
    target_compile_options( {name}_executable PRIVATE {indented_compile_flags} )
    target_link_libraries( {name}_executable {policy} {name} )
endif()
        """)
        
    return textwrap.dedent(contents.format(name=state['name'],
                                           driver=(state['driver'] if 'driver' in state else ''),
                                           language=language[state['language']],
                                           policy=policy,
                                           sources=sources,
                                           compile_flags=compile_flags,
                                           indented_compile_flags=compile_flags.replace('\n', '\n    '),
                                           link_flags=link_flags,
                                           include_path=state['include path'] if 'include path' in state else ''))


def add_tests(state):
    contents=''
    if not state['is external project']:
        if has_tests(state):
            name=state['name']
            contents=""" 
            if( NOT is_subproject )
                enable_testing() """
            if state['tests']:
                expression=test_flags_expression(state)
                contents += """
                if ( unit_tests )"""
                for test_name, sources in state['tests'].items():
                    executable_name=test_name + '.test'
                    directory=os.path.dirname(sources[0])
                    contents += """
                    add_subdirectory( {} )""".format(directory)
                    test_contents="""
add_executable( {executable_name}"""
                    split='\n                ' if len(sources) > 1 else ' '
                    test_contents += split + split.join([os.path.split(entry)[1] for entry in sources]) + ' )'
                    test_contents += """
target_compile_options( {executable_name} PRIVATE {expression} )
target_link_libraries( {executable_name} PUBLIC {name} ) """
                    if os.path.isdir(os.path.join(directory, 'resources')):
                        test_contents += """
file( GLOB resources "resources/*" )"""
                        test_contents += """
foreach( resource ${{resources}})"""
                        test_contents += """
    file( COPY "${{resource}}" DESTINATION "${{CMAKE_CURRENT_BINARY_DIR}}" )"""
                        test_contents += """
endforeach()"""
                        
                    test_contents += """
add_test( NAME {test_name} COMMAND {executable_name} )"""
                    test_contents=textwrap.dedent(test_contents).format(executable_name=executable_name,
                                                                        test_name=test_name,
                                                                        expression=expression,
                                                                           name=name)
                    with open(os.path.join(directory, 'CMakeLists.txt'), 'w') as TestCMakeFile:
                        TestCMakeFile.write(test_contents)
                    
                contents += """
                endif() """

            if os.path.isdir(os.path.join(os.getcwd(), 'tests')) :
                contents += """
                add_subdirectory( tests )"""

            contents += """
            endif()"""
            contents=textwrap.dedent(contents)
            
    return contents


def is_subdirectory(child, parent):
    true_child=os.path.realpath(child)
    true_parent=os.path.realpath(parent)
    while true_child:
        true_child=os.path.split(true_child)[0]
        if true_child == true_parent:
            return True

    return False


def install(state):
    targets=[]
    contents="\n"
    if has_library(state):
        targets.append("{name}".format(name=state['name']))

    if has_executable(state):
        targets.append("{name}_executable".format(name=state['name']))

    if targets:
        block = """
        install( TARGETS ${{installation_targets}} 
                 RUNTIME DESTINATION bin
                 LIBRARY DESTINATION lib
                 ARCHIVE DESTINATION lib
                 PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
                             GROUP_EXECUTE GROUP_READ 
                             WORLD_EXECUTE WORLD_READ"""
        if "group id" in state:
            block += """
                             SETGID {gid}"""

        block += """ )
        """
        
        if has_executable(state):
            if len(targets) > 1:
                contents += """
        set( installation_targets {0} )""".format(targets[0])
                contents += """
        if ( NOT is_subproject )
            list( APPEND installation_targets {0} )
        endif()
                """.format(targets[-1])
                contents += block
            else:
                contents += """
        if ( NOT is_subproject )
            list( APPEND installation_targets {0} )"""
                contents += block.replace('\n', '\n    ')
                contents += """
        endif()
                """
                


    regex=[]
    if 'include path' in state and is_subdirectory(state['include path'], os.getcwd()):
        if 'header files' in state['file extension']:
            for extension in state['file extension']['header files']:
                regex.append(".*\\.{0}$".format(extension).replace('+', '[+]'))

            contents += """
        install( DIRECTORY {include_path}/ DESTINATION include
                 FILE_PERMISSIONS OWNER_READ OWNER_WRITE 
                                  GROUP_READ 
                                  WORLD_READ
                 DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE 
                                       GROUP_READ 
                                       WORLD_READ"""
            if "group id" in state:
                contents += """
                SETGID {gid}"""


            contents +=     """
                 FILES_MATCHING REGEX "{regex}" """

            contents += """ )
                """

    if state['language'] == 'fortran':
        contents += """
        file( RELATIVE_PATH relative_fortran_module_files_path 
              "${{CMAKE_CURRENT_SOURCE_DIR}}" "${{CMAKE_Fortran_MODULE_DIRECTORY}}" )
        file( GLOB fortran_module_files 
              RELATIVE "${{relative_fortran_module_files_path}}"
              *.mod )
        install( FILES ${{fortran_module_files}} 
                 DESTINATION include
                 PERMISSIONS OWNER_READ OWNER_WRITE 
                             GROUP_READ 
                             WORLD_READ"""
        if "group id" in state:
            contents += """
                             SETGID {gid}"""

        contents += """ )
        """

    contents = contents.format(name=state['name'],
                               targets=targets,
                               include_path=state['include path'] if 'include path' in state else '',
                               regex='|'.join(regex),
                               gid=state['group id'] if 'group id' in state else '')

    return textwrap.dedent(contents)

def generate():
    state=description.deserialize()
    description.collect_subprojects(state, state['project path'])
    contents = \
        """
cmake_minimum_required( VERSION 3.2 ) 
set( CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Supported configuration types" FORCE )
        """
    contents += fetch_subprojects(state)
    contents += project_statement(state)
    if not state['is external project']:
        contents += compiler_minimum(state)
        contents += define_options(state)
        contents += make_aux_directories(state)
        contents += define_compiler_flags(state)
        contents += set_library_type(state)

    contents += traverse_subprojects(state)
    contents += collect_revision_info(state)
    contents += print_banner(state)
    if not state['is external project']:
        contents += add_targets(state)
        contents += add_tests(state) 
        contents += install(state)
        
    with open('CMakeLists.txt', 'w') as CMakeFile:
        CMakeFile.write(contents)
