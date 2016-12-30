languages = [ 'c', 'c++', 'fortran' ]

implementation_extensions = { 'c' : ['c'],
                              'c++' : ['cpp', 'cxx', 'cc'],
                              'fortran' : ['f', 'for', 'f90'] }

header_extensions = { 'c' : ['h'],
                      'c++' : ['hpp', 'hxx', 'hh', 'h'],
                      'fortran' : [] }

language_versions = { 'c' : ['C89', 'C99', 'C11'],
                      'c++' : ['C++98', 'C++11', 'C++14'],
                      'fortran' : ['Fortran77', 'Fortran90', 'Fortran95', 'Fortran2003', 'Fortran2008'] }

compilers = { 'c' : [ 'gcc', 'clang', 'apple clang' ],
              'c++' : [ 'g++', 'clang++', 'apple clang++' ],
              'fortran' : [ 'gfortran' ] }

revision_flags = { 'gcc' : { 'C89' : '-std=c90',
                             'C99' : '-std=c99',
                             'C11' : '-std=c11' },
                   'clang' : { 'C89' : '-std=c90',
                               'C99' : '-std=c99',
                               'C11' : '-std=c11' },
                   'apple clang' : { 'C89' : '-std=c90',
                                     'C99' : '-std=c99',
                                     'C11' : '-std=c11' },
                   'g++' : { 'C++98' : '-std=c++98',
                             'C++11' : '-std=c++11',
                             'C++14' : '-std=c++14' },
                   'clang++' : { 'C++98' : '-std=c++98',
                                 'C++11' : '-std=c++11',
                                 'C++14' : '-std=c++14' },
                   'apple clang++' : { 'C++98' : '-std=c++98',
                                       'C++11' : '-std=c++11',
                                       'C++14' : '-std=c++14' },
                   'gfortran' : { 'Fortran77' : '-std=legacy' ,
                                  'Fortran95' : '-std=f95',
                                  'Fortran2003' : '-std=f2003',
                                  'Fortran2008' : '-std=f2008' } }

def revision_flag( language, revision, compiler ):
    assert compiler in compilers[language]
    assert revision in language_versions[language]
    trial_revision = revision
    while( True ):
        try:
            return revision_flags[compiler][trial_revision]
        except KeyError:
            try:
                index = language_versions[language].index(trial_revision) + 1
                trial_revision = language_versions[language][index]
            except IndexError:
                raise RuntimeError(
                    '{} does not support {} version {} or later'.format(
                          compiler, language, revision) )

for language in languages:
    assert( language in compilers )

warning_flags = { 'gcc' : '-Wall -Wextra -Wpedantic -Werror',
                  'g++' : '-Wall -Wextra -Wpedantic -Werror',
                  'clang' : '-Wall -Wextra -Wpedantic -Werror',
                  'clang++' : '-Wall -Wextra -Wpedantic -Werror',
                  'apple clang' : '-Wall -Wextra -Wpedantic -Werror',
                  'apple clang++' : '-Wall -Wextra -Wpedantic -Werror',
                  'gfortran' : '-Wall -Wextra -Wpedantic -Werror' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in warning_flags )

debug_flags = { 'gcc' : '-g -gdwarf-3',
                'g++' : '-g -gdwarf-3',
                'clang' : '-g -gdwarf-3',
                'clang++' : '-g -gdwarf-3',
                'apple clang' : '-g -gdwarf-3',
                'apple clang++' : '-g -gdwarf-3',
                'gfortran' : '-g -gdwarf-3' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in debug_flags )

coverage_flags = { 'gcc' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'g++' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'clang' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'clang++' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'apple clang' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'apple clang++' : '-fprofile-arcs -ftest-coverage -fno-inline',
                   'gfortran' : '-fprofile-arcs -ftest-coverage -fno-inline' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in coverage_flags )

optimization_flags = { 'gcc' : '-O3 -DNDEBUG',
                       'g++' : '-O3 -DNDEBUG',
                       'clang' : '-O3 -DNDEBUG',
                       'clang++' : '-O3 -DNDEBUG',
                       'apple clang' : '-O3 -DNDEBUG',
                       'apple clang++' : '-O3 -DNDEBUG',
                       'gfortran' : '-O3 -DNDEBUG' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in optimization_flags )

native_flags = { 'gcc' : '-march=native',
                 'g++' : '-march=native',
                 'clang' : '-march=native',
                 'clang++' : '-march=native',
                 'apple clang' : '-march=native',
                 'apple clang++' : '-march=native',
                 'gfortran' : '-march=native' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in optimization_flags )

link_time_optimization_flags = { 'gcc' : '-flto',
                                 'g++' : '-flto',
                                 'clang' : '-flto',
                                 'clang++' : '-flto',
                                 'apple clang' : '-flto',
                                 'apple clang++' : '-flto',
                                 'gfortran' : '-flto' }

for language, compiler_list in compilers.items():
    for compiler in compiler_list:
        assert( compiler in link_time_optimization_flags )
