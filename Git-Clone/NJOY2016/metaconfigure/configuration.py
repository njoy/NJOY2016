import copy

languages = {'c' : {}, 'c++' : {}, 'fortran' : {} }

for langauge in languages.keys():
    languages[langauge]['file extension'] = {}
    languages[langauge]['compiler'] = {}

languages['c']['standards'] = ['c89', 'c99', 'c11']
languages['c++']['standards'] = ['c++98', 'c++11', 'c++14']
languages['fortran']['standards'] = ['fortran77', 'fortran90', 'fortran95', 'fortran2003', 'fortran2008']

languages['c']['file extension']['implementation files'] = ['c']
languages['c']['file extension']['header files'] = ['h']
languages['c++']['file extension']['implementation files'] = ['c++', 'cxx', 'cpp', 'cc']
languages['c++']['file extension']['header files'] = ['h++', 'hxx', 'hpp', 'hh', 'h']
languages['fortran']['file extension']['implementation files'] = ['f', 'for', 'f90']
    
languages['c']['compiler']['gcc'] = {}
languages['c']['compiler']['llvm clang'] = {}
languages['c']['compiler']['apple clang'] = {}
languages['c++']['compiler']['g++'] = {}
languages['c++']['compiler']['llvm clang++'] = {}
languages['c++']['compiler']['apple clang++'] = {}
languages['fortran']['compiler']['gfortran'] = {}

languages['c']['compiler']['gcc']['linux'] = {}
languages['c']['compiler']['gcc']['osx'] = languages['c']['compiler']['gcc']['linux']
languages['c']['compiler']['gcc']['windows'] = languages['c']['compiler']['gcc']['linux']
languages['c']['compiler']['llvm clang']['linux'] = languages['c']['compiler']['gcc']['linux']
languages['c']['compiler']['llvm clang']['osx'] = languages['c']['compiler']['gcc']['linux']
languages['c']['compiler']['llvm clang']['windows'] = languages['c']['compiler']['gcc']['linux']
languages['c']['compiler']['apple clang']['osx'] = languages['c']['compiler']['gcc']['linux']

languages['fortran']['compiler']['gfortran']['linux'] = {}
languages['fortran']['compiler']['gfortran']['osx'] = languages['fortran']['compiler']['gfortran']['linux']
languages['fortran']['compiler']['gfortran']['windows'] = languages['fortran']['compiler']['gfortran']['linux']

languages['c']['compiler']['gcc']['linux']['standard'] = {'c89' : '-std=c90',
                                                          'c99' : '-std=c99',
                                                          'c11' : '-std=c11'}

languages['c++']['compiler']['g++']['linux'] = {}
languages['c++']['compiler']['g++']['linux']['standard'] = {'c++98' : '-std=c++98',
                                                            'c++11' : '-std=c++11',
                                                            'c++14' : '-std=c++14'}

languages['fortran']['compiler']['gfortran']['linux']['standard'] = {'fortran77' : '-std=legacy' ,
                                                                     'fortran95' : '-std=f95',
                                                                     'fortran2003' : '-std=f2003',
                                                                     'fortran2008' : '-std=f2008' }

languages['c']['compiler']['gcc']['linux']['flags'] = {}
languages['fortran']['compiler']['gfortran']['linux']['flags'] = {}

languages['c']['compiler']['gcc']['linux']['flags']['common'] = []
languages['fortran']['compiler']['gfortran']['linux']['flags']['common'] = languages['c']['compiler']['gcc']['linux']['flags']['common']

languages['c']['compiler']['gcc']['linux']['flags']['subproject'] = []
languages['fortran']['compiler']['gfortran']['linux']['flags']['subproject'] = languages['c']['compiler']['gcc']['linux']['flags']['subproject']

languages['c']['compiler']['gcc']['linux']['flags']['base project'] = []
languages['fortran']['compiler']['gfortran']['linux']['flags']['base project'] = languages['c']['compiler']['gcc']['linux']['flags']['subproject']

languages['c']['compiler']['gcc']['linux']['flags']['warning'] = ['-Wall', '-Wextra', '-Wpedantic']
languages['fortran']['compiler']['gfortran']['linux']['flags']['warning'] = languages['c']['compiler']['gcc']['linux']['flags']['warning']

languages['c']['compiler']['gcc']['linux']['flags']['strict'] = ['-Werror']
languages['fortran']['compiler']['gfortran']['linux']['flags']['strict'] = languages['c']['compiler']['gcc']['linux']['flags']['strict']

languages['c']['compiler']['gcc']['linux']['flags']['debug'] = ['-O0', '-g', '-gdwarf-3', '-fsignaling-nans']
languages['fortran']['compiler']['gfortran']['linux']['flags']['debug'] = languages['c']['compiler']['gcc']['linux']['flags']['debug'] + ['-fcheck=all', '-ffpe-trap=invalid,zero,overflow']

languages['c']['compiler']['gcc']['linux']['flags']['coverage'] = ['--coverage']
languages['fortran']['compiler']['gfortran']['linux']['flags']['coverage'] = languages['c']['compiler']['gcc']['linux']['flags']['coverage']

languages['c']['compiler']['gcc']['linux']['flags']['static'] = ['-static']
languages['fortran']['compiler']['gfortran']['linux']['flags']['static'] = languages['c']['compiler']['gcc']['linux']['flags']['static']

languages['c']['compiler']['gcc']['linux']['flags']['optimization'] = ['-O3', '-DNDEBUG']
languages['fortran']['compiler']['gfortran']['linux']['flags']['optimization'] = languages['c']['compiler']['gcc']['linux']['flags']['optimization']

languages['c']['compiler']['gcc']['linux']['flags']['nonportable optimization'] = ['-march=native']
languages['fortran']['compiler']['gfortran']['linux']['flags']['nonportable optimization'] = languages['c']['compiler']['gcc']['linux']['flags']['nonportable optimization']

languages['c']['compiler']['gcc']['linux']['flags']['link time optimization'] = ['-flto']
languages['fortran']['compiler']['gfortran']['linux']['flags']['link time optimization'] = languages['c']['compiler']['gcc']['linux']['flags']['link time optimization']

languages['c']['compiler']['gcc']['linux']['flags']['profile generate'] = ["-fprofile-generate='{profile_path}'"]
languages['fortran']['compiler']['gfortran']['linux']['flags']['profile generate'] = languages['c']['compiler']['gcc']['linux']['flags']['profile generate']

languages['c']['compiler']['gcc']['linux']['flags']['profile use'] = ["-fprofile-use='{profile_path}'"]
languages['fortran']['compiler']['gfortran']['linux']['flags']['profile use'] = languages['c']['compiler']['gcc']['linux']['flags']['profile use']

languages['c++']['compiler']['g++']['linux']['flags'] = copy.deepcopy( languages['c']['compiler']['gcc']['linux']['flags'] )
languages['c++']['compiler']['g++']['osx'] = copy.deepcopy( languages['c++']['compiler']['g++']['linux'] )
languages['c++']['compiler']['g++']['windows'] = languages['c++']['compiler']['g++']['linux']

languages['c++']['compiler']['llvm clang++'] = copy.deepcopy( languages['c++']['compiler']['g++'])
languages['c++']['compiler']['llvm clang++']['linux']['flags']['common'] = ['-stdlib=libstdc++']
languages['c++']['compiler']['llvm clang++']['linux']['flags']['debug'].remove('-fsignaling-nans')

languages['c++']['compiler']['llvm clang++']['osx']['flags'] = copy.deepcopy(languages['c++']['compiler']['llvm clang++']['linux']['flags'])
languages['c++']['compiler']['llvm clang++']['osx']['flags']['common'] = ['-stdlib=libc++']
languages['c++']['compiler']['llvm clang++']['windows'] = languages['c++']['compiler']['llvm clang++']['osx']

languages['c++']['compiler']['apple clang++']['osx'] = languages['c++']['compiler']['llvm clang++']['osx']

languages['c']['compiler']['gcc']['minimum version'] = '6.2'
languages['c']['compiler']['llvm clang'] ['minimum version'] = '3.8'
languages['c']['compiler']['apple clang']['minimum version'] = '8.0'
languages['c++']['compiler']['g++']['minimum version'] = '6.2'
languages['c++']['compiler']['llvm clang++']['minimum version'] = '3.8'
languages['c++']['compiler']['apple clang++']['minimum version'] = '8.0'
languages['fortran']['compiler']['gfortran']['minimum version'] = '5.1'
