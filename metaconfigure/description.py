"""
A script to generate a serialized project description for build script
generation.
"""
import os
import glob
import json
import re
from metaconfigure import configuration


def serialize(state):
    with open("metaconfigure/description.json", "w") as json_file:
        project_path = state.pop('project path')
        subprojects = state.pop('subprojects')
        json_file.write(json.dumps(state, indent=2, sort_keys=True))
        state['subprojects'] = subprojects
        state['project path'] = project_path


def deserialize():
    with open("metaconfigure/description.json", "r") as json_file:
        state = json.loads(json_file.read())
        state['project path'] = os.getcwd()
        state['subprojects'] = {}
        return state


def collect_subprojects(state, root):
    if os.path.isdir(os.path.join(os.getcwd(), 'dependencies')):
        os.chdir('dependencies')
        anchor = os.getcwd()
        for name in os.listdir(os.getcwd()):
            if os.path.isdir(os.path.join(os.getcwd(), name)):
                os.chdir(os.path.join(root, 'subprojects', name))
                subproject = deserialize()
                collect_subprojects(subproject, root)
                state['subprojects'][name] = subproject
                os.chdir(anchor)

        os.chdir('..')


def reconstruct_dependency_graph(state):
    dependencies = {state['name']: state['subprojects'].keys()}
    for name in state['subprojects'].keys():
        dependencies.update(
            reconstruct_dependency_graph(state['subprojects'][name]))

    return dependencies


def reconstruct_build_queue(state):
    queue = []
    graph = reconstruct_dependency_graph(state)
    while len(graph) > 0:
        found_next = False
        next_node = None
        for node, edges in graph.items():
            queueable = True
            for edge in edges:
                if edge not in queue:
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


def relative_path(state):
    return os.path.relpath(os.getcwd(), state['project path'])


def evaluate_test_directory(state):
    path, name = os.path.split(os.path.dirname(relative_path(state)))
    if path:
        while True:
            path, directory = os.path.split(path)
            if not directory == "src":
                name = '.'.join([directory, name])
            else:
                break

    else:
        name = state['name']

    directory = relative_path(state)
    files = []
    for extension in state['file extension']['implementation files']:
        filenames = glob.glob('*.' + extension)
        file_paths = [os.path.join(directory, filename)
                      for filename in filenames]
        files.extend(file_paths)

    state['tests'][name] = files
    os.chdir('..')


def evaluate_directory(state):
    path = relative_path(state)
    for group in state['file extension']:
        files = []
        for extension in state['file extension'][group]:
            filenames = glob.glob('*.' + extension)
            file_paths = [os.path.join(path, filename)
                          for filename in filenames]
            files.extend(file_paths)

        state[group].extend(files)

    for name in os.listdir(os.getcwd()):
        if os.path.isdir(os.path.join(os.getcwd(), name)):
            os.chdir(name)
            if not re.match(state['ignore pattern'], name):
                if re.match(state['test pattern'], name):
                    evaluate_test_directory(state)
                else:
                    evaluate_directory(state)
            else:
                os.chdir('..')

    os.chdir('..')


def collect_driver(state):
    if 'driver' not in state:
        for extension in state['file extension']['implementation files']:
            possible_driver = 'src/main.' + extension
            if possible_driver in state['implementation files']:
                state['driver'] = possible_driver
                state['implementation files'].remove(possible_driver)
                return
        return


def generate(name, language, **kwargs):
    state = {'strict':True, 'test pattern':'test$', 'ignore pattern':'$^',
             'is external project':False}
    state.update({'name':name, 'language':language})
    state.update(kwargs)
    if 'initialized' not in state:
        state.update(configuration.languages[state['language']])
        if state['is external project']:
            state.pop('compiler')
        state['initialized'] = True

    state['project path'] = os.getcwd()
    state['tests'] = {}
    state['subprojects'] = {}
    if 'update' not in state:
        state['update'] = not state['is external project']

    for group in state['file extension']:
        state[group] = []

    os.chdir('src')
    evaluate_directory(state)
    collect_driver(state)
    serialize(state)
