#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 20:42:29 2025

@author: john
"""

import os
import numpy as np

from .mathematica_parsing import read_lines, get_terms





def generate_veff_module(args):
    """Creates the Veff module.
    """
    parent_dir = os.path.dirname(os.getcwd())
    data_dir   = os.path.join(parent_dir, 'src', 'ThreeHiggs')
    module_dir = os.path.join(parent_dir, 'src', 'ThreeHiggs', 'Veff')
    if not os.path.exists(module_dir):
        os.mkdir(module_dir)
    print("Generating Veff submodule")
    
    #================================== lo ===================================#
    name = 'lo'
    lo_file  = os.path.join(data_dir, args.loFile)
    filename = os.path.join(module_dir, 'lo.pyx')
    lo_params = generate_lo_submodule(name, filename, lo_file)
    
    #================================== nlo ==================================#
    name     = 'nlo'
    nlo_file = os.path.join(data_dir, args.nloFile)
    filename = os.path.join(module_dir, 'nlo.pyx')
    nlo_params = generate_lo_submodule(name, filename, nlo_file) if args.loopOrder > 0 else None
    
    #================================== nnlo =================================#
    name = 'nnlo'
    nnlo_file = os.path.join(data_dir, args.nnloFile)
    filename  = os.path.join(module_dir, 'nnlo.pyx')
    nnlo_params = generate_lo_submodule(name, filename, nnlo_file) if args.loopOrder > 1 else None
    
    #================================== Veff =================================#
    filename = os.path.join(module_dir, 'veff.py')
    generate_veff_submodule(filename, lo_params, nlo_params, nnlo_params)
    
    #================================ init file ==============================#
    with open(os.path.join(module_dir, '__init__.py'), 'w') as file:
        file.write("from .veff import *")
    
    #=============================== setup file ==============================#
    with open(os.path.join(module_dir, 'setup.py'), 'w') as file:
        from textwrap import dedent
        from jinja2 import Environment
        file.writelines(Environment().from_string(dedent("""\
            #!/usr/bin/env python3
            # -*- coding: utf-8 -*-
            from setuptools import setup, Extension
            from Cython.Build import cythonize
            
            extensions = [Extension("lo", ["lo.pyx"])]
            {% if args.loopOrder >= 1 %}
            extensions.append(Extension("nlo", ["nlo.pyx"]))
            {% endif %}
            {% if args.loopOrder >= 2 %}
            extensions.append(Extension("nnlo", ["nnlo.pyx"]))
            {% endif %}
            
            setup(
                name="Veff_cython",
                ext_modules=cythonize(
                    extensions, compiler_directives={"language_level": "3"}
                ),
            )
            """
        )).render(args = args))
        
def generate_veff_submodule(filename, lo_params, nlo_params, nnlo_params):
    """Creates a submodule with Veff and Veff_params functions (see below).
    """
    if os.path.exists(filename):
        os.remove(filename)
    
    with open(filename, 'w') as file:
        write_veff_function(file, lo_params, nlo_params, nnlo_params)
        file.write('\n\n')
        write_veff_params_function(file, lo_params, nlo_params, nnlo_params)



def write_veff_function(file, lo_params, nlo_params, nnlo_params):
    """Adds function called Veff to the given file. Veff imports lo, nlo and
    nnlo functions and evaluates them, returning the results in a tuple.
    """
    file.write('from .lo import lo\n')

    if nlo_params is not None:
        file.write('from .nlo import nlo\n')

    if nnlo_params is not None:
        file.write('from .nnlo import nnlo\n')

    file.write('\n')
    
    # Function name and input
    file.write('def Veff(\n')
    
    params = np.unique(np.concatenate((
        lo_params, 
        nlo_params if nlo_params is not None else [], 
        nnlo_params if nnlo_params is not None else [], 
    )))

    for param in params:
        param = convert_to_cython_syntax(param)
        file.write(f'    {param} = 1,\n')
    
    file.write('    ):\n')
    
    # Function body
    file.write('    val_lo = lo(\n')
    for param in lo_params:
        param = convert_to_cython_syntax(param)
        file.write(f'        {param},\n')
    file.write('    )\n')
    
    if nlo_params is not None:
        file.write('    val_nlo = nlo(\n')
        for param in nlo_params:
            param = convert_to_cython_syntax(param)
            file.write(f'        {param},\n')
        file.write('    )\n')
    else:
        file.write('    val_nlo = 0\n')

    if nnlo_params is not None:
        file.write('    val_nnlo = nnlo(\n')
        for param in nnlo_params:
            param = convert_to_cython_syntax(param)
            file.write(f'        {param},\n')
        file.write('    )\n')
    else:
        file.write('    val_nnlo = 0\n')
    
    file.write('    return (val_lo, val_nlo, val_nnlo)\n')
    


def write_veff_params_function(file, lo_params, nlo_params, nnlo_params):
    """Adds a function called Veff_params to the given file. Veff_params
    returns a tuple of parameters for the corresponding Veff function. The
    parameters are pulled from a provided parameter dictionary.
    """
    file.write('def Veff_params(params):\n')
    file.write('    return (\n')
    
    params = np.unique(np.concatenate((
        lo_params, 
        nlo_params if nlo_params is not None else [], 
        nnlo_params if nnlo_params is not None else [], 
    )))

    for param in params:
        param = convert_to_cython_syntax(param)
        file.write(f'        params["{param}"],\n')
    file.write('    )\n')

def generate_lo_submodule(name, filename, lo_file):
    """Creates a cython module with a function that evaluates an expression for
    Veff. 
    
    The expression is assumed to be broken into a list of `terms` to be added 
    together. The sign of each term should be in `signs` and `params` is an 
    array of parameters that appear in the expression.
    """
    if os.path.exists(filename):
        os.remove(filename)
        
    lines = read_lines(lo_file)
    params, signs, terms = get_terms(lines)
    
    with open(filename, 'w') as file:
        # Function imports used by Veff
        file.write('# cython: cdivision=True\n')
        file.write('from libc.complex cimport csqrt\n')
        file.write('from libc.complex cimport clog\n')
        file.write('from libc.math cimport M_PI\n')
        file.write('\n')
        
        # Function declaration and inputs
        file.write(f'cpdef double complex {name}(\n')
        
        for param in params:
            param = convert_to_cython_syntax(param)
            file.write(f'    double complex {param},\n')
        
        file.write('    ):\n')
        file.write(f'    return _{name}(\n')
        for param in params:
            param = convert_to_cython_syntax(param)
            file.write(f'        {param},\n')
        file.write('    )\n\n\n')
        
        file.write(f'cdef double complex _{name}(\n')
        
        for param in params:
            param = convert_to_cython_syntax(param)
            file.write(f'    double complex {param},\n')
        
        file.write('    ):\n')
        
        # Function body
        file.write('    cdef double complex a = 0.0\n')
        
        term = convert_to_cython_syntax(terms[0])
        file.write(f'    a += {term}\n')
        
        for sign, term in zip(signs, terms[1:]):
            term = convert_to_cython_syntax(term)
            if sign > 0:
                file.write(f'    a += {term}\n')
            else:
                file.write(f'    a -= {term}\n')
        
        file.write('    return a\n')
        
    return params
        


def convert_to_cython_syntax(term):
    term = term.replace('Sqrt', 'csqrt')
    term = term.replace('Pi', 'M_PI')
    term = term.replace('Log', 'clog')
    term = term.replace('[', '(')
    term = term.replace(']', ')')
    term = term.replace('^', '**')
    term = term.replace('λ', 'lam')
    term = term.replace('μ', 'mu')
    return term
