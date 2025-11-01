import os
from textwrap import dedent
from jinja2 import Environment
import numpy as np
import json
## Cannot import things from this module as gives cicular import
import Bloop.PythoniseMathematica as PythoniseMathematica


def generateModules(
    args, 
    allSymbols, 
    scalarMassMatrixFile,
    scalarMassNames,
    scalarPermutationMatrixFile,
    scalarRotationMatrixFile,
    vectorMasses,
    vectorShorthands,
    gccFlags
):
    
    parent_dir = os.path.dirname(os.getcwd())
    data_dir   = os.path.join(parent_dir, 'src', 'Bloop')
    module_dir = os.path.join(parent_dir, 'src', 'Bloop', 'CythonModules')
    
    if not os.path.exists(module_dir):
        os.mkdir(module_dir)
    if args.verbose:
        print("Generating Veff submodule")
    
    loopOrder = args.loopOrder 
    
    veffFilePaths   = [args.loFilePath, args.nloFilePath]
    veffNames = ["lo", "nlo"]
    if loopOrder >1:
        veffFilePaths.append(args.nnloFilePath)
        veffNames.append("nnlo")
        
    for idx, name in enumerate(veffNames):
        generateVeffSubModule(
            name, 
            os.path.join(module_dir, f"{name}.pyx"), 
            os.path.join(data_dir, veffFilePaths[idx]), 
            allSymbols
        )

    generateDiagonalizeSubModule(
        os.path.join(module_dir, "eigen.pyx"), 
        allSymbols,
        os.path.join(data_dir, scalarMassMatrixFile),
        scalarMassNames,
        os.path.join(data_dir, scalarPermutationMatrixFile),
        os.path.join(data_dir, scalarRotationMatrixFile),
        vectorMasses,
        vectorShorthands,
    )

    generateVeffModule(
        os.path.join(module_dir, 'VeffTotal.py'), 
        loopOrder, 
        allSymbols
    )
    
    #=============================== setup file ==============================#
    with open(os.path.join(module_dir, 'Setup.py'), 'w') as file:
        file.writelines(Environment().from_string(dedent("""\
            #!/usr/bin/env python3
            # -*- coding: utf-8 -*-
            from setuptools import setup, Extension
            from Cython.Build import cythonize
            
            extensions = [Extension("lo", ["lo.pyx"], extra_compile_args = {{gccFlags}})]
            extensions.append(Extension("nlo", ["nlo.pyx"], extra_compile_args = {{gccFlags}}))
            {% if args.loopOrder > 1 %}
            extensions.append(Extension("nnlo", ["nnlo.pyx"], extra_compile_args = {{gccFlags}}))
            {% endif %}
            extensions.append(Extension("eigen", ["eigen.pyx"], extra_compile_args = {{gccFlags}}))

            setup(
                name="Veff_cython",
                ext_modules=cythonize(
                    extensions, compiler_directives={"language_level": "3"}
                ),
            )
            """
        )).render(args = args, gccFlags = [f"-{flag}" for flag in gccFlags] ))
        
def generateVeffModule(filename, loopOrder, allSymbols):
    """Write a  function that imports veff submodules based on loopOrder,
    returns the evaluated submodules as a tuple.
    """
    with open(filename, 'w') as file:
        file.write(Environment().from_string(dedent(
        """\
        from Bloop.CythonModules.lo import lo 
        from Bloop.CythonModules.nlo import nlo 
        {%- if loopOrder > 1 %}
        from Bloop.CythonModules.nnlo import nnlo 
        {%- endif %}

        def veffTotal(
        {%- for symbol in allSymbols %}
            {{ symbol }} = 1,
        {%- endfor %}
            ):
            val_lo = lo(
        {%- for symbol in allSymbols %}
                {{ symbol }},
        {%- endfor %}
            )
            
            val_nlo = nlo(
        {%- for symbol in allSymbols %}
                {{ symbol }},
        {%- endfor %}
            )
            
        {%- if loopOrder > 1 %}
            val_nnlo = nnlo(
        {%- for symbol in allSymbols %}
                {{ symbol }},
        {%- endfor %}
            )
            return (val_lo, val_nlo, val_nnlo)
        
        {%- else %}
            return (val_lo, val_nlo)
        {%- endif %}
        """)).render(loopOrder=loopOrder, allSymbols=allSymbols))
     

def generateVeffSubModule(name, moduleName, veffFp, allSymbols):
    # Creates a cython module with that computes an order of Veff
    with open(moduleName, 'w') as file:
    
        file.write(Environment().from_string(dedent("""\
            #cython: cdivision=False
            from libc.complex cimport csqrt
            from libc.complex cimport clog
            
            cpdef double complex {{ name }}(
            {%- for symbol in allSymbols %}
                double complex {{ symbol }},
            {%- endfor %}
                ):
                ## Calling _name decreases compile time, maybe increases perfomance
                return _{{ name }}(
            {%- for symbol in allSymbols %}
                    {{ symbol }},
            {%- endfor %}
                )
            
            cdef double complex _{{ name }}(
            {%- for symbol in allSymbols %}
                double complex {{ symbol }},
            {%- endfor %}
                ):
                cdef double complex a = 0.0
            {%- for op, term in opsAndExpressions %}
                a {{ op }} {{ term }}
            {%- endfor %}
                return a
            """)).render(name=name, allSymbols=allSymbols, opsAndExpressions=np.transpose(mutliLineExpression(veffFp))))


def generateDiagonalizeSubModule(
    moduleName, 
    allSymbols, 
    scalarMassMatrixFile, 
    scalarMassNames,
    scalarPermutationMatrixFile,
    scalarRotationMatrixFile,
    vectorMasses,
    vectorShorthands,
):
    with open(scalarMassMatrixFile) as file:
        scalarMassMatrices = [convertMatrixToCythonSyntax(line) for line in file.readlines()]
    if scalarPermutationMatrixFile.lower() == "none":
        scalarPermutationMatrixFile = scalarPermutationMatrixFile.lower()
    else:
        with open(scalarPermutationMatrixFile) as file:
            scalarPermutationMatrix = convertMatrixToCythonSyntax(file.read())

    with open(scalarRotationMatrixFile) as file:
        scalarRotationMatrix = json.loads(file.read())

    # Creates a cython module with that computes an order of Veff
    with open(moduleName, 'w') as file:
        file.write(Environment().from_string(dedent("""\
            from scipy.linalg import lapack
            from scipy.linalg import block_diag
            from scipy.linalg.blas import dgemm
            from numpy import divide 
            from numpy import sqrt

            cpdef void eigen(complex [:] parameters):
            {%- for symbol in allSymbols %}
                cdef double complex * {{ symbol }} = &parameters[{{ loop.index0 }}]
            {%- endfor %}

                _eigen(
            {%- for symbol in allSymbols %}
                    {{ symbol }},
            {%- endfor %}
                )

            cdef void _eigen(
            {%- for symbol in allSymbols %}
                double complex * _{{ symbol }},
            {%- endfor %}
            ):
            {%- for symbol in allSymbols %}
                cdef double {{ symbol }} = _{{ symbol }}[0].real
            {%- endfor %}

            {%- for expression in vectorMasses %}
                {{ expression.identifier }} = {{ expression.expression }}
                _{{ expression.identifier }}[0] = {{ expression.identifier }}
            {%- endfor %}

            {%- for expression in vectorShorthands %}
                {{ expression.identifier }} = {{ expression.expression }}
                _{{ expression.identifier }}[0] = {{ expression.identifier }}
            {%- endfor %}

            {%- for scalarMassMatrix in scalarMassMatrices %}
                scalarMassMatrix{{ loop.index0 }} = divide({{ scalarMassMatrix -}}, (T ** 2))
                eigenValues{{ loop.index0 }}, eigenVectors{{ loop.index0 }}, _ = lapack.dsyevd(scalarMassMatrix{{ loop.index0 }}, compute_v = 1)
                eigenValues{{ loop.index0 }} *= (T ** 2)
            {%- endfor %}
            
                scalarPermutationMatrix = {{ scalarPermutationMatrix }}
                eigenVectors = block_diag(
            {%- for scalarMassMatrix in scalarMassMatrices %}
                    eigenVectors{{ loop.index0 }},
            {%- endfor %}
                )
                
                {%- if not scalarPermutationMatrixFile == none %}
                eigenVectors = dgemm(1, scalarPermutationMatrix, eigenVectors)
                {%- endif %}
                

            {%- for symbol, indices in scalarRotationMatrix.items() %}
                _{{ symbol }}[0] = eigenVectors[{{ indices[0] }}][{{ indices[1] }}]
            {%- endfor %}

            {% set scalarMassMatrixLength = (scalarMassNames | length) / (scalarMassMatrices | length) | int %}
            {%- for massSymbol in scalarMassNames %}
                _{{ massSymbol }}[0] = eigenValues{{ (loop.index0 / scalarMassMatrixLength) | int }}[{{ (loop.index0 % scalarMassMatrixLength) | int }}]
            {%- endfor %}
            """)).render(
                allSymbols=allSymbols, 
                scalarMassMatrices = scalarMassMatrices,
                scalarMassNames = scalarMassNames,
                scalarPermutationMatrix = scalarPermutationMatrix,
                scalarRotationMatrix = scalarRotationMatrix,
                vectorMasses = vectorMasses,
                vectorShorthands = vectorShorthands,
                scalarPermutationMatrixFile = scalarPermutationMatrixFile
            ))

def mutliLineExpression(filePointer):
    ## Takes an expressions and breaks it down into a mutli line expression
    ## (Cython seems to struggle with the one line NNLO veff)
    
    with open(filePointer, 'r') as file:
        veff = file.read()
    
    operations = ["+="]
    expressions = []
    
    netBrackets = 0
    start = 0
    
    for i, char in enumerate(veff):
        if char == '(':
            netBrackets += 1
        elif char == ')':
            netBrackets -= 1
        if char == ' ' and netBrackets == 0:
            ##+1 to catch space
            line = veff[start:i+1]
            if line in ["+ ", "- "]:
                operations.append("+=" if line == "+ " else "-=")
            else:
                expressions.append(convertToCythonSyntax(line))
            start = i + 1
    
    # Any remaining characters should just be expressions
    if start < len(veff):
        line = veff[start:]
        expressions.append(convertToCythonSyntax(line))
    return operations, expressions

def convertMatrixToCythonSyntax(term):
    term = convertToCythonSyntax(term)
    term = term.replace('{', '[')
    return term.replace('}', ']')
    
def convertToCythonSyntax(term):
    term = term.replace('Sqrt', 'csqrt')
    term = term.replace('Log', 'clog')
    term = term.replace('[', '(')
    term = term.replace(']', ')')
    term = term.replace('^', '**')
    term = PythoniseMathematica.replaceSymbolsConst(term)
    return PythoniseMathematica.replaceGreekSymbols(term)
