from textwrap import dedent
from jinja2 import Environment
import numpy as np
import json
from pathlib import Path
import os
import sys
import time
import subprocess
## Cannot import things from this module as gives cicular import
import PythoniseMathematica as PythoniseMathematica



def generateModules(
    args, 
    allSymbols, 
    scalarMassMatrixFilePath,
    scalarMassNames,
    scalarPermutationMatrixFilePath,
    scalarRotationMatrixFilePath,
    vectorMasses,
    vectorShorthands,
    gccFlags,
    fieldNames
):
    
    (CythonModulesDir := Path("../Build/CythonModules")).mkdir(
        exist_ok=True, 
        parents=True
    )   
    
    if args.verbose:
        print("Generating cython modules")
    
    loopOrder = args.loopOrder 
    
    veffFilePaths   = [args.loFilePath, args.nloFilePath]
    veffNames = ["lo", "nlo"]
    if loopOrder >1:
        veffFilePaths.append(args.nnloFilePath)
        veffNames.append("nnlo")
        
    veffSubModules = []
    for idx, name in enumerate(veffNames):
        veffSubModules.append(generateVeffSubModule(
            name, 
            veffFilePaths[idx], 
            allSymbols
        ))
    
    computeMassesModule = generateComputeMassesModule(
        allSymbols,
        scalarMassMatrixFilePath,
        scalarMassNames,
        scalarPermutationMatrixFilePath,
        scalarRotationMatrixFilePath,
        vectorMasses,
        vectorShorthands,
        loopOrder,
    )
    
    generateEvaluatePotentialModule(
        CythonModulesDir / "evaluatePotential.pyx", 
        loopOrder,
        allSymbols, 
        fieldNames,
        veffSubModules,
        computeMassesModule,
    )
    
    generateSetupFile(
        CythonModulesDir / "Setup.py",
        loopOrder, 
        gccFlags,
        args.profile
    )
    compileCythonModules(args.verbose, CythonModulesDir)
    
def generateSetupFile(
    fileName, 
    loopOrder, 
    gccFlags,
    profile
):  
    with open(fileName, 'w') as file:
        file.writelines(Environment().from_string(dedent("""\
            #!/usr/bin/env python3
            # -*- coding: utf-8 -*-
            from setuptools import setup, Extension
            from Cython.Build import cythonize
            extensions = [Extension("evaluatePotential", ["evaluatePotential.pyx"], extra_compile_args = {{gccFlags}})]

            setup(
                name="Veff_cython",
                ext_modules = cythonize(
                        extensions, 
                        compiler_directives={
                            "language_level": "3", 
                            "boundscheck": False,
                            "wraparound": False,
                            "profile": {{profile}},
                            }
                ),
            )
            """
        )).render(loopOrder = loopOrder, 
        gccFlags = [f"-{flag}" for flag in gccFlags],
        profile = profile))
    
def generateEvaluatePotentialModule(
    filename, 
    loopOrder, 
    allSymbols, 
    fieldNames, 
    veffSubModules, 
    computeMassesModule
):
    with open(filename, 'w') as file:
        file.write(Environment().from_string(dedent(
        """
        from libc.complex cimport csqrt
        from libc.complex cimport clog

        cpdef evaluatePotential(fields, double [:] parameters):
        
        {% for name in fieldNames %}
            parameters[{{ allSymbols.index(name) }}] = fields[{{ loop.index0 }}]
         {%- endfor %}
        
            computeMasses(parameters)
            
        {%- for symbol in allSymbols %}
            cdef double {{ symbol }} = parameters[{{ loop.index0 }}]
        {%- endfor %}
            valueLO = _lo(
        {%- for symbol in allSymbols %}
            {{ symbol }},
        {%- endfor %}
            )
            valueNLO = _nlo(
        {%- for symbol in allSymbols %}
            {{ symbol }},
        {%- endfor %}
            )
        {%- if loopOrder > 1 %}
            valueNNLO = _nnlo(
        {%- for symbol in allSymbols %}
            {{ symbol }},
        {%- endfor %}
                )
            return valueLO + valueNLO + valueNNLO
            
        {%- else %}
            return valueLO + valueNLO
        {%- endif %}
        
        {%- for veffSubModule in veffSubModules %}
        {{ veffSubModule }}
        {%- endfor %}
        
        {{computeMassesModule}}
        
        """)).render(
        loopOrder=loopOrder, 
        allSymbols=allSymbols, 
        fieldNames=fieldNames, 
        veffSubModules = veffSubModules, 
        computeMassesModule = computeMassesModule
        )
        )
def generateVeffSubModule(name, veffFp, allSymbols):
    # Creates a cython module with that computes an order of Veff
    ## NOTE this is the one thing the can return complex
    return Environment().from_string(dedent("""\
            cdef double complex _{{ name }}(
            {%- for symbol in allSymbols %}
                float {{ symbol }},
            {%- endfor %}
                ):
                cdef double complex a = 0.0
            {%- for op, term in opsAndExpressions %}
                a {{ op }} {{ term }}
            {%- endfor %}
                return a
            """)).render(name=name, allSymbols=allSymbols, opsAndExpressions=np.transpose(mutliLineExpression(veffFp)))


def generateComputeMassesModule(
    allSymbols, 
    scalarMassMatrixFile, 
    scalarMassNames,
    scalarPermutationMatrixFile,
    scalarRotationMatrixFile,
    vectorMasses,
    vectorShorthands,
    loopOrder
):
    with open(scalarMassMatrixFile) as file:
        scalarMassMatrices = [convertMatrixToCythonSyntax(line) for line in file.readlines()]
    if "none" in scalarPermutationMatrixFile.lower():
        scalarPermutationMatrix = None
    else:
        with open(scalarPermutationMatrixFile) as file:
            scalarPermutationMatrix = file.readline()

    with open(scalarRotationMatrixFile) as file:
        scalarRotationMatrix = json.loads(file.read())
    
    return Environment().from_string(dedent("""\
        from scipy.linalg import lapack, block_diag
        from numpy import array, sqrt
        from scipy.linalg.blas import dgemm
        
        cdef void computeMasses(double [:] params):
        {%- for symbol in allSymbols %}
            cdef double {{ symbol }} = params[{{ loop.index0 }}]
        {%- endfor %}

        {%- for expression in vectorMasses %}
            params[{{allSymbols.index(expression.identifier)}}] = {{ expression.expression }}
        {%- endfor %}

        {%- for expression in vectorShorthands %}
            params[{{allSymbols.index(expression.identifier)}}] = {{ expression.expression }}
        {%- endfor %}
            #cdef double Tsq = T*T
            #cdef double TsqInverse = 1.0 / Tsq
        {%- for scalarMassMatrix in scalarMassMatrices %}
            scalarMassMatrix{{ loop.index0 }} = array({{ scalarMassMatrix -}}, dtype=float, order='F')# * TsqInverse
            
            eigenValues{{ loop.index0 }}, eigenVectors{{ loop.index0 }}, _ = lapack.dsyevd(scalarMassMatrix{{ loop.index0 }}, compute_v = {{bEigenVectors}})
            #eigenValues{{ loop.index0 }} *= Tsq
        {%- endfor %}
        
        {%- if bEigenVectors %}
            eigenVectors = block_diag(
        {%- for scalarMassMatrix in scalarMassMatrices %}
                eigenVectors{{ loop.index0 }},
        {%- endfor %}
            )
        
        {%- if not scalarPermutationMatrix == none %}
            scalarPermutationMatrix = {{ scalarPermutationMatrix }}
            eigenVectors = dgemm(1,  scalarPermutationMatrix, eigenVectors)
        {%- endif %}
        
        {%- for symbol, indices in scalarRotationMatrix.items() %}
            params[{{allSymbols.index( symbol )}}] = eigenVectors[{{ indices[0] }}][{{ indices[1] }}]
        {%- endfor %}

        {%- endif %}
        {% set scalarMassMatrixLength = (scalarMassNames | length) / (scalarMassMatrices | length) | int %}
        {%- for symbol in scalarMassNames %}
            params[{{allSymbols.index( symbol )}}] = eigenValues{{ (loop.index0 / scalarMassMatrixLength) | int }}[{{ (loop.index0 % scalarMassMatrixLength) | int }}]
        {%- endfor %}

        """)).render(
            allSymbols=allSymbols, 
            scalarMassMatrices = scalarMassMatrices,
            scalarMassNames = scalarMassNames,
            scalarPermutationMatrix = scalarPermutationMatrix,
            scalarRotationMatrix = scalarRotationMatrix,
            vectorMasses = vectorMasses,
            vectorShorthands = vectorShorthands,
            bEigenVectors = 0 if loopOrder ==1 else 1,
        )

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

def compileCythonModules(verbose, cythonModuleDir):
    if not os.path.isfile(cythonModuleDir / "Setup.py"):
        raise FileNotFoundError(f"No Setup.py found in {cythonModuleDir}")
    
    if verbose:
        print("Compiling cython modules")
    
    ti = time.time()
    result = subprocess.run(
        [sys.executable, "Setup.py", "build_ext", "--inplace"],
        cwd=cythonModuleDir,
        capture_output=True,
        text=True,
    )
    tf = time.time()

    if result.returncode != 0:
        print("Compilation failed:")
        print(result.stderr)
        raise RuntimeError("Cython build failed")
    else:
        if verbose:        
            print("Cython compilation succeeded:")
            print(result.stdout)
            print(f'Compilation took {tf - ti} seconds.')
        
    # TODO: Add a clean up step to remove any compilation artifacts.
