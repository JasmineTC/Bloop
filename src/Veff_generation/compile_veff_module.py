#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 17:33:46 2025

@author: john
"""

import os
import sys
import subprocess


def compile_veff_submodule():
    parent_dir = os.path.dirname(os.getcwd())
    module_dir = os.path.join(parent_dir, 'src', 'ThreeHiggs', 'Veff')
    setup_path = os.path.join(module_dir, "setup.py")
    
    if not os.path.isfile(setup_path):
        raise FileNotFoundError(f"No setup.py found in {module_dir}")
    
    print("Compiling Veff submodule")
    result = subprocess.run(
        [sys.executable, "setup.py", "build_ext", "--inplace"],
        cwd=module_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print("Compilation failed:")
        print(result.stderr)
        raise RuntimeError("Cython build failed")
    else:
        print("Cython compilation succeeded:")
        print(result.stdout)
        
    # TODO: Add a clean up step to remove any compilation artifacts.