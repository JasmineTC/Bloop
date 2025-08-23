#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 20:40:32 2025

@author: john
"""

import re
import numpy as np


def read_lines(filename):
    """Reads the expression in the given file and breaks it into a list of 
    lines, each line representing a term in the expression to be summed.
    """
    num_unpaired_brackets = 0
    next_line = []
    lines = []
    
    with open(filename, 'r') as file:
        while True:    
            char = file.read(1)
            if not char:
                lines.append(''.join(next_line))
                break
            
            next_line.append(char)
            
            if char == '(':
                num_unpaired_brackets += 1
            elif char == ')':
                num_unpaired_brackets -= 1
                
            if char == ' ' and num_unpaired_brackets == 0:
                lines.append(''.join(next_line))
                next_line = []
    
    return lines



def find_parameters(expression):
    """Returns a list of parameter names in the given expression.
    """
    pattern = r"[a-zA-Z_\u0370-\u03FF][\w\u0370-\u03FF]*"
    tokens = re.findall(pattern, expression)
    # Exclude known functions or keywords (extend as needed)
    known_functions = {"Sqrt", "Log", "Pi"}  # Add others as necessary
    return [token for token in tokens if token not in known_functions]



def get_terms(lines):
    """Breaks the given list of lines (ie terms to be summed in an expression)
    into two lists: a list containing the terms to be summed and a list of 
    leading signs for each term (excluding the first term).
    
    A set of all parameters that appear in the expression is also returned.
    """
    params = set()
    terms = []
    signs = []
    
    for line in lines:
        if line in ["+ ", "- "]:
            signs.append(1 if line == "+ " else -1)
            continue
        
        for param in find_parameters(line):
            params.add(param)
        
        terms.append(line)
        
    params = np.sort(list(params))
    return params, signs, terms
