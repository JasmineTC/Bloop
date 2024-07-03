import numpy as np

from typing import Callable, Tuple

from .ParsedExpression import ParsedExpression
from .CommonUtils import replaceGreekSymbols

""" Notes.
Suppose you give DRalgo a model with a parameter called lam123. By default the parameter name after integrating out hard modes is then lam1233d, 
ie. DRalgo appends a '3d'. I don't want the '3d' part, but this is easy to remove after reading the file in a dict.  
"""

## TODO replace this with ParsedExpressionSystem. It does almost the same things already.

## This would make a good abstract class
class ParameterMatching:
    ## This list specifies names for parameters that enter the matching relations. 
    # Needs to match symbol names in the file that defines matching relations except that:
    # 1. "Unicode" symbols like λ are automatically converted to ANSI according to rules in CommonUtils.replaceGreekSymbols(). For example λ -> lam  
    # 2. You can include extra symbols that do not appear explicitly in the matching relations. Eg. "RGScale" 
    #parameterNames = [ 'T', 'Lb', 'Lf', 'g1', 'g2', 'g3',  ]  
    ## LN: not used currently since the symbols are read automatically from .txt

    def __init__(self, lines):
        ## Dict for storing ParsedExpression objects
        self.matchingRelations = {}
        self.parameterNames = [] ## Automatically find all symbols that appear in matching relations

        for line in lines:
            lhs, rhs = map(str.strip, line.split("->"))

            from ThreeHiggs.MathematicaParsers import parseExpression
            expr = ParsedExpression(parseExpression(rhs))

            ## Replace Greeks also on the LHS and store in dict as LHS : parsedRHS
            self.matchingRelations[ replaceGreekSymbols(lhs) ] = expr

            ## find symbols but store as string, not the sympy type  
            for symbol in expr.symbols:
                ## NOTE this conversion will cause issues with pathological symbols like "A B" 
                # https://stackoverflow.com/questions/59401738/convert-sympy-symbol-to-string-such-that-it-can-always-be-parsed
                symbol_str = str(symbol)
                if symbol_str not in self.parameterNames:
                    self.parameterNames.append(symbol_str)

        self.parameterNames.sort()

        """Modifies notation in matching relations so that the matched param dict does not use the "3d" suffix (comes from DRalgo by default)
        """
        return

        if not arg:
            return

        newDict = self.matchingRelations

        # DRalgo also gives gauge couplings as "g13d^2" etc, which is terrible => change to g1sq.
        # Crazy oneliner, creates a new dict where just the key names are different:
        newDict = { key.replace("^2", "sq") : value for key, value in newDict.items() }

        """ For ultrasoft theory DRalgo appends "US" => remove that too. Gauge couplings again need special treatment."""
        newDict = { key[:-len("US")] if key.endswith("US") else key : value for key, value in newDict.items() }
        newDict = { key.replace("USsq", "sq") if key.endswith("USsq") else key : value for key, value in newDict.items() }

        ## Remove "3d" suffix with even crazier oneliner (suffix meaning that it's removed only from end of the string)
        newDict = { key[:-len("3d")] if key.endswith("3d") else key : value for key, value in newDict.items() }

        ## Gauge couplings are originally of form g3d^2 so account for that too 
        newDict = { key.replace("3dsq", "sq") if key.endswith("3dsq") else key : value for key, value in newDict.items() }
    
        self.matchingRelations = newDict

    def __call__(self, inputParams):
        return {key: expr(inputParams) for key, expr in self.matchingRelations.items()}

