(* ::Package:: *)

(************************************************************************)
(*
This is a collection of functions that can be used to manipulate mathematica
functions into a form we can turn into excutable pythong.
Some functions have been written by Claude v3.5 and are denoted as such
*)
(************************************************************************)



exportUTF8[fileName_, expr_] := Export[fileName, expr, CharacterEncoding -> "UTF-8"\[NonBreakingSpace]];


(* ::Input::Initialization:: *)
(** Replace non-numeric matrix elements by simple symbols. eg. long expression at matrix[[i,j]] -> elementSymbolX
where X is a number starting from 0. Each element gets a unique symbol. Returns the symbolic matrix and list of shorthands. **)
toSymbolicMatrix[matrix_, elementSymbol_, bIsSymmetric_: False] := Block[
{i, j, count, symbolicMatrix, shape, rows, columns, shorthandDefinitions, tempMatrix},
	
	count = 0;
	shape = Dimensions[matrix];
	{rows, columns} = shape;
	
	(** don't replace numerical elements (esp. 0 or 1) **) 
	IsTrivialQ[el_] := Return[ NumericQ[el] ];
	
	(** Helper. This is messy but am lazy. SetAttributes is used in order to modify the shorthandList argument **)
	SetAttributes[AppendSymbolicShorthand, HoldAll];
	AppendSymbolicShorthand[el_, shorthandBase_, shorthandList_] := Block[{substRule, shorthand},
		If[ !IsTrivialQ[el],
			shorthand = ToExpression[ToString[shorthandBase]<>ToString[count]];
			substRule = shorthand -> el;
			AppendTo[shorthandList, substRule];
			count++;
			Return[shorthand];
			,
			(* else *)
			Return[el];
		];
	];
	
	shorthandDefinitions = {};
		
	If[!bIsSymmetric,
		symbolicMatrix = Table[ AppendSymbolicShorthand[ matrix[[i,j]], elementSymbol, shorthandDefinitions ], {i,1,rows},{j,1,columns}];
	,
	(* Symmetric matrix: fill in symbols in the upper right half only *)
		tempMatrix = Table[
			AppendSymbolicShorthand[ matrix[[i,j]], elementSymbol, shorthandDefinitions ],
		{i,1,rows},{j,1,i}];
		
		symbolicMatrix = Table[
			If[j<=i, tempMatrix[[i,j]], tempMatrix[[j,i]]], {i,1,rows}, {j,1,columns}
		];
	];
	
	Return[{symbolicMatrix, shorthandDefinitions}];
];


(* Written by Claude v4*)
extractSymbolsDetailed3[expr_, includeProtected_: False, asStrings_: True] := Module[{symbols, lhsSymbols, rhsSymbols, filterProtected, stringConverter},
  
  (* Helper function to filter protected symbols *)
  filterProtected = If[includeProtected, Identity, Select[#, !MemberQ[Attributes[#], Protected] &] &];
  
  (* Helper function to convert to strings if requested *)
  stringConverter = If[asStrings, 
    Select[#, Head[#] === Symbol &] /. s_Symbol :> SymbolName[s] &, 
    Identity
  ];
  
  Which[
    (* Handle list of substitution rules *)
    ListQ[expr] && AllTrue[expr, MatchQ[#, _Rule | _RuleDelayed] &],
    lhsSymbols = stringConverter@DeleteDuplicates@filterProtected@Cases[expr[[All, 1]], _Symbol, Infinity];
    rhsSymbols = stringConverter@DeleteDuplicates@filterProtected@Cases[expr[[All, 2]], _Symbol, Infinity];
    <|"LHS" -> lhsSymbols, "RHS" -> rhsSymbols|>,
    
    (* Handle SparseArray *)
    Head[expr] === SparseArray,
    stringConverter@DeleteDuplicates@filterProtected@Cases[expr["NonzeroValues"], _Symbol, Infinity],
    
    (* Handle regular expressions *)
    True,
    stringConverter@DeleteDuplicates@filterProtected@Cases[expr, _Symbol, Infinity]
  ]
]


(*Claude 3.5*)
(* Helper function to remove any suffix from a symbol *)
RemoveSymbolSuffix[expr_Symbol, suffix_String] := 
  If[StringEndsQ[ToString[expr], suffix],
     Symbol[StringReplace[ToString[expr], suffix -> ""]],
     expr]

(* Main function to recursively process expressions *)
RemoveExpressionSuffix[expr_, suffix_String] := 
  expr /. s_Symbol :> RemoveSymbolSuffix[s, suffix]

(* Process entire rule list for a single suffix *)
RemoveSingleSuffix[rules_List, suffix_String] := 
  Map[Rule[
    RemoveSymbolSuffix[#[[1]], suffix], 
    RemoveExpressionSuffix[#[[2]], suffix]
  ] &, rules]

(* Process rules with multiple suffixes *)
RemoveSuffixes[rules_List, suffixes_List] := 
  Fold[RemoveSingleSuffix[#1, #2] &, rules, suffixes]


(**  Here the lists are assumed to have elements of form symbol -> expr. 
This function joins substitution rules from list1 and list2 so that the final subst rule is symbol -> expr1 + expr2. **)
CombineSubstRules[list1_, list2_] := Block[{combinedList,groupedRules},

	(** Magic code written by ChatGPT. But it works **)
	
	combinedList = Join[list1, list2];
	(* Group the rules by their left-hand sides *)
	groupedRules = GroupBy[combinedList, #[[1]] &];

	(* Sum up the right-hand sides for each group *)
	resultList = Rule @@@ KeyValueMap[{#1, Total[#2[[All, 2]]]} &, groupedRules];
	Return[resultList];
];


(** For table building etc it's convenient to use indices with []. 
But for exporting let's get rid of the [] since those are actually function calls in Mathematica.
Instead, we will just use a symbol with numbers attached to denote the indices, but to guarantee each
combination of indices produces an unique symbol we have to pad the indices with zeros.
**)

(* Turns a number into string and pads it with leading zeros *)
toPaddedString[idx_, numZeros_Integer] := Block[{},
	If[numZeros < 1, ToString[idx], 
		StringJoin@ConstantArray["0", numZeros]<>ToString[idx]
	]
];

toIndexedSymbol[symbol_, idx_, minDigits_Integer: 1] := Block[{paddedIdx},
	paddedIdx = toPaddedString[idx, minDigits-1 - Floor[Log10[idx]]];
	
	ToExpression[ ToString[symbol]<>paddedIdx ]
];

toIndexedSymbol2[symbol_, idx1_, idx2_, minDigits_Integer: 1] := Block[{paddedIdx1, paddedIdx2},
	paddedIdx1 = toPaddedString[idx1, minDigits-1 - Floor[Log10[idx1]]];
	paddedIdx2 = toPaddedString[idx2, minDigits-1 - Floor[Log10[idx2]]];
	
	ToExpression[ ToString[symbol]<>paddedIdx1<>paddedIdx2 ]
];


(*(*Claude 3.5 *)
(*Function to extract symbols from both sides of substitutions*)
GetSubstitutionSymbols[substitutions_]:=Module[{leftSymbols={},rightSymbols={}},(*Handle single substitution or list of substitutions*)substList=If[Head[substitutions]===List,substitutions,{substitutions}];
(*Helper function to get symbols from an expression*)getSymbols[expr_]:=Union[(*Handle the case where expr is itself a symbol*)If[Head[expr]===Symbol,{expr},{}],(*Handle nested symbols*)Cases[expr,s_Symbol:>s,Infinity]];
(*Process each substitution rule*)Do[(*Extract both sides using pattern matching for Rule*)lhs=sub[[1]];(*First part of the Rule*)rhs=sub[[2]];(*Second part of the Rule*)(*Get all symbols from both sides*)leftSymbols=Union[leftSymbols,getSymbols[lhs]];
rightSymbols=Union[rightSymbols,getSymbols[rhs]];,{sub,substList}];
(*Return both lists*){leftSymbols,rightSymbols}]*)


(*extractSymbols[expr_] := Union[Cases[expr, s_Symbol /; 
    (* Exclude built-in symbols *)
    !MemberQ[Attributes[s], Protected] && 
    (* Exclude temporary pattern variables *)
    !StringMatchQ[SymbolName[s], "$" ~~ ___] &&
    (* Exclude context-specific symbols *)
    Context[s] =!= "System`", 
    Infinity]]*)


(*(* This is code generated by claude ai, V3.5, I do not pretend to understand it, but it seems to work (with some help with recursion).
Please double check it gives sensible answers!!!*) 
(* Function to remove numerical coefficients and split a single expression into its components *)
splitExpr[expr_] := Module[{terms},
  (* Convert expression to a list of terms split by Plus *)
  terms = List @@ If[Head[expr] === Plus, expr, {expr}];
  
  (* Handle negative terms and extract symbols *)
  DeleteDuplicates @ Flatten[
    terms //. {
      (* Handle division by sqrt and complex combinations *)
      Times[x_Symbol, rest___] :> x,
      (* Handle fractions *)
      x_ /; Head[x] === Times :> List @@ Cases[x, _Symbol],
      (* Handle division *)
      x_ /; Head[x] === Power :> Which[
        MatchQ[x, Power[_, Rational[1, 2]]], First[x],  (* Handle sqrt *)
        MatchQ[x, Power[_, Rational[-1, 2]]], First[x], (* Handle 1/sqrt *)
        True, First[x]
      ],
      (* Ignore pure numbers, I, and sqrt of numbers *)
      _?NumberQ | I | Sqrt[_?NumberQ] -> Sequence[]
    }
  ]
];
symbolsToStrings[symbols_List] := ToString /@ symbols*)
