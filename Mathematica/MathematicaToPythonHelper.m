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


(* Written by Claude v4, bug fixed by ChatGPT o4-mini-high*)
extractSymbols[expr_, includeProtected_: False, asStrings_: True] :=
 Module[{syms, lhs, rhs, keepQ, toStr},

  (* keep or drop protected symbols *)
  keepQ = If[includeProtected, True &, ! MemberQ[Attributes[#], Protected] &];

  (* convert ONLY non-System symbols to strings *)
  toStr[s_] := If[asStrings, SymbolName[s], s];

  (* core extractor *)
  symsFrom[e_] :=
    DeleteDuplicates @ Select[Cases[e, _Symbol, Infinity], keepQ];

  Which[
   (* list of rules *)
   ListQ[expr] && AllTrue[expr, MatchQ[#, _Rule | _RuleDelayed] &],
   lhs = toStr /@ symsFrom[expr[[All, 1]]];
   rhs = toStr /@ symsFrom[expr[[All, 2]]];
   <|"LHS" -> lhs, "RHS" -> rhs|>,

   (* SparseArray *)
   Head[expr] === SparseArray,
   toStr /@ symsFrom[expr["NonzeroValues"]],

   True,
   toStr /@ symsFrom[expr]
  ]
]


symbolsFromDict[rules_List] :=
  DeleteDuplicates @ Flatten @ Cases[
    rules,
    r : (_Rule | _RuleDelayed) /; FreeQ[r[[2]], _Rule | _RuleDelayed] :> r[[2]],
    {0, \[Infinity]}
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


(*Gauge couplings given as g_i^2 even though g_i is what is needed, so take sqrt. Fine to do so long as g_i >0*)
sqrtSubRules[ruleList_]:=Module[{newRules},
  newRules = ruleList /. (lhs_ -> expr_) /; MatchQ[lhs, _^2] :> (PowerExpand[Sqrt[lhs]] -> Sqrt[expr]);
  Return[newRules];
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
