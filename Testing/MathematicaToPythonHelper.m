(* ::Package:: *)

(************************************************************************)
(*
This is a collection of functions that can be used to manipulate mathematica
functions into a form we can turn into excutable pythong.
Some functions have been written by Claude v3.5 and are denoted as such
*)
(************************************************************************)



ExportUTF8[fileName_, expr_] := Export[fileName, expr, CharacterEncoding -> "UTF-8"\[NonBreakingSpace]];


(* ::Input::Initialization:: *)
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

