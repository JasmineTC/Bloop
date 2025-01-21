(* ::Package:: *)

(** This file defines a 3HDM and computes Veff with 3 background fields. We use permute the fields to bring the mass matrix
into a block diagonal form consisting of two 6x6 matrices. **)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
(* This is pointing to my DRalgo repo for easier updating *)
(*pathToDRalgo = "/home/lani/repos/DRalgo/DRalgo.m"*)
pathToDRalgo = "/home/jasmine/.Mathematica/Applications/DRalgo/DRalgo.m"
Get[pathToDRalgo]


(* ::Section::Closed:: *)
(*Helper functions*)


(** Export expression or list of expressions to file using UTF-8 (important for compatibility) **)
ExportUTF8[fileName_, expr_] := Export[fileName, expr, CharacterEncoding -> "UTF-8"\[NonBreakingSpace]];


(* ::Chapter:: *)
(*3HDM a la Venus*)


(*See 1909.09234 [hep-ph], eq (1) *)


(* ::Subsection:: *)
(*Specify file paths for exporting*)


(** All file paths are relative to the working directory (set above). **)
hardToSoftDirectory = "HardToSoft";
softToUltrasoftDirectory = "SoftToUltrasoft";
effectivePotentialDirectory = "EffectivePotential";
variables = "Variables";


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsDoublet2={{{0,0},{1},1/2},"C"};
HiggsDoublet3={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2,HiggsDoublet3};

ClearAll[g1, g2, g3];
GaugeCouplings={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,GaugeCouplings,RepFermion3Gen,RepScalar];
(** Note that AllocateTensors[] brings some GroupMath symbols to the global namespace. And these are only removed later when calling ImportModelDRalgo...
**)


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


(** Here just list all possible gauge-invariant operators containing 2 doublets **)
(** DRalgo notation is that \[Phi]1\[Phi]2^+ = \!\(
\*SubsuperscriptBox[\(\[Phi]\), \(2\), \(\[Dagger]\)]
\*SubscriptBox[\(\[Phi]\), \(1\)]\ in\ standard\ \(notation . \ \nSo\)\ careful\ here\ to\ make\ sure\ imaginary\ parts\ match\ to\ the\ potential\ in\ our\ draft\) **)
(** I have changed the notation a lot from DRalgo's example 3HDM file. 
My notation for doublet products is \[Phi]ij, where first index is the conjugated doublet **)

InputInv={{1,1},{True,False}}; (*\[Phi]1 \[Phi]1^+*)
\[Phi]11=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,2},{True,False}}; (*\[Phi]2 \[Phi]2^+*)
\[Phi]22=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,3},{True,False}}; (*\[Phi]3 \[Phi]3^+*)
\[Phi]33=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{1,2},{True,False}}; (*\[Phi]1\[Phi]2^+*)
\[Phi]21=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,1},{True,False}};(*\[Phi]2\[Phi]1^+*)
\[Phi]12=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{1,3},{True,False}}; (*\[Phi]1\[Phi]3^+*)
\[Phi]31=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,1},{True,False}};(*\[Phi]3\[Phi]1^+*)
\[Phi]13=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,3},{True,False}}; (*\[Phi]2\[Phi]3^+*)
\[Phi]32=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,2},{True,False}};(*\[Phi]3\[Phi]2^+*)
\[Phi]23=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


(* Define some shorthands for complex parameters. DRalgo performs better if real/imag parts are used separately. *)
\[Mu]12sq = \[Mu]12sqRe + I*\[Mu]12sqIm;
\[Mu]12sqConj = Conjugate[\[Mu]12sq]//ComplexExpand;

\[Lambda]1 = \[Lambda]1Re + I*\[Lambda]1Im;
\[Lambda]1Conj = Conjugate[\[Lambda]1]//ComplexExpand;
\[Lambda]2 = \[Lambda]2Re + I*\[Lambda]2Im;
\[Lambda]2Conj = Conjugate[\[Lambda]2]//ComplexExpand;
\[Lambda]3 = \[Lambda]3Re + I*\[Lambda]3Im;
\[Lambda]3Conj = Conjugate[\[Lambda]3]//ComplexExpand;


(* Quadratic terms. Careful with complex conjugates *)
VMass=-\[Mu]1sq*\[Phi]11 - \[Mu]2sq*\[Phi]22 - \[Mu]3sq*\[Phi]33 - \[Mu]12sq*\[Phi]12 - \[Mu]12sqConj*\[Phi]21 // Simplify // Expand; (* simplify to get rid of imag units *)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


(** Def. quartic terms. Careful with complex conjugates **) 
QuarticTerm1 = \[Lambda]11*\[Phi]11^2 + \[Lambda]22*\[Phi]22^2 + \[Lambda]33*\[Phi]33^2;
QuarticTerm2 = \[Lambda]12*\[Phi]11*\[Phi]22 + \[Lambda]23*\[Phi]22*\[Phi]33 + \[Lambda]31*\[Phi]33*\[Phi]11;
QuarticTerm3 = \[Lambda]12p*\[Phi]12*\[Phi]21 + \[Lambda]23p*\[Phi]23*\[Phi]32 + \[Lambda]31p*\[Phi]31*\[Phi]13;
QuarticTerm4 = \[Lambda]1*\[Phi]12^2 + \[Lambda]2*\[Phi]23^2 + \[Lambda]3*\[Phi]31^2;
(* Hermitian conjugate of QuarticTerm4. Just adding it as ConjugateTranspose[...] didn't work, DRalgo just seemed to get stuck. *)
QuarticTerm5 = \[Lambda]1Conj*\[Phi]21^2 + \[Lambda]2Conj*\[Phi]32^2 + \[Lambda]3Conj*\[Phi]13^2;


VQuartic=QuarticTerm1 + QuarticTerm2 + QuarticTerm3 + QuarticTerm4 + QuarticTerm5 // Simplify // Expand; (* simplify to get rid of imag units *)


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{2,1,2},{False,False,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{3,1,2},{False,False,True}}; 
YukawaDoublet3=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


(*Ysff=-GradYukawa[yt1*YukawaDoublet1+yt2*YukawaDoublet2+yt3*YukawaDoublet3];*)
Ysff=-GradYukawa[yt3*YukawaDoublet3];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt3>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


(* ::Text:: *)
(*Parametric accuracy goal of the EFT matchings need to be specified already in ImportModelDRalgo[] (Mode option). I will first do a LO matching and export that, then import the model again and repeat with order g^4 matching.*)
(*Mode -> 0 : Match couplings at tree level and masses at 1-loop (full g^2)*)
(*Mode -> 1 : Match everything at 1-loop (partial g^4)*)
(*Mode -> 2 : Match couplings at 1-loop and masses at 2-loop (full g^4) *)
(**)
(*However Mode->0 does not really work ATM,  it doesn't give couplings etc...*)


(* ::Subsection::Closed:: *)
(*Helper function for combining LO and NLO substitution rules*)


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


(* ::Subsection::Closed:: *)
(*LO matching. TODO: currently Mode -> 0 is broken...*)


(*
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False, Mode->0];
PerformDRhard[];
*)


(* ::Subsection:: *)
(*NLO matching, by which I mean Mode -> 2*)


(** Normalization4D flag = preserve 4D units so that the EFT path integral weight is e^{-S/T} **)
(** TODO well currently the Normalization4D flag does not work... I've made an issue on Github **)
(** AutoRG->True means that 3D running is built in to the matching. This is bad for automatization since 
the 3D masses become be functions of other 3D parameters. To dodge this we match with AutoRG->False
and do the RG running manually in an additional stage. **)
 
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False, Mode->2, Normalization4D->False, AutoRG->False];
PerformDRhard[];


(* This is code generated by claude ai, V3.5, I do not pretend to understand it, but it seems to work (with some help with recursion).
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
      (* Keep other symbols as is *)
      x_Symbol :> x,
      (* Ignore pure numbers, I, and sqrt of numbers *)
      _?NumberQ | I | Sqrt[_?NumberQ] -> Sequence[]
    }
  ]
];

(* Apply the function to each element in the list and flatten results *)
fourPointSymbols = DeleteDuplicates @ Flatten[splitExpr /@ DeleteDuplicates @ Flatten[splitExpr /@ \[CapitalLambda]4]]
twoPointSymbols = DeleteDuplicates @ Flatten[splitExpr /@ DeleteDuplicates @ Flatten[splitExpr /@ \[Mu]ij]]
yukawaSymbols = DeleteDuplicates @ Flatten[splitExpr /@ DeleteDuplicates @ Flatten[splitExpr /@ Ysff]]


ExportUTF8[variables<>"/scalar4PointSymbols.txt", fourPointSymbols];
ExportUTF8[variables<>"/scalarMassSymbols.txt", twoPointSymbols];
ExportUTF8[variables<>"/YukawaSymbols.txt", yukawaSymbols];
ExportUTF8[variables<>"/GaugeSymbols.txt", GaugeCouplings];


couplingsSoft = PrintCouplings[];
temporalScalarCouplings = PrintTemporalScalarCouplings[];
debyeMasses = PrintDebyeMass["LO"]; (** For Debyes we only take LO result, NLO not needed since we integrate these out anyway **)
scalarMasses = CombineSubstRules[PrintScalarMass["LO"], PrintScalarMass["NLO"]];
softScaleParamsIntermediate = Join[couplingsSoft, temporalScalarCouplings, debyeMasses, scalarMasses];


(* ::Text:: *)
(*DRalgo gives gauge coupling matchings g3d^2 -> ... which is terrible, and furthermore all expressions use g3d instead of its square. So let's change the matching relation so that we get g3d instead of g3d^2. No issue with sqrt signs as long as we calculate with g > 0.*)
(*Also, DRalgo uses \[Lambda]VL[1] etc to label some couplings. But [] means a function call so this is hard to parse as a symbol. So let's change the notation. Easiest way to do this is to just define \[Lambda]VL as a function returning a symbol.*)


gsqSymbols3D = Map[ToExpression[ToString[#]<>"3d"] &, GaugeCouplings];

ReplaceGaugeMatching[ruleList_]:=Module[{newRules},
  newRules = ruleList /. (lhs_ -> expr_) /; MatchQ[lhs, _^2] :> (PowerExpand[Sqrt[lhs]] -> Sqrt[expr]);
  Return[newRules];
];

\[Lambda]VL[i_]:=ToExpression["\[Lambda]VL"<>ToString[i]];
\[Lambda]VLL[i_]:=ToExpression["\[Lambda]VLL"<>ToString[i]];

softScaleParamsIntermediate = ReplaceGaugeMatching[softScaleParamsIntermediate];


softScaleParams = softScaleParamsIntermediate;
ExportUTF8[hardToSoftDirectory<>"/softScaleParams_NLO.txt", softScaleParams]


BetaFunctions4D[];
ExportUTF8[hardToSoftDirectory<>"/BetaFunctions4D[].txt", BetaFunctions4D[]];


(* 3D RG equations can be solved exactly, so do that here. We will export subst rules analogous to the matching relations:
	msq -> msq + \[Beta][msq] Log[\[Mu]3/\[Mu]] where RHS msq is the 3D mass at scale \[Mu] and LHS msq is the mass at scale \[Mu]3 *)
	
SolveRunning3D[betaFunctions_] := Block[{exprList},
	(* Extracting lhs and beta for each list element *)
	(* DRalgo gives the lhs WITHOUT 3D symbol so append it here for consistency *)
	exprList = {ToExpression[ToString[#[[1]]]<>"3d" ], #[[2]]} & /@ betaFunctions;

	(* Make new list with RGE solution on RHS *)
	newRulesList = (#1 -> #1 + #2*Log[goalScale/startScale]) & @@@ exprList;
	Return[newRulesList];
];

running3D = SolveRunning3D[BetaFunctions3DS[]];

ExportUTF8[hardToSoftDirectory<>"/softScaleRGE.txt", running3D]


(* ::Subsection:: *)
(*Soft -> Ultrasoft matching*)


PerformDRsoft[{}];


(** This now works properly as of DRalgo 2023/11/24 update **)
couplingsUS = PrintCouplingsUS[];
scalarMassesUS = CombineSubstRules[PrintScalarMassUS["LO"], PrintScalarMassUS["NLO"]];

ultrasoftScaleParamsIntermediate = Join[couplingsUS, scalarMassesUS] /. \[Mu]3->RGScale;

ultrasoftScaleParams = ReplaceGaugeMatching[ultrasoftScaleParamsIntermediate];


ExportUTF8[softToUltrasoftDirectory<>"/ultrasoftScaleParams_NLO.txt", ultrasoftScaleParams];


runningUS = SolveRunning3D[BetaFunctions3DUS[]];

ExportUTF8[softToUltrasoftDirectory<>"/ultrasoftScaleRGE.txt", runningUS]


(* ::Section:: *)
(*Effective potential*)


(* ::Text:: *)
(*We need to give DRalgo two things: *)
(*1) Rotation matrices for scalar and gauge fields that bring the original field vectors to mass eigenstate basis*)
(*2) Diagonal mass-squared matrices for both scalars and gauges. The ordering of masses needs to match the order to which the rotation matrix brings the fields.*)
(**)
(*We will export a lot of symbolical data for diagonalization, mass eigenvalues, shorthand symbols in rotation matrices etc. These can then be evaluated numerically in an external program. The order of evaluations should be:*)
(*1. Obtain action parameters*)
(*2. Obtain background field values*)
(*3. Solve diagonalization equations for scalars and vectors -> get angles or sines (cosines) of the angles.*)
(*For scalars, use any linear algebra library to diagonalize the mass matrix and get the diagonalizing rotation.*)
(*4. Evaluate masses and rotation matrix elements in the form in which they enter the Veff expression. *)
(*5. Plug all of the above in to the Veff expressions*)
(**)
(*In principle there could be additional shorthands computed between steps 2 and 3 if the diagonalization conditions themselves depend on some shorthand symbols, but we're not including this currently.*)


(* ::Subsection:: *)
(*Specify background fields and init DRalgo stuff*)


UseUltraSoftTheory[];
DefineNewTensorsUS[\[Mu]ij,\[CapitalLambda]4,\[CapitalLambda]3,gvss,gvvv]; 
PrintScalarRepPositions[](** This is supposed to tell which index is which field, but it's cryptic... **)

(** DRalgo ordering: real parts go before imag parts, and this is repeated 3 times (because we have 3 complex doublets). 
So the "usual" place for BG field is second index in each doublet**)

backgroundFieldsFull = {(*\[Phi]1*)0, v1, 0, 0,(*\[Phi]2*)0, v2, 0, 0, (*\[Phi]3*)0, v3, 0, 0}//SparseArray; 
DefineVEVS[backgroundFieldsFull];
(* store a list of all nonzero background field symbols *)
backgroundFields = {v1, v2, v3};


(* ::Subsection:: *)
(*Diagonalizing scalar mass matrix*)


(* ::Text:: *)
(*Finding the diagonalizing matrix analytically is too difficult, so we brute force this by giving DRalgo a symbolic 12x12 matrix with unknown symbols and compute the Veff in terms of those symbols. For masses we give an unknown diagonal matrix. These unknown are field-dependent so we must solve them numerically every time the potential is evaluated. *)
(**)
(*Note that using unknown symbols in the rotation matrix means that DRalgo will compute the effective potential in what is effectively a non-diagonal field basis, ie. there is quadratic mixing. But all the mixing effects vanish when numerical values are fixed by diagonalization conditions later on. (And I believe DRalgo ignores quadratic vertices anyway, but have not confirmed this). *)
(**)
(*Now, there are some optimizations that we can do. For example with the 3-field configuration defined above, the mass matrix can be brought into block-diagonal form by permuting the fields. This reduces the problem to diagonalization of two 6x6 matrices which is faster. We use this approach here; the required 12x12 rotation then has zeros in many places so DRalgo has easier time working with it.*)


scalarMM = PrintTensorsVEV[1]//Normal//Simplify; (* Scalar mass matrix, simplify to get rid of possible imaginary units *)


scalarMM//MatrixForm


(* ::Subsubsection::Closed:: *)
(*Helpers for indexed symbols*)


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


(* ::Subsubsection:: *)
(*Permute scalars to make mass matrix block-diagonal *)


scalarPermutationMatrix = {
{1,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,1,0},
{0,0,1,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,1,0,0,0},
{0,0,0,0,1,0,0,0,0,0,0,0},
{0,0,0,0,0,0,1,0,0,0,0,0},
{0,0,0,0,0,1,0,0,0,0,0,0},
{0,0,0,0,0,0,0,1,0,0,0,0},
{0,0,0,1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,1,0,0},
{0,1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,1}};
(* Permutation matrix swaps the following rows and colums: 2<->11, 4<->9,6<->7  to get a block diagonal matrix
Order is now 
Re\[Phi]1, Im\[Phi]3, Im\[Phi]1, Re\[Phi]3, Re\[Phi]2, Im\[Phi]2 (charged dof)
Re\[Phi]2, Im\[Phi]2, Im\[Phi]1, Re\[Phi]3, Re\[Phi]1, Im\[Phi]3  (neutral dof)*)
If[!OrthogonalMatrixQ[scalarPermutationMatrix], Print["Error, permutation matrix is not orthogonal"]];
(*Our case has Transpose[scalarPermutationMatrix] = scalarPermutationMatrix but taking transpose anyway for consistency/future proofing*)
blockDiagonalMM = Transpose[scalarPermutationMatrix] . scalarMM . scalarPermutationMatrix;
Print["Block diagonal mass matrix:"];
blockDiagonalMM//MatrixForm


(*Extract permutation matrix and do consistency check*)
upperLeftMM = Take[blockDiagonalMM,{1,6},{1,6}];
bottomRightMM = Take[blockDiagonalMM,{7,12},{7,12}];

If[!SymmetricMatrixQ[upperLeftMM] || !SymmetricMatrixQ[bottomRightMM], Print["Error, block not symmetric!"]];

ExportUTF8[effectivePotentialDirectory<>"/scalarPermutationMatrix.txt", scalarPermutationMatrix];


(* ::Subsubsection::Closed:: *)
(*Export scalar mass matrix*)


(* Simplify both blocks by introducing additional symbols, then extract them separately *)
{upperLeftMMSymbolic, upperLeftMMDefinitions} = toSymbolicMatrix[upperLeftMM, "MMUL", True]//Simplify;
{bottomRightMMSymbolic, bottomRightMMDefinitions} = toSymbolicMatrix[bottomRightMM, "MMBR", True]//Simplify;

(* Export expressions separately because we have the code to parse that *)
ExportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_upperLeft.txt", upperLeftMMSymbolic];
ExportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_upperLeft_definitions.txt", upperLeftMMDefinitions];

ExportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_bottomRight.txt", bottomRightMMSymbolic];
ExportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_bottomRight_definitions.txt", bottomRightMMDefinitions];


(* ::Subsubsection:: *)
(*Construct scalar rotation matrix *)


blockSize = 6;
(** Diagonalizing rotation, this will be SO(N) with N=12. But we know a permutation transformation to reduce it to two SO(6) matrices,
so we first construct the two SO(6) and apply the inverse permutation. 
There's no easy way of generating a symbolic orthogonal matrix so just use a generic 6x6
This was done before DRalgo's fast rotate mode -TODO investigate fast rotate **)
rotUpperLeft = Table[ toIndexedSymbol2[ "RUL", i, j, Total[DigitCount[blockSize]] ], {i, 1, blockSize}, {j, 1, blockSize}];
rotBottomRight = Table[ toIndexedSymbol2[ "RBR", i, j, Total[DigitCount[blockSize]] ], {i, 1, blockSize}, {j, 1, blockSize}];

DSRotBlock = ArrayFlatten[{
{rotUpperLeft, 0},
{0, rotBottomRight}
}];
(* V = \[Phi]^T.M.\[Phi] 
	 = \[Phi]^T.P.P^T.M.P.P^T.\[Phi] = \[Phi]^T.P.B.P^T.\[Phi], make the mass matrix block diagonal, with some permutation matrix P: P^T.M.P = B
	 = \[Phi]^T.P.S.S^T.B.S.S^T.P^T.\[Phi] = \[Phi]^T.P.S.B.S^T.P^T.\[Phi], make the block diagonal mass matrix diagonal, with some similarity transform S: S^T.B.P = D'
Note P = scalarPermutationMatrix, S = DSRotBlock
Since we give DRalgo an arbitrary diagonal matrix and rotation matrix we have
V = \[Phi]^T.M.\[Phi]
  = \[Phi]^T.R.R^T.M.R.R^T.\[Phi] = \[Phi]^T.R.D.R^T.\[Phi]
We impose D = D' so R = P.S
We compute D' and S in the python code numerically
*)

DSRot = scalarPermutationMatrix . DSRotBlock;
Print["Scalar diagonalizing rotation:"];
DSRot//MatrixForm

ExportUTF8[effectivePotentialDirectory<>"/scalarRotationMatrix.txt", DSRot];

(** Diagonal mass matrix, unknown symbols **)
ScalarMassDiag = DiagonalMatrix[ Table[toIndexedSymbol["MSsq", i, Total[DigitCount[12]]], {i, 1, 12}] ];


(* ::Subsection::Closed:: *)
(*Gauge field diagonalization*)


(* ::Text:: *)
(*This is mostly hardcoded these for now. The diagonalization is identical to the Standard Model case*)


VectorMassMatrix = PrintTensorsVEV[2]//Normal;
vectorN = 12; (* SU3 x SU2 x U1 *)

(** Take the only nontrivial 2x2 submatrix and diagonalize that **) 
VectorMassMatrixNontrivial = VectorMassMatrix[[11;;12,11;;12]];
{VectorEigenvalues, VectorEigenvectors} = Eigensystem[VectorMassMatrixNontrivial];
VectorEigenvectors = FullSimplify[ Normalize /@ VectorEigenvectors, Assumptions -> {g3>0, g2>0, g1>0}];

(** Diagonalizing rotation: **)
DVRot = ArrayFlatten[{
	{IdentityMatrix[10], 0, 0},
	{0, 0, VectorEigenvectors}
}];

(* Eigenvalues in correct order. Just pick diagonals from the original mass matrix, then replace the ones we had to diagonalize manually  *)
VectorMassDiag = Table[ If[i==j, VectorMassMatrix[[i,i]], 0], {i, 1, 12}, {j, 1, 12}] //Simplify;
VectorMassDiag[[11,11]] = VectorEigenvalues[[1]];
VectorMassDiag[[12,12]] = VectorEigenvalues[[2]];

Print["Diagonalized vector mass matrix:"];
VectorMassDiag // MatrixForm

(** Simplify with easier symbols **)
gaugeRotationSubst = {g1/Sqrt[g1^2+g2^2] -> stW, g2/Sqrt[g1^2+g2^2] -> ctW};
DVRot = DVRot /. gaugeRotationSubst;

vectorShorthands = {stW-> g1/Sqrt[g1^2+g2^2], ctW-> g2/Sqrt[g1^2+g2^2]};

(** Vector masses mVsq[i]. **)
{VectorMassDiagSimple, VectorMassExpressions} = toSymbolicMatrix[VectorMassDiag, mVsq];

ExportUTF8[effectivePotentialDirectory<>"/vectorMasses.txt", VectorMassExpressions];
ExportUTF8[effectivePotentialDirectory<>"/vectorShorthands.txt", vectorShorthands];


(* ::Subsection::Closed:: *)
(*Export Veff derivatives*)


(* ::Text:: *)
(*Tree level minimization: \[PartialD]V/\[PartialD]Subscript[v, i] == 0, one for each background field.*)


(** DRalgo has no option for printing just the tree level potential!!
=> need to compute the full Veff, which is slow
=> cheat here with the PerturbativeDiagonalization flag
**)

(* LN: This didn't work! The PerturbativeDiagonalization seems to push off-diagonals terms like \[Phi]1^\[Dagger]\[Phi]2 to higher orders,
so the tree-level derivative is wrong! So let's just do the derivatives with sympy instead. *)
(*
CalculatePotentialUS[PerturbativeDiagonalization->True]
VeffLO = PrintEffectivePotential["LO"]//Expand;

extremaPolynomials = Table[D[VeffLO, backgroundFields[[i]]],{i, 1, Length[backgroundFields]}];
(* probs need to expand these for hpcpack to work *)
extremaPolynomials = extremaPolynomials // Simplify // Expand;
ExportUTF8[effectivePotentialDirectory<>"/extremaPolynomials.txt", extremaPolynomials];
*)


(* ::Subsection:: *)
(*Calculating the effective potential*)


(*(** NB! RotateTensorsCustomMass[] is very very slow, this can run for hours!
It's because our scalar rotation matrix is so large. **)
AbsoluteTiming[
	(** Tell DRalgo to rotate the fields to mass diagonal basis **)
	RotateTensorsCustomMass[DSRot,DVRot,ScalarMassDiag,VectorMassDiagSimple];
	CalculatePotentialUS[]
]*)


(*VeffLO = PrintEffectivePotential["LO"]//Simplify; (* Simplify to get rid of possible imaginary units *)
VeffNLO = PrintEffectivePotential["NLO"]//Simplify; (* Simplify to factor 1/pi division for tiny speed up *)
VeffNNLO = PrintEffectivePotential["NNLO"]; (* NOT simplified as seems to change numerical result for unknown reasons *)

ExportUTF8[effectivePotentialDirectory<>"/Veff_LO.txt", VeffLO];
ExportUTF8[effectivePotentialDirectory<>"/Veff_NLO.txt", VeffNLO];
ExportUTF8[effectivePotentialDirectory<>"/Veff_NNLO.txt", VeffNNLO];*)


(*$Assumptions = _Symbol \[Element] Reals
VeffNNLOSIMPREAL = PrintEffectivePotential["NNLO"]//Simplify;
ExportUTF8[effectivePotentialDirectory<>"/Veff_NNLOSIMPREAL.txt", VeffNNLOSIMPREAL];*)
