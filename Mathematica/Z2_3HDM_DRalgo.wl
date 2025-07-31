(* ::Package:: *)

(* ::Subsection:: *)
(*Import DRalgo and Group Math*)


ClearAll[]
SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
(* This is pointing to my DRalgo repo for easier updating *)
(*pathToDRalgo = "/home/lani/repos/DRalgo/DRalgo.m"*)
pathToDRalgo = "/home/jasmine/.Mathematica/Applications/DRalgo/DRalgo.m"
Get[pathToDRalgo]


(* ::Subsection::Closed:: *)
(*Import helper functions *)


Get["MathematicaToPythonHelper.m"]


(* ::Subsection:: *)
(*Specify file paths for exporting*)


hardToSoftDirectory = "DRalgoOutput/Z2_3HDM/HardToSoft";
softToUltrasoftDirectory = "DRalgoOutput/Z2_3HDM/SoftToUltrasoft";
effectivePotentialDirectory = "DRalgoOutput/Z2_3HDM/EffectivePotential";
variables = "DRalgoOutput/Z2_3HDM/Variables";


(* ::Section::Closed:: *)
(*Model*)


(*See 1909.09234 [hep-ph], eq (1) *)
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
(** Note that AllocateTensors[] brings some GroupMath symbols to the global namespace. And these are only removed later when calling ImportModelDRalgo...**)


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


(** Here just list all possible gauge-invariant operators containing 2 doublets **)
(** DRalgo notation is that \[Phi]1\[Phi]2^+ = \!\(
\*SubsuperscriptBox[\(\[Phi]\), \(2\), \(\[Dagger]\)]
\*SubscriptBox[\(\[Phi]\), \(1\)]\ in\ standard\ \(notation . \ \nSo\)\ careful\ here\ to\ make\ sure\ imaginary\ parts\ match\ to\ the\ potential\ in\ our\ draft\) **)
(** I have changed the notation a lot from DRalgo's example 3HDM file. 
My notation for doublet products is \[Phi]ij, where first index is the conjugated doublet **)

InputInv={{1,1},{True,False}}; (*\[Phi]1 \[Phi]1^\[Dagger]*)
\[Phi]11=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,2},{True,False}}; (*\[Phi]2 \[Phi]2^\[Dagger]*)
\[Phi]22=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,3},{True,False}}; (*\[Phi]3 \[Phi]3^\[Dagger]*)
\[Phi]33=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{1,2},{True,False}}; (*\[Phi]1\[Phi]2^\[Dagger]*)
\[Phi]21=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,1},{True,False}};(*\[Phi]2\[Phi]1^\[Dagger]*)
\[Phi]12=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{1,3},{True,False}}; (*\[Phi]1\[Phi]3^\[Dagger]*)
\[Phi]31=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,1},{True,False}};(*\[Phi]3\[Phi]1^\[Dagger]*)
\[Phi]13=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{2,3},{True,False}}; (*\[Phi]2\[Phi]3^\[Dagger]*)
\[Phi]32=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;

InputInv={{3,2},{True,False}};(*\[Phi]3\[Phi]2^\[Dagger]*)
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
\[CapitalLambda]4=GradQuartic[VQuartic];


InputInv={{3,1,2},{False,False,True}}; 
YukawaDoublet3=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;

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


(* ::Subsection:: *)
(*NLO matching, by which I mean Mode -> 2*)


(** Normalization4D flag = preserve 4D units so that the EFT path integral weight is e^{-S/T} (didn't work at time of writing)
 utoRG->True means that 3D running is built in to the matching. This is bad for automatization since 
the 3D masses become be functions of other 3D parameters. To dodge this we match with AutoRG->False
and do the RG running manually in an additional stage. **)
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False, Mode->2, Normalization4D->False, AutoRG->False];
PerformDRhard[];


betaFunctions4DUnsquared = BetaFunctions4D[] /. {(x_^2 -> y_) :> (x -> y/(2*x))};
exportUTF8[hardToSoftDirectory<>"/BetaFunctions4D.txt", betaFunctions4DUnsquared];


couplingsSoft = PrintCouplings[];
temporalScalarCouplings = PrintTemporalScalarCouplings[];
debyeMasses = PrintDebyeMass["LO"]; (** For Debyes we only take LO result, NLO not needed since we integrate these out anyway **)
scalarMasses = CombineSubstRules[PrintScalarMass["LO"], PrintScalarMass["NLO"]];
allSoftScaleParams = Join[couplingsSoft, temporalScalarCouplings, debyeMasses, scalarMasses];


(*DRalgo gives temporal couplings with [] which is a function call which makes things awkward so remove the []*)
\[Lambda]VL[i_]:=ToExpression["\[Lambda]VL"<>ToString[i]];
\[Lambda]VLL[i_]:=ToExpression["\[Lambda]VLL"<>ToString[i]];


(*We want to do in place updating of parameters in the python code i.e. \[Lambda]14D gets updated to \[Lambda]13D which gets updated to \[Lambda]13DUS,
it's easier to do this if we remove the suffices so its the same variable name throughout*)
allSoftScaleParamsSqrtSuffixFree = RemoveSuffixes[sqrtSubRules[allSoftScaleParams], {"3d"}];
(*Sometimes compute these equations and put the results into a np.zeros,
without T->T etc we would lose what T is *)
allSoftScaleParamsSqrtSuffixFree = Join[allSoftScaleParamsSqrtSuffixFree, {T->T,RGScale->RGScale}];
exportUTF8[hardToSoftDirectory<>"/softScaleParams_NLO.txt", allSoftScaleParamsSqrtSuffixFree];


(* 3D RG equations can be solved exactly, so do that here. We will export subst rules analogous to the matching relations:
	msq -> msq + \[Beta][msq] Log[\[Mu]3/\[Mu]] where RHS msq is the 3D mass at scale \[Mu] and LHS msq is the mass at scale \[Mu]3 *)
	
SolveRunning3D[betaFunctions_] := Block[{exprList},
	(* Extracting lhs and beta for each list element *)
	exprList = {#[[1]], #[[2]]} & /@ betaFunctions;

	(* Make new list with RGE solution on RHS *)
	newRulesList = (#1 -> #1 + #2*Log[T/RGScale]) & @@@ exprList;
	Return[newRulesList];
];


running3DSoft = RemoveSuffixes[SolveRunning3D[BetaFunctions3DS[]],{"3d"}];
running3DSoft= Join[running3DSoft, {RGScale->T}];
exportUTF8[hardToSoftDirectory<>"/softScaleRGE.txt", running3DSoft];


(* ::Subsection:: *)
(*Soft -> Ultrasoft matching*)


PerformDRsoft[{}];
(** This now works properly as of DRalgo 2023/11/24 update **)
couplingsUS = PrintCouplingsUS[];
scalarMassesUS = CombineSubstRules[PrintScalarMassUS["LO"], PrintScalarMassUS["NLO"]];
(*Change \[Mu]3 to RGScale for easier arraynesss in python (we previously did this explicitly in python)*)
allUltrasoftScaleParams = Join[couplingsUS, scalarMassesUS] /. \[Mu]3->RGScale;


allUltrasoftScaleParamsSqrt = RemoveSuffixes[sqrtSubRules[allUltrasoftScaleParams], {"US", "3d"}];(*Some reduant sqrt operations here? e.g. g13dUS*)
allUltrasoftScaleParamsSqrt= Join[allUltrasoftScaleParamsSqrt, {mu3US -> T}];


exportUTF8[softToUltrasoftDirectory<>"/ultrasoftScaleParams_NLO.txt", allUltrasoftScaleParamsSqrt];


runningUS = RemoveSuffixes[SolveRunning3D[BetaFunctions3DUS[]],{"US", "3d"}];
exportUTF8[softToUltrasoftDirectory<>"/ultrasoftScaleRGE.txt", runningUS];


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


(* ::Subsection:: *)
(*Diagonalizing scalar mass matrix*)


(* ::Text:: *)
(*Finding the diagonalizing matrix analytically is likely impossible, so we brute force this by giving DRalgo a symbolic 12x12 matrix with unknown symbols and compute the Veff in terms of those symbols. For masses we give an unknown diagonal matrix. These unknown are field-dependent so we must solve them numerically every time the potential is evaluated. *)
(**)
(*Note that using unknown symbols in the rotation matrix means that DRalgo will compute the effective potential in what is effectively a non-diagonal field basis, ie. there is quadratic mixing. But all the mixing effects vanish when numerical values are fixed by diagonalization conditions later on. (And I believe DRalgo ignores quadratic vertices anyway, but have not confirmed this). *)
(**)
(*Now, there are some optimizations that we can do. For example with the 3-field configuration defined above, the mass matrix can be brought into block-diagonal form by permuting the fields. This reduces the problem to diagonalization of two 6x6 matrices which is faster. We use this approach here; the required 12x12 rotation then has zeros in many places so DRalgo has easier time working with it.*)


scalarMM = PrintTensorsVEV[1]//Normal//Simplify; (* Scalar mass matrix, simplify to get rid of possible imaginary units *)


(* ::Subsubsection:: *)
(*Permute scalars to make mass matrix block-diagonal *)


(* Permutation matrix swaps the following rows and colums: 2<->11, 4<->9,6<->7  to get a block diagonal matrix. 
These are the minimal swaps to get a block diagonal matrix, and conviently leave the permutation symmetric.  
Once permuted the order of the mass matrix will be: 
Re\[Phi]1, Im\[Phi]3, Im\[Phi]1, Re\[Phi]3, Re\[Phi]2, Im\[Phi]2 (charged dof)
Re\[Phi]2, Im\[Phi]2, Im\[Phi]1, Re\[Phi]3, Re\[Phi]1, Im\[Phi]3  (neutral dof)*)
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
If[!OrthogonalMatrixQ[scalarPermutationMatrix], Print["Error, permutation matrix is not orthogonal"]];
exportUTF8[effectivePotentialDirectory<>"/scalarPermutationMatrix.txt", StringReplace[ToString[scalarPermutationMatrix],{"{"->"[","}"->"]"}]];


(*Our casescalarPermutationMatrix is symmetric but taking transpose anyway for consistency/future proofing*)
blockDiagonalMM = Transpose[scalarPermutationMatrix] . scalarMM . scalarPermutationMatrix;
Print["Block diagonal mass matrix:"];
blockDiagonalMM//MatrixForm

(*Extract permutation matrix and do consistency check*)
upperLeftMM = Take[blockDiagonalMM,{1,6},{1,6}];
bottomRightMM = Take[blockDiagonalMM,{7,12},{7,12}];

If[!SymmetricMatrixQ[upperLeftMM] || !SymmetricMatrixQ[bottomRightMM], Print["Error, block not symmetric!"]];


(* ::Subsubsection:: *)
(*Export scalar mass matrix*)


(* Simplify both blocks by introducing additional symbols, then extract them separately *)
{upperLeftMMSymbolic, upperLeftMMDefinitions} = toSymbolicMatrix[upperLeftMM, "MMUL", True]//Simplify;
{bottomRightMMSymbolic, bottomRightMMDefinitions} = toSymbolicMatrix[bottomRightMM, "MMBR", True]//Simplify;

(* Export expressions separately because we have the code to parse that *)
exportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_upperLeft.txt", upperLeftMMSymbolic];
exportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_upperLeft_definitions.txt", upperLeftMMDefinitions];


exportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_bottomRight.txt", bottomRightMMSymbolic];
exportUTF8[effectivePotentialDirectory<>"/scalarMassMatrix_bottomRight_definitions.txt", bottomRightMMDefinitions];


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
DSRot//MatrixForm;

exportUTF8[effectivePotentialDirectory<>"/scalarRotationMatrix.txt", DSRot];

(** Diagonal mass matrix, unknown symbols **)
ScalarMassDiag = DiagonalMatrix[ Table[toIndexedSymbol["MSsq", i, Total[DigitCount[12]]], {i, 1, 12}] ];


(* ::Subsection:: *)
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
VectorMassDiag // MatrixForm;

(** Simplify with easier symbols **)
gaugeRotationSubst = {g1/Sqrt[g1^2+g2^2] -> stW, g2/Sqrt[g1^2+g2^2] -> ctW};
DVRot = DVRot /. gaugeRotationSubst;

vectorShorthands = {stW-> g1/Sqrt[g1^2+g2^2], ctW-> g2/Sqrt[g1^2+g2^2]};

(** Vector masses mVsq[i]. **)
{VectorMassDiagSimple, VectorMassExpressions} = toSymbolicMatrix[VectorMassDiag, mVsq];

exportUTF8[effectivePotentialDirectory<>"/vectorMasses.txt", VectorMassExpressions];
exportUTF8[effectivePotentialDirectory<>"/vectorShorthands.txt", vectorShorthands];


(* ::Subsection:: *)
(*Calculating the effective potential*)


(** NB! RotateTensorsCustomMass[] is very very slow, this can run for hours!
It's because our scalar rotation matrix is so large. **)
AbsoluteTiming[
	(** Tell DRalgo to rotate the fields to mass diagonal basis **)
	RotateTensorsCustomMass[DSRot,DVRot,ScalarMassDiag,VectorMassDiagSimple];
	CalculatePotentialUS[]
]


veffLO = PrintEffectivePotential["LO"]//Simplify; (* Simplify to get rid of possible imaginaryDetailed units *)
veffNLO = PrintEffectivePotential["NLO"]//Simplify; (* Simplify to factor 1/pi division for tiny speed up *)
veffNNLO = PrintEffectivePotential["NNLO"]; (* NOT simplified as seems to change numerical result for unknown reasons *)


(*Done for consistent in out structure for python*)
veffLOR = {LO -> veffLO};
veffNLOR = {NLO -> veffNLO};
veffNNLOR = {NNLO -> veffNNLO};
exportUTF8[effectivePotentialDirectory<>"/Veff_LO.txt", veffLOR];
exportUTF8[effectivePotentialDirectory<>"/Veff_NLO.txt", veffNLOR];
exportUTF8[effectivePotentialDirectory<>"/Veff_NNLO.txt", veffNNLOR];


exportUTF8[
	variables<>"/LagranianSymbols.json", 
	{"fourPointSymbols"-> extractSymbols[\[CapitalLambda]4],
	"threePointSymbols"-> extractSymbols[\[CapitalLambda]3],
	"twoPointSymbols"-> extractSymbols[\[Mu]ij],
	"gaugeSymbols"-> extractSymbols[GaugeCouplings],
	"yukawaSymbols" -> extractSymbols[Ysff],
	"fieldSymbols" -> extractSymbols[backgroundFieldsFull]}];


(*The scalar mass matrices are a bit hacky since they don't have a direct out and the block diagonal nature means they shouldn't share the same out*)
equationSymbols={
	"hardScaleRGE"->{
		"Out" -> extractSymbols[betaFunctions4DUnsquared]["LHS"],
		"In" -> extractSymbols[betaFunctions4DUnsquared]["RHS"]},
	"softScaleParams"->{
		"Out" -> extractSymbols[allSoftScaleParamsSqrtSuffixFree]["LHS"],
		"In" -> extractSymbols[allSoftScaleParamsSqrtSuffixFree]["RHS"]},
	"softScaleRGE"->{
		"Out" -> extractSymbols[running3DSoft]["LHS"],
		"In" -> extractSymbols[running3DSoft]["RHS"]},	
	"ultraSoftScaleParams"->{
		"Out" -> extractSymbols[allUltrasoftScaleParamsSqrt]["LHS"],
		"In" -> extractSymbols[allUltrasoftScaleParamsSqrt]["RHS"]},
	"ultraSoftScaleRGE"->{
		"Out" -> extractSymbols[runningUS]["LHS"],
		"In" -> extractSymbols[runningUS]["RHS"]},	
	"upperLeftMMDefinitions"->{
		"Out" -> extractSymbols[ScalarMassDiag],
		"In" -> extractSymbols[upperLeftMMDefinitions]["RHS"]},
	"bottomRightMMDefinitions"->{
		"Out" -> extractSymbols[ScalarMassDiag],
		"In" -> extractSymbols[bottomRightMMDefinitions]["RHS"]},
	"vectorMasses"->{
		"Out" -> extractSymbols[VectorMassDiagSimple],
		"In" -> extractSymbols[VectorMassExpressions]["RHS"]},
	"LO"->{
		"Out" -> extractSymbols[veffLOR]["LHS"],
		"In" -> extractSymbols[veffLOR]["RHS"]},
	"NLO"->{
		"Out" -> extractSymbols[veffNLOR]["LHS"],
		"In" -> extractSymbols[veffNLOR]["RHS"]},
	"NNLO"->{
		"Out" -> extractSymbols[veffNNLOR]["LHS"],
		"In" -> extractSymbols[veffNNLOR]["RHS"]},
	"rotationSymbols"->extractSymbols[DSRot],
	"scalarMassNames"->extractSymbols[ScalarMassDiag]};


exportUTF8[variables<>"/EquationSymbols.json", equationSymbols];


exportUTF8[variables<>"/allSymbols.json",symbolsFromDict[equationSymbols]];


extractSymbols[ScalarMassDiag]
