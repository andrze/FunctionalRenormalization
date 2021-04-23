(* ::Package:: *)

(* ::Section:: *)
(*Numerical Library*)


(* ::Subsection::Closed:: *)
(*Constants*)


(* ::Text:: *)
(*Size of \[Rho] grid*)


If[!NumberQ[gridSize],
	gridSize = 120];
	
SetAttributes[gridSize,{Protected,Constant}];
potentialMinimum = Floor[57/60*gridSize];
Protect[potentialMinimum];


(* ::Text:: *)
(*Constants related to momentum integration (can be altered in your source code)*)


IntegralConfigurations = {
{0.001, 0.1, 3},
{0.1, 1., 5},
{1., 2., 5},
{2., 3., 4},
{3., 5., 3}};

SetAttributes[IntegralConfigurations,{Protected,Constant}]


(* ::Text:: *)
(*Point of normalization of Z[\[Rho]] function (Z[zNormalization] := 1)*)


If[!NumberQ[zNormalization],
	If[LPA,
		zNormalization = potentialMinimum;
		zNormalization = 0];
];	
SetAttributes[zNormalization,{Protected,Constant}];


(* ::Subsection:: *)
(*Replacement rules*)


(* ::Subsubsection::Closed:: *)
(*Coefficients for finite difference derivatives at order O(h^4)*)


(* ::Text:: *)
(*Forward/BackwardCoefficients - derivative at the boundary*)
(*SemiForward/BackwardCoefficients - derivative one step from the boundary*)
(*CentralCoefficcients - derivatives away from boundary*)


Block[{n},
CentralCoefficients[n_] := Piecewise[{{{1/12,-2/3,0,2/3,-1/12},n==1}, {{-1/12,4/3,-5/2,4/3,-1/12},n==2},{{-(1/2),1,0,-1,1/2},n==3}},{0}];
ForwardCoefficients[n_] := Piecewise[{{{-25/12,4,-3,4/3,-1/4},n==1}, {{15/4,-77/6,107/6,-13,61/12,-5/6},n==2},{{-(3/2),5,-6,3,-(1/2)},n==3}},{0}];
SemiForwardCoefficients[n_] := Piecewise[{{{-1/4,-5/6,3/2,-1/2,1/12},n==1}, {{5/6,-5/4,-1/3,7/6,-1/2,1/12},n==2},{{-(5/2),9,-12,7,-(3/2)},n==3}},{0}];
BackwardCoefficients[n_] := Piecewise[{{-ForwardCoefficients[n],OddQ[n]}}, ForwardCoefficients[n]];
SemiBackwardCoefficients[n_] := Piecewise[{{-SemiForwardCoefficients[n],OddQ[n]}}, SemiForwardCoefficients[n]]];


(* ::Input:: *)
(*(*n=5;*)
(*d=3;*)
(*rep = s[j_]\[Rule]j-3;*)
(*stencil = Table[s[j],{j,1,n}]/. rep*)
(*mat = Table[s[j]^i,{i,0,n-1},{j,1,n}]/.rep;*)
(*vec = d! Table[KroneckerDelta[d,i],{i,0,n-1}];*)
(*Inverse[mat].vec*)*)


(* ::Subsubsection::Closed:: *)
(*Rules for replacing derivatives with finite differences*)


(* ::Text:: *)
(*repder - replacement of the derivatives away from the boundary*)


repder = {};
Block[{n, a},
For[n=1, n<=3, n++, a=CentralCoefficients[n];
	AppendTo[repder,\!\(\*SuperscriptBox[\(ff_\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[ix_]-> Sum[a[[k]]*ff[ix+k-3],{k,1,Length[a]}]];
];]


(* ::Text:: *)
(*repderborder - replacement near the boundary*)


Block[{n, a, a2, b, b2},
repderborder = {};
For[n=1, n<=3, n++,
	a = ForwardCoefficients[n]; a2 = SemiForwardCoefficients[n]; 
	b = BackwardCoefficients[n]; b2 = SemiBackwardCoefficients[n];
	AppendTo[repderborder,\!\(\*SuperscriptBox[\(ff_\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[ix_] /; ix==0 -> Sum[a[[k]]*ff[ix+k-1], {k, 1, Length[a]}]];
	AppendTo[repderborder,\!\(\*SuperscriptBox[\(ff_\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[ix_] /; ix==1 -> Sum[a2[[k]]*ff[ix+k-2], {k, 1, Length[a2]}]];
	AppendTo[repderborder,\!\(\*SuperscriptBox[\(ff_\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[ix_] /; ix==gridSize -> Sum[b[[k]]*ff[ix-k+1], {k, 1, Length[b]}]];
	AppendTo[repderborder,\!\(\*SuperscriptBox[\(ff_\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[ix_] /; ix==gridSize-1 -> Sum[b2[[k]]*ff[ix-k+2],{k, 1, Length[b2]}]];
];]
AppendTo[repderborder,\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[0]->(56 f[0]-333 f[1]+852 f[2]-1219 f[3]+1056 f[4]-555 f[5]+164 f[6]-21 f[7])/6];


(* ::Subsubsection::Closed:: *)
(*Shortcuts for replacement*)


(* ::Text:: *)
(*NumericDerivatives - replaces exact derivatives with finite differences*)
(*ONSymmetry - simplifies equations to ones for the O (N) model*)
(*LPArep - simplifies equations to LPA' approximation*)


ONSymmetry = {Subscript[w, i_][\[Rho]_] -> 0, \!\(\*SuperscriptBox[
SubscriptBox[\(w\), \(i_\)], 
TagBox[
RowBox[{"(", "m_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0};
NumericDerivatives := (# //. repderborder //.repder)&
LPArep = {z[\[Rho]_]->1, \!\(\*SuperscriptBox[\(z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0, yy[\[Rho]_]->Y, \!\(\*SuperscriptBox[\(yy\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0, \!\(\*SuperscriptBox[\(t\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0, t[\[Rho]_] -> T};
ConstRep[dim_, \[Rho]max_, alpha_:defaultAlpha] := 
	{d->N[dim], vd->1,Subscript[zp, 1][0]->Subscript[zs, 1][0],
	eps->\[Rho]max/gridSize, \[Alpha]->alpha, z[zNormalization]->1, zp[zNormalization]->1, zs[zNormalization]->1};


(* ::Subsubsection::Closed:: *)
(*MakeLPA*)


(* ::Text:: *)
(*Yields integrands in LPA' approximation*)


MakeLPA[integrands_] := Block[{approximate},

Approximate[func_]:= Block[{approx},
	If[MatchQ[func[[1]], f_[0]],
		Return[Null];];
	If[MatchQ[func[[1]], f_[\[Rho]/eps]],
		approx = func/.LPArep;
		If[approx[[1]] == 1, approx[[1]] = Z];
		Return[approx];
		];
	Return[func];
];
Return[DeleteCases[Map[Approximate,integrands],Null]];
];


(* ::Subsection::Closed:: *)
(*Gauss - Legendre integration*)


(* ::Text:: *)
(*Generates weights and evaluation points Gauss - Legendre integration for n \[Element] [1, MaxPointsLegendre].*)


Block[{i, legendrePoints,roots,weights,rootsLegendre,weightGL},
MAXPOINTSLEGENDRE = 30;
rootsLegendre[n_] := N[Solve[LegendreP[n, x] == 0, x]]; (* generates roots of the n-th Legendre poly. *)
weightGL[n_, x_] := N[2/((1-x^2) D[LegendreP[n, x], x]^2)];
For[legendrePoints=1, legendrePoints <= MAXPOINTSLEGENDRE, legendrePoints++,
	roots = rootsLegendre[legendrePoints];
	weights = weightGL[legendrePoints, x] /. roots;
	GLTable[legendrePoints] = Table[{w[i] -> Chop[N[weights[[i]]]], N[roots[[i,1]]]},
		{i, 1, Length[rootsLegendre[legendrePoints]]}]
	];
];


(* ::Text:: *)
(*GLIntegral - carries out Gauss - Legendre integration of function on the interval [a, b] with n = pointsLegendre*)
(*GetIntegral - performs Gauss - Legendre integration of function split into parts specified in IntegralConfigurations*)


GLIntegral[f_,a_,b_,pointsLegendre_]:= Module[{i},(b-a)/2 Sum[(f[(a+b)/2 + 1/2 (b-a) x] w[i] /. GLTable[pointsLegendre][[i]]),{i, 1, pointsLegendre}]]


GetIntegral[f_] := Module[{i},Sum[GLIntegral[f,IntegralConfigurations[[i,1]],IntegralConfigurations[[i,2]],IntegralConfigurations[[i,3]]],{i,1,Length[IntegralConfigurations]}]]


(* ::Subsection::Closed:: *)
(*Regulator functions*)


SelectRegulator[name_]:=Block[{},
isCutoffLitim = False;
If[ MemberQ[name, {"theta", "litim"}],
	regulatorReplacement = {r -> (\[Alpha] (1-#) HeavisideTheta[1-#] &)};
	defaultAlpha = 1.;
	isCutoffLitim = True;
	regulatorLabel = "Litim";
	Return[];];
If[ name == "exponential",
	regulatorReplacement = {r -> (\[Alpha] Exp[-#] &)};
	defaultAlpha = 1.;
	regulatorLabel = "Exponential";
	Return[];];
If[ name == "smooth",
	regulatorReplacement = {r -> (\[Alpha] #/(Exp[#]-1) &)};
	defaultAlpha = 2.;
	regulatorLabel = "Wetterich";
	Return[];];

Print["Select supported regulator function"];
Quit[];
];


(* ::Subsection::Closed:: *)
(*Other functions*)


(* ::Subsubsection::Closed:: *)
(*Transpose Association*)


(* ::Text:: *)
(*Exchanges primary and secondary keys in a nested dictionary*)


TransposeAssociation[association_] := 
Block[{i, j, transposed, primaryKeys, pk, secondaryKeys, sk},
transposed = Association[];
secondaryKeys = Keys[association];
primaryKeys = Keys[association[secondaryKeys[[1]] ]];
For[i=1, i<=Length[primaryKeys], i++,
	pk = primaryKeys[[i]];
	transposed[pk] = Association[];
	For[j=1, j<=Length[secondaryKeys], j++,
		sk = secondaryKeys[[j]];
		If[MemberQ[Keys[association[sk]],pk],
			transposed[pk][sk] = association[sk][pk]];
		];
	];
Return[transposed];
];


(* ::Subsubsection::Closed:: *)
(*OneLevelNest*)


(* ::Text:: *)
(*Transforms a dictionary with list - like keys into a nested dictionary. One level nest - extracts last value from a list-like key and organizes it into a nested dictionary*)


OneLevelNest[association_] := 
Block[{i, newEntry, nested, key, k, rest},
nested = Association[];
If[Length[Keys[association]]==0 || Length[Keys[association][[1]]] ==0,
	Return[association]]; 
For[i=1, i<=Length[association], i++,
	key = Keys[association][[i]];
	k = key[[-1]];
	rest = key[[1;;-2]];
	If[Length[rest] == 1, rest = rest[[1]]];
	newEntry = Association[k->association[[i]]];
	If[MemberQ[Keys[nested], rest],
		nested[rest] = KeySort[Join[nested[rest],newEntry],#1>#2&],
		nested[rest] = newEntry];
];

Return[nested];
];


(* ::Subsubsection::Closed:: *)
(*NestAssociation*)


(* ::Text:: *)
(*Transforms a dictionary with list - like keys into a nested dictionary*)


NestAssociation[association_] := 
Block[{i, newEntry, nested, key, k, rest},

nested = association;
While[Length[Keys[nested]]>0 && Length[Keys[nested][[1]]] > 0,
	nested = KeySort[OneLevelNest[nested],#1>#2&];
];

Return[nested];
];


(* ::Subsubsection::Closed:: *)
(*UnnestAssociation*)


(* ::Text:: *)
(*Combines keys of a nested dictionary into first-level dictionary*)


UnnestAssociation[association_] := 
Block[{i, j, keysLevels, level, newlevel, newkey, pk, sk},
keysLevels = {};

level = association;
While[AssociationQ[level[[1]]],
	newlevel = Association[];
	For[i=1, i<=Length[level], i++,
		pk = Keys[level][[i]];
		For[j=1, j<=Length[level[[i]]], j++,
			sk = Keys[level[pk]][[j]];
			newkey = Flatten[{pk,sk}];
			newlevel[newkey] = level[pk][sk];
		];
	];
	level = newlevel;
];

Return[newlevel];
];


(* ::Subsubsection::Closed:: *)
(*SelectClosest*)


(* ::Text:: *)
(*Selects value in association 'assoc' corresponding to the key with a value closest to that of 'key'. If no key in 'assoc' is within 0.05 from 'key' Null is returned.*)


SelectClosest[assoc_, key_] := Block[{sorted, keys, diffs, pos},
  sorted = KeySort[assoc];
  keys = Keys[sorted];
  diffs = keys - key;
  pos = Position[Abs[diffs], Min[Abs[diffs]]][[1, 1]];
  If[pos > 1 && Sign[diffs[[pos]]]*Sign[diffs[[pos - 1]]] == -1, Return[sorted[[pos]]];];
  If[pos < Length[sorted] && Sign[diffs[[pos]]]*Sign[diffs[[pos + 1]]] == -1, Return[sorted[[pos]]];];
  If[Abs[diffs[[pos]]] < 0.05, Return[sorted[[pos]]]];
  Return[Null];
  ]


(* ::Subsubsection::Closed:: *)
(*DropOutliers*)


(* ::Text:: *)
(*Drops keys associated with values outside the specified range*)


DropOutliers[assoc_,min_,max_] := Block[{filter, a=Association[]},
	filter[i_] := If[assoc[[i]]>=min && assoc[[i]]<=max, a[Keys[assoc][[i]]] = assoc[[i]]];
	Map[filter, Range[Length[assoc]]];
	Return[a];
];


(* ::Section:: *)
(*Renormalization Group*)


(* ::Subsection::Closed:: *)
(*Integration of the flow equations*)


(* ::Text:: *)
(*Initializes the list of flow equations and list of model parameters (without substituting the constants).*)


IntegrateEquations[integrandsList_]:=
Block[{fff,i, j, equationsList={}, equation, parameters={}, function, scaling, integrand, integrand0, minPoint, reppedEquation},
For[i=1, i<=Length[integrandsList], i++,
	(* We identify the components of integrandsList *)
	function = integrandsList[[i, 1]];
	scaling = integrandsList[[i, 2]];
	integrand[q_] = integrandsList[[i, 3]] /. regulatorReplacement /. y -> q;

	(* Full flow equation consists of scaling term and integrated term *)
	equation = scaling + GetIntegral[integrand];

	
	eqRep[j_] := Block[{},
		If[j==0,
			(* If j\[Equal]0 we load a separate \[Rho]=0 evaluation of the integrand *)
			integrand0[q_] = integrandsList[[i, 4]] /. regulatorReplacement /. y -> q;
			reppedEquation = (scaling /.{\[Rho]->0}) + GetIntegral[integrand0],
	
			(* If not we just substitute \[Rho] *)
			reppedEquation = equation /. \[Rho] -> j eps
		];

		(* We append equation to the list of equations and function to the list of parameters*)
		AppendTo[equationsList, reppedEquation];
		AppendTo[parameters, function/. \[Rho] -> j eps];
	];
	
	(* Iteration over \[Rho] grid *)
	(* If `function` is not a function of \[Rho] but constant then no iteration is performed *)
	If[MatchQ[function,fff_[\[Rho]_]],
		minPoint = 1;
		If[Length[integrandsList[[i]]]==4, minPoint=0];
		Map[eqRep,Range[minPoint,gridSize]],
		eqRep[zNormalization];
	];	
];
(* Replaces derivatives with finite differences *)
equationsList = NumericDerivatives[equationsList];
Return[equationsList];
];


(* ::Subsection:: *)
(*Fixed Point Search*)


(* ::Subsubsection::Closed:: *)
(*AdjustNumberOfPoints *)


(* ::Text:: *)
(*Returns a guess interpolated from supplied one onto grid with specified number of points together with new position for zNormalization point. *)


AdjustNumberOfPoints[guess_, numberOfPoints_, zNormalization_, \[Rho]Max_]:= Block[{guessRepar,newGuess, scale, functions, 
	guessLength, Interpolate,zDeviation, oldFunctions,interpolation,Reparametrize},
functions = {v,zs,zp};
oldFunctions = DeleteCases[DeleteDuplicates[ guess[[All,1,0]]],Symbol];
guessLength = Length[Position[guess,v[i_]]]-1;

guessRepar = guess;
If[oldFunctions == {v,z,yy},
	Reparametrize[f_] := Block[{j},
		If[MatchQ[f,v[i_]->c_?NumberQ], Return[f]];
		
		If[MatchQ[f,z[i_]->c_?NumberQ], Return[zs[f[[1,1]]]->f[[2]] ]];
		
		If[MatchQ[f,yy[i_]->c_?NumberQ], 
			j = f[[1,1]];
			Return[zp[j]->(z[j]-2 j eps yy[j])/. guess /. {z[i_]->1, eps->\[Rho]Max/guessLength}];
		];
		Return[f];];
	guessRepar = Map[Reparametrize,guess];
];

(* Function for linear interpolation onto new grid *)
scale = guessLength/numberOfPoints;

(* Construction of a new grid representation *)
interpolation := Map[Interpolation[Table[{j, #[j]},{j,0,guessLength}] /. guessRepar /. {zp[k_] ->1,zs[k_]->1},InterpolationOrder->4]&,functions];

newGuess = Table[functions[[i]][j] -> interpolation[[i]][j*scale]/(zs[zNormalization]/.guessRepar/.zs[k_]->1), {i,1,Length[functions]}, {j,0,numberOfPoints}];
newGuess = Flatten[newGuess,1];

(* Remove zNormalization from FP specification *)
newGuess = DeleteCases[newGuess, zs[zNormalization]->c_];
newGuess = DeleteCases[newGuess, zp[0]->c_];

(* Add anomalous dimension \[Eta] *)
AppendTo[newGuess,\[Eta]->(\[Eta]/.guess)];

Return[newGuess];
]


MatchQ[f->1,f_->c_?NumberQ]


fu = v[x]->1;


zs[fu[[1,1]]]->fu[[2]]


v[x][[1]]


(* ::Subsubsection::Closed:: *)
(*LPAGuess*)


LPAGuess[guess_] := Block[{},
Return[Join[guess[[1;;gridSize+1]],{\[Eta]->(\[Eta]/.guess)}]]];


(* ::Subsubsection::Closed:: *)
(*LPApGuess*)


LPApGuess[guess_] := Block[{},
Return[Join[guess[[1;;gridSize+1]],{Y->yy[zNormalization]/.guess,\[Eta]->(\[Eta]/.guess)}]]];


(* ::Subsubsection::Closed:: *)
(*ConvertGuess*)


(* ::Text:: *)
(*Converts guess from Z, Y parametrization to Zs, Zp parametrization*)


ConvertGuess[guess_, \[Rho]Max_]:=Block[{vs,zss,zps,newGuess},
vs = Table[v[i]->(v[i]/.guess),{i,0,gridSize}];
zss = Table[zs[i]->(z[i]/.guess),{i,1,gridSize}];
zps = Table[zp[i]->((z[i]-2 i eps yy[i])/.guess/.ConstRep[x,\[Rho]Max]),{i,1,gridSize}];
newGuess = Join[vs,zss,zps,{\[Eta]->(\[Eta]/.guess)}];
Return[newGuess];
];


(* ::Subsubsection::Closed:: *)
(*PlotFixedPoint*)


(* ::Text:: *)
(*Plots a ListPlot for a list of replacements*)


PlotFixedPoint[fixedPoint_] :=Block[{fp=fixedPoint},
If[fixedPoint[[0]] == Association,
	If[Length[fixedPoint]==1,
		Return[PlotFixedPoint[fixedPoint[[1]],Keys[fixedPoint][[1]]]],
		Print["Argument is not a valid model parametrization"];
		Return[Null];]];
If[MemberQ[Keys[fixedPoint],Y], AppendTo[fp,yy[k_]->(Y/.fp)];];
ListPlot[Transpose[Table[{v[i],z[i],yy[i]},{i,0,gridSize}]/.fp/.z[k_]->1], 
	PlotLegends->{"V","Z","Y"},PlotLabel->"\[Eta]="<>ToString[\[Eta]/.fp],AxesLabel->{"i","f[i]"}]];


PlotFixedPoint[guess_,\[Rho]Max_] := Block[{tabelize, functions= {v,zs,zp}, fp=guess, oldFunctions, guessLength, guessRepar,Reparametrize},
oldFunctions = DeleteCases[DeleteDuplicates[ guess[[All,1,0]]],Symbol];
guessLength = Length[Position[guess,v[i_]]]-1;

If[MemberQ[fp[[All,1]],Y], AppendTo[fp,yy[k_]->(Y/.fp)];];

guessRepar = guess;
If[oldFunctions == {v,z,yy},
	Reparametrize[f_] := Block[{j},
		If[MatchQ[f,v[i_]->c_?NumberQ], Return[f]];
		
		If[MatchQ[f,z[i_]->c_?NumberQ], Return[zs[f[[1,1]]]->f[[2]] ]];
		
		If[MatchQ[f,yy[i_]->c_?NumberQ], 
			j = f[[1,1]];
			Return[zp[j]->(z[j]-2 j eps yy[j])/. guess /. {z[i_]->1, eps->\[Rho]Max/guessLength}];
		];
		Return[f];];
	guessRepar = Map[Reparametrize,guess];
];

tabelize[f_] := Table[{i \[Rho]Max/guessLength, f[i]},{i,0,guessLength}]/.guessRepar/.{z[k_]->1,zs[k_]->1,zp[k_]->1};
functions = Map[tabelize,functions];
ListPlot[functions, PlotLegends->{"V",Subscript["Z","\[Sigma]"],Subscript["Z","\[Pi]"]},AxesLabel->{"\[Rho]","f[\[Rho]]"},PlotLabel->"\[Eta]="<>ToString[\[Eta]/.fp]]];


(* ::Subsubsection:: *)
(*FindFixedPoint*)


(* ::Text:: *)
(*Finds a fixed point for provided equations*)


FindFixedPoint[eqlist_, guess_, constants_,damping_:2]:=
Block[{solution, eqlistSubstituted, plot, check, guessrep},
FindFixedPoint::singularEq = "Flow equations singular due to too small regulator prefactor";

If[(\[Alpha] /.constants) < -(v[0]/.guess),
	Message[FindFixedPoint::singularEq];
	Return[guess],
	
	(* Substitutes for constants in the equations *)
	eqlistSubstituted = eqlist /. constants;

	guessrep = Table[{guess[[i,1]], guess[[i,2]]}, {i,1,Length[guess]}];

	(* Finds numerical solution *)
	solution = Quiet[FindRoot[eqlistSubstituted, guessrep,
					MaxIterations->10000, PrecisionGoal->5, DampingFactor->damping],
				{FindRoot::lstol}];
	
	Return[solution];
];
]


(* ::Subsubsection::Closed:: *)
(*IsFixedPoint *)


(* ::Text:: *)
(*Checks if guess is a fixed point for provided equations*)


IsFixedPoint[eqlist_, guess_, constants_]:=
Block[{eqlistSubstituted,check},
eqlistSubstituted = eqlist /. constants;
(* Finds numerical solution *)
check = eqlistSubstituted /. guess;
Print["Check for d="<>ToString[d/.constants]<>" Precision digits:"<>ToString[Floor[-Log10[Norm[check]]]]];
Return[Norm[check]<10^-5];
]


(* ::Subsubsection::Closed:: *)
(*IsWilsonFisher *)


(* ::Text:: *)
(*Checks if guess is proper Wilson-Fisher FP*)


IsWilsonFisher[guess_]:=
Block[{minmax},

minmax = MinimaMaxima[guess];
If[minmax["Minima"][[1]] != 1 || minmax["Maxima"][[1]] != 1, Return[False]];

If[(\[Eta]/.guess) <0, Return[False]];

Return[True];
]


(* ::Subsubsection::Closed:: *)
(*IsTricritical *)


(* ::Text:: *)
(*Checks if guess is a tricritical point*)


IsTricritical[guess_]:=
Block[{minmax},
minmax = MinimaMaxima[guess];
If[minmax["Minima"][[1]] != 2 || minmax["Maxima"][[1]] != 1, Return[False]];

If[(\[Eta]/.guess) <0, Return[False]];

Return[True];
]


(* ::Subsubsection::Closed:: *)
(*PotentialMinimum*)


(* ::Text:: *)
(*Finds the first minimum of the potential*)


PotentialMinimum[guess_]:=
Block[{i},
For[i=0,i<gridSize,i++,
	If[(v[i]/.guess) <0 && (v[i+1]/.guess) > 0,
		Return[i];
		];
	];
Return[Null];
];


(* ::Subsubsection::Closed:: *)
(*MinimaMaxima*)


(* ::Text:: *)
(*Finds the all minima and maxima of the potential (as function of order parameter)*)


MinimaMaxima[guess_]:=
Block[{u, i, min = {0, {}}, max = {0, {}}},
For[i=0, i<=gridSize, i++, u[i] = v[i]/.guess];

If[u[0] > 0, 
	min[[1]] = 1;
	AppendTo[min[[2]],0],
	max[[1]] = 1;
	AppendTo[max[[2]],0]];

For[i=1,i<gridSize,i++,
	If[u[i] <0 && u[i+1] > 0,
		min[[1]] += 1;
		AppendTo[min[[2]], i];
		];
	
	If[u[i] >0 && u[i+1] < 0,
		max[[1]] += 1;
		AppendTo[max[[2]], i];
		];
	];
Return[Association["Minima"->min, "Maxima"->max]];
];


(* ::Subsubsection::Closed:: *)
(*\[Rho]MaxScan*)


(* ::Text:: *)
(*Returns an association with entries \[Rho]max -> FP. Starts at initial\[Rho]max and scans with additive increments in \[Rho] up to maxIterations.*)
(*On approaching inproper FPs (see IsWIlsonFisher function) it discards it, but does not stop the iterations. Iterations are stopped after 5 consecutive errors.*)


\[Rho]MaxScan[equations_, initialFixedPoint_, dim_, initial\[Rho]Max_, alpha_, \[Rho]Increment_, maxIterations_] := Block[
{d\[Rho], min, j, guess, newFixedPoint, fixedPoints, fail, \[Rho]Max, desiredMinimumPosition=potentialMinimum, FPCheck=IsWilsonFisher, fpMsg, errorMsg},

If[searchTricritical, FPCheck=IsTricritical];

fpMsg[d_,\[Rho]_,a_,m_] := "FP at d="<>ToString[d]<>" \[Rho]max="<>ToString[\[Rho]]<>" \[Alpha]="<>ToString[a]<>" Potential minimum at "<>ToString[m]<>"/"<>ToString[gridSize];
errorMsg[d_, \[Rho]_] := "Error during calculation of a Fixed Point at d="<>ToString[d]<>" for \[Rho]max="<>ToString[\[Rho]];

\[Rho]Max = initial\[Rho]Max;

fail=False;
Check[
	newFixedPoint = FindFixedPoint[equations, initialFixedPoint, ConstRep[dim, \[Rho]Max, alpha]],
	fail = True;
	];
	
min = PotentialMinimum[newFixedPoint];
If[fail || !FPCheck[newFixedPoint],
	Print[errorMsg[dim, \[Rho]Max]];
	Return[{}],
	Print[fpMsg[dim, \[Rho]Max, alpha, min]];
	];

If[min<desiredMinimumPosition, 
	d\[Rho] = -Abs[\[Rho]Increment],
	d\[Rho] = Abs[\[Rho]Increment]];

\[Rho]Max += d\[Rho];
guess = newFixedPoint;
For[j=0, j<maxIterations, j++,
	fail=False;
	Check[
		newFixedPoint = FindFixedPoint[equations, guess, ConstRep[dim, \[Rho]Max, alpha]],
		fail = True;
		];
		
	fail = fail || !FPCheck[newFixedPoint];
		
	If[fail, 
		Print[errorMsg[dim,\[Rho]Max]];
		Return[{}],(*Return[{guess, \[Rho]Max-d\[Rho]}],*)
		
		min = MinimaMaxima[newFixedPoint]["Minima"][[2,-1]];
		Print[fpMsg[dim, \[Rho]Max, alpha, min]];
		If[d\[Rho]<0 && PotentialMinimum[newFixedPoint] >= desiredMinimumPosition, Return[{newFixedPoint, \[Rho]Max}]];
		If[d\[Rho]>0 && PotentialMinimum[newFixedPoint] <= desiredMinimumPosition, Return[{guess, \[Rho]Max-d\[Rho]}]];
		guess = newFixedPoint;
		];
	\[Rho]Max += d\[Rho];
	];

Return[{}];(*Return[{guess, \[Rho]Max-d\[Rho]}];*)
];


(* ::Subsubsection::Closed:: *)
(*FindFixedPointDict*)


(* ::Text:: *)
(*Creates a dictionary of dictionaries of fixed points.*)
(*The structure of outcome is like  (Replacement list FP) = dict[dimension][\[Rho]Max]*)


FindFixedPointDict[eqlist_, guess_, keys_, keyLabels_, constants_] := 
Block[{fixedPointDict, fixedPoints, fixedPoints2, params, equations, newFixedPoint, newGuess, fail, i, j, scan, 
	\[Rho]Max, oldPK, pkChanged=False, newSeries\[Rho]Max, newSeriesGuess},
fixedPointDict = Association[];
newGuess = guess;

If[searchTricritical,
	Print["Searching for tricritical fixed points"];
	Print["Searching for critical fixed points"]];

If[Length[keyLabels] > 1,
	oldPK = keys[[1,1]];
	];

For[i=1, i<=Length[keys], i++,
	params = constants;
	If[ValueQ[\[Rho]Max], params["\[Rho]Max"] = \[Rho]Max];
	
	If[Length[keyLabels] == 1,
		params[keyLabels[[1]]] = keys[[i]],
		
		pkChanged = (oldPK != keys[[i,1]]);
		For[j=1, j<=Length[keyLabels], j++,
			params[keyLabels[[j]]] = keys[[i,j]];
		];
		
		If[pkChanged,
			params["\[Rho]Max"] = newSeries\[Rho]Max;
			newGuess = newSeriesGuess;
		];
	];
	
	equations = eqlist;
	If[MemberQ[Keys[params],"N"], equations = equations /. n -> params["N"]];
	Print[params];
	
	scan = \[Rho]MaxScan[equations, newGuess, params["dim"], params["\[Rho]Max"], params["alpha"], 0.025*params["\[Rho]Max"], 20];
	
	If[Length[scan]<2 || !IsFixedPoint[equations, scan[[1]], ConstRep[params["dim"], scan[[2]], params["alpha"]]], 
		Print["No FP found for "<> ToString[params]];
		If[pkChanged || Length[keyLabels]==1,
			Break[],
			Continue[]];
		];
	newGuess = scan[[1]];
	\[Rho]Max = scan[[2]];

	If[pkChanged || i==1,
		newSeriesGuess = newGuess;
		newSeries\[Rho]Max = \[Rho]Max];

	If[Length[keyLabels] > 1, oldPK = keys[[i,1]]];
	
	If[Length[keyLabels] == 1,
		fixedPointDict[{keys[[i]],\[Rho]Max}] = newGuess,
		fixedPointDict[Join[keys[[i]],{\[Rho]Max}]] = newGuess;
		];
	
];
	
Return[fixedPointDict];
];


(* ::Subsubsection::Closed:: *)
(*FixedPointDimScan*)


(* ::Text:: *)
(*Creates a dictionary of dictionaries of fixed points.*)
(*The structure of outcome is like  (Replacement list FP) = dict[dimension][\[Rho]Max]*)


FixedPointDimScan[eqlist_, guess_, constants_, d_, initial\[Rho]Max_] := 
Block[{fixedPointDict=Association[], newGuess=guess, scan,
	\[Rho]Max=initial\[Rho]Max, step=0.05, steps=0, maxSteps=100, dim=d, factor=1.01},
If[searchTricritical,
	Print["Searching for tricritical fixed points"],
	Print["Searching for critical fixed points"]];

While[True,
	steps++;
	scan = \[Rho]MaxScan[eqlist, newGuess, dim, \[Rho]Max*factor, constants["alpha"], 0.0125*\[Rho]Max, 20];
	
	If[Length[scan]<2 || !IsFixedPoint[eqlist, scan[[1]], ConstRep[dim, scan[[2]], constants["alpha"]]], 
		Print["No FP found for d="<> ToString[dim]];
		If[step > 10^-4,
			step /= 2;
			dim += step;
			Continue[],
			Break[]];
		];
	factor = scan[[2]]/\[Rho]Max;
	newGuess = scan[[1]];
	\[Rho]Max = scan[[2]];
	fixedPointDict[{dim,\[Rho]Max}] = newGuess;
	dim -= step;
	If[steps>=maxSteps,Break[]];
];
	
Return[fixedPointDict];
];


(* ::Subsubsection::Closed:: *)
(*Rescale\[Rho]*)


Rescale\[Rho][fp_,\[Rho]Max_]:= Block[{min, new, newfp,funcs,points,Inter,Guess},
	points = Table[i \[Rho]Max/gridSize,{i,0,gridSize}];
	funcs = DeleteCases[DeleteDuplicates[ fp[[All,1,0]]],Symbol];
	Inter[f_] := Interpolation[Transpose[{points,Map[f,Range[0,gridSize]]}]/.fp/.h_[k_]->1.];
	
	Quiet[min = \[Rho]/. FindRoot[Inter[v][\[Rho]],{\[Rho],\[Rho]Max*potentialMinimum/gridSize}]];
	new = min*gridSize/potentialMinimum;
	Guess[f_] := Block[{int,start=1},
		int = Inter[f];
		If[f==v,start=0];
		Return[Table[f[i]->Quiet[int[i*new/gridSize]],{i,start,gridSize}]];
	];
	newfp = Flatten[Map[Guess,funcs]];
	Return[{Append[newfp,\[Eta]->(\[Eta]/.fp)],new}];
]


(* ::Subsubsection::Closed:: *)
(*MakeDStep*)


MakeDStep[eqs_, oldFP_, \[Rho]Max_, dim_, dd_, alpha_] := Block[{der, sm, dG, newFP, new\[Rho],zpart, zpartdG, bound, delta, exps},
	der = D[eqs,d]/.oldFP/.ConstRep[dim,\[Rho]Max,alpha];
	der = Drop[der,{gridSize+2}];
	sm = StabilityMatrix[integrandsListIsotropic,oldFP,ConstRep[dim,\[Rho]Max,alpha]];
	exps = ExtractExponents[sm];
	exps[\[Eta]] = \[Eta]/.oldFP;
	dG = -Quiet[Inverse[sm]].(der /. oldFP /. ConstRep[dim, \[Rho]Max, alpha]);
	AppendTo[dG,0.]; (* Append no change to \[Eta] *)
	
	
	newFP = oldFP; (* zpart includes also the last point of V *)
	zpart = oldFP[[gridSize+1;;-2]];
	zpartdG = dG[[gridSize+1;;-2]];
	
	(* delta must be small enough that at any point z is decreasing not more than 1/16. of its value *)
	bound = zpart/(16 zpartdG); 
	bound[[1]] = Abs[bound[[1]]]; 
	bound = Min[Select[bound, #>0&]];
	delta=dd;
	While[delta > bound, delta/=2];
	
	newFP[[All,2]] = oldFP[[All,2]] - delta dG;
	
	{newFP,new\[Rho]} = Rescale\[Rho][newFP,\[Rho]Max];
	
	newFP = FindFixedPoint[eqs, newFP, ConstRep[dim-delta, new\[Rho], alpha], 4];
	Return[{exps, newFP, new\[Rho], dim-delta}];
];


(* ::Subsubsection:: *)
(*FullDScan*)


FullDScan[eqlist_, guess_, d_, dd_, initial\[Rho]Max_, alpha_]:=
Block[{fixedPointDict=Association[], expDict=Association[], oldFP=guess, newFP, scan, exps,
	\[Rho]Max=initial\[Rho]Max, steps=0, maxSteps=300, newdim, dim=d, delta, factor=1.01, new\[Rho]},


For[steps=0, steps<maxSteps, steps++,
	If[dim>=1.35, delta=dd, delta=0.01];
	{exps, newFP, new\[Rho], newdim} = MakeDStep[eqlist, oldFP, \[Rho]Max, dim, delta, alpha];
	
	fixedPointDict[{dim,\[Rho]Max}] = oldFP;
	expDict[dim] = exps;
	If[!IsFixedPoint[eqlist,newFP,ConstRep[newdim,new\[Rho],alpha]],
		Break[]];
	
	oldFP=newFP;
	dim=newdim;
	\[Rho]Max=new\[Rho];
];
	
Return[{fixedPointDict, TransposeAssociation[expDict]}];
];


(* ::Subsection:: *)
(*Stability Matrix*)


(* ::Subsubsection::Closed:: *)
(*Stability Matrix*)


(* ::Text:: *)
(*Returns the stability matrix of a system. M_i,j is a derivative of the i-th flow equation with respect to j-th free parameter of the model (evaluated at the fixed point). *)


StabilityMatrix[integrandsList_, fixedPoint_, constants_] := 
Block[{params, fixedPointRep, anisotropyOrder, i, j,
	FPGrad, FPSeries, ZIndex, \[Eta]Flow, \[Eta]Equation, \[Eta]Integrand, \[Eta]Solution, \[Eta]Replacement, minPoint,
	function, integrand, scaling, scalingGradient, integrandPoint, stabilityMatrix, eqlist, check},

(* List of free parameters in a model *)
minPoint = Map[4-Length[#]&,integrandsList];
params = Flatten[DeleteDuplicates[Table[integrandsList[[i,1]]/.{\[Rho]->j eps},{i,1,Length[integrandsList]},{j,minPoint[[i]],gridSize}]]];

params = DeleteCases[params, z[zNormalization]];
params = DeleteCases[params, zp[zNormalization]];
params = DeleteCases[params, zs[zNormalization]];
params = DeleteCases[params, Z];

(* Fixed Point replacement including only free parameters and cutoff for expansion in anisotropic field *)
fixedPointRep = Table[params[[i]]->(params[[i]] /. {t[0]->NumericDerivatives[zs'[0]-zp'[0]] /(4 eps)}/. {t[k_] -> (zs[k]-zp[k])/(4 k eps)} /. fixedPoint /. 
{Subscript[f_, j_][k_] -> 0}/. constants ),
		{i,1,Length[params]}];
AppendTo[fixedPointRep, Subscript[f_, ii_][jj_]->0];
If[!MemberQ[params,t[0]],
	AppendTo[fixedPointRep, t[0]->(NumericDerivatives[zs'[0]-zp'[0]] /(4eps)/. fixedPoint/. constants)];
];

(* Gradient operator in model parameter space *)
FPGrad[f_] := Grad[f, params]/. fixedPointRep;
(* Linear expansion around the fixed point *)
FPSeries[f_] := (f /. fixedPointRep) + (params -(params/.fixedPointRep)).FPGrad[f];

(* Anomalous dimension in terms of other parameters *)
\[Eta]Replacement = {\[Eta]->0};
ZIndex = Position[integrandsList[[All, 1]], z[\[Rho]/eps]];
If[Length[ZIndex]==0,
	ZIndex = Position[integrandsList[[All, 1]], Z];];
If[Length[ZIndex]==0,
	ZIndex = Position[integrandsList[[All, 1]], zs[\[Rho]/eps]];];

If[Length[ZIndex]>0,
	ZIndex = ZIndex[[1,1]];
	\[Eta]Flow = integrandsList[[ZIndex]];
	If[zNormalization==0,
		\[Eta]Integrand[q_] = \[Eta]Flow[[4]] /. regulatorReplacement /. y -> q,
		\[Eta]Integrand[q_] = \[Eta]Flow[[3]] /. regulatorReplacement /. y -> q];
	\[Eta]Equation = \[Eta]Flow[[2]] + GetIntegral[\[Eta]Integrand];
	\[Eta]Equation = NumericDerivatives[\[Eta]Equation /. \[Rho] -> zNormalization eps];
	\[Eta]Solution = Solve[(\[Eta]Equation /. constants)==0, \[Eta]] // Flatten;
	\[Eta]Replacement = { \[Eta] -> FPSeries[\[Eta] /. \[Eta]Solution] };];

(* Derivatives of equations*)
stabilityMatrix = {};
For[i=1, i<=Length[integrandsList], i++,
	function = integrandsList[[i,1]];
	integrand[q_] = integrandsList[[i,3]] /. regulatorReplacement /. y -> q;
	scaling = integrandsList[[i,2]];

	AddRow[j_] := Block[{},
		If[!MemberQ[params,function/.\[Rho]->j eps], Return[]];
		(*If[Mod[j,10]==0,Print[function/.\[Rho]->j eps]];*)
		If[j==0,
			(* Y[\[Rho]] has separate formula for integrand at \[Rho]=0 *)
			integrandPoint[q_] = FPGrad[NumericDerivatives[integrandsList[[i,4]] 
				/. regulatorReplacement /. y -> q] /. \[Eta]Replacement /. constants],
			(* Takes gradient of the flow equation *)
			integrandPoint[q_] = FPGrad[NumericDerivatives[integrand[q] /. \[Rho] -> j eps]
				/. \[Eta]Replacement /. constants];
			];

		scalingGradient = NumericDerivatives[scaling/. \[Rho]->j eps] /. \[Eta]Replacement /. constants;

		AppendTo[stabilityMatrix, FPGrad[scalingGradient] + GetIntegral[integrandPoint]];
	];
	Map[AddRow,Range[0,gridSize]];
];
Return[stabilityMatrix];
];


(* ::Subsubsection::Closed:: *)
(*EigenSys*)


(* ::Text:: *)
(*Returns all real pairs {eigenvalue, eigenvector} from eigensystem of mat, sorted by eigenvalues in ascending order*)


EigenSys[mat_] := Block[{eigenSys},
eigenSys = Sort[Transpose[Eigensystem[mat]]];
(*Return[Re[eigenSys]]];*)
eigenSys = Sort[eigenSys,Re[#1]<Re[#2]&];
Return[eigenSys]];
(*Return[Re[Sort[Select[eigenSys, Abs[Im[#[[1]]]]<0.00001 &]]]]]*)


(* ::Subsubsection::Closed:: *)
(*ExtractExponents*)


ExtractExponents[stabilityMatrix_]:=Block[{eigenvalues,criticalExponents=<| |>},
eigenvalues = EigenSys[stabilityMatrix][[All,1]];

criticalExponents[e1] = -eigenvalues[[1]];
criticalExponents[e2] = -eigenvalues[[2]];
criticalExponents[e3] = -eigenvalues[[3]];
criticalExponents[e4] = -eigenvalues[[4]];
criticalExponents[e5] = -eigenvalues[[5]];
Return[criticalExponents];
];


(* ::Subsubsection::Closed:: *)
(*GetCriticalExponents*)


(* ::Text:: *)
(*Finds a dictionary of critical exponents \[Nu], \[CapitalOmega], \[Eta], y_p for a fixed point.*)


GetCriticalExponents[integrandsListIsotropic_, integrandsListAnisotropic_, fixedPoint_, constants_] :=
Block[{j, stabilityMatrixIsotropic, stabilityMatrixAnisotropic,
	eigenvaluesIsotropic, eigenvaluesAnisotropic, eqlist, check, criticalExponents = Association[]},

stabilityMatrixIsotropic = StabilityMatrix[integrandsListIsotropic, fixedPoint, constants];
eigenvaluesIsotropic = EigenSys[stabilityMatrixIsotropic][[All,1]];

criticalExponents[\[Nu]] = -1/eigenvaluesIsotropic[[1]];
criticalExponents[\[CapitalOmega]] = eigenvaluesIsotropic[[2]];
criticalExponents[\[Eta]] = \[Eta] /. fixedPoint;
criticalExponents[e1] = -eigenvaluesIsotropic[[1]];
criticalExponents[e2] = -eigenvaluesIsotropic[[2]];
criticalExponents[e3] = -eigenvaluesIsotropic[[3]];
criticalExponents[e4] = -eigenvaluesIsotropic[[4]];
criticalExponents[e5] = -eigenvaluesIsotropic[[5]];

If[Length[integrandsListAnisotropic] > 0, 
	stabilityMatrixAnisotropic = StabilityMatrix[integrandsListAnisotropic, fixedPoint, constants];
	eigenvaluesAnisotropic = EigenSys[stabilityMatrixAnisotropic][[All,1]];
	
	criticalExponents[y1] = -eigenvaluesAnisotropic[[1]];
	criticalExponents[y2] = -eigenvaluesAnisotropic[[2]];
	criticalExponents[y3] = -eigenvaluesAnisotropic[[3]];
	criticalExponents[y4] = -eigenvaluesAnisotropic[[4]];
	criticalExponents[y5] = -eigenvaluesAnisotropic[[5]];

	For[j=1, j<=Length[eigenvaluesIsotropic], j++,
		(* y is the lowest eigenvalue that appears in anisotropic model but not in isotropic *)
		If[Abs[eigenvaluesIsotropic[[j]] - eigenvaluesAnisotropic[[j]]] < 10^-3,
			Continue[]];
		criticalExponents[y] = eigenvaluesAnisotropic[[j]];
		Break[];
	];
];

Return[{criticalExponents, stabilityMatrixIsotropic}];
];


(* ::Subsubsection::Closed:: *)
(*FindExponents*)


(* ::Text:: *)
(*Finds a dictionary of dictionaries of critical exponents \[Nu], \[CapitalOmega], \[Eta], y_p .*)
(*Result maps an exponent label (e.g. \[Nu], \[CapitalOmega], \[Eta], y) to a dictionary of its values with the same keys as fixedPointDict, except for \[Rho]Max which is removed.*)


FindExponents[integrandsListIsotropic_, integrandsListAnisotropic_, fixedPointDict_, keyLabels_, constants_] := 
Block[{keys, key, i, j, fixedPoint, criticalExponents, exponents, exp, params, cRep, equationsIsotropic,matrices},

keys = Keys[fixedPointDict];

exponents = Association[];
matrices = Association[];

For[i=1, i<=Length[keys], i++,
	key = keys[[i]];
	fixedPoint = fixedPointDict[key];
	
	params = constants;
	For[j=1, j<=Length[key], j++,
		params[keyLabels[[j]]] = key[[j]];
	];
	
	If[! SubsetQ[Keys[params], {"dim", "\[Rho]Max", "alpha"}],
		Print["Missing model parameters"];
		Quit[];];
	
	cRep = ConstRep[params["dim"], params["\[Rho]Max"], params["alpha"]];

	equationsIsotropic = integrandsListIsotropic;
	If[MemberQ[Keys[params],"N"], 
		equationsIsotropic = equationsIsotropic /. n -> params["N"]];
	
	{exponents[key[[1;;-2]]], matrices[key[[1;;-2]]]} = GetCriticalExponents[equationsIsotropic, integrandsListAnisotropic, fixedPoint, cRep];
	Print["Exponents for "<>ToString[keyLabels[[1;;-2]]]<>" = "<>ToString[key[[1;;-2]]]<>": "<>ToString[exponents[key[[1;;-2]]]]];
	];
criticalExponents = TransposeAssociation[exponents];

(*If["alpha" == keyLabels[[-2]],
	criticalExponents = FindPMSValues[criticalExponents];
];*)

Return[criticalExponents];
];



(* ::Subsubsection::Closed:: *)
(*FindPMSValues*)


FindPMSValues[exponents_] := Block[{i, pms, k, j, expDict, PMSDict, key, keys, oldKey, values, derivative},
expDict = KeySort[exponents,!OrderedQ[{#1,#2}]&];
PMSDict = Association[];
	
keys = Keys[expDict];
	
oldKey = keys[[1, 1;;-2]];
values = {};

For[j=1, j<=Length[expDict], j++,
	key = keys[[j, 1;;-2]];
	If[key != oldKey || j==Length[expDict],
		If[Length[values]<4,
			Print["PMS value not found at "<>ToString[oldKey]],
			
			derivative = Join[Join[{values[[2]]-values[[1]]}, 
						   (values[[3;;-1]] - values[[1;;-3]])/2],
						   {values[[-1]]-values[[-2]]}];
		
			pms = Position[Abs[derivative], Min[Abs[derivative]]][[1,1]];
			(*If[Length[oldKey] == 1,
				PMSDict[oldKey[[1]]] = values[[pms]],
				PMSDict[oldKey] = values[[pms]]];*)
				For[k=1, k<Length[derivative], k++,
					If[Sign[derivative[[k]]] != Sign[derivative[[k+1]]],
						If[Length[oldKey] == 1,
							PMSDict[oldKey[[1]]] = values[[k]],
							PMSDict[oldKey] = values[[k]];
							];
					
						Break[];
					];
				];
				If[k==Length[derivative], 
					Print["PMS value not found at "<>ToString[oldKey]]];
			];
		oldKey = key;
		values = {};
	];
	AppendTo[values, expDict[[j]]];
];

Return[PMSDict];
];

