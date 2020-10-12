(* ::Package:: *)

(* ::Section:: *)
(*Numerical Library*)


(* ::Subsection::Closed:: *)
(*Constants*)


(* ::Text:: *)
(*Size of \[Rho] grid*)


If[!NumberQ[gridSize],
	gridSize = 120];


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
	zNormalization = 0];


(* ::Subsection:: *)
(*Replacement rules*)


(* ::Subsubsection::Closed:: *)
(*Coefficients for finite difference derivatives at order O(h^4)*)


(* ::Text:: *)
(*Forward/BackwardCoefficients - derivative at the boundary*)
(*SemiForward/BackwardCoefficients - derivative one step from the boundary*)
(*CentralCoefficcients - derivatives away from boundary*)


Block[{n},
CentralCoefficients[n_] := Piecewise[{{{1/12,-2/3,0,2/3,-1/12},n==1}, {{-1/12,4/3,-5/2,4/3,-1/12},n==2}},{0}];
ForwardCoefficients[n_] := Piecewise[{{{-25/12,4,-3,4/3,-1/4},n==1}, {{15/4,-77/6,107/6,-13,61/12,-5/6},n==2}},{0}];
SemiForwardCoefficients[n_] := Piecewise[{{{-1/4,-5/6,3/2,-1/2,1/12},n==1}, {{5/6,-5/4,-1/3,7/6,-1/2,1/12},n==2}},{0}];
BackwardCoefficients[n_] := Piecewise[{{-ForwardCoefficients[n],OddQ[n]}}, ForwardCoefficients[n]];
SemiBackwardCoefficients[n_] := Piecewise[{{-SemiForwardCoefficients[n],OddQ[n]}}, SemiForwardCoefficients[n]]];


(* ::Input:: *)
(*(*n=5;*)
(*d=2;*)
(*rep = s[j_]\[Rule]j-1;*)
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
For[n=1, n<=2, n++, a=CentralCoefficients[n];
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
For[n=1, n<=2, n++,
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


(* ::Subsubsection::Closed:: *)
(*Shortcuts for replacement*)


(* ::Text:: *)
(*NumericDerivatives - replaces exact derivatives with finite differences*)
(*ONSymmetry - simplifies equations to ones for the O (N) model*)
(*LPA - simplifies equations to LPA approximation*)


ONSymmetry = {Subscript[w, i_][\[Rho]_] -> 0, \!\(\*SuperscriptBox[
SubscriptBox[\(w\), \(i_\)], 
TagBox[
RowBox[{"(", "m_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0};
NumericDerivatives := (# //. repderborder //.repder)&
LPA = {z[\[Rho]_]->1, \!\(\*SuperscriptBox[\(z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0, yy[\[Rho]_]->0, \!\(\*SuperscriptBox[\(yy\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0, \[Eta]->0};
ConstRep[dim_, \[Rho]max_, alpha_:defaultAlpha] := 
	{d->N[dim], vd->N[(2.^(-1-dim) \[Pi]^(-dim/2))/Gamma[dim/2]],
	eps->\[Rho]max/gridSize, \[Alpha]->alpha, z[zNormalization]->1};


(* ::Subsection::Closed:: *)
(*Gauss - Legendre integration*)


(* ::Text:: *)
(*Generates weights and evaluation points Gauss - Legendre integration for n \[Element] [1, MaxPointsLegendre].*)


MAXPOINTSLEGENDRE = 30;
rootsLegendre[n_] := N[Solve[LegendreP[n, x] == 0, x]] (* generates roots of the n-th Legendre poly. *) 
weightGL[n_, x_] := N[2/((1-x^2) D[LegendreP[n, x], x]^2)];
For[legendrePoints=1, legendrePoints <= MAXPOINTSLEGENDRE, legendrePoints++,
	roots = rootsLegendre[legendrePoints];
	weights = weightGL[legendrePoints, x] /. roots;
	GLTable[legendrePoints] = Table[{w[i] -> Chop[N[weights[[i]]]], N[roots[[i,1]]]},
		{i, 1, Length[rootsLegendre[legendrePoints]]}]
	];


(* ::Text:: *)
(*GLIntegral - carries out Gauss - Legendre integration of function on the interval [a, b] with n = pointsLegendre*)
(*GetIntegral - performs Gauss - Legendre integration of function split into parts specified in IntegralConfigurations*)


GLIntegral[f_,a_,b_,pointsLegendre_]:= (b-a)/2 Sum[(f[(a+b)/2 + 1/2 (b-a) x] w[i] /. GLTable[pointsLegendre][[i]]),{i, 1, pointsLegendre}]


GetIntegral[f_] := Sum[GLIntegral[f,IntegralConfigurations[[i,1]],IntegralConfigurations[[i,2]],IntegralConfigurations[[i,3]]],{i,1,Length[IntegralConfigurations]}]


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
		transposed[pk][sk] = association[sk][pk];
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
		nested[rest] = Join[nested[rest],newEntry],
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
	nested = KeySort[OneLevelNest[nested],#1<#2&];
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
  If[Abs[diffs[[pos]]] < 0.05, Return[sorted[pos]]];
  Return[Null];
  ]


(* ::Section:: *)
(*Renormalization Group*)


(* ::Subsection::Closed:: *)
(*Integration of the flow equations*)


(* ::Text:: *)
(*Initializes the list of flow equations and list of model parameters (without substituting the constants).*)


IntegrateEquations[integrandsList_]:=
Block[{i, j, equationsList, equation, parameters, function, scaling, integrand, yintegrand, yscaling, reppedEquation},
parameters = {}; equationsList = {};
For[i=1, i<=Length[integrandsList], i++,
	(* We identify the components of integrandsList *)
	function = integrandsList[[i, 1]];
	If[MatchQ[function, f_[0]],Continue[]];
	scaling = integrandsList[[i, 2]];
	integrand[q_] = integrandsList[[i, 3]] /. regulatorReplacement /. y -> q;

	(* Full flow equation consists of scaling term and integrated term *)
	equation = scaling + GetIntegral[integrand];

	(* Iteration over \[Rho] grid *)
	For[j=0, j<=gridSize, j++,
		If[j==0 && MatchQ[function, yy[\[Rho]/eps]],
			(* If the function is Y[\[Rho]] we have separate expression for the integrand at \[Rho]=0 *)
			yintegrand[q_] = integrandsList[[i+1, 3]] /. regulatorReplacement /. y -> q;
			yscaling = integrandsList[[i+1, 2]];
			reppedEquation = yscaling + GetIntegral[yintegrand],
	
			(* If not we just substitute \[Rho] *)
			reppedEquation = equation /. \[Rho] -> j eps
		];

		(* We append equation to the list of equations and function to the list of parameters*)
		AppendTo[equationsList, reppedEquation];
		AppendTo[parameters, function/. \[Rho] -> j eps];
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


AdjustNumberOfPoints[guess_, numberOfPoints_,zNormalization_]:= Block[{newGuess, scale, functions, guessLength, Interpolate,zDeviation},
functions = {v,z,yy};

guessLength = Length[Position[guess,v[i_]]]-1;

(* Function for linear interpolation onto new grid *)
scale = guessLength/numberOfPoints;
(* Construction of a new grid representation *)
interpolation[i_] := Interpolation[Table[{j, functions[[i]][j]},{j,0,guessLength}] /.guess /. z[k_]->1,InterpolationOrder->4];

newGuess = Table[functions[[i]][j] -> interpolation[i][j*scale]/(z[zNormalization]/.guess/.z[k_]->1), {i,1,Length[functions]}, {j,0,numberOfPoints}];
newGuess = Flatten[newGuess,1];

(* Remove zNormalization from FP specification *)
newGuess = DeleteCases[newGuess, z[zNormalization]->c_];

(* Add anomalous dimension \[Eta] *)
AppendTo[newGuess,\[Eta]->(\[Eta]/.guess)];

Return[newGuess];
]


(* ::Subsubsection::Closed:: *)
(*PlotFixedPoint*)


(* ::Text:: *)
(*Plots a ListPlot for a list of replacements*)


PlotFixedPoint[fixedPoint_] :=ListPlot[Transpose[Table[{v[i],z[i],yy[i]},{i,0,gridSize}]/.fixedPoint/.z[k_]->1], 
	PlotLegends->{"V","Z","Y"},PlotLabel->"\[Eta]="<>ToString[\[Eta]/.fixedPoint],AxesLabel->{"i","f[i]"}];


PlotFixedPoint[fixedPoint_,\[Rho]Max_] := Block[{tabelize,functions},
tabelize[f_] := Table[{i \[Rho]Max/gridSize, f[i]},{i,0,gridSize}]/.fixedPoint/.z[k_]->1;
functions = Map[tabelize,{v,z,yy}];
functions[[3,All,2]] = functions[[2,All,2]]- 2 functions[[3,All,1]]*functions[[3,All,2]];
ListPlot[functions, PlotLegends->{"V",Subscript["Z","\[Sigma]"],Subscript["Z","\[Pi]"]},AxesLabel->{"\[Rho]","f[\[Rho]]"},PlotLabel->"\[Eta]="<>ToString[\[Eta]/.fixedPoint]]];


(* ::Subsubsection::Closed:: *)
(*FindFixedPoint*)


(* ::Text:: *)
(*Finds a fixed point for provided equations*)


FindFixedPoint[eqlist_, guess_, constants_]:=
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
					MaxIterations->1000, PrecisionGoal->5, DampingFactor->2],
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


(* ::Subsubsection:: *)
(*\[Rho]MaxScan*)


(* ::Text:: *)
(*Returns an association with entries \[Rho]max -> FP. Starts at initial\[Rho]max and scans with additive increments in \[Rho] up to maxIterations.*)
(*On approaching inproper FPs (see IsWIlsonFisher function) it discards it, but does not stop the iterations. Iterations are stopped after 5 consecutive errors.*)


\[Rho]MaxScan[equations_, initialFixedPoint_, dim_, initial\[Rho]Max_, alpha_, \[Rho]Increment_, maxIterations_] := Block[
{d\[Rho], min, j, guess, newFixedPoint, fixedPoints, fail, \[Rho]Max, desiredMinimumPosition=Floor[50*gridSize/60], FPCheck=IsWilsonFisher, fpMsg, errorMsg},

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


(* ::Subsubsection:: *)
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
	\[Rho]Max=initial\[Rho]Max, step=0.05, steps=0, maxSteps=400, dim=d},
If[searchTricritical,
	Print["Searching for tricritical fixed points"],
	Print["Searching for critical fixed points"]];

While[True,
	steps++;
	scan = \[Rho]MaxScan[eqlist, newGuess, dim, \[Rho]Max, constants["alpha"], 0.0125*\[Rho]Max, 20];
	
	If[Length[scan]<2 || !IsFixedPoint[eqlist, scan[[1]], ConstRep[dim, scan[[2]], constants["alpha"]]], 
		Print["No FP found for d="<> ToString[dim]];
		If[step > 10^-5,
			step /= 2;
			dim += step;
			Continue[],
			Break[]];
		];
	newGuess = scan[[1]];
	\[Rho]Max = scan[[2]];
	fixedPointDict[{dim,\[Rho]Max}] = newGuess;
	dim -= step;
	If[steps>=maxSteps,Break[]];
];
	
Return[fixedPointDict];
];


(* ::Subsection:: *)
(*Stability Matrix*)


(* ::Subsubsection::Closed:: *)
(*Stability Matrix*)


(* ::Text:: *)
(*Returns the stability matrix of a system. M_i,j is a derivative of the i-th flow equation with respect to j-th free parameter of the model (evaluated at the fixed point). *)


StabilityMatrix[integrandsList_, fixedPoint_, constants_] := 
Block[{params, fixedPointRep, anisotropyOrder, i, j,
	FPGrad, FPSeries, ZIndex, \[Eta]Flow, \[Eta]Equation, \[Eta]Integrand, \[Eta]Solution, \[Eta]Replacement,
	function, integrand, scaling, scalingGradient, integrandPoint, stabilityMatrix, eqlist, check},
(* List of free parameters in a model *)
params = DeleteDuplicates[Flatten[Transpose[Table[integrandsList[[All, 1]]/. {\[Rho] -> i eps}, {i, 0, gridSize}]]]];
params = DeleteCases[params, z[zNormalization]];

(* Fixed Point replacement including only free parameters and cutoff for expansion in anisotropic field *)
fixedPointRep = Table[params[[i]]->(params[[i]] /. {t[k_] -> yy[k]} /. fixedPoint /. z[zNormalization]->1 /. {Subscript[w, j_][k_] -> 0} ),
		{i,1,Length[params]}];
		
(* Gradient operator in model parameter space *)
FPGrad[f_] := Grad[f, params]/. fixedPointRep;
(* Linear expansion around the fixed point *)
FPSeries[f_] := (f /. fixedPointRep) + (params -(params/.fixedPointRep)).FPGrad[f];

(* Anomalous dimension in terms of other parameters *)
\[Eta]Replacement = {\[Eta]->0};
ZIndex = Position[integrandsList[[All, 1]], z[\[Rho]/eps]];
If[Length[ZIndex]>0,
	ZIndex = ZIndex[[1,1]];
	\[Eta]Flow = integrandsList[[ZIndex]];
	\[Eta]Integrand[q_] = \[Eta]Flow[[3]] /. regulatorReplacement /. y -> q;
	\[Eta]Equation = \[Eta]Flow[[2]] + GetIntegral[\[Eta]Integrand];
	\[Eta]Equation = NumericDerivatives[\[Eta]Equation /. \[Rho] -> zNormalization eps];
	\[Eta]Solution = Solve[(\[Eta]Equation /. constants)==0, \[Eta]] // Flatten;
	\[Eta]Replacement = { \[Eta] -> FPSeries[\[Eta] /. \[Eta]Solution] };];
	
(* Derivatives of equations*)
stabilityMatrix = {};
For[i=1, i<=Length[integrandsList], i++,
	function = integrandsList[[i,1]];
	If[MatchQ[function, f_[0]],Continue[]];

	integrand[q_] = integrandsList[[i,3]] /. regulatorReplacement /. y -> q;
	scaling = integrandsList[[i,2]];

	For[j=0, j<=gridSize, j++,
		(* Skips one iteration at Z[normalization] *)
		If[function==z[\[Rho]/eps] && j==zNormalization, Continue[]];
		
		If[j==0 && MatchQ[function, yy[\[Rho]/eps]],
			(* Y[\[Rho]] has separate formula for integrand at \[Rho]=0 *)
			integrandPoint[q_] = FPGrad[NumericDerivatives[integrandsList[[-1,3]] 
				/. regulatorReplacement /. y -> q] /. \[Eta]Replacement /. constants];,
			(* Takes gradient of the flow equation *)
			integrandPoint[q_] = FPGrad[NumericDerivatives[integrand[q] /. \[Rho] -> j eps]
				/. \[Eta]Replacement /. constants];
		];

		scalingGradient = NumericDerivatives[scaling/. \[Rho]->j eps] /. \[Eta]Replacement /. constants;

		AppendTo[stabilityMatrix, FPGrad[scalingGradient] + GetIntegral[integrandPoint]];
		];
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
Return[Re[Sort[Select[eigenSys, Abs[Im[#[[1]]]]<0.00001 &]]]]]


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

If[Length[integrandsListAnisotropic] > 0, 
	stabilityMatrixAnisotropic = StabilityMatrix[integrandsListAnisotropic, fixedPoint, constants];
	eigenvaluesAnisotropic = EigenSys[stabilityMatrixAnisotropic][[All,1]];

	For[j=1, j<=Length[eigenvaluesIsotropic], j++,
		(* y is the lowest eigenvalue that appears in anisotropic model but not in isotropic *)
		If[Abs[eigenvaluesIsotropic[[j]] - eigenvaluesAnisotropic[[j]]] < 10^-4,
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
	If[Length[integrandsListAnisotropic] == 0, 
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

