(* ::Package:: *)

(* ::Title:: *)
(*NPRG Procedure for O(2) Model with Z4-symmetric Anisotropy*)


(* ::Subtitle:: *)
(*Functional Approach within Local Potential Approximation*)


(* ::Section::Closed:: *)
(*Initialization*)


arguments = $CommandLine;

If[MemberQ[$CommandLine,"-script"], debug = False, debug = True];
SetAttributes[debug,Protected];

If[debug,
	(* Debugging configuration *)
	n\[Phi]=0; 
	n\[Tau]=0; 
	ON=True;
	LPA=False;
	strict=False;
	anisotropySymmetry="ON";
	(*n=2;*)
	SetAttributes[n,{Protected,Constant}],
		
	(* Run configuration serialization *)
	notebookPath = arguments[[3]];
	slashes = StringPosition[notebookPath,"/"];
	If[Length[slashes]==0,
		projectPath = "./",
		projectPath = StringTake[notebookPath, slashes[[-1,1]]];
	];

	Get[projectPath<>"ModuleNaming.wl"];

	n\[Phi] = Null;
	n\[Tau] = Null;
	ON = False;
	strict = False;
	anisotropySymmetry = Null;
	Block[{i},
	For[i = 1, i<=Length[arguments], i++,
		If[arguments[[i]] == "Zp", 
			If[arguments[[i+1]] == "ON",
				ON = True;
				anisotropySymmetry = "ON";
				task = "ON";
				LPA=False;
				n\[Tau] = 0,
				ON = False;
				LPA = False;
				n = 2;
				anisotropySymmetry = ToExpression[arguments[[i+1]]];
				task = "Z"<>ToString[anisotropySymmetry]]];
		If[arguments[[i]] == "phi", n\[Phi] = ToExpression[arguments[[i+1]]]];
		If[arguments[[i]] == "tau", n\[Tau] = ToExpression[arguments[[i+1]]]];
		If[arguments[[i]] == "strict", strict=True];
		];
	];

	SetAttributes[{n,p,q,r,s},{Protected}];

	fail = False;
	If[anisotropySymmetry==Null,
		Print["Specify symmetry of the anistropic field - keyword: Zp"];
		fail = True;
	];
	If[n\[Phi] == Null,
		Print["Specify order of the field expansion - keyword: phi (0 for functional)"];
		fail = True;
	];
	If[n\[Tau] == Null,
		Print["Specify order of the expansion in powers of anistropy \[Tau] - keyword: tau"];
		fail = True;
	];
	If[fail, Quit[]];
];


(* ::Section:: *)
(*Invariants, Effective Action and Code Organizing Functions*)


(* ::Subsection::Closed:: *)
(*Protected symbols*)


Protect[{p,q,r,s,n,\[Phi],Z,U,V,W,T,dim}];


(* ::Subsection::Closed:: *)
(*Invariants*)


If[ON,
\[Rho][m_, x_] := Sum[m[i][x]^2, {i, 1, 3}]/2;
\[Rho][\[Phi]] := Sum[\[Phi][i]^2, {i, 1, 3}]/2,
\[Rho][m_, x_] := Sum[m[i][x]^2, {i, 1, 2}]/2;
\[Rho][\[Phi]] := Sum[\[Phi][i]^2, {i, 1, 2}]/2];


(* Choice of anisotropic invariant *)
If[anisotropySymmetry == 4 || anisotropySymmetry == "ON", 
	\[Tau][m_, x_] := m[1][x]^2*m[2][x]^2/2; \[Tau][\[Phi]] := \[Phi][1]^2*\[Phi][2]^2/2; 
	VariablesRange = \[Rho] > 0 && \[Rho]^2 > 2\[Tau]>0 && \[Phi][1] >= \[Phi][2] && \[Phi][2] >= 0 && \[Phi][1] >= 0 && \[Tau] > 0; 
	
	(* Solution for \[Phi][1] and \[Phi][2] expressed by \[Rho] and \[Tau] *)
	disentangled = Simplify[Solve[{\[Rho] == \[Rho][\[Phi]], \[Tau] == \[Tau][\[Phi]]}, {\[Phi][1], \[Phi][2]}, Reals][[8]], Assumptions -> VariablesRange];
    
    (* Expressions appearing in the flow equations that need to be expressed in terms of invariants *)
    repExpressions = {\[Phi][1]^2*\[Phi][2]^2, \[Phi][1]^4*\[Phi][2]^4, \[Phi][1]^2 + \[Phi][2]^2, \[Phi][1]^4 + \[Phi][2]^4, (\[Phi][1]^2 - \[Phi][2]^2)^2, \[Phi][1]^4 - \[Phi][2]^4, 
    (\[Phi][1]^2-3 \[Phi][2]^2), (-29 \[Phi][1]^2 \[Phi][2]^2+10 \[Phi][1]^4+3 \[Phi][2]^4)}; 
    ]; 

If[anisotropySymmetry == 6, 
	\[Tau][m_, x_] := (m[2][x]^3 - 3*m[1][x]^2*m[2][x])^2/4; \[Tau][\[Phi]] := (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2/4; 
    VariablesRange = \[Rho] > 0 && 2\[Rho]^3 > \[Tau] && \[Tau] > 0 && \[Phi][1] >= \[Phi][2] && \[Phi][2] >= 0 && \[Phi][1] >= 0; 
    
    (* Solution for \[Phi][1] and \[Phi][2] expressed by \[Rho] and \[Tau] *)
    disentangled = Simplify[Solve[{\[Rho] == (\[Phi][1]^2 + \[Phi][2]^2)/2, \[Tau] == (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2/4}, {\[Phi][1], \[Phi][2]}, Reals][[12]], Assumptions -> VariablesRange];
    
    (* Expressions appearing in the flow equations that need to be expressed in terms of invariants *)
    repExpressions = {(\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2, (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^4, \[Phi][1]^2 + \[Phi][2]^2, \[Phi][1]^6 - 15*\[Phi][1]^4*\[Phi][2]^2 + 15*\[Phi][1]^2*\[Phi][2]^4 - \[Phi][2]^6, 
      (3*\[Phi][1]^5*\[Phi][2] - 10*\[Phi][1]^3*\[Phi][2]^3 + 3*\[Phi][1]*\[Phi][2]^5)^2, (3*\[Phi][1]^7*\[Phi][2] - 7*\[Phi][1]^5*\[Phi][2]^3 - 7*\[Phi][1]^3*\[Phi][2]^5 + 3*\[Phi][1]*\[Phi][2]^7)^2, 
      3*\[Phi][1]^6 - 36*\[Phi][1]^4*\[Phi][2]^2 + 39*\[Phi][1]^2*\[Phi][2]^4 - 2*\[Phi][2]^6, (-3*\[Phi][1]^4*\[Phi][2] - 2*\[Phi][1]^2*\[Phi][2]^3 + \[Phi][2]^5)^2, (\[Phi][1]^3 - 3*\[Phi][1]*\[Phi][2]^2)^2, 
      (-3*\[Phi][1]^4 - 2*\[Phi][1]^2*\[Phi][2]^2 + \[Phi][2]^4)^2*\[Phi][2]^2, -3*\[Phi][1]^6 + 36*\[Phi][1]^4*\[Phi][2]^2 - 39*\[Phi][1]^2*\[Phi][2]^4 + 2*\[Phi][2]^6}; 
    ]; 
$Assumptions = VariablesRange;
varRep = Table[repExpressions[[i]] -> FullSimplify[repExpressions[[i]] /. disentangled, Assumptions -> VariablesRange], {i, 1, Length[repExpressions]}]; 


(* ::Subsection::Closed:: *)
(*Effective action (function of \[Rho] and \[Tau])*)


If[ON,
\[CapitalGamma]:= int[U[\[Rho][x]] + Zp[\[Rho][x]]/2  Sum[D[m[i][x], x]^2, {i, 1, 3}] + (Zs[\[Rho][x]]-Zp[\[Rho][x]])/(4 \[Rho][m,x]) (Sum[m[i][x] D[m[i][x], x], {i, 1, 3}])^2, x],
(*\[CapitalGamma]:= int[U[\[Rho][x],\[Tau][x]] + (Z[\[Rho][x]]-2 \[Rho][m, x] Y[\[Rho][x]])/2  Sum[D[m[i][x], x]^2, {i, 1, 2}] 
	+ Y[\[Rho][x]]/2 (Sum[m[i][x] D[m[i][x], x], {i, 1, 2}])^2 + T[\[Rho][x]]m[1][x]m[2][x] D[m[1][x],x]D[m[2][x],x], x];*)
(*\[CapitalGamma]:= int[U[\[Rho][x],\[Tau][x]] + (Z[\[Rho][x]])/2  Sum[D[m[i][x], x]^2, {i, 1, 2}],x ],*)
\[CapitalGamma]:= int[U[\[Rho][x],\[Tau][x]] + (Zs[\[Rho][x],\[Tau][x]])/2  Sum[D[m[i][x], x]^2, {i, 1, 2}] - (Zs[\[Rho][x],\[Tau][x]]-Zp[\[Rho][x],\[Tau][x]])/(4 \[Rho][m,x]) (m[1][x]^2 D[m[2][x], x]^2+m[2][x]^2 D[m[1][x], x]^2) + 2 T[\[Rho][x]] m[1][x]m[2][x] D[m[1][x],x]D[m[2][x],x], x];
];



(*rot = {m[1] \[Rule] (Cos[\[Phi]]m[1][#] - Sin[\[Phi]]m[2][#]&),m[2] \[Rule] (Cos[\[Phi]]m[2][#] + Sin[\[Phi]]m[1][#]&)};*)


(* ::Subsection::Closed:: *)
(*Regulator derivatives*)


regulatorRep = {R[\[Lambda],k_] -> Z[\[Lambda]] \[Lambda]^2 r[k/\[Lambda]^2]};
Block[{i, j}, For[i=0, i<=1, i++,
	For[j=KroneckerDelta[i, 0], j<=2, j++,
		AppendTo[regulatorRep, \!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"i", ",", "j"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda]_, k_] ->Evaluate[D[D[Z[\[Lambda]]\[Lambda]^2 r[k/\[Lambda]^2], {\[Lambda], i}], {k, j}]]];
];];]
dimensionlessRep = {d[q, q]->y^2, q->y, \[Lambda]->1, Z'[\[Lambda]]->-\[Eta], Z[\[Lambda]]->1};


(* ::Subsection::Closed:: *)
(*Homogeneous configuration*)


homogeneous := {\[Rho][x_] -> \[Rho], \!\(\*SuperscriptBox[\(\[Rho]\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0, \[Tau][x_] -> \[Tau], \!\(\*SuperscriptBox[\(\[Tau]\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0, m[i_][x_] -> \[Phi][i], \!\(\*SuperscriptBox[\(m[i_]\), 
TagBox[
RowBox[{"(", "j_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0}
minimal := {m[i_][x_] -> KroneckerDelta[i,1] Sqrt[2\[Rho]], \!\(\*SuperscriptBox[\(m[i_]\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0, \[Phi][i_] -> KroneckerDelta[i,1] Sqrt[2\[Rho]]}


(* ::Subsection::Closed:: *)
(*Rules for finding dot products of momenta*)


allPairs[q_List] := Flatten[Table[q[[i]]q[[j]] -> d[q[[i]], q[[j]]], {i, 1, Length[q]}, {j, 1, i}], 1];
CollectMomenta[q_List] := Simplify[(Expand[#] /. allPairs[q])]&


(* ::Subsection::Closed:: *)
(*Organizing functions *)


PotentialDers = {\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_], \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho],\[Tau]], T[\[Rho]_], \!\(\*SuperscriptBox[\(T\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_], Z[\[Rho]], Z[\[Rho],\[Tau]_], \!\(\*SuperscriptBox[\(Z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], \!\(\*SuperscriptBox[\(Z\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho],\[Tau]], Y[\[Rho]],Y[\[Rho],\[Tau]], \!\(\*SuperscriptBox[\(Y\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], \!\(\*SuperscriptBox[\(Y\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho],\[Tau]] };
DersCollect := Collect[#, PotentialDers, FullSimplify[#, Assumptions -> VariablesRange]& ]& 
Clean[q_List] := CollectMomenta[q][#/. homogeneous]&
CleanInvariants := DersCollect[#/. homogeneous ] //. varRep& 
Minimal := Simplify[# /. { \[Phi][i_] -> KroneckerDelta[i,1] Sqrt[2\[Rho]], \[Tau] -> 0}]&
ONSym := Simplify[# /. {\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_] /; j>0 -> 0, Subscript[f_, j_]->(0&), 
T[0]-> (Zs'[0]-Zp'[0])/4, T'[0] -> (Zs''[0]-Zp''[0])/8, T''[0]->1/12(-\!\(\*SuperscriptBox[\(Zp\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[0]+\!\(\*SuperscriptBox[\(Zs\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[0]), 
T -> ((Zs[#]-Zp[#])/(4#)&), \[Tau] -> 0}]&
Dimensionless := 4 y^(d-1) vd (# /. regulatorRep /. dimensionlessRep)&
(*Dimensionless := 2 vd y^(d/2 - 1) # /. regulatorRep /. dimensionlessRep&*)


(* ::Section:: *)
(*Definitions of Functionals*)


(* ::Subsection::Closed:: *)
(*Dot product d[p_,q_]*)


(*Homogeneity*)
d[a_ p_, b_ q_] := a b d[p, q];
d[a_ p_, q_] := a d[p, q];
d[ p_, b_ q_] := b d[p, q];
d[-q_, p_] :=-d[q, p];
d[q_, -p_] :=-d[q, p];
d[-q_, -p_] :=d[q, p];
(* Linearity *)
d[p1_+p2_, q_]:=d[p1, q]+d[p2, q];
d[q_, p1_+p2_]:=d[q, p1]+d[q, p2];
d[0, q_] :=0;
d[q_, 0] :=0;

SetAttributes[d,{Orderless,Protected}];


(* ::Subsection::Closed:: *)
(*Gradient grad[f_, p_,q_]*)


(* ::Text:: *)
(*Returns derivative of f with respect to p at p=0 with q as the integration variable in the flow equation*)


grad[f_, p_, q_] := Module[{der, simple, res, repped, a,b,s,r,i,j},
repped= f /.{d[q,q]->q^2, d[p,p] -> p^2, d[p,q]->p q/Sqrt[d], (a_. p + b_. q)^2 -> a^2 p^2 + 2 a b p q/Sqrt[d] + b^2 q^2};
der = SeriesCoefficient[repped,{p, 0, 1}];
Return[der];
];


(* ::Subsection::Closed:: *)
(*Laplacian lap[f_,p_,q_]*)


(* ::Text:: *)
(*Returns 1/(2 d)(Subscript[\[CapitalDelta], p]f )\!\(\*SubscriptBox[\("|"\), \("p=0"\)]\), q is the integration variable in the flow equation*)


lap[f_, p_, q_] := Module[{der, simple, res, repped, a,b,s,r,i,j},
repped= f /.{d[q,q]->q^2, d[p,p] -> p^2, d[p,q]->p q/Sqrt[d], (a_. p + b_. q)^2 -> a^2 p^2 + 2 a b p q/Sqrt[d] + b^2 q^2};
der = SeriesCoefficient[repped,{p,0, 2}];
Return[der]
];


(* ::Subsection::Closed:: *)
(*Integral int[f_,x_]*)


(* Integral over dx^d *)
int[f_+g_, x_] := int[f, x]+int[g, x] (* additivity *)
int[f_List, x_] := Table[int[f[[i]], x],{i, 1, Length[f]}] (* additivity *)
int[a_(f_+g_),x_] := int[a f, x]+int[a g, x] (* linearity *)
int[w_*c_, x_] := c int[w, x] /; FreeQ[{c}, x] (* homogeneity *)
int[Exp[a_ x_ + b_]f_., x_] := Exp[b] int[Exp[a x]f, x] /; FreeQ[{b}, x] (* homogeneity *)
int[Exp[I(a_ x_ + b_)]f_., x_] := Exp[I b] int[Exp[I a x]f, x] /; FreeQ[{b}, x] (* homogeneity *)
int[w_[a_*x_ +b_], x_] := 1/a int[w[x], x] /; FreeQ[{a, b}, x]  (* change of variables with translational invariance *)
int[w_[a_*x_], x_] := 1/a int[w[x], x] /; FreeQ[{a}, x] (* change of variables *)
int[w_[xnormalization_ + b_], x_] := int[w[x], x] /; FreeQ[{b}, x] (* translational invariance *)
int[f_ Tr[g_], x_] := Tr[int[f g, x]]  (* commutation with trace *)
int[Sum[w_, {i_, 1, n_}], x_] := Sum[int[w, x], {i, 1, n}] (* commutation with summation *)
int[0, x_] := 0 (* integral of 0 *)
int[f_.*Exp[I x_ k_], x_] := f DDelta[k] /; FreeQ[f,x]
\!\(\*SuperscriptBox[\(int\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[f_, x_]:=1;
\!\(\*SuperscriptBox[\(int\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[f_, x_]:=0;


(* ::Subsection::Closed:: *)
(*Functional Derivative in Momentum Space (fder[f, m_[q_]])*)


(* Definition for magnetic field *)
fder[m_[i_][x_], m_[j_][q_]] := Exp[I q x] KroneckerDelta[i-j]
fder[m_[i_]'[x_], m_[j_][q_]] := I q Exp[I q x]KroneckerDelta[i-j] 
fder[m_[i_]''[x_], m_[j_][q_]] := -d[q, q] Exp[I q x] KroneckerDelta[i-j]


(* Functional derivative in direct space *)
fder[int[f_, x_], m_[y_]] := int[fder[f, m[y]], x]  (* Commutative with integration *)
fder[Sum[f_, {i_, 1, j_}], m_[y_]] := Sum[fder[f, m[y]], {i, 1, j}]  (* Commutative with summation *)
fder[Tr[f_], m_[y_]] := Tr[fder[f, m[y]]]  (* Commutative with summation *)
fder[f_List, m_[y_]] := Table[fder[f[[i]], m[y]], {i, 1, Length[f]}]  (* Commutative with summation *)
fder[f_^i_, m_[y_]]:=i f^(i-1) fder[f, m[y]] (* Power derivative *)
fder[a_ int[f_, x_], m_[y_]] := a fder[ int[f, x], m[y]] (* Homogeneity *)
fder[f_ + g_, m_[y_]] := fder[f, m[y]] + fder[g, m[y]] (* Linearity *)
fder[f_ g_, m_[y_]] := fder[f, m[y]] g + f fder[g, m[y]] (* Leibnitz rule *)
fder[f_ .g_, m_[y_]] := fder[f, m[y]] .g +f .fder[g, m[y]] (* Leibnitz rule *)
fder[d[f_, g_], m_[y_]] := d[fder[f, m[y]] , g] +d[f, fder[g, m[y]]] (* Leibnitz rule *)
fder[f_[m[1][x_], m[2][x_]], m_[i_][y_]] :=  D[f[m[1][y], m[2][y]], m[i][y]] DDelta[y-x] (* Chain rule - function of m[1] and m[2] *)
fder[f_[\[Rho][x_]],  m_[i_][y_]] := fder[\[Rho][m, x], m[i][y]] D[f[\[Rho][x]], \[Rho][x]]  (* Chain rule - function of rho *)
fder[f_[\[Rho][x_], \[Tau][x_]], m_[i_][y_]] := fder[\[Rho][m, x], m[i][y]] D[f[\[Rho][x], \[Tau][x]], \[Rho][x]] + fder[\[Tau][m, x], m[i][y]] D[f[\[Rho][x], \[Tau][x]], \[Tau][x]]   (* Chain rule - function of rho and tau *)


(* Functions with functional derivative equal 0 *)
fder[0, m_[y_]] := 0 
fder[a_, m_[y_]] := 0 /; NumericQ[a] 
fder[int[x_, x_], m_[y_]] := 0  
fder[p, m_[i_][y_]] := 0
fder[q, m_[i_][y_]] := 0
fder[r, m_[i_][y_]] := 0
fder[s, m_[i_][y_]] := 0
fder[t, m_[i_][y_]] := 0
fder[R[\[Lambda],q_], m_[y_]] := 0
fder[DDelta[x_], m_[y_]] := 0
fder[DDelta'[x_], m_[y_]] := 0  
fder[DDelta''[x_], m_[y_]] := 0
DDelta[0] = 1;


(* ::Section::Closed:: *)
(*Correlation Functions*)


(* ::Subsubsection::Closed:: *)
(*Superscript[\[CapitalGamma],(2)] and propagator functions (without DiracDelta)*)


If[ON, fieldIndices = 3, fieldIndices = 2];


raw\[CapitalGamma]2[q_,p_]:=Table[fder[fder[\[CapitalGamma],m[i][q]],m[j][p]],{i,1,fieldIndices},{j,1,fieldIndices}];
\[CapitalGamma]2[q_,p_]:=FullSimplify[Clean[{p,q}][raw\[CapitalGamma]2[q,p]]];
rawG[q_,p_] := Inverse[raw\[CapitalGamma]2[q,p]+DiagonalMatrix[Table[DDelta[p+q] R[\[Lambda],q^2],{i,1,fieldIndices}]]];
G[q_,p_] := Simplify[Inverse[\[CapitalGamma]2[q,p]+DiagonalMatrix[Table[DDelta[p+q] R[\[Lambda],q^2],{i,1,fieldIndices}]]]];


(* ::Subsubsection::Closed:: *)
(*Superscript[\[CapitalGamma],(3)] function*)


raw\[CapitalGamma]3[p_,q_,r_] :=Module[{n},Table[fder[raw\[CapitalGamma]2[q,r],m[n][p]],{n,1,fieldIndices}]];
\[CapitalGamma]3[p_,q_,r_]:=FullSimplify[Clean[{p,q,r}][raw\[CapitalGamma]3[p,q,r]]];


(* ::Subsubsection::Closed:: *)
(*Superscript[\[CapitalGamma],(4)]function*)


raw\[CapitalGamma]4[p_,q_,r_,s_] :=Module[{n},Table[fder[raw\[CapitalGamma]3[q,r,s],m[n][p]],{n,1,fieldIndices}]];
\[CapitalGamma]4[p_,q_,r_,s_]:=FullSimplify[Clean[{p,q,r,s}][raw\[CapitalGamma]4[p,q,r,s]]];


(* ::Subsubsection::Closed:: *)
(*Superscript[\[CapitalGamma],(5)]function*)


raw\[CapitalGamma]5[p_,q_,r_,s_,t_] :=Module[{n},Table[fder[raw\[CapitalGamma]4[q,r,s,t],m[n][p]],{n,1,fieldIndices}]];
\[CapitalGamma]5[p_,q_,r_,s_,t_]:=FullSimplify[Clean[{p,q,r,s,t}][raw\[CapitalGamma]t[p,q,r,s,t]]];


(* ::Section:: *)
(*Flow equations*)


(* ::Subsection::Closed:: *)
(*Propagators in minimal configuration*)


DerRep = {\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"i", "-", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, 0] /; i>0 -> \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"i", "-", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], 
	\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho], 0] /; j>n\[Tau] -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, 0] -> \!\(\*SuperscriptBox[
SubscriptBox[\(W\), \(j\)], 
TagBox[
RowBox[{"(", "i", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], 
	\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, 0]/;MemberQ[{Zs,Zp},f] -> \!\(\*SuperscriptBox[\(f\), 
TagBox[
RowBox[{"(", "i", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], \!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, 0] /;MemberQ[{Zs,Zp},f]&&j>1 -> 0, 
	\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, 0] /;MemberQ[{Zs,Zp},f] -> \!\(\*SuperscriptBox[
SubscriptBox[\(f\), \(j\)], 
TagBox[
RowBox[{"(", "i", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], f_[\[Rho]_,0]/;MemberQ[{Zs,Zp},f] -> f[\[Rho]]};
anisotropyExpansion = { 
	\!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho], \[Tau]_] /; j>n\[Tau] -> 0, \[Tau]^i_ /; i>n\[Tau] -> 0, \[Phi][2]^i_ /; i/anisotropySymmetry > n\[Tau] ->0, Subscript[W, j_]/;j>n\[Tau]->(0&),
	Subscript[Zs, j_]/;j>1->(0&),Subscript[Zp, j_]/;j>1->(0&)
	};
If[n\[Tau]>0, AppendTo[anisotropyExpansion, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "n\[Tau]"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, \[Tau]_] -> \!\(\*SuperscriptBox[
SubscriptBox[\(W\), \(n\[Tau]\)], 
TagBox[
RowBox[{"(", "i", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]]];
	AppendTo[anisotropyExpansion, \!\(\*SuperscriptBox[\(f_\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "n\[Tau]"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho], \[Tau]_] -> \!\(\*SuperscriptBox[
SubscriptBox[\(f\), \(n\[Tau]\)], 
TagBox[
RowBox[{"(", "i", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]]]];


Propagator[p_, q_] := Minimal[DDelta[p+q]^2 G[p, q]]/.DerRep
Propq = Propagator[p,q] /. {p -> -q} /. {DDelta[0] -> 1};
If[ON, Propq = ONSym[Propq]];


(* ::Text:: *)
(*Gs is the longitudinal (sigma) part of the propagator, Gp is the transverse part (pi)*)


Gs = Propq[[1, 1]];
Gp = Propq[[2, 2]];


(* ::Text:: *)
(*S is the operator G (Subscript["\[PartialD]", \[CapitalLambda]]R) G, it consists of Ss and Sp parts *)


Ss =  D[R[\[Lambda],q^2], \[Lambda]] Gs^2;
Sp =  D[R[\[Lambda],q^2], \[Lambda]] Gp^2;


(* ::Text:: *)
(*Replacement rules and functions for collecting powers of propagator*)


props = {Propq[[1,1]], Propq[[2,2]], 1/gs, 1/gp}; 
propRep = {1/props[[1]] -> gs, Expand[1/props[[1]]] -> gs, Simplify[1/props[[2]]] -> gp, Expand[1/props[[2]]] -> gp,
	Expand[\[Rho]/props[[2]]] -> \[Rho] gp, Simplify[(1/props[[1]] + 1/props[[2]])/2] -> (gp + gs)/2, Expand[(1/props[[1]] + 1/props[[2]])/2] -> (gp + gs)/2,
	FullSimplify[1/Propq[[1,1]]+1/Propq[[2,2]]]->gp+gs};
propRep = propRep/.anisotropyExpansion;
propRep= Join[propRep,(propRep/.d[q_,q_]->q^2)];
propRep = Join[propRep, propRep /. d[q_,p_] -> q p];
propReverse = {gs -> 1/props[[1]], gp -> 1/props[[2]], g-> (1/props[[1]] /.\[Rho]->0)}/.d[q,q]->q^2;

propReverseDimensionless = propReverse /. regulatorRep /. dimensionlessRep;

ReplaceProps := Collect[Collect[Collect[#//.propRep, props, Simplify] //. propRep, props, FullSimplify] //. propRep,
	{1/gs,1/gp},FullSimplify]&;


order = 3;
Map[({Subscript[gs,#],Subscript[gp,#]} = SeriesCoefficient[{gs,gp}/.propReverse,{\[Rho],0,#}])&,Range[order]];

EvAt\[Rho]0[f_] := (Collect[f /. {\[Tau] -> 0, 
 gs^i_. -> Series[gs[\[Rho]]^i,{\[Rho],0,order}]/.{gs[0]->g, Derivative[k_][gs][0]-> k! Subscript[gs,k]},
 gp^i_. -> Series[gp[\[Rho]]^i,{\[Rho],0,order}]/.{gp[0]->g, Derivative[k_][gp][0]-> k! Subscript[gp,k]}, 
    Subscript[Zs, 1]->(Subscript[Zs, 1][#]&),
	(*Subscript[Zs, 1]\[Rule]Subscript[Zp, 1][0]+# Subscript[Zs, 1]'[0]+#^2/2 Subscript[Zs, 1]''[0]&), *)
	Zs->(Zp[0]+# Zs'[0]+#^2/2 Zs''[0] + #^3/6 Zs'''[0]&)},{g}, Normal[Series[#,{\[Rho],0,0}]]&])


(* ::Subsection::Closed:: *)
(*Wetterich equation (Exact)*)


If[ON, 
	
	wetterich = 1/2 D[R[\[Lambda],q^2],\[Lambda]] (Gs + (n-1)Gp),
	
	(*We calculate trace of the modified propagator*)
	prop = G[q,p];
	traced = 1/2 Together[prop[[1]][[1]] + prop[[2]][[2]]] /.p->-q;

	(*We separately simplify the numerator and denominator of the trace.*)
	GNum = Numerator[traced];
	GDen = Denominator[traced];
	
	
	GNumCollected = FullSimplify[DersCollect[GNum] //. varRep] //. varRep//.disentangled;
	GDenCollected = FullSimplify[DersCollect[GDen] //. varRep] //. varRep//.disentangled;
	(*We recover full Wetterich equation at O(\[PartialD]^2)*)
	wetterich = D[R[\[Lambda],q^2],\[Lambda]]  FullSimplify[(GNumCollected/GDenCollected)/.anisotropyExpansion]/.anisotropyExpansion //. varRep;
];

wetterich = wetterich /. p->-q /. DDelta[0] -> 1;


(* ::Subsection::Closed:: *)
(*Flow equations for Local Potential Expanded at \[Tau]=0*)


flowEquations = {};


If[n\[Phi]==0,
VEquation =Collect[D[FullSimplify[wetterich /.\[Tau]->0], \[Rho]], GDenCollected/.\[Tau]->0, Simplify]/.DerRep;
AppendTo[flowEquations, {V[\[Rho]], (-2+\[Eta])V[\[Rho]]+(-2+d+\[Eta]) \[Rho] V'[\[Rho]], ReplaceProps[VEquation]}];

If[!ON,
Block[{kk},
For[kk=1, kk<=n\[Tau], kk++,
	WiEquation = kk! ReplaceProps[SeriesCoefficient[wetterich, {\[Tau], 0, kk}]/.DerRep];
	AppendTo[flowEquations,{Subscript[W, kk][\[Rho]],
		(anisotropySymmetry (d-2+\[Eta]) kk/2 - d) Subscript[W, kk][\[Rho]] + (-2+d+\[Eta]) \[Rho] Subscript[W, kk]'[\[Rho]],
		ReplaceProps[Collect[WiEquation,{gs,gp},FullSimplify] /.DerRep]}];
	]];
];
];



(* ::Subsection::Closed:: *)
(*Flow equations for Local Potential Expanded at \[Tau]=0 and \[Rho]=\[Kappa] at LPA'*)


If[n\[Phi]!=0,
LPArep = {Z[\[Rho]_] -> Z, \!\(\*SuperscriptBox[\(Z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0, Y[\[Rho]_] -> 0, \!\(\*SuperscriptBox[\(Y\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0};
If[n\[Tau]==0,
DerRep = Join[{\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] /; i<=1 -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]/; 2i>n\[Phi] -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->Subscript[U, i]}, LPArep],
DerRep = Join[{\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_] /; i<=1 -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_]/; 2i + anisotropySymmetry j>n\[Phi] || j > n\[Tau] -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i_", ",", "j_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_]->Subscript[U, i, j]}, LPArep];
];

wetterichLPA = wetterich /. LPArep;

wetterichTauExpansion = Table[ii!SeriesCoefficient[wetterichLPA, {\[Tau], 0, ii}], {ii, 0, n\[Tau]}];
wetterichFullExpansion = Dimensionless[Table[jj! SeriesCoefficient[wetterichTauExpansion[[ii+1]],{\[Rho],\[Kappa],jj}]/.DerRep ,{ii,0,n\[Tau]},{jj,0,(n\[Phi]-anisotropySymmetry ii)/2}]];

If[n\[Tau]==0,
dt\[Kappa] =- wetterichFullExpansion[[1]][[2]]/Subscript[U, 2],
dt\[Kappa] =- wetterichFullExpansion[[1]][[2]]/Subscript[U, 2,0]];

equationsExpanded = Table[If[i==0&&j==0, {0,0,0}, If[i==1&&j==0, {\[Kappa],2-d-\[Eta],dt\[Kappa]}, 
	{Subscript[u, i,j],Expand[(i +anisotropySymmetry/2 j)(d-2+\[Eta])-d], dt\[Kappa] Subscript[u, i+1,j] + wetterichFullExpansion[[j+1]][[i+1]]}]],
	{j,0,n\[Tau]},{i,0,(n\[Phi]-anisotropySymmetry j)/2}];
	
flowEquations = Drop[Flatten[equationsExpanded,1],1] /. {Subscript[u, i_,j_] /; 2i + anisotropySymmetry j>n\[Phi]->0};
];


(* ::Subsection:: *)
(*Equations for flow of fluctuation suppressors*)


(* ::Text:: *)
(*We assume \[PartialD]\[Tau] Z =0, so equations are evaluated at configuration \[Tau]=0, \[Phi][2]=0.*)


(* ::Subsubsection::Closed:: *)
(*Vertices in general*)


If[!ON,
	\[CapitalGamma]3 = Simplify[Simplify[Clean[{p,q,r}][raw\[CapitalGamma]3[p,q,r]]]/. anisotropyExpansion/. {r->-p-q}];
	\[CapitalGamma]4 = Simplify[Simplify[Clean[{p,q,r,s}][raw\[CapitalGamma]4[p,q,r,s]]]/.anisotropyExpansion/. {q->-p,r->q,s->-q}];
	Gq = G[p,q] /. {p->-q}/.anisotropyExpansion;]


(* ::Subsubsection::Closed:: *)
(*Flow equations from \[CapitalGamma](2) - for Z4 anisotropy*)


(* ::Text:: *)
(*General function for extracting flow equations for \[CapitalGamma](2) functions*)


If[!ON, 
	\[CapitalGamma]4exp = Normal[Series[\[CapitalGamma]4//.varRep/.disentangled,{\[Tau],0,Min[n\[Tau],1]}]]/.DerRep;
	\[CapitalGamma]3exp = Normal[Series[\[CapitalGamma]3//.varRep/.disentangled,{\[Tau],0,Min[n\[Tau],1]}]]/.DerRep;
	Gqexp = Normal[Series[Gq//.varRep/.disentangled,{\[Tau],0,Min[n\[Tau],1]}]]/.DerRep;
	
	propExpand[i_] := (Minimal[1/Gq[[i,i]]] + p (grad[Minimal[1/Gq[[i,i]]]/.{q->p+q},p,q]) + p^2 FullSimplify[lap[Minimal[1/Gq[[i,i]]]/.{q->p+q},p,q]]) /.DerRep/.propRep;
	propSeries = {gs[p+q] ->  propExpand[1], gp[p+q] -> propExpand[2]};

	G[r_] = Gqexp /. propRep /. q->r /. {gp->gp[r], gs->gs[r]};
	G[q] = G[q]/.{gp[q]->gp,gs[q]->gs};
	];


If[!ON,
	If[BooleanQ[strict] && strict,
	Apply[(\[CapitalGamma]3sq[#1,#2,#3,#4,#5,#6] = FullSimplify[ExpandAll[Series[\[CapitalGamma]3exp[[#1,#2,#3]]\[CapitalGamma]3exp[[#4,#5,#6]],{\[Tau],0,1}]] /.{d[x_,y_]^2->0,d[x_,y_] d[z_,v_]->0}])&,Tuples[Table[Range[2],{i,1,6}]],1],
	Apply[(\[CapitalGamma]3sq[#1,#2,#3,#4,#5,#6] = FullSimplify[Series[\[CapitalGamma]3exp[[#1,#2,#3]]\[CapitalGamma]3exp[[#4,#5,#6]],{\[Tau],0,1}]])&,Tuples[Table[Range[2],{i,1,6}]],1];
	];
];


G2flow[i0_,j0_, order_] := Block[{tadpole,Eval,bubble,res},
tadpole = SeriesCoefficient[-1/2 Tr[\[CapitalGamma]4exp[[i0,j0]].G[q].G[q]],{\[Tau],0,order/2}];
(*bubble = SeriesCoefficient[ Tr[\[CapitalGamma]3exp[[i0]].G[q+p].Transpose[\[CapitalGamma]3exp[[j0]]].G[q].G[q]],{\[Tau],0,order/2}];*)
bubble = SeriesCoefficient[ Sum[\[CapitalGamma]3sq[i0,i1,i2,j0,i4,i3] G[q+p][[i2,i3]]G[q][[i4,i5]]G[q][[i5,i1]],{i1,1,2},{i2,1,2},{i3,1,2},{i4,1,2},{i5,1,2}],{\[Tau],0,order/2}];
Eval[diag_] := SeriesCoefficient[diag /. propSeries/. {d[r_,r_]->r^2,d[p,q] -> p q/Sqrt[d]},{p,0,2}];

res = Eval[tadpole]+Eval[bubble];
If[order==1, res = res/Sqrt[8]];

Return[D[R[\[Lambda],q^2],\[Lambda]]Collect[res,{1/gs,1/gp},Simplify]];
];


If[!ON,
	flowZ1 = G2flow[1,1,0];
	AppendTo[flowEquations,{Zs[\[Rho]], \[Eta] Zs[\[Rho]]+(-2+d+\[Eta]) \[Rho] Zs'[\[Rho]], flowZ1}];
	
	If[!LPA, 
		flowZ2 = G2flow[2,2,0];
		AppendTo[flowEquations, {Zp[\[Rho]], \[Eta] Zp[\[Rho]]+(-2+d+\[Eta]) \[Rho] Zp'[\[Rho]], flowZ2}];];

	If[!LPA && n\[Tau]>0,
		yterm = (flowZ1-flowZ2)/(2\[Rho]^2);
		flowZ11 = G2flow[1,1,2] + yterm;
		AppendTo[flowEquations, {Subscript[Zs, 1][\[Rho]], (2d-4+3\[Eta]) Subscript[Zs, 1][\[Rho]] + (-2+d+\[Eta]) \[Rho] Subscript[Zs, 1]'[\[Rho]], flowZ11}];
		
		flowZ21 = G2flow[2,2,2] - yterm;
		AppendTo[flowEquations, {Subscript[Zp, 1][\[Rho]], (2d-4+3\[Eta]) Subscript[Zp, 1][\[Rho]] + (-2+d+\[Eta]) \[Rho] Subscript[Zp, 1]'[\[Rho]], flowZ21}];
		
		flowT = G2flow[1,2,1];
		AppendTo[flowEquations, {T[\[Rho]], (d-2+2\[Eta]) T[\[Rho]]+(-2+d+\[Eta]) \[Rho] T'[\[Rho]], flowT}];
	];
];


(*on2 = ONSym[flowZ11];
on3 = ONSym[flowZ21];
Collect[(on2 * gs^6*gp^6)/(R^(1,0))[\[Lambda],q^2]/. gs \[Rule] gp+q^2 (Zs[\[Rho]]-Zp[\[Rho]])+2 \[Rho] V'[\[Rho]],gp,Simplify]
Collect[(on3 * gs^6*gp^6)/(R^(1,0))[\[Lambda],q^2]/. gs \[Rule] gp+q^2 (Zs[\[Rho]]-Zp[\[Rho]])+2 \[Rho] V'[\[Rho]],gp,Simplify]*)


(* ::Subsubsection::Closed:: *)
(*Vertices in minimal configuration - for pure O(N) models only*)


If[ON,
	M\[CapitalGamma]3 = Simplify[Minimal[Clean[{p,q,r}][raw\[CapitalGamma]3[p,q,r]]]]/. r->-p-q;
	M\[CapitalGamma]3 = DersCollect[ONSym[M\[CapitalGamma]3/.DerRep]];
	
	M\[CapitalGamma]4 = Simplify[Minimal[Clean[{p,q,r,s}][raw\[CapitalGamma]4[p,q,r,s]]]]/. {q->-p, r->q, s->-q};
	M\[CapitalGamma]4 = DersCollect[ONSym[M\[CapitalGamma]4/.DerRep]];
];


(*If[!ON && !LPA,
	M\[CapitalGamma]5 = Simplify[Minimal[Clean[{p,q,r,s,t}][raw\[CapitalGamma]5[p,q,r,s,t]]]];
	M\[CapitalGamma]5 = DersCollect[M\[CapitalGamma]5];]*)


(* ::Text:: *)
(*For flows of both Z and Y we have two contributing diagrams: tadpole with one \[CapitalGamma]4 vertex and bubble with two \[CapitalGamma]3 vertices.*)


(* ::Subsubsection:: *)
(*Flow of Z[\[Rho]] - for pure O(N) models only*)


(* ::Text:: *)
(*The tadpole diagram*)


If[ON,
tadpole = Simplify[-1/2(gs^-2 M\[CapitalGamma]4[[1,1,1,1]] + (n-1) gp^-2 M\[CapitalGamma]4[[1,1,2,2]] )];
tadpoleDifferentiated = ReplaceProps[lap[tadpole, p, q]];

bubble = Simplify[M\[CapitalGamma]3[[1,1,1]]^2] gs^-2 (Gs /. q -> q+p) + (n-1) Simplify[M\[CapitalGamma]3[[1,2,2]]^2] gp^-2 (Gp /. q -> q+p);
bubbleDifferentiated = ReplaceProps[lap[bubble, p, q]];

Zflow = \!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda],q^2] Collect[(tadpoleDifferentiated + bubbleDifferentiated), {1/gs, 1/gp}, FullSimplify]/.propReverse;

If[n\[Phi]==0,
	ZEquation = {Zs[\[Rho]], \[Eta] Zs[\[Rho]]+(-2+d+\[Eta]) \[Rho] Zs'[\[Rho]], Zflow/.DerRep},
	ZEquation = {Z, \[Eta], Zflow /. DerRep /. \[Rho] -> \[Kappa]};
];
AppendTo[flowEquations, ZEquation];
];


(* ::Subsubsection:: *)
(*Flow equations for Y[\[Rho]] - for pure O(N) models only*)


If[n\[Phi]==0 && ON,
	Ytadpole = Collect[-1/2(M\[CapitalGamma]4[[2,2,1,1]] gs^-2 + M\[CapitalGamma]4[[2,2,2,2]] gp^-2 + (n-2) M\[CapitalGamma]4[[2,2,3,3]] gp^-2) /. propRep,{gs,gp}, Simplify];
	YtadpoleDifferentiated = ReplaceProps[lap[Ytadpole, p, q]];

	Ybubble = Simplify[ M\[CapitalGamma]3[[2,2,1]]^2] (gp^-2 (Gs /. q -> q+p)) + Simplify[M\[CapitalGamma]3[[2,1,2]]^2] (gs^-2(Gp /. q -> q+p));
	YbubbleDifferentiated = ReplaceProps[lap[Ybubble, p, q]];

	Yflow = YtadpoleDifferentiated + YbubbleDifferentiated;
	YflowSimplified = \!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda],q^2] Collect[Yflow, {1/gs, 1/gp}, Simplify]/. propReverse ;
	
	AppendTo[flowEquations, {Zp[\[Rho]], \[Eta] Zp[\[Rho]] + (-2+d+\[Eta]) \[Rho] Zp'[\[Rho]], YflowSimplified /. DerRep}];
];


(* ::Subsection::Closed:: *)
(*Apply chosen expansions*)


If[n\[Phi]==0,
	flowEquations = flowEquations /. {Subscript[W, i_]/;i>n\[Tau] -> (0&), 
		Subscript[Zs, i_]/;i>Min[n\[Tau],1] -> (0&),Subscript[Zp, i_]/;i>Min[n\[Tau],1] -> (0&)};
	];
If[n\[Tau]==0,
	flowEquations = flowEquations /. T->((Zs[#]-Zp[#])/(4#)&);];
If[LPA,
	flowEquations[[All,2;;3]] = flowEquations[[All,2;;3]] /. {Z->(1&)};
	flowEquations = flowEquations /. {T->(T&),Y->(Y&),Z->(Z&)}];


If[!LPA,
	Add\[Rho]0Eq[f_] := Block[{\[Rho]0eval, newf = f},
		newf[[3]] = Dimensionless[f[[3]]/.propReverse];
		If[MemberQ[{V[\[Rho]],Subscript[W, 1][\[Rho]],Subscript[W, 2][\[Rho]],Zs[\[Rho]], Subscript[Zs,1][\[Rho]],T[\[Rho]]},f[[1]]], 
			\[Rho]0eval = Dimensionless[EvAt\[Rho]0[f[[3]]]/.propReverse];
			AppendTo[newf, \[Rho]0eval];];
		Return[newf]];
	dimlessFlowEquations = Map[Add\[Rho]0Eq,flowEquations];];


(* ::Section:: *)
(*Cache the equations*)


(* ::Subsection:: *)
(*Consistency check for functional equations at O(\[Tau])*)


(* ::Text:: *)
(*Replacements for O(N) models without anisotropy*)


If[ON,
	repJ ={YY->Y, W[\[Rho]_]->0, W'[\[Rho]_]->0, W''[\[Rho]_]->0};
	repA0 = {y^p_. -> y^(p/2)};
	repA = {\!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"i", "+", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]],Subscript[W, i_][\[Rho]_]->0,\!\(\*SuperscriptBox[
SubscriptBox[\(W\), \(i_\)], 
TagBox[
RowBox[{"(", "m_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0,V[\[Rho]_]->U'[\[Rho]],Zs->(Z[#]+2 # Y[#]&),Zp->(Z[#]&),n->2, r[y]->y r[y], r'[y]->r[y] + y r'[y], r''[y]->2 r'[y]+y r''[y]};
	repProp = {y r[y]+y Z[\[Rho]_]+Derivative[1][U][\[Rho]_]->gp, y (r[y]+ Z[\[Rho]_])+Derivative[1][U][\[Rho]_]->gp,y r[y]+y (2 \[Rho] Y[\[Rho]_]+Z[\[Rho]_])+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs,y (r[y]+2 \[Rho] Y[\[Rho]_]+Z[\[Rho]_])+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs};
	];


(* ::Text:: *)
(*Replacements for O(2) models with anisotropy*)


If[!ON,
	repJ ={ YY->(Y[#]&)};
	repA0 = {y^p_. -> y^(p/2)};
	repA = {V->(U'[#]&),Subscript[W, i_]/;i==1->(W[#]&),Subscript[W, i_]/;i>1->(0&), 
		Subscript[Zp, i_]->(0&),Subscript[Zs, i_]->(0&), T->(Y[#]/2&),
		(2 r[y]-\[Eta] r[y]-2 y Derivative[1][r][y]) -> (-y \[Eta] r[y]-2 y^2 Derivative[1][r][y]), r->(# r[#]&),Zp->(Z[#]&), Zs->(Z[#] +2 # Y[#]&), dim->d};
	repProp = {y r[y]+y Z[\[Rho]]+Derivative[1][U][\[Rho]_]+ 2 \[Rho] W[\[Rho]]->gp,y (r[y]+ Z[\[Rho]]) +Derivative[1][U][\[Rho]_]+ 2 \[Rho] W[\[Rho]]->gp,
		y r[y]+y Z[\[Rho]]+ y 2 \[Rho] Y[\[Rho]]+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs, y r[y]+y (Z[\[Rho]]+  2 \[Rho] Y[\[Rho]])+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs, (y (r[y]+2 \[Rho] Y[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] (U^\[Prime]\[Prime])[\[Rho]])->gs};
	]


(* ::Text:: *)
(*Equations copied from Pawe\[LSlash]'s notebook*)


flowZJ = 2 vd y^(-1+d/2) (-y \[Eta] r[y]-2 y^2 Derivative[1][r][y]) ((4 \[Rho] Derivative[1][Z][\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])-(2 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2)/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+1/d y (2 \[Rho] Derivative[1][Z][\[Rho]]^2-(8 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) Derivative[1][Z][\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+2 \[Rho] ((4 (Z[\[Rho]]^2+2 Z[\[Rho]] (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2))/((y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]]) (y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2)+(4 \[Rho] YY[\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])-(2 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2)/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+1/d y (2 \[Rho] (-2 YY[\[Rho]]+Derivative[1][Z][\[Rho]])^2+(8 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (-2 YY[\[Rho]]+Derivative[1][Z][\[Rho]]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+2 \[Rho] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2 ((4 ((2 \[Rho] YY[\[Rho]]+Z[\[Rho]])^2+2 (2 \[Rho] YY[\[Rho]]+Z[\[Rho]]) (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]]))))/((y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2 (y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]]))+1/2 (-((2 YY[\[Rho]]+Derivative[1][Z][\[Rho]])/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2)-(Derivative[1][Z][\[Rho]]+2 \[Rho] Z''[\[Rho]])/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2));
flowYJ = 1/(2 \[Rho]) (-2 vd y^(-1+d/2) (-y \[Eta] r[y]-2 y^2 Derivative[1][r][y]) ((4 \[Rho] Derivative[1][Z][\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])-(2 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2)/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+1/d y (2 \[Rho] Derivative[1][Z][\[Rho]]^2-(8 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) Derivative[1][Z][\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+2 \[Rho] ((4 (Z[\[Rho]]^2+2 Z[\[Rho]] (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2))/((y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]]) (y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2)+(4 \[Rho] YY[\[Rho]] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])-(2 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2)/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+1/d y (2 \[Rho] (-2 YY[\[Rho]]+Derivative[1][Z][\[Rho]])^2+(8 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (-2 YY[\[Rho]]+Derivative[1][Z][\[Rho]]) (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+2 \[Rho] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2 ((4 ((2 \[Rho] YY[\[Rho]]+Z[\[Rho]])^2+2 (2 \[Rho] YY[\[Rho]]+Z[\[Rho]]) (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]]))))/((y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2 (y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]]))+1/2 (-((2 YY[\[Rho]]+Derivative[1][Z][\[Rho]])/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2)-(Derivative[1][Z][\[Rho]]+2 \[Rho] Z''[\[Rho]])/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2))+2 vd y^(-1+d/2) (-y \[Eta] r[y]-2 y^2 Derivative[1][r][y]) (1/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^3(4 \[Rho] YY[\[Rho]] (2 W[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]])-(2 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) (2 W[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]])^2)/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+1/d y (2 \[Rho] Derivative[1][Z][\[Rho]]^2-(8 \[Rho] (r[y]+Z[\[Rho]]+y Derivative[1][r][y]) Derivative[1][Z][\[Rho]] (2 W[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])+2 \[Rho] ((4 (Z[\[Rho]]^2+2 Z[\[Rho]] (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])) (2 W[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]])^2))+1/2 (-((2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])/(y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2)-(2 YY[\[Rho]]+10 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]]+4 \[Rho]^2 Y''[\[Rho]]+2 \[Rho] Z''[\[Rho]])/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2)+(4 \[Rho] (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]]) (y (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])+3 U''[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]])-(2 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (y (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])+3 U''[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]])^2)/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+1/d y (2 \[Rho] (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])^2-(8 \[Rho] (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]]+y Derivative[1][r][y]) (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]]) (y (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])+3 U''[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])+2 \[Rho] ((4 ((2 \[Rho] YY[\[Rho]]+Z[\[Rho]])^2+2 (2 \[Rho] YY[\[Rho]]+Z[\[Rho]]) (r[y]+y Derivative[1][r][y])+(r[y]+y Derivative[1][r][y])^2))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2-(2 (2 Derivative[1][r][y]+y r''[y]))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])) (y (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])+3 U''[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]])^2))/(y r[y]+y (2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^3));
flowY0J=1/(y r[y]+y Z[0]+Derivative[1][U][0])^3 2 vd y^(d/2) (\[Eta] r[y]+2 y Derivative[1][r][y]) (2 W[0] YY[0]+3 y r[y] Derivative[1][YY][0]+3 y Z[0] Derivative[1][YY][0]+3 Derivative[1][U][0] Derivative[1][YY][0]+4 W[0] Derivative[1][Z][0]-2 YY[0] U''[0]+2 Derivative[1][Z][0] U''[0]-(2 (r[y]+Z[0]+y Derivative[1][r][y]) (2 W[0]+y YY[0]+U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])+((r[y]+Z[0]+y Derivative[1][r][y]) (2 W[0]+y Derivative[1][Z][0]+U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])-2 (2 YY[0]+Derivative[1][Z][0]) (y (2 YY[0]+Derivative[1][Z][0])+3 U''[0])+((r[y]+Z[0]+y Derivative[1][r][y]) (y (2 YY[0]+Derivative[1][Z][0])+3 U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])+1/d y (Derivative[1][Z][0]^2-(4 (r[y]+Z[0]+y Derivative[1][r][y]) Derivative[1][Z][0] (2 W[0]+y YY[0]+U''[0]))/(y r[y]+y Z[0]+Derivative[1][U][0])+(2 (2 (r[y]+Z[0]+y Derivative[1][r][y])^2-(y r[y]+y Z[0]+Derivative[1][U][0]) (2 Derivative[1][r][y]+y r''[y])) (2 W[0]+y YY[0]+U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])^2)+1/d y ((-2 YY[0]+Derivative[1][Z][0])^2+(4 (r[y]+Z[0]+y Derivative[1][r][y]) (-2 YY[0]+Derivative[1][Z][0]) (2 W[0]+y YY[0]+U''[0]))/(y r[y]+y Z[0]+Derivative[1][U][0])+(2 (2 (r[y]+Z[0]+y Derivative[1][r][y])^2-(y r[y]+y Z[0]+Derivative[1][U][0]) (2 Derivative[1][r][y]+y r''[y])) (2 W[0]+y YY[0]+U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])^2)-1/d y (Derivative[1][Z][0]^2-(4 (r[y]+Z[0]+y Derivative[1][r][y]) Derivative[1][Z][0] (2 W[0]+y Derivative[1][Z][0]+U''[0]))/(y r[y]+y Z[0]+Derivative[1][U][0])+(2 (2 (r[y]+Z[0]+y Derivative[1][r][y])^2-(y r[y]+y Z[0]+Derivative[1][U][0]) (2 Derivative[1][r][y]+y r''[y])) (2 W[0]+y Derivative[1][Z][0]+U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])^2)-1/d y ((2 YY[0]+Derivative[1][Z][0])^2-(4 (r[y]+Z[0]+y Derivative[1][r][y]) (2 YY[0]+Derivative[1][Z][0]) (y (2 YY[0]+Derivative[1][Z][0])+3 U''[0]))/(y r[y]+y Z[0]+Derivative[1][U][0])+(2 (2 (r[y]+Z[0]+y Derivative[1][r][y])^2-(y r[y]+y Z[0]+Derivative[1][U][0]) (2 Derivative[1][r][y]+y r''[y])) (y (2 YY[0]+Derivative[1][Z][0])+3 U''[0])^2)/(y r[y]+y Z[0]+Derivative[1][U][0])^2));
flowW1J=-(1/(2 \[Rho]))vd y^(d/2) (\[Eta] r[y]+2 y Derivative[1][r][y]) ((2 W[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]])/(2 \[Rho] W[\[Rho]]+y (r[y]+Z[\[Rho]])+Derivative[1][U][\[Rho]])^2+(y (2 YY[\[Rho]]+2 \[Rho] Derivative[1][YY][\[Rho]]+Derivative[1][Z][\[Rho]])+3 U''[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]])/(y (r[y]+2 \[Rho] YY[\[Rho]]+Z[\[Rho]])+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2+2 (-((2 y YY[\[Rho]]+12 \[Rho] Derivative[1][W][\[Rho]]+y Derivative[1][Z][\[Rho]]+3 U''[\[Rho]])/(2 (y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2))+(4 \[Rho] (2 W[\[Rho]]+y YY[\[Rho]]+2 \[Rho] Derivative[1][W][\[Rho]]+U''[\[Rho]])^2 (y r[y]+\[Rho] W[\[Rho]]+y \[Rho] YY[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]]+\[Rho] U''[\[Rho]]))/((y r[y]+2 \[Rho] W[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]])^2 (y r[y]+2 y \[Rho] YY[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2)-(2 W[\[Rho]]+10 \[Rho] Derivative[1][W][\[Rho]]+2 y \[Rho] Derivative[1][YY][\[Rho]]+y Derivative[1][Z][\[Rho]]+U''[\[Rho]]+4 \[Rho]^2 (W^\[Prime]\[Prime])[\[Rho]]+2 \[Rho] \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]])/(2 (y r[y]+2 y \[Rho] YY[\[Rho]]+y Z[\[Rho]]+Derivative[1][U][\[Rho]]+2 \[Rho] U''[\[Rho]])^2)));


(* ::Text:: *)
(*Rules for simplification*)


If[debug,
fullA[f_]:=(f/(2y)/.repA0/.repA /. repProp);
fullJ[f_]:=(f//.repJ/.repProp);
If[ON,
	WflowA = 0;
	ZflowA = fullA[dimlessFlowEquations[[3,3]]];
	YflowA = fullA[(dimlessFlowEquations[[2,3]]-dimlessFlowEquations[[3,3]])/(2\[Rho])],
	WflowA = fullA[dimlessFlowEquations[[2,3]]];
	ZflowA = fullA[dimlessFlowEquations[[4,3]]];
	YflowA = fullA[(dimlessFlowEquations[[3,3]]-dimlessFlowEquations[[4,3]])/(2\[Rho])];
	WflowJ = fullJ[flowW1J];];
ZflowJ = fullJ[flowZJ];
YflowJ = fullJ[flowYJ];
gpsCoef[f_,i_,j_] := SeriesCoefficient[SeriesCoefficient[f,{gp,0,-i}],{gs,0,-j}];

If[!ON,
	Print[Simplify[FullSimplify[Simplify[(WflowA-WflowJ)/(2 vd y^(1/2 (d)) (\[Eta] r[y]+2 y Derivative[1][r][y]))]/.(propReverse/.repA/.q^2->y/.R[\[Lambda],y]->y r[y])]]]];
Print[FullSimplify[Collect[ZflowA-ZflowJ,{gs,gp},Simplify]/.(propReverse/.repA/.q^2->y/.R[\[Lambda],y]->y r[y])]];
Print[FullSimplify[Collect[YflowA-YflowJ,{gs,gp},Simplify]/.(propReverse/.repA/.q^2->y/.R[\[Lambda],y]->y r[y])]];]


(* ::Subsection::Closed:: *)
(*Export to LaTeX*)


If[n\[Phi]==0,
PropRep = {(r[y^2]+V[\[Rho]]+y^2 Z[\[Rho]]+2 \[Rho] Derivative[1][V][\[Rho]]) ->Subscript[G,\[Sigma]][\[Rho]]^-1,(r[y^2]+V[\[Rho]]+y^2 Z[\[Rho]]+2 \[Rho] Derivative[1][V][\[Rho]]) ->Subscript[G,\[Sigma]][\[Rho]],
	r[y^2]+V[\[Rho]]+y^2 (Z[\[Rho]]-2 \[Rho] Y[\[Rho]]) ->Subscript[G,\[Pi]][\[Rho]]^-1 , (r[y^2]+V[\[Rho]]+y^2 Z[\[Rho]] + y^2 \[Rho] Y[\[Rho]]+\[Rho] V'[\[Rho]]) ->1/2( Subscript[G,\[Pi]][\[Rho]]^-1 +Subscript[G,\[Sigma]][\[Rho]]^-1 ),
	(r[y^2]+V[\[Rho]]+y^2 (Z[\[Rho]] + \[Rho] Y[\[Rho]])+ \[Rho]( Subscript[W, 1][\[Rho]]+ V'[\[Rho]]))->1/2(Subscript[G,\[Pi]][\[Rho]]^-1 +Subscript[G,\[Sigma]][\[Rho]]^-1 ),(r[y^2]+V[0]+y^2 Z[0] )->g^-1};,
PropRep = {(r[y^2]+2\[Kappa] U2+y^2 Z)->gs^-1,(r[y^2]+y^2 Z)->gp^-1};
];


GCollect[F_] :=
Collect[F/.PropRep,{gs,gp},FullSimplify]


TeXReady[F_] := TeXForm[GCollect[GCollect[
F/((2 r[y^2]-\[Eta] r[y^2]-2 y^2 Derivative[1][r][y^2])vd y^(d-1))]]]


GCoef[F_,n_,m_] := FullSimplify[Coefficient[F,gp^n gs^m  ]/.{gs->0,gp->0}]


(*x=2;
Print[flowEquations[[x,1]]]
Print[TeXForm[flowEquations[[x,2]]]]
Print[GCollect[TeXReady[flowEquations[[x,3]]]]]
Remove[x];*)


(* ::Subsection::Closed:: *)
(*Export to C++ code*)


If[n\[Phi]==0,
CRepU = {V[\[Rho]_] ->v1,V'[\[Rho]_]->v2,V''[\[Rho]_]->v3,Z[\[Rho]_] ->z,Z'[\[Rho]_]->z1,Z''[\[Rho]_]->z2,Y[\[Rho]_] ->yy,Y'[\[Rho]_]->yy1,Y''[\[Rho]_]->yy2, Subscript[W, 1][\[Rho]_]->w, Subscript[W, 1]'[\[Rho]_]->w1, Subscript[W, 1]''[\[Rho]_]->w2,  Subscript[W, 2][\[Rho]_]->0},
CRepU = {Z->1,Z1->0,Z2->0,Y->0,Y1->0,Y2->0};
For[i=1,i<=n\[Phi]/2,i++, 
	AppendTo[CRepU, Subscript[U,i]->ToExpression["U"<>ToString[i]]];
	AppendTo[CRepU, Subscript[U,i,j_]->ToExpression["U"<>ToString[i]<>ToString[j]]]];];
CRepR = {r'[y_]->rp[y],r''[y_]->rp2[y], \[Rho]->rho};
CRep = Join[CRepU,CRepR];


If[n\[Phi]==0,
PropRep = {(r[y]+V[\[Rho]]+y Z[\[Rho]]+2 \[Rho] Derivative[1][V][\[Rho]]) ->gs^-1,(r[y]+V[\[Rho]]+y Z[\[Rho]]+2 \[Rho] Derivative[1][V][\[Rho]]) ->gs^-1,r[y]+V[\[Rho]]+y (Z[\[Rho]]-2 \[Rho] Y[\[Rho]]) ->gp^-1 , (r[y]+V[\[Rho]]+y Z[\[Rho]] + y \[Rho] Y[\[Rho]]+\[Rho] V'[\[Rho]]) ->1/2( gp^-1 + gs^-1),(r[y]+V[\[Rho]]+y (Z[\[Rho]] + \[Rho] Y[\[Rho]])+ \[Rho]( Subscript[W, 1][\[Rho]]+ V'[\[Rho]]))->1/2( gp^-1 + gs^-1),(r[y]+V[0]+y Z[0] )->g^-1};,
PropRep = {(r[y^2]+2\[Kappa] U2+y^2 Z)->gs^-1,(r[y^2]+y^2 Z)->gp^-1};
];


PropRep = PropRep/.CRep;


GCollect[F_] :=
Collect[F/.PropRep,{gs,gp},FullSimplify]


CReady[F_] := CForm[GCollect[GCollect[
F/((2 r[y^2]-\[Eta] r[y^2]-2 y^2 Derivative[1][r][y^2])vd y^(d-1))/.CRep]]/.{y^2->y2,\[Kappa]->K}]


GCoef[F_,n_,m_] := FullSimplify[Coefficient[F,gp^n gs^m  ]/.{gs->0,gp->0}]


(*x=3;
flowEquations[[x,1]]
GCollect[CReady[flowEquations[[x,3]]]]
Remove[x];*)


(* ::Subsection:: *)
(*Finite grid representation for the equations*)


(* ::Text:: *)
(*Formulates equations in finite grid representation in functional calculations*)


If[n\[Phi]==0,
Block[{k},
finiteGridReplacements = {V -> v, Subscript[W, i_] -> Subscript[w, i], Zs -> zs, Zp -> zp, Subscript[Zs, i_] ->Subscript[zs, i],Subscript[Zp, i_] ->Subscript[zp, i],
	Y -> yy,Subscript[yy, i]->Subscript[yy, i_], T ->t};
finiteGrid = {};
For[k=1, k<=Length[finiteGridReplacements], k++,
	AppendTo[finiteGrid, finiteGridReplacements[[k,1]][\[Rho]_] -> finiteGridReplacements[[k,2]][\[Rho]/eps]];
	AppendTo[finiteGrid, \!\(\*SuperscriptBox[\(finiteGridReplacements[\([k, 1]\)]\), 
TagBox[
RowBox[{"(", "j_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 1/eps^j \!\(\*SuperscriptBox[\(finiteGridReplacements[\([k, 2]\)]\), 
TagBox[
RowBox[{"(", "j", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]/eps]];
	];
];
gridFlowEquations = dimlessFlowEquations /. finiteGrid;
];


If[BooleanQ[LPA]&&LPA, comment="LPA"];
If[BooleanQ[strict]&&strict, comment="strict",comment=""];

filename = GetEquationsFileName[task, n\[Phi], n\[Tau], comment];
Export[filename, gridFlowEquations];
