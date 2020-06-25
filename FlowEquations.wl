(* ::Package:: *)

(* ::Title:: *)
(*NPRG Procedure for O(2) Model with Z4-symmetric Anisotropy*)


(* ::Subtitle:: *)
(*Functional Approach within Local Potential Approximation*)


(* ::Section:: *)
(*Initialization*)


arguments = $CommandLine;
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
anisotropySymmetry = Null;
Block[{i},
For[i = 1, i<=Length[arguments], i++,
	If[arguments[[i]] == "Zp", 
		If[arguments[[i+1]] == "ON",
			ON = True;
			anisotropySymmetry = "ON";
			n\[Tau] = 0,
			ON = False;
			anisotropySymmetry = ToExpression[arguments[[i+1]]]]];
	If[arguments[[i]] == "phi", n\[Phi] = ToExpression[arguments[[i+1]]]];
	If[arguments[[i]] == "tau", n\[Tau] = ToExpression[arguments[[i+1]]]];
	];
];

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


(*n\[Phi]=0; n\[Tau]=0; ON=True; anisotropySymmetry="ON";*)


(* ::Section:: *)
(*Invariants, Effective Action and Code Organizing Functions*)


(* ::Subsection::Closed:: *)
(*Invariants*)


If[ON,
\[Rho][m_, x_] := Sum[m[i][x]^2, {i, 1, 3}]/2;
\[Rho][\[Phi]] := Sum[\[Phi][i]^2, {i, 1, 3}]/2,
\[Rho][m_, x_] := Sum[m[i][x]^2, {i, 1, 2}]/2;
\[Rho][\[Phi]] := Sum[\[Phi][i]^2, {i, 1, 2}]/2];


(* Choice of anisotropic invariant *)
If[anisotropySymmetry == 4 || anisotropySymmetry == "ON", 
	\[Tau][m_, x_] := m[1][x]^2*(m[2][x]^2/2); \[Tau][\[Phi]] := \[Phi][1]^2*(\[Phi][2]^2/2); 
	VariablesRange = \[Rho] > 0 && \[Rho]^2 > 2*\[Tau] && \[Phi][1] >= \[Phi][2] && \[Phi][2] >= 0 && \[Phi][1] >= 0 && \[Tau] > 0; 
	
	(* Solution for \[Phi][1] and \[Phi][2] expressed by \[Rho] and \[Tau] *)
    disentangled = Simplify[Solve[{\[Rho] == (\[Phi][1]^2 + \[Phi][2]^2)/2, \[Tau] == (1/2)*\[Phi][1]^2*\[Phi][2]^2}, {\[Phi][1], \[Phi][2]}, Reals][[4]], Assumptions -> VariablesRange]; 
    
    (* Expressions appearing in the flow equations that need to be expressed in terms of invariants *)
    repExpressions = {\[Phi][1]^2*\[Phi][2]^2, \[Phi][1]^4*\[Phi][2]^4, \[Phi][1]^2 + \[Phi][2]^2, \[Phi][1]^4 + \[Phi][2]^4, (\[Phi][1]^2 - \[Phi][2]^2)^2, \[Phi][1]^4 - \[Phi][2]^4}; 
    ]; 

If[anisotropySymmetry == 6, 
	\[Tau][m_, x_] := (m[2][x]^3 - 3*m[1][x]^2*m[2][x])^2/4; \[Tau][\[Phi]] := (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2/4; 
    VariablesRange = \[Rho] > 0 && 2*\[Rho]^3 > \[Tau] && \[Tau] > 0 && \[Phi][1] >= \[Phi][2] && \[Phi][2] >= 0 && \[Phi][1] >= 0; 
    
    (* Solution for \[Phi][1] and \[Phi][2] expressed by \[Rho] and \[Tau] *)
    disentangled = Simplify[Solve[{\[Rho] == (\[Phi][1]^2 + \[Phi][2]^2)/2, \[Tau] == (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2/4}, {\[Phi][1], \[Phi][2]}, Reals][[12]], Assumptions -> VariablesRange];
    
    (* Expressions appearing in the flow equations that need to be expressed in terms of invariants *)
    repExpressions = {(\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^2, (\[Phi][2]^3 - 3*\[Phi][1]^2*\[Phi][2])^4, \[Phi][1]^2 + \[Phi][2]^2, \[Phi][1]^6 - 15*\[Phi][1]^4*\[Phi][2]^2 + 15*\[Phi][1]^2*\[Phi][2]^4 - \[Phi][2]^6, 
      (3*\[Phi][1]^5*\[Phi][2] - 10*\[Phi][1]^3*\[Phi][2]^3 + 3*\[Phi][1]*\[Phi][2]^5)^2, (3*\[Phi][1]^7*\[Phi][2] - 7*\[Phi][1]^5*\[Phi][2]^3 - 7*\[Phi][1]^3*\[Phi][2]^5 + 3*\[Phi][1]*\[Phi][2]^7)^2, 
      3*\[Phi][1]^6 - 36*\[Phi][1]^4*\[Phi][2]^2 + 39*\[Phi][1]^2*\[Phi][2]^4 - 2*\[Phi][2]^6, (-3*\[Phi][1]^4*\[Phi][2] - 2*\[Phi][1]^2*\[Phi][2]^3 + \[Phi][2]^5)^2, (\[Phi][1]^3 - 3*\[Phi][1]*\[Phi][2]^2)^2, 
      (-3*\[Phi][1]^4 - 2*\[Phi][1]^2*\[Phi][2]^2 + \[Phi][2]^4)^2*\[Phi][2]^2, -3*\[Phi][1]^6 + 36*\[Phi][1]^4*\[Phi][2]^2 - 39*\[Phi][1]^2*\[Phi][2]^4 + 2*\[Phi][2]^6}; 
    ]; 

varRep = Table[repExpressions[[i]] -> FullSimplify[repExpressions[[i]] /. disentangled, Assumptions -> VariablesRange], {i, 1, Length[repExpressions]}]; 


(* ::Subsection::Closed:: *)
(*Effective action (function of \[Rho] and \[Tau])*)


If[ON,
\[CapitalGamma]:= int[U[\[Rho][x],\[Tau][x]] + (Z[\[Rho][x]]-2 \[Rho][m, x] Y[\[Rho][x]])/2  Sum[D[m[i][x], x]^2, {i, 1, 3}] + Y[\[Rho][x]]/2 (Sum[m[i][x] D[m[i][x], x], {i, 1, 3}])^2, x],
\[CapitalGamma]:= int[U[\[Rho][x],\[Tau][x]] + (Z[\[Rho][x]]-2 \[Rho][m, x] Y[\[Rho][x]])/2  Sum[D[m[i][x], x]^2, {i, 1, 2}] + Y[\[Rho][x]]/2 (Sum[m[i][x] D[m[i][x], x], {i, 1, 2}])^2, x]
]


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
];];
dimensionlessRep = {d[q, q]->y, q^n_->y^(n/2), \[Lambda]->1, Z'[\[Lambda]]->-\[Eta], Z[\[Lambda]]->1};]


(* ::Subsection::Closed:: *)
(*Homogeneous configuration*)


homogeneous := {\[Rho][x_] -> \[Rho], \!\(\*SuperscriptBox[\(\[Rho]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0, \[Tau][x_] -> \[Tau], \!\(\*SuperscriptBox[\(\[Tau]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0, m[i_][x_] -> \[Phi][i], \!\(\*SuperscriptBox[\(m[i_]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[x_] -> 0}
minimal := {m[i_][x_] -> KroneckerDelta[i,1] Sqrt[2\[Rho]], \!\(\*SuperscriptBox[\(m[i_]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
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
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_], \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho],\[Tau]], \!\(\*SuperscriptBox[\(T\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_], Z[\[Rho]], \!\(\*SuperscriptBox[\(Z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], Y[\[Rho]], \!\(\*SuperscriptBox[\(Y\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]]};
DersCollect := Collect[#, PotentialDers, FullSimplify[#, Assumptions -> VariablesRange]& ]& 
Clean[q_List] := CollectMomenta[q][#/. homogeneous]&
CleanInvariants := DersCollect[#/. homogeneous ] //. varRep& 
Dimensionless := 2 vd y^(d/2 - 1)# /. regulatorRep /. dimensionlessRep&
Minimal := Simplify[# /. { \[Phi][i_] -> KroneckerDelta[i,1]Sqrt[2\[Rho]], \[Tau] -> 0}]&
ONSym := Simplify[# /. {\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_,\[Tau]_] /; n>0 -> 0, \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho],\[Tau]] /; n>0 -> 0, \!\(\*SuperscriptBox[\(T\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, \[Tau]_] -> 0, \[Tau] -> 0}]&
EvAt\[Rho]0 := (SeriesCoefficient[# /. \[Tau] -> 0, {\[Rho], 0, 0}])&


(* ::Section::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Laplacian lap[f_,p_,q_]*)


(* ::Text:: *)
(*Returns 1/(2 d)(Subscript[\[CapitalDelta], p]f )\!\(\*SubscriptBox[\("|"\), \("p=0"\)]\), q is the integration variable in the flow equation*)


$Assumptions = d>1;


lap[f_, p_, q_] := Module[{der, simple, res, repped},
repped= f /.{d[s_, r_] -> Sum[r[j]s[j], {j, 1, d}], R[\[Lambda],(a_. r_+b_. s_)^2] -> R[\[Lambda], Sum[(a r[i]+b s[i])^2, {i, 1, d}]],
	\!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"n_", ",", "m_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda], (a_. r_+b_. s_)^2] -> \!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"n", ",", "m"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda], Sum[(a r[i]+b s[i])^2, {i, 1, d}]]};
der = D[repped, {p[1], 2}] ;
simple = der  /. {p[i_] -> 0, p -> 0} /. {Sum[q[j_]^2, {j_, 1, d}] -> q^2};
res = 1/2(SeriesCoefficient[simple, {q[1], 0, 0}] + q^2 / d SeriesCoefficient[simple, {q[1], 0, 2}]);
Return[res]
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
(*Functional Derivative in Momentum Space (ffder[f, m_[q_]])*)


(* Definition for magnetic field *)
fder[m_[i_][x_], m_[j_][q_]] := Exp[I q x] KroneckerDelta[i-j]
fder[m_[i_]'[x_], m_[j_][q_]] := I q Exp[I q x]KroneckerDelta[i-j] 
fder[m_[i_]''[x_], m_[j_][q_]] := -d[q, q] Exp[I q x] KroneckerDelta[i-j]


(* Functional derivative in direct space *)
fder[int[f_, x_], m_[y_]] := int[fder[f, m[y]], x]  (* Commutative with integration *)
fder[Sum[f_, {i_, 1, n_}], m_[y_]] := Sum[fder[f, m[y]], {i, 1, n}]  (* Commutative with summation *)
fder[Tr[f_], m_[y_]] := Tr[fder[f, m[y]]]  (* Commutative with summation *)
fder[f_List, m_[y_]] := Table[fder[f[[i]], m[y]], {i, 1, Length[f]}]  (* Commutative with summation *)
fder[f_^n_, m_[y_]]:=n f^(n-1) fder[f, m[y]] (* Power derivative *)
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
fder[R[\[Lambda],q_], m_[y_]] := 0
fder[DDelta[x_], m_[y_]] := 0
fder[DDelta'[x_], m_[y_]] := 0  
fder[DDelta''[x_], m_[y_]] := 0  


(* ::Section::Closed:: *)
(*Correlation Functions*)


(* ::Subsubsection:: *)
(*Superscript["\[CapitalGamma]","(2)"] and propagator functions (without DiracDelta)*)


raw\[CapitalGamma]2[q_,p_]:=Table[fder[fder[\[CapitalGamma],m[i][q]],m[j][p]],{i,1,3},{j,1,3}];
\[CapitalGamma]2[q_,p_]:=FullSimplify[Clean[{p,q}][raw\[CapitalGamma]2[q,p]]];
rawG[q_,p_] := Inverse[raw\[CapitalGamma]2[q,p]+{{DDelta[p+q] R[\[Lambda],q^2],0,0},{0,DDelta[p+q] R[\[Lambda],q^2],0},{0,0,DDelta[p+q] R[\[Lambda],q^2]}}]
G[q_,p_] := Simplify[Inverse[\[CapitalGamma]2[q,p]+{{DDelta[p+q] R[\[Lambda],q^2],0,0},{0,DDelta[p+q] R[\[Lambda],q^2],0},{0,0,DDelta[p+q] R[\[Lambda],q^2]}}]]


(* ::Subsubsection:: *)
(*Superscript[\[CapitalGamma],(3)] function*)


raw\[CapitalGamma]3[p_,q_,r_] :=Module[{n},Table[fder[raw\[CapitalGamma]2[q,r],m[n][p]],{n,1,3}]];
\[CapitalGamma]3[p_,q_,r_]:=FullSimplify[Clean[{p,q,r}][raw\[CapitalGamma]3[p,q,r]]];


(* ::Subsubsection:: *)
(*Superscript[\[CapitalGamma],(4)]function*)


raw\[CapitalGamma]4[p_,q_,r_,s_] :=Module[{n},Table[fder[raw\[CapitalGamma]3[q,r,s],m[n][p]],{n,1,3}]];
\[CapitalGamma]4[p_,q_,r_,s_]:=FullSimplify[Clean[{p,q,r,s}][raw\[CapitalGamma]4[p,q,r,s]]];


(* ::Input:: *)
(*(*\[CapitalGamma]4 = Table[DD[{m,n,l,o}][U[\[Rho],\[Tau]]],{m,1,2},{n,1,2},{l,1,2},{o,1,2}] // AbsoluteTiming;*)*)


(* ::Section:: *)
(*Flow equations*)


(* ::Subsection:: *)
(*Propagators in minimal configuration*)


Propagator[p_, q_] := Minimal[DDelta[p+q]^2 G[p, q]]
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


props = {Propq[[1]][[1]], Propq[[2]][[2]], 1/gs, 1/gp}; 
propRep = {1/props[[1]] -> gs, Expand[1/props[[1]]] -> gs, Simplify[1/props[[2]]] -> gp, Expand[1/props[[2]]] -> gp,
	Expand[\[Rho]/props[[2]]] -> \[Rho] gp, Simplify[(1/props[[1]] + 1/props[[2]])/2] -> (gp + gs)/2, Expand[(1/props[[1]] + 1/props[[2]])/2] -> (gp + gs)/2};
propRep = Join[propRep, propRep /. d[q_,p_] -> q p];
propReverse = {gs -> 1/props[[1]], gp -> 1/props[[2]]};
ReplaceProps := Collect[Collect[#/.propRep, props, Simplify] /. propRep, props, FullSimplify] /. propRep&


(* ::Subsection:: *)
(*Wetterich equation (Exact)*)


If[ON, 
	
	wetterich = 1/2 D[R[\[Lambda],q^2],\[Lambda]] (Gs + (n-1)Gp),
	
	(*We calculate trace of the modified propagator*)
	prop = G[q,p];
	traced = 1/2 Together[prop[[1]][[1]] + prop[[2]][[2]]];

	(*We separately simplify the numerator and denominator of the trace.*)
	GNum = Numerator[traced];
	GDen = Denominator[traced];
	
	
	GNumCollected = FullSimplify[DersCollect[GNum] //. varRep] //. varRep;
	GDenCollected = FullSimplify[DersCollect[GDen] //. varRep] //. varRep;
	(*We recover full Wetterich equation at O(\[PartialD]^2)*)
	wetterich = D[R[\[Lambda],q^2],\[Lambda]]  Simplify[GNumCollected/GDenCollected] //. varRep;
];

wetterich = wetterich /. p->-q /. DDelta[0] -> 1;


(* ::Subsection:: *)
(*Flow equations for Local Potential Expanded at \[Tau]=0*)


flowEquations = {};


If[n\[Phi]==0,
DerRep = {\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, \[Tau]_] /; m>0 -> \!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", 
RowBox[{"m", "-", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]], \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, \[Tau]_] /; n>n\[Tau] -> 0, \!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"m_", ",", "n_"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_, \[Tau]_] -> \!\(\*SuperscriptBox[
SubscriptBox[\(W\), \(n\)], 
TagBox[
RowBox[{"(", "m", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]]};

VEquation = Collect[D[FullSimplify[wetterich /.\[Tau]->0], \[Rho]], GDenCollected, Simplify];
AppendTo[flowEquations, {V[\[Rho]], (-2+\[Eta])V[\[Rho]]+(-2+d+\[Eta]) \[Rho] V'[\[Rho]], Dimensionless[VEquation] /. DerRep}];

If[!ON,
For[ii=1, ii<=n\[Tau], ii++,
WiEquation = ii! ReplaceProps[SeriesCoefficient[wetterich, {\[Tau], 0, ii}]];
AppendTo[flowEquations,{Subscript[W, ii][\[Rho]],
	(anisotropySymmetry (d-2+\[Eta]) ii/2 - d) Subscript[W, ii][\[Rho]] + (-2+d+\[Eta]) \[Rho] Subscript[W, ii]'[\[Rho]],
	Dimensionless[WiEquation /. propReverse] /. DerRep}];
]];
];


(* ::Subsection::Closed:: *)
(*Flow equations for Local Potential Expanded at \[Tau]=0 and \[Rho]=\[Kappa] at LPA'*)


If[n\[Phi]!=0,
LPA = {Z[\[Rho]_] -> Z, \!\(\*SuperscriptBox[\(Z\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0, Y[\[Rho]_] -> 0, \!\(\*SuperscriptBox[\(Y\), 
TagBox[
RowBox[{"(", "i_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 0};
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
MultilineFunction->None]\)[\[Rho]_,\[Tau]_]->Subscript[u, i, j]}, LPA];

wetterichLPA = wetterich /. LPA;

wetterichTauExpansion = Table[ii!SeriesCoefficient[wetterichLPA, {\[Tau], 0, ii}], {ii, 0, n\[Tau]}];
wetterichFullExpansion =Table[Dimensionless[jj! SeriesCoefficient[wetterichTauExpansion[[ii+1]],{\[Rho],\[Kappa],jj}]/.DerRep ],{ii,0,n\[Tau]},{jj,0,(n\[Phi]-anisotropySymmetry ii)/2}];

dt\[Kappa] =- wetterichFullExpansion[[1]][[2]]/Subscript[u, 2,0];

equationsExpanded = Table[If[i==0&&j==0, {0,0,0}, If[i==1&&j==0, {\[Kappa],2-d-\[Eta],dt\[Kappa]}, 
	{Subscript[u, i,j],Expand[(i +anisotropySymmetry/2 j)(d-2+\[Eta])-d], dt\[Kappa] Subscript[u, i+1,j] + wetterichFullExpansion[[j+1]][[i+1]]}]],
	{j,0,n\[Tau]},{i,0,(n\[Phi]-anisotropySymmetry j)/2}];
	
flowEquations = Drop[Flatten[equationsExpanded,1],1] /. {Subscript[u, i_,j_] /; 2i + anisotropySymmetry j>n\[Phi]->0};
];


(* ::Subsection:: *)
(*Equation for flow of Z[\[Rho]] and Y[\[Rho]]*)


(* ::Text:: *)
(*We assume Subscript["\[PartialD]", \[Tau]]Z =0, so equations are evaluated at configuration \[Tau]=0, \[Phi][2]=0.*)


(* ::Subsubsection:: *)
(*Vertices in minimal configuration*)


M\[CapitalGamma]3 = Simplify[Minimal[Clean[{p,q,r}][raw\[CapitalGamma]3[p,q,r]]] /. r->-p-q /. {d[p,q] -> d[q,p], DDelta[0] -> 1}];
M\[CapitalGamma]3 = DersCollect[M\[CapitalGamma]3];


M\[CapitalGamma]4 = Simplify[Minimal[Clean[{p,q,r,s}][raw\[CapitalGamma]4[p,q,r,s]]] /. {q->-p, r->q, s->-q} /. {d[p,q] -> d[q,p], DDelta[0] -> 1}];
M\[CapitalGamma]4 = DersCollect[M\[CapitalGamma]4];


If[ON,
M\[CapitalGamma]3=ONSym[M\[CapitalGamma]3];
M\[CapitalGamma]4=ONSym[M\[CapitalGamma]4];
];


(* ::Text:: *)
(*For flows of both Z and Y we have two contributing diagrams: tadpole with one \[CapitalGamma]4 vertex and bubble with two \[CapitalGamma]3 vertices.*)


(* ::Subsubsection:: *)
(*Flow of Z[\[Rho]]*)


(* ::Text:: *)
(*The tadpole diagram*)


tadpole = Simplify[-1/2(gs^-2 M\[CapitalGamma]4[[1,1,1,1]] + (n-1) gp^-2 M\[CapitalGamma]4[[1,1,2,2]] )];
tadpoleDifferentiated = ReplaceProps[lap[tadpole, p, q]];


(* ::Text:: *)
(*The bubble diagram*)


bubble =(M\[CapitalGamma]3[[1,1,1]]^2 gs^-2 (Gs /. q -> q+p) + (n-1) M\[CapitalGamma]3[[1,2,2]]^2 gp^-2 (Gp /. q -> q+p));
bubbleDifferentiated = ReplaceProps[lap[bubble, p, q]] ;


(* ::Text:: *)
(*Flow equation for Z:*)


Zflow = Dimensionless[\!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda],q^2] Collect[(tadpoleDifferentiated + bubbleDifferentiated), {1/gs, 1/gp}, FullSimplify]/.propReverse];
If[n\[Phi]==0,
ZEquation = {Z[\[Rho]], \[Eta] Z[\[Rho]]+(-2+d+\[Eta]) \[Rho] Z'[\[Rho]], Zflow/.DerRep},
ZEquation = {Z, \[Eta], Zflow /. DerRep /. \[Rho] -> \[Kappa]}
];
AppendTo[flowEquations, ZEquation];


(* ::Subsubsection:: *)
(*Flow equations for Y[\[Rho]]*)


If[n\[Phi]==0,
Ytadpole = Simplify[(tadpole-(-1/2(M\[CapitalGamma]4[[2,2,1,1]] gs^-2 + M\[CapitalGamma]4[[2,2,2,2]] gp^-2 + (n-2) M\[CapitalGamma]4[[2,2,3,3]] gp^-2)))/(2\[Rho])];
YtadpoleDifferentiated = ReplaceProps[lap[Ytadpole, p, q]];

Ybubble = Simplify[(bubble- M\[CapitalGamma]3[[2,2,1]]^2 (gp^-2 (Gs /. q -> q+p)) - M\[CapitalGamma]3[[2,1,2]]^2 (gs^-2(Gp /. q -> q+p)))/(2\[Rho])];
YbubbleDifferentiated = ReplaceProps[lap[Ybubble, p, q]];

Yflow = ReplaceProps[YtadpoleDifferentiated + YbubbleDifferentiated];
YflowSimplified = Dimensionless[\!\(\*SuperscriptBox[\(R\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Lambda],q^2] Collect[Yflow, {1/gs, 1/gp}, Simplify]/. propReverse ];

AppendTo[flowEquations, {Y[\[Rho]], (d-2+2\[Eta]) Y[\[Rho]] + (-2+d+\[Eta]) \[Rho] Y'[\[Rho]], YflowSimplified  /. DerRep}];

flowY\[Rho]0 = Collect[EvAt\[Rho]0[YflowSimplified /.DerRep]/.{r[y] + V[0] + y Z[0] -> g}, {g^-1}, FullSimplify] /. g -> r[y] + V[0] + y Z[0];
AppendTo[flowEquations, {Y[0], (d-2+2\[Eta]) Y[0], flowY\[Rho]0/.DerRep}];
];


(* ::Subsection:: *)
(*Apply chosen expansion in powers of anisotropic coupling*)


If[n\[Phi]==0,
flowEquations = flowEquations /. {Subscript[W, i_][\[Rho]_]/;i>n\[Tau] -> 0, 
	\!\(\*SuperscriptBox[\(Subscript[W, \ i_]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]/;i>n\[Tau] -> 0}];
If[!ON,
flowEquations = flowEquations /. n->2;];


(* ::Subsection:: *)
(*Finite grid representation for the equations*)


(* ::Text:: *)
(*Formulates equations in finite grid representation in functional calculations*)


If[n\[Phi]==0,
Block[{k},
finiteGridReplacements = {V -> v, Subscript[W, i_] -> Subscript[w, i], Z -> z, Y -> yy};
finiteGrid = {};
For[k=1, k<=Length[finiteGridReplacements], k++,
	AppendTo[finiteGrid, finiteGridReplacements[[k,1]][\[Rho]_] -> finiteGridReplacements[[k,2]][\[Rho]/eps]];
	AppendTo[finiteGrid, \!\(\*SuperscriptBox[\(finiteGridReplacements[\([k, 1]\)]\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_] -> 1/eps^n \!\(\*SuperscriptBox[\(finiteGridReplacements[\([k, 2]\)]\), 
TagBox[
RowBox[{"(", "n", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]/eps]];
	];
];
flowEquations = flowEquations /. finiteGrid;
];


(* ::Section:: *)
(*Cache the equations*)


filename = GetEquationsFileName[anisotropySymmetry, n\[Phi], n\[Tau]];
Export[filename, flowEquations];


(* ::Subsection:: *)
(*Consistency check for functional equations at O(\[Tau])*)


repJ ={YY->Y, W[\[Rho]_]->0,\!\(\*SuperscriptBox[\(W\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0};
repA = {\!\(\*SuperscriptBox[\(V\), 
TagBox[
RowBox[{"(", "n_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->\!\(\*SuperscriptBox[\(U\), 
TagBox[
RowBox[{"(", 
RowBox[{"n", "+", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]],Subscript[W, n_][\[Rho]_]->0,\!\(\*SuperscriptBox[
SubscriptBox[\(W\), \(n_\)], 
TagBox[
RowBox[{"(", "m_", ")"}],
Derivative],
MultilineFunction->None]\)[\[Rho]_]->0,V[\[Rho]_]->U'[\[Rho]],Z[\[Rho]_]->Z[\[Rho]]+2 \[Rho] Y[\[Rho]], Z'[\[Rho]_]->2 Y[\[Rho]]+2 \[Rho] Derivative[1][Y][\[Rho]]+Derivative[1][Z][\[Rho]],Z''[\[Rho]_]->4 Derivative[1][Y][\[Rho]]+2 \[Rho] Y''[\[Rho]]+Z''[\[Rho]],n->2, r[y]->y r[y], r'[y]->r[y] + y r'[y], r''[y]->2 r'[y]+y r''[y]};
repProp = {y r[y]+y Z[\[Rho]_]+Derivative[1][U][\[Rho]_]->gp,y (r[y]+ Z[\[Rho]_])+Derivative[1][U][\[Rho]_]->gp,y r[y]+y (2 \[Rho] Y[\[Rho]_]+Z[\[Rho]_])+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs,y (r[y]+2 \[Rho] Y[\[Rho]_]+Z[\[Rho]_])+Derivative[1][U][\[Rho]_]+2 \[Rho] U''[\[Rho]_]->gs};


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


(*0!=Collect[((flowEquations[[2,3]]/.repA/.repProp)-(2\[Rho] flowEquations[[3,3]]/.repA/.repProp)-((flowZJ/.repJ))/.repProp),{gs,gp},Simplify]
0!=Collect[(flowEquations[[3,3]]/.repA/.repProp)-((flowYJ/.repJ))/.repProp,{gs,gp},Simplify]*)
