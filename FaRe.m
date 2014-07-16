BeginPackage["FaRe`"]

If[StringCases[ToString[$Packages],"HighEnergyPhysics"]=={},
	Print["FeynCalc is not loaded. Please be sure FeynCalc is loaded before running FaRe."];
	Abort[]
	];

Print["\!\(\*
StyleBox[\"FaRe\",\nFontSize->16,\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"v1\",\nFontSize->14]\)\!\(\*
StyleBox[\".0\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"-\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"package\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"for\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"tensor\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"reduction\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"of\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"Feynman\",\nFontSize->14]\)\!\(\*
StyleBox[\" \",\nFontSize->14]\)\!\(\*
StyleBox[\"integrals\",\nFontSize->14]\)\!\(\*
StyleBox[\".\",\nFontSize->14]\)"];


(* °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° *)

PQ::usage =
	"PQ[ loopNumber_Integer ]

	PQ returns expression of the Schwinger exponent P/Q after the momentum transformations.
	It has only one argument: the number of loops.
	The output is in function of coefficient ac_ij and d_i."

TRed::usage =
	"TRed[ D,num,den,loopMomenta]
	
	TRed performs the integral reduction of the input integral.
		\[Bullet] D: is an integer representing the dimension D in which the integral must be computed.
		\[Bullet] num is a list containing the loop momenta that appear in the numerator of the integral.
		\[Bullet] den is a list containing the arguments of the propagators, and their relative exponent.
		  Each propagator must be written without the square. 
		\[Bullet] loopMomenta is a list containing the loop momenta.
	"

FIREType::usage =
	"FIREType[expr,irr,zeros,fLetter]

	FIREType converts the output of TRed into a form suitable for further computation with FIRE.
		\[Bullet] expr is a linear combination of scalar integrals as obtained from TRed output. 
		  It is convenient to remove the tensor structures, even if it is not strictly necessary.
		\[Bullet] irr is the number of irreducible numerators treated by FIRE. 
		  Their exponents will be appended to the \[Nu]i\[CloseCurlyQuote]s list treated by TRed and will appear as zeros.
		\[Bullet] zeros is the list of positions where to insert a zero in the propagator exponent list. 
		  This can be used if FIRE asks for a list of exponents longer than that used by FaRe, 
		  for instance because of simplifications.
		\[Bullet] fLetter specifies the name used to represent the integrals in a suitable way for FIRE.
	"

(* °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° *)

Begin["`Private`"]

g[mu_,nu_]:=HighEnergyPhysics`FeynCalc`MetricTensor`MetricTensor[mu,nu];
fourv[x_,m_]:=HighEnergyPhysics`FeynCalc`FourVector`FourVector[x,m];
d[i_,m_]:=HighEnergyPhysics`FeynCalc`FourVector`FourVector["d"<>ToString[i],m];
contr[x_]:=HighEnergyPhysics`fctools`Contract`Contract[x];
ind[i_]:="m"<>ToString[i];

xi[nui_Integer,s_Integer,loops_Integer]:=Pochhammer[nui,s]*(-1)^s;

PQ[loops_Integer]:=Block[
	{ac,f,k,K,dd,
	P,Q,
	Ex1,Ex2,ck2,ck,subk,
	m,n,al,be,
	i,j},

	dd=Table["d"<>ToString[i],{i,1,loops}];
	ac=Table[Table[If[i<j,"ac"<>ToString[i]<>ToString[j],
						"ac"<>ToString[j]<>ToString[i]],{j,1,loops}],{i,1,loops}];
	k:=Table["k"<>ToString[i],{i,1,loops}];
	K=Table["K"<>ToString[i],{i,1,loops}];
	f="f";
	
	Ex1=Sum[Sum[ac[[i]][[j]]k[[i]]k[[j]],{i,1,loops}],{j,1,loops}]+
		Sum[2*k[[i]]dd[[i]],{i,1,loops}];

	For[i=1,i<=loops,i++,
		ck2=Coefficient[Ex1,k[[i]]^2];
		ck=Coefficient[Ex1,k[[i]]];
		subk=K[[i]]-ck/(2 ck2);
		Ex2=Ex1/.k[[i]]->subk;
		Ex1=Ex2;
	];
	Q=Ex1-Sum[Coefficient[Ex1,K[[i]]^2]*K[[i]]^2,{i,1,loops}]//FullSimplify;
	Return[Q];
];

PQ1[loops_Integer,dd_]:=Block[
	{ac,f,k,K,
	P,Q,
	Ex1,Ex2,ck2,ck,subk,
	m,n,al,be,
	i,j},

	ac=Table[Table[If[i<j,"ac"<>ToString[i]<>ToString[j],
						"ac"<>ToString[j]<>ToString[i]],{j,1,loops}],{i,1,loops}];
	k:=Table["k"<>ToString[i],{i,1,loops}];
	K=Table["K"<>ToString[i],{i,1,loops}];
	f="f";
	
	Ex1=Sum[Sum[ac[[i]][[j]]k[[i]]k[[j]],{i,1,loops}],{j,1,loops}]+
		Sum[2*k[[i]]dd[[i]],{i,1,loops}];

	For[i=1,i<=loops,i++,
		ck2=Coefficient[Ex1,k[[i]]^2];
		ck=Coefficient[Ex1,k[[i]]];
		subk=K[[i]]-ck/(2 ck2);
		Ex2=Ex1/.k[[i]]->subk;
		Ex1=Ex2;
	];
	Q=Ex1-Sum[Coefficient[Ex1,K[[i]]^2]*K[[i]]^2,{i,1,loops}];
	Return[Q];
];

NumV[Q_,num_,loops_Integer,dd_]:=Block[
	{du,i,j,Exp1,Exp2},	
	
	deriv[i_,m_,x_]:=HighEnergyPhysics`fctools`FourDivergence`FourDivergence[x,d[i,m]];
	
	Exp1=Exp[Q/"P"];
	For[i=1,i<=Length[num],i++,
		du=HighEnergyPhysics`FeynCalc`FourVector`FourVector["d"<>ToString[num[[i]]],"m"<>ToString[i]];
		Exp1=HighEnergyPhysics`fctools`FourDivergence`FourDivergence[Exp1,du]/2;
	];
	Return[Exp1/Exp[Q/"P"]//FullSimplify];
];

RepV[Q_,loops_Integer,dd_]:=Block[
	{Q1,Q2,i,j},
	(* d[m_]:=Table[HighEnergyPhysics`FeynCalc`FourVector`FourVector["d"<>ToString[i],m],{i,1,loops}];*)
	Q1=Expand[Q];

	For[i=1,i<=loops,i++,
		For[j=1,j<=loops,j++,
			Q2=Q1/.{dd[[i]]dd[[j]]->contr[d[i,m]d[j,n]g[m,n]]};
			Q1=Q2;
		];
	];

	Return[FullSimplify[Q1]];
];

Repdi[d_,pnames_,len_,index_]:=Block[
	{i,j,
	d1},

	d1=d;
	For[i=1,i<=len,i++,
		If[pnames[[i]]!=0,
			d1=d1/.{ToExpression[pnames[[i]]]->fourv[pnames[[i]],index]};
		];
	];
	Return[d1];
];

FIREType[str1_,irr_,zeros_,fletter_]:=Block[
	{str,letter,A1,A2,A3,x,y,i},
	letter="\[ScriptCapitalI]";

	str=StandardForm[str1];
	(*A1=str/.{Subsuperscript[letter,y_,x_]:>"G["<>ToString[Flatten[Append[Insert[y,0,zeros],Table[0,{i,1,Length[irr]}]]]<>"]"]}; *)
	If[Length[zeros]==1,
		A1=str/.{Subsuperscript[letter,y_,x_]:>ToString[fletter]<>"["<>ToString[Flatten[Append[If[zeros[[1]]!=0,Insert[y,0,zeros],y],Table[0,{i,1,irr}]]]]<>"]"},
		A1=str/.{Subsuperscript[letter,y_,x_]:>ToString[fletter]<>"["<>ToString[Flatten[Append[Insert[y,0,zeros],Table[0,{i,1,irr}]]]]<>"]"};
	];
	A2=ToString[A1];
	Return[A2];
];

TRed[D_Integer,nuu_,de_,loop_]:=Block[
	{i,j,l,m,
		num,num1,den,den1,den2,nu,nu1,A,numk,numi,numx,numxstr,deni,denj,denpw,
		ac,acNames,di,diNames,
		pNames,p1,pq,mlnumx,mli,sp,p,
		P,Q,Q1,R,n,n1,nRep,ml,res,resi,resj,inti,
		minpos,letters,loops,len,x,idx,pw,coeff,coeffden,dime},

(*	If[Length[loop]>3,
		Print["No more than 2 loop momenta, so far."];
		Abort[]];
	For[i=1,i<=Length[loop],i++,
		If[Length[Position[de,loop[[i]]]]==0,
			Print["No loop momentum found in denominator"];
			Abort[]];
	];
*)
	den=de;
	loops=Length[loop];
	len=Length[den];
	letters={};
	For[i=1,i<=loops,i++,
		letters=Append[letters,"k"<>ToString[i]];
	];

	num=nuu;
	For[i=1,i<=loops,i++,
		num=num/.{loop[[i]]->letters[[i]]};
	];
	num1=Table[0,{i,1,Length[num]}];
	For[i=1,i<=Length[num],i++,
		For[j=1,j<=loops,j++,
			If[num[[i]]==letters[[j]],num1[[i]]=j];
		];
	];
	den=de;
	For[i=1,i<=loops,i++,
		den=den/.{loop[[i]]->letters[[i]]};
	];
	p={};
	den1=den;
	For[j=1,j<=loops,j++,
			den1=den1/.{letters[[j]]->0};
	];
	pNames=Table[0,{i,1,len}];
	For[i=1,i<=len,i++,
		p1=ToString[den1[[i]][[1]]];
		If[!DigitQ[p1],
			minpos=StringPosition[p1,"-"];
			If[Length[minpos]!=0,
				pNames[[i]]=StringDrop[p1,minpos[[1]]],
				pNames[[i]]=p1;
			];
		];
	];

	den2=Table[{},{i,1,len}];
	For[i=1,i<=len,i++,
		den2[[i]]=den[[i]][[1]];
	];

(*
	a=Table[{},{i,1,len}];
	b=Table[{},{i,1,len}];
	c=Table[{},{i,1,len}];
	d=Table[{},{i,1,len}];	
	dd=Table[{},{i,1,len}];	
	e=Table[{},{i,1,len}];
	ee=Table[{},{i,1,len}];
	f=Table[{},{i,1,len}];
*)
	nu=Table[den[[i]][[2]],{i,1,len}];
	(* x=Table[Symbol["x"<>ToString[i]],{i,1,len}]; *)
	x=Table[Module[{x},x],{i,1,len}];
	A=Sum[Expand[den[[i]][[1]]^2]*x[[i]],{i,1,len}];
	ac=Table[Table[0,{j,1,loops}],{i,1,loops}];
	acNames=Table[Table[0,{j,1,loops}],{i,1,loops}];
	di=Table[0,{i,1,loops}];	
	diNames=Table[0,{i,1,loops}];	

	For[i=1,i<=loops,i++,
		For[j=1,j<=loops,j++,
			acNames[[i]][[j]]="ac"<>ToString[i]<>ToString[j];
			If[i==j,
				ac[[i]][[j]]=Coefficient[A,letters[[i]]*letters[[j]]],
				ac[[i]][[j]]=Coefficient[A,letters[[i]]*letters[[j]]]/2;
			];
		];
	];
	For[j=1,j<=loops,j++,
		di[[j]]=Coefficient[A,letters[[j]]]/2//FullSimplify;
		diNames[[j]]="d"<>ToString[j];
		For[i=1,i<=loops,i++,
			di[[j]]=di[[j]]/.{letters[[i]]->0};
		];
	(*	For[i=1,i<=len,i++,
			di[[j]]=di[[j]]/.{den1[[i]][[1]]->p[[i]]};
		];*)
	];
	pq=PQ1[loops,diNames]//FullSimplify;
	Q=Numerator[pq];
	R=RepV[Q,loops,diNames];
	n=NumV[R,num1,loops,diNames];
	For[m=1,m<=len,m++,
		n=n/.{x[[m]]->ToExpression[x[[m]]]};
	];
	For[i=1,i<=loops,i++,
		For[j=1,j<=loops,j++,
			n=n/.{acNames[[i]][[j]]->ac[[i]][[j]]};
		];
	];
	For[idx=1,idx<=Length[num],idx++,
		For[i=1,i<=loops,i++,
			n=n/.fourv[diNames[[i]],ind[idx]]:>Repdi[di[[i]],pNames,len,ind[idx]];
		];
	];
	ml=MonomialList[n]//FullSimplify;
	res=Table[0,{j,1,Length[ml]}];
	nu1=Table[0,{j,1,len}];

	For[i=1,i<=Length[ml],i++,
		mli=ml[[i]];
		numi=Numerator[mli];
		numk=numi;
		deni=Denominator[mli];
		dime=Exponent[deni,"P"];
		deni=deni/."P"->1;
		For[j=1,j<=len,j++,
			numk=numk/.{x[[j]]->1};
		];
		numx=Coefficient[numi,numk];
		coeff=1;
		For[j=1,j<=len,j++,
			pw=Exponent[numx,x[[j]]];
			nu1[[j]]=nu[[j]]+pw;
			coeff*=xi[nu[[j]],pw,loops];
		];
		inti=coeff*Subsuperscript["\[ScriptCapitalI]",nu1,D+2*dime];
		res[[i]]=numk/deni*inti;
	];

	Return[Sum[res[[i]],{i,1,Length[res]}]];
	
];

End[ ]

EndPackage[ ]


