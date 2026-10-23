(* v2 source-based comparison. No mathematical function bodies are rewritten.
   Serial/headless orchestration only. Requires a licensed Wolfram Language
   kernel and network access to the public, hash-checked StableWein notebook.
   Exact zero is not scored against StableWein's strict-positive convention.
   This script does NOT reproduce the Python timing environment. *)
ClearAll[runStableWeinV2];
runStableWeinV2[rows_List, outdir_String, modes_List:{4}, seconds_:30] := Module[
 {url,txt,nb,boxes,held,ids,hash,parse,toParameters,results={},row,p,par,m,flags,timed,rec,index,mode},
 url="https://raw.githubusercontent.com/jurciukonis/StableWein/main/StableWein.nb";
 hash="acf62cd0c9c4de76f2cef5e54b127071c085cf50f2cf5dcecac1c2aa91e8d098";
 txt=Import[url,"Text"];
 If[!StringQ[txt]||Hash[txt,"SHA256","HexString"]!=hash,Return[Failure["SourceChangedOrUnavailable",<|"url"->url|>]]];
 nb=ToExpression[txt,InputForm,HoldComplete];
 boxes=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];
 held=Map[Function[b,Quiet[If[ListQ[b],ToExpression[b,StandardForm,HoldComplete],{ToExpression[b,StandardForm,HoldComplete]}]]],boxes];
 ids={33,34,35,36,37,38,39,40,41,42,43,44,46,48,51,52};
 If[!FreeQ[held[[ids]],$Failed],Return[Failure["NotebookParse",<||>]]];
 wbLoadBinariesPac=False;wpCompTargetPac="WVM";
 Scan[ReleaseHold,Flatten[held[[ids]]]/.{ParallelMap->Map,ParallelTable->Table}];
 Clear[ShowStatus];ShowStatus[_String]:=Null;
 If[MemberQ[modes,3],{wbfaceListWein,wb\[Theta]ListWein,wb\[Theta]ListWein2}=wbWeinGridGen[{17,40,0}]];
 parse[s_String]:=If[StringMatchQ[s,RegularExpression["[+-]?[0-9]+(/[0-9]+)?"]],ToExpression[s,InputForm],Throw["Invalid exact coefficient"]];
 (* Allowed sign rephasings move the invariant product sign onto source edge 23.
    Arrays go from (12,23,31) to (23,31,12); source d magnitudes are nonnegative. *)
 toParameters[p_Association]:=Module[{a,b,c,d,phase},
   {a,b,c,d}=(parse/@p[#]&)/@{"a","b","c","d"};
   phase=If[Times@@d<0,Pi,0];N[Join[a,b[[{2,3,1}]],c[[{2,3,1}]],Abs[d[[{2,3,1}]]],{phase}]]];
 If[!DirectoryQ[outdir],CreateDirectory[outdir,CreateIntermediateDirectories->True]];
 Do[
  row=rows[[index]];p=row["couplings"];m=parse[row["minimum"]];par=toParameters[p];
  Do[
   SeedRandom[Lookup[row,"seed",2026101700+10 index+mode]];
   flags=<|"sufficient"->wbWeinSufCond[par],"NCL1"->wbWeinNecCond1[par],"NCL2"->wbWeinNecCond2[par]|>;
   timed=If[m==0,{0,"FLAT_NOT_SCORED"},AbsoluteTiming[TimeConstrained[
    Length[Last[wbWeinUniBFBcheck[{par},mode,"ImposeUNI"->False,"ShowMessages"->False]]],seconds,"TIMEOUT"]]];
   rec=Join[<|"id"->row["id"],"mode"->mode,"minimum"->row["minimum"],"seed"->Lookup[row,"seed",2026101700+10 index+mode],
    "seconds"->timed[[1]],"accepted_or_status"->timed[[2]],"expected_strict_positive"->(m>0),
    "boundary_sensitive"->(Abs[N[m]]/Max[1,Max[Abs[Take[par,12]]]]<10^-10),
    "numerical_stage_required_by_filters"->(mode==4&&!TrueQ[flags["sufficient"]]&&TrueQ[flags["NCL1"]]&&TrueQ[flags["NCL2"]]),
    "source_Wolfram_string_sha256"->hash,"kernel"->$Version|>,flags];
   Export[FileNameJoin[{outdir,ToString[index]<>"_mode"<>ToString[mode]<>".json"}],rec,"RawJSON"];
   AppendTo[results,rec],{mode,modes}],{index,Length[rows]}];
 <|"results"->results,"source"->url,"source_hash"->hash,"scope"->"Real coefficients only; no UNI cuts or neural networks. Timeouts and exact-flat cases are not passes."|>
];
(* Example, from the bundle root in a Wolfram session:
 Get["revision_v2/run_stablewein_benchmark.wl"];
 rows=Import["validation/revision_v2/stablewein_full_inputs.json","RawJSON"];
 report=runStableWeinV2[rows,"validation/stablewein_local",{3,4},30];
 Export["validation/stablewein_local/summary.json",report,"RawJSON"];
 This all-input run is provided for reproduction; it was NOT executed in v2.
*)
