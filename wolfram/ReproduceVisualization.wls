ClearAll[fail, styleMesh, renderMesh];
fail[code_Integer, message_String] := (Print["ERROR: " <> message]; Exit[code]);
styleMesh[mesh_] := mesh /. GraphicsComplex[points_, primitives_, options___] :>
  GraphicsComplex[points, {Directive[RGBColor[0.78, 0.64, 0.34],
    Specularity[White, 24]], EdgeForm[], primitives}, options];
renderMesh[mesh_, view_, size_: 900] := Rasterize[Show[styleMesh[mesh],
  BoxRatios -> Automatic, SphericalRegion -> True, Boxed -> False, Axes -> False,
  PlotRange -> All, ViewPoint -> view, ViewVertical -> {0, 0, 1},
  Lighting -> "Neutral", Background -> White, ImageSize -> size],
  "Image", ImageResolution -> 144];

root = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
arguments = Rest[$ScriptCommandLine];
outputRoot = If[arguments === {}, FileNameJoin[{root, "artifacts"}],
  ExpandFileName[First[arguments]]];
If[!DirectoryQ[outputRoot], CreateDirectory[outputRoot,
  CreateIntermediateDirectories -> True]];
Get[FileNameJoin[{root, "wolfram", "CPHeveaReducedSphere.wl"}]];
result = MyFun`CPHeveaReducedSphere["OutputDirectory" -> outputRoot,
  "TimeLimit" -> 3600];
If[!AssociationQ[result], fail[2, "CPHeveaReducedSphere 未返回 Association。"]];
If[!TrueQ[result["Succeeded"]], Print[result["Log"]];
  fail[3, "native 运行或数值验收失败。"]];
finalFile = SelectFirst[result["Files"], StringContainsQ[#, "stage=3"] &&
  StringContainsQ[#, "sampled_reduced_sphere"] &, Missing["NotFound"]];
If[MissingQ[finalFile], fail[4, "找不到最终 sampled VTK。"]];
mesh = Quiet@Check[Import[finalFile], $Failed];
If[Head[mesh] =!= Graphics3D, fail[5, "MMA 未能导入最终 VTK。"]];
image = renderMesh[mesh, {2.1, -2.5, 1.6}];
imageFile = FileNameJoin[{outputRoot, "reduced-sphere-v5.png"}];
Export[imageFile, image];
Print["reproduction=pass image=" <> imageFile <> " manifest=" <>
  ToString[result["ManifestFile"]]];
