ClearAll[fail, styleMesh, renderMesh, binaryMask, largestComponent,
  componentCount, maskAspect];

fail[code_Integer, message_String] := (Print["ERROR: " <> message]; Exit[code]);

styleMesh[mesh_] := mesh /. GraphicsComplex[points_, primitives_, options___] :>
  GraphicsComplex[points, {Directive[RGBColor[0.78, 0.64, 0.34],
      Specularity[White, 24]], EdgeForm[], primitives}, options];

renderMesh[mesh_, view_, size_: 800] := Rasterize[
  Show[styleMesh[mesh], BoxRatios -> Automatic, SphericalRegion -> True,
    Boxed -> False, Axes -> False, PlotRange -> All, ViewPoint -> view,
    ViewVertical -> {0, 0, 1}, Lighting -> "Neutral", Background -> White,
    ImageSize -> size], "Image", ImageResolution -> 144];

binaryMask[mesh_, view_] := Module[{silhouette},
  silhouette = mesh /. GraphicsComplex[points_, primitives_, options___] :>
    GraphicsComplex[points, {Black, EdgeForm[], primitives}, options];
  Binarize[ColorNegate@ColorConvert[Rasterize[
    Show[silhouette, BoxRatios -> Automatic, SphericalRegion -> True,
      Boxed -> False, Axes -> False, PlotRange -> All, ViewPoint -> view,
      ViewVertical -> {0, 0, 1}, Lighting -> {{"Ambient", White}},
      Background -> White, ImageSize -> 500], "Image"], "Grayscale"], 0.5]
];

largestComponent[image_] := Module[{labels, counts, label},
  labels = MorphologicalComponents[image];
  counts = KeyDrop[Counts@Flatten[labels], 0];
  If[Length[counts] == 0, Return[ConstantImage[0, ImageDimensions[image]]]];
  label = First@First@MaximalBy[Normal[counts], Last];
  Image[Map[Boole[# == label] &, labels, {2}], "Bit"]
];

componentCount[image_, minimumArea_: 50] := Module[{labels, counts},
  labels = MorphologicalComponents[image];
  counts = KeyDrop[Counts@Flatten[labels], 0];
  Count[Values[counts], area_ /; area >= minimumArea]
];

maskAspect[image_] := Module[{positions, ranges},
  positions = PixelValuePositions[image, 1];
  If[Length[positions] == 0, Return[Indeterminate]];
  ranges = MinMax /@ Transpose[positions];
  N[(Subtract @@ Reverse[ranges[[1]]] + 1)/(Subtract @@ Reverse[ranges[[2]]] + 1)]
];

root = DirectoryName[DirectoryName[$InputFileName]];
paperDirectoryOverride = Environment["HEVEA_FIGURE9_PAPER_DIRECTORY"];
renderSubdirectoryOverride = Environment["HEVEA_FIGURE9_RENDER_SUBDIRECTORY"];
diagnosticMode = Environment["HEVEA_FIGURE9_DIAGNOSTIC"] === "1";
skipConversion = Environment["HEVEA_FIGURE9_SKIP_CONVERSION"] === "1";
referenceMeshOverride = Environment["HEVEA_FIGURE9_REFERENCE_MESH"];
paperDirectory = If[StringQ[paperDirectoryOverride] &&
  StringLength[paperDirectoryOverride] > 0, paperDirectoryOverride,
  FileNameJoin[{root, "artifacts", "paper-grid"}]];
renderSubdirectory = If[StringQ[renderSubdirectoryOverride] &&
  StringLength[renderSubdirectoryOverride] > 0, renderSubdirectoryOverride,
  "figure9-mma-v5"];
renderDirectory = FileNameJoin[{paperDirectory, renderSubdirectory}];
If[!DirectoryQ[renderDirectory],
  CreateDirectory[renderDirectory, CreateIntermediateDirectories -> True]];

sourceManifestFile = FileNameJoin[{paperDirectory, "manifest.json"}];
converter = FileNameJoin[{root, "tools", "export_paper_meshes"}];
conversionArguments = If[diagnosticMode,
  {converter, "--diagnostic", sourceManifestFile, renderDirectory},
  {converter, sourceManifestFile, renderDirectory}];
If[!skipConversion,
  conversion = RunProcess[conversionArguments];
  If[conversion["ExitCode"] =!= 0,
    Print[conversion["StandardOutput"]]; Print[conversion["StandardError"]]; Exit[2]]
];

meshManifestFile = FileNameJoin[{renderDirectory, "figure9-mesh-manifest.json"}];
If[!FileExistsQ[meshManifestFile], fail[3, "C++ 导出器没有生成 mesh manifest。"]];
meshManifest = Quiet@Check[Import[meshManifestFile, "RawJSON"], $Failed];
If[!AssociationQ[meshManifest], fail[3, "无法读取 mesh manifest。"]];
envelopeTool = FileNameJoin[{root, "tools", "corrugation_envelope"}];
envelopeFile = FileNameJoin[{renderDirectory, "stage3-mode21-envelope.csv"}];
referenceEnvelopeFile = FileNameJoin[{root, "tests", "fixtures",
  "reference_stage1_mode21_envelope.csv"}];
envelopeExtraction = RunProcess[{envelopeTool, "--binary",
  meshManifest["full_resolution_final"], envelopeFile, "2", "20"}];
If[envelopeExtraction["ExitCode"] =!= 0,
  Print[envelopeExtraction["StandardOutput"]];
  Print[envelopeExtraction["StandardError"]];
  fail[4, "无法提取最终曲面的纬向 corrugation 包络。"]];
envelopeComparison = RunProcess[{envelopeTool, "--compare", envelopeFile,
  referenceEnvelopeFile, "1.10", ".10", ".15"}];
If[envelopeComparison["ExitCode"] =!= 0,
  Print[envelopeComparison["StandardOutput"]];
  Print[envelopeComparison["StandardError"]];
  fail[4, "中纬度至高纬度的 mode-21 褶皱包络与作者 WRL 不符。"]];
expectedRidges = meshManifest["ridges_expected"];
countedRidges = meshManifest["ridges_counted"];
directionRatios = meshManifest["direction_energy_ratios"];
If[countedRidges =!= expectedRidges || !AllTrue[directionRatios, # > 1 &],
  fail[4, "真实阶段差分的 ridge 数或方向检测失败。"]];

globalFile = FileNameJoin[{renderDirectory, "sampled-reduced-sphere-stage3.vtk"}];
detailFiles = FileNameJoin[{renderDirectory, #}] & /@
  {"detail-stage1-ridges21.vtk", "detail-stage2-ridges142.vtk",
   "detail-stage3-ridges997.vtk"};
If[!FileExistsQ[globalFile] || !AllTrue[detailFiles, FileExistsQ],
  fail[5, "全局网格或局部 detail patch 缺失。"]];

globalMesh = Quiet@Check[Import[globalFile], $Failed];
If[Head[globalMesh] =!= Graphics3D, fail[6, "MMA 未能把最终 VTK 导入为 Graphics3D。"]];
complexes = Cases[globalMesh,
  GraphicsComplex[points_, primitives_, ___] :> {points, primitives}, Infinity];
If[Length[complexes] =!= 1, fail[6, "最终 VTK 不是单一连通 GraphicsComplex。"]];
points = complexes[[1, 1]];
faces = FirstCase[complexes[[1, 2]], Polygon[index_] :> index, $Failed, Infinity];
If[!MatrixQ[points, NumericQ] || !ListQ[faces], fail[6, "导入网格缺少数值顶点或面。"]];
importedBounds = MinMax /@ Transpose[points];
manifestBounds = meshManifest["global_bounds"];
If[Max[Abs[Flatten[N[importedBounds] - N[manifestBounds]]]] > 10^-8,
  fail[7, "MMA 导入坐标范围与导出 manifest 不一致。"]];

topology = meshManifest["topology"];
faceSizes = Counts[Length /@ faces];
closedSphere = Length[points] == topology["point_count"] &&
  Length[faces] == topology["cell_count"] &&
  Lookup[faceSizes, 3, 0] == 2 First[meshManifest["global_grid"]] &&
  Lookup[faceSizes, 4, 0] == topology["cell_count"] -
    2 First[meshManifest["global_grid"]] && TrueQ[topology["closed"]] &&
  topology["connected_components"] == 1 && topology["euler_characteristic"] == 2;
If[!TrueQ[closedSphere], fail[8, "MMA 导入网格的结构计数不是闭合单球面。"]];

viewPoint = {2.1, -2.5, 1.6};
globalImage = renderMesh[globalMesh, viewPoint, 900];
globalPath = FileNameJoin[{renderDirectory, "global-fixed-camera.png"}];
Export[globalPath, globalImage];

globalMask = largestComponent@binaryMask[globalMesh, viewPoint];
connectedComponents = componentCount[globalMask];
holeMask = DeleteSmallComponents[
  Binarize@ImageSubtract[FillingTransform[globalMask], globalMask], 50];
holeCount = componentCount[holeMask];
globalAspect = maskAspect[globalMask];

If[StringQ[referenceMeshOverride] && StringLength[referenceMeshOverride] > 0,
  referenceMesh = Quiet@Check[Import[referenceMeshOverride], $Failed];
  If[Head[referenceMesh] =!= Graphics3D, fail[9, "无法导入作者 WRL 的诊断 VTK。"]];
  referenceCrop = renderMesh[referenceMesh, viewPoint, 900];
  referenceMask = largestComponent@binaryMask[referenceMesh, viewPoint];
  referenceRole = "author-wrl-fixed-camera",
  reference = Import[FileNameJoin[{root, "tests", "fixtures", "figure9-page-17.png"}]];
  referenceCrop = ImageTake[reference, {210, 620}, {220, 650}];
  referenceRawMask = DeleteSmallComponents[
    Binarize[ColorNegate@ColorConvert[referenceCrop, "Grayscale"], 0.08], 100];
  referenceMask = largestComponent[referenceRawMask];
  referenceRole = "paper-crop-only"
];
referenceAspect = maskAspect[referenceMask];
aspectRelativeError = Abs[globalAspect/referenceAspect - 1];
If[connectedComponents =!= 1 || holeCount =!= 0 ||
   !NumericQ[aspectRelativeError] || aspectRelativeError >= 0.10,
  fail[9, "固定相机轮廓的连通性、内部孔洞或宽高比验收失败。"]];

comparisonPath = FileNameJoin[{renderDirectory, "figure9-side-by-side.png"}];
Export[comparisonPath, GraphicsRow[{referenceCrop, globalImage},
  Spacings -> 10, ImageSize -> 1600]];

detailPaths = MapIndexed[Function[{file, index}, Module[{mesh, image, path},
  mesh = Import[file];
  image = renderMesh[mesh, viewPoint, 700];
  path = FileNameJoin[{renderDirectory,
    "detail-stage" <> ToString[First[index]] <> ".png"}];
  Export[path, image]; Clear[mesh, image]; path]], detailFiles];

globalGrid = meshManifest["global_grid"];
capRows = meshManifest["cap_rows"];
gridStart = 2 + (capRows - 1) First[globalGrid];
southSeam = points[[gridStart ;; gridStart + First[globalGrid] - 1]];
northStart = gridStart + (Last[globalGrid] - 1) First[globalGrid];
northSeam = points[[northStart ;; northStart + First[globalGrid] - 1]];
seamOverlay = Graphics3D[{Directive[Red, Opacity[0.9]],
  Tube[Append[southSeam, First[southSeam]], 0.0025],
  Tube[Append[northSeam, First[northSeam]], 0.0025]}];
seamImage = Rasterize[Show[styleMesh[globalMesh], seamOverlay,
  BoxRatios -> Automatic, SphericalRegion -> True, Boxed -> False, Axes -> False,
  PlotRange -> All, ViewPoint -> viewPoint, ViewVertical -> {0, 0, 1},
  Lighting -> "Neutral", Background -> White, ImageSize -> 800],
  "Image", ImageResolution -> 144];
seamPath = FileNameJoin[{renderDirectory, "cap-seams.png"}];
Export[seamPath, seamImage];

ballOverlay = Graphics3D[{Directive[RGBColor[0.1, 0.35, 0.9], Opacity[0.12]],
  Sphere[{0, 0, 0}, 0.52]}];
ballImage = Rasterize[Show[styleMesh[globalMesh], ballOverlay,
  BoxRatios -> Automatic, SphericalRegion -> True, Boxed -> False, Axes -> False,
  PlotRange -> All, ViewPoint -> viewPoint, ViewVertical -> {0, 0, 1},
  Lighting -> "Neutral", Background -> White, ImageSize -> 800],
  "Image", ImageResolution -> 144];
ballPath = FileNameJoin[{renderDirectory, "bounding-ball-overlay.png"}];
Export[ballPath, ballImage];

stageFiles = FileNameJoin[{renderDirectory, "global-stage" <> ToString[#] <> ".vtk"}] & /@
  Range[0, 3];
stageImages = Map[renderMesh[Import[#], viewPoint, 450] &, stageFiles];
stagesPath = FileNameJoin[{renderDirectory, "global-stages-0-through-3.png"}];
Export[stagesPath, GraphicsRow[stageImages, Spacings -> 5, ImageSize -> 1800]];

imageFiles = Join[{globalPath, comparisonPath, seamPath, ballPath, stagesPath}, detailPaths];
If[!AllTrue[imageFiles, FileExistsQ[#] && FileByteCount[#] > 10000 &],
  fail[10, "一个或多个视觉验收图缺失或为空。"]];
imageHashes = AssociationThread[imageFiles,
  IntegerString[FileHash[#, "SHA256"], 16, 64] & /@ imageFiles];

visualReport = <|
  "Figure9ReferenceRole" -> referenceRole,
  "ImportedBounds" -> importedBounds,
  "ManifestBounds" -> manifestBounds,
  "ClosedSphere" -> closedSphere,
  "ConnectedComponents" -> connectedComponents,
  "InternalHoles" -> holeCount,
  "GlobalAspect" -> globalAspect,
  "ReferenceAspect" -> referenceAspect,
  "AspectRelativeError" -> aspectRelativeError,
  "RidgesExpected" -> expectedRidges,
  "RidgesCounted" -> countedRidges,
  "DirectionEnergyRatios" -> directionRatios,
  "EnvelopeComparison" -> StringTrim[envelopeComparison["StandardOutput"]],
  "EnvelopeFile" -> envelopeFile,
  "ImageSHA256" -> imageHashes,
  "Files" -> imageFiles|>;
reportPath = FileNameJoin[{renderDirectory, "figure9-visual-report.json"}];
Export[reportPath, visualReport, "RawJSON"];

Print["figure9_render=pass mode=" <>
  If[diagnosticMode, "diagnostic", "certified"] <>
  " mma_import=True closed_sphere=" <>
  ToString[closedSphere] <> " connected_components=" <>
  ToString[connectedComponents] <> " internal_holes=" <> ToString[holeCount] <>
  " global_aspect=" <> ToString[N[globalAspect]] <> " reference_aspect=" <>
  ToString[N[referenceAspect]] <> " aspect_relative_error=" <>
  ToString[N[aspectRelativeError]] <> " bounds_delta=" <>
  ToString[Max[Abs[Flatten[N[importedBounds] - N[manifestBounds]]]]] <>
  " ridges_counted=" <> StringRiffle[ToString /@ countedRidges, ","] <>
  " direction_ratios=" <> StringRiffle[ToString /@ N[directionRatios], ","] <>
  " " <> StringTrim[envelopeComparison["StandardOutput"]] <>
  " report=" <> reportPath <> " files=" <> ToString[imageFiles, InputForm]];
