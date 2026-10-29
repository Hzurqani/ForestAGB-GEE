/**
 * ForestAGB-GEE v1.0.0
 * Reproducible multi-source forest aboveground biomass density modeling in GEE.
 *
 * Scientific companion:
 * Zurqani, H.A. (2025), Ecological Informatics 86, 103052.
 * DOI: 10.1016/j.ecoinf.2025.103052
 *
 * License: MIT. See LICENSE.txt and NOTICE.md.
 */

// -----------------------------------------------------------------------------
// 1. USER CONFIGURATION
// -----------------------------------------------------------------------------
var CONFIG = {
  // Approximate published GTR geographic extent for demonstration only.
  // Replace with your exact study boundary for analysis/reproduction.
  roi:  ee.Geometry.Polygon(
        [[[-91.63473534912838, 34.33618714567295],
          [-91.63473534912838, 34.29989179794097],
          [-91.57946038574947, 34.29989179794097],
          [-91.57946038574947, 34.33618714567295]]], null, false),

  startDate: '2021-01-01',
  endDate: '2022-01-01',

  // Optical preprocessing.
  cloudScoreBand: 'cs',
  cloudScoreThreshold: 0.50,
  s2Composite: 'median',             // 'median' or 'mean'

  // SAR preprocessing.
  s1Composite: 'mean',               // 'median' or 'mean'
  s1OrbitPass:  null ,               // e.g., 'ASCENDING', 'DESCENDING', or null
  sarIndexMode: 'linear_power', 
  includeS1Angle: false,

  // Terrain.
  demSource: 'USGS_3DEP_1m',

  // Canopy height.
  useCanopyHeight: true,
  canopyHeightAsset: 'projects/meta-forest-monitoring-okw37/assets/CanopyHeight',
  canopyMinMeters: 1,
  maskBelowCanopyMin: true,

  // GEDI quality control. ee.Terrain.slope() is in degrees.
  gediRelativeSeMax: 0.50,
  gediMaxSlopeDeg: 16.699244234,     // 30% grade expressed in degrees

  // Predictor profile.
  // Paper-consistent Sentinel-2 bands exclude B1 aerosol and B9 water vapor.
  s2Bands: ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'],
  useS2Indices: true,
  useS1: true,
  useTerrain: true,

  // Sampling/modeling.
  sampleScale: 25,
  modelScale: 20,
  trainFraction: 0.70,
  randomSeed: 42,
  kFolds: 10,
  cvRepeats: 1,
  modelName: 'RF',                    // 'RF', 'CART', 'SVM', 'GTB'

  rf: {
    numberOfTrees: 120,
    variablesPerSplit: 7,
    minLeafPopulation: 5,
    bagFraction: 0.60,
    maxNodes: 50,
    seed: 1
  },
  cart: {
    maxNodes: 40,
    minLeafPopulation: 10
  },
  svm: {
    svmType: 'EPSILON_SVR',
    kernelType: 'RBF',
    gamma: 0.1,
    cost: 10,
    lossEpsilon: 0.1,
    shrinking: true
  },
  gtb: {
    numberOfTrees: 150,
    shrinkage: 0.05,
    samplingRate: 0.70,
    maxNodes: 40,
    loss: 'Huber',
    seed: 42
  },

  // Export options.
  exportFolder: 'ForestAGB_GEE',
  exportPrefix: 'ForestAGB_GEE_2021',
  exportScale: 20,
  maxPixels: 1e13
};

var roi = CONFIG.roi;
var startDate = ee.Date(CONFIG.startDate);
var endDate = ee.Date(CONFIG.endDate);
Map.centerObject(roi, 12);
Map.addLayer(roi, {color: 'ffffff'}, 'ROI', true);

// -----------------------------------------------------------------------------
// 2. SENTINEL-2 PREPROCESSING AND SOURCE-CONSISTENT SPECTRAL TRANSFORMS
// -----------------------------------------------------------------------------
function scaleS2(image) {
  return image.multiply(0.0001)
    .copyProperties(image, ['system:time_start', 'system:index']);
}

function maskS2CloudScore(image) {
  var mask = image.select(CONFIG.cloudScoreBand).gte(CONFIG.cloudScoreThreshold);
  return image.updateMask(mask);
}

function addS2Indices(image) {
  var B2 = image.select('B2');
  var B3 = image.select('B3');
  var B4 = image.select('B4');
  var B5 = image.select('B5');
  var B8 = image.select('B8');
  var B11 = image.select('B11');

  var s2RVI = image.expression('NIR / RED', {NIR: B8, RED: B4}).rename('S2_RVI');
  var ndvi = image.expression('(NIR - RED) / (NIR + RED)', {NIR: B8, RED: B4}).rename('NDVI');
  var ndre = image.expression('(NIR - RE) / (NIR + RE)', {NIR: B8, RE: B5}).rename('NDRE');
  var gndvi = image.expression('(NIR - G) / (NIR + G)', {NIR: B8, G: B3}).rename('GNDVI');
  var datt = image.expression('RED / (G * RE)', {RED: B4, G: B3, RE: B5}).rename('DATT');
  var npci = image.expression('(RED - BLUE) / (RED + BLUE)', {RED: B4, BLUE: B2}).rename('NPCI');
  var mtci = image.expression('(NIR - RE) / (NIR - RED)', {NIR: B8, RE: B5, RED: B4}).rename('MTCI');
  var sipi = image.expression('(NIR - BLUE) / (NIR - RED)', {NIR: B8, BLUE: B2, RED: B4}).rename('SIPI');
  var psri = image.expression('(RED - BLUE) / NIR', {RED: B4, BLUE: B2, NIR: B8}).rename('PSRI');
  var evi = image.expression('2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)',
    {NIR: B8, RED: B4, BLUE: B2}).rename('EVI');
  var msr = image.expression('NIR / RED - 1 / sqrt(NIR / RED + 1)', {NIR: B8, RED: B4}).rename('MSR');
  var tcari = image.expression('3 * ((RE - RED) - 0.2 * (RE - G) * (RE / RED))',
    {RE: B5, RED: B4, G: B3}).rename('TCARI');
  var ndbi = image.normalizedDifference(['B11','B8']).rename('S2_NDBI');
  var gli = image.expression('(2 * G - RED - BLUE) / (2 * G + RED + BLUE)',
    {G: B3, RED: B4, BLUE: B2}).rename('GLI');
  // Reflectance is already scaled to 0-1; do NOT divide by 10000 again.
  var avi = image.expression('NIR * pow((1 - RED) * (NIR - RED), 1.0/3.0)',
    {NIR: B8, RED: B4}).float().rename('AVI');

  return image.addBands([
    s2RVI, ndvi, ndre, gndvi, datt, npci, mtci, sipi, psri,
    evi, msr, tcari, ndbi, gli, avi
  ]);
}

function reduceCollection(collection, method) {
  if (method === 'mean') return collection.mean();
  return collection.median();
}

function buildS2() {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterDate(startDate, endDate)
    .filterBounds(roi);
  var s2Projection = ee.Image(s2.first()).select('B4').projection();

  var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
  var csPlusBands = csPlus.first().bandNames();
  var linked = s2.linkCollection(csPlus, csPlusBands)
    .map(maskS2CloudScore)
    .select('B.*')
    .map(scaleS2);

  if (CONFIG.useS2Indices) linked = linked.map(addS2Indices);
  return reduceCollection(linked, CONFIG.s2Composite)
    .setDefaultProjection(s2Projection);
}

// -----------------------------------------------------------------------------
// 3. SENTINEL-1
// -----------------------------------------------------------------------------
function dbToPower(db) {
  return ee.Image(10).pow(db.divide(10));
}

function addS1Indices(image) {
  var vv = image.select('VV');
  var vh = image.select('VH');
  var vvCalc = vv;
  var vhCalc = vh;
  if (CONFIG.sarIndexMode === 'linear_power') {
    vvCalc = dbToPower(vv);
    vhCalc = dbToPower(vh);
  }
  var denom = vvCalc.add(vhCalc);
  var s1Ndbi = vvCalc.subtract(vhCalc).divide(denom).rename('S1_NDBI');
  var s1Rvi = vhCalc.multiply(4).divide(denom).rename('S1_RVI');
  return image.addBands([s1Ndbi, s1Rvi]);
}

function buildS1() {
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
  if (CONFIG.s1OrbitPass) {
    s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', CONFIG.s1OrbitPass));
  }
  s1 = s1.map(addS1Indices);
  return reduceCollection(s1, CONFIG.s1Composite);
}

// -----------------------------------------------------------------------------
// 4. TERRAIN AND CANOPY HEIGHT
// -----------------------------------------------------------------------------
function buildTerrain() {
  var elevation;
  if (CONFIG.demSource === 'USGS_3DEP_1m') {
    var dep = ee.ImageCollection('USGS/3DEP/1m').filterBounds(roi).select('elevation');
    var depProj = ee.Image(dep.first()).select('elevation').projection();
    elevation = dep.mosaic().rename('dem').setDefaultProjection(depProj);
  } else {
    // Compatibility option for the supplied research script. The Earth Engine
    // catalog currently labels COPERNICUS/DEM/GLO30 as deprecated; use 3DEP for
    // the U.S. reference application or replace with an appropriate current DEM.
    var glo = ee.ImageCollection('COPERNICUS/DEM/GLO30').filterBounds(roi).select('DEM');
    var proj = ee.Image(glo.first()).select(0).projection();
    elevation = glo.mosaic().rename('dem').setDefaultProjection(proj);
  }
  var slope = ee.Terrain.slope(elevation).rename('slope_deg');
  return elevation.addBands(slope);
}

function buildCanopyHeight() {
  var ch = ee.ImageCollection(CONFIG.canopyHeightAsset)
    .filterBounds(roi)
    .mosaic()
    .select([0])
    .rename('canopy_height');
  if (CONFIG.maskBelowCanopyMin) {
    ch = ch.updateMask(ch.gte(CONFIG.canopyMinMeters));
  }
  return ch;
}

// -----------------------------------------------------------------------------
// 5. GEDI L4A QUALITY CONTROL
// -----------------------------------------------------------------------------
function buildGedi(slopeImage) {
  var gedi = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
    .filterDate(startDate, endDate)
    .filterBounds(roi);

  var qualityMask = function(image) {
    return image.updateMask(image.select('l4_quality_flag').eq(1))
      .updateMask(image.select('degrade_flag').eq(0));
  };
  var errorMask = function(image) {
    var relativeSe = image.select('agbd_se').divide(image.select('agbd'));
    return image.updateMask(relativeSe.lte(CONFIG.gediRelativeSeMax));
  };
  var slopeMask = function(image) {
    return image.updateMask(slopeImage.lt(CONFIG.gediMaxSlopeDeg));
  };

  var proj = ee.Image(gedi.first()).select('agbd').projection();
  return gedi.map(qualityMask).map(errorMask).map(slopeMask)
    .mosaic().select('agbd').rename('agbd').setDefaultProjection(proj);
}

// -----------------------------------------------------------------------------
// 6. MODEL CONSTRUCTORS, METRICS, AND CROSS-VALIDATION
// -----------------------------------------------------------------------------
function makeRegressor(name) {
  if (name === 'CART') {
    return ee.Classifier.smileCart({
      maxNodes: CONFIG.cart.maxNodes,
      minLeafPopulation: CONFIG.cart.minLeafPopulation
    }).setOutputMode('REGRESSION');
  }
  if (name === 'SVM') {
    return ee.Classifier.libsvm({
      svmType: CONFIG.svm.svmType,
      kernelType: CONFIG.svm.kernelType,
      gamma: CONFIG.svm.gamma,
      cost: CONFIG.svm.cost,
      lossEpsilon: CONFIG.svm.lossEpsilon,
      shrinking: CONFIG.svm.shrinking
    }).setOutputMode('REGRESSION');
  }
  if (name === 'GTB') {
    return ee.Classifier.smileGradientTreeBoost({
      numberOfTrees: CONFIG.gtb.numberOfTrees,
      shrinkage: CONFIG.gtb.shrinkage,
      samplingRate: CONFIG.gtb.samplingRate,
      maxNodes: CONFIG.gtb.maxNodes,
      loss: CONFIG.gtb.loss,
      seed: CONFIG.gtb.seed
    }).setOutputMode('REGRESSION');
  }
  return ee.Classifier.smileRandomForest({
    numberOfTrees: CONFIG.rf.numberOfTrees,
    variablesPerSplit: CONFIG.rf.variablesPerSplit,
    minLeafPopulation: CONFIG.rf.minLeafPopulation,
    bagFraction: CONFIG.rf.bagFraction,
    maxNodes: CONFIG.rf.maxNodes,
    seed: CONFIG.rf.seed
  }).setOutputMode('REGRESSION');
}

function scoreMetrics(scored, observedName, predictedName) {
  var withErrors = scored.map(function(f) {
    var o = ee.Number(f.get(observedName));
    var p = ee.Number(f.get(predictedName));
    var e = p.subtract(o);
    return f.set({sqerr: e.pow(2), abserr: e.abs()});
  });
  var meanObs = ee.Number(withErrors.aggregate_mean(observedName));
  var withDev = withErrors.map(function(f) {
    return f.set('sqdev', ee.Number(f.get(observedName)).subtract(meanObs).pow(2));
  });
  var sse = ee.Number(withDev.aggregate_sum('sqerr'));
  var sst = ee.Number(withDev.aggregate_sum('sqdev'));
  var r2 = ee.Number(1).subtract(sse.divide(sst));
  var rmse = ee.Number(withDev.aggregate_mean('sqerr')).sqrt();
  var mae = ee.Number(withDev.aggregate_mean('abserr'));
  return ee.Dictionary({R2: r2, RMSE: rmse, MAE: mae, n: withDev.size()});
}

function kFoldCV(training, predictors, target) {
  var all = ee.FeatureCollection([]);
  for (var repeat = 0; repeat < CONFIG.cvRepeats; repeat++) {
    var seed = CONFIG.randomSeed + repeat * 1000;
    var folded = training.randomColumn('cv_random', seed).map(function(f) {
      var fold = ee.Number(f.get('cv_random')).multiply(CONFIG.kFolds).floor().int();
      return f.set('cv_fold', fold);
    });
    var foldFeatures = ee.List.sequence(0, CONFIG.kFolds - 1).map(function(i) {
      i = ee.Number(i);
      var trainFold = folded.filter(ee.Filter.neq('cv_fold', i));
      var validFold = folded.filter(ee.Filter.eq('cv_fold', i));
      var model = makeRegressor(CONFIG.modelName).train({
        features: trainFold,
        classProperty: target,
        inputProperties: predictors
      });
      var scored = validFold.classify(model, 'prediction');
      var m = scoreMetrics(scored, target, 'prediction');
      return ee.Feature(null, m).set({fold: i, repeat: repeat + 1, model: CONFIG.modelName});
    });
    all = all.merge(ee.FeatureCollection(foldFeatures));
  }
  return all;
}

function importanceToTable(model) {
  var explanation = ee.Dictionary(model.explain());
  var importance = ee.Dictionary(ee.Algorithms.If(
    explanation.contains('importance'), explanation.get('importance'), ee.Dictionary({})
  ));
  return ee.FeatureCollection(importance.keys().map(function(k) {
    return ee.Feature(null, {variable: k, importance: importance.get(k)});
  }));
}

// -----------------------------------------------------------------------------
// 7. BUILD MULTI-SOURCE STACK AND REFERENCE SAMPLES
// -----------------------------------------------------------------------------
var s2Composite = buildS2();
var terrain = buildTerrain();
var slopeDeg = terrain.select('slope_deg');
var gediAgbd = buildGedi(slopeDeg);

var predictorsImage = s2Composite.select(CONFIG.s2Bands);
if (CONFIG.useS2Indices) {
  predictorsImage = predictorsImage.addBands(s2Composite.select([
    'S2_RVI','NDVI','NDRE','GNDVI','DATT','NPCI','MTCI','SIPI','PSRI',
    'EVI','MSR','TCARI','S2_NDBI','GLI','AVI'
  ]));
}
if (CONFIG.useS1) {
  var s1Composite = buildS1();
  var s1Bands = ['VV','VH','S1_NDBI','S1_RVI'];
  if (CONFIG.includeS1Angle) s1Bands.push('angle');
  predictorsImage = predictorsImage.addBands(s1Composite.select(s1Bands));
}
if (CONFIG.useTerrain) predictorsImage = predictorsImage.addBands(terrain);
if (CONFIG.useCanopyHeight) predictorsImage = predictorsImage.addBands(buildCanopyHeight());

var predictorNames = predictorsImage.bandNames();
var trainingTable = predictorsImage.addBands(gediAgbd).sample({
  region: roi,
  scale: CONFIG.sampleScale,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('ForestAGB-GEE version', '1.0.0');
print('Predictor bands', predictorNames);
print('Number of retained GEDI/predictor samples', trainingTable.size());

// -----------------------------------------------------------------------------
// 8. TRAIN / HOLDOUT EVALUATION / K-FOLD CV
// -----------------------------------------------------------------------------
var splitTable = trainingTable.randomColumn('split_random', CONFIG.randomSeed);
var train = splitTable.filter(ee.Filter.lt('split_random', CONFIG.trainFraction));
var test = splitTable.filter(ee.Filter.gte('split_random', CONFIG.trainFraction));

var trainedModel = makeRegressor(CONFIG.modelName).train({
  features: train,
  classProperty: 'agbd',
  inputProperties: predictorNames
});

var holdoutPredictions = test.classify(trainedModel, 'prediction');
var holdoutMetrics = scoreMetrics(holdoutPredictions, 'agbd', 'prediction');
var cvMetrics = kFoldCV(train, predictorNames, 'agbd');
var importanceTable = importanceToTable(trainedModel);

print('Training n', train.size());
print('Holdout n', test.size());
print('Holdout metrics', holdoutMetrics);
print('Cross-validation metrics', cvMetrics);
print('Model explanation', trainedModel.explain());

// -----------------------------------------------------------------------------
// 9. SPATIAL PREDICTION AND VISUALIZATION
// -----------------------------------------------------------------------------
var predictedAgbd = predictorsImage.classify(trainedModel, 'agbd_predicted');

Map.addLayer(s2Composite.clip(roi), {bands:['B4','B3','B2'], min:0, max:0.3, gamma:1.2}, 'Sentinel-2', false);
Map.addLayer(gediAgbd.clip(roi), {min:0, max:200, palette:['edf8fb','b2e2e2','66c2a4','2ca25f','006d2c']}, 'Filtered GEDI AGBD', true);
Map.addLayer(predictedAgbd.clip(roi), {min:0, max:200, palette:['ffffcc','a1dab4','41b6c4','2c7fb8','253494']}, 'Predicted AGBD', true);

// -----------------------------------------------------------------------------
// 10. EXPORTS
// -----------------------------------------------------------------------------
var metricsFeature = ee.Feature(null, holdoutMetrics)
  .set({model: CONFIG.modelName, software_version: '1.0.0'});

Export.image.toDrive({
  image: predictedAgbd.clip(roi),
  description: CONFIG.exportPrefix + '_Predicted_AGBD',
  folder: CONFIG.exportFolder,
  fileNamePrefix: CONFIG.exportPrefix + '_predicted_agbd',
  region: roi,
  scale: CONFIG.exportScale,
  maxPixels: CONFIG.maxPixels
});

Export.table.toDrive({
  collection: holdoutPredictions,
  description: CONFIG.exportPrefix + '_Holdout_Predictions',
  folder: CONFIG.exportFolder,
  fileNamePrefix: CONFIG.exportPrefix + '_holdout_predictions',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: ee.FeatureCollection([metricsFeature]),
  description: CONFIG.exportPrefix + '_Holdout_Metrics',
  folder: CONFIG.exportFolder,
  fileNamePrefix: CONFIG.exportPrefix + '_holdout_metrics',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: cvMetrics,
  description: CONFIG.exportPrefix + '_CV_Metrics',
  folder: CONFIG.exportFolder,
  fileNamePrefix: CONFIG.exportPrefix + '_cv_metrics',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: importanceTable,
  description: CONFIG.exportPrefix + '_Variable_Importance',
  folder: CONFIG.exportFolder,
  fileNamePrefix: CONFIG.exportPrefix + '_variable_importance',
  fileFormat: 'CSV'
});
