// ForestAGB-GEE Green Tree Reservoir 2021 demonstration configuration.
// Copy these values into the CONFIG block in src/ForestAGB_GEE.js.
// IMPORTANT: the rectangle below is an approximate extent derived from coordinates
// reported in Zurqani (2025), not the exact original GTR boundary.

var GTR_DEMO_CONFIG = {
  roi: ee.Geometry.Rectangle([-91.6333, 34.3000, -91.5833, 34.3333], null, false),
  startDate: '2021-01-01',
  endDate: '2022-01-01',
  cloudScoreThreshold: 0.50,
  s2Composite: 'median',
  s1Composite: 'mean',
  demSource: 'USGS_3DEP_1m',
  useCanopyHeight: true,
  gediRelativeSeMax: 0.50,
  gediMaxSlopeDeg: 16.699244234,
  trainFraction: 0.70,
  randomSeed: 42,
  kFolds: 10,
  modelName: 'RF'
};

// The paper reports a 30% slope threshold. Earth Engine's ee.Terrain.slope()
// returns degrees, so the reference profile uses atan(0.30) = 16.699 degrees.
// To reproduce the supplied legacy script logic instead, use 30.0 degrees.
