# ForestAGB-GEE

**ForestAGB-GEE** is a reproducible Google Earth Engine (GEE) workflow for multi-source forest aboveground biomass density (AGBD) modeling using GEDI L4A, Sentinel-2, Sentinel-1, terrain, canopy-height information, and machine-learning regression.

This repository is the companion to:

> Zurqani, H.A. (2025). *A multi-source approach combining GEDI LiDAR, satellite data, and machine learning algorithms for estimating forest aboveground biomass on Google Earth Engine platform*. Ecological Informatics, 86, 103052. https://doi.org/10.1016/j.ecoinf.2025.103052

## What is included

- `src/ForestAGB_GEE.js` — standalone v1.0 workflow intended for public use.
- `examples/GTR_2021_demo.js` — a compact example configuration for the 2021 Green Tree Reservoir reference setting. The included rectangle is an **approximate published geographic extent**, not the exact original reservoir boundary.
- `docs/USER_GUIDE.md` — setup, execution, outputs, and troubleshooting.
- `results/reference_benchmark_metrics.csv` — benchmark metrics reported in the 2025 companion paper.
- `figures/` — figures used in the SoftwareX manuscript.
- `LICENSE.txt`, and `NOTICE.md`.

## Quick start

1. Sign in to the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Open `src/ForestAGB_GEE.js` and paste it into a new Code Editor script.
3. Replace `CONFIG.roi` with your own `ee.Geometry` or `ee.FeatureCollection`, or retain the approximate demonstration extent.
4. Review `CONFIG` at the top of the script, especially dates, DEM source, GEDI filters, model, and export settings.
5. Run the script. Review Console outputs and map layers.
6. Start the generated export tasks from the Earth Engine **Tasks** tab.

## Requirements

- Google Earth Engine account/project with Code Editor access.
- Modern web browser.
- Access to the public GEE datasets used by the selected configuration.
- If canopy height is enabled, access to `projects/meta-forest-monitoring-okw37/assets/CanopyHeight` is expected; otherwise set `useCanopyHeight: false`.

No local compilation is required.

## Citation

Please cite both the softwareX release and the scientific companion paper.

## License

MIT License. See `LICENSE.txt` and `NOTICE.md`.
