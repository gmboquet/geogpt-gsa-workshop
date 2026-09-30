# Course input provenance

Prepared by the course author from retained source files, preparation code and receipts; delivered to the hosted workspace after the dogfood source audit identified missing preparation documentation. This is supplied provenance, not evidence independently rediscovered by a language model.

## Fossils

The original M0027A occurrence table is PANGAEA doi:10.1594/PANGAEA.779542, with its metadata preserved. A nonblank species code is treated as a recorded occurrence; a blank is not established biological absence. The method selects the deepest recorded P. sicana and shallowest recorded C. dissimilis within its stated source selection. These are occurrence constraints, not directly dated local event horizons.

The calculation explicitly uses the GTS2004-labelled column of IODP Expedition 342 Methods Table T2. Ages in legacy occurrence comments must not silently replace that selected calibration. Taxonomy comes from the attributed Treatise reference, not an identification of unseen core specimens.

The twelve Sr-context rows were transcribed from Browning et al. (2013), Table 1, M27 column. Their depth is published mcd; fossil depths are mbsf. No conversion or depth tie has been established. The integrated published chronology is not an independent validation of its own fossil inputs.

## Resolution

F02-1 logs and the inline 362 export were acquired from the SEG 2017 colored-inversion tutorial at commit b38f798866dd4e1a38c78cfd45edb7b2cfd3d785. The tutorial manuscript describes a dip-steered median-filtered stacked F3 inline. This course does not possess the full original unfiltered seismic volume or a processing flow sufficient to reproduce the filtering.

The LAS header declares depth in m, RHOB in kg/m3, DT in microseconds/m and NULL=-999.25. The course retains rows with positive density and slowness. The two half-open depth windows are 1000 <= depth < 1010 m and 1010 <= depth < 1020 m. The existing loader constructs the exported section time axis at 0.004 s; this is a documented input convention in the supplied method, not newly inferred acquisition metadata.

The supplied time-depth model and checkshots are retained as separate source files. Their presence does not independently validate a seismic-to-well tie. The exercise uses measured window properties with an assumed Ricker wavelet, phase and bed geometry; its forward model is distinct from the SEG colored-inversion algorithm.

## Spatial

The course preparer used prepare.py with USGS/ScienceBase Roback_Nepal_final_files.zip (doi:10.5066/F7DZ06F9). The source archive MD5 was matched to publisher metadata. Layer names: Source20170209, MappingExtent20170209, ObscuredAreas20170209 and ImageQuality20170209. Invalid source and obscuration geometries were repaired with make_valid before spatial operations; the retained preparation receipt records 53 and 12 repairs respectively.

The prepared coordinate system is WGS84 / UTM zone 45N, EPSG:32645; x_m and y_m are cell-centre coordinates in metres. A central subset bounded from 85.08–85.92 degrees E and 27.55–28.55 degrees N was snapped to a 500 m grid. Eligible cells must be fully covered by the mapping extent, not intersect an obscured-area polygon, and have finite slope, relief and PGA. This produced 27413 retained cells from 36573 before fitting; 9160 were excluded by extent/obscuration, and none by nonfinite predictors.

The response is 1 when a 500 m cell intersects any mapped landslide-source polygon, not a landslide centroid, full runout footprint, fractional affected area or future failure probability. The source layer contains 24795 polygons. Retained presence count is 4098. Nonmapped cells are not individually verified stable slopes. The inventory concerns the earthquake sequence; the ShakeMap predictor is the mainshock product.

Terrain uses Mapzen Skadi tiles N27E085 and N28E085, an SRTM-based compilation with fills, not independently certified here as pure pre-event SRTM. Elevation was bilinearly reprojected to 100 m EPSG:32645 pixels. Slope is atan of the gradient magnitude in degrees; each 500 m cell averages 25 slope/elevation pixels, and relief is maximum minus minimum elevation over those pixels.

PGA comes from the archived USGS ShakeMap grid.xml product us20002926/atlas/1594162031303, with units percent g. Regular-grid interpolation samples it at cell centres. Image quality copies the source img_qualit value only for cells fully covered by a quality polygon; -1 means no such assignment. Class counts are 0:5742, 1:21234, -1:437. The 0/1 quality ranking has no verified codebook here, so the core exercise uses no quality filter.

## Source records

https://doi.pangaea.de/10.1594/PANGAEA.779542

https://publications.iodp.org/proceedings/342/102/102_t2.htm

https://doi.org/10.1130/GES00857.1

https://github.com/seg/tutorials-2017/tree/b38f798866dd4e1a38c78cfd45edb7b2cfd3d785/1710_Colored_inversion

https://doi.org/10.5066/F7DZ06F9

https://earthquake.usgs.gov/product/shakemap/us20002926/atlas/1594162031303/download/grid.xml

## Reproducibility identities

prepare.py: SHA-256 2f3b3a95990d9d8934371738aa7fc30d891cea79af5d46790c3030e2504b8287

workshop.py: SHA-256 e036a37fe513805a8e09631bde0913be7f53f2563c65c14dccb7cb7e22b9b5de

classroom.py: SHA-256 d0aa610b6776b3d97507afcccc92f970914dc1e8b6e1a9de80fd75432cfb4417

Input identities are in participant-inputs.json. The preparation script and original archive/receipt remain instructor evidence; only the frozen analysis inputs are needed for the workshop.
