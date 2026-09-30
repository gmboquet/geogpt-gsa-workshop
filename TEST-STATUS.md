# Public-link setup rehearsal

The GitHub-to-Research workflow completed in Safari on the international Zero2x Research service, using a new project in an existing account.

## Verified by actual execution and downloaded artifacts

- Research retrieved the public `workshop-runtime.zip`: HTTP 200, 2,438,141 bytes, SHA-256 `01ef78539d632bb2c6e2e5516613523bbe934563722d6215698a81a7163d0b4b`.
- Safe extraction produced thirteen files. All twelve manifest entries passed verification before and after execution. No packages or replacement scientific methods were installed.
- Hosted Python 3.12.3 imported numpy 2.4.6, pandas 3.0.5, matplotlib 3.11.1, scipy 1.18.0 and scikit-learn 1.9.0.
- The seismic run used 55 Hz, 30-degree phase and a 12 m target bed. The spatial run used the north holdout, 3 km buffer and no image-quality filter.
- Both runs completed through the unchanged `classroom.run` interface, produced their figures and tables, and saved research memos that Research read back.
- The project was exported through the Research browser interface. The downloaded export retained all twelve source hashes. Six output CSVs—including 6,920 held-out predictions—agreed with the separate local reference within numerical precision (largest absolute difference approximately 1.4e-14). Elapsed times were excluded from numerical equality checks.
- Recorded calculation times were 0.72 seconds for seismic and 2.55 seconds for spatial. These are calculation times, not full conversation or lesson times. The package-setup response displayed 1 minute 14 seconds.
- Both downloaded PNGs were inspected: the seismic spectrum/wedge/bias plot and the geographic holdout/prediction/performance plot contain the expected results.

## Interpretation remains part of the exercise

Successful execution does not validate every sentence of a model-written memo. The rehearsal caught unsupported statements about a measured seismic noise floor, attributing a difference to phase without a control, treating Brier score as calibration alone, and ranking undocumented image-quality codes. These require evidence review; they are not supplied answers for students to copy.

## Scope and remaining checks

This test establishes public-file retrieval, unchanged hosted execution, artifact generation, memo saving and export. It does not establish permissions for a newly registered student account, class-scale capacity or participant pacing. The fossil exercise uses the embedded CT viewer and GeoGPT Chat; it does not use the quantitative runtime package and was not re-run through Chat during this setup test.

The landing page and copy control loaded publicly. The first full-handout download was interrupted by a local connection reset; the page now explains its loading state. Complete the handout download before class and retain a saved copy.
