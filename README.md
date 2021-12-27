# ekg_tda
This code uses tools from persistent homology to identify optimal representative 1-cycles as P,Q,S, and T-waves. These cycles are then used to measure the PR-interval, QT-interval, ST-segment, QRS-duration, P-wave duration, and T-wave duration based off of the upper and lower time axis bounds of these representative cycles. Briefly, P,Q,S, and T waves are characterized using dimension one homological features with a persistence and centroid of representative cycle within certain ranges depending on the specific waveform, i.e. a P,Q,S, or T-wave. For better understanding and intuition about what the topological invariants of the data are and how they're used to measure various intervals of interest, see "ekg_feature_extraction_intro.docx".

Here is a brief description of the code:
ekg_sim.py is used to analyze simultated ECG signals
ekg_real_data.py is used to analyze real ECG signals
cycles.py contains functions 
intervals.py contains functions
processing.py contains functions

![cycle_example](https://user-images.githubusercontent.com/73852653/147366648-d563e3a3-68db-4663-a7d3-add220ce05e1.png)
