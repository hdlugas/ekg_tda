# ekg_tda
This code uses tools from persistent homology to identify optimal representative 1-cycles as P,Q,S, and T-waves of an electrocardiogram (ECG) signal. These cycles are then used to measure the PR-interval, QT-interval, ST-segment, QRS-duration, P-wave duration, and T-wave duration based off of the upper and lower time axis bounds of these representative cycles. Briefly, P,Q,S, and T waves are characterized using dimension one homological features with a persistence and centroid of representative cycle within certain ranges depending on the specific waveform, i.e. a P,Q,S, or T-wave. For better understanding and intuition about what the topological invariants of the data are and how they're used to measure various intervals of interest, see "ekg_feature_extraction_intro.docx".

Here is a brief description of the code:

ekg_sim.py is used to analyze simulated ECG signals

ekg_real_data.py is used to analyze real ECG signals

cycles.py contains functions to compute the centroid of boundary points of optimal cycles, compute the onset and offset of the intervals of interest, and draw optimal 1-cycles identified as P,Q,S, and T-waves

intervals.py contains functions to measure intervals of interest

processing.py contains functions to process an ECG signal prior to computing its persistent homology

ekg_example_persistent_homology.R generates teh figures shown in the document 'ekg_feature_extraction_intro.docx'


Here is an image of a simulated ECG signal with area-optimal 1-cycles with certain properties depending on their birth filtration, persistence, and centroid identified as P,Q,S, and T-waves:
![cycle_example](https://user-images.githubusercontent.com/73852653/147366648-d563e3a3-68db-4663-a7d3-add220ce05e1.png)
