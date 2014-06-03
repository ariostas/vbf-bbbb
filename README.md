Here you can find all my code used for the vbf-bbbb analysis.

The main files are the ones called vbf-bbbb_X.C, all of which run over the vbf-bbbb and ttbar samples. All of these perform the analysis in a different way.

vbf-bbbb_0.C: selects 4 bjets and 2 light jets.

vbf-bbbb_1.C: selects 3 bjets and 3 light jets.

vbf-bbbb_2.C: selects 6 jets (whatever they are) and chooses the light jets to be the outermost on opposite sides.

vbf-bbbb_3.C: selects 6 jets (whatever they are) and chooses the light jets to be the two jets with highest pt.

The output histograms of these macros can be found at the Histograms folder, and the event yields can be found at Results.txt
