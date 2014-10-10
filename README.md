Here you can find all my code used for the vbf-bbbb analysis.

The main file is the one named vbf-bbbb.C, which runs over the samples given as a parameter. It runs over the small ntuples saved in /afs/cern.ch/work/a/ariostas/public/vbf-bbbb/ . The small ntuples were created with the code found on the Selection folder.

To run use root vbf-bbbb.C+ to run over all samples or vbf-bbbb.C+(\"name of background sample\") to use a specific background sample.

There is also an experimental selection code that uses a similar method to the one used in gf-bbbb. It constructs di-bjets, makes sure the reconstructed object looks like a Higgs, and then constructs di-light-jets with high separation such that the two di-bjets lie between them. The small ntuples are saved in /afs/cern.ch/work/a/ariostas/public/vbf-bbbb_test/ . To run over these ntuples use vbf-bbbb_experimental.C+ . It works pretty much the same as vbf-bbbb.C . Right now, the signal yields are fairly low, but I think it is much more reliable than the normal method, so it would be intresting to see if it can be improved.
