# PathHomologyDataAndScripts
For Chowdhury, Huntsman, and Yutin: "Path homologies of motifs and temporal network representations." Cf. https://arxiv.org/abs/2008.11885 

The scripts are not polished, but merely intended to make the paper's results reproducible with a modicum of (versus zero) effort. To run scripts, save the data files to an appropriate location and edit the directories in the scripts. For apnsMLPrepresentations.m, run interactively one code cell at a time for the appropriate data set, editing any filenames/locations as needed. 

Data:
  - sx-mathoverflow-a2q.txt
  - email-Eu-core-temporal.txt
  - out.facebook-wosn-wall.txt

Scripts: all M-files except for temporaldigraph.m (for which see below). Scripts ending in "1" are for the aggregated representations; the script apnsMLPrepresentations.m is for the MLP/layered representation; and the scripts ending in "TD" are for the temporal digraph representations. For path homology MATLAB code called by these scripts, see https://github.com/SteveHuntsmanBAESystems/PerformantPathHomology

Code: temporaldigraph.m produces a temporal digraph from a DCN. Lower-level/faster approaches are possible, but this exploits MATLAB's digraph functions. 
