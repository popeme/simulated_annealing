%an example of running the annealing function

%say I want to anneal for negative O-information on the HCP data... start
%by loading in the HCP data

load grandaverage_HCP.mat
dataset = FC; %this is the dataset I'm optimizing from

%because this is a matrix, and I want to index both rows and columns for the cost function, I set flag_mat equal to 1. 
%If it were a time-series (index only on rows), I would set flag_mat equal to zero
%I also don't need the within-trial output.
flag_mat = 1;
flag_out = 0;

%then set the annealing parameters to reasonable values-- play with these to make the algorithm converge better
H = 1e04; %maximum number of annealing steps
T_slope = 10; %the higher this number, the faster the temperature decreases
H_cut = 5000; %cut off after this number of steps with no better solution
T0 = 1; %initial temperature

%I would like to anneal for subsets of size 10
setsize = 10;

%and minimize the O-information (my cost function)
minmax = -1;

%starting with 1000 trials and recording only the first trial

trials = 1000;
record_tc = 1; %record the progress of this many trials (don't pick too many-- memory intensive!)

%so now I just call the annealing function with my cost function, which is
%calcO_logdet2-- this function just outputs the O-information. Switch in
%your own function just by typing its name after the '@'.
[optvarall,inivarall,optcostall,inicostall,opttc] = fcn_anneal(@calcO_logdet2,dataset,setsize,H,T_slope,H_cut,T0,trials,minmax,flag_out,flag_mat,record_tc);

%see whether or not the algorithm converged on the first trial
figure; plot(opttc)

%see the nodes chosen by the algorithm
figure; imagesc(optvarall)

