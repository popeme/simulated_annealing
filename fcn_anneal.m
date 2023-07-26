function [optvarall,inivarall,optcostall,inicostall,opttc] = fcn_anneal(cost_func,dataset,setsize,H,T_slope,H_cut,T0,trials,minmax,flag_out,flag_mat,record_tc,varargin)

%this function performs simulated annealing for a number of trials on the dataset
%provided, minimizing or maximizing the cost function provided. Simulated
%annealing works by becoming increasingly less stochastic as the trial
%progresses, i.e. it becomes more and more like a greedy optimization algorithm. 
% This progression is governed by the temperature function and annealing
%parameters. The higher the initial temperature(T0), the more random the process is in the
%beginning, and the higher the temperature's slope (T_slope), the more rapidly
%it becomes like a greedy optimization. For unknown cost function
%landscapes, I recommend a low slope with a high total number of annealing
%steps (H) (with H_cut=H), until you get a better idea of how easy it is for the algorithm
%to converge. Note that increasing H_cut will save computation time, but can decrease
%convergence.

%Inputs:
%  cost_func: a function performed on the dataset provided. It must output
%       a single number, which is the property being optimized.
%  dataset: an array with size N x data, where N are the variables subject to the
%       optimization, and data can be operated on by the cost_func
%  setsize: the number of variables to be considered for annealing
%  H: maximum number of annealing steps the optimization will take
%  T_slope: governs the slope of the temperature decrease: higher =
%       decreasing faster
%  H_cut: number of steps to cut off after if have found no better value
%  T0: initial temperature
%  trials: number of trials of annealing to be performed
%  minmax: set to +1 to maximize the cost_func or -1 to minimize the
%       cost_func
%  flag_out: if true, then the annealer will display the cost as annealing
%       proceeds, set to false to suppress this output.
%  flag_mat: if true, then the annealer assumes that the data it is
%       annealing is a matrix-- subsets of variables have to be indexed on
%       both rows and columns
%  record_tc: number of trials to record optimization time course for
%  varargin: passes any additional inputs to the cost function
%
%Outputs:
%  optvarall: logical, NxT, stores the chosen optimized variables for all
%       trials
%  inivarall: logical, NxT, stores the initial variables for all trials
%  optcostall: double, Tx1, stores the optimized value of the cost function to all
%       trials
%  inicostall: double, Tx1, stores the initial value of the cost function
%       for all trials
%  opttc: double, Hx1, stores the optimization time course for record_tc
%       trials

%Assembled by: Maria Pope, Indiana University, July 2023
%Based on code written by Olaf Sporns, Indiana University, a long time ago


% set annealing variables
Texp = 1-T_slope/H;  
Tall = T0.*Texp.^[1:1:H];
minmax = minmax.*-1; %I have done this to make the input more intuitive: minimizing means inputting a negative 1, 
% rather than minimizing means inputting a positive 1

%get total number of variables
N = size(dataset,1);

%initialize storage variables
optvarall = logical(zeros(N,trials));
inivarall = logical(zeros(N,trials));
optcostall = zeros(trials,1);
inicostall = zeros(trials,1);
opttc = zeros(H,record_tc);

% loop over trials
for t=1:trials
    
    disp(strcat('Now beginning trial ',num2str(t),' of ', num2str(trials)))

    h = 0; hcnt = 0;

    % select initial set of variables
    rp = randperm(N);
    inibag = rp(1:setsize);
    minbag = inibag;
    minbagglobal = minbag;
    
    %calculate initial cost
    if flag_mat==1
        cost = minmax.*cost_func(dataset(inibag,inibag),varargin{:});
    else
        cost = minmax.*cost_func(dataset(inibag,:),varargin{:});
    end
    
    %store initial cost
    inicost = cost; 
    mincost = inicost; 
    lowcost = inicost;


    costh = zeros(H,1);

    while (h<H)&&(hcnt<H_cut)

        h = h+1; hcnt = hcnt+1;

        % current temperature
        Tc = T0*Texp^h;

        if ((mod(h,H/10)==0)&&(flag_out==1))
            disp(['at step ',num2str(h),' - elapsed time ',num2str(toc),' - lowest cost = ',num2str(mincost)]);
        end

        % swap in new variables
        newbag = minbag;
        % swap how many?
        swap = min(setsize,ceil(abs(randn))); % no more than setsize
        candidates = setdiff(1:N,newbag);
        newbag(randperm(setsize,swap)) = candidates(randperm(length(candidates),swap));

        %calculate the cost of the new variables
        if flag_mat==1
            cost = minmax.*cost_func(dataset(newbag,newbag),varargin{:});
        else
            cost = minmax.*cost_func(dataset(newbag,:),varargin{:});
        end
        costnew = cost;

        % annealing
        randcon = rand < exp(-(costnew-lowcost)/Tc);
        if (costnew < lowcost) || randcon
            minbag = newbag;
            lowcost = costnew;
            % is this a new absolute best?
            if (lowcost < mincost)
                minbagglobal = minbag;
                mincost = lowcost;
                hcnt = 0;
            end
        end

        % time course of optimization
        costh(h) = lowcost;

    end
    %package the data
    optvarall(minbagglobal,t) = 1;
    inivarall(inibag,t) = 1;
    optcostall(t) = mincost;
    inicostall(t) = inicost;
    if t<=record_tc;opttc(:,t) = costh;end
    
    
end
end