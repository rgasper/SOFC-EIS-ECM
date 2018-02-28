%A wee example script for how to use the EIS effective circuit model fitting function
%Ray Gasper, 2018, UMass Amherst
clear;clc

%% simplest possible run- no initial guess, no constraints
%If you run this multiple times you'll notice the ECM can change
%significantly, this means there are many local minima in the
%error:parameter space- only if you're lucky will this fit be realistic
%sometimes it won't converge at all, often (not always) meaning the random 
%initial guess is in a particularly bad region of parameter space
figure(1)
[fit_1, err_1] = fit_eis_dat('exp_data_fine.csv',[],[],[]);
title('Random Initial Guess & Default UB, LB')

%% with initial guess
%this is using a pretty good initial guess.
%Notice with re-running there's some, but little change- we're close to a good
%minima, with these values listed corresponding to realistic ones
guess=[0.01,0.01,0.75,0.15,0.25,0.75,0.05,0.01,0.05,0.05,0.5,0.3];
figure(2)
[fit_2, err_2] = fit_eis_dat('exp_data_fine.csv',guess,[],[]);
title('Good Initial Guess & Default UB, LB')

%% with bad UB and LB, no initial guess
% you need very good upper bounds and lower bounds, see the result of a poor set here
%just arbitrarily setting upper bounds on all the elements almost always
%causes failure to produce a good set
ub=[0.1,0.1,0.99,0.5,0.5,0.99,0.5,0.5,0.5,0.5,0.99,0.99];
lb=[0,0,0,0,0,0,0,0,0,0,0,0];
figure(3)
[fit_3, err_3] = fit_eis_dat('exp_data_fine.csv',[],ub,lb);
title('Random Initial Guess & BAD Custom UB, LB')

%% with good UB, no initial guess
ub = [1 1 1 1 1 1 1 1 0.2 1 1 1]; %We're limiting the resistance of the GE element
figure(4)
[fit_4, err_4] = fit_eis_dat('exp_data_fine.csv',[],ub,lb);
title('Random Initial Guess & Good Custom UB, Default LB')
