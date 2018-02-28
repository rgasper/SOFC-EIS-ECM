function [fitted_params, final_error] = fit_eis_dat(exp_dat, params, ub, lb)
%% Setup
%Ray Gasper, 2018, UMass Amherst
%fits an EIS curve using a predefined effective circuit model that is
%typical for solid-oxide fuel cells
%Produces a Nyquist plot with the experimental data, total ECM fit, and ECM
%split up element-wise
%Currently ignores positive imaginary parts of the EIS, thus also the L
%element
%ECM structure: R-RQ1-RQ2-GE-FLW
%R: Resistor; parameters: Resistance
%RQ: RQ element; parameters: Yq, nq, Resistance
%GE: Gerischer element; parameters: Tc, Resistance
%FLW: Generalized Finite-Length Warburg element; parameters: Tw, Nw, Resistance

%I advise using a relatively good initial guess or reasonable
%constraints in order to ensure the realism of the fit

%exp_dat should be a string containing the address of a csv with EIS data
%in it.
% EIS data structure:
% Frequencies, Real, Imaginary
%   dat      , dat , dat
%   dat      , dat , dat
%   dat      , dat , dat
%   dat      , dat , dat
%   dat      , dat , dat
%   etc

%initial guess, upper bound, and lower bound are all 1x12 vectors
%corresponding to the paramater list:
%  1      2       3      4     5       6       7       8     9      10       11    12
%R1(R),RQ1(Yq),RQ1(nq),RQ1(R),RQ2(Yq),RQ2(nq),RQ2(R),GE(Tc),GE(R),FLW(Tw),FLW(nw),FLW(R)

%read the experimental data
dat = csvread(exp_dat,1,0);
exp_dat = clean_eis(dat);
global Freq exp_i exp_r
Freq = exp_dat(:,1);
exp_r = exp_dat(:,2);
exp_i = exp_dat(:,3);

%initialize fmincon inputs
if isempty(params)
    %if there's no initial guess use a random one
    %for the SOFCS initial parameters, ub, and lb all 0<x<1 work alright
    params = rand(1,12);
end

A = [];
b = [];
Aeq = [];
beq = [];
if isempty(lb)
    lb = [0 0 0 0 0 0 0 0 0 0 0 0];
end
if isempty(ub)
    ub = [1 1 1 1 1 1 1 1 1 1 1 1];
end

%% Fit the model
[fitted_params, final_error] = ...
    fmincon(@ecm_min_fit, params, A, b, Aeq, beq, lb, ub);

%% Plot
R_shift = R_element(Freq,[],1); %unit horizontal shift element for fancy plotting

% to modify the ECM structure, please modify in ecm_min_fit first then copy here, to ensure consistency
R = R_element(Freq,[], fitted_params(1));
RQ1 = RQ_element(Freq, fitted_params(2:3), fitted_params(4)); % Yq, nq
RQ2 = RQ_element(Freq, fitted_params(5:6), fitted_params(7)); % Yq, nq
GE = GE_element(Freq, fitted_params(8), fitted_params(9)); % Tc
FLW = FLW_element(Freq, fitted_params(10:11), fitted_params(12));% Tw, nw

Rr = (:,2);           Ri = R(:,3);
RQr1 = R1(:,2);       RQi1 = RQ1(:,3);
RQr2 = RQ2(:,2);      RQi2 = RQ2(:,3);
GEr = GE(:,2);        GEi = GE(:,3);
FLWr = FLW(:,2);      FLWi = FLW(:,3);
Sr = R_shift(:,2);    Si = R_shift(:,3);

sim_r = Rr + RQr1 + RQr2 + GEr + FLWr;
sim_i = Ri + RQi1 + RQi2 + GEi + FLWi;

plot(exp_r,exp_i,'ok','LineWidth',1)
xs = xlim;
ylim([-xs(2),0])
set(gca,'Ydir','reverse')
title('EIS Fit with ECM')
xlabel('Z_{real}(\Omega cm^2)')
ylabel('Z_{imaginary}(\Omega cm^2)')
axis manual
hold on
plot(sim_r,sim_i,'b-','LineWidth',1)
plot(RQr1+Sr*(fitted_params(1)), RQi1,'-r')
plot(RQr2+Sr*(fitted_params(1)+fitted_params(4)), RQi2,'-g')
plot(GEr+Sr*(fitted_params(1)+fitted_params(4)+fitted_params(7)), GEi,'-m')
plot(FLWr+Sr*(fitted_params(1)+fitted_params(4)+fitted_params(7)+fitted_params(9)), FLWi,'c-')
hold off
legend('Exp. Data','ECM fit','RQ1','RQ2','GE','FLW')

%% function to clean data by excluding positive imaginary components
function dat = clean_eis(exp_dat)
 % cleans EIS data of Zi > 0 elements
 omega = exp_dat(:,1);
 Zr = exp_dat(:,2);
 Zi = exp_dat(:,3);
 
 bad = any(Zi>0,2);
 if sum(bad) > length(Zi)/2
    Zi = Zi.*-1;
    bad = any(Zi>0,2);
 end
 
 c_omega = omega(~bad,:);
 c_Zr = Zr(~bad,:);
 c_Zi = Zi(~bad,:);
 
 dat = [c_omega, c_Zr, c_Zi]; 
end
%% Function to calculate total squared error of a particular ECM parameter set
function err = ecm_min_fit(params)
% ecm fit function with defined circuit for implementation in matlab's
% fminsearch funtion
% note that this is equal to the complex conjugate of the true complex
% error value

%If you want to modify the ECM structure, make your changes here then copy to the plotting section
Rs = R_element(Freq,[], params(1));
RQ1s = RQ_element(Freq, params(2:3), params(4)); % Yq, nq
RQ2s = RQ_element(Freq, params(5:6), params(7)); % Yq, nq
GEs = GE_element(Freq, params(8), params(9)); % Tc
FLWs = FLW_element(Freq, params(10:11), params(12));% Tw, nw

fit_r = Rs(:,2) + RQ1s(:,2) + RQ2s(:,2) + GEs(:,2) + FLWs(:,2);
fit_i = Rs(:,3) + RQ1s(:,3) + RQ2s(:,3) + GEs(:,3) + FLWs(:,3);
% fit = [Freq, fit_r, fit_i];

err_r = sum((exp_r - fit_r).^2);
err_i = sum((exp_i - fit_i).^2);

err = err_r + err_i;
end
%% Element functions
function dat = R_element(omega, params, R)
 % from a given set of frequencies generates the Zr and Zi for a resistor

 Z = R;
 
 dat(:,1) = omega;
 dat(:,2) = real(Z);
 dat(:,3) = imag(Z);
end
function dat = RQ_element(omega, params, R)
 % from a given set of frequencies generates the Zr and Zi for a RQ
 % element
 Y_q = params(1);
 n_q = params(2);
 
 Q = 1 ./ (Y_q .* (omega.*1i).^n_q);
 Z = R ./ ( 1 + R.*(Q.^-1) );
 
 dat(:,1) = omega;
 dat(:,2) = real(Z);
 dat(:,3) = imag(Z); 
end
function dat = GE_element(omega, params, R)
 % from a given set of frequencies generates the Zr and Zi for a
 % Gerischer Element
 t_c = params(1);
  
 Z = R ./ sqrt(1+omega.*t_c.*1i); 
 
 dat(:,1) = omega;
 dat(:,2) = real(Z);
 dat(:,3) = imag(Z); 
end
function dat = FLW_element(omega, params, R)
 % from a given set of frequencies generates the Zr and Zi for a
 % Finite-Length Warburg Element
 T_w = params(1);
 n_w = params(2);
 
 Z = R .* ( tanh( (omega.*T_w.*1i).^n_w ) )./...
              ( (omega.*T_w.*1i ).^n_w ) ;
 
 dat(:,1) = omega;
 dat(:,2) = real(Z);
 dat(:,3) = imag(Z); 
end
end
