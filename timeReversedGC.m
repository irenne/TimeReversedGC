function [h_diffTRGC, stats] = timeReversedGC(x,y, varargin)
% returns the test decision of the difference-based time-reversed Granger-
% causality test and other Granger causality variants, as described in [1].
%
% IN   x, y         - two one-dimensional signals of the same length
% OUT  h_diffTRGC   - test decision for the difference-based time-reversed
%                     Granger-causality (diff-TRGC) 
%                       1  : causality x => y
%                       0  : no causality
%                       -1 : causality y => x    
%      stats        - statistics and test decisions for several variants
%       stats.lagOrder      
%                  lagOrder at which all tests were computed
%       stats.p_GCxtoy, stats.h_GCxtoy   
%                  p-value and test decision for the null hypothesis that
%                  x Granger causes y  (standard F-Test)
%       stats.p_GCytox, stats.h_GCytox   
%                  p-value and test decision for the null hypothesis that
%                  y Granger causes x  (standard F-Test)
%       stats.h_diffTRGC, stats.ci_diffTRGC
%                  decision and bootstrap confidence interval of the 
%                  diff-TRGC statistic 
%       stats.h_GCnet, stats.ci_GCnet
%                  decision and bootstrap confidence interval of the net-GC
%                  statistic (1 : x => y; 0: no causality, -1: y => x) 
%       stats.h_GCnet_fut, stats.ci_GCnet_fut
%                  decision and bootstrap confidence interval of the net-GC 
%                  statistic on the reversed time series 
%                  (1 : x_reversed => y_reversed; 0: no causality,
%                  -1: y_reversed => x_reversed) 
%       stats.h_conjTRGC
%                  decision of conjunction-based time-reversed Granger
%                  causality (conj-TRGC)
%                    1 : x => y        if h_GCnet = 1 and h_GCnet_fut = -1 
%                    -1: y => x        if h_GCnet = -1 and h_GCnet_fut = 1 
%                     0: no causality  otherwise
%
% OPTIONAL PARAMETERS
%   [h_diffTRGC, stats] = timeReversedGC(...,'lagOrder', P) 
%           considers P time lags for all tests
%           default: P = []; P is determined using the Bayesian Information 
%                    Criterion (BIC) 
%   [h_diffTRGC, stats] = timeReversedGC(...,'alpha', ALPHA) 
%           performs all tests at the significance level ALPHA. 
%           default: ALPHA = 0.05
%   [h_diffTRGC, stats] = timeReversedGC(...,'nBootstrapIterations', NBOOT) 
%           uses NBOOT bootstrap iterations 
%           default: NBOOT = 100 
%
% References: 
%   [1] I. Winkler, D. Panknin, D. Bartz, K.-R. Müller, S. Haufe. Validity 
%   of time reversal for testing Granger causality (2016). IEEE 
%   Transactions on Signal Processing, 64(11), 2746-2760.
%   [2] S. Haufe, V.V. Nikulin, K.-R. Müller, K. R., G. Nolte (2013). A 
%   critical assessment of connectivity measures for EEG data: a simulation 
%   study. Neuroimage, 64, 120-133.

 
%% Parse optional inputs
parser = inputParser; 
parser.addParamValue('lagOrder', []); %lag order, empty matrix: determined by BIC
parser.addParamValue('nBootstrapIterations', 100); % number of bootstrap samples
parser.addParamValue('alpha', 0.05); %significance level
parser.parse(varargin{:})
options = parser.Results; 

%% Standardise input time series
if (length(x) ~= length(y))
    error('x and y must be the same length');
end

x = (x - mean(x))/std(x);   
y = (y - mean(y))/std(y); 
x = x(:)'; 
y = y(:)'; 

%% Determine lag order by BIC if requested
if isempty(options.lagOrder)    
    p = getLagOrderByBIC(x, y);   
else 
    p = options.lagOrder; 
end

%% Prepare dependent and independent variables for Granger causality regressions
%dependent variables: cut x and y at the start and the end
y_dep = y(p+1:end-p);
x_dep = x(p+1:end-p); 

%independent variables: create matrices of past and future x and y
Y_past = zeros(p, length(x_dep));
X_past = zeros(p, length(x_dep));
Y_fut = zeros(p, length(x_dep));
X_fut = zeros(p, length(x_dep));
for i = 1:p   
    Y_past(i,:) = y(p+1-i:end-p-i);
    X_past(i,:) = x(p+1-i:end-p-i);
    Y_fut(i,:)  = y(p+1+i:end-p+i);
    X_fut(i,:)  = x(p+1+i:end-p+i);
end
Z_past = [Y_past; X_past]; 

%% Compute standard GC-tests: F-tests of regressions on past variables
%get RSS (Residual Sum of Square)
RSSy_res = norm( y_dep - y_dep * Y_past' * inv(Y_past * Y_past') * Y_past)^2;   
RSSy_full = norm( y_dep - y_dep * Z_past' * inv(Z_past * Z_past') * Z_past)^2; 
RSSx_res = norm( x_dep - x_dep * X_past' * inv(X_past * X_past') * X_past)^2;   
RSSx_full = norm( x_dep - x_dep * Z_past' * inv(Z_past * Z_past') * Z_past)^2; 

%compute standard GC F-test
df1 = p;  % cost in d.f (number of new variables added) 
df2 = length(x_dep) - 2*p; % degress of freedom remaining 
Fstaty = ((RSSy_res - RSSy_full) / df1 ) / (RSSy_full / df2);
Fstatx = ((RSSx_res - RSSx_full) / df1 ) / (RSSx_full / df2); 
stats.lagOrder = p; 
stats.p_GCxtoy = 1 - fcdf(Fstaty,df1,df2); 
stats.h_GCxtoy = (stats.p_GCxtoy < options.alpha); 
stats.p_GCytox = 1 - fcdf(Fstatx,df1,df2);
stats.h_GCytox = (stats.p_GCytox < options.alpha); 

     
%% Bootstrap net-GC and Diff-TRGC
% We bootstrap the residuals of regressions on both past and future values
Z_pastfut = [Z_past; X_fut; Y_fut]; 
x_fit =  x_dep * Z_pastfut' * inv(Z_pastfut * Z_pastfut') * Z_pastfut;
y_fit = y_dep * Z_pastfut' * inv(Z_pastfut * Z_pastfut') * Z_pastfut;
residuals = [ x_dep - x_fit; y_dep - y_fit]';
ci = bootci(options.nBootstrapIterations, {@(e)calcDiffNetGC(x_fit'+e(:,1),y_fit'+e(:,2), ...
   X_past', Y_past', X_fut', Y_fut'),residuals}, 'type', 'per', 'alpha', options.alpha); 
ci_GCnet = ci(:,1); 
ci_GCnet_fut = ci(:,2); 
ci_diffTRGC = ci(:,3); 

%collect results: net Granger causality
if ci_GCnet(1) > 0 
    stats.h_GCnet = 1; 
elseif ci_GCnet(2) < 0
    stats.h_GCnet = -1; 
else
    stats.h_GCnet = 0; 
end
stats.ci_GCnet = ci_GCnet;  

%collect results: net Granger causality on time-reversed signals
if ci_GCnet_fut(1) > 0 
    stats.h_GCnet_fut = 1; 
elseif ci_GCnet_fut(2) < 0
    stats.h_GCnet_fut = -1; 
else
    stats.h_GCnet_fut = 0; 
end
stats.ci_GCnet_fut = ci_GCnet_fut; 

%collect results: difference-based time-reversed Granger causality
if ci_diffTRGC(1) > 0 
    h_diffTRGC = 1; 
elseif ci_diffTRGC(2) < 0
    h_diffTRGC = -1; 
else
    h_diffTRGC = 0; 
end
stats.h_diffTRGC = h_diffTRGC; 
stats.ci_diffTRGC = ci_diffTRGC; 

%infer conjunction-based time-reversed Granger causality
if stats.h_GCnet == 1 && stats.h_GCnet_fut == -1
    stats.h_conjTRGC = 1; 
elseif stats.h_GCnet == -1 && stats.h_GCnet_fut == 1
    stats.h_conjTRGC = -1; 
else
    stats.h_conjTRGC = 0; 
end

return

%Determine the number of lags using the Bayesian Information Criterion 
function P = getLagOrderByBIC(x, y)
    pmax = min(120, floor(length(y)/4) -1);   
    
    y_dep = y(pmax+1:end); 
    Z_past = zeros(2*pmax, length(y_dep));
    for i = 1:pmax  
        Z_past(2*i-1, :) = y(pmax+1-i:end-i); 
        Z_past(2*i, :) = x(pmax+1-i:end-i);
    end
    
    [C,R] = qr(sparse(Z_past'),sparse(y_dep'));
    C = full(C); 
    residuals = C(3:end);
    RSSperP = zeros(1, pmax);
    for i = 1:pmax
        RSSperP(i) = norm(residuals(2*i-1:end))^2;
    end
    
    n = length(y_dep); 
    PredictionTerm = n*log( RSSperP./n); 
    ComplexityTerm = (1:pmax).*log(n); 
    BIC = PredictionTerm + ComplexityTerm; 
    [foo, P] = min(BIC);  
    return
    
%Given the prepared regression matrices, calculate netGC, netGC on the 
%reversed series, and the difference
function GCnets = calcDiffNetGC(x_dep, y_dep, X_past, Y_past, X_fut, Y_fut)
    x_dep = x_dep';
    y_dep = y_dep'; 
    X_past = X_past'; 
    Y_past = Y_past'; 
    Z_past = [X_past; Y_past]; 
    X_fut = X_fut'; 
    Y_fut = Y_fut'; 
    Z_fut = [X_fut; Y_fut]; 
    
    RSSy_res = norm( y_dep - y_dep * Y_past' * inv(Y_past * Y_past') * Y_past)^2;   
    RSSy_full = norm( y_dep - y_dep * Z_past' * inv(Z_past * Z_past') * Z_past)^2; 
    RSSx_res = norm( x_dep - x_dep * X_past' * inv(X_past * X_past') * X_past)^2;   
    RSSx_full = norm( x_dep - x_dep * Z_past' * inv(Z_past * Z_past') * Z_past)^2; 
    GCnet_xtoy = log(RSSy_res / RSSy_full) - log(RSSx_res / RSSx_full);
    
    RSSy_res_fut = norm( y_dep - y_dep * Y_fut' * inv(Y_fut * Y_fut') * Y_fut)^2;   
    RSSy_full_fut = norm( y_dep - y_dep * Z_fut' * inv(Z_fut * Z_fut') * Z_fut)^2; 
    RSSx_res_fut = norm( x_dep - x_dep * X_fut' * inv(X_fut * X_fut') * X_fut)^2;   
    RSSx_full_fut = norm( x_dep - x_dep * Z_fut' * inv(Z_fut * Z_fut') * Z_fut)^2; 
    GCnet_xtoy_fut = log(RSSy_res_fut / RSSy_full_fut) - log(RSSx_res_fut / RSSx_full_fut);
    
    diff = GCnet_xtoy - GCnet_xtoy_fut;
    
    GCnets = [GCnet_xtoy, GCnet_xtoy_fut, diff]; 
    return