% Simulation of unidirectional causal interaction in the presence of noise,
% which reproduces the result presented in Figure 1b in [1]. 
% 
% Two-dimensional auto-regressive time series 
%   z_latent = [x_latent; y_latent] and z_noise = [x_noise; y_noise]
% with a causal interation x_latent => y_latent and no causal interaction
% between x_noise and y_noise are superimposed as: 
%   z = (1 - gamma) * z_latent + gamma * z_noise
% The signal-to-noise ratio gamma is varied. 
% 
% Standard Granger-causality (GC), net Granger causality (net-GC), d
% difference based time-reversed Granger causality (diff-TRGC) and 
% conjunction-based TRGC (conj-TRGC) are computed on z. 
% A summary plot is generated. 
% 
% References: 
%   [1] I. Winkler, D. Panknin, D. Bartz, K.-R. Müller, S. Haufe. Validity 
%   of time reversal for testing Granger causality (2016). IEEE 
%   Transactions on Signal Processing, 64(11), 2746-2760.
%   [2] S. Haufe, V.V. Nikulin, K.-R. Müller, K. R., G. Nolte (2013). A 
%   critical assessment of connectivity measures for EEG data: a simulation 
%   study. Neuroimage, 64, 120-133.

 

clear

nRepetitions = 100; %number of repetitions of the simulation 

T = 2000; %length of the generated time series 
sigmaA = 0.2; %standard deviation of the AR-coefficients
p_latent = 5; %lag order of z_latent
p_noise = 5; %lag order of z_noise (set to 0 for non auto-correlated noise)
isNoiseLinearlyMixed = true; %if false, x_noise and y_noise are independent

gammas = [0, 0.25, 0.5, 0.75, 0.9, 1]; %signal-to-noise ratios


h_diffTRGC = zeros(length(gammas),nRepetitions); 
h_conjTRGC = zeros(length(gammas),nRepetitions); 
h_netGC = zeros(length(gammas),nRepetitions); 
h_GCxtoy = zeros(length(gammas),nRepetitions); 
h_GCytox = zeros(length(gammas),nRepetitions); 
fprintf('Run %d repetitions  ', nRepetitions)
for n = 1:nRepetitions
    fprintf('.')
    %generate z_latent and z_noise
    isCausalInteraction = true; 
    z_latent = create2dVAR(T, p_latent, sigmaA, isCausalInteraction); 
    
    isCausalInteraction = false; 
    z_noise = create2dVAR(T, p_noise, sigmaA, isCausalInteraction); 
    if isNoiseLinearlyMixed 
        B = randn(2,2); %mixing matrix
        B = B/sqrt(abs(det(B))); %make sure that det(B) = +-1
        z_noise = B * z_noise; 
    end
    
    %for each signal-to-noise ratio gamma, superimpose z_latent and z_noise
    %and evaluate causality methods
    for idxGamma = 1:length(gammas)  
        gamma = gammas(idxGamma); 
        z = (1-gamma)* z_latent + gamma * z_noise; 
        
        [h_diffTRGC(idxGamma, n),  stats] = timeReversedGC(z(1,:),z(2,:)); 
        h_conjTRGC(idxGamma, n) = stats.h_conjTRGC; 
        h_netGC(idxGamma,n) = stats.h_GCnet; 
        h_GCxtoy(idxGamma, n) = stats.h_GCxtoy; 
        h_GCytox(idxGamma, n) = stats.h_GCytox;
    end
end
fprintf('\n')


%Plot results
figure
subplot(1,2,1)
plot(gammas, mean(h_GCxtoy == 1,2),'h-', 'Color', [0.5, 0.5, 0.5])
hold on
plot(gammas, mean(h_netGC == 1,2),'+k-')
plot(gammas, mean(h_diffTRGC == 1,2),'dr-')
plot(gammas, mean(h_conjTRGC == 1,2),'o-', 'Color',  [0 0.45 0.75])

set(gca, 'Xtick', gammas, 'YGrid', 'on')
xlim([0, 1]); ylim([0 1]); 
xlabel('Noise level')
ylabel('True Positives')

subplot(1,2,2)
plot(gammas, mean(h_GCytox==1,2),'h-', 'Color', [0.5, 0.5, 0.5])
hold on
plot(gammas, mean(h_netGC==-1,2),'+k-')
plot(gammas, mean(h_diffTRGC == -1,2),'dr-')
plot(gammas, mean(h_conjTRGC == -1,2),'o-', 'Color',[0 0.45 0.75] )

set(gca, 'Xtick', gammas, 'YGrid', 'on')
xlim([0, 1]); ylim([0 1])
xlabel('Noise level')
ylabel('False Positives')

legend('Standard GC', 'Net GC', 'Diff-TRGC', 'Conj-TRGT')