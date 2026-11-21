%% N-QPCA BASED QUANTUM NOISE SIGNATURE ANALYSIS
% Main project entry point.
%
% The central algorithm is a 15-dimensional QPCA simulation.
% The 15x15 covariance-derived density matrix is embedded into
% a 16-dimensional Hilbert space, represented by 4 system qubits.
%
% Run from the project root:
%   addpath("src")
%   main_NQPCA

clear; clc; close all;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(projectRoot,"src"));

dataFile = fullfile(projectRoot,"data","noise_dataset.csv");
resultDir = fullfile(projectRoot,"results");
figDir = fullfile(resultDir,"figures");

if ~exist(resultDir,"dir"), mkdir(resultDir); end
if ~exist(figDir,"dir"), mkdir(figDir); end

fprintf("\n============================================================\n");
fprintf(" N-QPCA QUANTUM NOISE SIGNATURE ANALYSIS\n");
fprintf("============================================================\n");

%% 1. Load data
T = readtable(dataFile);
featureNames = T.Properties.VariableNames(3:end);

X = T{:,3:end};
noiseType = string(T.NoiseType);
noiseStrength = T.NoiseStrength;

fprintf("Rows                 : %d\n", size(X,1));
fprintf("Feature dimension N  : %d\n", size(X,2));
fprintf("Noise conditions     : %d\n", numel(unique(noiseType)));

%% 2. Standardize the 15-dimensional feature matrix
mu = mean(X,1);
sigma = std(X,0,1);
sigma(sigma == 0) = 1;
Xz = (X-mu)./sigma;

%% 3. Build covariance and density matrix
C15 = cov(Xz);
rho15 = C15/trace(C15);

% Embed N=15 into 2^4=16 dimensional Hilbert space.
rho16 = zeros(16,16);
rho16(1:15,1:15) = rho15;
rho16 = (rho16+rho16')/2;
rho16 = rho16/trace(rho16);

fprintf("Covariance size      : %d x %d\n",size(C15,1),size(C15,2));
fprintf("Quantum matrix size  : %d x %d\n",size(rho16,1),size(rho16,2));
fprintf("System qubits        : 4\n");

%% 4. Classical numerical reference ONLY
[V,D] = eig(rho16);
lambdaClassical = real(diag(D));
[lambdaClassical,order] = sort(lambdaClassical,"descend");
V = V(:,order);

lambdaClassical(lambdaClassical < 0) = 0;
lambdaClassical = lambdaClassical/sum(lambdaClassical);

explainedClassical = 100*lambdaClassical;
cumExplained = cumsum(explainedClassical);

fprintf("\nClassical reference top eigenvalues:\n");
disp(lambdaClassical(1:5)');

%% 5. THE MAIN N-QPCA COMPUTATION
% QPE is run on a purification of rho16.
qpeBits = 8;
fprintf("\nRunning N-QPCA / QPE with %d ancilla qubits...\n",qpeBits);

qpca = n_qpca_qpe(rho16,qpeBits,3);

fprintf("QPE probability sum : %.12f\n",sum(qpca.probability));
fprintf("QPE resolution       : 1/%d = %.8f\n",2^qpeBits,1/2^qpeBits);

%% 6. Compare quantum QPCA estimates with classical reference
qLambda = qpca.estimatedEigenvalues(:);
cLambda = lambdaClassical(1:numel(qLambda));

absErr = abs(qLambda-cLambda);
relErr = absErr./max(cLambda,eps);

fprintf("\n------------------------------------------------------------\n");
fprintf("GLOBAL N-QPCA EVALUATION\n");
fprintf("------------------------------------------------------------\n");
fprintf("Component    Classical       QPCA/QPE        Abs Error\n");
for k = 1:numel(qLambda)
    fprintf("PC%-9d %.8f       %.8f       %.8f\n", ...
        k,cLambda(k),qLambda(k),absErr(k));
end
fprintf("Mean absolute error: %.8f\n",mean(absErr));

%% 7. QPCA scores using the dominant reference modes
nPC = 3;
scores = Xz*V(1:15,1:nPC);
Xrec = scores*V(1:15,1:nPC)';
reconstructionError = norm(Xz-Xrec,"fro")/norm(Xz,"fro");

fprintf("\nTop-%d classical reference explained variance: %.4f %%\n", ...
    nPC,cumExplained(nPC));
fprintf("Top-%d reconstruction error              : %.6f\n", ...
    nPC,reconstructionError);

%% 8. Noise-type evaluation in the QPCA score space
types = unique(noiseType,"stable");
typeResults = table('Size',[numel(types),7], ...
    'VariableTypes',{'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'NoiseType','MeanStrength','MeanPC1','MeanPC2','MeanPC3','MeanAbsPC1','PC1StrengthCorrelation'});

for i = 1:numel(types)
    idx = noiseType == types(i);
    typeResults.NoiseType(i) = types(i);
    typeResults.MeanStrength(i) = mean(noiseStrength(idx));
    typeResults.MeanPC1(i) = mean(scores(idx,1));
    typeResults.MeanPC2(i) = mean(scores(idx,2));
    typeResults.MeanPC3(i) = mean(scores(idx,3));
    typeResults.MeanAbsPC1(i) = mean(abs(scores(idx,1)));
    typeResults.PC1StrengthCorrelation(i) = simple_corr(noiseStrength(idx),abs(scores(idx,1)));
end

%% 9. Nearest-centroid validation in QPCA space
centroids = zeros(numel(types),nPC);
for i = 1:numel(types)
    centroids(i,:) = mean(scores(noiseType==types(i),:),1);
end

predicted = strings(size(noiseType));
for r = 1:size(scores,1)
    d = sum((centroids-scores(r,:)).^2,2);
    [~,j] = min(d);
    predicted(r) = types(j);
end
centroidAccuracy = mean(predicted==noiseType)*100;

fprintf("Nearest-centroid QPCA-space accuracy: %.2f %%\n",centroidAccuracy);

%% 10. Per-noise N-QPCA signatures
% Each noise condition is analyzed using the SAME global standardization.
perNoise = table('Size',[numel(types),8], ...
    'VariableTypes',{'string','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'NoiseType','Lambda1_Classical','Lambda2_Classical','Lambda3_Classical', ...
    'Lambda1_QPCA','Lambda2_QPCA','Lambda3_QPCA','Top3_QPCA_ExplainedPercent'});

for i = 1:numel(types)
    idx = noiseType == types(i);
    Cg = cov(Xz(idx,:));
    rg = Cg/trace(Cg);
    r16 = zeros(16,16);
    r16(1:15,1:15) = rg;
    r16 = (r16+r16')/2;
    r16 = r16/trace(r16);

    [~,dg] = eig(r16);
    lg = sort(real(diag(dg)),"descend");
    lg(lg<0)=0;
    lg = lg/sum(lg);

    qg = n_qpca_qpe(r16,qpeBits,3);

    perNoise.NoiseType(i) = types(i);
    perNoise.Lambda1_Classical(i) = lg(1);
    perNoise.Lambda2_Classical(i) = lg(2);
    perNoise.Lambda3_Classical(i) = lg(3);
    perNoise.Lambda1_QPCA(i) = qg.estimatedEigenvalues(1);
    perNoise.Lambda2_QPCA(i) = qg.estimatedEigenvalues(2);
    perNoise.Lambda3_QPCA(i) = qg.estimatedEigenvalues(3);
    perNoise.Top3_QPCA_ExplainedPercent(i) = 100*sum(qg.estimatedEigenvalues);
end

%% 11. Save numerical results
globalResults = table((1:3)',cLambda,qLambda,absErr,relErr, ...
    'VariableNames',{'Component','ClassicalEigenvalue','QPCALambda','AbsoluteError','RelativeError'});

strengthResults = table(noiseStrength,abs(scores(:,1)),noiseType, ...
    'VariableNames',{'NoiseStrength','AbsQPCA_PC1','NoiseType'});

writetable(globalResults,fullfile(resultDir,"nqpca_global_spectrum.csv"));
writetable(typeResults,fullfile(resultDir,"nqpca_noise_type_results.csv"));
writetable(perNoise,fullfile(resultDir,"nqpca_per_noise_signatures.csv"));
writetable(strengthResults,fullfile(resultDir,"nqpca_noise_strength_results.csv"));

qpeDistribution = table((0:length(qpca.probability)-1)', ...
    (0:length(qpca.probability)-1)'/2^qpeBits, qpca.probability, ...
    'VariableNames',{'Bin','EstimatedEigenvalue','Probability'});
writetable(qpeDistribution,fullfile(resultDir,"nqpca_qpe_distribution.csv"));

save(fullfile(resultDir,"nqpca_main_results.mat"), ...
    "T","X","Xz","mu","sigma","C15","rho15","rho16", ...
    "V","lambdaClassical","explainedClassical","cumExplained", ...
    "qpca","scores","Xrec","reconstructionError", ...
    "typeResults","perNoise","centroidAccuracy", ...
    "noiseType","noiseStrength","featureNames");

%% 12. Figures
figure('Name','N-QPCA Spectrum');
stem(1:15,lambdaClassical(1:15),'filled');
hold on;
stem(1:3,qLambda,'LineWidth',1.5);
grid on;
xlabel('Principal component');
ylabel('Eigenvalue');
title('15-Dimensional N-QPCA Spectrum');
legend('Classical reference','QPCA/QPE estimate','Location','northeast');
saveas(gcf,fullfile(figDir,"01_NQPCA_spectrum.png"));

figure('Name','QPE Distribution');
bar(qpeDistribution.EstimatedEigenvalue,qpeDistribution.Probability);
grid on;
xlabel('Estimated eigenvalue');
ylabel('Probability');
title('N-QPCA Quantum Phase Estimation Distribution');
xlim([0 1]);
saveas(gcf,fullfile(figDir,"02_QPE_distribution.png"));

figure('Name','Explained Variance');
bar(1:15,explainedClassical(1:15));
hold on;
plot(1:15,cumExplained(1:15),'o-','LineWidth',1.5);
grid on;
xlabel('Principal component');
ylabel('Percent');
title('N-QPCA Principal Spectrum and Cumulative Variance');
legend('Individual','Cumulative','Location','best');
saveas(gcf,fullfile(figDir,"03_explained_variance.png"));

figure('Name','QPCA Scores');
hold on;
for ii = 1:numel(types)
    idx = noiseType == types(ii);
    scatter(scores(idx,1),scores(idx,2),18,'filled');
end
grid on;
xlabel('QPCA PC1');
ylabel('QPCA PC2');
title('Noise Signatures in N-QPCA Score Space');
legend(types,'Location','best');
hold off;
saveas(gcf,fullfile(figDir,"04_noise_signatures.png"));

figure('Name','Noise Strength');
hold on;
for ii = 1:numel(types)
    idx = noiseType == types(ii);
    scatter(noiseStrength(idx),abs(scores(idx,1)),18,'filled');
end
grid on;
xlabel('Noise strength');
ylabel('|QPCA PC1|');
title('Noise Strength vs QPCA PC1 Magnitude');
legend(types,'Location','best');
hold off;
saveas(gcf,fullfile(figDir,"05_strength_vs_qpca.png"));

figure('Name','Reconstruction');
bar(1,reconstructionError);
ylim([0 max(0.3,reconstructionError*1.2)]);
grid on;
xticks(1); xticklabels({'Top 3 QPCA'});
ylabel('Normalized reconstruction error');
title('Reconstruction Using Dominant Principal Modes');
saveas(gcf,fullfile(figDir,"06_reconstruction_error.png"));

fprintf("\n============================================================\n");
fprintf(" N-QPCA PROJECT COMPLETE\n");
fprintf("============================================================\n");
fprintf("Results saved to: %s\n",resultDir);
fprintf("Main MAT file   : nqpca_main_results.mat\n");
fprintf("Figures folder  : results/figures/\n");
fprintf("============================================================\n");


function r = simple_corr(a,b)
a = a(:); b = b(:);
a = a-mean(a);
b = b-mean(b);
den = sqrt(sum(a.^2)*sum(b.^2));
if den == 0
    r = 0;
else
    r = sum(a.*b)/den;
end
end
