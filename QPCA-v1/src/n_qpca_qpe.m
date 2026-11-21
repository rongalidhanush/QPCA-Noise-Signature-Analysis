function result = n_qpca_qpe(rho,qpeBits,numComponents)
%N_QPCA_QPE  N-dimensional QPCA spectral estimation by QPE.
%
% INPUT
%   rho            - density matrix, dimension d x d
%   qpeBits        - number of phase-estimation ancilla qubits
%   numComponents  - number of dominant eigenvalues to report
%
% OUTPUT
%   result.probability
%   result.phaseGrid
%   result.estimatedEigenvalues
%   result.peakBins
%   result.unitary
%   result.purification
%
% The implementation is a state-vector simulation of the QPCA/QPE
% spectral step. The input state is a canonical purification of rho:
%
%       |Psi_rho> = sum_j sqrt(lambda_j)|u_j>|u_j^*>
%
% so QPE samples eigencomponents according to their density-matrix
% weights rather than uniformly.
%
% The phase unitary is:
%
%       U = exp(2*pi*i*rho)
%
% Hence an eigenvalue lambda gives an eigenphase lambda modulo 1.

rho = (rho+rho')/2;
rho = rho/trace(rho);

d = size(rho,1);
D = 2^qpeBits;

if D < d
    error("Number of QPE ancilla states must be >= density-matrix dimension.");
end

%% Construct U = exp(2*pi*i*rho)
U = expm(2*pi*1i*rho);

%% Construct a stable canonical purification from sqrt(rho)
[V,L] = eig(rho);
lam = real(diag(L));
lam(lam<0) = 0;
sqrtRho = V*diag(sqrt(lam))*V';
sqrtRho = (sqrtRho+sqrtRho')/2;

% Normalization should be one because trace(rho)=1.
sqrtRho = sqrtRho/sqrt(real(trace(sqrtRho*sqrtRho')));

%% State register:
%   ancilla: D amplitudes
%   system : d
%   reference : d
%
% Initial ancilla state is |+>^m.
psi = repmat(sqrtRho/sqrt(D),1,1,D);
psi = permute(reshape(psi,[d,d,D]),[3 1 2]);

%% Controlled-U^(2^k)
% Ancilla index uses the least-significant-bit convention.
Up = U;

for k = 0:qpeBits-1
    power = 2^k;
    if k == 0
        Up = U;
    else
        Up = Up*Up;
    end

    for a = 0:D-1
        if bitand(a,power) ~= 0
            page = squeeze(psi(a+1,:,:));
            page = Up*page;
            psi(a+1,:,:) = reshape(page,1,d,d);
        end
    end
end

%% Inverse QFT on the ancilla register
Fdag = inverse_qft_matrix(D);
psi2 = reshape(psi,[D,d*d]);
psi2 = Fdag*psi2;
psi = reshape(psi2,[D,d,d]);

%% Measurement probability distribution
probability = squeeze(sum(sum(abs(psi).^2,3),2));
probability = real(probability);
probability(probability<0)=0;
probability = probability/sum(probability);

phaseGrid = (0:D-1)'/D;

%% Extract dominant QPCA eigenvalue peaks
peakBins = select_qpca_peaks(probability,numComponents);
estimatedEigenvalues = phaseGrid(peakBins);

% QPE peak heights are affected by finite phase resolution and spectral
% leakage.  Therefore peak probability is NOT used to rank principal
% components.  QPCA components are ordered by their estimated eigenvalue.
[estimatedEigenvalues,ord] = sort(estimatedEigenvalues,"descend");
peakBins = peakBins(ord);

result.probability = probability;
result.phaseGrid = phaseGrid;
result.estimatedEigenvalues = estimatedEigenvalues;
result.peakBins = peakBins;
result.unitary = U;
result.purification = sqrtRho;
result.qpeBits = qpeBits;
result.phaseResolution = 1/D;
result.stateDimension = d;

end
