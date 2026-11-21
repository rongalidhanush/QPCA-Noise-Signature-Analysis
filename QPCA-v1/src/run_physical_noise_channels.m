function run_physical_noise_channels(resultDir)
%RUN_PHYSICAL_NOISE_CHANNELS
% Small physics sanity check for the five noise mechanisms represented
% in the project. This is supplementary; N-QPCA remains the main method.

if nargin<1
    resultDir="results";
end

fprintf("\nPhysical noise-channel sanity check\n");

channels = ["bitflip","phaseflip","depolarizing","amplitude_damping","correlated"];
probability = 0.05;

results = table('Size',[numel(channels),4], ...
    'VariableTypes',{'string','double','double','double'}, ...
    'VariableNames',{'Channel','TracePreservingError','HermiticityError','MinimumEigenvalue'});

for k=1:numel(channels)
    rho = eye(2)/2;
    out = apply_single_channel(rho,channels(k),probability);
    results.Channel(k)=channels(k);
    results.TracePreservingError(k)=abs(trace(out)-1);
    results.HermiticityError(k)=norm(out-out','fro');
    results.MinimumEigenvalue(k)=min(real(eig((out+out')/2)));
end

writetable(results,fullfile(resultDir,"physical_channel_sanity.csv"));
disp(results);

end

function out=apply_single_channel(rho,name,p)
I=eye(2);
X=[0 1;1 0];
Z=[1 0;0 -1];
Y=[0 -1i;1i 0];

switch lower(name)
    case "bitflip"
        K={sqrt(1-p)*I,sqrt(p)*X};
    case "phaseflip"
        K={sqrt(1-p)*I,sqrt(p)*Z};
    case "depolarizing"
        K={sqrt(1-p)*I,sqrt(p/3)*X,sqrt(p/3)*Y,sqrt(p/3)*Z};
    case "amplitude_damping"
        K={ [1 0;0 sqrt(1-p)], [0 sqrt(p);0 0] };
    otherwise
        % A simple correlated two-qubit Z⊗Z channel is represented
        % separately from the one-qubit sanity check.
        ZZ=kron(Z,Z);
        r2=kron(rho,rho);
        out=(1-p)*r2+p*ZZ*r2*ZZ';
        return;
end

out=zeros(size(rho));
for i=1:numel(K)
    out=out+K{i}*rho*K{i}';
end
end
