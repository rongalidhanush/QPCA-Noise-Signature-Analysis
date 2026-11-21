function Fdag = inverse_qft_matrix(N)
%INVERSE_QFT_MATRIX Construct the inverse quantum Fourier transform.
%
% [Fdag]_(a,b) = exp(-2*pi*i*a*b/N)/sqrt(N)

[a,b] = ndgrid(0:N-1,0:N-1);
Fdag = exp(-2*pi*1i.*a.*b/N)/sqrt(N);

end
