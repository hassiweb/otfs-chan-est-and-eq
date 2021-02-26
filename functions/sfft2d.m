% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

function X = sfft2d(x) % Symplectic FFT
[N,M,numLayers] = size(x);
X = zeros(M,N,numLayers);
for iter=1:numLayers
    X(:,:,iter) = fft(ifft(x(:,:,iter)).'); 
end
% a = M/sqrt(M*N);
a = N/sqrt(M*N);
X = a * X;     
end