% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

function x = isfft2d(X)  % Inverse Symplectic FFT
[M,N,numLayers] = size(X);
x = zeros(N,M,numLayers);
for iter=1:numLayers
    x(:,:,iter) = ifft(fft(X(:,:,iter)).'); 
end
a = N/sqrt(M*N);
x = a * x;          
end