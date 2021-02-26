% Copyright (c) 2021, KDDI Research, Inc. and KDDI Corp. All rights reserved.

function out = invcramer(in)

M = size(in,1);
N = size(in,2);
adj = zeros(M,N); 
for m = 1:M
    for n = 1:N
        adj(m,n) = (-1)^(m+n)*det(in([1:m-1 m+1:M],[1:n-1 n+1:N]));
    end
end
out = adj.'/det(in);

end
