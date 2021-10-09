function A = gena( N, n, sb )
%GENA Generates the BTTB matrix bigA, as well as the Toeplitz matrix A.
%   This function is separate since many unneeded smaller Toeplitz
%   matrices, as well as many vectors, are created in the process of
%   creating bigA.

%bigQ = cell(N,1);
Q = cell(2*n,1);
for d_ly = 0:(2*n-1)
    %bigrow = exp(-((0:(N-1)).^2+d_ly^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);
    row = exp(-((0:(2*n-1)).^2+d_ly^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);
    %bigQ{d_ly+1} = toeplitz(bigrow);
    %if d_ly < 2*n
    Q{d_ly+1} = toeplitz(row);
    %end
end
%bigA = cell2mat(bigQ(toeplitz(1:N)));
A = cell2mat(Q(toeplitz(1:2*n)));

end

