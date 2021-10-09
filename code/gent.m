function [ bigt, t ] = gent( N, n, sb )
%GENT Generates the BTTB array for the matrix bigA
%   For the BTTB matrix bigA with size(bigA)=[N^2,N^2] and block size
%   [N,N], this is the size [2*N-1,2*N-1] array such that the ith diagonal
%   of the jth block of A is filled with the value t(i,j).
%   Similarly, for the BTTB matrix A with size(A)=[2*n^2,4*n^2] and block
%   size [2*n,2*n], this is the [4*n-1,4*n-1] size array with the same
%   property.
%   The matrices bigA and A themselves do not need to be created.

[DLX,DLY] = meshgrid((1-N):(N-1));
%[DPX,DPY] = meshgrid((1-2*ell*n):(2*ell*n-1));
[DLamX,DLamY] = meshgrid((1-2*n):(2*n-1));
bigt = exp(-(DLX.^2+DLY.^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);
%medt = exp(-(DPX.^2+DPY.^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);
t = exp(-(DLamX.^2+DLamY.^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);

end

