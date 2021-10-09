function [ g, g_ext ] = bttbmvpcbcedfft( t, f )
%BTTBMVPCBCEDFFT Block Toeplitz with Toeplitz Blocks Matrix-Vector Computation by Block Circulant Extension and Discrete Fast Fourier Transforms
%   This is based on Algorithm 5.2.1 from Curtis R. Vogel, "Computational
%   Methods for Inverse Problems", 2002, page 74.
% INPUTS
% t: For a BTTB matrix A with size(A)=[N^2,N^2] and block size [N,N], this
%    is the size [2*N-1,2*N-1] array such that the ith diagonal of the jth
%    block of A is filled with the value t(i,j).
% f: Matrix of size [N,N], which leads to a vector f(:) of size [N^2,1]
%    whose product with A we wish to calculate.
% OUTPUTS
% g: Matrix of size [N,N], related to the BTTB matrix-vector product A*f(:)
%    by g = reshape(A*f(:), [N,N]).

N = size(f,1);
% If dimensions are wrong, throw error before doing calculations
assert(N == size(f,2));
assert(2*N-1 == size(t,1) && 2*N-1 == size(t,2));

t_tild = padarray(t,[1,1],'pre');
c_ext = [t_tild(N+1:2*N,N+1:2*N),t_tild(N+1:2*N,1:N);t_tild(1:N,N+1:2*N),t_tild(1:N,1:N)];
f_ext = padarray(f,[N,N],'post');
g_ext = ifft2(fft2(c_ext).*fft2(f_ext));
g = g_ext(1:N,1:N);


end

