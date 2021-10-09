function normmat = gennormmat_p( A, Siginv )
%GENNORMMAT_P Summary of this function goes here
%   Detailed explanation goes here

normmat = A.'*Siginv*A;

end

