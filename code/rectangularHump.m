function u = rectangularHump( N, a, b )
%RECTANGULARHUMP Generates a rectangular hump in two dimensions
%   a,b are the beginning,end of hump in either axial direction

[X,Y] = meshgrid(1/N:1/N:1);
u = (a<X).*(X<b).*(a<Y).*(Y<b);

end

