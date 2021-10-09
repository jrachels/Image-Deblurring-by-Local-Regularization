function bmat = genb( N, n, ell, sb )
%GENB Generates B matrix

bmat = zeros((2*n)^2,(2*ell*n)^2);
for px=1:2*ell*n
    for py=1:2*ell*n
        if px<(ell-1)*n+1 || px>(ell+1)*n || py<(ell-1)*n+1 || py>(ell+1)*n
            for lamx=1:2*n
                for lamy=1:2*n
                    bmat(lamx+2*n*(lamy-1),px+2*ell*n*(py-1)) = exp(-((px+(ell+1)*n-lamx)^2+(py+(ell+1)*n-lamy)^2)/(2*N^2*sb^2))/(2*pi*N^2*sb^2);
                end
            end
        end
    end
end

end

