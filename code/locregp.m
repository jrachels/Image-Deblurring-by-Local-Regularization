function [ u_d, time ] = locregp( N, n, alpha, sn, bigt, t, A, bmat, f_d, c, number )
%LOCREG_P 2D local regularization with Poisson or mixed noise

% N = 256;
% n = 4;
% % Slightly smaller than this, probably
% alpha = 0.000015;
% ell = 1;
% 
% u = 10*rectangularHump(N,0.2,0.8) + 220*rectangularHump(N,0.22,0.78) - 80*rectangularHump(N,0.32,0.68);
% %u = 10 + 222*verticals(N,3);
% 
% sb = 0.045;
% sn = 0;
% 
% [bigt,t] = gent(N,n,sb);
% 
% A = gena(N,n,sb);
% bmat = genb(N,n,ell,sb);
% 
% f = bttbmvpcbcedfft(bigt,u);
% f_d = poissrnd(f) + sn*randn(N);
% 
% figure(1)
% imagesc(u,[0,255]);
% drawnow
% colormap hot
% figure(2)
% imagesc(f_d,[0,255]);
% colormap hot
% drawnow
% 
% c = zeros((N-2*n)^2,(2*n)^2);
% for lx = 1:2*n
%     for ly = 1:2*n
%         for x = 1:N-2*n
%             for y = 1:N-2*n
%                 c(x+(N-2*n)*(y-1),lx+2*n*(ly-1)) = f_d(x+lx,y+ly);
%             end
%         end
%     end
% end
u_d = padarray(reshape(c(:,n*(2*n-1)),[N-2*n,N-2*n]),[n,n]);
% 
% figure(3)
% imagesc(u_d,[0,255]);
% colormap hot
% drawnow

% max_iter = 1;

eps = zeros(N-2*n);
% iter = 0;
% while iter==0 || iter<max_iter
    tic;
    eps = 0*eps;
    bigApadcbarold = bttbmvpcbcedfft(bigt,u_d);
    for x=1:N-2*n
        for y=1:N-2*n
            L_padcbarold = u_d(x+1:x+2*n,y+1:y+2*n);
            L_bigApadcbarold = bigApadcbarold(x+1:x+2*n,y+1:y+2*n);
            L_f = f_d(x+1:x+2*n,y+1:y+2*n);
            % padLeps = padarray(eps,[ell*n,ell*n]);
            padLeps = padarray(eps,[n,n]);
            % L_padLeps = padLeps(x+1:x+2*ell*n,y+1:y+2*ell*n);
            L_padLeps = padLeps(x+1:x+2*n,y+1:y+2*n);
            g = L_f - L_bigApadcbarold + bttbmvpcbcedfft(t,L_padcbarold) - reshape(bmat*L_padLeps(:),[2*n,2*n]);
            Siginv = diag((L_bigApadcbarold(:) + sn^2).^(-1));
            mat = gennormmat_p(A,Siginv) + alpha*eye((2*n)^2);
            ATSiginvg = bttbmvpcbcedfft(t,reshape(Siginv*g(:),[2*n,2*n]));
            c(x+(N-2*n)*(y-1),:) = transpose(mat\(ATSiginvg(:)));
            eps(x,y) = c(x+(N-2*n)*(y-1),n*(2*n-1)) - u_d(x+n,y+n);
        end
    end
    u_d = padarray(reshape(c(:,n*(2*n-1)),[N-2*n,N-2*n]),[n,n]);
    time = toc;
%     figure(number)
%     imagesc(u_d,[0,255]);
%     colormap hot
%     drawnow
%     iter = iter + 1;
% end


end

