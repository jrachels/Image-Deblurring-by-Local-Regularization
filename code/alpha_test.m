function [ u, f_d, recons ] = alpha_test(  )
%ALPHA_TEST Summary of this function goes here
%   Detailed explanation goes here

N = 512;
n = 16;
%alpha = 0.000015;
ell = 1;

u = 10*rectangularHump(N,0.2,0.8) + 220*rectangularHump(N,0.22,0.78) - 80*rectangularHump(N,0.32,0.68);
%u = 10 + 222*verticals(N,3);

sb = 0.045;
sn = 0;

[bigt,t] = gent(N,n,sb);

A = gena(N,n,sb);
bmat = genb(N,n,ell,sb);

f = bttbmvpcbcedfft(bigt,u);
f_d = poissrnd(f) + sn*randn(N);

figure(1)
imagesc(u,[0,255]);
drawnow
colormap hot
figure(2)
imagesc(f_d,[0,255]);
colormap hot
drawnow

c = zeros((N-2*n)^2,(2*n)^2);
for lx = 1:2*n
    for ly = 1:2*n
        for x = 1:N-2*n
            for y = 1:N-2*n
                c(x+(N-2*n)*(y-1),lx+2*n*(ly-1)) = f_d(x+lx,y+ly);
            end
        end
    end
end

recons = cell(10,2);

alpha = [0.0003] ;

for number=1:5
    alpha = alpha/1.27;
    [u_d,time]=locregp(N,n,alpha,sn,bigt,t,A,bmat,f_d,c,number);
    figure(number+2)
    imagesc(u_d,[0,255]);
    colormap hot
    drawnow
    recons{number,1} = u_d;
    recons{number,2} = alpha;
    fprintf('time = %d\n', time);
end


end

