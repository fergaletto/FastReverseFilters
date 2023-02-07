function [x1,pc,mseE] = cheby_v2(x,b,g,T, alpha,MODE)

%Chebyshev iteration for T-method
% "Chebyshev Inertial Iteration for Accelerating Fixed-Point Iterations"
%Tadashi Wadayama? and Satoshi Takabe, arXiv:2001.03280v1
%Fig.10 and section V.D

%input:
%x -- original image used to calculate psnr in each iteration
%b -- image to be reversed
%g -- blackbox filter, defined as an inline function
%T -- length one period for the chebyshev coefficient omega
%     if T is a scalar, 500 iteration is used, number of period is 500/T.
%     if T is a vector, it specifies the period in each iteration
%MODE -- 1 = T, 2 = TDA, other = P

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index
%mseE -- mean square difference of b and g(x0) as a function of iteration
%        index, an indication of quality in each iteration
maxX = max(x(:));
maxPeriod = length(T);
if maxPeriod == 1
    maxIter = 200;
    maxPeriod = ceil(maxIter/T); %default number of period
    T = T*ones(maxPeriod,1);
end
maxIter = sum(T);
x0 = b;
pc = zeros(maxIter,1);
mseE = pc;


pc(1) = psnr(b,x,maxX);
n = 2;

switch MODE
    case 1 %T-method
        disp('T-method')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = min(alpha,2/(1 + cos(pi*(k+0.5)/T(numPeriod))));
                %omega_k = 1+180*exp(-(T-1-k)^2/1.2);
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                %this is the T update
                x1 = x0 + omega_k*e;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
            end
            alpha = 0.75*alpha;
        end
    case 2 %TDA
        disp('TDA-method')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = min(alpha,2/(1 + cos(pi*(k+0.5)/T(numPeriod))));
                %omega_k = 1+180*exp(-(T-1-k)^2/1.2);
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                %this is the TDA update
                x1 = x0 + omega_k*(g(x0+e)-y);
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
            end
            alpha = 0.75*alpha;
        end
        
    case 3 %p-method
        disp('P-method with norm')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = min(alpha,2/(1 + cos(pi*(k+0.5)/T(numPeriod))));
                %omega_k = 1+180*exp(-(T-1-k)^2/1.2);
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                %this is the P update
                d = g(x0+e) - g(x0-e);  
                Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
                x1 = x0 + omega_k*Ax;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
            end
            alpha = 0.75*alpha;
        end
    otherwise
        disp('P-method without norm')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = min(alpha,2/(1 + cos(pi*(k+0.5)/T(numPeriod))));
                %omega_k = 1+180*exp(-(T-1-k)^2/1.2);
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                %this is the P update
                 %d = g(x0+e) - g(x0-e);  
                 %Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
                Ax = (g(x0+e) - g(x0-e))/2;
                x1 = x0 + omega_k*Ax;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
            end
            alpha = 0.75*alpha;
        end
    end

%end of function
end

function r = norm3(x,p)
 % TO calculate the norm, can handle both color and grayscale images. 
    r = [];
    for n=1: size(x, 3)
        y = norm(x(:,:,n),p);
        r=cat(3, r,y);
    end

end