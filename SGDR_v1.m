function [x1,pc,mseE] = SGDR_v1(x,b,g,T,nmin, nmax, MODE)

% SGDR for T-method
% @article{loshchilov2016sgdr,
%   title={Sgdr: Stochastic gradient descent with warm restarts},
%   author={Loshchilov, Ilya and Hutter, Frank},
%   journal={arXiv preprint arXiv:1608.03983},
%   year={2016}
% }

%input:
%x -- original image used to calculate psnr in each iteration
%b -- image to be reversed
%g -- blackbox filter, defined as an inline function
%T -- number of iterations until restart. 
%     if T is a scalar, 500 iteration is used, number of period is 500/T.
%     if T is a vector, it specifies the period in each iteration
%MODE -- 1 = T, 2 = TDA, 3 = P-method with norm, other P-method without
%norm. 


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
                omega_k = nmin+0.5*(nmax-nmin)*(1+cos(k/T(numPeriod)*pi));
                
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
    
                %Ax = g(x0+e)-y;
                x1 = x0 + omega_k*e;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
    
            end
        end
    case 2 %TDA
        disp('TDA-method')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = nmin+0.5*(nmax-nmin)*(1+cos(k/T(numPeriod)*pi));
                
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
    
                Ax = g(x0+e)-y;
                x1 = x0 + omega_k*Ax;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;
            end
        end
        
    case 3 %p-method w norm
        disp('P-method with norm')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = nmin+0.5*(nmax-nmin)*(1+cos(k/T(numPeriod)*pi));
                
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                d = g(x0+e) - g(x0-e);  
                Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
                %Ax = d/2;
                x1 = x0 + omega_k*Ax;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;

            end
        end
    otherwise
        disp('P-method without norm')
        for numPeriod = 1 : maxPeriod
            for k = 0 : T(numPeriod)-1
                omega_k = nmin+0.5*(nmax-nmin)*(1+cos(k/T(numPeriod)*pi));
                
                y = g(x0);
                e = b - y;
                mseE(n) = mean(e(:).^2);
                d = g(x0+e) - g(x0-e);  
                %Ax = 4*norm(h,2).^2./(norm(d,2)+eps).^2.*d/2;
                Ax = d/2;
                x1 = x0 + omega_k*Ax;
                pc(n) = psnr(x1,x,maxX);
                x0 = x1;
                n = n + 1;

            end
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
