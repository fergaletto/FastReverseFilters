function [x0,pc,mseE] = MGD(x,g,maxIter,Beta,alp,MODE)

%Momentum

%input:
%x -- original image used to calculate psnr in each iteration
%g -- blackbox filter, defined as an inline function
%maxIter -- max number of iterations
%Beta -- parameter of momentum, can set it to 0.99
%alp -- parameter of momentum, can set it to 1
%MODE -- 1 = T, 2 = TDA, other = P

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index
%mseE -- mean square difference of b and g(x0) as a function of iteration
%        index, an indication of quality in each iteration
pc = zeros(maxIter+1,1);
mseE = zeros(maxIter,1);
maxX = max(x(:));
b = g(x);
pc(1) = psnr(b,x,maxX);

% Set method's parameters and Initialization
x0 = b;
v = 0;
%Beta = 0.99;
%alp = 1;

switch MODE
    case 1 %T-method
        disp('T-method')
        
        for k = 1 : maxIter
            h = b - g(x0);
            v = Beta.*v + alp.*h;
            x0 = x0 + v;
            pc(k+1) = psnr(x0,x, maxX);
            mseE(k) = mean(h(:).^2);
        end
    case 2 %TDA
        disp('TDA-method')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            h = g(x0+e) - y;
            v = Beta.*v + alp.*h;
            x0 = x0 + v;
            pc(k+1) = psnr(x0,x,maxX);
            mseE(k) = mean(e(:).^2);
        end
    case 3 %p-method
        disp('P-method with norm')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            d = g(x0+e) - g(x0-e);
            h = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
            v = Beta.*v + alp.*h;
            x0 = x0 + v;
            pc(k+1) = psnr(x0,x,maxX);
            mseE(k) = mean(e(:).^2);
        end
    otherwise
        disp('P-method without norm')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            h = g(x0+e) - g(x0-e);

            v = Beta.*v + alp.*h/2;
            x0 = x0 + v;
            pc(k+1) = psnr(x0,x,maxX);
            mseE(k) = mean(e(:).^2);
        end
    end
end


function r = norm3(x,p)
 % TO calculate the norm, can handle both color and grayscale images. 
    r = [];
    for n=1: size(x, 3)
        y = norm(x(:,:,n),p);
        r=cat(3, r,y);
    end

end


