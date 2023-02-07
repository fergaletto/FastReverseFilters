function [x1,pc,mseE] = Aitken_v3(x,b,g,maxIter, MODE)

% This code implements Aitken's delta-squared process to improve the
% convergence of three reverse image filtering methods. 

%input:
%x -- original image used to calculate psnr in each iteration
%b -- image to be reversed
%g -- blackbox filter, defined as an inline function
% maxIter -- number of iterations. 
%MODE -- 1 = T, 2 = TDA, other = P

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index
%mseE -- mean square difference of b and g(x0) as a function of iteration

maxX = max(x(:));
x0 = b;
pc = zeros(maxIter,1);
mseE = pc;

alpha = 1; % increment multiplier. 

pc(1) = psnr(b,x,maxX);
epsilon = 0.01;
switch MODE
    case 1 %T-method
    disp('T-method')
    f = @(in) in + b - g(in); 
    case 2 %TDA
        disp('TDA-method')
        f = @(in) in + alpha*(g(in+b-g(in))-g(in));    
    case 3 %p-method
            disp('P-method w/norm')
            ee = @(in) b-g(in);
            d = @(in,ei) g(in+ei) - g(in-ei); 
            P = @(ein, din) 4*norm3(ein,2).^2./(norm3(din,2)+eps).^2.*din/2; %Original Version
            f = @(in) in+P(ee(in), d(in, ee(in)));     
    case 4
            disp('P-method w/o norm') % 
            ee = @(in) b-g(in);
            d = @(in,ei) g(in+ei) - g(in-ei); 
            P = @(ein, din) din/2; % Dennis version Performs as TDA method. 
            f = @(in) in+P(ee(in), d(in, ee(in)));  
 
end

for k = 1 : maxIter 

    
    x1 = f(x0);
    x2 = f(x1);
    aitkenX = x2;
    denom = x2 - 2*x1 + x0;
    idx = find(abs(denom)> epsilon);
    
    aitkenX(idx) = x2(idx) - ((x2(idx)-x1(idx)).^2)./denom(idx);
    x0 = max(min(1,aitkenX),0);
    
    e = b - g(x0);
    mseE(k+1) = mean(e(:).^2);
    pc(k+1) = psnr(x0,x,maxX);
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

