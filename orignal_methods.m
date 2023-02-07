function [x1,pc,mseE] = orignal_methods(x,b,g,alpha,maxIter, MODE)

% This code implements the original methods for reverse image filering
% called T-method, TDA-method and P-method. 

%input:
%x -- original image used to calculate psnr in each iteration
%b -- image to be reversed
%g -- blackbox filter, defined as an inline function
%alpha -- Increment multiplier. 
% maxIter -- number of iterations. 
%MODE -- 1 = T, 2 = TDA, 3 = P, other = P without norm. 

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index
%mseE -- mean square difference of b and g(x0) as a function of iteration

maxX = max(x(:));
x0 = b;
pc = zeros(maxIter,1);
mseE = pc;

pc(1) = psnr(b,x,maxX);
n = 2;
switch MODE
    case 1  %T-method
        disp('T-method')
            for k = 1 : maxIter           
                y = g(x0);
                e = b - y;
                mseE(k+1) = mean(e(:).^2);
                x1 = x0 + alpha*e;
                pc(k+1) = psnr(x1,x,maxX);
                x0 = x1;
            end
    case 2
        disp('TDA-method')
        for k=1:maxIter
            
            y = g(x0);
            e = b - y;
            mseE(k+1) = mean(e(:).^2);
            Ax = g(x0+e)-y;
            x1 = x0 + alpha*Ax;
            pc(k+1) = psnr(x1,x,maxX);
            x0 = x1;          
        end
    case 3%p-method
        disp('P-method with norm')
        
        for k=1:maxIter                
            y = g(x0);
            e = b - y;
            mseE(k+1) = mean(e(:).^2);
            d = g(x0+e) - g(x0-e);  
            Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
            x1 = x0 + alpha*Ax;
            pc(k+1) = psnr(x1,x,maxX);
            x0 = x1;
        end
    otherwise
        disp('P-method without norm')
        
        for k=1:maxIter                
            y = g(x0);
            e = b - y;
            mseE(k+1) = mean(e(:).^2);
            d = g(x0+e) - g(x0-e);  
            Ax = d/2;
            x1 = x0 + alpha*Ax;
            pc(k+1) = psnr(x1,x,maxX);
            x0 = x1;
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

