function [x0,pc,mseE] = Adam(x,g,maxIter,MODE)

%Nestorov

%input:
%x -- original image used to calculate psnr in each iteration
%g -- blackbox filter, defined as an inline function
%maxIter -- max number of iterations
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
s = 0;
r = 0;
%Beta = 0.99;
%alp = 1;
epsilon = 0.01;
rho1 = 0.9;
rho2 = 0.999;


switch MODE
    case 1 %T-method
        disp('T-method')
        
        for k = 1 : maxIter
            %gradient
            h = g(x0) - b; %for proper gradient descent
            %update 1st moment
            s = rho1*s + (1-rho1)*h;
            %update 2nd moment
            r = rho2*r + (1-rho2)*(h.^2);
            %correct bias 1st moment
            sHat = s/(1-rho1);
            %correct bias 2nd moment
            rHat = r/(1-rho2);
            %apply update
            x0 = x0 - epsilon*(sHat./(sqrt(rHat) + 1e-8));
            
            pc(k+1) = psnr(x0,x, maxX);
            mseE(k) = mean(h(:).^2);
        end
    case 2 %TDA
        disp('TDA-method')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            %gradient
            h = y - g(x0+ e); %for proper gradient descent
            %update 1st moment
            s = rho1*s + (1-rho1)*h;
            %update 2nd moment
            r = rho2*r + (1-rho2)*(h.^2);
            %correct bias 1st moment
            sHat = s/(1-rho1);
            %correct bias 2nd moment
            rHat = r/(1-rho2);
            %apply update
            x0 = x0 - epsilon*(sHat./(sqrt(rHat) + 1e-8));
            pc(k+1) = psnr(x0,x,maxX);
            mseE(k) = mean(e(:).^2);
        end
    case 3 
        disp('P-method with norm')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            %gradient
            d = g(x0 - e) - g(x0 + e); %for proper gradient descent
            h = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;

            %update 1st moment
            s = rho1*s + (1-rho1)*h;
            %update 2nd moment
            r = rho2*r + (1-rho2)*(h.^2);
            %correct bias 1st moment
            sHat = s/(1-rho1);
            %correct bias 2nd moment
            rHat = r/(1-rho2);
            %apply update
            x0 = x0 - epsilon*(sHat./(sqrt(rHat) + 1e-8));
            pc(k+1) = psnr(x0,x,maxX);
            mseE(k) = mean(e(:).^2);
        end
    otherwise %p-method
        disp('P-method without norm')
        
        for k = 1 : maxIter
            y = g(x0);
            e = b - y;
            %gradient
            h = (g(x0 - e) - g(x0 + e))/2; %for proper gradient descent
            %update 1st moment
            s = rho1*s + (1-rho1)*h;
            %update 2nd moment
            r = rho2*r + (1-rho2)*(h.^2);
            %correct bias 1st moment
            sHat = s/(1-rho1);
            %correct bias 2nd moment
            rHat = r/(1-rho2);
            %apply update
            x0 = x0 - epsilon*(sHat./(sqrt(rHat) + 1e-8));
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




