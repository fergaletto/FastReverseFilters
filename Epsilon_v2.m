function [x0,pc] = Epsilon_v2(x,g,maxIter, MODE)

% This code implements epsilon algorithm for vector sequences. 
% Souce xIterativeresidual-basedvectormethodstoacceleratefixedpointiterations
% Eq. 21

%input:
%x -- original image used to calculate psnr in each iteration
%b -- image to be reversed
%g -- blackbox filter, defined as an inline function
% maxIter -- number of iterations. 
%MODE -- 1 = T, 2 = TDA, 3 = P w/norm , other = P w/o norm

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index
%mseE -- mean square difference of b and g(x0) as a function of iteration

maxX = max(x(:));
b = g(x);
x0 = b;
pc = zeros(maxIter,1);
mseE = pc;

alpha = 1; % increment multiplier. 

pc(1) = psnr(b,x,maxX);
switch MODE
    case 1 %T-method
        disp('T-method')
        f = @(in) in + b - g(in); 
    case 2 %TDA
        disp('TDA-method')
        f = @(in) TDAmethod(in, b, g);
    case 3 %p-method
        disp('P-method w/norm')
        f = @(in) Pmethod(in, b, g);   
    case 4
        disp('P-method w/o norm') % 
        f = @(in) P2method(in, b, g); 
 
end
for k = 1 : maxIter 

        gx = f(x0);
        ggx = f(gx);
        Ax = gx-x0;
        Agx = ggx-gx;
        A2x = Agx-Ax;


    x0 = gx + (Agx.* norm3(Ax,2).^2 - Ax .* norm3(Agx,2).^2)./norm3(A2x,2).^2;
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
function x1 = Pmethod(x0, b, g)
    % This function does one iteration of the P method and returns the
    % increment. 
    % X0 previous iteration value. 
    % b input blurred image. 
    % g Filter to reverse (Function handler).
            y = g(x0);
            e = b - y;
            d = g(x0+e) - g(x0-e);  
            Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;
            x1 = x0 + Ax;


end

function x1 = P2method(x0, b, g)
    % This function does one iteration of the new P method and returns the
    % increment. 
    % X0 previous iteration value. 
    % b input blurred image. 
    % g Filter to reverse (Function handler).
            y = g(x0);
            e = b - y;
            d = g(x0+e) - g(x0-e);  
            Ax = d/2;
            x1 = x0 + Ax;

end

function x1 = TDAmethod(x0, b, g)
    % This function does one iteration of the new P method and returns the
    % increment. 
    % X0 previous iteration value. 
    % b input blurred image. 
    % g Filter to reverse (Function handler).

            y = g(x0);
            e = b - y;
            Ax = g(x0+e)-y;
            x1 = x0 + Ax;
             



end


