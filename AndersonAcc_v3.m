function [ua, pc] = AndersonAcc_v3(x,g,maxIter, MODE)

%AndersonAccKelley4AnyFilter -- inverse filtering
%Anderson implementation of page 23 of Kelley's talk
%faster version by reuse results

%input:
%x -- original image used to calculate psnr in each iteration
%g -- blackbox filter, defined as an inline function
% maxIter -- number of iterations. 
%MODE -- 1 = T, 2 = TDA, 3 = P other = new P

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index

maxX = max(x(:));
b=g(x);

pc = zeros(maxIter,1);

% Select the method 
switch MODE 
    case 1 %T-method
    disp('T-method')
    F = @(in)  b - g(in); 
    case 2 %TDA
        disp('TDA-method')
        F = @(in) (g(in+b-g(in))-g(in));    
    case 3 %p-method
            disp('P-method w/norm')
             F = @(in) Pmethod(in, b, g);   
    case 4
            disp('P-method w/o norm') % 
            F = @(in) P2method(in, b, g); 
 
end


%% Apply Anderson Acc. 

%G and F function as in Kelley's presentation
G = @(in) in + F(in);
 % because F = @(x) G(x) - x;

%initialization
u0 = b;         pc(1) = psnr(b,x,maxX);
u1 = G(u0);     pc(2) = psnr(u1,x,maxX);
u2 = G(u1);     pc(3) = psnr(u2,x,maxX);


Gu0 = u1;
Gu1 = u2;
Gu2 = G(u2);

Fu1 = Gu1(:) - u1(:);
v0 = Fu1(:) - (Gu0(:) - u0(:));
v2 = Gu2(:) - u2(:); %same as Fu2
v1 = v2(:) - Fu1(:);


for k = 3 : maxIter 
%step 1   
    D = [v0'*v0 v0'*v1;v0'*v1 v1'*v1];
    den = [v0'*v2;v1'*v2];
    a = D\den;
    %an interesting note
    %in theory a = D\d  is equvalent to a = [v0 v1]\v2
    %but this implementation is 8 times slower!
    
    %Step 2
    ua = Gu2 - a(1)*(Gu1 - Gu0) - a(2)*(Gu2 - Gu1);
    %update variable and reuse results
    v0 = v1;
    Fu1 = v2; % because Fu2 = v2
    Gu0 = Gu1;
    Gu1 = Gu2;
    u2 = ua;
    Gu2 = G(u2); %call the filter once
    v2 = Gu2(:) - u2(:);
    v1 = v2 - Fu1;

    pc(k+1) = psnr(ua,x,maxX);
end     
%end of function
end


%% Auxiliary functions. 

function r = norm3(x,p)
 % TO calculate the norm, can handle both color and grayscale images.
    r = [];
    for n=1: size(x, 3)
        y = norm(x(:,:,n),p);
        r=cat(3, r,y);
    end

end


function Ax = Pmethod(x0, b, g)
    % This function does one iteration of the P method and returns the
    % increment. 
    % X0 previous iteration value. 
    % b input blurred image. 
    % g Filter to reverse (Function handler).
            y = g(x0);
            e = b - y;
            d = g(x0+e) - g(x0-e);  
            Ax = 4*norm3(e,2).^2./(norm3(d,2)+eps).^2.*d/2;

end

function Ax = P2method(x0, b, g)
    % This function does one iteration of the new P method and returns the
    % increment. 
    % X0 previous iteration value. 
    % b input blurred image. 
    % g Filter to reverse (Function handler).
            y = g(x0);
            e = b - y;
            d = g(x0+e) - g(x0-e);  
            Ax = d/2;

end

