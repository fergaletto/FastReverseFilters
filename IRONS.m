function [x0, pc] = IRONS(x, g, maxIter, MODE)
% This code implements Irons's Algorithm from the paper: "A version of the
% Aitken acceleration for computer iteration" 1969.
% Source: https://pdf.sciencedirectassets.com/271503/1-s2.0-S0898122115X00204/1-s2.0-S0898122115004046/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEO7%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQCAJd%2FtfloY7Ccgiq96N7eZF6cnmFSr279c6DVfzmgxiQIhAOuLhBfQMZKNPOxfhQ2ibL6BlhnjQSRogUTZc4KYJe6qKvoDCGYQBBoMMDU5MDAzNTQ2ODY1IgyQv1PYJTcSGS4bnSEq1wO4cjChgQjnRI1x91ZP0FV2FT2L9wKaYUufJx%2BdyIF1wsM%2FgoGenZKclHdMG8awYPT05ci2ExjbX3G9pObw8S%2F3SsAtB9ghWe8tKfcB5GuNSk8I9WASYGXhqFhlCwz5wqMk3Om%2BI%2F3sCg9L5oc%2BL6VC2ZKz5zEWNjHt0Wsv7wTF0kwHxOX2MkDSSh5%2BrGySWLm5wU%2FtKk4rgtygZDRiDeRmkOI98z04kc1s8u1MKEkk0ieeoBcCMhpsgCU0V1GlSi5lgSpQWPFjxpgtODmIAfFfAQ97zQZLu%2FC93Rp7uClmAr6%2F0Jr6%2FGOw%2FDHE41%2FrVog0q659V97UY8z3eaQkTcE7k7VTOrL0P3OZbx2PjILTKlDBGs22ITgihQOg3F%2Fh5e3WtM0Zzs7t1RWztO1X2XeOCPuhHXi5oboNjLK5LwfMiO%2BOHkci4ZJ8bWJ%2FGFhvxkvOfWe7fbSr3AnLfuRyP0x2J5SD7HyTrdcUiTEFLLkQ0fsUZ3PtPTy7tc0gdCSlU1W7z%2B1ESZXqPV2yrCyFV9hOk3ANQkGaY4sJhMFCKD%2BQQ4mX1UEqZ88u2lq8ROJsIFoIbLHmWKWOZOVRd2wefLQ8S%2F7NeyluOh44wdbY7A%2FfrFIcA6XLjjww796%2BkQY6pAGUdCy%2BOatKIwMTk0bqKrCbkFYyitPWNL9yfbR4VBcGJsnSDaRUgssLOFvIf1QoN0yNXjYu8v%2F7OR0krRxaL1QbdVYm5TKvpCm8S6KbPGmBpXpQgx%2Bg3VzdUXQLm448A%2FlxuDYDnHVu%2BngglSsMoBy1v8IORDxrGUtKKLWk7c5hxpXyKSSwOU83pLE5ejzXB1suFeVehdtNlOCvLS%2BRPbfKJCI5UA%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220314T223729Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYS2IXMOJZ%2F20220314%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=5a838ad3fb0c29b36cff7aafb92cd9428a6198f85d0b2955ba6dff5773314224&hash=48747a3fdb8ac7da94fa9e9e543403392507ecd507fb322eb455e776f43a533e&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0898122115004046&tid=spdf-2030f236-8366-4161-abe8-cc7775329a36&sid=af3075bb457e59402b7a433250812678d226gxrqa&type=client&download=true&ua=5906045752065d575e&rr=6ec07e829cdf3772
% page 2212
%
%input:
%x -- original image used to calculate psnr in each iteration
%g -- blackbox filter, defined as an inline function
% maxIter -- number of iterations. 
%MODE -- 1 = T, 2 = TDA, 3 = P w/norm , other = P w/o norm

%output
%x1 -- result at the last iteration
%pc -- psnr as a function of iteration index

b = g(x);
pc = zeros(maxIter,1);
maxX = max(x(:));
pc(1) = psnr(b,x,maxX);

% Select the function f (depends on the method to reverse)

switch MODE
    case 1 %T-method
        disp('T-method')
        f = @(in) in + b - g(in); 
    case 2 %TDA
        disp('TDA-method')
        f = @(in) TDAmethod(in, b, g); % See function definition below
    case 3 %p-method
        disp('P-method w/norm')
        f = @(in) Pmethod(in, b, g); % See function definition below  
    case 4
        disp('P-method w/o norm') % 
        f = @(in) P2method(in, b, g); % See function definition below
 
end

%% The method. 
    x0=b; % initialize the variable. 
    for k = 1:maxIter
        gx = f(x0);
        ggx = f(gx);
        Ax = gx-x0;
        Agx = ggx-gx;
        A2x = Agx-Ax;

        x0 = ggx -Agx.*A2x.*Agx./norm3(A2x,2).^2;   
        pc(k+1) = psnr(x0,x,maxX);
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


