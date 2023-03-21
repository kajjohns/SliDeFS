function M0 = Kostrov(X,Y,Exx,Exy,Eyy,H,mu,tri_areas)



% First, determine which function to call
if nargin < 8, tri = false; else tri = true; end % assume quadrilateral elements by default
if nargin < 7, mu = 30e9; end % assume default shear modulus


eps = [Exx(:)'; Eyy(:)'; Exy(:)']; 

if tri
    A = tri_areas;
else
    A = area_quad(X, Y);
end

epsmax = getEpsMax(eps);
maxGam = getMaxShearStrain(eps);

M0 = getM0(epsmax, A, H, mu);
%M0 = getM0(maxGam, A, H, mu);

end


function A = area_quad(X,Y)

% compute grid spacing and the area of each gridded block in meters 
h1 = abs(X(1,1) - X(2,1)); 
if h1==0
    h1 = abs(X(1,1) - X(1,2));
    if h1==0
        error('strain_quad: h1 is zero')
    end
end
h2 = abs(Y(1,1) - Y(2,1));
if h2 == 0
    h2 = abs(Y(1,1) - Y(1,2));
    if h2==0
        error('strain_quad: h2 is zero')
    end
end
h1 = h1*1000; 
h2 = h2*1000; 
A = h1*h2; 

end

function [epsmax] = getEpsMax(epsall)
% Compute the maximum principal strain
dV1dx1  = epsall(1,:); 
dV2dx2  = epsall(2,:); 
e12     = epsall(3,:); 

% compute max(|e_1|, |e_2|, |e_1 + e_2|) for the moment calculation
T = dV1dx1(:) + dV2dx2(:); 
D = dV1dx1(:).*dV2dx2(:) - (e12(:).^2); 

lam1 = T/2 + sqrt((T.^2)/4 - D);
lam2 = T/2 - sqrt((T.^2)/4 - D); 

L = [abs(lam1'); abs(lam2'); abs(lam1' + lam2')]; 
%L = [abs(lam1'); abs(lam2')]; 
epsmax = max(L,[],1); 
epsmax(isnan(epsmax))=0;
end

function [M0] = getM0(epsmax, A, H, mu)
% getM0: convert maximum principal strain rate to moment using the Savage &
% Simpson 1994 method.
% Inputs: H (km), mu (Pa), default 30e9
    if nargin < 3, mu=30e9; end
    % convert H to meters
    H = H*1000; 
    if length(A)==1
        M0 = 2*sum(mu*H*A*epsmax); 
    else
        M0 = 2*mu*H*(A(:)'*epsmax(:));
    end
end

function [maxGam] = getMaxShearStrain(eps)
    maxGam = sqrt(((eps(1,:) - eps(2,:))/2).^2 + eps(3,:).^2);
end

