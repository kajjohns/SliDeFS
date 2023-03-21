function [displacement,gradient,Stress,flag]= DisloOkada3mode(m,coordinates,shear_m,poission_ratio)
%This function is to derive displacements, displacement gradients and stresses anywhere in a half space.
%It is based on solutions provided by Y. Okada (1992)
%in April, 1992. Bulletin of the Seismological Society of America. Vol.
%82, No.2 pp. 1018-1040

%auther: Wen-Jeng Huang
%date: Nov. 9th. 2006
%most recent modification:  Nov. 12th. 2006

      %parameters for DisloOkada3mode.m is equivalent to disloc3dR12
    %North_position and east_position refers to position of center of
    %   bottom edge of dislocation.
    %depth referes to depth to bottom edge of dislocation
    %m=[length width depth dip strike east_position north_position strike-slip dip-slip tensile]
    %strike-slip (positive): left-lateral;
    %dip-slip (positive): reverse fault; 
    %tensile (positive): opening:
    %[displacement]=DisloOkadaDisplacement(m,[x; y; z],G,Nu);
    % x,y,z: cooridnates
    % G: shear modulus
    % Nu: poission ratio
    %results: U(1,:) = displacement in east(+); U(2,:) = displacement in north(+); U(3,:) = displacement upward(+); 

    flag = [];
    
X = coordinates(1,:);
Y = coordinates(2,:);
Z = coordinates(3,:);


  L = m(1);
  W = m(2);
  D = m(3);
 angle = m(4); % dip in degree
 
 if D >= W*sin(angle/180*pi) && 1 == isequal(Z,-abs(Z))% A1
strikeAngle = m(5); %closewise is positive
 Xc = m(6); 
 Yc = m(7);
 slip_str = m(8);
 slip_dip = m(9);
 tensile = m(10);
 Gshear = shear_m;
 nu = poission_ratio;

 % unit transformation and constants
delta = angle/180 *pi;
angle_Str = -strikeAngle/180 * pi;
if nu == 0.5
    nu = 0.4999;
    youngs = 2 * Gshear * (1 +nu);
    lambda = nu *youngs / ( (1 + nu) *(1 - 2*nu));   
else
   youngs = 2 * Gshear * (1 +nu);
    lambda = nu *youngs / ( (1 + nu) *(1 - 2*nu));
end
    mu = Gshear;
c = D;
% coordinate transformation
x =  -sin(angle_Str) * (X - Xc) +  cos(angle_Str) * (Y - Yc ) +1/2 *L;
y =  -cos(angle_Str) * (X - Xc) -  sin(angle_Str) * (Y - Yc);
z = Z;
%%%%%%%%%%%%%%%%%
%integating 
d = c - z; 
p = y * cos(delta) + d * sin(delta);
   
            xi = [x;x;x-L;x-L]; % 
            eta = [p;p-W;p;p-W];%

   q = ones(4,1)*y * sin(delta) - ones(4,1)*d * cos(delta);
   alpha = (lambda + mu) / (lambda + 2* mu);
   %alpha = 2/3;
   Rsquare = xi.^2 + eta.^2 + q.^2;
   R = sqrt(Rsquare);
   y_ = eta * cos(delta) + q * sin(delta);
   d_ = eta * sin(delta) - q * cos(delta);
   c_ = d_ + ones(4,1)*z;
   %%%%%% For displacement
   X11 = 1./(R .* (R + xi));
   X32 = (2 .* R + xi )./ ( R.^3 .* ( R + xi ).^2 );
   X53 = (8*R.^2 + 9*R.*xi + 3*xi.^2) ./ (R.^5 .*( R + xi).^3 );
   Y11 = 1 ./ ( R.* ( R + eta ) );
   Y32 = (2 * R + eta ) ./ (R.^3 .* (R + eta).^2 );
   Y53 = (8*R.^2 + 9*R.*eta + 3*eta.^2) ./ (R.^5 .*( R + eta).^3 );
   h = q*cos(delta) - ones(4,1)*z;
   Z32 = sin(delta)./R.^3 - h.*Y32;
   Z53 = 3*sin(delta)./R.^5 - h.*Y53;
   Y0 = Y11 - xi.^2.*Y32;
   Z0 = Z32 - xi.^2.*Z53;
   %%%
   %%% Selecting a right root for theta is the most tricky problem
   qsign = sign(q);
   theta = atan2(xi .* eta, abs(q) .* R);
   theta = qsign.* theta;

   X = sqrt( xi.^2 + q.^2);
    if abs(cos(delta))< 0.000001
        I3 = 1/2 * ( eta ./ (R + d_) + y_.*q ./( (R + d_).^2 ) - log(R + eta) );
        I4 = 1/2 * ( xi.* y_ ./( (R + d_).^2) );
    else
    I3 = 1/cos(delta) .* y_ ./ ( R + d_ ) - 1 ./cos(delta).^2 .*( log(R + eta) - sin(delta) .* log(R + d_));
    I4 = sin(delta)./cos(delta) .* xi./(R + d_) + 2./(cos(delta).^2) * atan2(eta.*(X + q.*cos(delta) )+ ...
        X .* (R + X) * sin(delta), xi.* (R + X)* cos(delta));
    end
    I1 = -(xi./(R +d_))*cos(delta) - I4 *sin(delta);
    I2 = log(R + d_) + I3* sin(delta);
    %%% ji
    D11 = 1./(R.*(R + d_));
    if abs(cos(delta))< 0.000001
        K1 = (xi .*  q) ./ (R + d_) .* D11;
        K3 = sin(delta)./ (R + d_) .* (xi.^2 .* D11 - 1);
    else
         K1 = xi/cos(delta) .* (D11 - Y11*sin(delta));
         K3 = 1/cos(delta) .* (q .* Y11 - y_ .* D11);
    end
    K2 = 1./R + K3*sin(delta);
    K4 = xi .* Y11 * cos(delta) - K1*sin(delta);
    
    J5 = -(d_ + y_.^2 ./ (R +d_) ) .* D11;
    J2 = xi .* y_ ./ (R + d_) .* D11;
    if abs(cos(delta))< 0.000001
        J6 = - y_ ./ (R + d_).^2 .* (xi.^2 .* D11 - 1/2);
        J3 = - xi ./(R + d_).^2 .* (q.^2 .* D11 - 1/2);
    else
        J6 = 1./cos(delta).*(K3 - J5*sin(delta));
        J3 = 1./cos(delta).*(K1 - J2*sin(delta));
    end
    J1 = J5*cos(delta) - J6*sin(delta);
    J4 = -xi.*Y11 - J2*cos(delta) + J3*sin(delta);
    %%% ki
    E = sin(delta) ./ R  - y_ .* q ./ R.^3;
    F = d_ ./ R.^3 + xi.^2 .* Y32 * sin(delta);
    G = 2 * X11 * sin(delta) - y_.* q.* X32;
    H = d_ .* q .* X32 + xi .* q .* Y32 * sin(delta);
    P = cos(delta) ./ R.^3 + q.* Y32 * sin(delta);
    Q = 3 * c_ .* d_ ./ R.^5 - (  ones(4,1)*z.* Y32 + Z32 + Z0)* sin(delta);
    % li
    E_ = cos(delta)./R + d_ .* q ./ R.^3;
    F_ = y_ ./ R.^3 + xi.^2 .* Y32 * cos(delta);
    G_ = 2 * X11.* cos(delta) + d_.* q.* X32;
    H_ = y_.* q.* X32 + xi.* q .* Y32 * cos(delta);
    P_ = sin(delta) ./ R.^3 - q.* Y32 * cos(delta);
    Q_ = (3 * c_ .* y_ )./ R.^5 + q.* Y32 - (  ones(4,1)* z.* Y32 + Z32 + Z0) * cos(delta);
if slip_str ~= 0
    %%%%%% strike-slip
    %%% displacement
    %%uA
    Su1A = theta/2 + alpha / 2 * xi .* q .* Y11;
    Su2A = alpha / 2 * q ./ R;
    Su3A = (1 - alpha) / 2 * log( R + eta) - alpha / 2 * q.^2 .* Y11;
    %%
    %uB
    Su1B = - xi .* q .* Y11 - theta - (1 - alpha) / alpha .* I1 * sin(delta);
    Su2B = -q ./ R  +  (1 - alpha) / alpha * y_ ./ ( R + d_ ) * sin(delta);
    Su3B = q.^2 .* Y11 - (1 - alpha) / alpha *  I2 * sin(delta);
    
    %uC
    Su1C = (1 - alpha) * xi .* Y11 * cos(delta) - alpha * xi .* q .* Z32;
    Su2C = (1 - alpha) * ( cos(delta) ./ R + 2 * q .* Y11 * sin(delta) ) - alpha * c_ .* q ./ R.^3;
    Su3C = (1 - alpha) * q .* Y11 * cos(delta) - alpha * ( c_ .* eta ./ R.^3 - ones(4,1)*z .* Y11 + xi.^2 .* Z32);
    %%%%%% strike-slip
    %%% displacement gradient
    %jA
    Sj1A = -(1 -alpha)/2 *q.*Y11 - alpha/2*(xi.^2).*q.*Y32;
    Sj2A = - alpha/2 * xi.*q./R.^3;
    Sj3A = (1 -alpha)/2 * xi.*Y11 + alpha/2 * xi.*q.^2.*Y32;
    
    %jB
    Sj1B = xi.^2 .*q.* Y32 -(1-alpha)/alpha * J1* sin(delta);
    Sj2B = xi.*q./R.^3 - (1-alpha)/alpha * J2 *sin(delta);
    Sj3B = -xi.*q.^2.*Y32 -(1-alpha)/alpha * J3* sin(delta); 
    
    %jC
    Sj1C = (1 -alpha)* Y0*cos(delta) - alpha * q .* Z0;
    Sj2C = -(1 -alpha)* xi.*(cos(delta)./R.^3 + 2*q.*Y32 *sin(delta)) + alpha * (3 *  c_ .*xi .*q) ./ R.^5 ;
    Sj3C = -(1 - alpha) * xi .* q .* Y32 * cos(delta) + alpha * xi.* ( (3 * c_ .* eta)./R.^5 - ones(4,1)*z.* Y32 - Z32 -Z0 );
    
    %kA
    Sk1A = (1 -alpha)/2 * xi.* Y11* sin(delta) + d_ / 2 .* X11 + alpha/2 * xi .* F; 
    Sk2A = alpha /2 * E;
    Sk3A = (1 -alpha)/2 *( cos(delta)./R + q.* Y11* sin(delta) ) - alpha/2 * q .* F;
    
    %kB
    Sk1B =  -xi.* F - d_ .* X11 + (1 -alpha)/alpha .*( xi.* Y11 + J4) * sin(delta);
    Sk2B =  -E + (1 -alpha)/alpha * (1./R + J5).* sin(delta);
    Sk3B = q.* F - (1 -alpha) /alpha * (q.*Y11 - J6) * sin(delta);
    
    %kC
    Sk1C = -(1 - alpha).* xi.* P * cos(delta) - alpha.* xi.* Q;
    Sk2C = 2 * (1 -alpha) * (d_./R.^3 - Y0 * sin(delta) ) * sin(delta) - y_ ./ R.^3 * cos(delta) - alpha* ( (c_ + d_)./R.^3 *sin(delta) - eta./R.^3 - 3*c_.* y_.* q./R.^5);
    Sk3C = -(1 -alpha)* q./R.^3 + ( y_ ./ R.^3 - Y0*cos(delta))* sin(delta) + alpha * ( (c_ + d_)./ R.^3 * cos(delta) + 3 * c_ .* d_ .* q ./R.^5 - (Y0 * cos(delta) + q .* Z0)* sin(delta));
    
    %lA
    Sl1A = (1 -alpha)/2 * xi .* Y11 * cos(delta) + y_/2 .* X11 + alpha/2 * xi.* F_;
    Sl2A = alpha/2 * E_;
    Sl3A =  - (1 - alpha)/2 * (sin(delta)./R - q.* Y11 * cos(delta)) - alpha/2 * q .* F_;
    
    %lB
    Sl1B = -xi.* F_ - y_ .* X11 + (1 -alpha) / alpha * K1 * sin(delta);
    Sl2B = - E_ + (1 -alpha)/alpha * y_ .* D11 * sin(delta);
    Sl3B = q.* F_ + (1 -alpha)/alpha * K2 * sin(delta);
    
    %lC
    Sl1C = (1 - alpha) * xi .* P_ * cos(delta) - alpha * xi.* Q_;
    Sl2C = 2 * (1 -alpha) * ( y_ ./ R.^3 - Y0 * cos(delta) ) * sin(delta) + d_ ./ R.^3 * cos(delta) - alpha * ( (c_ + d_)./R.^3 * cos(delta) +3 *c_.*d_.*q./R.^5); 
    Sl3C = (y_./R.^3 - Y0 * cos(delta) ) * cos(delta) - alpha * ( (c_ + d_)./R.^3 * sin(delta) - 3*c_.*y_.*q./R.^5 - Y0* sin(delta)^2 + q.* Z0*cos(delta));
end  
if slip_dip ~= 0 % B2
    %%% %%% dip-slip
    %%%%%%% displacement
        %%uA
    Du1A = alpha / 2 * q ./ R; 
    Du2A = theta/ 2 + alpha /2 * eta .* q .* X11;
    Du3A = (1-alpha) /2 * log(R + xi) - alpha/2 * q.^2 .* X11;
    %%
    %uB
    Du1B = - q./R + (1 -alpha)/alpha * I3 * sin(delta) * cos(delta);
    Du2B = -eta .* q .* X11 - theta - (1 -alpha)/ alpha * xi ./ (R + d_) * sin(delta) * cos(delta);
    Du3B = q.^2 .* X11 + (1 - alpha)/alpha * I4 * sin(delta) * cos(delta);
    
    %uC
    Du1C = (1 - alpha)* cos(delta) ./ R - q .* Y11 * sin(delta) - alpha * c_  .* q ./ R.^3; 
    Du2C = (1 - alpha)* y_ .* X11 - alpha * c_ .* eta .* q .* X32;
    Du3C = - d_ .* X11 - xi .* Y11 * sin(delta) - alpha * c_ .* ( X11 - q.^2 .* X32);
    %%%%%%%%%%%% displacement gradient
    %%% dip-slip
    %jA
    Dj1A = -alpha/2 * xi .* q ./ R.^3;
    Dj2A = - q/2 .* Y11 - alpha/2 * eta .* q ./ R.^3;
    Dj3A = (1 - alpha)/2 * 1./ R + alpha/2 .* q.^2 ./ R.^3;
    
    %jB
    Dj1B = xi .* q ./ R.^3 + (1 - alpha) / alpha * J4 * sin(delta) * cos(delta);
    Dj2B = eta .* q ./ R.^3 + q.* Y11 + ( 1 -alpha) / alpha * J5 * sin(delta) * cos(delta);
    Dj3B = - q.^2 ./ R.^3 + (1 - alpha) / alpha * J6 * sin(delta)* cos(delta);
    
    %jC
    Dj1C = - (1 - alpha) .* xi./R.^3 * cos(delta) + xi .* q .* Y32 * sin(delta) + alpha.* (3*c_.*xi.*q./R.^5);
    Dj2C = - (1 - alpha) .* y_ ./ R.^3 + alpha * 3 * c_ .* eta .* q ./ R.^5;
    Dj3C =  d_ ./ R.^3 - Y0 .* sin(delta) + alpha * c_ ./ R.^3 .* ( 1 - 3 * q.^2 ./ R.^2);

    %kA
    Dk1A = alpha / 2 * E;
    Dk2A = (1 -alpha) / 2 * d_ .* X11 + xi/ 2 .* Y11 * sin(delta) + alpha/2 * eta .* G; 
    Dk3A = (1 -alpha) / 2 * y_ .* X11 - alpha/2 .* q .* G;
    
    %kB
    Dk1B = -E + (1 - alpha)./ alpha .* J1 * sin(delta) * cos(delta);
    Dk2B = -eta .* G - xi.* Y11 * sin(delta) + (1 - alpha) /alpha *J2 * sin(delta) * cos(delta);
    Dk3B = q .* G + (1 - alpha)./ alpha .* J3 * sin(delta) * cos(delta);
    
    %kC
    Dk1C = -(1 -alpha).* eta ./ R.^3 + Y0 * sin(delta).^2 - alpha * ( (c_ + d_) ./ R.^3 * sin(delta) - 3 * c_ .*  y_ .* q ./ R.^5);
    Dk2C = (1 - alpha) * (X11 - y_.^2 .* X32) - alpha* c_.* ( (d_ + 2* q *cos(delta) ).* X32 - y_ .* eta.* q.* X53);
    Dk3C = xi .* P * sin(delta) + y_ .* d_.* X32 + alpha * c_ .* ( ( y_ + 2 * q * sin(delta) ).* X32 - y_ .* q.^2.* X53);
    
    %lA
    Dl1A = alpha/2 * E_;
    Dl2A = (1 - alpha)/2 .* y_ .* X11 + xi/2 .* Y11 * cos(delta) + alpha/2 * eta.* G_;
    Dl3A = - (1 - alpha)/ 2 * d_ .* X11 - alpha /2 * q .* G_;
    %lB
    Dl1B =  - E_ - (1 - alpha)/ alpha .* K3 .* sin(delta) .* cos(delta);
    Dl2B = - eta .* G_ - xi.* Y11* cos(delta) - (1 -alpha)./alpha .* xi .* D11 * sin(delta) * cos(delta);
    Dl3B = q.*G_ - (1 -alpha)/alpha * K4 .* sin(delta) * cos(delta);
    %lB
    Dl1C = -q./ R.^3 + Y0 * sin(delta) * cos(delta) - alpha * ( ( c_ + d_ )./ R.^3 * cos(delta)  + 3* c_ .* d_ .* q ./ R.^5);
    Dl2C = (1 - alpha) * y_ .* d_ .* X32 - alpha * c_.* ( ( y_ -2* q * sin(delta) ).* X32 + d_ .* eta .* q .* X53 ) ;
    Dl3C = - xi.* P_ * sin(delta) + X11 - d_.^2.* X32 - alpha* c_ .* ( ( d_ - 2*q * cos(delta) ) .* X32 - d_ .* q.^2 .* X53);
end % B2
if tensile ~= 0
    %%%tensile
    %displacement
    %%uA
    Tu1A = -(1 - alpha) / 2 * log(R + eta) - alpha / 2 * q.^2 .* Y11;  
    Tu2A = -(1 - alpha) / 2 * log(R + xi) - alpha / 2 * q.^2 .* X11;
    Tu3A = theta / 2 - alpha / 2 * q .* (eta .* X11 + xi .* Y11); 
    %%
    %uB
    Tu1B = q.^2 .* Y11 - (1 - alpha)/ alpha * I3 * sin(delta)^2;
    Tu2B = q.^2 .* X11 + (1 - alpha)/ alpha * xi ./ (R + d_) * sin(delta)^2;
    Tu3B = q .* (eta .* X11 + xi .* Y11) - theta - (1 - alpha) /alpha * I4 * sin(delta)^2;
    
    %uC
    Tu1C = -(1 - alpha) * ( sin(delta)./ R + q .* Y11 * cos(delta) ) - alpha * ( ones(4,1)*z .* Y11 - q.^2 .* Z32); 
    Tu2C =  (1 - alpha) * 2 * xi .* Y11 * sin(delta) + d_ .* X11 - alpha * c_ .* ( X11 - q.^2 .* X32);
    Tu3C =  (1 - alpha) * (y_ .* X11 + xi .* Y11 * cos(delta)) + alpha * q .* (c_ .* eta .* X32 + xi .* Z32);
    %%%%%%%%%%%% displacement gradient
    %%% tensile
    %jA
    Tj1A = - (1- alpha)/2 * xi .* Y11 + alpha/2 * xi.*q.^2.* Y32;
    Tj2A = -(1-alpha)/2 * 1./R + alpha/2 * q.^2 ./R.^3;
    Tj3A = -(1 -alpha)/2 * q.* Y11 - alpha/2 *q.^3.*Y32;
    
    %jB
    Tj1B = -xi.* q.^2.* Y32 - (1-alpha)/alpha*J4*sin(delta)^2;
    Tj2B = -q.^2./R.^3 -(1-alpha)/alpha*J5*sin(delta)^2;
    Tj3B = q.^3.*Y32 - (1 -alpha)/alpha*J6*sin(delta)^2;
    
    %jC
    Tj1C = (1 -alpha)*xi./R.^3* sin(delta) + xi.*q.*Y32*cos(delta) + alpha * xi.*(3*c_ .*eta./R.^5 - 2 * Z32 - Z0);
    Tj2C = (1 -alpha).* 2 * Y0 * sin(delta) - d_./R.^3 + alpha.* c_./R.^3.*(1 - 3*q.^2./R.^2);
    Tj3C = -(1 -alpha).*( y_./R.^3 - Y0*cos(delta)) - alpha* (3*c_.*eta.*q./R.^5 - q.* Z0);
    %%%tensile
    %kA
    Tk1A = -(1 - alpha)/2 * (cos(delta)./R + q.*Y11*sin(delta)) - alpha/2 * q.* F;
    Tk2A = -(1 - alpha)/2 * y_.*X11 - alpha/2 * q .* G ;
    Tk3A =  (1 - alpha)/2 * (d_.*X11 + xi.* Y11 * sin(delta) ) + alpha /2 * q .* H;
    
    %kB
    Tk1B =  q.*F - (1 -alpha)/alpha * J1 * sin(delta)^2;
    Tk2B =  q.*G - (1 -alpha)/alpha * J2 * sin(delta)^2;
    Tk3B = -q.*H - (1 -alpha)/alpha * J3 * sin(delta)^2;
    
    %kC
    Tk1C = (1 -alpha)* (q./R.^3 + Y0 * sin(delta) * cos(delta)) + alpha *( ones(4,1)*z./ R.^3 * cos(delta) + 3 * c_ .* d_ .* q ./ R.^5 - q.* Z0 * sin(delta) );
    Tk2C = -(1 -alpha) * 2 * xi .* P * sin(delta) - y_.* d_ .* X32 + alpha * c_ .* ( ( y_ +2 *q*sin(delta)).*X32 - y_.* q.^2 .*X53);
    Tk3C = -(1 - alpha) * ( xi.* P * cos(delta) - X11 + y_.^2 .*X32) + alpha* c_ .* ( ( d_ + 2 * q * cos(delta)) .* X32 - y_.* eta .* q .* X53 ) + alpha * xi.* Q;
    
    %%%% tensile
    Tl1A = (1 - alpha)/2 *( sin(delta)./R - q.* Y11 * cos(delta) ) - alpha/2 *q .* F_;
    Tl2A = (1 - alpha)/2 * d_.* X11 - alpha/2 * q .* G_;
    Tl3A = (1 - alpha)/2 * (y_ .* X11 + xi.* Y11 * cos(delta)) + alpha/2 * q .* H_;
    
    Tl1B =  q.* F_ + (1 -alpha)/alpha * K3 * sin(delta)^2;
    Tl2B =  q.* G_ + (1 -alpha)/alpha * xi .* D11 * sin(delta)^2;
    Tl3B = -q.* H_ + (1 -alpha)/alpha * K4 * sin(delta)^2;
    
    Tl1C = - eta./R.^3 + Y0* cos(delta)^2 - alpha *(ones(4,1)*z./R.^3 * sin(delta) - 3*c_.* y_.* q./R.^5 - Y0* sin(delta)^2 + q.*Z0 * cos(delta) );
    Tl2C =  (1 - alpha)* 2* xi.* P_ *sin(delta) - X11 + d_.^2.*X32 -alpha * c_ .*( (d_ -2*q*cos(delta)).*X32 -d_.*q.^2.*X53);
    Tl3C = (1 -alpha)* (xi.* P_ * cos(delta) + y_.* d_.* X32) + alpha.* c_ .* ((y_ - 2 * q * sin(delta) ).* X32 + d_.*eta.*q.*X53) + alpha*xi.* Q_;
end % B3
    %%%%%%%%%%%%
    %%% for u1A_, 
    Z_ = -z;
    d = c - Z_;
    p = y * cos(delta) + d * sin(delta);
   
            xi = [x;x;x-L;x-L]; % 
            eta = [p;p-W;p;p-W];%
    
    q = ones(4,1)*y * sin(delta) - ones(4,1)*d * cos(delta);
   %alpha = (lambda + mu) / (lambda + 2* mu);
   %alpha = 2/3;
   Rsquare = xi.^2 + eta.^2 + q.^2;
   R = sqrt(Rsquare);
   y_ = eta * cos(delta) + q * sin(delta);
   d_ = eta * sin(delta) - q * cos(delta);
   %%
   X11 = 1./(R .* (R + xi));
   X32 = (2 * R + xi )./ ( R.^3 .* ( R + xi ).^2 );
   %X53 = (8*R.^2 + 9*R.*xi + 3*xi.^2) ./ (R.^5 *( R + xi).^3 );
   Y11 = 1 ./ ( R.* ( R + eta ) );
   Y32 = (2 * R + eta ) ./ (R.^3 .* (R + eta).^2 );
   %%%
   %%% Selecting a right root for theta is the most tricky problem
   qsign = sign(q);
   theta = atan2(xi .* eta, abs(q) .* R);
   theta = qsign.* theta;

    X = sqrt( xi.^2 + q.^2);
    if abs(cos(delta))< 0.000001
        I3 = 1/2 * ( eta ./ (R + d_) + y_.*q ./( (R + d_).^2 ) - log(R + eta) );
        I4 = 1/2 * ( xi.* y_ ./( (R + d_).^2) );
    else
    I3 = 1/cos(delta)* y_ ./ ( R + d_ ) - 1 /cos(delta)^2 *( log(R + eta) - sin(delta)* log(R + d_));
    I4 = sin(delta)/cos(delta) * xi./(R + d_) + 2/(cos(delta)^2) * atan2(eta.*(X + q*cos(delta) )...
        + X .* (R + X) * sin(delta), xi.* (R + X)* cos(delta));
    end
    
    %%% ki
    E = sin(delta) ./ R  - y_ .* q ./ R.^3;
    F = d_ ./ R.^3 + xi.^2 .* Y32 * sin(delta);
    G = 2 * X11 * sin(delta) - y_.* q.* X32;
    H = d_ .* q .* X32 + xi .* q .* Y32 * sin(delta);
    
    % li
    E_ = cos(delta)./R + d_ .* q ./ R.^3;
    F_ = y_ ./ R.^3 + xi.^2 .* Y32 * cos(delta);
    G_ = 2 * X11.* cos(delta) + d_.* q.* X32;
    H_ = y_.* q.* X32 + xi.* q .* Y32 * cos(delta);

    if slip_str ~= 0
    %%%%%% strike-slip
    %displacement
     % u1A_
    Su1A_ = theta/2 + alpha / 2 * xi .* q .* Y11;
    Su2A_ = alpha / 2 * q ./ R;
    Su3A_ = (1 - alpha) / 2 * log( R + eta) - alpha / 2 * q.^2 .* Y11;
    %%%%%%%%% for displacement gradient
    % jA_
    Sj1A_ = -(1 -alpha)/2 *q.*Y11 - alpha/2*xi.^2.*q.*Y32;
    Sj2A_ = - alpha/2 * xi.*q./R.^3;
    Sj3A_ = (1 -alpha)/2 .* xi.*Y11 + alpha/2 * xi.*q.^2.*Y32;
    % kA
    Sk1A_ =(1 -alpha)/2 * xi.* Y11* sin(delta) + d_ / 2 .* X11 + alpha/2 * xi .* F; 
    Sk2A_ = alpha /2 * E; 
    Sk3A_ = (1 -alpha)/2 *( cos(delta)./R + q.* Y11* sin(delta) ) - alpha/2 * q .* F;
    %lA
    Sl1A_ = (1 -alpha)/2 * xi .* Y11 * cos(delta) + y_/2 .* X11 + alpha/2 * xi.* F_;
    Sl2A_ =  alpha/2 * E_;
    Sl3A_ = - (1 - alpha)/2 * (sin(delta)./R - q.* Y11 * cos(delta)) - alpha/2 * q .* F_;
    end % B4
if slip_dip ~= 0
    %%%%% dip-slip
    %%% displacement
    % u1A_
    Du1A_ = alpha / 2 * q ./ R; 
    Du2A_ = theta/ 2 + alpha /2 * eta .* q .* X11;
    Du3A_ = (1-alpha) /2 * log(R + xi) - alpha/2 * q.^2 .* X11;
    %%% displacement gradient
    % jA_
    Dj1A_ = -alpha/2 * xi .* q ./ R.^3;
    Dj2A_ = - q/2 .* Y11 - alpha/2 * eta .* q ./ R.^3;
    Dj3A_ = (1 - alpha)/2 * 1./ R + alpha/2 .* q.^2 ./ R.^3;
    %kA
    Dk1A_ = alpha / 2 .* E;
    Dk2A_ = (1 -alpha) / 2 * d_ .* X11 + xi/ 2 .* Y11 * sin(delta) + alpha/2 * eta .* G; 
    Dk3A_ = (1 -alpha) / 2 * y_ .* X11 - alpha/2 .* q .* G;
    %lA
    Dl1A_ = alpha/2 * E_;
    Dl2A_ = (1 - alpha)/2 .* y_ .* X11 + xi/2 .* Y11 * cos(delta) + alpha/2 * eta.* G_;
    Dl3A_ = - (1 - alpha)/ 2 * d_ .* X11 - alpha /2 * q .* G_;
end % B5
if tensile ~= 0 % B6
    %%% tensile for u1A_
    %%%%%%%%%%%%%% displacement
    Tu1A_ = -(1 - alpha) / 2 * log(R + eta) - alpha / 2 * q.^2 .* Y11;  
    Tu2A_ = -(1 - alpha) / 2 * log(R + xi) - alpha / 2 * q.^2 .* X11;
    Tu3A_ = theta / 2 - alpha / 2 * q .* (eta .* X11 + xi .* Y11); 
    %%%%%%%%% displacement gradient
    %%%tensile 
    % jA_
    Tj1A_ = - (1- alpha)/2 * xi .* Y11 + alpha/2 * xi.*q.^2.* Y32;
    Tj2A_ =-(1-alpha)/2 * 1./R + alpha/2 * q.^2 ./R.^3;
    Tj3A_ = -(1 -alpha)/2 * q.* Y11 - alpha/2 *q.^3.*Y32;
    %kA
    Tk1A_ =-(1 - alpha)/2 * (cos(delta)./R + q.*Y11*sin(delta)) - alpha/2 * q.* F;
    Tk2A_ =-(1 - alpha)/2 * y_.*X11 - alpha/2 * q .* G ;
    Tk3A_ =  (1 - alpha)/2 * (d_.*X11 + xi.* Y11 * sin(delta) ) + alpha /2 * q .* H;
    %lA
    Tl1A_ =(1 - alpha)/2 *( sin(delta)./R - q.* Y11 * cos(delta) ) - alpha/2 *q .* F_;
    Tl2A_ = (1 - alpha)/2 * d_.* X11 - alpha/2 * q .* G_;
    Tl3A_ =(1 - alpha)/2 * (y_ .* X11 + xi.* Y11 * cos(delta)) + alpha/2 * q .* H_;
end % B6    
if slip_str ~= 0 % B7
    %%%%%% strike-slip
    %% displacement
    Sux = 1 / (2 * pi) * slip_str * (Su1A - Su1A_ + Su1B + ones(4,1)*z.* Su1C  );
    Suy = 1 / (2 * pi) * slip_str * ( (Su2A - Su2A_ + Su2B + ones(4,1)*z.* Su2C) * cos(delta) - (Su3A - Su3A_ + Su3B + ones(4,1)*z .* Su3C) * sin(delta) );
    Suz = 1 / (2 * pi) * slip_str * ( (Su2A - Su2A_ + Su2B - ones(4,1)*z.* Su2C) * sin(delta) + (Su3A - Su3A_ + Su3B - ones(4,1)*z .* Su3C) * cos(delta) );

    % displacement gradients
    Sduxdx = 1 / (2 * pi) * slip_str * (Sj1A - Sj1A_ + Sj1B + ones(4,1)*z.* Sj1C  );
    Sduydx = 1 / (2 * pi) * slip_str * ( (Sj2A - Sj2A_ + Sj2B + ones(4,1)*z.* Sj2C) * cos(delta) - (Sj3A - Sj3A_ + Sj3B + ones(4,1)*z .* Sj3C) * sin(delta) );
    Sduzdx = 1 / (2 * pi) * slip_str * ( (Sj2A - Sj2A_ + Sj2B - ones(4,1)*z.* Sj2C) * sin(delta) + (Sj3A - Sj3A_ + Sj3B - ones(4,1)*z .* Sj3C) * cos(delta) );

    Sduxdy = 1 / (2 * pi) * slip_str * (Sk1A - Sk1A_ + Sk1B + ones(4,1)*z.* Sk1C  );
    Sduydy = 1 / (2 * pi) * slip_str * ( (Sk2A - Sk2A_ + Sk2B + ones(4,1)*z.* Sk2C) * cos(delta) - (Sk3A - Sk3A_ + Sk3B + ones(4,1)*z .* Sk3C) * sin(delta) );
    Sduzdy = 1 / (2 * pi) * slip_str * ( (Sk2A - Sk2A_ + Sk2B - ones(4,1)*z.* Sk2C) * sin(delta) + (Sk3A - Sk3A_ + Sk3B - ones(4,1)*z .* Sk3C) * cos(delta) );

    Sduxdz = 1 / (2 * pi) * slip_str * (  Sl1A + Sl1A_ + Sl1B + Su1C + ones(4,1)*z.* Sl1C  );
    Sduydz = 1 / (2 * pi) * slip_str * ( (Sl2A + Sl2A_ + Sl2B + Su2C + ones(4,1)*z.* Sl2C) * cos(delta) - (Sl3A + Sl3A_ + Sl3B + Su3C + ones(4,1)*z .* Sl3C) * sin(delta) );
    Sduzdz = 1 / (2 * pi) * slip_str * ( (Sl2A + Sl2A_ + Sl2B - Su2C - ones(4,1)*z.* Sl2C) * sin(delta) + (Sl3A + Sl3A_ + Sl3B - Su3C - ones(4,1)*z .* Sl3C) * cos(delta) );
else % B7
    Sux = 0; Suy = 0; Suz = 0;
    Sduxdx = 0; Sduydx = 0; Sduzdx = 0; 
    Sduxdy = 0; Sduydy = 0; Sduzdy = 0;
    Sduxdz = 0; Sduydz = 0; Sduzdz = 0;
end
    %%%%% dip-slip
    %%% displacement
if slip_dip ~= 0 % B8
    Dux = 1 / (2 * pi) * slip_dip * (Du1A - Du1A_ + Du1B + ones(4,1)*z.* Du1C  );
    Duy = 1 / (2 * pi) * slip_dip * ( (Du2A - Du2A_ + Du2B + ones(4,1)*z.* Du2C) * cos(delta) - (Du3A - Du3A_ + Du3B + ones(4,1)*z .* Du3C) * sin(delta) );
    Duz = 1 / (2 * pi) * slip_dip * ( (Du2A - Du2A_ + Du2B - ones(4,1)*z.* Du2C) * sin(delta) + (Du3A - Du3A_ + Du3B - ones(4,1)*z .* Du3C) * cos(delta) );

    % displacement gradients
    Dduxdx = 1 / (2 * pi) * slip_dip * (Dj1A - Dj1A_ + Dj1B + ones(4,1)*z.* Dj1C  );
    Dduydx = 1 / (2 * pi) * slip_dip * ( (Dj2A - Dj2A_ + Dj2B + ones(4,1)*z.* Dj2C) * cos(delta) - (Dj3A - Dj3A_ + Dj3B + ones(4,1)*z .* Dj3C) * sin(delta) );
    Dduzdx = 1 / (2 * pi) * slip_dip * ( (Dj2A - Dj2A_ + Dj2B - ones(4,1)*z.* Dj2C) * sin(delta) + (Dj3A - Dj3A_ + Dj3B - ones(4,1)*z .* Dj3C) * cos(delta) );

    Dduxdy = 1 / (2 * pi) * slip_dip * (Dk1A - Dk1A_ + Dk1B + ones(4,1)*z.* Dk1C  );
    Dduydy = 1 / (2 * pi) * slip_dip * ( (Dk2A - Dk2A_ + Dk2B + ones(4,1)*z.* Dk2C) * cos(delta) - (Dk3A - Dk3A_ + Dk3B + ones(4,1)*z .* Dk3C) * sin(delta) );
    Dduzdy = 1 / (2 * pi) * slip_dip * ( (Dk2A - Dk2A_ + Dk2B - ones(4,1)*z.* Dk2C) * sin(delta) + (Dk3A - Dk3A_ + Dk3B - ones(4,1)*z .* Dk3C) * cos(delta) );

    Dduxdz = 1 / (2 * pi) * slip_dip * (  Dl1A + Dl1A_ + Dl1B + Du1C + ones(4,1)*z.* Dl1C  );
    Dduydz = 1 / (2 * pi) * slip_dip * ( (Dl2A + Dl2A_ + Dl2B + Du2C + ones(4,1)*z.* Dl2C) * cos(delta) - (Dl3A + Dl3A_ + Dl3B + Du3C + ones(4,1)*z .* Dl3C) * sin(delta) );
    Dduzdz = 1 / (2 * pi) * slip_dip * ( (Dl2A + Dl2A_ + Dl2B - Du2C - ones(4,1)*z.* Dl2C) * sin(delta) + (Dl3A + Dl3A_ + Dl3B - Du3C - ones(4,1)*z .* Dl3C) * cos(delta) );
else % B8
    Dux = 0; Duy = 0; Duz = 0;
    Dduxdx = 0; Dduydx = 0; Dduzdx = 0;
    Dduxdy = 0; Dduydy = 0; Dduzdy = 0;
    Dduxdz = 0; Dduydz = 0; Dduzdz = 0;
end  % B8
    if tensile ~= 0
    % tensile
    %%%%%%%%%%% displacement
    Tux = 1 / (2 * pi) * tensile * (Tu1A - Tu1A_ + Tu1B + ones(4,1)*z.* Tu1C  );
    Tuy = 1 / (2 * pi) * tensile * ( (Tu2A - Tu2A_ + Tu2B + ones(4,1)*z.* Tu2C) * cos(delta) - (Tu3A - Tu3A_ + Tu3B + ones(4,1)*z .* Tu3C) * sin(delta) );
    Tuz = 1 / (2 * pi) * tensile * ( (Tu2A - Tu2A_ + Tu2B - ones(4,1)*z.* Tu2C) * sin(delta) + (Tu3A - Tu3A_ + Tu3B - ones(4,1)*z .* Tu3C) * cos(delta) );

    %%%%%%%%%% displacement gradients
    Tduxdx = 1 / (2 * pi) * tensile * (Tj1A - Tj1A_ + Tj1B + ones(4,1)*z.* Tj1C  );
    Tduydx = 1 / (2 * pi) * tensile * ( (Tj2A - Tj2A_ + Tj2B + ones(4,1)*z.* Tj2C) * cos(delta) - (Tj3A - Tj3A_ + Tj3B + ones(4,1)*z .* Tj3C) * sin(delta) );
    Tduzdx = 1 / (2 * pi) * tensile * ( (Tj2A - Tj2A_ + Tj2B - ones(4,1)*z.* Tj2C) * sin(delta) + (Tj3A - Tj3A_ + Tj3B - ones(4,1)*z .* Tj3C) * cos(delta) );

    Tduxdy = 1 / (2 * pi) * tensile * (Tk1A - Tk1A_ + Tk1B + ones(4,1)*z.* Tk1C  );
    Tduydy = 1 / (2 * pi) * tensile * ( (Tk2A - Tk2A_ + Tk2B + ones(4,1)*z.* Tk2C) * cos(delta) - (Tk3A - Tk3A_ + Tk3B + ones(4,1)*z .* Tk3C) * sin(delta) );
    Tduzdy = 1 / (2 * pi) * tensile * ( (Tk2A - Tk2A_ + Tk2B - ones(4,1)*z.* Tk2C) * sin(delta) + (Tk3A - Tk3A_ + Tk3B - ones(4,1)*z .* Tk3C) * cos(delta) );

    Tduxdz = 1 / (2 * pi) * tensile * (  Tl1A + Tl1A_ + Tl1B + Tu1C + ones(4,1)*z.* Tl1C  );
    Tduydz = 1 / (2 * pi) * tensile * ( (Tl2A + Tl2A_ + Tl2B + Tu2C + ones(4,1)*z.* Tl2C) * cos(delta) - (Tl3A + Tl3A_ + Tl3B + Tu3C + ones(4,1)*z .* Tl3C) * sin(delta) );
    Tduzdz = 1 / (2 * pi) * tensile * ( (Tl2A + Tl2A_ + Tl2B - Tu2C - ones(4,1)*z.* Tl2C) * sin(delta) + (Tl3A + Tl3A_ + Tl3B - Tu3C - ones(4,1)*z .* Tl3C) * cos(delta) );
    else % B8
    Tux = 0;  Tuy = 0;  Tuz = 0;
    Tduxdx = 0; Tduydx = 0; Tduzdx = 0;
    Tduxdy = 0; Tduydy = 0; Tduzdy = 0;
    Tduxdz = 0; Tduydz = 0; Tduzdz = 0;
    end % B8
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    factor =[ones(1,size(xi,2)); -ones(1,size(xi,2)); -ones(1,size(xi,2)); ones(1,size(xi,2))];
    G1 =  sum(factor.* (Sux + Dux + Tux));
    G2 =  sum(factor.* (Suy + Duy + Tuy));
    G3 =  sum(factor.* (Suz + Duz + Tuz));

    Dg11 = sum(factor.*(Sduxdx + Dduxdx + Tduxdx));
    Dg12 = sum(factor.*(Sduxdy + Dduxdy + Tduxdy));
    Dg13 = sum(factor.*(Sduxdz + Dduxdz + Tduxdz));

    Dg21 = sum(factor.*(Sduydx + Dduydx + Tduydx));
    Dg22 = sum(factor.*(Sduydy + Dduydy + Tduydy));
    Dg23 = sum(factor.*(Sduydz + Dduydz + Tduydz));

    Dg31 = sum(factor.*(Sduzdx + Dduzdx + Tduzdx));
    Dg32 = sum(factor.*(Sduzdy + Dduzdy + Tduzdy));
    Dg33 = sum(factor.*(Sduzdz + Dduzdz + Tduzdz));

%%%%%%%%%%%%%%%%%%%%
% coordinate transformation
%Ux_y_ =[cos(angle_Str), -sin(angle_Str); sin(angle_Str),cos(angle_Str)] * U_;
Gx = cos(angle_Str)*(-G2) - sin(angle_Str)*G1;
Gy = sin(angle_Str)*(-G2) + cos(angle_Str)*G1;
Gz = G3;

displacement=[Gx;Gy;Gz];
%tempGrad1 = [Dg11,Dg12,Dg13;Dg21,Dg22,Dg23;Dg31,Dg32,Dg33]
Dg11_ =   Dg22;
Dg12_ =  -Dg21;
Dg13_ =  -Dg23;

Dg21_ =  -Dg12;
Dg22_ =   Dg11;
Dg23_ =   Dg13;

Dg31_ =  -Dg32;
Dg32_ =   Dg31;
Dg33_ =   Dg33;

%tempGrad2 = [Dg11_, Dg12_, Dg13_;Dg21_, Dg22_, Dg23_;Dg31_,Dg32_,Dg33_]
    
Dgxx = (cos(angle_Str)*Dg11_ - sin(angle_Str)*Dg21_) * cos(angle_Str)...
    + (cos(angle_Str)* Dg12_ - sin(angle_Str)*Dg22_ )* (-sin(angle_Str));
Dgyx = (sin(angle_Str)*Dg11_ + cos(angle_Str)*Dg21_ ) * ( cos(angle_Str))...
    - (sin(angle_Str)*Dg12_ + cos(angle_Str)*Dg22_ ) * sin(angle_Str);
Dgzx = Dg31_ * cos(angle_Str) - Dg32_ * sin(angle_Str);

Dgxy = (cos(angle_Str)*Dg11_ - sin(angle_Str)*Dg21_) * sin(angle_Str)...
    + (cos(angle_Str) * Dg12_ - sin(angle_Str) * Dg22_) * (cos(angle_Str));
Dgyy = ( sin(angle_Str)*Dg11_ + cos(angle_Str)*Dg21_) * (sin(angle_Str))...
    + (sin(angle_Str)*Dg12_ + cos(angle_Str)*Dg22_) * cos(angle_Str);
Dgzy =  sin(angle_Str) *Dg31_ + cos(angle_Str) * Dg32_;

Dgxz =  cos(angle_Str)*Dg13_ - sin(angle_Str)* Dg23_;
Dgyz = sin(angle_Str)* Dg13_ + cos(angle_Str)* Dg23_ ;
Dgzz =  Dg33_;

gradient =[ Dgxx; Dgxy;Dgxz; Dgyx;Dgyy;Dgyz; Dgzx; Dgzy;Dgzz];

Ex = Dgxx;
Ey = Dgyy;
Ez = Dgzz;
Exy = 1/2*(Dgyx + Dgxy);
Eyz = 1/2*(Dgyz + Dgzy);
Ezx = 1/2*(Dgzx + Dgxz);

Sx = youngs/((1 + nu)*(1 -2*nu))*(Ex + nu *(Ey + Ez - Ex));
Sy = youngs/((1 + nu)*(1 -2*nu))*(Ey + nu *(Ex + Ez - Ey));
Sz = youngs/((1 + nu)*(1 -2*nu))*(Ez + nu *(Ey + Ex - Ez));
Sxy = 2*Gshear*Exy;
Syz = 2*Gshear*Eyz;
Szx = 2*Gshear*Ezx;

Stress =[ Sx; Sxy; Szx; Sy; Syz; Sz];

else % A1
    if  1 == isequal(Z,-abs(Z))% A1
       fprintf('warning: physically impossible')     
    else
       fprintf('warning: All z should be negtive.')    
    end
        displacement = nan;
        gradient = nan;
        Stress = nan;
end  % A1