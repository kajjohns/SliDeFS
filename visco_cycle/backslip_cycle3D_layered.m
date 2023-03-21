function U=backslip_cycle3D(m,xloc,H1,H2,mu,lam,t,tR1,tR2,T,Nterms)


%%INPUTS:
%fault model (standard Okada paramterization -- must be column vector):
% m(1) = length (km)
% m(2) = width (km)
% m(3) = depth to down-dip edge (km)
% m(4) = dip (degrees)
% m(5) = strike (degrees)
% m(6) = east position of down-dip edge (km)
% m(7) = north position (km)
% m(8) = strike-slip (m)
% m(9) = dip-slip (m)
% m(10) = tensile (m)
% d = 1x(M-1) vector of depths to bottom of layer -- M is number of layers
%      note: bottom layer is halfspace, thus no depth is specified
% mu,lam = 1xM vector of Lame constants -- mu=lam for poisson ratio = 0.25
%        -- M is number of layers
% xloc = 2xn matrix of n station coordinates (km) -- first row x, second row y
%%OUTPUTS:
%U is a 3xn matrix of VELOCITIES (first row east, second row north, third row up) 

if nargin<11
    Nterms=[];
end

xy=xloc'; 
%initialize basis vectors to zero
G1=zeros(3*size(xy,1),1);
G2=zeros(3*size(xy,1),1);
G3=zeros(3*size(xy,1),1);

%define components of slip
ss=m(8);
ds=m(9);
ten=m(10);

%convert parameters to be consistent with Okada
if m(4)<=90 & m(4)>=-90
   dipdir = (m(5) + 90)*pi/180;
else
   dipdir = (m(5) - 90)*pi/180;
end

offset = abs(m(2).*cos(m(4)*pi/180));
EastOff = offset.*sin(dipdir);
NorthOff = offset.*cos(dipdir);
xy(:,1)=xy(:,1)-m(6); %+EastOff;
xy(:,2)=xy(:,2)-m(7); %+NorthOff;

L=m(1);
W=m(2);
dip=m(4);
strike=m(5);
D=m(3)-W*sin(dip*pi/180);


Xcoord=xy(:,2);	%switch x and y
Ycoord=xy(:,1);

bulk=lam+2*mu/3;   %bulk modulus
g=lam+2*mu;
beta1 = 1/tR1;
beta2 = 1/tR2;


%position of point sources
%setup spacing of point sources for Gaussian quadrature
%number of point sources near ground surface
   

%NL=ceil(.0075*m(1));
NL=ceil(.025*m(1));
NW=ceil(.05*m(2));



u=(1:NL-1)./sqrt((2*(1:NL-1)).^2-1);
	[vc,bp]=eig(diag(u,-1)+diag(u,1));
	[bp,k]=sort(diag(bp));
	a=-L/2;b=L/2;
	wL=(2*vc(1,k).^2)*(b-a)/2;
	xp=(a+b)/2+(b-a)/2*bp;

	u=(1:NW-1)./sqrt((2*(1:NW-1)).^2-1);
	[vc,bp]=eig(diag(u,-1)+diag(u,1));
	[bp,k]=sort(diag(bp));
   
   a=0;b=W;


   wW=(2*vc(1,k).^2)*(b-a)/2;
	yp=(a+b)/2+(b-a)/2*bp;
   
[xp,yp]=meshgrid(xp,yp);   
  
xp=xp(:)';
yp=yp(:)';
zp=zeros(size(xp));


wW=repmat(wW,size(xy,1),1);
wW=repmat(wW,[1 1 NL]);
wL=reshape(wL,[1 1 NL]);
wL=repmat(wL,size(xy,1),NW);

%mesh points
%rotate into true position
%rotate plane to true dip
R=[1 0 0;0 cos(dip*pi/180) -sin(dip*pi/180);0 sin(dip*pi/180) cos(dip*pi/180)];	%rotation matrix x-axis
Xp=R*[xp;yp;zp];
%rotate to true strike
R=[cos(strike*pi/180) -sin(strike*pi/180) 0;sin(strike*pi/180) cos(strike*pi/180) 0;0 0 1];	%rotation matrix about z-axis
Xp=R*Xp;
xpos=Xp(1,:);
ypos=Xp(2,:);
zpos=Xp(3,:);
xpos=reshape(xpos,NW,NL);
ypos=reshape(ypos,NW,NL);
zpos=zpos(1:NW);


for psW=1:NW

      
zs=D+zpos(psW);		%depth of point source

ZS=zs;
if zs<0
   disp('Warning: fault extends above ground surface')
end

if isempty(Nterms)
    NN=100;
    N=100;
else
    NN=Nterms;
    N=Nterms;
end

%kmax=.0625;
kmax=.25;
k=linspace(0.000001,kmax,NN);
K=k;

kstore{psW}=K;
Nstore{psW}=NN;

[M1,M2,M3]=momtensor_inverse(strike,dip,lam,mu);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate Bessel coefficients

for j=1:N
   
[sourceP4x4,sourceP2x2,halfspaceP4x4,halfspaceP2x2]=getprop(H1,k,j,mu,lam,g,zs);
 
Ph2=halfspaceP2x2;
Ph4=halfspaceP4x4;

%% make source vectors and get bessel coefficients

if ss~=0
   [UR_0_1(:,j),US_0_1(:,j), ...
         UR_n1_1(:,j),US_n1_1(:,j),UT_n1_1(:,j), ...
         UR_p1_1(:,j),US_p1_1(:,j),UT_p1_1(:,j), ...
         UR_n2_1(:,j),US_n2_1(:,j),UT_n2_1(:,j), ...
         UR_p2_1(:,j),US_p2_1(:,j),UT_p2_1(:,j)]=getbesselco(Ph2,Ph4,sourceP4x4,...
         sourceP2x2,M1,g,k,bulk,j,lam,mu,beta1,beta2,H1,H2,t,T);

end

if ds~=0
[UR_0_2(:,j),US_0_2(:,j), ...
         UR_n1_2(:,j),US_n1_2(:,j),UT_n1_2(:,j), ...
         UR_p1_2(:,j),US_p1_2(:,j),UT_p1_2(:,j), ...
         UR_n2_2(:,j),US_n2_2(:,j),UT_n2_2(:,j), ...
         UR_p2_2(:,j),US_p2_2(:,j),UT_p2_2(:,j)]=getbesselco(Ph2,Ph4,sourceP4x4,...
         sourceP2x2,M2,g,k,bulk,j,lam,mu,beta1,beta2,H1,H2,t,T);
end

if ten~=0
[UR_0_3(:,j),US_0_3(:,j), ...
         UR_n1_3(:,j),US_n1_3(:,j),UT_n1_3(:,j), ...
         UR_p1_3(:,j),US_p1_3(:,j),UT_p1_3(:,j), ...
         UR_n2_3(:,j),US_n2_3(:,j),UT_n2_3(:,j), ...
         UR_p2_3(:,j),US_p2_3(:,j),UT_p2_3(:,j)]=getbesselco(Ph2,Ph4,sourceP4x4,...
         sourceP2x2,M3,g,k,bulk,j,lam,mu,beta1,beta2,H1,H2,t,T);
end

end %j

if ss~=0
	   UR_0_1s{psW}=UR_0_1;
      US_0_1s{psW}=US_0_1;
      
  	   UR_n1_1s{psW}=UR_n1_1;
      US_n1_1s{psW}=US_n1_1;
      UT_n1_1s{psW}=UT_n1_1;
      
  	   UR_p1_1s{psW}=UR_p1_1;
      US_p1_1s{psW}=US_p1_1;
      UT_p1_1s{psW}=UT_p1_1;
      
  	   UR_n2_1s{psW}=UR_n2_1;
      US_n2_1s{psW}=US_n2_1;
      UT_n2_1s{psW}=UT_n2_1;
      
  	   UR_p2_1s{psW}=UR_p2_1;
      US_p2_1s{psW}=US_p2_1;
      UT_p2_1s{psW}=UT_p2_1;
end
   
if ds~=0
	   UR_0_2s{psW}=UR_0_2;
      US_0_2s{psW}=US_0_2;
      
  	   UR_n1_2s{psW}=UR_n1_2;
      US_n1_2s{psW}=US_n1_2;
      UT_n1_2s{psW}=UT_n1_2;
      
  	   UR_p1_2s{psW}=UR_p1_2;
      US_p1_2s{psW}=US_p1_2;
      UT_p1_2s{psW}=UT_p1_2;
      
  	   UR_n2_2s{psW}=UR_n2_2;
      US_n2_2s{psW}=US_n2_2;
      UT_n2_2s{psW}=UT_n2_2;
      
  	   UR_p2_2s{psW}=UR_p2_2;
      US_p2_2s{psW}=US_p2_2;
      UT_p2_2s{psW}=UT_p2_2;
end

if ten~=0
	   UR_0_3s{psW}=UR_0_3;
      US_0_3s{psW}=US_0_3;
      
  	   UR_n1_3s{psW}=UR_n1_3;
      US_n1_3s{psW}=US_n1_3;
      UT_n1_3s{psW}=UT_n1_3;
      
  	   UR_p1_3s{psW}=UR_p1_3;
      US_p1_3s{psW}=US_p1_3;
      UT_p1_3s{psW}=UT_p1_3;
      
  	   UR_n2_3s{psW}=UR_n2_3;
      US_n2_3s{psW}=US_n2_3;
      UT_n2_3s{psW}=UT_n2_3;
      
  	   UR_p2_3s{psW}=UR_p2_3;
      US_p2_3s{psW}=US_p2_3;
      UT_p2_3s{psW}=UT_p2_3;
end


end %psW


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inverse HT into physical space

for psW=1:NW

   for psL=1:NL


xshift=repmat(xpos(psW,psL),length(Xcoord),1);
yshift=repmat(ypos(psW,psL),length(Xcoord),1);

X=Xcoord-xshift;
Y=Ycoord-yshift;

r=sqrt(X.^2+Y.^2);
theta=atan2(Y,X);


k=kstore{psW};
N=Nstore{psW};     




rr=r;
kk=k;
r=repmat(r,1,length(k));
k=repmat(k,length(rr),1);
kr=k.*r;
theta_orig=theta;
theta=repmat(theta,1,N);

if ss~=0
UR_0_1=repmat(UR_0_1s{psW},length(rr),1);
US_0_1=repmat(US_0_1s{psW},length(rr),1);
UR_n1_1=repmat(UR_n1_1s{psW},length(rr),1);
US_n1_1=repmat(US_n1_1s{psW},length(rr),1);
UT_n1_1=repmat(UT_n1_1s{psW},length(rr),1);

UR_p1_1=repmat(UR_p1_1s{psW},length(rr),1);
US_p1_1=repmat(US_p1_1s{psW},length(rr),1);
UT_p1_1=repmat(UT_p1_1s{psW},length(rr),1);

UR_n2_1=repmat(UR_n2_1s{psW},length(rr),1);
US_n2_1=repmat(US_n2_1s{psW},length(rr),1);
UT_n2_1=repmat(UT_n2_1s{psW},length(rr),1);

UR_p2_1=repmat(UR_p2_1s{psW},length(rr),1);
US_p2_1=repmat(US_p2_1s{psW},length(rr),1);
UT_p2_1=repmat(UT_p2_1s{psW},length(rr),1);
end


if ds~=0
UR_0_2=repmat(UR_0_2s{psW},length(rr),1);
US_0_2=repmat(US_0_2s{psW},length(rr),1);
UR_n1_2=repmat(UR_n1_2s{psW},length(rr),1);
US_n1_2=repmat(US_n1_2s{psW},length(rr),1);
UT_n1_2=repmat(UT_n1_2s{psW},length(rr),1);

UR_p1_2=repmat(UR_p1_2s{psW},length(rr),1);
US_p1_2=repmat(US_p1_2s{psW},length(rr),1);
UT_p1_2=repmat(UT_p1_2s{psW},length(rr),1);

UR_n2_2=repmat(UR_n2_2s{psW},length(rr),1);
US_n2_2=repmat(US_n2_2s{psW},length(rr),1);
UT_n2_2=repmat(UT_n2_2s{psW},length(rr),1);

UR_p2_2=repmat(UR_p2_2s{psW},length(rr),1);
US_p2_2=repmat(US_p2_2s{psW},length(rr),1);
UT_p2_2=repmat(UT_p2_2s{psW},length(rr),1);
end

if ten~=0
UR_0_3=repmat(UR_0_3s{psW},length(rr),1);
US_0_3=repmat(US_0_3s{psW},length(rr),1);
UR_n1_3=repmat(UR_n1_3s{psW},length(rr),1);
US_n1_3=repmat(US_n1_3s{psW},length(rr),1);
UT_n1_3=repmat(UT_n1_3s{psW},length(rr),1);

UR_p1_3=repmat(UR_p1_3s{psW},length(rr),1);
US_p1_3=repmat(US_p1_3s{psW},length(rr),1);
UT_p1_3=repmat(UT_p1_3s{psW},length(rr),1);

UR_n2_3=repmat(UR_n2_3s{psW},length(rr),1);
US_n2_3=repmat(US_n2_3s{psW},length(rr),1);
UT_n2_3=repmat(UT_n2_3s{psW},length(rr),1);

UR_p2_3=repmat(UR_p2_3s{psW},length(rr),1);
US_p2_3=repmat(US_p2_3s{psW},length(rr),1);
UT_p2_3=repmat(UT_p2_3s{psW},length(rr),1);
end


M=0;

b0=besselj(0,kr);
b1=besselj(1,kr);
b2=-b0+2./kr.*b1;
b3=-b1+4./kr.*b2;

DrJm=-k.*b1;

   if ss~=0
	Uz_integrand_1=k.*UR_0_1.*b0;
	Ur_integrand_1=US_0_1.*DrJm;
	end
   
   if ds~=0
   Uz_integrand_2=k.*UR_0_2.*b0;
   Ur_integrand_2=US_0_2.*DrJm;
	end
   
   if ten~=0
   Uz_integrand_3=k.*UR_0_3.*b0;
   Ur_integrand_3=US_0_3.*DrJm;
	end


if ss~=0
Uz_1(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,1)=zeros(length(rr),1);
end

if ds~=0
Uz_2(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,1)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,1)=zeros(length(rr),1);
end

M=-1;
DrJm=.5*k.*(b2-b0);
   
   if ss~=0
	Uz_integrand_1=-k.*UR_n1_1.*b1.*exp(i*M*theta);
	Ur_integrand_1=(US_n1_1.*DrJm - i*M*UT_n1_1.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=-k.*UR_n1_2.*b1.*exp(i*M*theta);
	Ur_integrand_2=(US_n1_2.*DrJm - i*M*UT_n1_2.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ten~=0
 	Uz_integrand_3=-k.*UR_n1_3.*b1.*exp(i*M*theta);
	Ur_integrand_3=(US_n1_3.*DrJm - i*M*UT_n1_3.*(1./r).*b1).*exp(i*M*theta);
	end


if ss~=0
Uz_1(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
end

if ds~=0
Uz_2(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,2)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,2)=zeros(length(rr),1);
end

M=1;
DrJm=.5*k.*(b0-b2);

   if ss~=0
	Uz_integrand_1=k.*UR_p1_1.*b1.*exp(i*M*theta);
	Ur_integrand_1=(US_p1_1.*DrJm + i*M*UT_p1_1.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_p1_2.*b1.*exp(i*M*theta);
	Ur_integrand_2=(US_p1_2.*DrJm + i*M*UT_p1_2.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ten~=0
	Uz_integrand_3=k.*UR_p1_3.*b1.*exp(i*M*theta);
	Ur_integrand_3=(US_p1_3.*DrJm + i*M*UT_p1_3.*(1./r).*b1).*exp(i*M*theta);
	end


if ss~=0
Uz_1(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,3)=zeros(length(rr),1);
end

if ds~=0
Uz_2(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,3)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,3)=zeros(length(rr),1);
end


M=-2;
DrJm=.5*k.*(b1-b3);
   
   if ss~=0
	Uz_integrand_1=k.*UR_n2_1.*b2.*exp(i*M*theta);
	Ur_integrand_1=(US_n2_1.*DrJm + i*M*UT_n2_1.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_1=(i*M*US_n2_1.*(1./r).*b2 - UT_n2_1.*DrJm).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_n2_2.*b2.*exp(i*M*theta);
	Ur_integrand_2=(US_n2_2.*DrJm + i*M*UT_n2_2.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_2=(i*M*US_n2_2.*(1./r).*b2 - UT_n2_2.*DrJm).*exp(i*M*theta);
  	end  
  
   if ten~=0
	Uz_integrand_3=k.*UR_n2_3.*b2.*exp(i*M*theta);
	Ur_integrand_3=(US_n2_3.*DrJm + i*M*UT_n2_3.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_3=(i*M*US_n2_3.*(1./r).*b2 - UT_n2_3.*DrJm).*exp(i*M*theta);
	end  

if ss~=0
Uz_1(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,4)=real((1/(2*pi))*trapz(kk,Utheta_integrand_1,2));
end

if ds~=0
Uz_2(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,4)=real((1/(2*pi))*trapz(kk,Utheta_integrand_2,2));
end

if ten~=0
Uz_3(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,4)=real((1/(2*pi))*trapz(kk,Utheta_integrand_3,2));
end


M=2;
DrJm=.5*k.*(b1-b3);
   
   if ss~=0
	Uz_integrand_1=k.*UR_p2_1.*b2.*exp(i*M*theta);
	Ur_integrand_1=(US_p2_1.*DrJm + i*M*UT_p2_1.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_1=(i*M*US_p2_1.*(1./r).*b2 - UT_p2_1.*DrJm).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_p2_2.*b2.*exp(i*M*theta);
	Ur_integrand_2=(US_p2_2.*DrJm + i*M*UT_p2_2.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_2=(i*M*US_p2_2.*(1./r).*b2 - UT_p2_2.*DrJm).*exp(i*M*theta);
  	end  
  
   if ten~=0
	Uz_integrand_3=k.*UR_p2_3.*b2.*exp(i*M*theta);
	Ur_integrand_3=(US_p2_3.*DrJm + i*M*UT_p2_3.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_3=(i*M*US_p2_3.*(1./r).*b2 - UT_p2_3.*DrJm).*exp(i*M*theta);
	end  

if ss~=0
Uz_1(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,5)=real((1/(2*pi))*trapz(kk,Utheta_integrand_1,2));
end

if ds~=0
Uz_2(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,5)=real((1/(2*pi))*trapz(kk,Utheta_integrand_2,2));
end

if ten~=0
Uz_3(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,5)=real((1/(2*pi))*trapz(kk,Utheta_integrand_3,2));
end



if ss~=0
Gz_1(:,psW,psL)=sum(Uz_1,2);
UR=sum(Ur_1,2);
UTHETA=sum(Utheta_1,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_1(:,psW,psL)=Urx+Uthetax;
Gy_1(:,psW,psL)=Ury+Uthetay;
end


if ds~=0
Gz_2(:,psW,psL)=sum(Uz_2,2);
UR=sum(Ur_2,2);
UTHETA=sum(Utheta_2,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_2(:,psW,psL)=Urx+Uthetax;
Gy_2(:,psW,psL)=Ury+Uthetay;
end

if ten~=0
Gz_3(:,psW,psL)=sum(Uz_3,2);
UR=sum(Ur_3,2);
UTHETA=sum(Utheta_3,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_3(:,psW,psL)=Urx+Uthetax;
Gy_3(:,psW,psL)=Ury+Uthetay;
end


end %psW
end %psL


if ss~=0
   
   Gx_1=wW.*wL.*Gx_1;
   Gy_1=wW.*wL.*Gy_1;
   Gz_1=wW.*wL.*Gz_1;
   
   Gx_1=sum(sum(Gx_1,2),3);
	Gy_1=sum(sum(Gy_1,2),3);
	Gz_1=sum(sum(Gz_1,2),3);

	G1=zeros(3*length(Gy_1),1);
	G1(1:3:end)=Gy_1;
	G1(2:3:end)=Gx_1;
	G1(3:3:end)=-Gz_1;
end

if ds~=0
   
   Gx_2=wW.*wL.*Gx_2;
   Gy_2=wW.*wL.*Gy_2;
   Gz_2=wW.*wL.*Gz_2;
   
   Gx_2=sum(sum(Gx_2,2),3);
	Gy_2=sum(sum(Gy_2,2),3);
	Gz_2=sum(sum(Gz_2,2),3);

	G2=zeros(3*length(Gy_2),1);
	G2(1:3:end)=Gy_2;
	G2(2:3:end)=Gx_2;
   G2(3:3:end)=-Gz_2;
end

if ten~=0
   
   Gx_3=wW.*wL.*Gx_3;
   Gy_3=wW.*wL.*Gy_3;
   Gz_3=wW.*wL.*Gz_3;
   
   Gx_3=sum(sum(Gx_3,2),3);
	Gy_3=sum(sum(Gy_3,2),3);
	Gz_3=sum(sum(Gz_3,2),3);

	G3=zeros(3*length(Gy_3),1);
	G3(1:3:end)=Gy_3;
	G3(2:3:end)=Gx_3;
   G3(3:3:end)=-Gz_3;
end


if ss==0
   G1=zeros(size(xloc,2)*3,1);
end
if ds==0
   G2=zeros(size(xloc,2)*3,1);
end
if ten==0
   G3=zeros(size(xloc,2)*3,1);
end


%change sign to be consistent with Okada's definition of slip
G2=-G2;

%%output displacments in matrix U
U=zeros(3,length(G1)/3);
%scale displacements with slip
temp=ss*G1+ds*G2+ten*G3;
U(1,:)=temp(1:3:end)';
U(2,:)=temp(2:3:end)';
U(3,:)=temp(3:3:end)';      %U is 3xn matrix of displacements


%add in elastic part for back-slip
[Ue,D,S,flag]=disloc3d(m',xloc,1,.25);%
%[Ue] = disloc(m',xloc(1:2,:),0.3);


U=U-Ue; 
          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M1,M2,M3]=momtensor_inverse(strike,dip,lam,mu)

strike=strike*pi/180;
dip=dip*pi/180;

Vs=[cos(strike) sin(strike) 0]';
Vd=[-cos(dip)*sin(strike) cos(dip)*cos(strike) sin(dip)]';
Vnorm=cross(Vd,Vs);

%moment tensor
%strike slip component
M1(1,1)=Vs(1)*Vnorm(1)*(lam+2*mu)+Vs(2)*Vnorm(2)*lam+Vs(3)*Vnorm(3)*lam;
M1(2,2)=Vs(1)*Vnorm(1)*lam+Vs(2)*Vnorm(2)*(lam+2*mu)+Vs(3)*Vnorm(3)*lam;
M1(3,3)=Vs(1)*Vnorm(1)*lam+Vs(2)*Vnorm(2)*lam+Vs(3)*Vnorm(3)*(lam+2*mu);
M1(1,2)=Vs(1)*Vnorm(2)*mu+Vs(2)*Vnorm(1)*mu;
M1(1,3)=Vs(1)*Vnorm(3)*mu+Vs(3)*Vnorm(1)*mu;
M1(2,3)=Vs(2)*Vnorm(3)*mu+Vs(3)*Vnorm(2)*mu;
M1(2,1)=M1(1,2);M1(3,1)=M1(1,3);M1(3,2)=M1(2,3);
%dip slip component
M2(1,1)=Vd(1)*Vnorm(1)*(lam+2*mu)+Vd(2)*Vnorm(2)*lam+Vd(3)*Vnorm(3)*lam;
M2(2,2)=Vd(1)*Vnorm(1)*lam+Vd(2)*Vnorm(2)*(lam+2*mu)+Vd(3)*Vnorm(3)*lam;
M2(3,3)=Vd(1)*Vnorm(1)*lam+Vd(2)*Vnorm(2)*lam+Vd(3)*Vnorm(3)*(lam+2*mu);
M2(1,2)=Vd(1)*Vnorm(2)*mu+Vd(2)*Vnorm(1)*mu;
M2(1,3)=Vd(1)*Vnorm(3)*mu+Vd(3)*Vnorm(1)*mu;
M2(2,3)=Vd(2)*Vnorm(3)*mu+Vd(3)*Vnorm(2)*mu;
M2(2,1)=M2(1,2);M2(3,1)=M2(1,3);M2(3,2)=M2(2,3);
%tensile component
M3(1,1)=Vnorm(1)*Vnorm(1)*(lam+2*mu)+Vnorm(2)*Vnorm(2)*lam+Vnorm(3)*Vnorm(3)*lam;
M3(2,2)=Vnorm(1)*Vnorm(1)*lam+Vnorm(2)*Vnorm(2)*(lam+2*mu)+Vnorm(3)*Vnorm(3)*lam;
M3(3,3)=Vnorm(1)*Vnorm(1)*lam+Vnorm(2)*Vnorm(2)*lam+Vnorm(3)*Vnorm(3)*(lam+2*mu);
M3(1,2)=Vnorm(1)*Vnorm(2)*mu+Vnorm(2)*Vnorm(1)*mu;
M3(1,3)=Vnorm(1)*Vnorm(3)*mu+Vnorm(3)*Vnorm(1)*mu;
M3(2,3)=Vnorm(2)*Vnorm(3)*mu+Vnorm(3)*Vnorm(2)*mu;
M3(2,1)=M3(1,2);M3(3,1)=M3(1,3);M3(3,2)=M3(2,3);

temp=M1./max(max(abs(M1)));
index=logical(abs(temp)<10^-13);
M1(index)=0;

temp=M2./max(max(abs(M2)));
index=logical(abs(temp)<10^-13);
M2(index)=0;

temp=M3./max(max(abs(M3)));
index=logical(abs(temp)<10^-13);
M3(index)=0;


function [sourceP4x4,sourceP2x2,halfspaceP4x4,halfspaceP2x2]=getprop(H,k,j,mu,lam,g,zs)
%subroutine to generate propagator matrices

A4x4=[0 k(j) 1/mu 0;-k(j)*lam/g 0 0 1/g;4*k(j)^2*mu*(lam+mu)/g 0 0 k(j)*lam/g;0 0 -k(j) 0];

%need to check the sign

global gravity
if isempty(gravity)
    gravity=1;
end
rg=gravity*3*10^-3;

G=[1 0 0 0;0 1 0 0;0 0 1 0;0 -rg 0 1];

z=0;
z0=H;
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;
P4x4(:,:)=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);
P2x2(1,:)=[cosh((z-z0)*abs(k(j))) (1/(mu*abs(k(j))))*sinh((z-z0)*abs(k(j)))];
P2x2(2,:)=[mu*abs(k(j))*sinh((z-z0)*abs(k(j))) cosh((z-z0)*abs(k(j)))];

halfspaceP4x4=G*P4x4;
halfspaceP2x2=P2x2;

%propagator from source to top of layer containing source
z=0;
z0=zs;
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;
A4x4=[0 k(j) 1/mu 0;-k(j)*lam/g 0 0 1/g;4*k(j)^2*mu*(lam+mu)/g 0 0 k(j)*lam/g;0 0 -k(j) 0];
P4x4zs=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);
P2x2zs(1,:)=[cosh((z-z0)*abs(k(j))) (1/(mu*abs(k(j))))*sinh((z-z0)*abs(k(j)))];
P2x2zs(2,:)=[mu*abs(k(j))*sinh((z-z0)*abs(k(j))) cosh((z-z0)*abs(k(j)))];

sourceP4x4=G*P4x4zs;
sourceP2x2=P2x2zs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%subroutine to generate bessel coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [UR_0,US_0,UR_n1,US_n1,UT_n1,UR_p1,US_p1,UT_p1,...
      UR_n2,US_n2,UT_n2,UR_p2,US_p2,UT_p2]=getbesselco(Ph2,Ph4,sourceP4x4,...
      sourceP2x2,M,g,k,K,j,lam,mu,B1,B2,H1,H2,t,T)

        
Pzs4=sourceP4x4;
Pzs2=sourceP2x2;

Mxx=M(1,1);
Myy=M(2,2);
Mzz=M(3,3);
Mxy=M(1,2);
Mxz=M(1,3);
Myz=M(2,3);

%point source
%order = 0
F1_0=0;
F2_0=Mzz/g;
F3_0=k(j)*((-Mxx-Myy)/2+lam*Mzz/g);

%order = -1
F1_n1=-1/2*(Mxz+i*Myz)/mu;
F2_n1=0;
F3_n1=0;
f1_n1=1/2*(Myz-i*Mxz)/mu;
f2_n1=0;

%order = 1
F1_p1=-1/2*(-Mxz+i*Myz)/mu;
F2_p1=0;
F3_p1=0;
f1_p1=-1/2*(Myz+i*Mxz)/mu;
f2_p1=0;

%order = -2
F1_n2=0;
F2_n2=0;
F3_n2=k(j)*(-(Myy-Mxx)/4+i*Mxy/2);
f1_n2=0;
f2_n2=k(j)*(i*(Mxx-Myy)/4-1/2*Mxy);  

%order = -2
F1_p2=0;
F2_p2=0;
F3_p2=k(j)*(-(Myy-Mxx)/4-i*Mxy/2);
f1_p2=0;
f2_p2=k(j)*(-i*(Mxx-Myy)/4-1/2*Mxy); 


   %order = 0
   F4x4=[F1_0;F2_0;F3_0;0];
    Pzsf=  Pzs4(3:4,:)*F4x4;
	
    	%compute roots
	[Nu11,Nu21,R,a] = roots_2visc([1 0],B1,B2,K,k(j),H1-H2,Ph4,mu);
    [Nu12,Nu22,R,a] = roots_2visc([0 1],B1,B2,K,k(j),H1-H2,Ph4,mu);
 
    Nu1 = Nu11*Pzsf(1) + Nu21*Pzsf(2);
    Nu2 = Nu12*Pzsf(1) + Nu22*Pzsf(2);
    
   
    s=R(2);
    Den = a*s*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term11 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term21 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(3);
    Den = a*s*(s-R(2))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term12 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term22 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(4);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term13 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term23 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(5);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term14 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term24 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(6);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term15 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term25 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(7);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term16 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term26 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(8);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term17 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term27 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
   u1 = term11 + term12 + term13 + term14+ term15+ term16+ term17;
   u2 = term21 + term22 + term23 + term24+ term25+ term26+ term27;
   
   US_0 = u1;
   UR_0 = -u2;
  
   %order = -1
    F4x4=[F1_n1;F2_n1;F3_n1;0];
    Pzsf=  Pzs4(3:4,:)*F4x4;
	
    Nu1 = Nu11*Pzsf(1) + Nu21*Pzsf(2);
    Nu2 = Nu12*Pzsf(1) + Nu22*Pzsf(2);
     
	s=R(2);
    Den = a*s*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term11 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term21 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(3);
    Den = a*s*(s-R(2))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term12 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term22 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(4);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term13 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term23 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(5);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term14 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term24 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(6);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term15 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term25 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(7);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term16 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term26 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(8);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term17 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term27 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
   u1 = term11 + term12 + term13 + term14+ term15+ term16+ term17;
   u2 = term21 + term22 + term23 + term24+ term25+ term26+ term27;
	   
   US_n1 = u1;
   UR_n1 = -u2;
  
  
%2x2 system
F2x2=[f1_n1;f2_n1];
Pzsf=Pzs2*F2x2;
%[N1,N2,N3,R2,a2] = roots_2visc_antiplane(Pzsf,B1,B2,k(j),H1-H2,Ph2,mu);  %negative sign becaue this is for k>0
[N11,N21,N31,R2,a2] = roots_2visc_antiplane([1 0],B1,B2,k(j),H1-H2,Ph2,mu);  %negative sign becaue this is for k>0
[N12,N22,N32,R2,a2] = roots_2visc_antiplane([0 1],B1,B2,k(j),H1-H2,Ph2,mu);  %negative sign becaue this is for k>0

N1 = N11*Pzsf(1) + N12*Pzsf(2);
N2 = N21*Pzsf(1) + N22*Pzsf(2);
N3 = N31*Pzsf(1) + N32*Pzsf(2);


    s=R2(2);
    Den = a2*s*(s-R2(3));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term1 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
     
    s=R2(3);
    Den = a2*s*(s-R2(2));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term2 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
    
    

UT_n1= -(term1 + term2);
 

   
   %order = 1
   F4x4=[F1_p1;F2_p1;F3_p1;0];
    Pzsf=  Pzs4(3:4,:)*F4x4;

    Nu1 = Nu11*Pzsf(1) + Nu21*Pzsf(2);
    Nu2 = Nu12*Pzsf(1) + Nu22*Pzsf(2);
        
	    s=R(2);
    Den = a*s*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term11 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term21 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(3);
    Den = a*s*(s-R(2))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term12 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term22 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(4);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term13 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term23 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(5);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term14 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term24 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(6);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term15 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term25 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(7);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term16 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term26 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(8);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term17 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term27 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
   u1 = term11 + term12 + term13 + term14+ term15+ term16+ term17;
   u2 = term21 + term22 + term23 + term24+ term25+ term26+ term27;
	   
   US_p1 = u1;
   UR_p1 = -u2;
   
   
%2x2 system
F2x2=[f1_p1;f2_p1];
Pzsf=Pzs2*F2x2;
%[N1,N2,N3,R2,a2] = roots_2visc_antiplane(Pzsf,B1,B2,k(j),H1-H2,Ph2,mu);

N1 = N11*Pzsf(1) + N12*Pzsf(2);
N2 = N21*Pzsf(1) + N22*Pzsf(2);
N3 = N31*Pzsf(1) + N32*Pzsf(2);

    s=R2(2);
    Den = a2*s*(s-R2(3));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term1 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
     
    s=R2(3);
    Den = a2*s*(s-R2(2));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term2 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;

UT_p1= -(term1 + term2);



   %order = -2
   F4x4=[F1_n2;F2_n2;F3_n2;0];
    Pzsf=  Pzs4(3:4,:)*F4x4;
	
    Nu1 = Nu11*Pzsf(1) + Nu21*Pzsf(2);
    Nu2 = Nu12*Pzsf(1) + Nu22*Pzsf(2);
    
	s=R(2);
    Den = a*s*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term11 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term21 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(3);
    Den = a*s*(s-R(2))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term12 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term22 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(4);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term13 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term23 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(5);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term14 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term24 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(6);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term15 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term25 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(7);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term16 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term26 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(8);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term17 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term27 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
   u1 = term11 + term12 + term13 + term14+ term15+ term16+ term17;
   u2 = term21 + term22 + term23 + term24+ term25+ term26+ term27;
	   
   US_n2 = u1;
   UR_n2 = -u2;
   
   
%2x2 system
F2x2=[f1_n2;f2_n2];
Pzsf=Pzs2*F2x2;
%[N1,N2,N3,R2,a2] = roots_2visc_antiplane(Pzsf,B1,B2,k(j),H1-H2,Ph2,mu);
 
N1 = N11*Pzsf(1) + N12*Pzsf(2);
N2 = N21*Pzsf(1) + N22*Pzsf(2);
N3 = N31*Pzsf(1) + N32*Pzsf(2);


    s=R2(2);
    Den = a2*s*(s-R2(3));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term1 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
     
    s=R2(3);
    Den = a2*s*(s-R2(2));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term2 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
      

UT_n2= -(term1 + term2);

    
   %order = 2
   F4x4=[F1_p2;F2_p2;F3_p2;0];
    Pzsf=  Pzs4(3:4,:)*F4x4;
	
    Nu1 = Nu11*Pzsf(1) + Nu21*Pzsf(2);
    Nu2 = Nu12*Pzsf(1) + Nu22*Pzsf(2);
    
	    s=R(2);
    Den = a*s*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term11 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term21 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(3);
    Den = a*s*(s-R(2))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term12 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term22 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(4);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(5))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term13 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term23 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(5);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(6))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term14 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term24 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
    
    s=R(6);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(7))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term15 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term25 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(7);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(8));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term16 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term26 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
    s=R(8);
    Den = a*s*(s-R(2))*(s-R(3))*(s-R(4))*(s-R(5))*(s-R(6))*(s-R(7));
    Num1 = -( Nu1(1) + Nu1(2)*s + Nu1(3)*s^2 + Nu1(4)*s^3 + Nu1(5)*s^4+ Nu1(6)*s^5+ Nu1(7)*s^6+ Nu1(8)*s^7 );
    Num2 = -( Nu2(1) + Nu2(2)*s + Nu2(3)*s^2 + Nu2(4)*s^3 + Nu2(5)*s^4+ Nu2(6)*s^5+ Nu2(7)*s^6+ Nu2(8)*s^7 );
    Tterm=1/(1-exp(s*T));
    term17 = T*Num1*s*exp(s*t)*Tterm/Den + Num1/Den;
    term27 = T*Num2*s*exp(s*t)*Tterm/Den + Num2/Den;
 
   u1 = term11 + term12 + term13 + term14+ term15+ term16+ term17;
   u2 = term21 + term22 + term23 + term24+ term25+ term26+ term27;
	   
   US_p2 = u1;
   UR_p2 = -u2;
  
%2x2 system
F2x2=[f1_p2;f2_p2];
Pzsf=Pzs2*F2x2;
%[N1,N2,N3,R2,a2] = roots_2visc_antiplane(Pzsf,B1,B2,k(j),H1-H2,Ph2,mu);
N1 = N11*Pzsf(1) + N12*Pzsf(2);
N2 = N21*Pzsf(1) + N22*Pzsf(2);
N3 = N31*Pzsf(1) + N32*Pzsf(2);
    
    s=R2(2);
    Den = a2*s*(s-R2(3));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term1 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;
     
    s=R2(3);
    Den = a2*s*(s-R2(2));
    Num = -( N1 + N2*s + N3*s^2);
    Tterm=1/(1-exp(s*T));
    term2 = T*Num*s*exp(s*t)*Tterm/Den + Num/Den ;

UT_p2= -(term1 + term2);

