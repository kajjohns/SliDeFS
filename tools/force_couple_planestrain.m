 function [df1dx,df2dx,df1dz,df2dz] = force_couple_planestrain(x1,x2,xi1,xi2,nu,mu)    

%greens function -- concentrated point force, free surface, Segall chapter 3.3
  
%xi1,xi2 are source coordinates

%x2 is positive up
%need to shift coordinates down, far from free surface
x2 = x2-10^4;
xi2 = xi2-10^4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dx1 couples
dx1=0.001;

%first force in couple
xi1 = xi1+dx1;

    x2minus = x2-xi2;
    x2plus = x2+xi2;
    x1minus = x1-xi1;
    
    r1 = sqrt(x1minus.^2 + x2minus.^2);
    r2 = sqrt(x1minus.^2 + x2plus.^2);
    theta2 = atan(x1minus./x2plus);
        
    
 %line force in 1 direction   
        term1 = (3-4*nu)/4.*log(r1);
        term2 = (8*nu^2-12*nu+5)/4.*log(r2);
        term3 =  x2minus.^2./(4*r1.^2);
        term4 = ((3-4*nu).*x2plus.^2 + 2*xi2.*x2plus - 2*xi2.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x2plus.^2./r2.^4;
   g11_1 = -1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
        term1 = (1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = -xi2.*x2.*x1minus.*x2plus./r2.^4;
   g21_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3+term4);
    

        term1 = x1minus.^3./r1.^4;
        term2 = x1minus.*(x1minus.^2 - 4*xi2.*x2 - 2*xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x1minus./r1.^2 + 3*x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s11_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = x2plus.*(2*xi2.*x2 + x1minus.^2)./r2.^4;
        term3 = -8*xi2.*x2.*x1minus.^2.*x2plus./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x2plus.^2./r2.^4 );
   s12_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2.*x1minus./r1.^4;
        term2 = -x1minus.*(xi2.^2 - x2.^2 + 6*xi2.*x2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.^3./r2.^6;
        term4 = -(1-2*nu)/2*( x1minus./r1.^2 - x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s22_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

 %line force in 2 direction       
        term1 = -(1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = xi2.*x2.*x1minus.*x2plus./r2.^4;
    g12_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = -(3-4*nu)/4.*log(r1);
        term2 = -(8*nu^2-12*nu+5)/4.*log(r2);
        term3 = -x1minus.^2./(4*r1.^2);
        term4 = (2*xi2.*x2-(3-4*nu).*x1minus.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x1minus.^2./r2.^4;
    g22_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
    
    
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = (x2plus.*(x1minus.^2+2*xi2.^2)-2*xi2.*x1minus.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( -x2minus./r1.^2 + (3*xi2+x2)./r2.^2 + 4*x2.*x1minus.^2./r2.^4 );
   s11_2_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2./r1.^4;
        term2 = (x2.^2 -2*xi2.*x2 - xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( 1./r1.^2 - 1./r2.^2 + 4*x2.*x2plus./r2.^4 );
   s12_2_1 = -x1minus./(2*pi*(1-nu)).*(term1 + term2 + term3 + term4);
        term1 = x2minus.^3./r1.^4;
        term2 = x2plus.*(x2plus.^2 + 2*xi2.*x2)./r2.^4;
        term3 = -8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x1minus.^2./r2.^4 );
   s22_2_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

        
%second force in couple   
   
xi1 = xi1-dx1-dx1;

    x2minus = x2-xi2;
    x2plus = x2+xi2;
    x1minus = x1-xi1;
    
    r1 = sqrt(x1minus.^2 + x2minus.^2);
    r2 = sqrt(x1minus.^2 + x2plus.^2);
    theta2 = atan(x1minus./x2plus);
         
    
 %line force in 1 direction   
        term1 = (3-4*nu)/4.*log(r1);
        term2 = (8*nu^2-12*nu+5)/4.*log(r2);
        term3 =  x2minus.^2./(4*r1.^2);
        term4 = ((3-4*nu).*x2plus.^2 + 2*xi2.*x2plus - 2*xi2.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x2plus.^2./r2.^4;
   g11_2 = -1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
        term1 = (1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = -xi2.*x2.*x1minus.*x2plus./r2.^4;
   g21_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3+term4);
    

        term1 = x1minus.^3./r1.^4;
        term2 = x1minus.*(x1minus.^2 - 4*xi2.*x2 - 2*xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x1minus./r1.^2 + 3*x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s11_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = x2plus.*(2*xi2.*x2 + x1minus.^2)./r2.^4;
        term3 = -8*xi2.*x2.*x1minus.^2.*x2plus./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x2plus.^2./r2.^4 );
   s12_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2.*x1minus./r1.^4;
        term2 = -x1minus.*(xi2.^2 - x2.^2 + 6*xi2.*x2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.^3./r2.^6;
        term4 = -(1-2*nu)/2*( x1minus./r1.^2 - x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s22_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

 %line force in 2 direction       
        term1 = -(1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = xi2.*x2.*x1minus.*x2plus./r2.^4;
    g12_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = -(3-4*nu)/4.*log(r1);
        term2 = -(8*nu^2-12*nu+5)/4.*log(r2);
        term3 = -x1minus.^2./(4*r1.^2);
        term4 = (2*xi2.*x2-(3-4*nu).*x1minus.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x1minus.^2./r2.^4;
    g22_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
    
    
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = (x2plus.*(x1minus.^2+2*xi2.^2)-2*xi2.*x1minus.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( -x2minus./r1.^2 + (3*xi2+x2)./r2.^2 + 4*x2.*x1minus.^2./r2.^4 );
   s11_2_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2./r1.^4;
        term2 = (x2.^2 -2*xi2.*x2 - xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( 1./r1.^2 - 1./r2.^2 + 4*x2.*x2plus./r2.^4 );
   s12_2_2 = -x1minus./(2*pi*(1-nu)).*(term1 + term2 + term3 + term4);
        term1 = x2minus.^3./r1.^4;
        term2 = x2plus.*(x2plus.^2 + 2*xi2.*x2)./r2.^4;
        term3 = -8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x1minus.^2./r2.^4 );
   s22_2_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

   
   df1dx.s11 = (s11_1_1 - s11_1_2)/(2*dx1);
   df1dx.s12 = (s12_1_1 - s12_1_2)/(2*dx1);
   df1dx.s22 = (s22_1_1 - s22_1_2)/(2*dx1);
   df1dx.u1 = (g11_1 - g11_2)/(2*dx1);
   df1dx.u2 = (g21_1 - g21_2)/(2*dx1);
  
   df2dx.s11 = (s11_2_1 - s11_2_2)/(2*dx1);
   df2dx.s12 = (s12_2_1 - s12_2_2)/(2*dx1);
   df2dx.s22 = (s22_2_1 - s22_2_2)/(2*dx1);
   df2dx.u1 = (g12_1 - g12_2)/(2*dx1);
   df2dx.u2 = (g22_1 - g22_2)/(2*dx1);
  
   
   
   
xi1 = xi1+dx1+dx1;

        
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dx2 couples
dx2=0.001;

%first force in couple
xi2 = xi2+dx2;

    x2minus = x2-xi2;
    x2plus = x2+xi2;
    x1minus = x1-xi1;
    
    r1 = sqrt(x1minus.^2 + x2minus.^2);
    r2 = sqrt(x1minus.^2 + x2plus.^2);
    theta2 = atan(x1minus./x2plus);
        
    
 %line force in 1 direction   
        term1 = (3-4*nu)/4.*log(r1);
        term2 = (8*nu^2-12*nu+5)/4.*log(r2);
        term3 =  x2minus.^2./(4*r1.^2);
        term4 = ((3-4*nu).*x2plus.^2 + 2*xi2.*x2plus - 2*xi2.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x2plus.^2./r2.^4;
   g11_1 = -1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
        term1 = (1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = -xi2.*x2.*x1minus.*x2plus./r2.^4;
   g21_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3+term4);
    

        term1 = x1minus.^3./r1.^4;
        term2 = x1minus.*(x1minus.^2 - 4*xi2.*x2 - 2*xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x1minus./r1.^2 + 3*x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s11_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = x2plus.*(2*xi2.*x2 + x1minus.^2)./r2.^4;
        term3 = -8*xi2.*x2.*x1minus.^2.*x2plus./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x2plus.^2./r2.^4 );
   s12_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2.*x1minus./r1.^4;
        term2 = -x1minus.*(xi2.^2 - x2.^2 + 6*xi2.*x2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.^3./r2.^6;
        term4 = -(1-2*nu)/2*( x1minus./r1.^2 - x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s22_1_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

 %line force in 2 direction       
        term1 = -(1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = xi2.*x2.*x1minus.*x2plus./r2.^4;
    g12_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = -(3-4*nu)/4.*log(r1);
        term2 = -(8*nu^2-12*nu+5)/4.*log(r2);
        term3 = -x1minus.^2./(4*r1.^2);
        term4 = (2*xi2.*x2-(3-4*nu).*x1minus.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x1minus.^2./r2.^4;
    g22_1 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
    
    
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = (x2plus.*(x1minus.^2+2*xi2.^2)-2*xi2.*x1minus.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( -x2minus./r1.^2 + (3*xi2+x2)./r2.^2 + 4*x2.*x1minus.^2./r2.^4 );
   s11_2_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2./r1.^4;
        term2 = (x2.^2 -2*xi2.*x2 - xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( 1./r1.^2 - 1./r2.^2 + 4*x2.*x2plus./r2.^4 );
   s12_2_1 = -x1minus./(2*pi*(1-nu)).*(term1 + term2 + term3 + term4);
        term1 = x2minus.^3./r1.^4;
        term2 = x2plus.*(x2plus.^2 + 2*xi2.*x2)./r2.^4;
        term3 = -8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x1minus.^2./r2.^4 );
   s22_2_1 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

        
%second force in couple   
   
xi2 = xi2-dx2-dx2;

    x2minus = x2-xi2;
    x2plus = x2+xi2;
    x1minus = x1-xi1;
    
    r1 = sqrt(x1minus.^2 + x2minus.^2);
    r2 = sqrt(x1minus.^2 + x2plus.^2);
    theta2 = atan(x1minus./x2plus);
        
    
 %line force in 1 direction   
        term1 = (3-4*nu)/4.*log(r1);
        term2 = (8*nu^2-12*nu+5)/4.*log(r2);
        term3 =  x2minus.^2./(4*r1.^2);
        term4 = ((3-4*nu).*x2plus.^2 + 2*xi2.*x2plus - 2*xi2.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x2plus.^2./r2.^4;
   g11_2 = -1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
        term1 = (1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = -xi2.*x2.*x1minus.*x2plus./r2.^4;
   g21_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3+term4);
    

        term1 = x1minus.^3./r1.^4;
        term2 = x1minus.*(x1minus.^2 - 4*xi2.*x2 - 2*xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x1minus./r1.^2 + 3*x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s11_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = x2plus.*(2*xi2.*x2 + x1minus.^2)./r2.^4;
        term3 = -8*xi2.*x2.*x1minus.^2.*x2plus./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x2plus.^2./r2.^4 );
   s12_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2.*x1minus./r1.^4;
        term2 = -x1minus.*(xi2.^2 - x2.^2 + 6*xi2.*x2)./r2.^4;
        term3 = 8*xi2.*x2.*x1minus.^3./r2.^6;
        term4 = -(1-2*nu)/2*( x1minus./r1.^2 - x1minus./r2.^2 - 4*x2.*x1minus.*x2plus./r2.^4 );
   s22_1_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

 %line force in 2 direction       
        term1 = -(1-2*nu)*(1-nu).*theta2;
        term2 = x2minus.*x1minus./(4*r1.^2);
        term3 = (3-4*nu).*x2minus.*x1minus./(4*r2.^2);
        term4 = xi2.*x2.*x1minus.*x2plus./r2.^4;
    g12_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = -(3-4*nu)/4.*log(r1);
        term2 = -(8*nu^2-12*nu+5)/4.*log(r2);
        term3 = -x1minus.^2./(4*r1.^2);
        term4 = (2*xi2.*x2-(3-4*nu).*x1minus.^2)./(4*r2.^2);
        term5 = -xi2.*x2.*x1minus.^2./r2.^4;
    g22_2 = 1/(2*pi*mu*(1-nu))*(term1 + term2 + term3 + term4 + term5);
    
    
        term1 = x2minus.*x1minus.^2./r1.^4;
        term2 = (x2plus.*(x1minus.^2+2*xi2.^2)-2*xi2.*x1minus.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( -x2minus./r1.^2 + (3*xi2+x2)./r2.^2 + 4*x2.*x1minus.^2./r2.^4 );
   s11_2_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);
        term1 = x2minus.^2./r1.^4;
        term2 = (x2.^2 -2*xi2.*x2 - xi2.^2)./r2.^4;
        term3 = 8*xi2.*x2.*x2plus.^2./r2.^6;
        term4 = (1-2*nu)/2*( 1./r1.^2 - 1./r2.^2 + 4*x2.*x2plus./r2.^4 );
   s12_2_2 = -x1minus./(2*pi*(1-nu)).*(term1 + term2 + term3 + term4);
        term1 = x2minus.^3./r1.^4;
        term2 = x2plus.*(x2plus.^2 + 2*xi2.*x2)./r2.^4;
        term3 = -8*xi2.*x2.*x2plus.*x1minus.^2./r2.^6;
        term4 = (1-2*nu)/2*( x2minus./r1.^2 + (3*x2+xi2)./r2.^2 - 4*x2.*x1minus.^2./r2.^4 );
   s22_2_2 = -1/(2*pi*(1-nu))*(term1 + term2 + term3 + term4);

   
   df1dz.s11 = (s11_1_1 - s11_1_2)/(2*dx2);
   df1dz.s12 = (s12_1_1 - s12_1_2)/(2*dx2);
   df1dz.s22 = (s22_1_1 - s22_1_2)/(2*dx2);
   df1dz.u1 = (g11_1 - g11_2)/(2*dx2);
   df1dz.u2 = (g21_1 - g21_2)/(2*dx2);
  
   df2dz.s11 = (s11_2_1 - s11_2_2)/(2*dx2);
   df2dz.s12 = (s12_2_1 - s12_2_2)/(2*dx2);
   df2dz.s22 = (s22_2_1 - s22_2_2)/(2*dx2);
   df2dz.u1 = (g12_1 - g12_2)/(2*dx2);
   df2dz.u2 = (g22_1 - g22_2)/(2*dx2);
  
        
  
    
    