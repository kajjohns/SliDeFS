function [N1,N2,N3,D1,D2,D3] = antiplane_2visc_terms(PzsF,B1,B2,k,d1,P0,u)

Ph11=P0(1,1);Ph12=P0(1,2);
Ph21=P0(2,1);Ph22=P0(2,2);

C1=cosh(abs(k)*d1); Z1=sinh(abs(k)*d1);

Pf2=PzsF(2);

N1= -4*B1*B2*C1*Pf2*Ph11 + 4*B1^2*Pf2*Ph11*Z1;

N2= -2*B1*C1*Pf2*Ph11 - 2*B2*C1*Pf2*Ph11 + 4*B1*Pf2*Ph11*Z1 + 2*B1*C1*Pf2*Ph12*u*abs(k)...
    - 2*B2*Pf2*Ph12*u*Z1*abs(k);

N3= (-C1)*Pf2*Ph11 + Pf2*Ph11*Z1 + C1*Pf2*Ph12*u*abs(k) - Pf2*Ph12*u*Z1*abs(k);


D1= 4*B1*B2*C1*Ph21 - 4*B1^2*Ph21*Z1;

D2= 2*B1*C1*Ph21 + 2*B2*C1*Ph21 - 4*B1*Ph21*Z1 - 2*B1*C1*Ph22*u*abs(k) + 2*B2*Ph22*u*Z1*abs(k);

D3= C1*Ph21 - Ph21*Z1 - C1*Ph22*u*abs(k) + Ph22*u*Z1*abs(k);