function [N1,N2,N3,R,D3] = roots_2visc_antiplane(PzsF,B1,B2,k,d1,P0,u)



[N1,N2,N3,D1,D2,D3] = antiplane_2visc_terms(PzsF,B1,B2,k,d1,P0,u);

R=roots([D3,D2,D1,0]);


