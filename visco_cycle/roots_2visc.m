function [N1,N2,R,D8] = roots_2visc(PzsF,B1,B2,K,k,d1,P0,u)


[D1,D2,D3,D4,D5,D6,D7,D8] = Denom_2visc(B1,B2,K,k,d1,P0,u);

[N1_1,N1_2,N1_3,N1_4,N1_5,N1_6,N1_7,N1_8] = Num1_2visc(PzsF,B1,B2,K,k,d1,P0,u);
N1=[N1_1,N1_2,N1_3,N1_4,N1_5,N1_6,N1_7,N1_8];

[N2_1,N2_2,N2_3,N2_4,N2_5,N2_6,N2_7,N2_8] = Num2_2visc(PzsF,B1,B2,K,k,d1,P0,u);
N2=[N2_1,N2_2,N2_3,N2_4,N2_5,N2_6,N2_7,N2_8];

R=roots([D8,D7,D6,D5,D4,D3,D2,D1,0]);

%[R1,P1,K1]=residue([N1_8 N1_7 N1_6 N1_5 N1_4 N1_3 N1_2 N1_1],[D8,D7,D6,D5,D4,D3,D2,D1,0]);

