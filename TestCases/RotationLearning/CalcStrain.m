function [STR]=CalcStrain(F)


FP=[F(1) F(2) F(3);
    F(4) F(5) F(6);
    F(7) F(8) F(9)];

STRP=0.5.*[FP.'*FP-eye(3,3)];

STR=[STRP(1,1);STRP(1,2);STRP(1,3);
    STRP(2,1);STRP(2,2);STRP(2,3);
    STRP(3,1);STRP(3,2);STRP(3,3)];


end