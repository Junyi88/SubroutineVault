function [STRESSR]=RotateStress(Ro,Stress)

R=[Ro(1), Ro(2),Ro(3);
    Ro(4), Ro(5),Ro(6);
    Ro(7), Ro(8),Ro(9)];

ST=[Stress(1), Stress(2),Stress(3);
    Stress(4), Stress(5),Stress(6);
    Stress(7), Stress(8),Stress(9)];

STR=R.'*(ST*R);

STRESSR=[STR(1,1); STR(1,2); STR(1,3);
    STR(2,1); STR(2,2); STR(2,3);
    STR(3,1); STR(3,2); STR(3,3)];
    
end