clear

tx=35;
ty=20;
tz=16;

Rx=[1 0 0;0 cosd(tx) -sind(tx);0 sind(tx) cosd(tx)];
Ry=[cosd(ty) 0 sind(ty);0 1 0;-sind(ty) 0 cosd(ty)];
Rz=[cosd(tz) -sind(tz) 0;sind(tz) cosd(tz) 0;0 0 1];

R=Rx*Ry*Rz;

% R=R.';

G=eye(3)+R;
Gi=inv(G);
DSx=(Gi*(R-eye(3)));
X=R-eye(3);
DS=2.*DSx;

RR=inv(eye(3)-0.5.*DS)*(eye(3)+0.5*DS);