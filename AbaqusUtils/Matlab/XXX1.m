clear;

x0=[-1 1;-1 1];
y0=[1 1;-1 -1];
z0=[0 0;0 0];

figure(1);
clf;
surf(x0,y0,z0);
xlabel('x')
ylabel('y')

theta=20;
Rot=[cosd(theta) -sind(theta);sind(theta) cosd(theta)];

x=x0;
y=y0;

for ny=1:2
    for nx=1:2
        G=Rot*[x0(ny,nx);y0(ny,nx)];
        x(ny,nx)=G(1);
         y(ny,nx)=G(2);
    end
end

figure(2);
clf;
surf(x,y,z0);
xlabel('x')
ylabel('y')

F=