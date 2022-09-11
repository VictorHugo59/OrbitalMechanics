%Gravitational Constant, we are working with international units. 
G=6.67E-11;
%mass of Sun and Earth
S=1.989E30;
E=5.972e34;
%we create the force field.
[x,y]=meshgrid(-2:0.5:2);
u = -G.*((S.*E)./(sqrt(x.^2+y.^2).^2)).*(x./(sqrt(x.^2+y.^2)));
v = -G.*((S.*E)./(sqrt(x.^2+y.^2).^2)).*(y./(sqrt(x.^2+y.^2)));
quiver(x,y,u,v)
hold on

%in order to solve the differential equation, will use the position verlet
%method, euler method isnt accurate enough
%One year time
tf=365.25*24*3600;
%time step
dt=360;
t=0:dt:tf;
vx=zeros(1,length(t));
vy=zeros(1,length(t));
%we create a x2 matrix because in position verlet we have 1/2 values.
y=zeros(1,length(t)*2);
x=zeros(1,length(t)*2);
%we can create an array that saves the radius
radio=zeros(1,length(t));
t(1)=0;
vx(1)=0;
vy(1)=30300; %Earth's velocity at perihelion
x(1)=147100632000; %Earth's distance at perihelion;
y(1)=0.000000000000;
radio(1)=147105052000;
vtot(1)=sqrt(vx(1).^2+vy(1).^2);
for i=1:length(t)
    x(2*i)=x(2*i-1)+vx(i).*0.5.*dt;
    y(2*i)=y(2*i-1)+vy(i).*0.5.*dt;
    vx(i+1)=vx(i)+(-G.*((S.*E)./(sqrt(x(2*i).^2+y(2*i).^2).^2)).*(x(2*i)./(sqrt(x(2*i).^2+y(2*i).^2)))./E)*dt;
    vy(i+1)=vy(i)+(-G.*((S.*E)./(sqrt(x(2*i).^2+y(2*i).^2).^2)).*(y(2*i)./(sqrt(x(2*i).^2+y(2*i).^2)))./E)*dt;
    x(2*i+1)=x(2*i)+vx(i+1).*0.5.*dt;
    y(2*i+1)=y(2*i)+vy(i+1).*0.5.*dt;
    %we calculate the radius of the orbit each iteration
    radio(i)=(sqrt(x(i).^2+y(i).^2));
    %velocity magnitude
    vtot(i+1)=sqrt(vx(i).^2+vy(i).^2);
end
hold on
plot(x,y,'color','b')
title("Trayectoria")
%major semiaxis (aphelion)
a=max(radio)
%minor semiaxis (perihelion)
b=min(radio)
%calculate excentricity
excentricidad = (max(radio)-min(radio))/(max(radio)+min(radio))
%we can create a sun for visualization
%position is accurate, but size has been modified for better visualization
rectangle('Position',[-696.34e6 -696.34e6 7e9 7e9],'Curvature',[1 1],'FaceColor',[1 1 0])