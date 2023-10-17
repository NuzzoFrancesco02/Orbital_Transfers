x = [1,0,0]; v = [0,10,0];
state_to_orbital(x,v,10);
teta = linspace(0,2*pi,500);
p = 7160; e = 0.0722;
r = p./(1+e.*cos(teta));
%plot3(r.*cos(teta),r.*sin(teta),r.*0)

%grid on;
%%
teta = 0:0.0001:2*pi;
a = 20000;

[r,v] = orbital_to_state(20000,0.6661,deg2rad(90),deg2rad(45),deg2rad(30),teta,398600);
t = (r(1,2:end)-r(1,1:end-1))./(v(1,2:end)-v(1,1:end-1));
h = cross(r,v);
e = cross(v,h)./398600-r./norm(r);


%%
[a,e,i,W,w,teta,v,e_vet] = state_to_orbital(r, v, 398600);
plot(1:size(v,2),e_vet)
a(1)
e(1)
i = rad2deg(i);
W = rad2deg(W);
w = rad2deg(w);
teta = rad2deg(teta);
%cross(v(:,1),h(:,j))/mu - r(:,j)./norm(r(:,j))

%%
%
plot3(r(1,:).*1e3,r(2,:).*1e3,r(3,:).*1e3,'LineWidth',2,'Color','r')
hold on
%Terra3d(6378)
%background('Milky Way');
hold on
planet3D('Earth');
light('Position',[1,-1,0]);
grid on;
figure()
grid on
subplot(3,1,1)
plot3(r(1,:),r(2,:),v(1,:),'LineWidth',2)
grid on
subplot(3,1,2)
plot3(r(1,:),r(2,:),v(2,:),'LineWidth',2)
grid on
subplot(3,1,3)
plot3(r(1,:),r(2,:),v(3,:),'LineWidth',2)
grid on
figure()
plot3(r(1,:).*1e3,r(2,:).*1e3,r(3,:).*1e3,'LineWidth',2,'Color','r')
hold on
planet3D('Earth');
light('Position',[1,-1,0]);
grid on
quiver3(r(1,1:10:end).*1e3,r(2,1:10:end).*1e3,r(3,1:10:end).*1e3,v(1,1:10:end).*1e2,v(2,1:10:end).*1e2,v(3,1:10:end),scale)
%%
r = [10000; 20000; 10000]; v = [-2.5; -2.5; 3];
mu = 398600;
[a,e,i,W,w,teta,v,e_vet] = state_to_orbital(r,v,mu);
a
e
rad2deg(i)
rad2deg(W)
rad2deg(w)
rad2deg(teta)
%%
teta = 0:pi/180:2*pi;
[r,v] = orbital_to_state(a,e,i,W,w,teta,mu);
figure()
planet3D('Earth');
light('Position',[1,-1,0]);
hold on
plot3(r(1,i).*1e3,r(2,i).*1e3,r(3,i).*1e3,'LineWidth',2,'Color','none')
for i = 1 : length(r(1,:))
    plot3(r(1,1:i).*1e3,r(2,1:i).*1e3,r(3,1:i).*1e3,'LineWidth',2,'Color','r')
    pause(0.1)
end

grid on
%%
a = 15000; e = 0.1; i = 15; W = 45; w = 30; teta = 180; mu = 398600;
[r,v] = orbital_to_state(a,e,deg2rad(i),deg2rad(W),deg2rad(w),deg2rad(teta),mu);
r,v;
p = a*(1-e^2);
r = p/(1+e*cos(deg2rad(teta)))
v_r = sqrt(mu/p)*(1+e*cos(deg2rad(teta)))
v_t = sqrt(mu/p)*(e*sin(deg2rad(teta)))
%%
mu = 398600.433;
r = [10000 20000 10000]'; v = [-2.5 -2.5 3]';
[a, e, i, OM, om, th] = car2kep(r,v,mu);
a 
e 

rad2deg(i)


rad2deg(OM)
rad2deg(om)
rad2deg(th)
[r,v] = kep2car(a, e, i, OM, om, th, mu);
r
v
th = 0:0.01:2*pi;
plotOrbit([a e i OM om th], mu);
%%
mu = 398600.433;
th = 0:0.01:2*pi;
r= kep2car(4e4, 0.8, 0, 0, 0, th, mu);
pt = interparc(1e2,r(1,:),r(2,:),r(3,:),'spline');
th = cart2pol(r(1,:),r(2,:),r(3,:));
pt = pt';
plot(pt(1,:),pt(2,:),'*')
axis('equal')
 %%
 n_orbit = 3;
 step_animation = 2;
 r= kep2car(4e4, 0.8, 0, 0, 0, th, mu);
 pt = pt';
 X = pt(1,:); Y = pt(2,:); Z = pt(3,:);
 Terra3d;
% Plot the 3D satellite orbit
plot3(X,Y,Z,'Color',"none");
grid on;

% Define an indefinite plot
v_view = [45 30 30];
h = plot3(nan,nan,nan,'or','MarkerFaceColor',"#77AC30",'MarkerEdgeColor',"#77AC30",'MarkerSize',10);
hold on;
traj = plot3(nan,nan,nan,'Color',"#A2142F");
view(v_view);
  e = 0.8
    t = (1/(1-e^2).^(3/2)).*(2.*atan(sqrt((1-e)./(1+e)).*tan(th.*0.5))-e.*sqrt(1-e^2).*sin(th)./(1+e.*cos(th)));
    t = [t(1) t];
    delta_t = (t(2:end)-t(1:end-1))*1e2;
    guess = find(delta_t<0);
    delta_t(guess) = delta_t(guess-1);
    [~,guess] = min(delta_t);
    delta_t(guess) = delta_t(guess+1);
    delta_t = delta_t.*1e-3;
for j = 1 :step_animation: n_orbit
    for i = 1:length(th)
        set(traj,'XData',X(1:i),'YData',Y(1:i),'ZData',Z(1:i));
        hold on;
        set(h,'XData',X(i),'YData',Y(i),'ZData',Z(i));
        drawnow
        
        pause(delta_t*1e2);
        
    end
    hold off
    
end
 
%%
teta = 0:0.1:2*pi;
for i = 1:length(teta)
    plot(teta(1:i),sin(teta(1:i)));
    pause(0.1)
end
%%
[r,v] = plotOrbit([100000 1.7 deg2rad(90) 0 0 0], 398600.433);
%%
%%
%plotOrbit([10000 0.3 90 0 0 0], 398600.433);
[r,v] = orbital_to_state(10000,0,deg2rad(90),0,0,180,398600.433);
r
%%
a = -100;
b = 100;
i = -100:100;
xi = 1./i(end).*i;
 x_nod = (a+b)/2 + (b-a)/2.*xi;
plot(x_nod,f(x_nod),'+')
%%
e = 0.8; p = 100; t = 0:0.01:360;
r  = p./(1+e.*cos(t));
plot(r.*cos(t),r.*sin(t),'*')
%%
R = 6378; ha = 278; hp = -196;
ra = ha + R; rp = hp + R; p = 2*rp*ra/(rp+ra); e = (ra-rp)/(ra+rp);
a = (rp+ra)/2;
mu = 398600.4418;
bitangent(6e4, 0.5, 0, 0, 0, 0, 8e4, 0.8, pi, pi, mu);
%%
