function [D_v1, D_v2, df, Dt] = changeOrbitShape2(ai, ei, i, OM, om_i, th_i, af, ef, om_f,th_f, mu)
R3 = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R1 = [cos(om_i) sin(om_i) 0; -sin(om_i) cos(om_i) 0; 0 0 1];
R = R1*R2*R3;

p1 = ai*(1-ei^2);
p2 = af*(1-ef^2);
rp1 = p1/(1+ei); ra1 = p1/(1-ei);
rp2 = p2/(1+ef); ra2 = p2/(1-ef);
at = (rp1+ra2)/2; et = (ra2-rp1)/(ra2+rp1); pt = at*(1-et^2);
D_v1 = sqrt(2*mu*(1/rp1-1/(2*at)))-sqrt(2*mu*(1/rp1-1/(2*ai)));
D_v2 = sqrt(2*mu*(1/ra2-1/(2*af)))-sqrt(2*mu*(1/ra2-1/(2*at)));

stepTh = pi/180;

%th_t0 = acos(p1/(pt*et)-1+(pi*ei)/(pt*et)*cos(th_i));
th_t0 = acos((pt/p1-1+pt/p1*ei*cos(th_i))/et);
%th_tf = acos((pt/p2-1+pt/p2*ef*cos(th_f))/et);
[r1,v1] = kep2car(ai,ei,i,OM,om_i,0:stepTh:th_i,mu);
D_v_car1 = D_v1 .* R'* [-sin(th_i) cos(th_i) 0]';
v1(:,end) = v1(:,end) + D_v_car1;
[rt,vt] = kep2car(at,et,i,OM,om_f,th_i:stepTh:th_f,mu);
D_v_car2 = D_v2 .* R'* [-sin(th_f) cos(th_f) 0]';
vt(:,end) = vt(:,end) + D_v_car2;
[rf,vf] = kep2car(af,ef,i,OM,om_f,th_f:stepTh:2*pi+th_f,mu);
r = [r1 rt rf];
for i = 1 : length(r(1,:))
    plot3(r(1,:),r(2,:),r(3,:))
end
plot3(r1(1,:),r1(2,:),r1(3,:));
grid on;
hold on
Terra3d;
plot3(rt(1,:),rt(2,:),rt(3,:));
plot3(rf(1,:),rf(2,:),rf(3,:));
