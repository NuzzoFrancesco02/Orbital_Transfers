function [r,v] = orbital_to_state(a,e,i,W,w,teta,mu)
R1 = [cos(w), sin(w), 0;-sin(w),cos(w),0;0,0,1];
R2 = [1, 0, 0; 0, cos(i), sin(i); 0 -sin(i) cos(i)];
R3 = [cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
R = R1*R2*R3;
p = a*(1-e^2);
vr = sqrt(mu/p).*e.*sin(teta);
vr_cart = [cos(teta); sin(teta)].*vr;
vtan = sqrt(mu/p).*(1+e.*cos(teta));
vtan_cart = [cos(teta+pi/2); sin(teta+pi/2)].*vtan;
v_segn = vr_cart+vtan_cart;
r = []; v = [];
for u = 1: length(teta)
    t = teta(u);
    r_new = [cos(t); sin(t); 0].*p./(1+e.*cos(t));
    r = [r R'*r_new];
    v = [v R'*[v_segn(1,u); v_segn(2,u);0]];
end