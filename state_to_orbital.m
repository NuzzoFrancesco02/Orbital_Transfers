function [a,e,i,W,w,teta,v,e_vet] = state_to_orbital(r, v, mu)
if size(r,1) ~= 3
    r = r'; %voglio i vari istanti di tempo lungo le colonne
end
if size(v,1) ~= 3
    v = v'; % stessa cosa
end

if size(r,2) == size(v,2)
    l = size(r,2);
else
    error('r e v hanno dimensioni diverse');
end

x = [1 0 0]';
y = [0 1 0]';
z = [0 0 1]';

e_vet = []; h = []; i = []; e = []; W = []; w = []; teta = []; v_sol = []; N = []; a = []; p = [];
for j = 1 : l
    % Calcolo momento della quantità di moto
    h = [h cross(r(:,j),v(:,j))];
    % Calcolo del diedro
    i = [i acos(dot(h(:,j),z)/norm(h(:,j)))];
    % Calcolo dell'eccentricità
    e_vet = [e_vet cross(v(:,j),h(:,j))/mu - r(:,j)./norm(r(:,j))];
    e = [e norm(e_vet(:,j))];
    % Calcolo del semiasse maggiore
    a = [a -mu./(2.*((norm(v(:,j))^2/2)-mu/norm(r(:,j))))];
    % Calcolo dell'asse dei nodi
    if norm(cross(z,h(:,j)))~=0
        N = [N cross(z,h(:,j))./norm(cross(z,h(:,j)))];
    else 
        N = [N [0 0 0]'];
    end
    
    if i(j) == 0
        W = [W 0];
        w = [w 0];
    else
        % Calcolo l'ascensione retta
        if dot(N(:,j),y)>=0
            W = [W acos(dot(x,N(:,j)))];
        elseif dot(N(:,j),y)<0
            W = [W 2*pi-acos(dot(x,N(:,j)))];
        end
        % Calcolo dell'anomalia del pericentro
        if dot(e_vet(:,j),z)>=0
        if norm(N(:,j)) ~= 0
            w = [w acos((dot(e_vet(:,j),N(:,j)))/(e(j)*norm(N(:,j))))];
        else
            w = [w 0];
        end
    elseif dot(e_vet(:,j),z)<0
        if norm(N(:,j)) ~= 0
            w = [w 2.*pi - acos((dot(e_vet(:,j),N(:,j)))/(e(j)*norm(N(:,j))))];
        else
            w = [w 0];
        end
    end
        
    
    
    
    end
    % Calcolo dell'anomalia vera
    if norm(r(:,j)) ~= 0
        vr = dot(v(:,j),r(:,j)./norm(r(:,j)));
    else
        vr = 0;
    end
    
    if vr >= 0
        if norm(r(:,j)) ~= 0
            teta = [teta acos((dot(e_vet(:,j),r(:,j)))/(e(j)*norm(r(:,j))))];
        else
            teta = [teta 0];
        end
    else
        if norm(r(:,j)) ~= 0
            teta = [teta 2*pi - acos((dot(e_vet(:,j),r(:,j)))/(e(j)*norm(r(:,j))))];
        else
            teta = [teta 0];
        end
    end
    p = [p a(j)*(1-e(j)^2)];
    v_r = sqrt(mu/p(j))*e(j)*sin(teta(j));
    v_teta = sqrt(mu/p(j))*(1+e(j)*cos(teta(j)));
    v_sol = [v_sol [v_r; v_teta]];
end

