function [dudn] = normalder(source, point, normal)


% source e point sono le coordinate rispettivamente di (vettori) sorgente e osservazione
% normal è il vettore contente le componenti della normale nel punto dove è calcolata la derivata


%r = norm(source-point);
pi2_r_2 = 2*pi*sum((source-point).^2);

% implementazione del calcolo della funzione derivata di U in 2D


dudx = (source(1)-point(1))/(pi2_r_2);
dudy = (source(2)-point(2))/(pi2_r_2);


% implementazione del calcolo della funzione derivata di U in 3D

%dudx = (source(1)-point(1))/(4*pi*r^3);
%dudy = (source(2)-point(2))/(4*pi*r^3);

dudn = normal(1)*dudx+normal(2)*dudy;