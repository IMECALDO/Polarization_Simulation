function [Ex_t,Ey_t] = jones(Ex,Ey,mat,theta)
Ex_s = size(Ex);
Ey_s = size(Ex);

Ex_t = zeros(Ex_s(1),Ex_s(2),length(theta));
Ey_t = zeros(Ey_s(1),Ey_s(2),length(theta));

for i = 1:length(theta)
    mat_t = mat(theta(i));
    Ex_t(:,:,i) = mat_t(1)*Ex + mat_t(3)*Ey;
    Ey_t(:,:,i) = mat_t(2)*Ex + mat_t(4)*Ey;
end

