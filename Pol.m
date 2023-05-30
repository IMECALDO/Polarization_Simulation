clc, clear, close all

% resolution
res = 100;
% lim
lim = 2;
% X 
x = linspace(-lim,lim,res);
% Y
y = linspace(-lim,lim,res);

% waist
w0 = 1;
% topological charge
m = 2;

[X,Y] = meshgrid(x,y);

% E
E = (sqrt(X.^2+Y.^2)/w0).^abs(m) .* exp((-(X.^2 + Y.^2))/(w0^2)) .* exp(1i*m*atan2(Y,X));

% Steps in the rotation of the waveplates
iterations_waveplates = 7;

% thetas variable
theta = linspace(0,2*pi,100);
theta_waveplates = linspace(0,pi/2,iterations_waveplates);

% Arrays to save Intensities in 
I1_array = zeros(length(theta),1);
I2_array = zeros(length(theta),2);
I3_array = zeros(iterations_waveplates,length(theta));
I4_array = zeros(iterations_waveplates,length(theta));

figure('Name','Beam','NumberTitle','off');
title("Beam")
imagesc(abs(E))
axis("square")

figure('Name','Fase of beam','NumberTitle','off');
title("Fase of beam")
imagesc(angle(E))
axis("square")

% functions for the matrices of Polarizer, quarter lambda waveplates, 
% & half lambda waveplates
Pol_mat = @(theta) [cos(theta).^2,sin(theta).*cos(theta);sin(theta).*cos(theta),sin(theta).^2];
Wplate_matqrt = @(theta) [1+1i*cos(2*theta),1i*sin(2*theta);1i*sin(2*theta),1-1i*cos(2*theta)]/sqrt(2);
Wplate_mathlf = @(theta) [cos(2*theta),sin(2*theta),sin(2*theta),-cos(2*theta)];

% Polarizes light horizontally
Polx = 1;
Poly = 0;

Ex = E.*Polx;
Ey = E.*Poly;

NormE = sum(abs(Ex(:)).^2+abs(Ey(:)).^2);
Ex = Ex/sqrt(NormE);
Ey = Ey/sqrt(NormE);

% different rotations of a linear polarizer
[Ex_t,Ey_t] = jones(Ex,Ey,Pol_mat,theta);
E_t = abs(Ex_t).^2 + abs(Ey_t).^2;
for i = 1:length(theta)
    I1_array(i) = sum(E_t(:,:,i),"all");
end

% Analitic mallus law
I0 = sum(abs(E),"all");

% Compare with analitic mallus law
figure('Name','Analitic mallus law and calculated','NumberTitle','off');
ax1 = subplot(2,1,1);
plot(theta,I1_array)
ylim([0,1])
title("Calculated")

ax2 = subplot(2,1,2);
plot(theta,I0*cos(theta).^2)
title("Analitic mallus law")


% plot format
xlim(ax1,[0,2*pi])
xlim(ax2,[0,2*pi])

set(ax1,'xticklabel',[]);
set(gca,'xtick',0:pi/2:2*pi)
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})

% Polarized beam through rotations of waveplates

[Ex_t,Ey_t] = jones(Ex,Ey,Wplate_mathlf,theta);
E_t = abs(Ex_t).^2 + abs(Ey_t).^2;
for i = 1:length(theta)
    I2_array(i,1) = sum(E_t(:,:,i),"all");
end

[Ex_t,Ey_t] = jones(Ex,Ey,Wplate_matqrt,theta);
E_t = abs(Ex_t).^2 + abs(Ey_t).^2;
for i = 1:length(theta)
    I2_array(i,2) = sum(E_t(:,:,i),"all");
end

% Intensities of polarized beam after different rotations of waveplates
figure('Name','Intensities of polarized beam after different rotations of waveplates','NumberTitle','off');

ax1 = subplot(2,1,1);
plot(theta,I2_array(:,1))
title("half lambda")

ax2 = subplot(2,1,2);
plot(theta,I2_array(:,2))
title("quarter lambda")

% plot format
xlim(ax1,[0,2*pi])
xlim(ax2,[0,2*pi])

set(ax1,'xticklabel',[]);
set(gca,'xtick',0:pi/2:2*pi)
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})

% Polarization of a polarized beam through diff rotations of waveplates of
% quarter lambda

figure('Name','Polarized beam after certain different quarter lambda waveplates rotation and after different rotations of polarization','NumberTitle','off');

for i = 1:length(theta_waveplates)
    % State of polarized beam after a waveplate
    [Ex_t,Ey_t] = jones(Ex,Ey,Wplate_matqrt,theta_waveplates(i));
    % Intensity of the beam after a several rotations of a polarizer
    [Ex_t_t,Ey_t_t] = jones(Ex_t,Ey_t,Pol_mat,theta);

    E_t = abs(Ex_t_t).^2 + abs(Ey_t_t).^2;
    for ii = 1:length(theta)
        I3_array(i,ii) = sum(E_t(:,:,ii),"all");
    end
    
    % Plot for each rotation of waveplate
    ax = subplot(iterations_waveplates,1,i);
    plot(theta,I3_array(i,:))
%     title("With " + string(theta_waveplates(i) + " rotation of quarter lambda waveplate"))
    title("Con " + string(theta_waveplates(i) + " radianes de la rotación del retardador de lambda/4"))

    set(ax,'xticklabel',[]);
    xlim([0,2*pi])
    ylim([0,1])
    ylabel("Intensidad")
end
xlabel("Rotación de polarizador lineal")

% plot format
set(gca,'xtick',0:pi/2:2*pi)
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})

% Polarization of a polarized beam through diff rotations of waveplates of
% half lambda

figure('Name','Polarized beam after certain different half lambda waveplates rotation and after different rotations of polarization','NumberTitle','off');

for i = 1:length(theta_waveplates)
    % State of polarized beam after a waveplate
    [Ex_t,Ey_t] = jones(Ex,Ey,Wplate_mathlf,theta_waveplates(i));
    % Intensity of the beam after a several rotations of a polarizer
    [Ex_t_t,Ey_t_t] = jones(Ex_t,Ey_t,Pol_mat,theta);

    E_t = abs(Ex_t_t).^2 + abs(Ey_t_t).^2;
    for ii = 1:length(theta)
        I4_array(i,ii) = sum(E_t(:,:,ii),"all");
    end
    
    % Plot for each rotation of waveplate
    ax = subplot(iterations_waveplates,1,i);
    plot(theta,I4_array(i,:))
%     title("With " + string(theta_waveplates(i) + " rotation of half lambda waveplate"))
    title("Con " + string(theta_waveplates(i) + " radianes de la rotación del retardador de lambda/2"))

    set(ax,'xticklabel',[]);
    xlim([0,2*pi])
    ylim([0,1])
    ylabel("Intensidad")
end
xlabel("Rotación de polarizador lineal")

% plot format
set(gca,'xtick',0:pi/2:2*pi)
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})

% Calculate stokes vectors of a right hand circle

% Create a left hand circle
[Ex_t,Ey_t] = jones(Ex,Ey,Wplate_matqrt,pi/4);
% Invert it 
[Ex_t,Ey_t] = jones(Ex_t,Ey_t,Wplate_mathlf,0);

%S0
[Ex_t_1,Ey_t_1] = jones(Ex_t,Ey_t,Pol_mat,0);
[Ex_t_2,Ey_t_2] = jones(Ex_t,Ey_t,Pol_mat,pi/2);
[Ex_t_3,Ey_t_3] = jones(Ex_t,Ey_t,Pol_mat,pi/4);
[Ex_t_4,Ey_t_4] = jones(Ex_t,Ey_t,Pol_mat,-pi/4);
[Ex_t_t,Ey_t_t] = jones(Ex_t,Ey_t,Wplate_matqrt,pi/4);
[Ex_t_5,Ey_t_5] = jones(Ex_t_t,Ey_t_t,Pol_mat,0);
[Ex_t_6,Ey_t_6] = jones(Ex_t_t,Ey_t_t,Pol_mat,pi/2);


I_S0 = sum(abs(Ex_t_1).^2 + abs(Ey_t_1).^2,"all") + sum(abs(Ex_t_2).^2 + abs(Ey_t_2).^2,"all");

I_S1 = sum(abs(Ex_t_1).^2 + abs(Ey_t_1).^2,"all") - sum(abs(Ex_t_2).^2 + abs(Ey_t_2).^2,"all");

I_S2 = sum(abs(Ex_t_3).^2 + abs(Ey_t_3).^2,"all") - sum(abs(Ex_t_4).^2 + abs(Ey_t_4).^2,"all");

I_S3 = sum(abs(Ex_t_5).^2 + abs(Ey_t_5).^2,"all") - sum(abs(Ex_t_6).^2 + abs(Ey_t_6).^2,"all");

P = sqrt(I_S1^2 + I_S2^2 + I_S3^2)/(I_S0);






