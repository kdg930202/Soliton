
%%% Multi-track simulation
%%% We will turn off all disp in GP2DF!!!!
% clearvars
% % close all

function Particle_num = threshold(ratio)
format compact
close all
tic


disp(date);
%%%%%% Define Constants
disp('Defining Constants');

%%% Physical Constants
hbar = 0.6582;	% Plank's Constant (meV.ps)
c=3*10^14*10^-12;	% Speed of light (um/ps)
me=511*10^6/c^2;	% Free electron mass (meV/c^2)
mp=(7*10^-5)*me; % mass of Cavity photon
%%% Calculation Contronl
tstart = 1;	% Initial Dynamics time (ps)
tfinish = 400;	% Finish time (ps)     % reservoir steady state time needed: 100 ps; +20 ps;
tp = tfinish/2;
tstep=0.1;%0.25;	% Frequency to stroe data
Nt=(tfinish-tstart)/tstep;
turn_Pr = 0;
turn_Pn = 1;
turn_defect = 0;
%%%%%% End

%%%%%% Grid of energy, momentum and real space
%%% For energy
Ns=400;	% Sub-steps;
% Although data will be store for every tstep, it will be calculate over smaller intervals, Ns sets the number of intervals for every tstep;
sstep=tstep/Ns;	% This time step in time used in solving the differential equations; (ps)


uc=1; % Unit cell size (um)
nuc=128; % Number of Unit Cell in one direction from 0-end


%%% For real space
num_p = 2;
Nx=nuc*num_p;	% Number of points in x-axis
Ny=nuc*num_p;	% Number of points in y-axis
xmax=nuc*uc;	% Set x range (um)
ymax=nuc*uc;	% set y range (um)
x=((1:Nx)'-Nx/2-1)*2*xmax/Nx';	% Grid in x
y=((1:Ny)'-Ny/2-1)*2*ymax/Ny';  % Grid in y
dx=abs(x(4)-x(3));	% dx (um)
dy=abs(y(4)-y(3));      % dy
x0=Nx/2+1;  % x=0 index
y0=Ny/2+1;
[X,Y]=meshgrid(x,y);

%%% For momentum space
kxmax=pi/dx; % maxmium of k valus (um^-1)
kymax=pi/dy;
kx=((1:Nx)'-Nx/2-1)*pi/xmax;	% Grid in kx
ky=((1:Ny)'-Ny/2-1)*pi/ymax;	% Grid in ky
dkx=2*kxmax/Nx;	% dkx
dky=2*kymax/Ny;	% dky
[KX,KY]=meshgrid(kx,ky);
%%%%%% END

%%%%%% Dispersion
EC=(hbar^2)*(KX.^2+KY.^2)/(2*mp);

%%%%% Nonlinear interaction %%%%%%%%
%%%%% Parameters controled by experimantal
alpha=2*10^-3;	%meV*um^2
Gin = 0.02;   % Gain rate 0.02 ps^-1 um^2
beta= 1*alpha;


%%%%%% Pump parameters: Laguerre Gaussian pumping
%%% Pump Rate Intensity no dxdy thing!!!!



omega = 0.15/hbar; %ps^-1
tauP = 14;   % ps
tauR = 120; % ps
Gamma = 1*hbar/(2*tauP);
% Gamma = 0;



% ratio = 10;
Inten_Pn = ratio*1/(Gin*tauP*tauR);
% Inten_Pr = 0.05*20;
Inten_Pr = 0.05*20;
% Inten_Pr = 0.05*100;
sigr = 6; % size
sigr_Poten = 2;
sigr_non = 10;

locate_pump = -50;
locate_poten = 10;
locate_non = 0;

RL = sqrt(X.^2+(Y-locate_pump).^2);
RL_non = sqrt(X.^2+(Y+locate_non).^2);
RL_poten = sqrt(X.^2+(Y+locate_poten).^2);
Phase = atan2(Y,X);

%Gaussian
r = sqrt(X.^2 + Y.^2);

kky = -0.5; %um^-1
% kky = -3.4;
Pr = turn_Pr*Inten_Pr*exp(-(RL/sigr).^2).*exp(-1i*kky*Y);
Pn = turn_Pn*Inten_Pn*exp(-(RL_non/sigr_non).^2);
% Pr(Y<locate_pump) = 0;
% Pr = Inten_Pr*exp(-(RL/sigr).^2);
PE = 10; %meV;
Poten = turn_defect*PE*exp(-(RL_poten/sigr_Poten).^2);


% system_name = sprintf('locatePn = %0.2f, locatePr = %0.2f, locateDefect = %0.2f.jpg',-locate_non,locate_pump,-locate_poten);
% figure()
% pcolor(x,y,abs(Pr));
% % pcolor(x,y,ss);
% hold on
% axis(2*[-32 32 -32 32])
% % pcolor(x,y,abs(Pn));
% viscircles([0,-locate_non],sigr_non,'Color','w');
% viscircles([0,-locate_poten],sigr_Poten);
% % title("Resonant","FontSize",20)
% axis on
% colorbar
% saveas(gcf,system_name)



%%


%%%%%% Define Fehlberg coefficients (for fifth order Runge-Kutta procedure)
ai(1)=16/135;
ai(2)=0;
ai(3)=6656/12825;
ai(4)=28561/56430;
ai(5)=-9/50;
ai(6)=2/55; %#ok<NASGU>

aibi(1)=1/360;
aibi(2)=0;
aibi(3)=-128/4275;
aibi(4)=-2197/75240;
aibi(5)=1/50;
aibi(6)=2/55;

dij(1,1)=0;
dij(2,1)=1/4;
dij(3,1)=3/32;          dij(3,2)=9/32;
dij(4,1)=1932/2197;     dij(4,2)=-7200/2197;    dij(4,3)=7296/2197;
dij(5,1)=439/216;       dij(5,2)=-8;            dij(5,3)=3680/513;
dij(5,4)=-845/4104;
dij(6,1)=-8/27;         dij(6,2)=2;             dij(6,3)=-3544/2565;
dij(6,4)=1859/4104;     dij(6,5)=-11/40;
%%%%%% END

% Load initial condition
%load psi_601.mat
%load res_601.mat
%psi_ini = psi5;
%rev_ini = rev5;
%clear psi5 rev5

%%%%%% Define the initial state

t=0;	% Real time (ps)
%%% Take the initial wavefunction
% set the random seed
rng('shuffle'); % make the initial seed depned on current time
%rng('default'); % make sure the same random initial condition


psi5=10e-10*(randn(Ny,Nx)+1i*randn(Ny,Nx));
% psi5 = 10e-1*(zeros(Ny,Nx)+1i*zeros(Ny,Nx));
rev5 = 10e-1*(zeros(Ny,Nx)+1i*zeros(Ny,Nx));
% noise = 10e-10*(randn(Ny,Nx)+1i*randn(Ny,Nx));
%%%%%% END






para = struct('hbar',hbar,'sstep',sstep,'alpha',alpha,'tauP',tauP, ...
    'Pr',Pr,'Pn',Pn,'omega',omega,'tp',tp,'Gamma',Gamma,'EC',EC,'Gin',Gin,'beta',beta,'tauR',tauR,'Poten',Poten);
%%%%%% 5th order Runge-Kutta Algorithm
for si=0:1:4
    t=si*sstep;

    % Label theta at next, current and 4 sub-timesteps behind
    psin1=['psi',sprintf('%g',si+5+1)];	% Label of next sub-timestep
    psin = ['psi',sprintf('%g',si+5)];      % Label of current sub-timestep
    psin_4 = ['psi',sprintf('%g',si+5-4)];  % Label of four sub-timesteps back

    % Label theta at next, current and 4 sub-timesteps behind
    revn1=['rev',sprintf('%g',si+5+1)];	% Label of next sub-timestep
    revn = ['rev',sprintf('%g',si+5)];      % Label of current sub-timestep
    revn_4 = ['rev',sprintf('%g',si+5-4)];  % Label of four sub-timesteps back


    % Calculate gradients
    [f1psi,f1rev] = fcal2D(...
        eval(psin),...
        eval(revn),...
        t,para);

    [f2psi,f2rev] = fcal2D(...
        eval(psin)+dij(2,1)*f1psi,...
        eval(revn)+dij(2,1)*f1rev,...
        t,para);

    [f3psi,f3rev] = fcal2D(...
        eval(psin)+dij(3,1)*f1psi+dij(3,2)*f2psi,...
        eval(revn)+dij(3,1)*f1rev+dij(3,2)*f2rev,...
        t,para);

    [f4psi,f4rev] = fcal2D(...
        eval(psin)+dij(4,1)*f1psi+dij(4,2)*f2psi+dij(4,3)*f3psi,...
        eval(revn)+dij(4,1)*f1rev+dij(4,2)*f2rev+dij(4,3)*f3rev,...
        t,para);

    [f5psi,f5rev] = fcal2D(...
        eval(psin)+dij(5,1)*f1psi+dij(5,2)*f2psi+dij(5,3)*f3psi+dij(5,4)*f4psi,...
        eval(revn)+dij(5,1)*f1rev+dij(5,2)*f2rev+dij(5,3)*f3rev+dij(5,4)*f4rev,...
        t,para);

    [f6psi,f6rev] = fcal2D(...
        eval(psin)+dij(6,1)*f1psi+dij(6,2)*f2psi+dij(6,3)*f3psi+dij(6,4)*f4psi+dij(6,5)*f5psi,...
        eval(revn)+dij(6,1)*f1rev+dij(6,2)*f2rev+dij(6,3)*f3rev+dij(6,4)*f4rev+dij(6,5)*f5rev,...
        t,para);

    % Calculate fields at next sub-timestep
    eval([psin1 '=eval(psin)+ai(1)*f1psi+ai(2)*f2psi+ai(3)*f3psi+ai(4)*f4psi+ai(5)*f5psi+ai(6)*f6psi;']);
    eval([revn1 '=eval(revn)+ai(1)*f1rev+ai(2)*f2rev+ai(3)*f3rev+ai(4)*f4rev+ai(5)*f5rev+ai(6)*f6rev;']);

    % Calculate errors
    errorpsi=sum(sum(abs(...
        aibi(1)*f1psi+aibi(2)*f2psi+aibi(3)*f3psi...
        +aibi(4)*f4psi+aibi(5)*f5psi+aibi(6)*f6psi)))...
        /(sum(sum(abs(eval(psin1)))));

    % Clear psi at sub-timestep 6 substeps behind next substep
    clear(psin_4);
    clear(revn_4);

end
clear f1psi f2psi f3psi f4psi f5psi f6psi
clear f1rev f2rev f3rev f4rev f5rev f6rev

clear ai aibi dij

%%%%%% END

%%%%%% Prepare for 5th order Adams-Bashforth-Moulton Algorithm
% This part of the code prepares for the 5th order Adams-Bashforth-Moulton Algorithm. From the stored values of psi from the Runge-Kutta procedure it calculates the gradients, fpsi=dpsi/dt, for each sub-timestep. The results are the arrays: % fpsi6, fpsi7, fpsi8, fpsi9 We also keep the last value of the field: (psi10, theta10) This part of the procedure uses the function, fcal2D(psi,...,t), which gives the value of dpsi/dt
% Increment substep by one
si=si+1;        %si=5 now

% Label psi at previous four sub-timesteps
psin_1 = ['psi',sprintf('%g',si+5-1)];
psin_2 = ['psi',sprintf('%g',si+5-2)];
psin_3 = ['psi',sprintf('%g',si+5-3)];
psin_4 = ['psi',sprintf('%g',si+5-4)];

% Label f at previous four sub-timesteps
fpsin_1 = ['fpsi',sprintf('%g',si+5-1)];
fpsin_2 = ['fpsi',sprintf('%g',si+5-2)];
fpsin_3 = ['fpsi',sprintf('%g',si+5-3)];
fpsin_4 = ['fpsi',sprintf('%g',si+5-4)];

% Label phi at previous four sub-timesteps
revn_1 = ['rev',sprintf('%g',si+5-1)];
revn_2 = ['rev',sprintf('%g',si+5-2)];
revn_3 = ['rev',sprintf('%g',si+5-3)];
revn_4 = ['rev',sprintf('%g',si+5-4)];

% Label f at previous four sub-timesteps
frevn_1 = ['frev',sprintf('%g',si+5-1)];
frevn_2 = ['frev',sprintf('%g',si+5-2)];
frevn_3 = ['frev',sprintf('%g',si+5-3)];
frevn_4 = ['frev',sprintf('%g',si+5-4)];


% Calculate f at previous four sub-timesteps
eval([ '[' fpsin_1  ',' frevn_1 ']=fcal2D(' ...
    psin_1  ',' revn_1 ','...
    't,para);']);

eval([ '[' fpsin_2  ',' frevn_2 ']=fcal2D(' ...
    psin_2 ',' revn_2 ','...
    't,para);']);

eval([ '[' fpsin_3  ',' frevn_3 ']=fcal2D(' ...
    psin_3 ',' revn_3 ','...
    't,para);']);

eval([ '[' fpsin_4  ',' frevn_4 ']=fcal2D(' ...
    psin_4 ',' revn_4 ','...
    't,para);']);

% Clear psi at sub-timesteps earlier than last sub-timestep
clear(psin_1); clear('psin_1');
clear(psin_2); clear('psin_2');
clear(psin_3); clear('psin_3');
clear(psin_4); clear('psin_4');

clear(revn_1); clear('revn_1');
clear(revn_2); clear('revn_2');
clear(revn_3); clear('revn_3');
clear(revn_4); clear('revn_4');


%%%%%% END

%%%%%% Start main 5th order Adams-Bashforth-Moulton Algorithm for all remaining timesteps This is the main part of the algorithm that calculates psi and theta for every sub-timestep.  For each timestep after the initial equilibrium time set by E_tstart plots are made for the first Monte Carlo realization if allow_plots==true For each value in the sample time grid, defined by E_Samples, E_start, and tfinish the field is saved in the variable psi_Ei where Ei is an index The time dynamics of the first Monte Carlo realization is also stored in the variable psit (if allow_plots==1) This part of the procedure uses the function, fcalc2DIS(psi,...,t), which gives the value of dpsi/dt
psit=zeros(Ny,Nx);
revt=zeros(Ny,Nx);


for ti=0:1:tfinish/tstep

    t=ti*tstep+si*sstep;
    plot_t(ti+1) = t;
    disp([sprintf('%1.4f',t),' ps', sprintf(' \t %g',ti),' step ']);
    display(['|psi|^2 : ' sprintf('%1.4g',sum(sum(abs(psit).^2)))]);	% Why 1e8?
    display(['|rev|^2 : ' sprintf('%1.4g',sum(sum(revt)))]);	% Why 1e8?
    sum_psi(ti+1) = sum(sum(abs(psit).^2));



    while si<Ns
        t=ti*tstep+si*sstep;

        % Label f at next and previous 5 sub-timesteps
        fpsin1 = ['fpsi',sprintf('%g',si+5+1)];
        fpsin = ['fpsi',sprintf('%g',si+5)];
        fpsin_1 = ['fpsi',sprintf('%g',si+5-1)];
        fpsin_2 = ['fpsi',sprintf('%g',si+5-2)];
        fpsin_3 = ['fpsi',sprintf('%g',si+5-3)];
        fpsin_4 = ['fpsi',sprintf('%g',si+5-4)];

        % Label f at next and previous 5 sub-timesteps
        frevn1 = ['frev',sprintf('%g',si+5+1)];
        frevn = ['frev',sprintf('%g',si+5)];
        frevn_1 = ['frev',sprintf('%g',si+5-1)];
        frevn_2 = ['frev',sprintf('%g',si+5-2)];
        frevn_3 = ['frev',sprintf('%g',si+5-3)];
        frevn_4 = ['frev',sprintf('%g',si+5-4)];



        % Label psi at next and last sub-timestep
        psin1 = ['psi',sprintf('%g',si+5+1)];
        psin = ['psi',sprintf('%g',si+5)];

        revn1 = ['rev',sprintf('%g',si+5+1)];
        revn = ['rev',sprintf('%g',si+5)];


        % Calculate f at last sub-timestep
        eval([ '[' fpsin  ',' frevn ']=fcal2D(eval(psin),eval(revn),t,para);']);

        % Predict psi at next sub-timestep (Adams-Bashforth Formula)
        psiapprox=eval(psin)+(1/720)*(...
            1901*eval(fpsin)-2774*eval(fpsin_1)...
            +2616*eval(fpsin_2)-1274*eval(fpsin_3)...
            +251*eval(fpsin_4));

        revapprox=eval(revn)+(1/720)*(...
            1901*eval(frevn)-2774*eval(frevn_1)...
            +2616*eval(frevn_2)-1274*eval(frevn_3)...
            +251*eval(frevn_4));

        clear(fpsin_4);
        clear(frevn_4);

        % Predict f at next sub-timestep

        eval([ '[' fpsin1  ',' frevn1 ']=fcal2D(psiapprox,revapprox,t,para);']);

        % Correct psi at next sub-timestep (Adams-Moulton Formula)
        eval([ psin1 '=eval(psin)+(1/720)*',...
            '(251*eval(fpsin1)+646*eval(fpsin)-264*eval(fpsin_1)',...
            '+106*eval(fpsin_2)-19*eval(fpsin_3));']);
         % Correct rev at next sub-timestep (Adams-Moulton Formula)
        eval([ revn1 '=eval(revn)+(1/720)*',...
            '(251*eval(frevn1)+646*eval(frevn)-264*eval(frevn_1)',...
            '+106*eval(frevn_2)-19*eval(frevn_3));']);



        clear(psin);
        clear(revn);

        % Increment sub-timestep
        si=si+1;
    end

    if eval(['isnan(sum(sum(' psin1 ')))'])
        calcsuccess=false;
        disp('  Calculation Stopped: Time step too large')
        return
    end

    % Prepare for next timestep
    psi5=eval(psin1); clear(psin1);
    rev5=eval(revn1); clear(revn1);

    clear(fpsin1);
    clear(frevn1);


    fpsi4=eval(fpsin); clear(fpsin);
    fpsi3=eval(fpsin_1); clear(fpsin_1);
    fpsi2=eval(fpsin_2); clear(fpsin_2);
    fpsi1=eval(fpsin_3); clear(fpsin_3);

    frev4=eval(frevn); clear(frevn);
    frev3=eval(frevn_1); clear(frevn_1);
    frev2=eval(frevn_2); clear(frevn_2);
    frev1=eval(frevn_3); clear(frevn_3);


    si=0;

    if ti>=tstart/tstep
        % store the result
        %         psiall=sum(sum(abs(psi5).^2));
        %         phiall=sum(sum(abs(phi5).^2));
        %         resall=sum(sum(rev5));
        save(['./result/psi_' sprintf('%g',ti-tstart/tstep+1)],'psi5','t');
        save(['./result/rev_' sprintf('%g',ti-tstart/tstep+1)],'rev5','t');

    end
    psit=psi5;
    revt=rev5;
    M(:,:,ti+1) = psit;


end
%%
Psi_M = abs(M).^2;

len_Psi = size(Psi_M);

for i = 1:len_Psi(3)
    Particle_num(i) = sum(sum(Psi_M(:,:,i)));
end


% figure()
% plot(Particle_num)
% saveas(gcf,'particle_num with ratio=%d.jpg', ratio)


% time = 0:tstep:tfinish;
% fig_name =  sprintf('Soliton : locatePn = %0.2f, locatePr = %0.2f, locateDefect = %0.2f.jpg',-locate_non,locate_pump,-locate_poten);
% figure()
% pcolor(X,Y,normal_max(:,:,end));
% viscircles([0,-locate_non],sigr_non,'Color','w');
% viscircles([0,-locate_poten],sigr_Poten);
% colorbar
% colormap hot
% shading flat
% axis equal
% title(['time : ' sprintf('%g',time(end)),'ps'],"FontSize",20);
% axis([-32 32 -32 32]);
% saveas(gcf,fig_name)
% 
% fig_name2 =  sprintf('Phase : locatePn = %0.2f, locatePr = %0.2f, locateDefect = %0.2f.jpg',-locate_non,locate_pump,-locate_poten);
% 
% figure()
% pcolor(X,Y,angle(M(:,:,end)));
% viscircles([0,-locate_non],sigr_non,'Color','w');
% viscircles([0,-locate_poten],sigr_Poten);
% colormap gray
% shading flat
% axis equal
% title(['time : ' sprintf('%g',time(end)),'ps'],"FontSize",20);
% axis([-32 32 -32 32]);
% 
% saveas(gcf,fig_name2)

end
