clear all
close all

%% input parameters

h=26.0;%DTu 10MW
g=9.81;

dt=0.1;
rho=1025;
CM=2.0;
Cd=1.0;
D=9.0; %Transition piece 8.5m, pile 9.0 
Nz=40;

%% Regular or irregular waves
if 0                                       % change to 0 for irregular waves
    % regular waves
    H=6; %was 1 changed to 6
    T=12;%was 11 changed to 12
    f=1/T;
    a=H/2;
    epsilon=0;
    TDur=2*T;
else
    % make a directional JONSWAP spectrum
    Hs=3.331; %DTU 10MW reference for 25m/s
    Tp=8.273; %DTU 10MW reference for 25m/s
    fp=1/Tp; 
    'Peak frequency of the spectra (Hz)', [fp]
    fHighCut=0.50;
    gamma_interval=Tp/sqrt(Hs);
    
    gamma = exp(5.75-1.15*gamma_interval); % 3.6<=gamma_test<=5.0
    %gamma=1 %test PM spectra
    %
    'Interval 3.6<=gamma_interval<5.0?', [gamma_interval]
    'Peak enhancement factor (gamma)', [gamma]
    
    
    TDur=600;
    df=1/TDur;
    f=[df:df:fHighCut]';
    
    % code the JONSWAP spectrum here. Right now it is just a uniform
    % spectrum.
    sigma = ones(size(f))*0.07;%TO DO need to change this to 0.09 for f>fp
    
    for i=1:size(f)
        if f(i)>fp 
            sigma(i)=0.09;
        end
    end
        
    %PM Spectra
    S_a = (0.3125*Hs^2*Tp).*(f./fp).^-5;
    S_b = exp(-1.25.*(f./fp).^-4);
    S_c = 1-0.287*log(gamma);
    S_d = exp(-0.5.*((f./fp-1)./sigma).^2);
    S = S_a.*S_b.*S_c.*gamma.^S_d;
    
    % scale spectrum so Hm0=4*std(eta)
    m0=trapz(f,S); % this is equal to sigma(eta)^2 which must be equal to (Hm0/4)^2
    S=S*((Hs/4)^2)/m0;
    
    % plot spectrum
    figure(1), clf, hold on
    plot(f,S)
    xlabel('f [Hz]')
    ylabel('S_{\eta} [m^2/Hz]');
    str = ['Wave Spectum, fs= ',num2str(fp),', Hs= ',num2str(Hs)];
    title(str)
    grid on
        
    % calc a and assign random phases
    a=sqrt(2*S*df);
    epsilon=rand(size(a))*2*pi;
end    

%% Solve the dispersion relation to get the wave numbers k

N=length(f);
omega=2*pi*f;
k=zeros(size(f));

myFun = @(k,omega,g,h) omega^2-g*k*tanh(k*h);
kGuess=omega(1)/sqrt(g*h); % we use shallow water limit as guess for first k
for j=1:N
    k(j) = fzero(@(x) myFun(x,omega(j),g,h) , kGuess);
    kGuess = k(j);         % we use solution for neighbooring frequency for all others 
end

%% Calculate free surface elevation and force in a loop

t = [0:dt:TDur]';
t=t(1:end-1);

Eta=zeros(size(t));
Fz = zeros(size(t));
u_bed = zeros(size(t));
u_free = zeros(size(t));
ut_bed = zeros(size(t));
ut_free = zeros(size(t));
% loop for calculation of free surface elevation and inline force
for j = 1:length(t)
    Eta(j) = sum(    a.*cos(omega.*t(j)+epsilon)    ) ;    
    zPhys=linspace(-h,Eta(j),Nz);
    if 0 % no Wheeler stretching. We use the physical z-values for calculation of kinematics. Change to 0 for Wheeler stretching.
        zCalc = zPhys;
    else % Wheeler stretching.              
        zCalc = h.*(zPhys-Eta)./(h+Eta);                       % Change zCalc here for question 5
    end
    % loop over zPhys. Calculate u, u_t and local force
    for l=1:length(zPhys)
        U= sum(a.*omega.*cosh(k.*(zCalc(l)+h))./sinh(k*h).*cos( omega*t(j)+epsilon));        
        Ut= sum(-a.*(omega.^2).*cosh(k.*(zCalc(l)+h))./sinh(k*h).*sin( omega*t(j)+epsilon)); % add acceleration here for question 4.3
        
        if(l==1)%first element in zPhys is at the seabed
            u_bed(j)=U;
            ut_bed(j)=Ut;
        end 
        z_free = length(zPhys);
        if(l==z_free)%last element in zPhys is at the free surface
            u_free(j)=U;
            ut_free(j)=Ut;
        end 
        dfx(l)=0.5*rho*Cd*D*abs(U)*U + rho*CM*pi*(D/2)^2*Ut;      % add inertia force here for question 4.3
        %dmy(l)= ...                       % add moment here for question 11
    end
    Fx(j)=trapz(zPhys,dfx);
    % My(j)= ...                           % complete this one for question 11 
end

%% Make some plots

figure(2), clf
subplot(2,1,1), hold on
plot(t,Eta,'b')
xlabel('t [s]')
ylabel('\eta [m]')

str = ['Free surface eleveation, fs= ',num2str(fp),', Hs= ',num2str(Hs)];
title(str)
grid on

if false
    subplot(3,1,2), hold on
    plot(t,u_bed,'b')
    xlabel('t [s]')
    ylabel('ubed [m/s]')
    title('Horizontal speed at seabed')
    axis([0 25 -3 3])
    grid on
end


if true
    subplot(2,1,2), hold on
    plot(t,Fx,'k')
    xlabel('t [s]')
    ylabel('Force [N]')
    str = ['Inline force, fs= ',num2str(fp),', Hs= ',num2str(Hs)];
    title(str)
    grid on
end

if false
subplot(3,1,3), hold on
plot(t,ut_bed,'b')
xlabel('t [s]')
ylabel('utbed [m/s^2]')
title('Horizontal acceleration at seabed')
axis([0 25 -1.2 1.2])
grid on

figure(3), clf
    subplot(3,1,1), hold on
    plot(t,Eta,'b')
    xlabel('t [s]')
    ylabel('\eta [m]')
    title('Free surface eleveation')
    grid on

    subplot(3,1,2), hold on
    plot(t,u_free,'b')
    xlabel('t [s]')
    ylabel('ufree [m/s]')
    title('Horizontal speed at free surface')
    axis([0 25 -3 3])
    grid on

    subplot(3,1,3), hold on
    plot(t,ut_free,'b')

    xlabel('t [s]')
    ylabel('utfree [m/s^2]')    
    title('Horizontal acceleration at free surface')
    axis([0 25 -1.2 1.2])
    grid on
end



%subplot(3,1,3), hold on
%plot(t,My,'r')
%xlabel('t [s]')
%ylabel('Moment [Nm]')
%grid on

%% Make some statistics
'eta Fx My'

'mean   ', [mean(Eta) mean(Fx) ]
'std    ', [std(Eta)  std(Fx) ]
'min    ', [min(Eta)  min(Fx) ]
'max    ', [max(Eta)  max(Fx) ]


'Hs and 4*std(Eta)', [Hs 4*std(Eta) ]

