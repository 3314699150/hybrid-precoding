% Modeling and Analysis of Reconfigurable Intelligent Surfaces for Indoor and Outdoor Applications in Future Wireless Networks

clear all
close all

Nsym = 10000;   % ���г���
N = [16 64 128 256 512];   % RIS��Ԫ��
Dh = 50:20:250;     % ���Ͷ���RIS��X�����ϵľ���
Dv = 10;        % ���Ͷ���RIS��Y�����ϵľ���
Dsr = sqrt(Dh.^2 + Dv^2);    % ���Ͷ���RIS�ľ���
Dsd = 4.*Dh;      % ���Ͷ�����ն˵ľ���
Drd = sqrt((Dsd - Dh).^2 + Dv^2);    % RIS����ն˵ľ���
K1 = 1;  % ���Ͷ���RIS֮��·������˹����
K2 = 1;  % RIS����ն�֮��·������˹����
Gi = 1;  % RIS�����䲨�����ϵ�����
Gr = 1;  % RIS�ڷ��䲨�����ϵ�����
ep = 1;  % RIS�����źŹ���������źŹ��ʵı�ֵ
f1 = 2.4e9;   % ����Ƶ��Ϊ2.4GHz
% f2 = 28e9;   % ����Ƶ��Ϊ28GHz
N0 = 10^(-9.5)*1e-3;   % ��˹�������ĵ��߹������ܶ�
pt = 5;    % �����źŵĹ���


% RIS�������ŵ�
for k = 1:length(N)
    
    C_RIS = zeros(1,length(Dh));    % �ŵ�����
    % ���Ͷ���RIS��X�����ϵľ��벻ͬ
    for j = 1:length(Dh)
        
        % ���Ͷ���RIS֮����ŵ�
        Hsr_LOS = exp(-1i*2*pi*f1/3e8*Dsr(j))*ones(Nsym,N(k));      % ����Dsr���ȴ�������λ�仯��ȷ����ֵ
        Hsr_NLOS = sqrt(1/2)*(randn(Nsym,N(k)) + 1i*randn(Nsym,N(k)));     % ����˹�ֲ�CN~(0,1),N��RIS
        Hsr = sqrt(K1/(K1+1))*Hsr_LOS + sqrt(1/(K1+1))*Hsr_NLOS;     % ��˹�ֲ�
        angle_Hsr = angle(Hsr);     % �����Ͷ���RIS֮����ŵ������λ�仯
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS����ն�֮����ŵ�
        Hrd_LOS = exp(-1i*2*pi*f1/3e8*Drd(j))*ones(Nsym,N(k));      % ����Drd���ȴ�������λ�仯��ȷ����ֵ
        Hrd_NLOS = sqrt(1/2)*(randn(Nsym,N(k)) + 1i*randn(Nsym,N(k)));  % ����˹�ֲ�CN~(0,1),N��RIS
        Hrd = sqrt(K2/(K2+1))*Hrd_LOS + sqrt(1/(K2+1))*Hrd_NLOS;     % ��˹�ֲ�
        angle_Hrd = angle(Hrd);     % ��RIS����ն�֮����ŵ������λ�仯
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLrd = 10^(PLrd_dB/10);
        
        % RIS
        u = - angle_Hsr - angle_Hrd;  % ˲ʱ�����������λ����,(Nsym*N)����
        w = exp(u*1i);  % RIS��λ����
        
        % ��˹������
        n = sqrt(N0/2)*(randn(1,Nsym) + 1i*randn(1,Nsym));
        
        % ·�����
        PL = PLsr*PLrd;
        
        % ���յ����ź�
        receive = sqrt(pt)*sqrt(PL^-1)*sum(Hsr.*w.*Hrd,2)' + n;
        
        % �����
        r = (abs(sqrt(pt)*sqrt(PL^-1)*sum(Hsr.*w.*Hrd,2)')).^2/N0;
        
        % �ŵ�����
        c_temp = log2(1 + r);
        
        % �ŵ������ľ�ֵ
        C_RIS(1,j) = sum(c_temp)/Nsym;
    end
    % ��ͼ
    plot(Dh,C_RIS(1,:),'->r');
    hold on
end

K = 1;
% RIS�������ŵ��������Ͻ�
for k = 1:length(N)
    
    C_RIS_upper = zeros(1,length(Dh));    % �ŵ�����
    % ���Ͷ���RIS��X�����ϵľ��벻ͬ
    for j = 1:length(Dh)
        
        % ���Ͷ���RIS֮����ŵ�
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS����ն�֮����ŵ�
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLrd = 10^(PLrd_dB/10);
        
        % ·�����
        PL = PLsr*PLrd;
        
        % ����Ⱦ�ֵ
        E_A = N(k)*sqrt(PL^-1)*pi*(laguerreL(0.5,-K^2/(K+1)))^2/(4*(K+1));
        D_A = N(k)*PL^(-1)-N(k)*PL^(-1)*pi^2*(laguerreL(0.5,-K^2/(K+1)))^4/(16*(K+1)^2);
        r = pt/N0*(D_A+E_A^2);
        
        % �ŵ��������Ͻ�
        C_RIS_upper(1,j) = log2(1 + r);
        
    end
    % ��ͼ
    plot(Dh,C_RIS_upper(1,:),'-*g');
    hold on
end

K = 1;
% RIS�������ŵ��������Ͻ�
for k = 1:length(N)
    
    C_RIS_upper1 = zeros(1,length(Dh));    % �ŵ�����
    % ���Ͷ���RIS��X�����ϵľ��벻ͬ
    for j = 1:length(Dh)
        
        % ���Ͷ���RIS֮����ŵ�
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS����ն�֮����ŵ�
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % ·�����ģ��LOS��d�ĵ�λΪm��fc�ĵ�λΪGHz
        PLrd = 10^(PLrd_dB/10);
        
        % ·�����
        PL = PLsr*PLrd;
        
        % ����Ⱦ�ֵ
        E_a = sqrt(pi/(4*(K+1)))*exp(-K/2)*((1+2*K/2)*besseli(0,K/2)+2*K/2*besseli(1,K/2));
        E_b = sqrt(pi/(4*(K+1)))*exp(-K/2)*((1+2*K/2)*besseli(0,K/2)+2*K/2*besseli(1,K/2));
        E_a2 = K/(K+1)+1/(K+1);
        E_b2 = K/(K+1)+1/(K+1);
        E_A = sqrt(PL^-1)*N(k)*E_a*E_b;
        D_A = PL^-1*N(k)*(E_a2*E_b2-E_a^2*E_b^2);
        r = pt/N0*(D_A+E_A^2);
        % https://everything2.com/title/complex+Gaussian+distribution
        % https://onlinelibrary.wiley.com/doi/pdf/10.1002/(SICI)1522-2594(199903)41:3%3C614::AID-MRM26%3E3.0.CO;2-1
        
        % �ŵ��������Ͻ�
        C_RIS_upper1(1,j) = log2(1 + r);
        
    end
    % ��ͼ
    plot(Dh,C_RIS_upper1(1,:),'-.b');
    hold on
end

% ��ͼ����
axis([50 250 0 15])
ylabel('Achievable Rate [bit/s/Hz]')
xlabel('Distance Dh [m]')