% Modeling and Analysis of Reconfigurable Intelligent Surfaces for Indoor and Outdoor Applications in Future Wireless Networks

clear all
close all

Nsym = 10000;   % 序列长度
N = [16 64 128 256 512];   % RIS单元数
Dh = 50:20:250;     % 发送端与RIS在X方向上的距离
Dv = 10;        % 发送端与RIS在Y方向上的距离
Dsr = sqrt(Dh.^2 + Dv^2);    % 发送端与RIS的距离
Dsd = 4.*Dh;      % 发送端与接收端的距离
Drd = sqrt((Dsd - Dh).^2 + Dv^2);    % RIS与接收端的距离
K1 = 1;  % 发送端与RIS之间路径的莱斯因子
K2 = 1;  % RIS与接收端之间路径的莱斯因子
Gi = 1;  % RIS在入射波方向上的增益
Gr = 1;  % RIS在反射波方向上的增益
ep = 1;  % RIS发射信号功率与接收信号功率的比值
f1 = 2.4e9;   % 中心频率为2.4GHz
% f2 = 28e9;   % 中心频率为28GHz
N0 = 10^(-9.5)*1e-3;   % 高斯白噪声的单边功率谱密度
pt = 5;    % 发送信号的功率


% RIS辅助的信道
for k = 1:length(N)
    
    C_RIS = zeros(1,length(Dh));    % 信道容量
    % 发送端与RIS在X方向上的距离不同
    for j = 1:length(Dh)
        
        % 发送端与RIS之间的信道
        Hsr_LOS = exp(-1i*2*pi*f1/3e8*Dsr(j))*ones(Nsym,N(k));      % 传播Dsr长度带来的相位变化，确定的值
        Hsr_NLOS = sqrt(1/2)*(randn(Nsym,N(k)) + 1i*randn(Nsym,N(k)));     % 复高斯分布CN~(0,1),N个RIS
        Hsr = sqrt(K1/(K1+1))*Hsr_LOS + sqrt(1/(K1+1))*Hsr_NLOS;     % 莱斯分布
        angle_Hsr = angle(Hsr);     % 经发送端与RIS之间的信道后的相位变化
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS与接收端之间的信道
        Hrd_LOS = exp(-1i*2*pi*f1/3e8*Drd(j))*ones(Nsym,N(k));      % 传播Drd长度带来的相位变化，确定的值
        Hrd_NLOS = sqrt(1/2)*(randn(Nsym,N(k)) + 1i*randn(Nsym,N(k)));  % 复高斯分布CN~(0,1),N个RIS
        Hrd = sqrt(K2/(K2+1))*Hrd_LOS + sqrt(1/(K2+1))*Hrd_NLOS;     % 莱斯分布
        angle_Hrd = angle(Hrd);     % 经RIS与接收端之间的信道后的相位变化
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLrd = 10^(PLrd_dB/10);
        
        % RIS
        u = - angle_Hsr - angle_Hrd;  % 瞬时信噪比最大的相位条件,(Nsym*N)矩阵
        w = exp(u*1i);  % RIS相位控制
        
        % 高斯白噪声
        n = sqrt(N0/2)*(randn(1,Nsym) + 1i*randn(1,Nsym));
        
        % 路径损耗
        PL = PLsr*PLrd;
        
        % 接收到的信号
        receive = sqrt(pt)*sqrt(PL^-1)*sum(Hsr.*w.*Hrd,2)' + n;
        
        % 信噪比
        r = (abs(sqrt(pt)*sqrt(PL^-1)*sum(Hsr.*w.*Hrd,2)')).^2/N0;
        
        % 信道容量
        c_temp = log2(1 + r);
        
        % 信道容量的均值
        C_RIS(1,j) = sum(c_temp)/Nsym;
    end
    % 绘图
    plot(Dh,C_RIS(1,:),'->r');
    hold on
end

K = 1;
% RIS辅助的信道的容量上界
for k = 1:length(N)
    
    C_RIS_upper = zeros(1,length(Dh));    % 信道容量
    % 发送端与RIS在X方向上的距离不同
    for j = 1:length(Dh)
        
        % 发送端与RIS之间的信道
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS与接收端之间的信道
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLrd = 10^(PLrd_dB/10);
        
        % 路径损耗
        PL = PLsr*PLrd;
        
        % 信噪比均值
        E_A = N(k)*sqrt(PL^-1)*pi*(laguerreL(0.5,-K^2/(K+1)))^2/(4*(K+1));
        D_A = N(k)*PL^(-1)-N(k)*PL^(-1)*pi^2*(laguerreL(0.5,-K^2/(K+1)))^4/(16*(K+1)^2);
        r = pt/N0*(D_A+E_A^2);
        
        % 信道容量的上界
        C_RIS_upper(1,j) = log2(1 + r);
        
    end
    % 绘图
    plot(Dh,C_RIS_upper(1,:),'-*g');
    hold on
end

K = 1;
% RIS辅助的信道的容量上界
for k = 1:length(N)
    
    C_RIS_upper1 = zeros(1,length(Dh));    % 信道容量
    % 发送端与RIS在X方向上的距离不同
    for j = 1:length(Dh)
        
        % 发送端与RIS之间的信道
        PLsr_dB = 22*log10(Dsr(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLsr = 10^(PLsr_dB/10);
        
        % RIS与接收端之间的信道
        PLrd_dB = 22*log10(Drd(j)) + 28 + 20*log10(f1/1e9);    % 路径损耗模型LOS，d的单位为m，fc的单位为GHz
        PLrd = 10^(PLrd_dB/10);
        
        % 路径损耗
        PL = PLsr*PLrd;
        
        % 信噪比均值
        E_a = sqrt(pi/(4*(K+1)))*exp(-K/2)*((1+2*K/2)*besseli(0,K/2)+2*K/2*besseli(1,K/2));
        E_b = sqrt(pi/(4*(K+1)))*exp(-K/2)*((1+2*K/2)*besseli(0,K/2)+2*K/2*besseli(1,K/2));
        E_a2 = K/(K+1)+1/(K+1);
        E_b2 = K/(K+1)+1/(K+1);
        E_A = sqrt(PL^-1)*N(k)*E_a*E_b;
        D_A = PL^-1*N(k)*(E_a2*E_b2-E_a^2*E_b^2);
        r = pt/N0*(D_A+E_A^2);
        % https://everything2.com/title/complex+Gaussian+distribution
        % https://onlinelibrary.wiley.com/doi/pdf/10.1002/(SICI)1522-2594(199903)41:3%3C614::AID-MRM26%3E3.0.CO;2-1
        
        % 信道容量的上界
        C_RIS_upper1(1,j) = log2(1 + r);
        
    end
    % 绘图
    plot(Dh,C_RIS_upper1(1,:),'-.b');
    hold on
end

% 绘图设置
axis([50 250 0 15])
ylabel('Achievable Rate [bit/s/Hz]')
xlabel('Distance Dh [m]')