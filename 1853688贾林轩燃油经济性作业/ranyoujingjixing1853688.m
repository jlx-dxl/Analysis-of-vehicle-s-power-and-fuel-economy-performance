%% info
% 1853688 ������
% ����ȼ�;�������ҵ
%%
clear all
clc

generate_fig_one();
[uaIVnV,b] = generate_fig_two();
calculate(uaIVnV,b)
%% const 
function const = generate_const()
    const.n = 600:1:4000; %������ת��
    const.m = 3880;%������
    const.g = 9.8;%�������ٶ�
    const.r = 0.367;%���ְ뾶
    const.nT = 0.85;%����ϵ��еЧ��
    const.f = 0.013;%��������ϵ��
    const.CDA = 2.77;%��������ϵ��*ӭ�����
    const.i0 = 5.83;%��������������
    const.If = 0.218;%����ת������
    const.Iw1 = 1.798;%��ǰ��ת������
    const.Iw2 = 3.598;%�ĺ���ת������
    const.Iw = const.Iw1+const.Iw2;
    const.ig = [5.56,2.769,1.644,1.00,0.793];%������������
    const.L = 3.2;%���
    const.a = 1.947;%������ǰ��ľ���
    const.hg = 0.9;%���ĸ�
    const.rou = 0.74;%ȡ�����ܶ�0.74g/ml
    const.Qid = 0.299;%���ٺ�����
    const.G = const.m*const.g;
    const.delta = 1+1/const.m*const.Iw/const.r^2+1/const.m*(const.If*const.ig.^2*const.i0^2*const.nT)/(const.r^2);%������ת��������ϵ��
    const.B = [815,1326.8,-416.46,72.379,-5.8629,0.17768;
        1207,1354.7,-303.98,36.657,-2.0553,0.043072;
        1614,1284.4,-189.75,14.524,-0.51184,0.0068164;
        2012,1122.9,-121.59,7.0035,-0.18517,0.0018555;
        2603,1141.0,-98.893,4.4763,-0.091077,0.00068906;
        3006,1051.2,-73.714,2.8593,-0.05138,0.00035032;
        3403,1233.9,-84.478,2.9788,-0.047449,0.00028230;
        3804,1129.7,-45.291,0.71113,-0.00075215,-0.000038568];%���ϵ��
end



%% ��������ƽ��ͼ����
function generate_fig_one()
    config = generate_const();
    config.Ttq = -19.313+295.27*(config.n/1000)-165.44*(config.n/1000).^2+40.874*(config.n/1000).^3-3.8445*(config.n/1000).^4;%������ת��
    ig_length = length(config.ig);
    n_length = length(config.n);
    ua = zeros(ig_length,n_length);
    Ft = zeros(ig_length,n_length);
    for i = 1:ig_length
        ua(i,:)=0.377*config.r*config.n/(config.ig(i)*config.i0);
        Ft(i,:) = config.Ttq.*config.ig(i)*config.i0*config.nT/config.r;
    end%����λʱ�ټ�������
    Pe = Ft.*ua/(config.nT*3600);%�������������
    figure()
    title('��������ƽ��ͼ')
    for i = 1:ig_length
        plot(ua(i,:),Pe(i,:),'-')
        hold on
    end
    ua1 = 0:0.1:120;
    Ffw1 = config.CDA*ua1.^2/21.15+config.G*config.f;%��ʻ����
    Pfw1 = Ffw1.*ua1/(3600*config.nT);%��ʻ�����ڷ�������������ĵĹ���
    plot(ua1,Pfw1)
    grid on
    legend('I��','II��','III��','IV��','V��','(Pf+Pw)/nT')
    xlabel('ua/(km/h)')
    ylabel('P/(kW)')

end

%% ��ߵ��ʹθߵ��ٹ����ͺ�
function [uaIVnV,b] = generate_fig_two()
    config = generate_const();
    [c,~] = size(config.B);
    uaIVnV = zeros(c,2);
    PeIVnV = zeros(c,2);
    b = zeros(c,2);
    uaT = (25:10:85)';
    PeT = zeros(length(uaT),2);
    for i = [length(config.ig)-1,length(config.ig)]
        uaIVnV(:,i-length(config.ig)+2)=0.377*config.r*config.B(:,1)/(config.ig(i)*config.i0);
        PeIVnV(:,i-length(config.ig)+2) = (config.CDA*uaIVnV(:,i-length(config.ig)+2).^2/21.15+config.G*config.f).*uaIVnV(:,i-length(config.ig)+2)/(config.nT*3600);
        b(:,i-length(config.ig)+2) = config.B(:,2)+config.B(:,3).*PeIVnV(:,i-length(config.ig)+2)+config.B(:,4).*PeIVnV(:,i-length(config.ig)+2).^2+config.B(:,5).*PeIVnV(:,i-length(config.ig)+2).^3+config.B(:,6).*PeIVnV(:,i-length(config.ig)+2).^4;
        PeT(:,i-length(config.ig)+2) = (config.CDA*uaT.^2/21.15+config.G*config.f).*uaT/(config.nT*3600);
    end%��ת�ٶ�Ӧʱ�ټ��������Ĺ���

    bT = zeros(length(uaT),2);
    for i = 1:length(uaT)
        for j = [1,2]
            index = find(uaIVnV(:,j)-uaT(i)<0, 1, 'last' );
            bT(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(uaT(i)-uaIVnV(index,j))+b(index,j);
        end
    end%��ֵ�õ�ϵ�г����µĴθߵ�����ߵ���ȼ��������
    Qs = zeros(length(uaT),2);
    for i = [1,2]
        Qs(:,i) = PeT(:,i).*bT(:,i)./(1.02*uaT*config.rou*config.g);
    end%��λת��
    figure()
    title('�θߵ�����ߵ��ٹ����ͺ�')
    plot(uaT,Qs(:,1),'-s')
    hold on
    plot(uaT,Qs(:,2),'->')
    legend('�θߵ�','��ߵ�')
    grid on
    xlabel('ua/(km/h)')
    ylabel('Qs/(L/100km)')
    
end

%% ��JB3352-83�涨��������ѭ����ʻ����ٹ����ͺ�
function calculate(uaIVnV,b)
    config = generate_const();
    ua_acc_1 = 25:0.001:40;%��һ�μ���
    Pe_acc_1 = zeros(2,length(ua_acc_1));
    acc_1 = 0.25;%��һ�μ��ٶ�
    for i = [length(config.ig)-1,length(config.ig)]
        Pe_acc_1(i-length(config.ig)+2,:) = (config.CDA*ua_acc_1.^2/21.15+config.G*config.f+config.delta(i)*config.m*acc_1).*ua_acc_1/(config.nT*3600);
    end%��һ�μ��ٹ����еĹ���
    b_acc_1 = zeros(length(ua_acc_1),2);%��һ�μ��ٹ����еĺ�����
    for i = 1:length(ua_acc_1)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_acc_1(i)<0, 1, 'last' );
            b_acc_1(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_acc_1(i)-uaIVnV(index,j))+b(index,j);
        end
    end%��ֵ�õ�ϵ�г����µĴθߵ�����ߵ���ȼ��������
    Qs = Pe_acc_1.*b_acc_1'./(367.1*config.rou*config.g);
    Q_acc_1 = sum(Qs,2)*(0.001)/(3.6*acc_1);%��һ�μ���ʱIV��V���ĺ�����

    ua_acc_2 = 40:0.001:50;%�ڶ��μ���
    Pe_acc_2 = zeros(2,length(ua_acc_2));
    acc_2 = 0.2;%�ڶ��μ��ٶ�
    for i = [length(config.ig)-1,length(config.ig)]
        Pe_acc_2(i-length(config.ig)+2,:) = (config.CDA*ua_acc_2.^2/21.15+config.G*config.f+config.delta(i)*config.m*acc_1).*ua_acc_2/(config.nT*3600);
    end%�ڶ��μ��ٹ����еĹ���
    b_acc_2 = zeros(length(ua_acc_2),2);%�ڶ��μ��ٹ����еĺ�����
    for i = 1:length(ua_acc_2)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_acc_2(i)<0, 1, 'last' );
            b_acc_2(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_acc_2(i)-uaIVnV(index,j))+b(index,j);
        end
    end%��ֵ�õ�ϵ�г����µĴθߵ�����ߵ���ȼ��������
    Qs = Pe_acc_2.*b_acc_2'./(367.1*config.rou*config.g);
    Q_acc_2 = sum(Qs,2)*(0.001)/(3.6*acc_2);%�ڶ��μ���ʱIV��V���ĺ�����

    dua_acc_3 = 25;%�����μ���
    acc_3 = 0.36;%�����μ��ٶ�
    Q_acc_3 = config.Qid*dua_acc_3/(3.6*0.36);%����ʱIV��V���ĺ�����

    ua_e = [25,40,50];%����
    b_e = zeros(length(ua_e),2);%���ٹ����еĺ�����
    for i = 1:length(ua_e)
        for j = [1,2]
            index = find(uaIVnV(:,j)-ua_e(i)<0, 1, 'last' );
            b_e(i,j) = (b(index+1,j)-b(index,j))/(uaIVnV(index+1,j)-uaIVnV(index,j))*(ua_e(i)-uaIVnV(index,j))+b(index,j);
        end
    end%��ֵ�õ�ϵ�г����µĴθߵ�����ߵ���ȼ��������
    Pe_e = (config.CDA*ua_e.^2/21.15+config.G*config.f).*ua_e/(config.nT*3600);%���ٹ����еĹ���
    t_e = [7.2,22.5,18]';%������ʻ��ʱ��
    Q_e = sum(Pe_e'.*b_e.*t_e/(367.1*config.rou*config.g))';%���ٹ����еĺ�����

    Q = Q_acc_1+Q_acc_2+Q_acc_3+Q_e;%���������ͺ�

    Qs = 100*Q/1075;
    fprintf('IV���İٹ����ͺ�Ϊ %f��V���İٹ����ͺ�Ϊ %f',Qs(1),Qs(2));

end







