function praktika7()
%% Created by L_DelOff
global type_of_noise time A f A_n1 A_n2 T N k
%% Пункт №1

%% Основные параметры
% сигнал
A=1;        % амплитуда сигнала[В]
f=50000;       % частота сигнала[Гц]
% модель
T=0.00001;    % Период дискретизации
time=10/f;    % всё время моделирования[с]
% когерентный накопитель
N=8;        % количество каналов при подсчёте ДПФ
k=4;        % интересующий нас канал = k+1
% шум
% максимальная амплитуда шума[В]
A_n1=0.01; %равномерного
A_n2=1000; %для нормального закона это множитель, который надо подгонять под цифру выше (шум 2)
type_of_noise=2; % выбор типа шума (1 или 2)

[s,n,t]=punkt1();
% s - массив значений сигнала
% t - массив значений времени
%% Пункт 2
y=punkt2(s,n,t);
% y = 
% [ выход ЦФ, если подать сигнал]
%%
%% Пункт 3
[report_all,report_sn]=punkt3(s,n,y)
% report_all =
% среднее(мат ожидание) СКО дисперсия
% x                     x   x
% ...                   ... ...
save('report_all.mat','report_all');
%% Внимание, эти строки нужны для построения зависимости SNR от амплитуды сигнала
%%1. Выполнить эти две строки для создания/обнуления таблицы
        %report_new=[]
        %save('report.mat','report_new');
%%2. Закоментировать 2 строки выше и раскоментировать дальше
    %первая строка - амплитуды
    %вторая строка - значение отношения сигнал шум
    load('report.mat','report_new');
    report_new(1,end+1)=A;
    report_new(2,end)=report_sn;
    save('report.mat','report_new');
%%
a=1;
fprintf('Конец');
end

function [s,n,t]=punkt1()
global type_of_noise time A f A_n1 A_n2 T
%% Задаю наблюдаемый промежуток времени
t=0:T:time;
%% Создаю комплексный гармонический сигнал
s=A*exp(1i*2*pi*f*t);
%% Создаю шум
n=[];
%n=[случайные числа (только rand())         ... это шум по равномерному распределению ...] 
%  [случайные числа (сложил много раз rand) ... это шум по нормальному распределению ... ]                               ]
for i=0:T:time
    n(1,end+1)=A_n1*((2*rand()-1)+1i*(2*rand()-1)); %Шум с равномерным распределением
    n(2,end)=0;
    for j=1:1000
        n(2,end)=n(2,end)+A_n1*((2*rand()-1)+1i*(2*rand()-1)); %Шум с нормальным распределением
    end
    n(2,end)=n(2,end)/1000*A_n2;
end
%% Формирование сигнала(графики)
    function grafiki1(s,n,t)
        figure(11)
        %% Чистый сигнал (действительная часть)
        subplot(2,1,1)
        plot(t,real(s))
        grid on
        title('s=Re[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Чистый сигнал (мнимая часть)
        subplot(2,1,2)
        plot(t,imag(s))
        grid on
        title('s=Im[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Шум (действительная часть)
        figure(12)
        subplot(2,2,1)
        plot(t,real(n(type_of_noise,:)))
        grid on
        title('Re[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        ylim([-1 1])
        %% Шум (мнимая часть)
        subplot(2,2,2)
        plot(t,imag(n(type_of_noise,:)))
        grid on
        title('Im[n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        ylim([-1 1])
        %% Полный сигнал (действительная часть) 
        subplot(2,2,3)
        plot(t,real(s+n(type_of_noise,:)))
        grid on
        title('Re[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Полный сигнал (мнимая часть) 
        subplot(2,2,4)
        plot(t,imag(s+n(type_of_noise,:)))
        grid on
        title('Im[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Распределение шума
        figure(13)
        hist(real(n(1,:)),20)
        grid on
        title('Шум(равномерный закон)')
        figure(14)
        hist(real(n(2,:)),20)
        grid on
        title('Шум(нормальный закон)')        
    end
%% Раскоментировать, если нужны графики
%grafiki1(s,n,t)
end

function y=punkt2(s,noise,t)
global T N k f A type_of_noise
%% Разностное уравнение фильтра
    function y=filter(x)
        sum=0;
        deltaw=2*pi/(N*T);
        eta=2*pi*f/deltaw;
        for n=0:N-1
            %sum=sum+x*exp(-1i*deltaw*T*n*k);
            sum=sum+A*exp(-1i*deltaw*T*n*(eta-k));
        end
        y=sum;        
    end
%% ЧХ фильтра
    function cfr1()
        figure(22)
        %% КЧХ
        deltaw=2*pi/(N*T);
        H=[];
        freq=1:100:1/T;
        for ii=freq
            eta=2*pi*ii/deltaw;
            fi=deltaw*T*(eta-k);
            H=[H sin(N*fi/2)/sin(fi/2)*exp(1j*(fi*(N-1))/2)];                
        end
        %% АЧХ
        subplot(2,1,1)
        plot(freq,abs(H));
        title('АЧХ');
        xlabel('Частота, Гц')
        ylabel('|H|, В/Гц')
        grid on
        %% ФЧХ
        subplot(2,1,2)
        plot(freq,angle(H));
        title('ФЧХ');
        xlabel('Частота, Гц')
        ylabel('\phi, \circ')
        grid on        
    end
%% Раскоментировать, если нужны графики
cfr1();
%% прохождение сингала через цифровой фильтр
y_s=[];
y_n=[];
y_sn=[];

for i=1:length(s)
    y_s=[y_s filter(s(i))]; % прогоняю чистый сигнал через ЦФ
    y_n=[y_n filter(noise(type_of_noise,i))]; % прогоняю чисто шум через ЦФ
    y_sn=[y_sn filter((s(i)+noise(type_of_noise,i)))]; % прогоняю смесь через ЦФ
end
y(1,:)=y_s;
y(2,:)=y_n;
y(3,:)=y_sn;
%% Временные диаграммы
    function grafiki2(s,t,y_s) 
        %% Аналогично реальные и мнимые части на входе и выходе ЦФ
        % Пределы по оси Y
        ylimit=[-1.5*A 1.5*A];        
        %% s(t) -----------------------------------
        figure(24)
        %%
        subplot(2,2,1)
        plot(t,real(s))
        grid on
        title('s=Re[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        ylim(ylimit);
        %%
        subplot(2,2,2)
        plot(t,imag(s))
        grid on
        title('s=Im[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        ylim(ylimit);
        %%
        subplot(2,2,3)
        plot(t,real(y_s))
        grid on
        title('Сигнал на выходе фильтра(действительная часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %ylim(ylimit);
        %%
        subplot(2,2,4)
        plot(t,imag(y_s))
        grid on
        title('Сигнал на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %ylim(ylimit); 
        %% --------------------------------------
        %% n(t)
        figure(25)
        %%
        subplot(2,2,1)
        plot(t,real(noise(type_of_noise,:)))
        grid on
        title('Re[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(noise(type_of_noise,:)))
        grid on
        title('Im[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        %%
        subplot(2,2,3)
        plot(t,real(y_n))
        grid on
        title('Шум на выходе фильтра(действительная часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %%
        subplot(2,2,4)
        plot(t,imag(y_n))
        grid on
        title('Шум на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %% --------------------------------------
        %% s(t)+n(t)
        figure(26)
        %%
        subplot(2,2,1)
        plot(t,real(s+noise(type_of_noise,:)))
        grid on
        title('Re[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(s+noise(type_of_noise,:)))
        grid on
        title('Im[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,3)
        plot(t,real(y_sn))
        grid on
        title('Сигнал на выходе фильтра(действительная часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %%
        subplot(2,2,4)
        plot(t,imag(y_sn))
        grid on
        title('Сигнал на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
    end
%% Раскоментировать, если нужны графики
grafiki2(s,t,y_s);
end

function [report_all,report_sn]=punkt3(s,noise,y)
global type_of_noise
    function report=izmerenie(x)
    %% Измерение 
    % Момент первого порядка - это математическое ожидание
    % Момент второго порядка - это дисперсия(квадрат среднеквадратического)
        fprintf('Математическое ожидание - ');
        MEAN=mean(x);
        fprintf(num2str(MEAN));
        fprintf(' В.\n');
        
        fprintf('Среднеквадратическое отклонение - ');
        RMS=rms(x);
        fprintf(num2str(RMS));
        fprintf(' В.\n');
        
        fprintf('Дисперсия - ');
        fprintf(num2str(RMS^2));
        fprintf(' В^2.\n\n');  
        
        report=[MEAN;RMS;RMS^2];
    end
report=[];
%% измерение отношения сигнал шум на выходе
ds=izmerenie(real(y(1,:)));
dn=izmerenie(real(y(2,:)));
report_sn=ds(2)/dn(2);
%% Эти строки для измерения параметров
%% --------
fprintf('Измерение параметров входного сигнала\n');
%% сигнал(вход)
fprintf('Чистый сигнал(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(s));
fprintf('Чистый сигнал(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(s));
%% шум(вход)
fprintf('Шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(noise(type_of_noise,:)));
fprintf('Шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(noise(type_of_noise,:)));
%% сигнал+шум(вход)
fprintf('Сигнал+шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(s+noise(type_of_noise,:)));
fprintf('Сигнал+шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(s+noise(type_of_noise,:)));
%% -------------
fprintf('Измерение параметров выходного сигнала\n');
%% сигнал(выход)
fprintf('Чистый сигнал(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(y(1,:)));
fprintf('Чистый сигнал(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(y(1,:)));
%% шум(выход)
fprintf('Шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(y(2,:)));
fprintf('Шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(y(2,:)));
%% сигнал+шум(выход)
fprintf('Сигнал+шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(y(3,:)));
fprintf('Сигнал+шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(y(3,:)));
report_all=report;
end


