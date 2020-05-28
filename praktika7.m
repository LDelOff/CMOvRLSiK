function praktika7()
%% Created by L_DelOff
global time A f T N k
%% Пункт №1
T=0.00001;    % Период дискретизации
%% Основные параметры
A=1;        % амплитуда сигнала[В]
f=1000;       % частота сигнала[Гц]
time=10/f;    % всё время моделирования[с]
N=8;        % количество каналов при подсчёте ДПФ
k=5;        % интересующий нас канал

[s,t]=punkt1();
% s - массив значений сигнала
% t - массив значений времени
%% Пункт 2
y=punkt2(s,t);
% y = 
% [ выход ЦФ, если подать сигнал      ]
%%
%% Пункт 3
report=punkt3(s,y,t);
% report =
% среднее(мат ожидание) СКО дисперсия
% x                     x   x
% ...                   ... ...

%% Внимание, эти строки нужны для построения зависимости SNR от амплитуды сигнала
%%Выполнить эти две строки для создания/обнуления таблицы
%report_new=[]
%save('report.mat','report_new');
%%Эти строки ломают код, но создают в папке таблицу(раскоментировать вместе со строками в п.3)
%первая строка - амплитуды
%вторая строка - значение отношения сигнал шум
%load('report.mat','report_new');
%report_new(1,end+1)=A;
%report_new(2,end)=report;
%save('report.mat','report_new');
%%
a=1;
fprintf('Конец');
end

function [s,t]=punkt1()
global time A f T
%% Задаю наблюдаемый промежуток времени
t=0:T:time;
%% Создаю комплексный гармонический сигнал
s=A*exp(1i*2*pi*f*t);
%% Формирование сигнала(графики)
    function grafiki1(s,t)
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
    end
%% Раскоментировать, если нужны графики
grafiki1(s,t)
end

function y=punkt2(s,t)
global T N k
%% Разностное уравнение фильтра
    function y=filter(x)
        sum=0;
        deltaw=2*pi/(N*T);
        for n=0:N-1
            sum=sum+x*exp(-1i*deltaw*T*n*k);
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
y_s=0;
for i=2:length(s)
    y_s=[y_s filter(s(i))]; % прогоняю чистый сигнал через ЦФ
end
y(1,:)=y_s;
%% Временные диаграммы
    function grafiki2(s,t,y_s) 
        %% Аналогично реальные и мнимые части на входе и выходе ЦФ
        % Пределы по оси Y
        ylimit=[-2 2];        
        %% s(t)
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
        ylim(ylimit);
        %%
        subplot(2,2,4)
        plot(t,imag(y_s))
        grid on
        title('Сигнал на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        ylim(ylimit);        
    end
%% Раскоментировать, если нужны графики
%grafiki2(s,t,y_s);

end

function report=punkt3(s,noise,y,t)
global type_of_noise
report=[]
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
%% измерение отношения сигнал шум на выходе  (либо это)
%ds=izmerenie(real(y(1,:)));
%dn=izmerenie(real(y(2,:)));
%report=ds(2)/dn(2);
%% Эти строки для измерения параметров (либо это, иначе может не заработать)
%{
%% --------
fprintf('Измерение параметров входного сигнала\n');
%% сигнал(вход)
fprintf('Чистый сигнал(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(s));
fprintf('Чистый сигнал(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(s));
%% шум(вход)
fprintf('Шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(n(type_of_noise,:)));
fprintf('Шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(n(type_of_noise,:)));
%% сигнал+шум(вход)
fprintf('Сигнал+шум(действительная составляющая): \n');
report(end+1,:)=izmerenie(real(s+n(type_of_noise,:)));
fprintf('Сигнал+шум(мнимая составляющая): \n');
report(end+1,:)=izmerenie(imag(s+n(type_of_noise,:)));
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
%}
end


