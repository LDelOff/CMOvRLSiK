function praktika()
%% Created by L_DelOff
global type_of_noise time A f A_n1 A_n2 T
%% Пункт №1 Задание:
%создаётся мод закона (шум 1)
A_n2=1000; %для нормального закона это множитель, который надо подгонять под цифру выше (шум 2)
type_of_noise=2; % выбор типа шума (1 или 2)

T=0.001;    % Период дискретизации
%% Выполнение пункта 1ель сигнала “шума” (распределение Рэлея), смесь ”сигнал+шум” (распределение Рэлея - Райса), параметры (проходили на практике) СКО шума, ампл. сигнала А.
%Среднее для Ш в квадратурах равно 0. Сигнал гармонический: частота (для одной точке по частоте, объём выборки по времени задаётся 100, 1000), фаза пост. начальная и по частоте изменяется в соответствии с заданной частотой.

%% Основные параметры
A=1;        % амплитуда сигнала[В]
f=10;       % частота сигнала[Гц]
time=10/f;    % всё время моделирования[с]

% максимальная амплитуда шума[В]
A_n1=0.01; %равномерного
[s,n,t]=punkt1();
% s - массив значений сигнала
% n - массив значений шума
% t - массив значений времени
%%
%% Пункт 2
%№2 – Алг., разрабатывается алгоритм (структура)- дам сегодня;
    %Приложение 1 (обобщённая схема на стр. 306).
    %а0=1, а1 = 1 для 0217;
    %а0=1, а1 = -1 для 0317.
    %Измеряется к-т передачи Ш и С+Ш (АЧХ, ФЧХ)
%% Выполняется пункт 2
y=punkt2(s,n,t);
% y = 
% [ выход ЦФ, если подать сигнал      ]
% [ выход ЦФ, если подать шум         ]
% [ выход ЦФ, если подать сигнал+шум  ]
%%
%% Пункт 3
%№3 – Изм., пишется Программа: измерение моментов: среднего, среднего квадрата, дисперсии.
%Измеренное значение определяются как предел,
%соответствующих выборочных моментов (1 – го, 2 –го порядка).
%% Выполняется пункт 3
report=punkt3(s,n,y,t);
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

function y=punkt2(s,n,t)
global type_of_noise f T
%% Параметры фильтра
a0=1;
a1=1;
%% Разностное уравнение фильтра
    function y=filter(a0,a1,x0,x1)
        % x0=x_i, x1=x_i-1
        % фильтр берёт предыдущие значения
        % если дать первое значение массива, то 
        % прога даст ошибку, что нет элемента с идексом 0
        % поэтому те несуществующие значения заменяем на нули
        if isnan(x1) x1=0; end 
        y=a0*x0+a1*x1;        
    end
%% ИХ фильтра
    function tr1(a0,a1)
        figure(21)
        % имитация единичного импульса, подавать буду начиная со второго
        % элемента, т.е. для фильтра первый ноль не учитывается
        % и сразу начнется с единицы
        delta=[0 1 0 0 0 0 0 0 0 0 0 0];
        h=[]; % импульсная характеристика
        for j=2:length(delta)
            h=[h filter(a0,a1,delta(j),delta(j-1))];
        end
        j=2:length(delta);
        stem(j-2,h,'filled','LineWidth',2);
        grid on
        hold on
        title('Импульсная характеристика ЦФ')
        xlabel('Отсчёты, бины')
        ylabel('h(t),В')
        xlim([-0.5 length(delta)-1.5]);
    end
%% Раскоментировать, если нужны графики
%tr1(a0,a1);
%% ЧХ фильтра
    function cfr1(a0,a1,T)
        figure(22)
        f=0:10^3;
        %% КЧХ - выражение берем из книжки
        H=a0+a1*exp(-1i*(2*pi*f)*T);
        %% АЧХ
        subplot(2,1,1)
        plot(f,abs(H));
        title('АЧХ');
        xlabel('Частота, Гц')
        ylabel('|H|, В/Гц')
        grid on
        %% ФЧХ
        subplot(2,1,2)
        plot(f,angle(H));
        title('ФЧХ');
        xlabel('Частота, Гц')
        ylabel('\phi, \circ')
        grid on
        %% Коэффициент передачи по мощности
        figure(23)
        plot(f,abs(H).^2);
        title('Коэффициент передачи');
        xlabel('Частота, Гц')
        ylabel('')
        grid on
    end
%% Раскоментировать, если нужны графики
%cfr1(a0,a1,T);

%% прохождение сингала через цифровой фильтр
y_s=0;
y_n=0;
y_sn=0;
for i=2:length(s)
    y_s=[y_s filter(a0,a1,s(i),s(i-1))]; % прогоняю чистый сигнал через ЦФ
    y_n=[y_n filter(a0,a1,n(type_of_noise,i),n(type_of_noise,i-1))]; % прогоняю чисто шум через ЦФ
    y_sn=[y_sn filter(a0,a1,(s(i)+n(type_of_noise,i)),(s(i-1)+n(type_of_noise,i-1)))]; % прогоняю смесь через ЦФ
end
y(1,:)=y_s;
y(2,:)=y_n;
y(3,:)=y_sn;

%% Временные диаграммы
    function grafiki2(s,n,t,y_s,y_n,y_sn) 
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
        %% --------------------------------------
        %% n(t)
        figure(25)
        %%
        subplot(2,2,1)
        plot(t,real(n(type_of_noise,:)))
        grid on
        title('Re[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(n(type_of_noise,:)))
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
        plot(t,real(s+n(type_of_noise,:)))
        grid on
        title('Re[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(s+n(type_of_noise,:)))
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
grafiki2(s,n,t,y_s,y_n,y_sn);

end

function report=punkt3(s,n,y,t)
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