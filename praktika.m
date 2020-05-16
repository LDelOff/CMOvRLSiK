function praktika()
%Цифровые методы обработки в радиолокационных системах и комплексах
%% Пункт №1
%создаётся модель сигнала “шума” (распределение Рэлея), смесь ”сигнал+шум” (распределение Рэлея - Райса), параметры (проходили на практике) СКО шума, ампл. сигнала А.
%Среднее для Ш в квадратурах равно 0. Сигнал гармонический: частота (для одной точке по частоте, объём выборки по времени задаётся 100, 1000), фаза пост. начальная и по частоте изменяется в соответствии с заданной частотой.


T_all=1;

A=1;
f=10;
A_n=0.1;
T=0.001;

%[s,n]=punkt1(время наблюдения[c], период дискретизации[c], амплитуда сигнала[В], частота сигнала[Гц], максимальная амплитуда шума[В])
[s,n,t]=punkt1(T_all,T, A, f,A_n);

%% Пункт 2
%№2 – Алг., разрабатывается алгоритм (структура)- дам сегодня;
    %Приложение 1 (обобщённая схема на стр. 306).
    %а0=1, а1 = 1 для 0217;
    %а0=1, а1 = -1 для 0317.
    %Измеряется к-т передачи Ш и С+Ш (АЧХ, ФЧХ)
punkt2(s,n,t,T);


%№3 – Изм., пишется Программа: измерение моментов: среднего, среднего квадрата, дисперсии.
%Измеренное значение определяются как предел,
%соответствующих выборочных моментов (1 – го, 2 –го порядка).

end

function [s,n,t]=punkt1(time,T,A,f, A_n)
%% Задаю наблюдаемый промежуток времени
t=0:T:time;
%% Создаю сигнал
s=A*exp(1i*2*pi*f*t);
%% Шум
n=[];
for i=0:time/1000:time
    n(1,end+1)=A_n*((2*rand()-1)+1i*(2*rand()-1)); %Шум с равномерным распределением
    n(2,end)=0;
    for j=1:1000
        n(2,end)=n(2,end)+A_n*((2*rand()-1)+1i*(2*rand()-1)); %Шум с нормальным распределением
    end
    n(2,end)=n(2,end)/1000;
end
%% Формирование сигнала(графики)

    function grafiki1(s,n,t)
        figure(11)
        %% 
        subplot(2,1,1)
        plot(t,real(s))
        grid on
        title('s=Re[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% 
        subplot(2,1,2)
        plot(t,imag(s))
        grid on
        title('s=Im[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Шум (графики)
        figure(12)
        subplot(2,2,1)
        plot(t,real(n))
        grid on
        title('Re[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        ylim([-1 1])
        %% 
        subplot(2,2,2)
        plot(t,imag(n))
        grid on
        title('Im[n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        ylim([-1 1])
        %% сигнал 
        subplot(2,2,3)
        plot(t,real(s+n))
        grid on
        title('Im[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% 
        subplot(2,2,4)
        plot(t,imag(s+n))
        grid on
        title('Im[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %% Распределение шума
        figure(13)
        hist(real(n(1,:)),20)
        grid on
        title('Шум 1')
        figure(14)
        hist(real(n(2,:)),20)
        grid on
        title('Шум 2')
        
    end
grafiki1(s,n,t)
end

function punkt2(s,n,t,T)
%% Параметры фильтра
a0=1;
a1=1;
%% Разностное уравнение фильтра
    function y=filter(a0,a1,x0,x1)
        % x0=x_i, x1=x_i-1
        if isnan(x1) x1=0; end
        y=a0*x0+a1*x1;        
    end
%% ИХ фильтра
    function tr1(a0,a1)
        figure(21)
        delta=[0 1 0 0 0 0 0 0 0 0 0 0];
        h=[];
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
%tr1(a0,a1);
%% ЧХ фильтра
    function cfr1(a0,a1,T)
        figure(22)
        f=0:10^3;
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
        %% Коэффициент передачи (неправильно)
        %figure(5)
        %plot(f,abs(H).^2);
        %title('Коэффициент передачи');
        %xlabel('Частота, Гц')
        %ylabel('')
        %grid on
    end
%cfr1(a0,a1,T);

%% прохождение сингала через фильтр
y_s=0;
y_n=0;
y_sn=0;
for i=2:length(s)
    y_s=[y_s filter(a0,a1,s(i),s(i-1))];
    y_n=[y_n filter(a0,a1,n(i),n(i-1))];
    y_sn=[y_sn filter(a0,a1,(s(i)+n(i)),(s(i-1)+n(i-1)))];
end
%% Временные диаграммы
    function grafiki2(s,n,t,y_s,y_n,y_sn) 
        %% -----------------------------------------------------
        %% s(t)
        figure(24)
        %%
        subplot(2,2,1)
        plot(t,real(s))
        grid on
        title('s=Re[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(s))
        grid on
        title('s=Im[A*e^{j*2\pift}]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,3)
        plot(t,real(y_s))
        grid on
        title('Сигнал на выходе фильтра(действительная часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %%
        subplot(2,2,4)
        plot(t,imag(y_s))
        grid on
        title('Сигнал на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %% --------------------------------------
        %% n(t)
        figure(25)
        %%
        subplot(2,2,1)
        plot(t,real(n))
        grid on
        title('Re[n(t)]')
        xlabel('Время, с')
        ylabel('n(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(n))
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
        plot(t,real(s+n))
        grid on
        title('Re[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,2)
        plot(t,imag(s+n))
        grid on
        title('Im[s(t)+n(t)]')
        xlabel('Время, с')
        ylabel('s(t), В')
        %%
        subplot(2,2,3)
        plot(t,real(s+n))
        grid on
        title('Сигнал на выходе фильтра(действительная часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
        %%
        subplot(2,2,4)
        plot(t,imag(s+n))
        grid on
        title('Сигнал на выходе фильтра(мнимая часть)')
        xlabel('Время, с')
        ylabel('y(t), В')
    end
%grafiki2(s,n,t,y_s,y_n,y_sn);

end