# Использование SIMD-инструкций для ускорения вычислений
## Цель
Цель данной практической работы - исследовать влияние использования SIMD(Single Instruction, Multiple Data) инструкций на производительность программы на примере построения [множества Мандельброта](https://ru.wikipedia.org/wiki/%D0%9C%D0%BD%D0%BE%D0%B6%D0%B5%D1%81%D1%82%D0%B2%D0%BE_%D0%9C%D0%B0%D0%BD%D0%B4%D0%B5%D0%BB%D1%8C%D0%B1%D1%80%D0%BE%D1%82%D0%B0). 


## Множество Мандельброта
Множество Мандельброта строится следующим образом:
1. Задаём радиус вспомогательной окружности $r_{max}$. 
2. Выбираем точку с координатами $(x_0, y_0)$.
3. Рассматриваем точки с координатами $(x_n, y_n)$, которые определяются формулами 
$$x_n = x_{n-1}^2 - y_{n-1}^2 + x_0,$$
$$y_n = 2·x_{n-1}·y_{n-1} + y_0.$$
4. Пусть на N-ой итераци точка $(x_N, y_N)$ вышла за пределы окружности радиуса $r_{max}$. Окрасим точку $(x_0, y_0)$ в цвет, соответствующий номеру N.

_Примечание 1_. Согласно свойствам множества Мандельброта, если точка вышла за пределы окружности, наследующих итерациях она не попадёт внутрь этой окружности. Значит, N определяется однозначно.

_Примечание 2_. Согласно свойствам множества Мандельброта, если точка не вышла за пределы окружности за достаточно большое количество итераций (в нашей работе применим ограничение $N_{max} = 256$), то на каждой следующей итерации она такж останется внутри окружности.

Визуализация построенного таким образом множества Мандельброта:

<img src = "Mandelbrot_picture.png" width="800" height="600">

## Реализация
Для построения множества Мандельброта были написаны 4 программы:
1. "Наивная" реализация: расчёт номера N производится для каждого пикселя отдельно;
2. "Векторная" реализация: расчёт производится для массивов пикселей, идущих подряд (в одном массиве 8 пикселей);
3. "Векторная" реализация с применением прописанных вручную inline функций, аналогичных SIMD-инструкциям; 
4. Реализация с использовнияем SIMD-инструкций.

Каждая программа может быть запущена в одном из двух режимов:
1. Вывод изображения с указанием количества кадров в секунду (FPS) в левом верхнем углу;
2. Расчёт времени работы программы. В этом режиме изображение не генерируется с целью подсчёта времени, затраченного только на расчёт номера N для каждой точки. Измерение проводится 100 раз, затем среднее время расчётов выводится на экран.

Расчёт времени работы программы производится с помощью ассемблерной инструкции `rdtsc`, возвращающей количество тактов процессора, прошедших с последнего сброса. Разность значений этой инструкции до и после измеряемой программы - количество тактов процессора, затраченных на работу этой программы. Так как в данной практической работе сравнивается время работы разных программ, для чего достаточно найти отношение времени их работы, все результаты приводятся в количестве тактов процессора и не переводятся в другие единицы измерения.

### 1. "Наивная" реализация
В данной версии расчёты производятся для каждого пикселя отдельно. Сравним время работы программы с разными уровнями оптимизации.

Обозначения: $t_i$ - время работы программы в i-ом эксперименте, $t_{ср}$ - среднее время работы программы

| Уровень оптимизации | $t_1$     | $t_2$     | $t_3$     | $t_{ср}$  |
|---------------------|-----------|-----------|-----------|-----------|
| -O0                 | 894402050 | 884516047 | 904521722 | 894479940 |
| -O3                 | 420541364 | 425600682 | 430086117 | 425409387 |

Таким образом, применение ключа оптимизации `-O3` ускоряет время работы программы в $\frac{894479940}{425409387} ≈ 2,10$ раз.

### 2. "Векторная" реализация
В данной версии расчёты производятся для массивов из восьми идущих подряд пикселей. Используемый уровень оптимизации - `-O3`. Ожидается, что данная версия будет быстрее первой за счёт параллельной обработки нескольких значений.

Обозначения: $t_i$ - время работы программы в i-ом эксперименте, $t_{ср}$ - среднее время работы программы

| $t_1$     | $t_2$     | $t_3$     | $t_{ср}$  |
|-----------|-----------|-----------|-----------|
| 190844396 | 205729103 | 189575441 | 195382980 |

Относительно первой версии с тем же уровнем оптимизации получаем ускорение в $\frac{425409387}{195382980} ≈ 2,17$ раз.

### 3. "Векторная" реализация с добавлением функций
Эта версия отличается от предыдущей использованием inline функций, имитирующих SIMD-инструкции. Используемый уровень оптимизации - `-O3`. Ожидается, что данная версия будет быстрее первой, но аналогична по времени второй.

Обозначения: $t_i$ - время работы программы в i-ом эксперименте, $t_{ср}$ - среднее время работы программы

| $t_1$     | $t_2$     | $t_3$     | $t_{ср}$  |
|-----------|-----------|-----------|-----------|
| 256830101 | 258743844 | 265300915 | 260291620 |

Относительно первой версии с тем же уровнем оптимизации получаем ускорение в $\frac{425409387}{260291620} ≈ 1, 63$ раза. Заметим, что эта версия получилась медленнее второй.

### 4. Реализация с использованием SIMD-инструкций
Процессор, на котором выполнялась работа, поддерживает AVX-инструкции, то есть может работать с векторами в 256 бит. Это 8 элементов типа `folat`. Таким образом, ожидается, что в этом эксперименте скорость работы программы вырастет в 8 раз. Используемый уровень оптимизации - `-O3`. Также используется флаг компиляции `-mavx2` для успешного применения AVX-инструкций. 

Обозначения: $t_i$ - время работы программы в i-ом эксперименте, $t_{ср}$ - среднее время работы программы

| $t_1$    | $t_2$    | $t_3$    | $t_{ср}$ |
|----------|----------|----------|----------|
| 59526824 | 58985164 | 58953974 | 59155321 |

Относительно первой версии с тем же уровнем оптимизации получаем ускорение в $\frac{425409387}{59155321} ≈ 7, 19$ раз.

