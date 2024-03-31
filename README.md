# Использование SIMD-инструкций для ускорения вычислений
## Цель
Цель данной практической работы - исследовать влияние использования SIMD(Single Instruction, Multiple Data) инструкций на производительность программы на примере построения [множества Мандельброта](https://ru.wikipedia.org/wiki/%D0%9C%D0%BD%D0%BE%D0%B6%D0%B5%D1%81%D1%82%D0%B2%D0%BE_%D0%9C%D0%B0%D0%BD%D0%B4%D0%B5%D0%BB%D1%8C%D0%B1%D1%80%D0%BE%D1%82%D0%B0). Для этого были реализованы 3 (4?) варианта построения множества Мандельброта:
1. Примитивная (наивная?) реализация: расчёт производится для каждого пикселя поотдельности;
2. "Векторная" реализация: расчёт производится для массивов пикселей, идущих подряд (в одном массиве 8 пикселей);
3. Реализация с использовнияем SIMD-инструкций.


## Множество Мандельброта
Множество Мандельброта строится следующим образом:
1. Задаём радиус вспомогательной окружности $r_{max}$. 
2. Выбираем точку с координатами $(x_0, y_0)$.
3. Рассматриваем точки с координатами $(x_n, y_n)$, которые определяются формулами 
$$x_n = x_{n-1}^2 - y_{n-1}^2 + x_0$$
$$y_n = 2\dotx_{n-1}\doty_{n-1} + y_0$$
4. Пусть на N-ой итераци точка $(x_N, y_N)$ вышла за пределы окружности радиуса $r_{max}$. Окрасим точку $(x_0, y_0)$ в цвет, соответствующий номеру N.


##