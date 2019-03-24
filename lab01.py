def func(x, u):
    return x ** 2  + u ** 2

def euler(n, h, x, y):
    y_out = []
    for i in range(n):
        y += h * func(x, y)
        x += h
        y_out.append(y)
    return y_out

def un_euler(n, h, x, y):
    y_out = []
    for i in range(n):
        y += h * (func(x, y) + func(x+h, y + h * func(x, y)))/2  # 
        x += h
        y_out.append(y)
    return y_out

def runge_kutta_2(n, h, x, y):
    y_out = []
    for i in range(n):
        y += h * func(x+h/2, y+h/2*func(x,y))
        x += h
        y_out.append(y)
    return y_out

def picar(n, h, x):
    # Производные для метода Пикара
    def f1(a):
        return a ** 3 / 3
    def f2(a):
        return f1(a) + a ** 7 / 63
    def f3(a):
        return f2(a) +  a ** 14 * 2 / 2079 + a ** 15 / 59535


    y_out = [0]
    for i in range(n-1):
        x += h
        y_out.append(f3(x))
    return y_out
        

n = 2*10 ** 6
h = 10 ** -6
x = 0
y0 = 0
x_arr = [x + h*i for i in range(n)]
y1 = euler(n, h, x, y0)
y2 = un_euler(n, h, x, y0)
y3 = runge_kutta_2(n, h, x, y0)
y4 = picar(n, h, x)

print("|  x  |      Пикара   |      Эйлера    |   Неявный Эйлера  |    Рунге-Кутты 2го порядка|")
print("-"*85)
for i in range(len(y1)):
    print("|{:.2f} |  {:.8f}   |   {:.8f}   |    {:.8f}     |        {:.8f}         |".format(x_arr[i],y4[i],y1[i],y2[i],y3[i]))
print("-"*76)
