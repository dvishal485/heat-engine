import tkinter as tk
from thermodynamics_logic import *

r = tk.Tk()
dest = []
globalDest = []
r.title('Thermodynamic Process')
initial = [tk.Entry(r), tk.Entry(r), tk.Entry(r)]
final = [tk.Entry(r), tk.Entry(r), tk.Entry(r)]
y_entry, k_entry = (tk.Entry(r), tk.Entry(r))


def clear():
    for i in globalDest:
        try:
            i.grid_remove()
        except:
            None
    initial[0].delete(0, tk.END)
    initial[1].delete(0, tk.END)
    initial[2].delete(0, tk.END)
    final[0].delete(0, tk.END)
    final[1].delete(0, tk.END)
    y_entry.delete(0, tk.END)
    k_entry.delete(0, tk.END)
    final[2].delete(0, tk.END)


def prcs(*args):
    for i in dest:
        try:
            i.grid_remove()
        except:
            None
    dest.clear()
    asp = process.get()[0]
    if asp == 'I':
        k = tk.Label(r, text='k')
        k.grid(row=12, column=1)
        k_entry.grid(row=12, column=2)
    elif asp == 'A':
        y = tk.Label(r, text='gamma')
        y.grid(row=11, column=1)
        y_entry.grid(row=11, column=2)
        k = tk.Label(r, text='k')
        k.grid(row=12, column=1)
        k_entry.grid(row=12, column=2)
    elif asp == 'P':
        y = tk.Label(r, text='y')
        y.grid(row=11, column=1)
        y_entry.grid(row=11, column=2)
        k = tk.Label(r, text='k')
        k.grid(row=12, column=1)
        k_entry.grid(row=12, column=2)
    dest.append(k)
    dest.append(k_entry)
    try:
        dest.append(y)
        dest.append(y_entry)
    except:
        None


def submit():
    alpUnit = []
    for u in [u1, u2]:
        pu = u[0].get()
        if pu == 'atm':
            alpUnit.append(101325)
        elif pu == 'bar':
            alpUnit.append(10**5)
        elif pu == 'Pa':
            alpUnit.append(1)
        elif pu == 'kPa':
            alpUnit.append(10**3)
        elif pu == 'MPa':
            alpUnit.append(10**6)
        pu = u[1].get()
        if pu == 'cubic meter':
            alpUnit.append(1)
        elif pu == 'cc' or pu == 'mL':
            alpUnit.append(10**-6)
        elif pu == 'L':
            alpUnit.append(10**-3)
    p, v, t = (initial[0].get(), initial[1].get(), initial[2].get())
    p = None if p == "" else float(p)*alpUnit[0]
    v = None if v == "" else float(v)*alpUnit[1]
    t = None if t == '' else (float(t)+273.15) if u1[2].get() == 'Celcius' else (
        float(t) if u1[2].get() == 'Kelvin' else (float(t)-32)*5/9+273.15)
    a = StateVariable(p, v, t)
    p, v, t = (final[0].get(), final[1].get(), final[2].get())
    p = None if p == "" else float(p)*alpUnit[2]
    v = None if v == "" else float(v)*alpUnit[3]
    t = None if t == '' else (float(t)+273.15) if u2[2].get() == 'Celcius' else (
        float(t) if u2[2].get() == 'Kelvin' else (float(t)-32)*5/9+273.15)
    k = k_entry.get()
    k = None if k == '' else float(k)
    try:
        y = y_entry.get()
        y = None if k == '' else float(y)
    except:
        y = None
    b = StateVariable(p, v, t)
    p = process.get()[0]
    clear()
    if p == 'I':
        pr = TherodynamicProcess.IsothermalReversible(a, b)
        k_entry.insert(0, pr.k)
    elif p == 'A':
        if y != None:
            pr = TherodynamicProcess.AdiabaticReversible(a, b, y, k)
        else:
            pr = TherodynamicProcess.AdiabaticReversible(a, b, k=k)
        k_entry.insert(0, pr.k)
        y_entry.insert(0, pr.gamma)
    elif p == 'P':
        pr = TherodynamicProcess.PolyIsotropicProcess(a, b)
        k_entry.insert(0, pr.k)
        y_entry.insert(0, pr.y)
    initial[0].insert(0, str(a.pressure/alpUnit[0]))
    initial[1].insert(0, str(a.volume/alpUnit[1]))
    t = a.temperature
    initial[2].insert(0, str((t-273.15) if u1[2].get() == 'Celcius' else (
        float(t) if u1[2].get() == 'Kelvin' else (t-273.15)*9/5+32)))
    t = b.temperature
    final[2].insert(0, str((t-273.15) if u2[2].get() == 'Celcius' else (float(t)
                    if u2[2].get() == 'Kelvin' else (t-273.15)*9/5+32)))
    final[0].insert(0, str(b.pressure/alpUnit[2]))
    final[1].insert(0, str(b.volume/alpUnit[3]))
    rs = tk.Label(r, text='Result')
    rs.grid(row=14, padx=5, pady=10)
    stats = pr.stats()
    ie = tk.Label(r, text='Internal Energy')
    ie.grid(row=15, column=1, padx=5, pady=5)
    wd = tk.Label(r, text='Work done by system')
    wd.grid(row=16, column=1, padx=5, pady=5)
    hr = tk.Label(r, text='Heat released by system')
    hr.grid(row=17, column=1, padx=5, pady=5)
    ier = tk.Label(r, text=f'{stats["u"]/1000:.3f} kJ')
    ier.grid(row=15, column=2, padx=5, pady=5)
    wdr = tk.Label(r, text=f'{-stats["w"]/1000:.3f} kJ')
    wdr.grid(row=16, column=2, padx=5, pady=5)
    hrd = tk.Label(r, text=f'{stats["q"]/1000:.3f} kJ')
    hrd.grid(row=17, column=2, padx=5, pady=5)
    globalDest.extend([rs, wd, hr, ie, ier, hrd, wdr])


p_units = ['Pa', 'kPa', 'MPa', 'atm', 'bar']
v_units = ['cubic meter', 'L', 'cc', 'mL']
t_units = ['Kelvin', 'Celcius', 'Faranheit']
units = [p_units, v_units, t_units]

tk.Label(r, text='Initial State').grid(row=0, column=0, padx=5, pady=5)
tk.Label(r, text='Final State').grid(row=4, column=0, padx=5, pady=5)
tk.Label(r, text='Process Type').grid(row=8, column=0, padx=5, pady=5)
for i in [[1, 2, 3], [5, 6, 7]]:
    tk.Label(r, text='Pressure').grid(row=i[0], column=1)
    tk.Label(r, text='Volume').grid(row=i[1], column=1)
    tk.Label(r, text='Temperature').grid(row=i[2], column=1)
u1 = [tk.StringVar(r), tk.StringVar(r), tk.StringVar(r)]
u2 = [tk.StringVar(r), tk.StringVar(r), tk.StringVar(r)]
for i in [u1, u2]:
    for k, x in enumerate(i):
        u1[k].set(units[k][0])
        u2[k].set(units[k][0])

for k, x in enumerate([1, 2, 3]):
    initial[k].grid(row=x, column=2, padx=15, pady=5)
    op = tk.OptionMenu(r, u1[k], *units[k])
    op.grid(row=x, column=3)
    op.config(width=12)
for k, x in enumerate([5, 6, 7]):
    final[k].grid(row=x, column=2, padx=15, pady=5)
    op = tk.OptionMenu(r, u2[k], *units[k])
    op.grid(row=x, column=3)
    op.config(width=12)

options = [
    "Isothermal Reversible Process",
    "Adiabatic Reversible Process",
    "Polyisotropic Process"
]
process = tk.StringVar(r)
process.trace("w", prcs)
process.set(options[0])
op = tk.OptionMenu(r, process, *options)
op.grid(row=9, column=2)
op.config(width=30)
tk.Button(r, text='Submit', command=submit).grid(
    row=13, column=0, padx=15, pady=10)
tk.Button(r, text='Clear', command=clear).grid(
    row=13, column=1, padx=15, pady=10)

r.mainloop()