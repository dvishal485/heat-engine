import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad as integrate

n = 1
R = 8.314
gamma = 5/3

Cv = R/(gamma-1)
Cp = Cv + R


class StateVariable:
    def __init__(self, pressure: float = None, volume: float = None, temperature: float = None) -> None:
        if pressure == None and volume != None and temperature != None:
            pressure = n * R * temperature / volume
        if pressure != None and volume == None and temperature != None:
            volume = n * R * temperature / pressure
        if pressure != None and volume != None and temperature == None:
            temperature = pressure * volume / (n * R)
        self.pressure = pressure
        self.volume = volume
        self.temperature = temperature


class TherodynamicProcess:
    def assign(variable: StateVariable, value: StateVariable, inplace: bool = True):
        '''Rewrites a state variable according to given information'''
        T = variable.temperature
        P = variable.pressure
        V = variable.volume
        if T == None:
            T = value.temperature
        if V == None:
            V = value.volume
        if P == None:
            P = value.pressure
        if inplace:
            variable.volume = V
            variable.pressure = P
            variable.temperature = T
            return variable
        else:
            return StateVariable(P, V, T)

    class IsothermalReversible:
        def __init__(self, a: StateVariable, b: StateVariable):
            TherodynamicProcess.assign(
                a, StateVariable(temperature=b.temperature))
            TherodynamicProcess.assign(
                b, StateVariable(temperature=a.temperature))
            k = None

            vars = [a, b]
            # PV = constant (say k) = nRT
            for i in vars:
                try:
                    k = i.pressure * i.volume
                except:
                    try:
                        k = n * R * i.temperature
                    except:
                        None
            if k == None:
                raise ValueError(
                    "Not enough information provided to process isothermal process")
            else:
                k: float
            for i in vars:
                if i.temperature == None:
                    try:
                        i.temperature = k/(n*R)
                    except:
                        None
                if i.pressure == None:
                    try:
                        TherodynamicProcess.assign(i, StateVariable(
                            pressure=k/i.volume))
                    except:
                        None
                if i.volume == None:
                    try:
                        TherodynamicProcess.assign(i, StateVariable(
                            volume=k/i.pressure))
                    except:
                        None

            TherodynamicProcess.assign(
                a, StateVariable(temperature=b.temperature))
            TherodynamicProcess.assign(
                b, StateVariable(temperature=a.temperature))

            if a.temperature != b.temperature:
                raise ArithmeticError(
                    'Temperature is not same at start and end point of isothermal process')
            self.a: StateVariable = a
            self.b: StateVariable = b
            self.k: float = k

        def getStateViaVolume(self, volume: float):
            '''Returns state for isothermal process for a given volume'''
            return StateVariable(pressure=self.k/volume, volume=volume)

        def getStateViaPressure(self, pressure: float):
            '''Returns state for isothermal process for a given pressure'''
            return StateVariable(volume=self.k/pressure, pressure=pressure)

        def stats(self):
            '''Returns certains quantities related to the process'''
            u = TherodynamicProcess.changeInInternalEnergy(self.a, self.b)
            w = float(-n*R*self.a.temperature *
                      np.log(self.b.volume/self.a.volume))
            return {
                'u': u,
                'w': w,
                'q': u - w
            }

        def coordinates(self):
            ll = np.min([self.a.volume, self.b.volume])
            l = np.max([self.a.volume, self.b.volume])
            volumes = np.linspace(ll, l, 100)
            p = []
            v = []
            t = []
            for volume in volumes:
                state = self.getStateViaVolume(volume)
                if state.pressure != None:
                    p.append(state.pressure)
                    v.append(state.volume)
                    t.append(state.temperature)
            return {'p': p, 'v': v, 't': t}

    def changeInInternalEnergy(a: StateVariable, b: StateVariable) -> float:
        '''Returns change in internal energy in process `a` to `b`'''
        return n * Cv * (b.temperature - a.temperature)

    class AdiabaticReversible:
        def __init__(self, a: StateVariable, b: StateVariable, gamma: float = gamma, k: float = None):
            '''	
            PV^&#120574; = constant
            '''
            vars = [a, b]
            # PV^gamma = constant (say k)
            for i in vars:
                try:
                    k = i.pressure * (i.volume ** gamma)
                except:
                    try:
                        k = n * R * i.temperature * (i.volume ** (gamma - 1))
                    except:
                        try:
                            k = i.pressure * \
                                ((n * R * i.temperature)/(i.pressure))**gamma
                        except:
                            None
            if k == None:
                raise ValueError(
                    "Not enough information provided to process adiabatic process")
            else:
                k: float
            for i in vars:
                if i.pressure == None:
                    try:
                        TherodynamicProcess.assign(
                            i, StateVariable(pressure=k/(i.volume**gamma)))
                    except:
                        TherodynamicProcess.assign(
                            i, StateVariable(pressure=(k/(n*R*i.temperature)**gamma)**(1/(1-gamma))))
                if i.volume == None:
                    try:
                        TherodynamicProcess.assign(i, StateVariable(
                            volume=(k/i.pressure)**(1/gamma)))
                    except:
                        TherodynamicProcess.assign(i, StateVariable(
                            volume=(k/n*R*i.temperature)**(1/(gamma-1))))
                if i.temperature == None:
                    try:
                        TherodynamicProcess.assign(i, StateVariable(
                            temperature=k/(n*R*(i.volume**(gamma-1)))))
                    except:
                        TherodynamicProcess.assign(i, StateVariable(
                            temperature=(((k/i.pressure)**(1/gamma))*i.pressure/(n*R))))
            self.a: StateVariable = a
            self.b: StateVariable = b
            self.gamma: float = gamma
            self.k: float = k

        def getStateViaVolume(self, volume: float):
            '''Returns state for adiabatic process for a given volume'''
            return StateVariable(pressure=self.k/(volume**self.gamma), volume=volume)

        def getStateViaPressure(self, pressure: float):
            '''Returns state for adiabatic process for a given pressure'''
            return StateVariable(pressure=pressure, volume=(self.k/pressure)**(1/self.gamma))

        def stats(self):
            u = TherodynamicProcess.changeInInternalEnergy(self.a, self.b)
            try:
                if self.a.temperature != None and self.b.temperature != None:
                    w = n * R * (self.b.temperature -
                                 self.a.temperature)/(self.gamma-1)
                else:
                    w = (self.b.pressure * self.b.volume -
                         self.a.pressure * self.a.volume)/(self.gamma-1)
            except:
                try:
                    w = (n * R * (self.b.temperature) -
                         self.a.pressure * self.a.volume)/(self.gamma-1)
                except:
                    try:
                        w = (self.b.pressure * self.b.volume - n *
                             R * (self.a.temperature))/(self.gamma-1)
                    except:
                        raise ValueError(
                            'Failed to process work done in the process')
            w: float
            return {
                'u': u,
                'w': w,
                'q': u-w
            }

        def coordinates(self):
            ll = np.min([self.a.volume, self.b.volume])
            l = np.max([self.a.volume, self.b.volume])
            volumes = np.linspace(ll, l, 100)
            p = []
            v = []
            t = []
            for volume in volumes:
                state = self.getStateViaVolume(volume)
                if state.pressure != None:
                    p.append(state.pressure)
                    v.append(state.volume)
                    t.append(state.temperature)
            return {'p': p, 'v': v, 't': t}

    class PolyIsotropicProcess:
        def __init__(self, a: StateVariable, b: StateVariable, y: float, k: float = None):
            '''
            Process satisfying `PV^y = constant`
            ---
            Parameters
            - `a` : Initial `StateVariable`
            - `b` : Final `StateVariable`
            - `y` : Value of power of Volume in the expression `PV^y = constant`
            '''
            process = TherodynamicProcess.AdiabaticReversible(
                a, b, gamma=y, k=k)
            self.process = process
            self.a = a
            self.b = b
            self.k = process.k
            self.y = y

        def getStateViaVolume(self, volume: float):
            return self.process.getStateViaVolume(volume)

        def getStateViaPressure(self, pressure: float):
            return self.process.getStateViaPressure(pressure)

        def stats(self):
            return self.process.stats()

        def coordinates(self):
            ll = np.min([self.a.volume, self.b.volume])
            l = np.max([self.a.volume, self.b.volume])
            volumes = np.linspace(ll, l, 100)
            p = []
            v = []
            t = []
            for volume in volumes:
                state = self.getStateViaVolume(volume)
                if state.pressure != None:
                    p.append(state.pressure)
                    v.append(state.volume)
                    t.append(state.temperature)
            return {'p': p, 'v': v, 't': t}

    class DefineRulePV:
        def __init__(self, relation) -> None:
            '''
            Create a custom P-V relation
            ---
            Parameters
            - `relation` : A function such that on calling `relation(volume:float)`
            the function should return value of pressure as a `float`
            '''
            self.relation = relation

        def getStateViaVolume(self, volume: float):
            return StateVariable(pressure=self.relation(volume), volume=volume)

        def stats(self, a: StateVariable, b: StateVariable):
            a.pressure = self.relation(a.volume)
            b.pressure = self.relation(b.volume)
            alpha = TherodynamicProcess.assign(StateVariable(), a, False)
            beta = TherodynamicProcess.assign(StateVariable(), b, False)
            TherodynamicProcess.assign(a, alpha)
            TherodynamicProcess.assign(b, beta)
            w = integrate(self.relation, a.volume, b.volume)[0]
            u = TherodynamicProcess.changeInInternalEnergy(a, b)
            return {
                'u': u, 'w': w, 'q': u - w
            }

        def coordinates(self, a: StateVariable, b: StateVariable):
            ll = np.min([a.volume, b.volume])
            l = np.max([a.volume, b.volume])
            TherodynamicProcess.assign(a, StateVariable(
                volume=a.volume, pressure=self.relation(a.volume)))
            TherodynamicProcess.assign(b, StateVariable(
                volume=b.volume, pressure=self.relation(b.volume)))
            volumes = np.linspace(ll, l, 100)
            p = []
            v = []
            t = []
            for volume in volumes:
                state = self.getStateViaVolume(volume)
                if state.pressure != None:
                    p.append(state.pressure)
                    v.append(state.volume)
                    t.append(state.temperature)
            return {'p': p, 'v': v, 't': t}

    def carnotEnginePlot(a: StateVariable, b: StateVariable, c: StateVariable, d: StateVariable, plot: bool = True, save_name: str = None):
        '''
        Plots the heat enigine's P-V curve
        ---
        Parameters

        - `a`, `b`, `c`, `d` : Class `StateVariables` containing end-points
        of the heat-engine's thermodynamic processes
        ---
        Process

        - `a` to `b` and `c` to `d` will be treated as Isothermal process
        - `b` to `c` and `d` to `a` will be treated as Adiabatic process
        '''
        trial = 0
        while(trial < 5):
            trial = trial + 1
            try:
                process1 = TherodynamicProcess.IsothermalReversible(a, b)
            except:
                None
            try:
                process2 = TherodynamicProcess.AdiabaticReversible(b, c)
            except:
                None
            try:
                process3 = TherodynamicProcess.IsothermalReversible(c, d)
            except:
                None
            try:
                process4 = TherodynamicProcess.AdiabaticReversible(d, a)
            except:
                None

        isothermalProcesses: list[TherodynamicProcess.IsothermalReversible] = [
            process1, process3]
        adiabaticProcesses: list[TherodynamicProcess.AdiabaticReversible] = [
            process2, process4]
        energy = []
        work = []
        heat = []
        if plot:
            plt.style.use('seaborn')
            fig, ax = plt.subplots(figsize=(13, 5), ncols=2)
        for process in isothermalProcesses:
            work.append(process.stats()['w'])
            energy.append(process.stats()['u'])
            heat.append(process.stats()['q'])
            pvt = process.coordinates()
            p1 = pvt['p']
            t1 = pvt['t']
            v1 = pvt['v']
            if plot:
                ax[0].plot(v1, p1, label="Isothermal Process")
                ax[1].plot(t1, v1, label="Isothermal Process")

        for process in adiabaticProcesses:
            work.append(process.stats()['w'])
            energy.append(process.stats()['u'])
            heat.append(process.stats()['q'])
            pvt = process.coordinates()
            p2 = pvt['p']
            t2 = pvt['t']
            v2 = pvt['v']
            if plot:
                ax[0].plot(v2, p2, label="Adiabatic Process")
                ax[1].plot(t2, v2, label="Adiabatic Process")
        w, q, u = np.array(work).sum(), np.array(
            heat).sum(), np.array(energy).sum()
        heat = np.array(heat)
        # work done by gas / heat absorbed by gas (effective)
        efficiency = -w/np.array((heat[heat > 0])).sum()
        if plot:
            ax[0].set_ylabel("Pressure", fontsize=13)
            ax[0].set_xlabel("Volume", fontsize=13)
            ax[1].set_ylabel("Volume", fontsize=13)
            ax[1].set_xlabel("Temperature", fontsize=13)
            vars = [a, b, c, d]
            names = ['A', 'B', 'C', 'D']
            vg = []
            pg = []
            tg = []
            for n, i in enumerate(vars):
                ax[0].annotate(
                    f"{names[n]} ({i.pressure:.1f},\n    {i.volume:.1f}, {i.temperature:.1f})", (i.volume, i.pressure + (ax[0].get_ylim()[1]-ax[0].get_ylim()[0])/100*5), backgroundcolor="w")
                ax[1].annotate(
                    f"{names[n]} ({i.pressure:.1f},\n    {i.volume:.1f},{i.temperature:.1f})", (i.temperature, i.volume+(ax[1].get_ylim()[1]-ax[1].get_ylim()[0])/100*5),  backgroundcolor="w")
                vg.append(i.volume)
                pg.append(i.pressure)
                tg.append(i.temperature)
            ax[0].set_ylim(ax[0].get_ylim()[0], ax[0].get_ylim()[1]*1.1)
            ax[1].set_ylim(ax[1].get_ylim()[0], ax[1].get_ylim()[1]*1.1)
            ax[0].set_xlim(ax[0].get_xlim()[0], ax[0].get_xlim()[1]*1.08)
            ax[1].set_xlim(ax[1].get_xlim()[0], ax[1].get_xlim()[1]*1.08)
            ax[0].scatter(vg, pg)
            ax[0].legend()
            ax[1].scatter(tg, vg)
            fig.suptitle(
                f'\u0394U = {u:.2f} J, w = {w:.2f} J, q = {q:.2f} J\n\u03B7 = {efficiency*100:.2f}%')
            ax[1].legend()
            fig.tight_layout()
        if save_name != None and plot:
            fig.savefig(f'{save_name}.png',
                        bbox_inches='tight',
                        transparent=True,)
        return {'u': u, 'w': w, 'q': q, 'e': efficiency}

    def plot(pvt: dict):
        '''Plots a P-V and V-T curve for given set of P, V, T Coordinates'''
        plt.style.use('seaborn')
        fig, ax = plt.subplots(figsize=(15, 5), ncols=2)
        v = pvt['v']
        p = pvt['p']
        t = pvt['t']
        l = len(v)-1
        ax[0].plot(v, p)
        ax[1].plot(t, v)
        ax[0].set_xlabel('Volume')
        ax[0].set_ylabel('Pressure')
        ax[1].set_xlabel('Temperature')
        ax[1].set_ylabel('Volume')
        ax[0].annotate(f"A", (v[0], p[0] + (ax[0].get_ylim()[1] -
                       ax[0].get_ylim()[0])/100*5), backgroundcolor="w")
        ax[1].annotate(f"A", (t[0], v[0]+(ax[1].get_ylim()[1] -
                       ax[1].get_ylim()[0])/100*5),  backgroundcolor="w")
        ax[0].annotate(f"B", (v[l], p[l] + (ax[0].get_ylim()[1] -
                       ax[0].get_ylim()[0])/100*5), backgroundcolor="w")
        ax[1].annotate(f"B", (t[l], v[l]+(ax[1].get_ylim()[1] -
                       ax[1].get_ylim()[0])/100*5),  backgroundcolor="w")
        ax[0].set_ylim(ax[0].get_ylim()[0], ax[0].get_ylim()[1]*1.05)
        ax[1].set_ylim(ax[1].get_ylim()[0], ax[1].get_ylim()[1]*1.05)
        ax[0].set_xlim(ax[0].get_xlim()[0], ax[0].get_xlim()[1]*1.05)
        ax[1].set_xlim(ax[1].get_xlim()[0], ax[1].get_xlim()[1]*1.05)
        ax[0].scatter([v[0], v[l]], [p[0], p[l]], c='red')
        ax[1].scatter([t[0], t[l]], [v[0], v[l]], c='red')
