import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from fluids import Reynolds, Prandtl

class Bioreactor:
    """Bioreaktor mit dynamischem h aus Drehzahl."""
    def __init__(self,
                 volumen=100,              # L
                 spez_c=4180,              # J/(kg·K)
                 dichte=1000,              # kg/m³
                 starttemperatur=20,       # °C
                 rpm=100,                  # U/min
                 impeller_d=0.1):          # m
        self.volumen = volumen
        self.spez_c = spez_c
        self.dichte = dichte
        self.ambient_temp = starttemperatur
        self.temperature = starttemperatur
        self.biological_heat = 0
        self.rpm = rpm
        self.impeller_d = impeller_d
        # geschätzte Oberfläche (Zylinder-Äquivalent)
        self.surface_area = 2 * np.pi * (volumen/1000)**(2/3)

    def calc_heat_transfer_coeff(self):
        """h über Dittus–Boelter aus rpm und Impellerdurchmesser."""
        # Umrechnung auf SI
        T_K = self.temperature + 273.15
        # Thermo-Eigenschaften bei 1 bar
        cp  = PropsSI('Cpmass','T',T_K,'P',101325,'Water')
        k   = PropsSI('CONDUCTIVITY','T',T_K,'P',101325,'Water')
        mu  = PropsSI('VISCOSITY','T',T_K,'P',101325,'Water')
        rho = PropsSI('Dmass','T',T_K,'P',101325,'Water')
        # Umfangsgeschwindigkeit am Propellerrand
        v = (self.rpm/60) * np.pi * self.impeller_d
        Re = Reynolds(rho, v, self.impeller_d, mu)
        Pr = Prandtl(cp, mu, k)
        Nu = (0.037 * Re ** 0.8 * Pr) / (1 + 2.443 * Re ** (- 0.1) * (Pr ** (2 / 3) - 1))
        h = Nu * k / self.impeller_d
        return h

    def update_temperature(self, heating_power, dt=1.0):
        """Bilanz mit dyn. h."""
        mass = self.volumen * self.dichte / 1000
        h = self.calc_heat_transfer_coeff()
        # Wärmeverlust
        q_loss = h * self.surface_area * (self.temperature - self.ambient_temp)
        # biologische Wärme
        q_bio  = self.biological_heat * mass
        # Netto
        q_net  = heating_power + q_bio - q_loss
        dT     = q_net * dt / (mass * self.spez_c)
        self.temperature += dT
        return self.temperature

class PIDController:
    def __init__(self, kp=1000, ki=10, kd=50):
        self.kp, self.ki, self.kd = kp, ki, kd
        self.setpoint = 37.0
        self.integral = 0
        self.last_error = 0
        self.dt = 1.0

    def update(self, current_temp):
        error = self.setpoint - current_temp
        p = self.kp * error
        self.integral += error * self.dt
        i = self.ki * self.integral
        d = self.kd * (error - self.last_error)/self.dt
        self.last_error = error
        u = p + i + d
        return max(0, min(5000, u)), p, i, d

def simulate_system(reactor, controller, duration, with_disturbance=False):
    times, temps, setpts, powers = [], [], [], []
    ps, is_, ds = [], [], []
    for t in range(duration):
        if with_disturbance and t==duration//2:
            reactor.biological_heat = 20
        u, p, i, d = controller.update(reactor.temperature)
        reactor.update_temperature(u)
        times.append(t); temps.append(reactor.temperature)
        setpts.append(controller.setpoint); powers.append(u)
        ps.append(p); is_.append(i); ds.append(d)
    return pd.DataFrame({
        'Zeit': times,
        'Temperatur': temps,
        'Sollwert': setpts,
        'Heizleistung': powers,
        'P': ps, 'I': is_, 'D': ds
    })

def main():
    st.title("🧪 Bioreaktor-PID mit drehzahlabhängigem h")
    st.sidebar.header("Parameter")
    # Reaktor
    vol = st.sidebar.slider("Volumen (L)", 50, 500, 100)
    Tamb = st.sidebar.slider("Umgebungstemperatur (°C)", 15, 30, 20)
    # Rührer
    rpm = st.sidebar.slider("Rührer-Drehzahl (U/min)", 0, 500, 100)
    d_imp = st.sidebar.slider("Impeller-D (cm)", 5, 30, 10)/100
    # PID
    kp = st.sidebar.slider("Kp", 100, 2000, 1000)
    ki = st.sidebar.slider("Ki", 1, 50, 10)
    kd = st.sidebar.slider("Kd", 10, 100, 50)
    sp = st.sidebar.slider("Soll-T (°C)", 30, 45, 37)
    # Simulation
    dur = st.sidebar.slider("Dauer (min)", 30, 180, 60)
    disturb = st.sidebar.checkbox("Störung", False)

    if st.button("🚀 Simulieren"):
        reactor   = Bioreactor(vol, 4180, 1000, Tamb, rpm, d_imp)
        controller= PIDController(kp, ki, kd)
        controller.setpoint = sp
        df = simulate_system(reactor, controller, dur, disturb)

        # Grafik wie gehabt ...
        fig, axes = plt.subplots(3,1,figsize=(10,8))
        axes[0].plot(df.Zeit, df.Temperatur,'b',label='T')
        axes[0].plot(df.Zeit, df.Sollwert,'r--',label='SP')
        axes[0].legend(); axes[0].set_ylabel("°C")
        axes[1].plot(df.Zeit, df.Heizleistung,'g'); axes[1].set_ylabel("W")
        axes[2].plot(df.Zeit, df.P,'r',alpha=0.7,label='P')
        axes[2].plot(df.Zeit, df.I,'y',alpha=0.7,label='I')
        axes[2].plot(df.Zeit, df.D,'c',alpha=0.7,label='D')
        axes[2].legend(); axes[2].set_ylabel("W")
        axes[2].set_xlabel("Zeit (min)")
        plt.tight_layout()
        st.pyplot(fig)

        st.dataframe(df.head())

if __name__=="__main__":
    main()
