#### Projekt 9: Simulation einer Temperaturregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
# Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen=10, start_temp=20, umg_temp=20, rpm=100, wandmaterial='Stahl', wandstaerke=0.005):
        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        temp_k = start_temp + 273.15
        self.volumen = volumen / 1000  # l -> m³
        self.spez_c = PropsSI("Cpmass", "T", temp_k, "P", 101325, "Water")
        self.dichte = PropsSI("Dmass", "T", temp_k, "P", 101325, "Water")
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)
        self.hoehe = 2 * self.radius
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe
        self.wandstaerke = wandstaerke
        self.d_ruehrer = (2 * self.radius) / 3
        self.umg_temp = umg_temp
        self.ist_temp = start_temp
        self.start_temp = start_temp
        self.rpm = rpm
        self.waerme_h = self.berechnung_h()
        self.lambda_wand = self.get_waermeleitfaehigkeit(wandmaterial)

    def berechnung_h(self):
        t_k = self.ist_temp + 273.15
        k = PropsSI("CONDUCTIVITY", "T", t_k, "P", 101325, "Water")
        mu = PropsSI("VISCOSITY", "T", t_k, "P", 101325, "Water")
        Re = ((self.rpm / 60) * (self.d_ruehrer ** 2) * self.dichte) / mu
        Pr = Prandtl(self.spez_c, mu, k)
        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:
            Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)
        elif Re >= 1e4 and Pr >= 0.6:
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
        else:
            Nu = 3.66
        return Nu * k / self.d_ruehrer

    def get_waermeleitfaehigkeit(self, material):
        material = str(material).lower()
        material_db = {
            'stahl': 21.0,
            'glas': 1.4,
        }
        return material_db.get(material, 16.0)

    def berechnung_waermeverlust(self):
        delta_t = self.ist_temp - self.umg_temp
        q_leitung = self.lambda_wand * (delta_t * self.flaeche / self.wandstaerke)
        q_konvektion = self.waerme_h * self.flaeche * delta_t
        q_verlust = q_leitung + q_konvektion
        return q_verlust

    def update_temperature(self, leistung, zeitintervall=1):
        masse = self.dichte * self.volumen
        waermeverlust = self.berechnung_waermeverlust()
        energie_netto = (leistung - waermeverlust) * zeitintervall
        temp_delta = energie_netto / (masse * self.spez_c)
        self.ist_temp += temp_delta
        return self.ist_temp

    def reset(self, start_temp=None):
        if start_temp is None:
            self.ist_temp = self.start_temp
        else:
            self.ist_temp = start_temp

## PID-Regler-Klasse
# Implementieren Sie einen PID-Regler, um eine Zieltemperatur zu halten

class PID:
    """
    PID-Regler zur Regelung der Temperatur eines Bioreaktors.
    """
    def __init__(self, dt=0.0, kp=0.0, ki=0.0, kd=0.0, output_min=0, output_max=1000, bioreaktor=None):
        self.dt = dt
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.output_min = output_min
        self.output_max = output_max
        self.fehler_vor = 0.0
        self.dt_vor = -1e-6
        self.integral = 0.0
        self.bioreaktor = bioreaktor  # Optional Referenz auf Bioreaktor

    def run(self, sollwert, istwert, offset=0.0):
        # Fehlerberechnung
        fehler = sollwert - istwert
        # P-Anteil
        p_anteil = self.kp * fehler
        # I-Anteil
        self.integral += fehler * self.dt
        i_anteil = self.ki * self.integral
        # D-Anteil
        delta_t = self.dt if self.dt - self.dt_vor > 0 else 1e-6
        d_anteil = self.kd * (fehler - self.fehler_vor) / delta_t
        # Ausgabe berechnen
        stellgroesse = offset + p_anteil + i_anteil + d_anteil
        # Begrenzung und Anti-Windup
        output = stellgroesse
        if output > self.output_max:
            output = self.output_max
            self.integral -= fehler * self.dt
        elif output < self.output_min:
            output = self.output_min
            self.integral -= fehler * self.dt
        self.dt_vor = self.dt
        self.fehler_vor = fehler
        return output

    def reset(self):
        self.fehler_vor = 0.0
        self.integral = 0.0

## Streamlit-Anwendung und Simulation

st.set_page_config(page_title="Bioreaktor Temperaturregelung", layout="centered")
st.title("Simulation einer Temperaturregelung von Bioreaktoren")

# Sidebar für Simulationseinstellungen
st.sidebar.title("Simulationseinstellungen")
st.sidebar.markdown("""
Hier können Sie die Parameter für die Simulation der Temperaturregelung eines Bioreaktors anpassen.
Die PID-Regelung wird verwendet, um die Temperatur des Reaktors auf einem Sollwert zu halten.
Die Auswirkungen der Regelparameter auf die Temperaturregelung können direkt beobachtet werden.
""")

with st.sidebar:
    st.header("Simulationsparameter")    
    volumen = st.slider("Reaktorvolumen (l)", 1, 10, 5)
    start_temp = st.slider("Starttemperatur (°C)", 0, 40, 20)
    umg_temp = st.slider("Umgebungstemperatur (°C)", 0, 40, 20)
    wandmaterial = st.selectbox("Wandmaterial", ["Stahl", "Glas"])                                                                         
    wandstärke = st.slider("Wandstärke (mm)", 1, 20, 5) / 1000                          
    soll_temp = st.slider("Solltemperatur (°C)", 20, 80, 37)
    simdauer = st.slider("Simulationsdauer (min)", 1, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 60, 5)
    rpm = st.slider("Rührerdrehzahl (1/min)", 1, 120, 60)
    st.header("PID-Parameter")
    kp = st.slider("Kp (proportional)", 0.0, 5.0, 1.0)
    ki = st.slider("Ki (integral)", 0.0, 1.0, 0.1)
    kd = st.slider("Kd (differential)", 0.00, 1.00, 0.10)
    st.markdown("---")
    st.caption("Tipp: Variiere die Reglerparameter und beobachte die Auswirkungen auf die Temperaturregelung.")

# Initialisierung der Bioreaktor- und PID-Regler-Objekte
reaktor_pid = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial, wandstaerke=wandstärke)
pid = PID(kp=kp, ki=ki, kd=kd, dt=dt, bioreaktor=reaktor_pid)
reaktor_ungeregelt = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial, wandstaerke=wandstärke)

# --- Simulation ---
n_steps = int(simdauer * 60 // dt)
zeiten = np.arange(0, n_steps * dt, dt) / 60
temps_pid, temps_offen, leistungen = [], [], []

reaktor_pid.reset(start_temp)
pid.reset()
reaktor_ungeregelt.reset(start_temp)

for t in range(n_steps):
    # PID-geregeltes System
    offset = reaktor_pid.berechnung_waermeverlust()
    leistung = pid.run(soll_temp, reaktor_pid.ist_temp, offset)
    temp_pid = reaktor_pid.update_temperature(leistung, dt)
    # Ungeregeltes System
    temp_offen = reaktor_ungeregelt.update_temperature(0, dt)
    # Speichern
    temps_pid.append(temp_pid)
    temps_offen.append(temp_offen)
    leistungen.append(leistung)

# --- Visualisierung ---
st.subheader("Temperaturverlauf im Bioreaktor")
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(zeiten, temps_pid, label="Mit PID-Regelung", color="tab:blue")
ax.plot(zeiten, temps_offen, "--", label="Ohne Regelung", color="tab:orange")
ax.axhline(soll_temp, color="grey", linestyle=":", label="Solltemperatur")
ax.set_xlabel("Zeit [min]")
ax.set_ylabel("Temperatur [°C]")
ax.legend()
ax.grid()
st.pyplot(fig)

st.subheader("Stellgröße (Heizleistung)")
fig2, ax2 = plt.subplots(figsize=(8, 2.5))
ax2.plot(zeiten, leistungen, color="tab:red")
ax2.set_xlabel("Zeit [min]")
ax2.set_ylabel("Leistung [W]")
ax2.grid()
st.pyplot(fig2)

with st.expander("Erklärung & Hinweise"):
    st.markdown("""
    - **Simulationseinstellungen**: Passen Sie die Parameter in der Sidebar an, um die Auswirkungen auf die Temperaturregelung zu beobachten.
    - **PID-Parameter**: Variieren Sie Kp, Ki und Kd, um das Verhalten des Reglers zu beeinflussen.
    - **Temperaturverlauf**: Der Graph zeigt den Temperaturverlauf im Bioreaktor mit und ohne PID-Regelung.
    - **Stellgröße**: Die Heizleistung des PID-Reglers wird im zweiten Graphen angezeigt.
    
    **Tipp**: Experimentieren Sie mit verschiedenen Reglerparametern, um ein besseres Verständnis für die Regelung von Bioreaktoren zu entwickeln.
    """)