import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# ------------------ Bioreaktor Klasse ------------------
class Bioreaktor:
    def __init__(self, volumen=10, start_temp=20, umg_temp=20, wandmaterial='stahl', wandstaerke=0.005):
        self.volumen = volumen / 1000  # l -> m³
        self.c = 4180                  # J/(kg*K) (angenommen Wasser)
        self.dichte = 1000             # kg/m³
        self.umg_temp = umg_temp
        self.ist_temp = start_temp
        self.wandmaterial = wandmaterial
        self.wandstaerke = wandstaerke
        self.radius = (self.volumen / (2 * np.pi)) ** (1 / 3)
        self.hoehe = 2 * self.radius
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe
        self.lambda_wand = {"stahl": 21.0, "glas": 1.4, "aluminium": 230.0}.get(wandmaterial.lower(), 21.0)
        self.h = 10.0  # W/(m²*K), typischer Wert für Wasser

    def berechnung_waermeverlust(self):
        delta_t = self.ist_temp - self.umg_temp
        q_leitung = self.lambda_wand * (delta_t * self.flaeche / self.wandstaerke)
        q_konvektion = self.h * self.flaeche * delta_t
        return q_leitung + q_konvektion

    def update_temperature(self, leistung, dt=1):
        masse = self.dichte * self.volumen
        q_verlust = self.berechnung_waermeverlust()
        energie = (leistung - q_verlust) * dt
        delta_temp = energie / (masse * self.c)
        self.ist_temp += delta_temp
        return self.ist_temp

    def reset(self, temp):
        self.ist_temp = temp

# ------------------ PID-Regler Klasse ------------------
class PID:
    def __init__(self, kp=500, ki=0.1, kd=200, dt=1, output_min=0, output_max=2000):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.dt = dt
        self.output_min = output_min
        self.output_max = output_max
        self.integral = 0.0
        self.fehler_vor = 0.0

    def run(self, soll, ist, offset=0):
        fehler = soll - ist
        self.integral += fehler * self.dt
        diff = (fehler - self.fehler_vor) / self.dt
        p = self.kp * fehler
        i = self.ki * self.integral
        d = self.kd * diff
        leistung = offset + p + i + d
        leistung = np.clip(leistung, self.output_min, self.output_max)
        if leistung in [self.output_min, self.output_max]:
            self.integral -= fehler * self.dt  # Anti-Windup
        self.fehler_vor = fehler
        return leistung

    def reset(self):
        self.integral = 0.0
        self.fehler_vor = 0.0

# ------------------ Streamlit-UI ------------------
st.set_page_config(page_title="Bioreaktor Temperaturregelung", layout="centered")
st.title("Simulation: Temperaturregelung im Bioreaktor")

# --- Sidebar: Parameter ---
with st.sidebar:
    st.header("Simulationsparameter")
    volumen = st.slider("Reaktorvolumen (l)", 1, 10, 5)
    start_temp = st.slider("Starttemperatur (°C)", 5, 35, 20)
    umg_temp = st.slider("Umgebungstemperatur (°C)", 5, 40, 25)
    wandmaterial = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Aluminium"])
    soll_temp = st.slider("Solltemperatur (°C)", 25, 45, 37)
    simdauer = st.slider("Simulationsdauer (min)", 1, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 60, 5)
    st.header("PID-Parameter")
    kp = st.slider("Kp (proportional)", 0.0, 2000.0, 500.0)
    ki = st.slider("Ki (integral)", 0.0, 5.0, 0.1)
    kd = st.slider("Kd (differential)", 0.0, 1000.0, 200.0)
    st.markdown("---")
    st.caption("Tipp: Variiere die Reglerparameter und beobachte die Auswirkungen auf die Temperaturregelung.")

# --- Initialisierung ---
reaktor_pid = Bioreaktor(volumen, start_temp, umg_temp, wandmaterial)
pid = PID(kp=kp, ki=ki, kd=kd, dt=dt)
reaktor_offen = Bioreaktor(volumen, start_temp, umg_temp, wandmaterial)

# --- Simulation ---
n_steps = int(simdauer * 60 // dt)
zeiten = np.arange(0, n_steps * dt, dt) / 60  # Minuten
temps_pid, temps_offen, leistungen = [], [], []

reaktor_pid.reset(start_temp)
pid.reset()
reaktor_offen.reset(start_temp)

for t in range(n_steps):
    # PID-geregeltes System
    offset = reaktor_pid.berechnung_waermeverlust()
    leistung = pid.run(soll_temp, reaktor_pid.ist_temp, offset)
    temp_pid = reaktor_pid.update_temperature(leistung, dt)
    # Ungeregeltes System
    temp_offen = reaktor_offen.update_temperature(0, dt)
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
    **Was passiert hier?**  
    Der PID-Regler sorgt dafür, dass die Bioreaktortemperatur möglichst konstant auf dem Sollwert bleibt.  
    Ohne Regelung kühlt der Reaktor durch Wärmeverluste ab.

    **PID-Parameter:**  
    - `Kp` bestimmt, wie stark auf die Abweichung reagiert wird.
    - `Ki` eliminiert bleibende Regelabweichungen (z.B. Offset).
    - `Kd` dämpft schnelle Änderungen (überschwingen).
    """)