import streamlit as st
import matplotlib.pyplot as plt

# --- Bioreaktor-Klasse ---
class Bioreaktor:
    def __init__(self, start_temp=25.0, umgeb_temp=20.0, volumen=1.0, k_wärme=0.1):
        self.temp = start_temp
        self.umgeb_temp = umgeb_temp
        self.volumen = volumen
        self.k_wärme = k_wärme

    def update(self, leistung, dt):
        dT = (leistung - self.k_wärme * (self.temp - self.umgeb_temp)) * dt / self.volumen
        self.temp += dT
        return self.temp

# --- Verbesserte PID-Regler-Klasse ---
class PIDRegler:
    def __init__(self, dt, output_max, output_min, kp, kd, ki):
        self.dt = dt
        self.output_max = output_max
        self.output_min = output_min
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.last_error = 0.0
        self.integral = 0.0

    def run(self, setpoint, actual):
        error = setpoint - actual
        P = self.kp * error
        self.integral += error * self.dt
        I = self.ki * self.integral
        D = self.kd * (error - self.last_error) / self.dt
        output = P + I + D
        output = max(self.output_min, min(self.output_max, output))
        self.last_error = error
        return output

# --- Simulationsfunktion ---
def simuliere_bioreaktor(zeit=100, dt=0.5, mit_regler=True, pid_param=(2.0,0.2,0.05), ziel_temp=37.0, start_temp=25.0, umgeb_temp=20.0, volumen=2.0, k_wärme=0.1):
    reaktor = Bioreaktor(start_temp=start_temp, umgeb_temp=umgeb_temp, volumen=volumen, k_wärme=k_wärme)
    pid = PIDRegler(dt=dt, output_max=10, output_min=0, kp=pid_param[0], ki=pid_param[1], kd=pid_param[2])
    temps, zeiten, leistungen = [], [], []
    for t in range(int(zeit / dt)):
        ist_temp = reaktor.temp
        if mit_regler:
            leistung = pid.run(ziel_temp, ist_temp)
        else:
            leistung = 5 if ist_temp < ziel_temp else 0
        reaktor.update(leistung, dt)
        temps.append(reaktor.temp)
        zeiten.append(t * dt)
        leistungen.append(leistung)
    return zeiten, temps, leistungen

# --- Streamlit-Interface ---
st.title("PID-Regelung für einen Bioreaktor")
st.markdown("Simuliere die Temperaturregelung mit einem verbesserten PID-Regler.")

# Seitenleiste für Parameter
st.sidebar.header("Simulationsparameter")
zeit = st.sidebar.slider("Simulationszeit [min]", 20, 200, 100, 5)
dt = st.sidebar.slider("Zeitschritt dt [min]", 1, 20, 5, 1) / 10
ziel_temp = st.sidebar.slider("Solltemperatur [°C]", 25, 50, 37, 1)
start_temp = st.sidebar.slider("Starttemperatur [°C]", 10, 45, 25, 1)
umgeb_temp = st.sidebar.slider("Umgebungstemperatur [°C]", 0, 40, 20, 1)
volumen = st.sidebar.slider("Reaktorvolumen [L]", 1, 10, 2, 1)
k_wärme = st.sidebar.slider("Wärmeverlust-Koeffizient", 1, 20, 10, 1) / 100

st.sidebar.header("PID-Regler Parameter")
Kp = st.sidebar.slider("Kp (Proportional)", 0.1, 10.0, 2.0, 0.1)
Ki = st.sidebar.slider("Ki (Integral)", 0.0, 1.0, 0.2, 0.01)
Kd = st.sidebar.slider("Kd (Differential)", 0.0, 1.0, 0.05, 0.01)

# Simulation und Plot
ez1, t1, l1 = simuliere_bioreaktor(zeit, dt, mit_regler=False, ziel_temp=ziel_temp, start_temp=start_temp, umgeb_temp=umgeb_temp, volumen=volumen, k_wärme=k_wärme, pid_param=(Kp, Ki, Kd))
ez2, t2, l2 = simuliere_bioreaktor(zeit, dt, mit_regler=True, ziel_temp=ziel_temp, start_temp=start_temp, umgeb_temp=umgeb_temp, volumen=volumen, k_wärme=k_wärme, pid_param=(Kp, Ki, Kd))

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(ez1, t1, label="Ohne Regelung", color="red", linestyle="--")
ax.plot(ez2, t2, label="Mit PID-Regelung", color="green")
ax.axhline(ziel_temp, color="blue", linestyle=":", label="Sollwert")
ax.set_xlabel("Zeit [min]")
ax.set_ylabel("Temperatur [°C]")
ax.set_title("Temperaturverlauf im Bioreaktor")
ax.legend()
ax.grid()
st.pyplot(fig)

st.markdown("""
**Hinweise:**  
- Die PID-Parameter beeinflussen, wie schnell und stabil die Regelung reagiert.  
- Teste verschiedene Szenarien – z.B. niedrige/heißere Start- oder Umgebungstemperatur, anderes Volumen, verschiedene Kp/Ki/Kd.  
- Der Plot zeigt den Vergleich zwischen einer ungeregelten ("nur Heizung an/aus") und einer PID-geregelten Temperatursteuerung.
""")
