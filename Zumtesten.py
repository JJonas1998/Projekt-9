##### Hochschule Weihenstephan-Triesdorf
##### Programmierung für Datenanalyse, Bildverarbeitung und Simulation 
##### Betreuerin: Prof. Dr. Kristina Eisen

#### Projekt 9: Simulation einer Temperaturregelung von Bioreaktoren 

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas
### Datum: 13.06.2025

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
class Bioreaktor:
    """Simuliert die Temperaturentwicklung in einem Bioreaktor."""

    # Materialdatenbank als Klassenattribut
    material_db = {
        'stahl': 46.0,  
        'edelstahl': 21.0,  
        'glas': 1.4,
        'kunststoff': 0.3,
        'aluminium': 230.0
    }

    def __init__(self, reaktor_vol=10, t_start=20, t_umgebung=20, drehz=100, wand_mat='stahl', wand_stk=5):
        self.t_start = t_start
        self.reaktor_vol = reaktor_vol / 1000                   # L → m³
        self.wand_stk = wand_stk / 1000                         # mm → m
        self.r_i = (self.reaktor_vol / (2 * np.pi)) ** (1 / 3)  # Innenradius in m 
        self.r_a = self.r_i + self.wand_stk                     # Außenradius
        self.h_i = 3 * self.r_i                                 # Innenhöhe
        self.h_a = self.h_i + 2 * self.wand_stk                 # Außenhöhe
        self.flaeche_i = 2 * np.pi * self.r_i**2 + 2 * np.pi * self.r_i * self.h_i
        self.flaeche_a = 2 * np.pi * self.r_a**2 + 2 * np.pi * self.r_a * self.h_a
        self.ruehrer_d = (2 * self.r_i) / 3                     # Rührerdurchmesser
        self.t_umgebung = t_umgebung
        self.t_ist = t_start
        self.drehz = drehz / 60                                 # 1/min → 1/s
        self.lambda_wand = self.berech_lambda(wand_mat)
        self.h_ext = 35                                         # Externer Wärmeübergangskoeffizient (W/(m²·K))
        self.max_leistung = 5000
        self.min_leistung = -2000
        self.update_stoffwerte(self.t_ist)
        self.h_int = self.berech_h_int()

    def update_stoffwerte(self, t):
        """Aktualisiert die thermophysikalischen Eigenschaften des Mediums (Wasser)."""
        if 0 < t < 100:  
            t_kelvin = t + 273.15   
            self.spez_c = PropsSI("Cpmass", "T", t_kelvin, "P", 101325, "Water")
            self.dichte = PropsSI("Dmass", "T", t_kelvin, "P", 101325, "Water")
            self.k = PropsSI("CONDUCTIVITY", "T", t_kelvin, "P", 101325, "Water")
            self.mu = PropsSI("VISCOSITY", "T", t_kelvin, "P", 101325, "Water")
        else:
            self.spez_c = 4186.0
            self.dichte = 997.0
            self.k = 0.606
            self.mu = 0.001002

    def berech_h_int(self):
        """Berechnet den internen Wärmeübergangskoeffizienten h."""
        try:
            Re = self.drehz * (self.ruehrer_d ** 2) * self.dichte / self.mu
            Pr = Prandtl(self.spez_c, self.mu, self.k)
            if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:
                Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)
            elif Re >= 1e4 and Pr >= 0.6:
                Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
            else:
                Nu = 3.66
            h = Nu * self.k / self.ruehrer_d
            return max(h, 150)
        except Exception:
            return 500

    @staticmethod
    def berech_lambda(wand_mat):
        return Bioreaktor.material_db.get(str(wand_mat).lower(), 46.0)

    def t_verlust(self):
        """Berechnet den Wärmeverlust des Reaktors (W)."""
        delta_t = self.t_ist - self.t_umgebung
        if abs(delta_t) < 0.01:
            return 0.0
        # Fläche für Wandmittelwert
        A_mittel = (self.flaeche_i + self.flaeche_a) / 2
        R_int = 1.0 / (self.h_int * self.flaeche_i)
        R_cond = self.wand_stk / (self.lambda_wand * A_mittel)
        R_ext = 1.0 / (self.h_ext * self.flaeche_a)
        R_gesamt = R_int + R_cond + R_ext
        q_verlust = delta_t / R_gesamt
        return q_verlust

    def update_temperatur(self, leistung, zeitintervall=1):
        """Aktualisiert die Reaktortemperatur."""
        leistung = np.clip(leistung, self.min_leistung, self.max_leistung)
        self.update_stoffwerte(self.t_ist)
        self.h_int = self.berech_h_int()
        masse = self.dichte * self.reaktor_vol
        energie_netto = (leistung - self.t_verlust()) * zeitintervall
        temp_delta = energie_netto / (masse * self.spez_c)
        self.t_ist += temp_delta
        self.t_ist = min(max(self.t_ist, 4), 100)
        return self.t_ist

    def reset(self, t_start=None):
        """Setzt die Ist-Temperatur auf die Starttemperatur zurück."""
        self.t_ist = self.t_start if t_start is None else t_start
        self.update_stoffwerte(self.t_ist)
        self.h_int = self.berech_h_int()

## PID-Regler-Klasse
class PID:
    """PID-Regler zur Regelung der Temperatur eines Bioreaktors."""
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, dt=1.0, output_min=-2000, output_max=5000):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.dt = dt
        self.output_min = output_min
        self.output_max = output_max
        self.reset()

    def run(self, sollwert, istwert):
        fehler = sollwert - istwert
        p_anteil = self.kp * fehler
        self.integral += fehler * self.dt
        self.integral = np.clip(self.integral, -self.integral_max, self.integral_max)
        i_anteil = self.ki * self.integral
        d_anteil = self.kd * (fehler - self.fehler_vor) / self.dt
        stellgroesse = p_anteil + i_anteil + d_anteil
        output = np.clip(stellgroesse, self.output_min, self.output_max)
        # Anti-Windup
        if self.ki > 0:
            if output != stellgroesse:
                self.integral -= (output - stellgroesse) / self.ki
        self.fehler_vor = fehler
        return output

    def reset(self):
        self.fehler_vor = 0.0
        self.integral = 0.0
        self.integral_max = abs(self.output_max - self.output_min) / max(self.ki, 1e-6)

# Streamlit-Anwendung und Simulation
st.set_page_config(page_title="Bioreaktor Temperaturregelung", layout="wide")  
st.title("Simulation einer Temperaturregelung von Bioreaktoren")  

# Sidebar für Einstellungen
st.sidebar.title("Simulationseinstellungen")  
st.sidebar.markdown("Passen Sie die Parameter an und beobachten Sie die Auswirkungen auf die Temperaturregelung.")

with st.sidebar:
    st.header("Reaktorparameter")    
    params = dict(
        reaktor_vol = st.slider("Reaktorvolumen (L)", 1, 100, 10),
        t_start = st.slider("Starttemperatur (°C)", 5, 40, 20),
        t_umgebung = st.slider("Umgebungstemperatur (°C)", 5, 40, 20),
        wand_mat = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Edelstahl", "Aluminium"]),
        wand_stk = st.slider("Wandstärke (mm)", 1, 50, 5),
        drehz = st.slider("Rührerdrehzahl (1/min)", 50, 180, 100)
    )

    st.header("Sollwert & Simulation")
    soll_temp = st.slider("Solltemperatur (°C)", 25, 80, 37)
    simdauer = st.slider("Simulationsdauer (min)", 10, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 30, 5)

    st.header("PID-Parameter")
    kp = st.slider("Kp (P-Anteil)", 0.0, 200.0, 50.0, 1.0, help = "Höher = schneller, aber instabiler")
    ki = st.slider("Ki (I-Anteil)", 0.0, 10.0, 1.0, 0.1, help = "Eliminiert bleibende Regelabweichung")
    kd = st.slider("Kd (D-Anteil)", 0.0, 50.0, 5.0, 0.5, help = "Dämpft Schwingungen")

# Tabs für Simulation und Analyse
tab1, tab2 = st.tabs(["📊 Simulation", "📋 Analyse"])

with tab1:
    col1, col2 = st.columns([2, 1])

    with col2:
        st.subheader("Reaktor-Eigenschaften")
        reaktor_info = Bioreaktor(**params)
        st.metric("Reaktorvolumen", f"{params['reaktor_vol']} L")
        st.metric("Reaktorradius", f"{reaktor_info.r_i*100:.1f} cm")
        st.metric("Reaktorhöhe", f"{reaktor_info.h_i*100:.1f} cm")
        st.metric("Wärmeübertr.-Fläche", f"{reaktor_info.flaeche_i:.2f} m²")
        st.metric("Wärmeleitfähigkeit", f"{reaktor_info.lambda_wand} W/(m·K)")

    with col1:
        try:
            # Initialisierung
            reaktor_pid = Bioreaktor(**params)
            reaktor_ungeregelt = Bioreaktor(**params)
            pid = PID(kp=kp, ki=ki, kd=kd, dt=dt, output_min=-2000, output_max=5000)
            n_steps = int(simdauer * 60 // dt)
            zeiten = np.arange(0, n_steps * dt, dt) / 60
            temps_pid, temps_offen, leistungen, sollwerte = [], [], [], []

            reaktor_pid.reset(params['t_start'])
            pid.reset()
            reaktor_ungeregelt.reset(params['t_start'])

            for i in range(n_steps):
                current_soll = soll_temp
                leistung = pid.run(current_soll, reaktor_pid.t_ist)
                temp_pid = reaktor_pid.update_temperatur(leistung, dt)
                temp_offen = reaktor_ungeregelt.update_temperatur(0, dt)
                temps_pid.append(temp_pid)
                temps_offen.append(temp_offen)
                leistungen.append(leistung)
                sollwerte.append(current_soll)

            # Visualisierung
            st.subheader("📈 Temperaturverlauf")
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
            ax1.plot(zeiten, temps_pid, label="Mit PID-Regelung", color="tab:blue", linewidth=2)
            ax1.plot(zeiten, temps_offen, "--", label="Ohne Regelung", color="tab:orange", linewidth=2)
            ax1.plot(zeiten, sollwerte, ":", label="Solltemperatur", color="tab:red", linewidth=2)
            ax1.set_ylabel("Temperatur [°C]")
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            ax1.set_title("Temperaturregelung im Bioreaktor")
            ax2.plot(zeiten, leistungen, color="tab:green", linewidth=2)
            ax2.set_xlabel("Zeit [min]")
            ax2.set_ylabel("Heizleistung [W]")
            ax2.grid(True, alpha=0.3)
            ax2.set_title("Stellgröße (Heizleistung)")
            plt.tight_layout()
            st.pyplot(fig)

            # Kennzahlen
            if temps_pid:
                steady_state_error = abs(temps_pid[-1] - soll_temp)
                overshoot = max(0, max(temps_pid) - soll_temp)
                col1_metric, col2_metric, col3_metric, col4_metric = st.columns(4)
                col1_metric.metric("Endtemperatur", f"{temps_pid[-1]:.1f}°C")
                col2_metric.metric("Regelabweichung", f"{steady_state_error:.2f}°C")
                col3_metric.metric("Überschwingen", f"{overshoot:.2f}°C")
                col4_metric.metric("Max. Heizleistung", f"{max(leistungen):.0f}W")
            # Hinweis bei Default-Werten
            if not (0 < temps_pid[-1] < 100):
                st.info("Achtung: Für Temperaturen außerhalb 0–100°C werden Default-Werte für Wasser verwendet.")

        except Exception as e:
            st.error(f"Fehler bei der Simulation: {str(e)}")
            st.info("Bitte überprüfen Sie die Eingabeparameter.")

with tab2:
    st.subheader("📊 Detaillierte Systemanalyse")
    try:
        st.write("**Regelgüte-Kennzahlen:**")
        # Einschwingzeit (Zeit bis 95% des Sollwerts erreicht)
        settling_time = "Nicht erreicht"
        for i, temp in enumerate(temps_pid):
            if abs(temp - soll_temp) <= 0.05 * soll_temp:
                settling_time = f"{zeiten[i]:.1f} min"
                break
        # Anstiegszeit (Zeit bis 90% des Sollwerts erreicht)
        rise_time = "Nicht erreicht"
        for i, temp in enumerate(temps_pid):
            if temp >= 0.9 * soll_temp:
                rise_time = f"{zeiten[i]:.1f} min"
                break
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Einschwingzeit (95%)", settling_time)
            st.metric("Anstiegszeit (90%)", rise_time)
            st.metric("Mittlere Abweichung", f"{np.mean(np.abs(np.array(temps_pid) - soll_temp)):.2f}°C")
        with col2:
            st.metric("Standardabweichung", f"{np.std(temps_pid):.2f}°C")
            st.metric("Energieverbrauch", f"{np.mean(np.maximum(leistungen, 0))*simdauer*60/1000:.1f}kJ")
            st.metric("Kühlenergie", f"{abs(np.mean(np.minimum(leistungen, 0)))*simdauer*60/1000:.1f}kJ")
        # Verlauf der Regelabweichung
        st.subheader("📈 Regelabweichung über Zeit")
        fig_error, ax_error = plt.subplots(figsize=(10, 4))
        fehler = np.array(temps_pid) - soll_temp
        ax_error.plot(zeiten, fehler, color="red", linewidth=2)
        ax_error.axhline(0, color="black", linestyle="--", alpha=0.5)
        ax_error.fill_between(zeiten, fehler, alpha=0.3, color="red")
        ax_error.set_xlabel("Zeit [min]")
        ax_error.set_ylabel("Regelabweichung [°C]")
        ax_error.set_title("Abweichung von der Solltemperatur")
        ax_error.grid(True, alpha=0.3)
        st.pyplot(fig_error)
    except Exception:
        st.info("Simulation wird geladen...")

st.sidebar.markdown("---")
st.sidebar.caption("🔄 **Live-Modus:** Ergebnisse werden automatisch aktualisiert!")