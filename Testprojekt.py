##### Hochschule Weihenstephan-Triesdorf
##### Programmierung für Datenanalyse, Bildverarbeitung und Simulation 
##### Betreuerin: Prof. Dr. Kristina Eisen

#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren 

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas
### Datum: 13.06.2025

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
class Bioreaktor:
    """Simuliert die Temperaturentwicklung in einem Bioreaktor.

    Die Temperatur wird beeinflusst durch Heizen, Kühlen und Umwelteinflüsse.
    """
    def __init__(
            self, 
            reaktor_vol = 10, 
            t_start = 20, 
            t_umgebung = 20, 
            drehz = 100, 
            wand_mat = 'stahl', 
            wand_stk = 5):
        """Initialisiert den Bioreaktor mit Volumen, Temperatur und Materialeigenschaften."""

        # Speichert die Starttemperatur für einen späteren Reset
        self.t_start = t_start
        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        self.update_stoffwerte(t_start)
        # Betriebsbedingungen
        self.t_umgebung = t_umgebung    # Umgebungstemperatur in °C 
        self.t_ist = t_start            # Aktuelle Temperatur in °C
        self.drehz = drehz / 60         # Rührerdrehzahl in 1/s (Umrechnung von 1/min in 1/s)
        # Geometrie des Bioreaktors 
        self.geometrie_daten(reaktor_vol, wand_stk)                
        # Wärmeübertragung                                                              
        self.h_int = self.berech_h_int()                 # Interner Wärmeübergangskoeffizient in W/(m²·K)
        self.h_ext = 35                                  # Externer Wärmeübergangskoeffizient in W/(m²·K) (naturliche Konvektion)
        self.lambda_wand = self.berech_lambda(wand_mat)  # Wärmeleitfähigkeit des Wandmaterials in W/(m·K)
        # Physikalische Grenzen
        self.max_leistung = 5000   # maximale Heizleistung in W
        self.min_leistung = -2000  # maximale Kühlleistung in W

    def geometrie_daten(
            self, 
            reaktor_vol,
            wand_stk):
        """Berechnet und speichert alle Geometriedaten des Bioreaktors."""

        # Umrechnungen
        self.reaktor_vol = reaktor_vol / 1000                   # Volumen in m³ (Umrechnung von L in m³)
        self.wand_stk = wand_stk / 1000                         # Wandstärke in m (Umrechnung von mm in m)
        # Radien des Bioreaktors (m)
        self.r_i = (self.reaktor_vol / (2 * np.pi)) ** (1 / 3)  # Innenradius (Annahme: Zylinderform) 
        self.r_a = self.r_i + self.wand_stk                     # Außenradius 
        self.h_i = 3 * self.r_i                                 # Innenhöhe (Annahme: Annahme: H = 3 * r)
        self.h_a = self.h_i + 2 * self.wand_stk                 # Außenhöhe
        # Flächen des Bioraktors (m²)
        self.flaeche_i = 2 * np.pi * self.r_i**2 + 2 * np.pi * self.r_i * self.h_i  # Innenfläche 
        self.flaeche_a = 2 * np.pi * self.r_a**2 + 2 * np.pi * self.r_a * self.h_a  # Außenfläche
        # Rührerdurchmesses (m)
        self.ruehrer_d = (2 * self.r_i) / 3                     # Annahme: 1/3 des Innendurchmessers

    def update_stoffwerte(self,t):
        """Aktualisiert die thermophysikalischen Eigenschaften des Mediums (Wasser)."""
            
        if 0 < t < 100:  
            # Umrechnung der Temperatur in Kelvin für CoolProp
            t_kelvin = t + 273.15   
            # Stoffwerte für Wasser 
            self.spez_c = PropsSI("Cpmass", "T", t_kelvin, "P", 101325, "Water")   # Spezifische Wärmekapazität in J/(kg·K)
            self.dichte = PropsSI("Dmass", "T", t_kelvin, "P", 101325, "Water")    # Dichte in kg/m³
            self.k = PropsSI("CONDUCTIVITY", "T", t_kelvin, "P", 101325, "Water")  # Wärmeleitfähigkeit in W/(m·K)
            self.mu = PropsSI("VISCOSITY", "T", t_kelvin, "P", 101325, "Water")    # Dynamische Viskosität in Pa·s
        
        else:   
            # Fallback-Werte für ungültige Temperaturen
            self.spez_c = 4186.0    # Spezifische Wärmekapazität in J/(kg·K) (Wasser bei 20°C)
            self.dichte = 997.0     # Dichte in kg/m³ (Wasser bei 20°C)
            self.k = 0.606          # Wärmeleitfähigkeit in W/(m·K) (Wasser bei 20°C)
            self.mu = 0.001002      # Dynamische Viskosität in Pa·s (Wasser bei 20°C)

     def berech_lambda(self, wand_mat):
        """Gibt die Wärmeleitfähigkeit eines Wandmaterials in W/(m·K) zurück."""

        wand_mat = str(wand_mat).lower()  # Normalisierung für Lookup

        # Datenbank für Wärmeleitfähigkeiten in W/(m·K)
        material_db = {
            'stahl': 46.0,  
            'edelstahl': 21.0,  
            'glas': 1.4,
            'kunststoff': 0.3,
            'aluminium': 230.0
        }
        return material_db.get(wand_mat, 46.0) # Default: Stahl

    def berech_h_int(self):
        """Berechnet den internen Wärmeübergangskoeffizienten h."""

        # Dimensionslose Kennzahlen
        Re = self.drehz * (self.ruehrer_d ** 2) * self.dichte / self.mu  # Reynoldszahl 
        Pr = Prandtl(self.spez_c, self.mu, self.k)                       # Prandtl-Zahl
        # Nusselt-Zahl basierend auf Reynolds- und Prandtl-Zahl
        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:  
            Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)  # Übergangsbereich                                                                               
        elif Re >= 1e4 and Pr >= 0.6:  
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)      # Turbulente Strömung                                                                                                 
        else:  
            Nu = 3.66                                   # Laminare Strömung 
        # Berechnung des Wärmeübergangskoeffizienten h
        h = Nu * self.k / self.ruehrer_d  # Wärmeübergangskoeffizient in W/(m²·K)

        return max(h, 150)  # Mindest-Wärmeübergangskoeffizient
    



    
    def t_verlust(self):
        """Berechnet den Wärmeverlust des Reaktors basierend auf Temperaturdifferenz und Wärmeübergangskoeffizienten.
        """
        
        #
        delta_t = self.t_ist - self.t_umgebung 
        
        # Vermeidung von Division durch Null
        if abs(delta_t) < 0.01:  
            return 0.0
        
        # Berechnung der Wärmeverluste
        R_int = 1.0 / (self.h_int * self.flaeche_i)  # Interner Widerstand in K/W
        R_cond = self.wand_stk / (self.lambda_wand * self.flaeche_i)  # Widerstand der Wand in K/W
        R_ext = 1.0 / (self.h_ext * self.flaeche_a)  # Externer Widerstand in K/W
        
        # Gesamtwiderstand
        R_gesamt = R_int + R_cond + R_ext  
        q_verlust = delta_t / R_gesamt  # Wärmeverlust in W
        
        return q_verlust
    
    def update_temperatur(
                         self, 
                         leistung, 
                         zeitintervall = 1):
        """Aktualisiert die Reaktortemperatur basierend auf der Heizleistung, Kühlleistung, Wärmeverlust und Zeitintervall.
        """

        # Leistungsbegrenzung
        leistung = np.clip(leistung, self.min_leistung, self.max_leistung)  
        
        # Eigenschaften aktualisieren
        self.update_stoffwerte(self.t_ist)
        self.h_int = self.berech_h_int()
        
        # Energiebilanz
        masse = self.dichte * self.reaktor_vol  # Masse des Mediums in kg
        energie_netto = (leistung - self.t_verlust()) * zeitintervall  # Nettoenergie in J
        temp_delta = energie_netto / (masse * self.spez_c)  # Temperaturänderung in K
        
        self.t_ist += temp_delta  
        
        # Physikalische Grenzen
        self.t_ist = max(self.t_ist, 4)     # Praktische untere Grenze
        self.t_ist = min(self.t_ist, 100)   # Praktische obere Grenze
        
        return self.t_ist
    
    def reset(self, 
        t_start = None):
        """"Setzt die Ist-Temperatur auf die Starttemperatur zurück und aktualisiert die Stoffwerte.
        """
        if t_start is None:
            self.t_ist = self.t_start
        else:
            self.t_ist = t_start
        
        # Aktualisiere Stoffwerte und interne Wärmeübergangskoeffizienten
        self.update_stoffwerte(self.t_ist)

## PID-Regler-Klasse
class PID:
    """PID-Regler zur Regelung der Temperatur eines Bioreaktors.
    """
    def __init__(
                self, 
                kp = 0.0, 
                ki = 0.0, 
                kd = 0.0, 
                dt = 1.0, 
                output_min = -2000, 
                output_max = 5000):
    
        self.kp = kp  # Proportionaler Anteil
        self.ki = ki  # Integraler Anteil
        self.kd = kd  # Differentialer Anteil
        self.dt = dt  # Zeitschritt in Sekunden
        self.output_min = output_min  # Minimale Ausgangsleistung
        self.output_max = output_max  # Maximale Ausgangsleistung
        
        # Interne Zustandsvariablen
        self.fehler_vor = 0.0  # Vorheriger Fehler
        self.integral = 0.0  # Integral des Fehlers
        self.output_vor = 0.0  # Vorherige Ausgangsleistung
        
        # Anti-Windup Parameter
        self.integral_max = abs(output_max - output_min) / max(ki, 1e-6)  # Maximales Integral zur Vermeidung von Windup
        
    def run(
            self, 
            sollwert, 
            istwert):
        """Führt einen PID-Regelschritt aus und berechnet die Stellgröße.
        """

        # Fehlerberechnung
        fehler = sollwert - istwert  # Regelabweichung
        
        # P-Anteil
        p_anteil = self.kp * fehler 
        
        # I-Anteil mit Anti-Windup
        self.integral += fehler * self.dt  # Integral des Fehlers
        self.integral = np.clip(self.integral, -self.integral_max, self.integral_max)  # Begrenzung des Integrals
        i_anteil = self.ki * self.integral  
        
        # D-Anteil 
        d_anteil = self.kd * (fehler - self.fehler_vor) / self.dt  
        
        # Gesamtausgabe
        stellgroesse = p_anteil + i_anteil + d_anteil
        
        # Ausgangsbegrenzung mit Anti-Windup
        if stellgroesse > self.output_max:
            output = self.output_max  
            excess = stellgroesse - self.output_max  
            self.integral -= excess / max(self.ki, 1e-6)
        elif stellgroesse < self.output_min:
            output = self.output_min  
            excess = self.output_min - stellgroesse  
            self.integral += excess / max(self.ki, 1e-6)
        else:
            output = stellgroesse
        
        # Zustand speichern
        self.fehler_vor = fehler  
        self.output_vor = output
        
        return output
    
    def reset(self):
        """Regler zurücksetzen."""

        self.fehler_vor = 0.0
        self.integral = 0.0
        self.output_vor = 0.0






## Streamlit-Anwendung und Simulation
st.set_page_config(page_title = "Bioreaktor Temperaturregelung", layout = "wide")  
st.title("Simulation einer Temperaturregelung von Bioreaktoren")  

# Sidebar für Einstellungen
st.sidebar.title("Simulationseinstellungen")  
st.sidebar.markdown("""Passen Sie die Parameter an und beobachten Sie die Auswirkungen auf die Temperaturregelung.""")

with st.sidebar:
    st.header("Reaktorparameter")    
    reaktor_vol = st.slider("Reaktorvolumen (L)", 1, 100, 10)
    t_start = st.slider("Starttemperatur (°C)", 5, 40, 20)
    t_umgebung = st.slider("Umgebungstemperatur (°C)", 5, 40, 20)
    wand_mat = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Edelstahl", "Aluminium"])
    wand_stk = st.slider("Wandstärke (mm)", 1, 50, 5)
    drehz = st.slider("Rührerdrehzahl (1/min)", 50, 180, 100)
    
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
        
        # Erstelle Bioreaktor-Instanz mit den aktuellen Parametern
        reaktor_info = Bioreaktor(reaktor_vol, t_start, t_umgebung, drehz, wand_mat.lower(), wand_stk)
        
        st.metric("Reaktorvolumen", f"{reaktor_vol} L")
        st.metric("Reaktorradius", f"{reaktor_info.r_i*100:.1f} cm")
        st.metric("Reaktorhöhe", f"{reaktor_info.h_i*100:.1f} cm")
        st.metric("Wärmeübertr.-Fläche", f"{reaktor_info.flaeche_i:.2f} m²")
        st.metric("Wärmeleitfähigkeit", f"{reaktor_info.lambda_wand} W/(m·K)")

with col1:
        # Automatische Simulation (läuft bei jeder Parameter-Änderung)
        try:
            # Initialisierung
            reaktor_pid = Bioreaktor(reaktor_vol, t_start, t_umgebung, drehz, wand_mat.lower(), wand_stk)
            pid = PID(kp = kp, ki = ki, kd = kd, dt = dt, output_min = -2000, output_max = 5000)
            reaktor_ungeregelt = Bioreaktor(reaktor_vol, t_start, t_umgebung, drehz, wand_mat.lower(), wand_stk)
            
            # Simulation
            n_steps = int(simdauer * 60 // dt)
            zeiten = np.arange(0, n_steps * dt, dt) / 60
            
            temps_pid, temps_offen, leistungen = [], [], []
            sollwerte = []
            
            reaktor_pid.reset(t_start)
            pid.reset()
            reaktor_ungeregelt.reset(t_start)
            
            # Simulationsschleife
            for i in range(n_steps):
                t_min = zeiten[i]
                current_soll = soll_temp
            
                # PID-geregeltes System
                leistung = pid.run(current_soll, reaktor_pid.t_ist)
                temp_pid = reaktor_pid.update_temperatur(leistung, dt)
                
                # Ungeregeltes System
                temp_offen = reaktor_ungeregelt.update_temperatur(0, dt)
                
                # Speichern
                temps_pid.append(temp_pid)
                temps_offen.append(temp_offen)
                leistungen.append(leistung)
                sollwerte.append(current_soll)
            
            # Visualisierung
            st.subheader("📈 Temperaturverlauf")
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
            
            # Temperaturplot
            ax1.plot(zeiten, temps_pid, label="Mit PID-Regelung", color="tab:blue", linewidth=2)
            ax1.plot(zeiten, temps_offen, "--", label="Ohne Regelung", color="tab:orange", linewidth=2)
            ax1.plot(zeiten, sollwerte, ":", label="Solltemperatur", color="tab:red", linewidth=2)
            
            ax1.set_ylabel("Temperatur [°C]")
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            ax1.set_title("Temperaturregelung im Bioreaktor")
            
            # Leistungsplot
            ax2.plot(zeiten, leistungen, color="tab:green", linewidth=2)
            ax2.set_xlabel("Zeit [min]")
            ax2.set_ylabel("Heizleistung [W]")
            ax2.grid(True, alpha=0.3)
            ax2.set_title("Stellgröße (Heizleistung)")
            
            plt.tight_layout()
            st.pyplot(fig)
            
            # Kennzahlen
            if len(temps_pid) > 0:
                steady_state_error = abs(temps_pid[-1] - soll_temp)
                overshoot = max(0, max(temps_pid) - soll_temp)
                
                col1_metric, col2_metric, col3_metric, col4_metric = st.columns(4)
                col1_metric.metric("Endtemperatur", f"{temps_pid[-1]:.1f}°C")
                col2_metric.metric("Regelabweichung", f"{steady_state_error:.2f}°C")
                col3_metric.metric("Überschwingen", f"{overshoot:.2f}°C")
                col4_metric.metric("Max. Heizleistung", f"{max(leistungen):.0f}W")
                
        except Exception as e:
            st.error(f"Fehler bei der Simulation: {str(e)}")
            st.info("Bitte überprüfen Sie die Eingabeparameter.")
    
with tab2:
    st.subheader("📊 Detaillierte Systemanalyse")
    
    try:
        # Performance-Metriken basierend auf der aktuellen Simulation
        st.write("**Regelgüte-Kennzahlen:**")
        
        # Einschwingzeit (Zeit bis 95% des Sollwerts erreicht)
        settling_time = "Nicht erreicht"
        for i, temp in enumerate(temps_pid):
            if abs(temp - soll_temp) <= 0.05 * soll_temp:
                settling_time = f"{zeiten[i]:.1f} min"
                break
        
        # Zusätzliche Metriken
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
        
    except:
        st.info("Simulation wird geladen...")

st.sidebar.markdown("---")
st.sidebar.caption("🔄 **Live-Modus:** Ergebnisse werden automatisch aktualisiert!")

