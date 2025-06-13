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
        self.h_ext = 35                                  # Externer Wärmeübergangskoeffizient in W/(m²·K) (natürliche Konvektion)
        self.lambda_wand = self.berech_lambda(wand_mat)  # Wärmeleitfähigkeit des Wandmaterials in W/(m·K)

        # Physikalische Grenzen
        self.max_leistung = 2000   # maximale Heizleistung in W
        self.min_leistung = -1000  # maximale Kühlleistung in W

    def geometrie_daten(
            self, 
            reaktor_vol,
            wand_stk):
        """Berechnet und speichert alle Geometriedaten des Bioreaktors."""

        # Umrechnungen
        self.reaktor_vol = reaktor_vol / 1000                   # Volumen in m³ (Umrechnung von L in m³)
        self.wand_stk = wand_stk / 1000                         # Wandstärke in m (Umrechnung von mm in m)

        # Radien des Bioreaktors (m)
        self.r_i = (self.reaktor_vol / 2 * np.pi)**(1/3)        # Innenradius (Annahme: Zylinderform) 
        self.r_a = self.r_i + self.wand_stk                     # Außenradius 
        self.h_i = 4 * self.r_i                                 # Innenhöhe (Annahme: H = 4 * r)
        self.h_a = self.h_i + 2 * self.wand_stk                 # Außenhöhe

        # Flächen des Bioreaktors (m²)
        self.flaeche_i = 2 * np.pi * self.r_i**2 + 2 * np.pi * self.r_i * self.h_i  # Innenfläche 
        self.flaeche_a = 2 * np.pi * self.r_a**2 + 2 * np.pi * self.r_a * self.h_a  # Außenfläche

        # Rührerdurchmesser (m)
        self.ruehrer_d = (2 * self.r_i) / 3  # Annahme: 1/3 des Innendurchmessers

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
            self.spez_c = 4185.0    # Spezifische Wärmekapazität in J/(kg·K) (Wasser bei 20°C)
            self.dichte = 998.0     # Dichte in kg/m³ (Wasser bei 20°C)
            self.k = 0.599          # Wärmeleitfähigkeit in W/(m·K) (Wasser bei 20°C)
            self.mu = 0.001002      # Dynamische Viskosität in Pa·s (Wasser bei 20°C)

    def berech_lambda(self, wand_mat):
        """Gibt die Wärmeleitfähigkeit eines Wandmaterials in W/(m·K) zurück."""

        wand_mat = str(wand_mat).lower()  # Normalisierung für Lookup

        # Datenbank für Wärmeleitfähigkeiten in W/(m·K)
        material_db = {
            'stahl': 50.0,  
            'edelstahl': 21.0,  
            'glas': 1.4,
            'kunststoff': 0.3,
            'aluminium': 237.0
        }
        return material_db.get(wand_mat, 50.0)  # Default: Stahl

    def berech_h_int(self):
        """Berechnet den internen Wärmeübergangskoeffizienten h."""

        # Dimensionslose Kennzahlen
        Re = (self.drehz * self.ruehrer_d**2 * self.dichte) / self.mu  # Reynoldszahl 
        Pr = Prandtl(self.spez_c, self.mu, self.k)                     # Prandtl-Zahl

        # Nusselt-Zahl basierend auf Reynolds- und Prandtl-Zahl
        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:  
            Nu = 0.354 * Re**0.714 * Pr**0.260  # Übergangsbereich                                                                               
        elif Re >= 1e4 and Pr >= 0.6:  
            Nu = 0.023 * Re**0.8 * Pr**0.4      # Turbulente Strömung                                                                                                 
        else:  
            Nu = 3.66                           # Laminare Strömung 

        # Berechnung des Wärmeübergangskoeffizienten h
        h = (Nu * self.k) / self.ruehrer_d  # Wärmeübergangskoeffizient in W/(m²·K)

        return max(h, 150)  # Mindest-Wärmeübergangskoeffizient

    def t_verlust(self):
        """Berechnet den Wärmeverlust des Bioreaktors."""
        
        # Temperaturdifferenz zwischen Reaktorinhalt und Umgebung
        delta_t = self.t_ist - self.t_umgebung 
    
        # Vermeidung von Division durch nahezu Null
        if abs(delta_t) < 0.01:  
            return 0.0
        
        # Mittlere Fläche für Wärmeleitung durch die Wand
        a_mittel = (self.flaeche_i + self.flaeche_a) / 2
        
        # Berechnung der Wärmeverluste
        r_int = 1.0 / (self.h_int * self.flaeche_i)             # Interner Widerstand in K/W (Konvektion)
        r_cond = self.wand_stk / (self.lambda_wand * a_mittel)  # Widerstand der Wand in K/W (Leitung)
        r_ext = 1.0 / (self.h_ext * self.flaeche_a)             # Externer Widerstand in K/W (Konvektion)
        
        # Gesamtwiderstand
        r_gesamt = r_int + r_cond + r_ext

        # Wärmeverlust in W
        q_verlust = delta_t / r_gesamt  
        
        return q_verlust
    
    def update_temperatur(
            self, 
            leistung, 
            zeitintervall = 1):
        """Aktualisiert die Reaktortemperatur.

        Die Temperaturänderung ergibt sich aus Heizleistung, Kühlleistung,
        Wärmeverlust sowie dem betrachteten Zeitintervall.
        """
        # Leistungsbegrenzung
        leistung = np.clip(leistung, self.min_leistung, self.max_leistung)  
        
        # Eigenschaften aktualisieren
        self.update_stoffwerte(self.t_ist)
        self.h_int = self.berech_h_int()
        
        # Energiebilanz
        masse = self.dichte * self.reaktor_vol                         # Masse des Mediums in kg
        energie_netto = (leistung - self.t_verlust()) * zeitintervall  # Nettoenergie in J
        temp_delta = energie_netto / (masse * self.spez_c)             # Temperaturänderung in K
        
        self.t_ist += temp_delta  
        
        # Physikalische Grenzen
        self.t_ist = max(self.t_ist, 4)     # Praktische untere Grenze
        self.t_ist = min(self.t_ist, 100)   # Praktische obere Grenze
        
        return self.t_ist
    
    def reset(
            self, 
            t_start = None):
        """Zurücksetzen der Temperatur und Aktualisierung der Stoffwerte."""

        # Setzt die aktuelle Temperatur zurück
        if t_start is None:
            self.t_ist = self.t_start
        else:
            self.t_ist = t_start
        
        # Aktualisiert temperaturabhängige Stoffwerte und Wärmeübergangskoeffizienten
        self.update_stoffwerte(self.t_ist)


## PID-Regler-Klasse
class PID:
    """PID-Regler zur Regelung der Temperatur eines Bioreaktors."""

    def __init__(
            self, 
            kp = 0.0, 
            ki = 0.0, 
            kd = 0.0, 
            dt = 1.0, 
            output_min = -1000, 
            output_max = 2000):
        """Initialisiert die PID-Reglerparameter und interne Zustände."""

        # PID-Regelparameter
        self.kp = kp  # Proportionalanteil
        self.ki = ki  # Integralanteil
        self.kd = kd  # Differentialanteil
        self.dt = dt  # Zeitschritt in s

        # Begrenzung der Reglerausgabe
        self.output_min = output_min  # Minimale Ausgangsleistung
        self.output_max = output_max  # Maximale Ausgangsleistung
        
        # Interne Zustandsvariablen
        self.fehler_vor = 0.0  # Vorheriger Fehler
        self.integral = 0.0    # Integral des Fehlers
        self.output_vor = 0.0  # Vorherige Ausgangsleistung
        
        # Anti-Windup Parameter
        self.integral_max = abs(output_max - output_min) / max(ki, 1e-6)  # Maximales Integral zur Vermeidung von Windup
        
    def run(
            self, 
            sollwert, 
            istwert):
        """Berechnet die Stellgröße anhand eines PID-Regelschritts."""

        # Fehlerberechnung (Regelabweichung)
        fehler = sollwert - istwert  
        
       # Proportionalanteil
        p_anteil = self.kp * fehler
        
        # Integralanteil mit Anti-Windup
        self.integral += fehler * self.dt  
        self.integral = np.clip(self.integral, -self.integral_max, self.integral_max)  
        i_anteil = self.ki * self.integral  
        
        # Integralanteil mit Anti-Windup
        d_anteil = self.kd * (fehler - self.fehler_vor) / self.dt  
        
        # PID-Gesamtausgabe
        stellgroesse = p_anteil + i_anteil + d_anteil
        
        # Ausgangsbegrenzung mit Anti-Windup
        if stellgroesse > self.output_max:
            output = self.output_max  
            abw_output = stellgroesse - self.output_max  
            self.integral -= abw_output / max(self.ki, 1e-6)
        elif stellgroesse < self.output_min:
            output = self.output_min  
            abw_output = self.output_min - stellgroesse  
            self.integral += abw_output / max(self.ki, 1e-6)
        else:
            output = stellgroesse
        
       # Zustand für nächste Iteration speichern
        self.fehler_vor = fehler  
        self.output_vor = output
        
        return output
    
    def reset(self):
        """Setzt den PID-Regler in den Ausgangszustand zurück."""

        self.fehler_vor = 0.0
        self.integral = 0.0
        self.output_vor = 0.0


## Streamlit-Anwendung und Simulation

# Streamlit konfigurieren
st.set_page_config(
    page_title = "Bioreaktor Temperaturregelung", 
    layout = "wide") 
st.image(
    "https://upload.wikimedia.org/wikipedia/commons/8/8b/HSWT_Logo_gruen.png",
    width = 500)                                                     
st.title("Simulation einer Temperaturregelung von Bioreaktoren") 

# Sidebar für Einstellungen
st.sidebar.title("Simulationseinstellungen")  
st.sidebar.markdown(
    "Passen Sie die Parameter an, " \
    "um das Verhalten der Temperaturregelung zu analysieren.")

# Sidebar: Eingabe der Simulationsparameter
with st.sidebar:

    # Reaktor-Konfiguration: Initiale Betriebsparameter
    st.header("Reaktorparameter") 
    reaktor_vol = st.slider("Reaktorvolumen (L)", 1, 20, 5)
    t_start = st.slider("Starttemperatur (°C)", 5, 40, 20)
    t_umgebung = st.slider("Umgebungstemperatur (°C)", 5, 40, 20)
    wand_mat = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Edelstahl", "Aluminium"])
    wand_stk = st.slider("Wandstärke (mm)", 1, 20, 5)
    drehz = st.slider("Rührerdrehzahl (1/min)", 1, 800, 100)
    
    # Sollwert und Simulationssteuerung
    st.header("Sollwert & Simulation")
    soll_temp = st.slider("Solltemperatur (°C)", 5, 60, 37)
    simdauer = st.slider("Simulationsdauer (min)", 10, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 30, 5)
    
    # PID-Reglerparameter zur Steuerung der Temperaturregelung
    st.header("PID-Parameter")
    kp = st.slider("Kp (P-Anteil)", 0.0, 200.0, 50.0, 1.0,
        help = "Kp × e(t): Lineare Reaktion auf Regelfehler e(t)")
    ki = st.slider("Ki (I-Anteil)", 0.0, 10.0, 1.0, 0.1, 
        help = "Ki × ∫e(τ)dτ: Korrektur über aufsummierten Fehlerverlauf")
    kd = st.slider("Kd (D-Anteil)", 0.0, 50.0, 5.0, 0.5, 
        help = "Kp × e(t): Lineare Reaktion auf Regelfehler e(t)")
    
# Tabs für Simulation und Analyse (zur strukturierten Darstellung von Inhalten)
tab1, tab2 = st.tabs(["📊 Simulation", "📋 Analyse"])

# Tab 1: Interaktive Eingabe und Visualisierung der Reaktorsimulation
with tab1:
    # Zwei-Spalten-Layout 
    col1, col2 = st.columns([2, 1])  # linke Spalte für Visualisierungen, rechte Spalte für Parameteranzeige
    
    with col2:
        st.subheader("Reaktor-Eigenschaften")
        
        # Erstellt eine Bioreaktor-Instanz mit den aktuellen Einstellungen (Sidebar)
        reaktor_info = Bioreaktor(
            reaktor_vol, t_start, t_umgebung, drehz, 
            wand_mat.lower(), wand_stk) 
        
        # Anzeige der wichtigsten Reaktorparameter als Kennzahlen (metrische Kacheln)
        st.metric("Reaktorvolumen", f"{reaktor_vol} L")                 # Eingestelltes Reaktorvolumen in Litern
        st.metric("Reaktorradius", f"{reaktor_info.r_i*100:.1f} cm")    # Berechneter Innenradius in cm (Umrechnung aus m)
        st.metric("Reaktorhöhe", f"{reaktor_info.h_i*100:.1f} cm")      # Berechnete Innenhöhe in cm (Umrechnung aus m)
        st.metric("Wärmeübertragungsfläche", f"{reaktor_info.flaeche_i:.2f} m²")    # Berechnete innere Wärmeübertragungsfläche in m²
        st.metric("Wärmeleitfähigkeit", f"{reaktor_info.lambda_wand} W/(m·K)")      # Wärmeleitfähigkeit des Wandmaterials in W/(m·K)

    with col1:
        # Automatische Simulation (läuft bei jeder Parameter-Änderung

        # Initialisierung des geregelten Reaktors (kontrolliert)
        reaktor_pid = Bioreaktor(#
            reaktor_vol, t_start, t_umgebung, 
            drehz, wand_mat.lower(), wand_stk)
            
        # Initialisierung des PID-Reglers
        pid = PID(kp = kp, ki = ki, kd = kd, dt = dt, 
            output_min = -1000, output_max = 2000)
            
        # Initialisierung des ungeregelten Reaktors (unkontrolliert)
        reaktor_ungeregelt = Bioreaktor(
            reaktor_vol, t_start, t_umgebung, 
            drehz, wand_mat.lower(), wand_stk)
            
        # Vorbereitung der Simulation: Zeitbasis und Initialisierung
        n_steps = int(simdauer * 60 // dt)            
        zeiten = np.arange(0, n_steps * dt, dt) / 60  # Zeitachse in Minuten
            
        # Leere Listen zur Speicherung der Simulationsverläufe
        temps_pid, temps_offen, leistungen = [], [], [] # Temperaturverlauf und Heiz-/Kühlleistung
        sollwerte = []  # Solltemperaturverlauf
            
        # Startwerte setzen
        reaktor_pid.reset(t_start)
        pid.reset()
        reaktor_ungeregelt.reset(t_start)
            
        # Simulationsschleife: berechnet Temperaturverläufe und Regelgröße
        for i in range(n_steps):
            zeit_min = zeiten[i]

            # Aktueller Sollwert (später ggf. zeitabhängig)
            t_soll = soll_temp  
            
            # PID-geregeltes System: Temperaturänderung mit Regler
            leistung = pid.run(t_soll, reaktor_pid.t_ist)
            t_pid = reaktor_pid.update_temperatur(leistung, dt)
                
            # Ungeregeltes System: Temperaturänderung ohne Regler
            t_offen = reaktor_ungeregelt.update_temperatur(0, dt)
                
            # Ergebnisse speichern: Temperaturverläufe, Leistung und Sollwert
            temps_pid.append(t_pid)
            temps_offen.append(t_offen)
            leistungen.append(leistung)
            sollwerte.append(t_soll)
            
        # Visualisierung: Temperaturverlauf und Heizleistung
        st.subheader("📈 Temperaturverlauf")
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12, 8))  # zwei vertiale Subplots
            
        # Temperaturplot: Temperaturverlauf in beiden Systemen
        ax1.plot(
            zeiten, temps_pid, 
            label = "Mit PID-Regelung", 
            color = "tab:blue", 
            linewidth = 2)
        ax1.plot(
            zeiten, temps_offen,
            "--", label = "Ohne Regelung", 
            color = "tab:orange", 
            linewidth = 2)
        ax1.plot(
            zeiten, sollwerte, 
            ":", label = "Solltemperatur", 
            color = "tab:red", 
            linewidth = 2)
            
        ax1.set_ylabel("Temperatur [°C]")
        ax1.legend()                 # Legende für die Temperaturkurven 
        ax1.grid(True, alpha = 0.3)  # Dezentes Gitter für bessere Lesbarkeit
        ax1.set_title("Temperaturregelung im Bioreaktor")
            
        # Leistungsplot: Heiz- und Kühlleistung über die Zeit
        ax2.plot(zeiten, leistungen, color = "tab:green", linewidth = 2)
        ax2.set_xlabel("Zeit [min]")
        ax2.set_ylabel("Heiz-/Kühlleistung [W]")
        ax2.grid(True, alpha = 0.3)
        ax2.set_title("Stellgröße (Heiz-/Kühlleistung)")
            
        plt.tight_layout()
        st.pyplot(fig)  # Anzeige des Plots in Streamlit
            
        # Kennzahlen berechnen und anzeigen, falls Daten vorhanden
        if len(temps_pid) > 0:
            stat_Fehler = abs(temps_pid[-1] - soll_temp)     # Regelabweichung am Ende
            ue_schwing = max(0, max(temps_pid) - soll_temp)  # Überschwingen über Sollwert
                
            # Streamlit-Metriken: Darstellung der wichtigsten Regelgüte-Kennzahlen   
            col1_metric, col2_metric, col3_metric, col4_metric = st.columns(4)

            col1_metric.metric("Endtemperatur", f"{temps_pid[-1]:.1f}°C")
            col2_metric.metric("Regelabweichung", f"{stat_Fehler:.2f}°C")
            col3_metric.metric("Überschwingen", f"{ue_schwing:.2f}°C")
            col4_metric.metric("Max. Heizleistung", f"{max(leistungen):.0f}W")
    
# Tab 2: Darstellung detaillierter Analysemetriken der Regelgüte    
with tab2:
    st.subheader("Detaillierte Systemanalyse")
    
    try:
        # Performance-Metriken auf Basis der aktuellen Simulation
        st.write("**Regelgüte-Kennzahlen:**")
        
        # Einschwingzeit: Zeit bis die Temperatur innerhalb von ±5 % des Sollwerts liegt
        schw_zeit = "Nicht erreicht"
        for i, temp in enumerate(temps_pid):
            if abs(temp - soll_temp) <= 0.05 * soll_temp:
                schw_zeit = f"{zeiten[i]:.1f} min"
                break
        
        # Anstiegszeit: Zeit bis 90 % des Sollwerts erreicht sind
        anstieg_zeit = "Nicht erreicht"
        for i, temp in enumerate(temps_pid):
            if temp >= 0.9 * soll_temp:
                anstieg_zeit = f"{zeiten[i]:.1f} min"
                break
        
        # Anzeige von Regelgüte-Metriken in zwei Spalten
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Einschwingzeit (95%)", schw_zeit)   # Zeit bis Temperatur 95 % des Sollwerts erreicht
            st.metric("Anstiegszeit (90%)", anstieg_zeit)  # Zeit bis 90 % des Sollwerts erreicht werden
            st.metric(
                "Mittlere Abweichung", 
                f"{np.mean(np.abs(np.array(temps_pid) - soll_temp)):.2f} °C")  # Durchschnittlicher Regelabweichungsbetrag
        
        with col2:
            st.metric("Standardabweichung", f"{np.std(temps_pid):.2f} °C")  # Schwankung um den Mittelwert
            st.metric(
                "Energieverbrauch", 
                f"{np.mean(np.maximum(leistungen, 0)) * simdauer * 60/1000:.1f} kJ")       # Gesamtenergie für Heizen
            st.metric(
                "Kühlenergie", 
                f"{abs(np.mean(np.minimum(leistungen, 0))) * simdauer * 60/1000:.1f} kJ")  # Gesamtenergie für Kühlen
        
       # Visualisierung der Regelabweichung über die Zeit
        st.subheader("📈 Regelabweichung über Zeit")

        fig_error, ax_error = plt.subplots(figsize=(10, 4))
        fehler = np.array(temps_pid) - soll_temp  # Abweichung Ist-Temperatur von Sollwert

        # Regelabweichung plotten
        ax_error.plot(zeiten, fehler, color = "red", linewidth = 2)
        ax_error.axhline(0, color = "black", linestyle = "--", alpha = 0.5)
        ax_error.fill_between(zeiten, fehler, alpha = 0.3, color = "red")
        ax_error.set_xlabel("Zeit [min]")
        ax_error.set_ylabel("Regelabweichung [°C]")
        ax_error.set_title("Abweichung von der Solltemperatur")
        ax_error.grid(True, alpha = 0.3)
        st.pyplot(fig_error)
        
    except Exception:
        st.info("Simulation wird geladen...")

st.sidebar.markdown("---")
st.sidebar.caption(
    "🔄 **Live-Modus:** Ergebnisse werden automatisch aktualisiert!")

