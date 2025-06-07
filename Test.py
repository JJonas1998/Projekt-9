#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
 # Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen = 10, start_temp = 20, umg_temp = 20, rpm = 100, wandmaterial = 'stahl', wandstaerke = 5):

        # Speichere die Starttemperatur für reset()
        self.start_temp = start_temp

        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        temp_k = start_temp + 273.15                                                                # Temperatur in Kelvin
        self.spez_c = PropsSI("Cpmass", "T", (temp_k), "P", 101325, "Water")                        # Spezifische Wärmekapazität des Mediums(Wasser) in J/(kg*K) 
        self.dichte = PropsSI("Dmass" ,"T", (temp_k), "P", 101325, "Water")                         # Dichte des Mediums (Wasser) in kg/m³

        # Geometrie des Bioreaktors (Zylinder, Annahme: H = 2 * r)
        self.volumen = volumen / 1000                                                                # Volumen des Bioreaktors in m³ (Liter → m³)
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)                                         # Radius des Bioreaktors in m 
        self.hoehe  = 2 * self.radius                                                               # Höhe des Bioreaktors in m 
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe        # Innenfläche des Bioreaktor in m² 
        self.wandstaerke = wandstaerke / 1000                                                       # Wandstärke des Bioreaktors in m
        self.d_ruehrer = (2 * self.radius) / 3                                                      # Durchmesser des Rührers in m (Annahme: Rührerdurchmesser 1/3 des Reaktordurchmessers) 

        # Betriebsbedingungen
        self.umg_temp = umg_temp                                                                    # Umgebungstemperatur in °C
        self.ist_temp = start_temp                                                                  # Reaktor-Innentemperatur in °C (Starttemperatur)
        self.rpm = rpm                                                                              # Drehzahl des Rührers in 1/min

        # Wärmeübertragung                                                              
        self.h_intern = self.berechnung_h()                                                         # Berechnung des Wärmeübergangskoeffizienten in W/(m²*K)   
        self.h_extern = 1.0                                                                        # Externer Wärmeübergangskoeffizient (Luft) in W/(m²*K)  
        self.lambda_wand = self.get_waermeleitfaehigkeit(wandmaterial)                              # Wärmeleitfähigkeit des Wandmaterials in W/(m*K)                   
        
    def berechnung_h(self):
        """
        Berechnung des Wärmeübergangskoeffizienten h anhand von Impeller-Drehzahl (rpm) und Rührerdurchmesser.
        """
        t_k = self.ist_temp + 273.15                                                                # Umrechnung der Innentemperatur in Kelvin   
        
        # Thermophysikalische Eigenschaften des Fluids                                                                                              
        k   = PropsSI("CONDUCTIVITY", "T", t_k, "P", 101325, "Water")                               # Wärmeleitfähigkeit des Mediums (Wasser) in W/(m*K)    
        mu  = PropsSI("VISCOSITY" ,"T", t_k, "P", 101325, "Water")                                  # Dynamische Viskosität des Mediums (Wasser) in Pa*s    

        # Dimensionlose Kennzahlen 
        Re = ((self.rpm / 60) * (self.d_ruehrer ** 2) * self.dichte) / mu                           # Reynolds-Zahl des Mediums (Wasser) im Bioreaktor    
        Pr = Prandtl(self.spez_c, mu, k)                                                            # Prandtl-Zahl des Mediums (Wasser) im Bioreaktor    

        # Nusselt-Zahl und Wärmeübergangskoeffizient
        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:
            Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)                                              # Nusselt-Zahl für turbulente Strömung (Impeller)                                       
        elif Re >= 1e4 and Pr >= 0.6:
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)                                                  # Nusselt-Zahl für turbulente Strömung (Dittus-Boelter)                                                  
        else:
            Nu = 3.66                                                                               # Nusselt-Zahl für laminare Strömung       
                                                                                                                                                                              
        return Nu * k / self.d_ruehrer    

    def get_waermeleitfaehigkeit(self, material):
        """
        Gibt die Wärmeleitfähigkeit [W/mK] für ein gegebenes Wandmaterial zurück.
        """
        material = str(material).lower()
        material_db = {
            'stahl': 21.0,                                                                          # Edelstahl V2A
            'glas': 1.4,                                                                            # Quarzglas     
        }                                                       
        return material_db.get(material, 16.0)                                                      # default: Stahl     

    def berechnung_waermeverlust(self):

        delta_t = self.ist_temp - self.umg_temp

        A = self.flaeche
        d = self.wandstaerke
        λ = self.lambda_wand
        h_int = self.h_intern                                                                       # Interner Wärmeübergangskoeffizient in W/(m²*K)
        h_ext = self.h_extern                                                                       # Externer Wärmeübergangskoeffizient in W/(m²*K)

        # Richtige Widerstände: 1/(h*A) bzw. d/(λ*A)
        R_int  = 1.0 / (h_int * A)
        R_cond = d / (λ * A)
        R_ext  = 1.0 / (h_ext * A)

        R_gesamt = R_int + R_cond + R_ext
        q_verlust = delta_t / R_gesamt

        return q_verlust
    
    def update_temperature(self, leistung, zeitintervall = 1):
        """
        Aktualisiert die Innentemperatur des Reaktors auf Basis der Energiebilanz:
        ΔT = (Q_zufuhr - Q_verlust) * dt / (m * c)
        """
        # Aktualisierung der temperaturabhängigen Eigenschaften
        temp_k = self.ist_temp + 273.15
        self.spez_c = PropsSI("Cpmass", "T", temp_k, "P", 101325, "Water")
        self.dichte = PropsSI("Dmass", "T", temp_k, "P", 101325, "Water")
        self.h_intern = self.berechnung_h()                                                         # Aktualisierung des internen Wärmeübergangskoeffizienten in W/(m²*K)
        
        # Berechnung der Masse des Reaktorinhalts und Wärmeverlust
        masse = self.dichte * self.volumen                                                          # Masse des Reaktorinhalts in kg 
        waermeverlust = self.berechnung_waermeverlust()                                             # Gesamtwärmeverlust in W 
        energie_netto = (leistung - waermeverlust) * zeitintervall                                  # Netto-Energiezufuhr in J (Leistung - Wärmeverlust) * Zeitintervall in s
        temp_delta = energie_netto / (masse * self.spez_c)                                          # Temperaturänderung in K 
        self.ist_temp += temp_delta                                                                 # Aktualisierung der Innentemperatur des Reaktors in °C

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
    def __init__(self, kp = 0.0, ki = 0.0, kd = 0.0, dt = 0.0, output_min = -1000, output_max = 1000):
        self.dt = dt                                                                                # Zeitschritt in Sekunden
        self.kp = kp                                                                                # Proportionalanteil  
        self.ki = ki                                                                                # Integralanteil    
        self.kd = kd                                                                                # Differentialanteil
        self.output_min = output_min                                                                # Minimalwert der Stellgröße (Heizleistung)
        self.output_max = output_max                                                                # Maximalwert der Stellgröße (Heizleistung)
        self.fehler_vor = 0.0                                                                       # Vorheriger Fehler für den D-Anteil                                                                          
        self.integral = 0.0                                                                         # Integralwert für den I-Anteil
        
    def run(self, sollwert, istwert):
        """
        Berechnet PID-Ausgabe
        
        Args:
            sollwert: Zieltemperatur
            istwert: Aktuelle Temperatur
            
        Returns:
            Stellgröße (Heizleistung)
        """
        # Fehlerberechnung
        fehler = sollwert - istwert                                                                  # Berechnung des Fehlers (Sollwert - Istwert)
        
        # P-Anteil
        p_anteil = self.kp * fehler                                                                  # Proportionalanteil
        
        # I-Anteil
        self.integral += fehler * self.dt                                                            # Integralwert aktualisieren
        i_anteil = self.ki * self.integral                                                           # Integralanteil  
        
        # D-Anteil
        if self.dt > 0:
            d_anteil = self.kd * (fehler - self.fehler_vor) / self.dt                               # Differentialanteil
            
        # Ausgabe berechnen
        stellgroesse = p_anteil + i_anteil + d_anteil                                               # PID-Ausgabe (Stellgröße)
        
        # Begrenzung mit Anti-Windup
        if stellgroesse > self.output_max:
            output = self.output_max
            # Anti-Windup: Integral zurücksetzen wenn Sättigung erreicht
            if fehler > 0:  # Nur bei positivem Fehler
                self.integral -= fehler * self.dt
        elif stellgroesse < self.output_min:
            output = self.output_min
            # Anti-Windup: Integral zurücksetzen wenn Sättigung erreicht
            if fehler < 0:  # Nur bei negativem Fehler
                self.integral -= fehler * self.dt
        else:
            output = stellgroesse                                                                    # Stellgröße innerhalb der Grenzen

        self.fehler_vor = fehler                                                                     # Aktuellen Fehler speichern
        
        return output
    
    def reset(self):
        """Regler zurücksetzen"""
        self.fehler_vor = 0.0
        self.integral = 0.0

## Streamlit-Anwendung und Simulation
 # Streamlit-UI

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
    volumen = st.slider("Reaktorvolumen (l)", 1, 100, 5)
    start_temp = st.slider("Starttemperatur (°C)", 0, 40, 20)
    umg_temp = st.slider("Umgebungstemperatur (°C)", 0, 40, 20)
    wandmaterial = st.selectbox("Wandmaterial", ["Stahl", "Glas"])                                                                         
    wandstaerke = st.slider("Wandstärke (mm)", 1, 100, 5)                         
    soll_temp = st.slider("Solltemperatur (°C)", 20, 80, 37)
    simdauer = st.slider("Simulationsdauer (min)", 1, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 60, 5)
    rpm = st.slider("Rührerdrehzahl (1/min)", 1, 240, 60)
    st.header("PID-Parameter")
    kp = st.slider("Kp (proportional)", 0.0, 100.0, 10.0, 0.1)
    ki = st.slider("Ki (integral)", 0.0, 5.0, 0.5, 0.01)
    kd = st.slider("Kd (differential)", 0.0, 100.0, 1.0, 0.01)
    st.markdown("---")
    st.caption("Tipp: Variiere die Reglerparameter und beobachte die Auswirkungen auf die Temperaturregelung.")

# Initialisierung der Bioreaktor- und PID-Regler-Objekte
reaktor_pid = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial, wandstaerke)
pid = PID(kp=kp, ki=ki, kd=kd, dt=dt)
reaktor_ungeregelt = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial, wandstaerke)  # Ungeregelter Reaktor für Vergleich

# --- Simulation ---
n_steps = int(simdauer * 60 // dt)
zeiten = np.arange(0, n_steps * dt, dt) / 60  
temps_pid, temps_offen, leistungen = [], [], []

reaktor_pid.reset(start_temp) 
pid.reset()
reaktor_ungeregelt.reset(start_temp)  # Reset für unregulierten Reaktor

for t in range(n_steps):
    # PID-geregeltes System
    leistung = pid.run(soll_temp, reaktor_pid.ist_temp)
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