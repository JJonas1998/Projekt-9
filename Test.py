#### Projekt 9: Simulation einer Temperaturregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

## Bioreaktor-Klasse
class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen=10, start_temp=20, umg_temp=20, rpm=100, wandmaterial='stahl', wandstaerke=0.005):

        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        temp_k = start_temp + 273.15
        self.volumen = volumen / 1000
        self.spez_c = PropsSI("Cpmass", "T", (temp_k), "P", 101325, "Water")
        self.dichte = PropsSI("Dmass" ,"T", (temp_k), "P", 101325, "Water")

        # Geometrie des Bioreaktors (Zylinder, Annahme: H = 2 * r)
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)
        self.hoehe  = 2 * self.radius
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe
        self.wandstaerke = wandstaerke
        self.d_ruehrer = (2 * self.radius) / 3

        # Betriebsbedingungen
        self.umg_temp = umg_temp
        self.ist_temp = start_temp
        self.rpm = rpm

        # Wärmeübertragung
        self.waerme_h = self.berechnung_h()
        self.lambda_wand = self.get_waermeleitfaehigkeit(wandmaterial)
        
    def berechnung_h(self):
        """
        Berechnung des Wärmeübergangskoeffizienten h anhand von Impeller-Drehzahl (rpm) und Rührerdurchmesser.
        """
        t_k = self.ist_temp + 273.15
        
        # Thermophysikalische Eigenschaften des Fluids
        k = PropsSI("CONDUCTIVITY", "T", t_k, "P", 101325, "Water")
        mu = PropsSI("VISCOSITY" ,"T", t_k, "P", 101325, "Water")

        # Dimensionlose Kennzahlen 
        Re = ((self.rpm / 60) * (self.d_ruehrer ** 2) * self.dichte) / mu
        Pr = Prandtl(self.spez_c, mu, k)

        # Nusselt-Zahl und Wärmeübergangskoeffizient
        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:
            Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)
        elif Re >= 1e4 and Pr >= 0.6:
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
        else:
            Nu = 3.66
                                                                                                                                                                              
        return Nu * k / self.d_ruehrer    

    def get_waermeleitfaehigkeit(self, material):
        """
        Gibt die Wärmeleitfähigkeit [W/mK] für ein gegebenes Wandmaterial zurück.
        """
        material_db = {
            'stahl': 21.0,
            'glas': 1.4,
            'aluminium': 230.0
        }
        return material_db.get(material.lower(), 16.0)

    def berechnung_waermeverlust(self):
        """
        Berechnung des gesamten Wärmeverlusts: 
        - durch die Wand (Fourier)
        - durch Konvektion mit h (z. B. Luft/Wasserbad)
        """
        # Temperaturdifferenz zwischen Innen- und Außentemperatur
        delta_T = self.ist_temp - self.umg_temp

        # 1. Leitung durch Wandmaterial
        q_leitung = self.lambda_wand * (delta_T * self.flaeche / self.wandstaerke)

        # 2. Konvektiver Verlust an Umgebung
        q_konvektion = self.waerme_h * self.flaeche * delta_T

        # Gesamtverlust
        q_verlust = q_leitung + q_konvektion

        return q_verlust

    def update_temperature(self, leistung, dt=1):
        """
        Aktualisiert die Innentemperatur des Reaktors auf Basis der Energiebilanz:
        ΔT = (Q_zufuhr - Q_verlust) * dt / (m * c)
        """
        # Berechnung der Masse des Reaktorinhalts und Wärmeverlust
        masse = self.dichte * self.volumen
        waermeverlust = self.berechnung_waermeverlust()
        energie_netto = (leistung - waermeverlust) * dt
        temp_delta = energie_netto / (masse * self.spez_c)
        self.ist_temp += temp_delta

        return self.ist_temp

## PID-Regler-Klasse
class PIDRegler:
    """
    PID-Regler zur Regelung der Temperatur eines Bioreaktors.
    """
    
    def __init__(self, Kp=100.0, Ki=10.0, Kd=5.0, sollwert=37.0, max_output=5000.0, min_output=-2000.0):
        # PID-Parameter
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        
        # Sollwert und Ausgangsbegrenzung
        self.sollwert = sollwert
        self.max_output = max_output
        self.min_output = min_output
        
        # Interne Zustandsvariablen
        self.fehler_integral = 0.0
        self.letzter_fehler = 0.0
        self.letzter_istwert = None
        
        # Anti-Windup Parameter
        self.integral_max = 1000.0
        self.integral_min = -1000.0
        
    def set_sollwert(self, neuer_sollwert):
        """
        Setzt einen neuen Sollwert für die Temperaturregelung.
        """
        self.sollwert = neuer_sollwert
        
    def set_pid_parameter(self, Kp=None, Ki=None, Kd=None):
        """
        Aktualisiert die PID-Parameter.
        """
        if Kp is not None:
            self.Kp = Kp
        if Ki is not None:
            self.Ki = Ki
        if Kd is not None:
            self.Kd = Kd
            
    def reset(self):
        """
        Setzt die internen Zustandsvariablen zurück.
        """
        self.fehler_integral = 0.0
        self.letzter_fehler = 0.0
        self.letzter_istwert = None
        
    def berechne_stellgroesse(self, istwert, dt=1.0):
        """
        Berechnet die Stellgröße (Heiz-/Kühlleistung) basierend auf der aktuellen Temperatur.
        """
        # Regelabweichung berechnen
        fehler = self.sollwert - istwert
        
        # Initialisierung beim ersten Aufruf
        if self.letzter_istwert is None:
            self.letzter_istwert = istwert
            self.letzter_fehler = fehler
            
        # P-Anteil
        p_anteil = self.Kp * fehler
        
        # I-Anteil
        self.fehler_integral += fehler * dt
        self.fehler_integral = max(self.integral_min, 
                                 min(self.integral_max, self.fehler_integral))
        i_anteil = self.Ki * self.fehler_integral
        
        # D-Anteil
        if dt > 0:
            fehler_ableitung = (fehler - self.letzter_fehler) / dt
        else:
            fehler_ableitung = 0.0
        d_anteil = self.Kd * fehler_ableitung
        
        # Gesamte Stellgröße
        stellgroesse = p_anteil + i_anteil + d_anteil
        stellgroesse = max(self.min_output, min(self.max_output, stellgroesse))
        
        # Zustandsvariablen aktualisieren
        self.letzter_fehler = fehler
        self.letzter_istwert = istwert
        
        return stellgroesse
    
    def get_pid_anteile(self, istwert, dt=1.0):
        """
        Gibt die einzelnen PID-Anteile zurück.
        """
        fehler = self.sollwert - istwert
        
        if self.letzter_istwert is None:
            return {'fehler': fehler, 'p_anteil': 0, 'i_anteil': 0, 'd_anteil': 0, 'stellgroesse': 0}
            
        p_anteil = self.Kp * fehler
        i_anteil = self.Ki * self.fehler_integral
        
        if dt > 0:
            fehler_ableitung = (fehler - self.letzter_fehler) / dt
        else:
            fehler_ableitung = 0.0
        d_anteil = self.Kd * fehler_ableitung
        
        stellgroesse_unbegrenzt = p_anteil + i_anteil + d_anteil
        stellgroesse = max(self.min_output, min(self.max_output, stellgroesse_unbegrenzt))
        
        return {
            'fehler': fehler,
            'p_anteil': p_anteil,
            'i_anteil': i_anteil,
            'd_anteil': d_anteil,
            'stellgroesse_unbegrenzt': stellgroesse_unbegrenzt,
            'stellgroesse': stellgroesse,
            'integral': self.fehler_integral
        }

## Streamlit Anwendung
def main():
    st.set_page_config(page_title="Bioreaktor Temperaturregelung", layout="wide")
    
    st.title("🧪 Bioreaktor Temperaturregelung Simulation")
    st.markdown("---")
    
    # Sidebar für Parameter
    st.sidebar.header("⚙️ Reaktor-Parameter")
    
    # Reaktor-Konfiguration
    volumen = st.sidebar.slider("Reaktorvolumen [L]", 1, 100, 10)
    start_temp = st.sidebar.slider("Starttemperatur [°C]", 15, 40, 20)
    umg_temp = st.sidebar.slider("Umgebungstemperatur [°C]", 10, 35, 20)
    rpm = st.sidebar.slider("Rührerdrehzahl [rpm]", 50, 500, 100)
    wandmaterial = st.sidebar.selectbox("Wandmaterial", ["stahl", "glas", "aluminium"])
    wandstaerke = st.sidebar.slider("Wandstärke [mm]", 1, 20, 5) / 1000
    
    # PID-Parameter
    st.sidebar.header("🎛️ PID-Parameter")
    kp = st.sidebar.slider("Kp (Proportional)", 0.1, 500.0, 100.0, 0.1)
    ki = st.sidebar.slider("Ki (Integral)", 0.0, 100.0, 10.0, 0.1)
    kd = st.sidebar.slider("Kd (Differential)", 0.0, 50.0, 5.0, 0.1)
    sollwert = st.sidebar.slider("Solltemperatur [°C]", 20, 80, 37)
    
    # Simulationsparameter
    st.sidebar.header("⏱️ Simulation")
    sim_dauer = st.sidebar.slider("Simulationsdauer [min]", 10, 180, 60)
    dt = st.sidebar.slider("Zeitschritt [s]", 1, 10, 1)
    
    # Störungen
    st.sidebar.header("🔧 Störungen")
    stoerung_aktiviert = st.sidebar.checkbox("Temperaturstörung aktivieren")
    stoerung_zeit = st.sidebar.slider("Störung nach [min]", 5, sim_dauer-10, 30)
    stoerung_groesse = st.sidebar.slider("Störungsgröße [°C]", -10, 10, -5)
    
    # Simulation starten
    if st.sidebar.button("🚀 Simulation starten", type="primary"):
        
        # Initialisierung
        reaktor = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial, wandstaerke)
        regler = PIDRegler(kp, ki, kd, sollwert)
        
        # Datensammlung
        zeiten = []
        temperaturen = []
        sollwerte = []
        stellgroessen = []
        p_anteile = []
        i_anteile = []
        d_anteile = []
        waermeverluste = []
        
        # Progress Bar
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Simulation
        t = 0
        while t <= sim_dauer * 60:  # Umrechnung in Sekunden
            
            # Störung anwenden
            if stoerung_aktiviert and abs(t - stoerung_zeit * 60) < dt:
                reaktor.ist_temp += stoerung_groesse
                st.sidebar.success(f"Störung angewendet bei t={t/60:.1f} min")
            
            # PID-Regelung
            stellgroesse = regler.berechne_stellgroesse(reaktor.ist_temp, dt)
            
            # Reaktor aktualisieren
            neue_temp = reaktor.update_temperature(stellgroesse, dt)
            
            # PID-Anteile für Analyse
            pid_info = regler.get_pid_anteile(reaktor.ist_temp, dt)
            
            # Daten sammeln
            zeiten.append(t / 60)  # Zeit in Minuten
            temperaturen.append(reaktor.ist_temp)
            sollwerte.append(sollwert)
            stellgroessen.append(stellgroesse)
            p_anteile.append(pid_info['p_anteil'])
            i_anteile.append(pid_info['i_anteil'])
            d_anteile.append(pid_info['d_anteil'])
            waermeverluste.append(reaktor.berechnung_waermeverlust())
            
            # Progress aktualisieren
            progress = min(t / (sim_dauer * 60), 1.0)
            progress_bar.progress(progress)
            status_text.text(f"Zeit: {t/60:.1f}/{sim_dauer} min | Temp: {reaktor.ist_temp:.1f}°C")
            
            t += dt
        
        progress_bar.empty()
        status_text.empty()
        
        # Ergebnisse anzeigen
        st.success("✅ Simulation abgeschlossen!")
        
        # Hauptdiagramm mit Subplots
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=('Temperaturverlauf', 'Stellgröße', 
                          'PID-Anteile', 'Wärmeverlust',
                          'Regelabweichung', 'Systemparameter'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"type": "table"}]]
        )
        
        # Temperaturverlauf
        fig.add_trace(go.Scatter(x=zeiten, y=temperaturen, name='Ist-Temperatur', 
                                line=dict(color='blue', width=2)), row=1, col=1)
        fig.add_trace(go.Scatter(x=zeiten, y=sollwerte, name='Soll-Temperatur', 
                                line=dict(color='red', dash='dash')), row=1, col=1)
        
        # Stellgröße
        fig.add_trace(go.Scatter(x=zeiten, y=stellgroessen, name='Stellgröße', 
                                line=dict(color='green')), row=1, col=2)
        
        # PID-Anteile
        fig.add_trace(go.Scatter(x=zeiten, y=p_anteile, name='P-Anteil', 
                                line=dict(color='orange')), row=2, col=1)
        fig.add_trace(go.Scatter(x=zeiten, y=i_anteile, name='I-Anteil', 
                                line=dict(color='purple')), row=2, col=1)
        fig.add_trace(go.Scatter(x=zeiten, y=d_anteile, name='D-Anteil', 
                                line=dict(color='brown')), row=2, col=1)
        
        # Wärmeverlust
        fig.add_trace(go.Scatter(x=zeiten, y=waermeverluste, name='Wärmeverlust', 
                                line=dict(color='cyan')), row=2, col=2)
        
        # Regelabweichung
        regelabweichung = [soll - ist for soll, ist in zip(sollwerte, temperaturen)]
        fig.add_trace(go.Scatter(x=zeiten, y=regelabweichung, name='Regelabweichung', 
                                line=dict(color='magenta')), row=3, col=1)
        
        # Layout anpassen
        fig.update_xaxes(title_text="Zeit [min]")
        fig.update_yaxes(title_text="Temperatur [°C]", row=1, col=1)
        fig.update_yaxes(title_text="Leistung [W]", row=1, col=2)
        fig.update_yaxes(title_text="PID-Werte", row=2, col=1)
        fig.update_yaxes(title_text="Wärmeverlust [W]", row=2, col=2)
        fig.update_yaxes(title_text="Abweichung [°C]", row=3, col=1)
        
        fig.update_layout(height=1000, showlegend=True, 
                         title_text="Bioreaktor Temperaturregelung - Simulationsergebnisse")
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Kennwerte anzeigen
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Endtemperatur", f"{temperaturen[-1]:.1f}°C")
            
        with col2:
            max_abw = max([abs(x) for x in regelabweichung])
            st.metric("Max. Abweichung", f"{max_abw:.1f}°C")
            
        with col3:
            einschw_zeit = None
            toleranz = 0.5
            for i, abw in enumerate(regelabweichung):
                if abs(abw) <= toleranz:
                    einschw_zeit = zeiten[i]
                    break
            st.metric("Einschwingzeit", f"{einschw_zeit:.1f} min" if einschw_zeit else "Nicht erreicht")
            
        with col4:
            mittlere_leistung = np.mean([abs(x) for x in stellgroessen])
            st.metric("Mittlere Leistung", f"{mittlere_leistung:.0f} W")
        
        # Datenexport
        if st.button("📊 Daten als CSV exportieren"):
            df = pd.DataFrame({
                'Zeit [min]': zeiten,
                'Temperatur [°C]': temperaturen,
                'Sollwert [°C]': sollwerte,
                'Stellgröße [W]': stellgroessen,
                'P-Anteil': p_anteile,
                'I-Anteil': i_anteile,
                'D-Anteil': d_anteile,
                'Wärmeverlust [W]': waermeverluste,
                'Regelabweichung [°C]': regelabweichung
            })
            
            csv = df.to_csv(index=False)
            st.download_button(
                label="CSV herunterladen",
                data=csv,
                file_name='bioreaktor_simulation.csv',
                mime='text/csv'
            )

if __name__ == "__main__":
    main()