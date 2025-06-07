#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren 
### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl
import logging

# Logging für Debugging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen=10, start_temp=20, umg_temp=20, rpm=100, wandmaterial='stahl', wandstaerke=5):
        # Input validation
        if volumen <= 0 or start_temp < -273.15 or rpm < 0 or wandstaerke <= 0:
            raise ValueError("Ungültige Eingabeparameter")
            
        # Speichere die Starttemperatur für reset()
        self.start_temp = start_temp

        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        self._update_fluid_properties(start_temp)

        # Geometrie des Bioreaktors (Zylinder, Annahme: H = 2 * r)
        self.volumen = volumen / 1000  # Volumen in m³ (Liter → m³)
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)  # Radius in m 
        self.hoehe = 2 * self.radius  # Höhe in m 
        
        # Korrekte Flächenberechnung für Wärmeverlust (Außenfläche)
        r_außen = self.radius + wandstaerke/1000
        h_außen = self.hoehe + 2*wandstaerke/1000
        self.flaeche = 2 * np.pi * (r_außen ** 2) + 2 * np.pi * r_außen * h_außen
        
        self.wandstaerke = wandstaerke / 1000  # Wandstärke in m
        self.d_ruehrer = (2 * self.radius) / 3  # Rührerdurchmesser in m

        # Betriebsbedingungen
        self.umg_temp = umg_temp
        self.ist_temp = start_temp
        self.rpm = rpm

        # Wärmeübertragung                                                              
        self.h_intern = self.berechnung_h()
        self.h_extern = self._berechnung_h_extern()  # Verbesserte externe Berechnung
        self.lambda_wand = self.get_waermeleitfaehigkeit(wandmaterial)
        
        # Zusätzliche Eigenschaften für bessere Simulation
        self.max_heizleistung = 5000  # W - Realistische Begrenzung
        self.min_heizleistung = -2000  # W - Kühlleistung
        
    def _update_fluid_properties(self, temp):
        """Aktualisiert temperaturabhängige Fluideigenschaften"""
        try:
            temp_k = temp + 273.15
            self.spez_c = PropsSI("Cpmass", "T", temp_k, "P", 101325, "Water")
            self.dichte = PropsSI("Dmass", "T", temp_k, "P", 101325, "Water")
        except Exception as e:
            logger.warning(f"CoolProp Fehler bei T={temp}°C: {e}")
            # Fallback-Werte für Wasser bei Raumtemperatur
            self.spez_c = 4186  # J/(kg·K)
            self.dichte = 1000  # kg/m³

    def _berechnung_h_extern(self):
        """Verbesserte Berechnung des externen Wärmeübergangskoeffizienten"""
        # Natürliche Konvektion an senkrechter Wand (vereinfacht)
        delta_t = abs(self.ist_temp - self.umg_temp)
        if delta_t > 0:
            # Vereinfachte Korrelation für natürliche Konvektion
            h_ext = 5.0 + 3.0 * (delta_t ** 0.25)
        else:
            h_ext = 5.0  # Minimum für ruhende Luft
        return min(h_ext, 25.0)  # Begrenzung auf realistische Werte

    def berechnung_h(self):
        """Berechnung des internen Wärmeübergangskoeffizienten"""
        try:
            t_k = self.ist_temp + 273.15
            
            k = PropsSI("CONDUCTIVITY", "T", t_k, "P", 101325, "Water")
            mu = PropsSI("VISCOSITY", "T", t_k, "P", 101325, "Water")
            
            Re = ((self.rpm / 60) * (self.d_ruehrer ** 2) * self.dichte) / mu
            Pr = Prandtl(self.spez_c, mu, k)
            
            # Verbesserte Korrelationen mit Gültigkeitsbereichen
            if Re < 10:
                Nu = 2.0  # Minimum für ruhende Flüssigkeit
            elif 10 <= Re < 4.5e3:
                Nu = 0.664 * (Re ** 0.5) * (Pr ** (1/3))  # Laminare Strömung
            elif 4.5e3 <= Re < 1e4 and 0.6 < Pr < 160:
                Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)
            elif Re >= 1e4 and Pr >= 0.6:
                Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)
            else:
                Nu = 3.66
                
            h = Nu * k / self.d_ruehrer
            return max(h, 100)  # Mindest-Wärmeübergangskoeffizient
            
        except Exception as e:
            logger.warning(f"Fehler bei h-Berechnung: {e}")
            return 500  # Fallback-Wert

    def get_waermeleitfaehigkeit(self, material):
        """Wärmeleitfähigkeit für Wandmaterialien"""
        material = str(material).lower()
        material_db = {
            'stahl': 21.0,
            'edelstahl': 16.0,
            'glas': 1.4,
            'kunststoff': 0.2,
            'aluminium': 237.0
        }
        return material_db.get(material, 16.0)

    def berechnung_waermeverlust(self):
        """Berechnung des Wärmeverlustes mit verbesserter Fehlerbehandlung"""
        delta_t = self.ist_temp - self.umg_temp
        
        if abs(delta_t) < 0.01:  # Vermeidung von Division durch sehr kleine Zahlen
            return 0.0
        
        A = self.flaeche
        d = self.wandstaerke
        λ = self.lambda_wand
        h_int = self.h_intern
        h_ext = self.h_extern
        
        # Wärmewiderstände
        R_int = 1.0 / (h_int * A)
        R_cond = d / (λ * A)
        R_ext = 1.0 / (h_ext * A)
        
        R_gesamt = R_int + R_cond + R_ext
        q_verlust = delta_t / R_gesamt
        
        return q_verlust
    
    def update_temperature(self, leistung, zeitintervall=1):
        """Aktualisiert die Reaktortemperatur"""
        # Leistungsbegrenzung
        leistung = np.clip(leistung, self.min_heizleistung, self.max_heizleistung)
        
        # Eigenschaften aktualisieren
        self._update_fluid_properties(self.ist_temp)
        self.h_intern = self.berechnung_h()
        self.h_extern = self._berechnung_h_extern()
        
        # Energiebilanz
        masse = self.dichte * self.volumen
        waermeverlust = self.berechnung_waermeverlust()
        energie_netto = (leistung - waermeverlust) * zeitintervall
        temp_delta = energie_netto / (masse * self.spez_c)
        
        self.ist_temp += temp_delta
        
        # Physikalische Grenzen
        self.ist_temp = max(self.ist_temp, -50)  # Praktische untere Grenze
        self.ist_temp = min(self.ist_temp, 150)  # Praktische obere Grenze
        
        return self.ist_temp
    
    def reset(self, start_temp=None):
        """Reset des Reaktors"""
        if start_temp is None:
            self.ist_temp = self.start_temp
        else:
            self.ist_temp = start_temp
        self._update_fluid_properties(self.ist_temp)

class PID:
    """Verbesserter PID-Regler mit robustem Anti-Windup"""
    
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, dt=1.0, output_min=-1000, output_max=1000):
        if dt <= 0:
            raise ValueError("Zeitschritt dt muss positiv sein")
            
        self.dt = dt
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.output_min = output_min
        self.output_max = output_max
        
        # Interne Zustandsvariablen
        self.fehler_vor = 0.0
        self.integral = 0.0
        self.output_vor = 0.0
        
        # Anti-Windup Parameter
        self.integral_max = abs(output_max - output_min) / max(ki, 1e-6)
        
    def run(self, sollwert, istwert):
        """Berechnet PID-Ausgabe mit verbessertem Anti-Windup"""
        fehler = sollwert - istwert
        
        # P-Anteil
        p_anteil = self.kp * fehler
        
        # I-Anteil mit Anti-Windup
        self.integral += fehler * self.dt
        # Begrenzung des Integrals
        self.integral = np.clip(self.integral, -self.integral_max, self.integral_max)
        i_anteil = self.ki * self.integral
        
        # D-Anteil (jetzt korrekt implementiert)
        d_anteil = self.kd * (fehler - self.fehler_vor) / self.dt
        
        # Gesamtausgabe
        stellgroesse = p_anteil + i_anteil + d_anteil
        
        # Ausgangsbegrenzung mit Anti-Windup
        if stellgroesse > self.output_max:
            output = self.output_max
            # Integral-Korrektur bei Sättigung
            excess = stellgroesse - self.output_max
            self.integral -= excess / max(self.ki, 1e-6)
        elif stellgroesse < self.output_min:
            output = self.output_min
            # Integral-Korrektur bei Sättigung
            excess = self.output_min - stellgroesse
            self.integral += excess / max(self.ki, 1e-6)
        else:
            output = stellgroesse
        
        # Zustand speichern
        self.fehler_vor = fehler
        self.output_vor = output
        
        return output
    
    def reset(self):
        """Regler zurücksetzen"""
        self.fehler_vor = 0.0
        self.integral = 0.0
        self.output_vor = 0.0
    
    def get_status(self):
        """Debug-Information für den Regler"""
        return {
            'integral': self.integral,
            'fehler_vor': self.fehler_vor,
            'output_vor': self.output_vor
        }

# Streamlit-Anwendung
st.set_page_config(page_title="Bioreaktor Temperaturregelung", layout="wide")
st.title("🧪 Simulation einer Temperaturregelung von Bioreaktoren")

# Info-Box am Anfang
st.info("🔄 **Live-Simulation:** Die Ergebnisse werden automatisch bei Parameteränderungen aktualisiert!")

# Sidebar für Einstellungen
st.sidebar.title("⚙️ Simulationseinstellungen")
st.sidebar.markdown("""
**Bioreaktor-Simulation mit PID-Regelung**

Passen Sie die Parameter an und beobachten Sie die Auswirkungen auf die Temperaturregelung.
""")

with st.sidebar:
    st.header("🔬 Reaktorparameter")    
    volumen = st.slider("Reaktorvolumen (L)", 1, 100, 10, help="Größeres Volumen = träger")
    start_temp = st.slider("Starttemperatur (°C)", 5, 40, 20)
    umg_temp = st.slider("Umgebungstemperatur (°C)", 5, 40, 20)
    wandmaterial = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Edelstahl", "Aluminium"])
    wandstaerke = st.slider("Wandstärke (mm)", 1, 20, 5)
    rpm = st.slider("Rührerdrehzahl (1/min)", 50, 500, 100, help="Höhere Drehzahl = bessere Wärmeübertragung")
    
    st.header("🎯 Sollwert & Simulation")
    soll_temp = st.slider("Solltemperatur (°C)", 25, 80, 37)
    simdauer = st.slider("Simulationsdauer (min)", 10, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 30, 5)
    
    st.header("🎛️ PID-Parameter")
    col1, col2 = st.columns(2)
    with col1:
        kp = st.slider("Kp (P-Anteil)", 0.0, 200.0, 50.0, 1.0, help="Höher = schneller, aber instabiler")
        ki = st.slider("Ki (I-Anteil)", 0.0, 10.0, 1.0, 0.1, help="Eliminiert bleibende Regelabweichung")
    with col2:
        kd = st.slider("Kd (D-Anteil)", 0.0, 50.0, 5.0, 0.5, help="Dämpft Schwingungen")
        
    st.header("🔧 Erweiterte Einstellungen")
    show_details = st.checkbox("Detaillierte Analyse anzeigen", False)
    add_disturbance = st.checkbox("Störung hinzufügen", False)
    if add_disturbance:
        störung_zeit = st.slider("Störungszeitpunkt (min)", 5, simdauer-5, simdauer//2)
        störung_größe = st.slider("Störungsgröße (°C)", -10, 10, -5)

# Hauptbereich mit Tabs
tab1, tab2, tab3 = st.tabs(["📊 Simulation", "📋 Analyse", "ℹ️ Info"])

with tab1:
    col1, col2 = st.columns([2, 1])
    
    with col2:
        st.subheader("Reaktor-Eigenschaften")
        
        # Erstelle Reaktor für Eigenschaften-Anzeige
        reaktor_info = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial.lower(), wandstaerke)
        
        st.metric("Reaktorvolumen", f"{volumen} L")
        st.metric("Reaktorradius", f"{reaktor_info.radius*100:.1f} cm")
        st.metric("Reaktorhöhe", f"{reaktor_info.hoehe*100:.1f} cm")
        st.metric("Wärmeübertr.-Fläche", f"{reaktor_info.flaeche:.2f} m²")
        st.metric("Wärmeleitfähigkeit", f"{reaktor_info.lambda_wand} W/(m·K)")
    
    with col1:
        # Automatische Simulation (läuft bei jeder Parameter-Änderung)
        try:
            # Initialisierung
            reaktor_pid = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial.lower(), wandstaerke)
            pid = PID(kp=kp, ki=ki, kd=kd, dt=dt, output_min=-2000, output_max=5000)
            reaktor_ungeregelt = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial.lower(), wandstaerke)
            
            # Simulation
            n_steps = int(simdauer * 60 // dt)
            zeiten = np.arange(0, n_steps * dt, dt) / 60
            
            temps_pid, temps_offen, leistungen = [], [], []
            sollwerte = []
            
            reaktor_pid.reset(start_temp)
            pid.reset()
            reaktor_ungeregelt.reset(start_temp)
            
            # Simulationsschleife
            for i in range(n_steps):
                t_min = zeiten[i]
                current_soll = soll_temp
                
                # Störung hinzufügen
                if add_disturbance and abs(t_min - störung_zeit) < dt/120:
                    reaktor_pid.ist_temp += störung_größe
                    reaktor_ungeregelt.ist_temp += störung_größe
                
                # PID-geregeltes System
                leistung = pid.run(current_soll, reaktor_pid.ist_temp)
                temp_pid = reaktor_pid.update_temperature(leistung, dt)
                
                # Ungeregeltes System
                temp_offen = reaktor_ungeregelt.update_temperature(0, dt)
                
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
            
            if add_disturbance:
                ax1.axvline(störung_zeit, color="red", linestyle="--", alpha=0.7, label="Störung")
            
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

with tab3:
    st.subheader("ℹ️ Informationen zur Simulation")
    
    st.markdown("""
    ### 🧪 Bioreaktor-Modell
    
    Das Modell simuliert einen zylindrischen Bioreaktor mit:
    - Realistischen thermodynamischen Eigenschaften (CoolProp)
    - Wärmeübertragung durch Konvektion und Leitung
    - Temperaturabhängigen Fluideigenschaften
    - Rührwerk für bessere Durchmischung
    
    ### 🎛️ PID-Regelung
    
    **P-Anteil (Kp):** Proportional zum aktuellen Fehler
    - Höhere Werte: Schnellere Reaktion, aber möglicherweise instabil
    
    **I-Anteil (Ki):** Eliminiert bleibende Regelabweichung
    - Zu hoch: Langsam und schwingend
    - Zu niedrig: Dauerhafte Abweichung
    
    **D-Anteil (Kd):** Dämpft Änderungen und Schwingungen
    - Verbessert Stabilität bei schnellen Änderungen
    
    ### 📊 Tipps für gute Regelung
    
    1. **Start konservativ:** Beginnen Sie mit niedrigen Werten
    2. **P-Anteil zuerst:** Erhöhen Sie Kp bis leichte Schwingungen auftreten
    3. **I-Anteil hinzufügen:** Reduziert bleibende Abweichung
    4. **D-Anteil feintunen:** Dämpft Schwingungen
    
    ### ⚠️ Physikalische Grenzen
    
    - Maximale Heizleistung: 5000 W
    - Maximale Kühlleistung: -2000 W
    - Temperaturbereich: -50°C bis 150°C
    """)
    
    with st.expander("🔧 Technische Details"):
        st.markdown("""
        **Wärmeübertragungskoeffizienten:**
        - Intern: Berechnet über Nusselt-Korrelationen
        - Extern: Natürliche Konvektion an Luft
        
        **Verwendete Korrelationen:**
        - Laminare Strömung: Nu = 0.664 × Re^0.5 × Pr^(1/3)
        - Turbulente Strömung: Nu = 0.023 × Re^0.8 × Pr^0.4
        - Rührkessel: Nu = 0.354 × Re^0.714 × Pr^0.260
        
        **Materialeigenschaften:**
        - Stahl: λ = 21 W/(m·K)
        - Glas: λ = 1.4 W/(m·K)
        - Edelstahl: λ = 16 W/(m·K)
        - Aluminium: λ = 237 W/(m·K)
        """)

st.sidebar.markdown("---")
st.sidebar.caption("🔄 **Live-Modus:** Ergebnisse werden automatisch aktualisiert!")
st.sidebar.caption("💡 **Tipp:** Experimentieren Sie mit verschiedenen Parametern!")