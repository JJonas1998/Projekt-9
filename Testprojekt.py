#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren 

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

### Datum: 13.06.2025

# Importiere benötigte Bibliotheken
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
class Bioreaktor:
    """Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen, Kühle und Umgebungseinflüsse verändert.
    """
    def __init__(self, 
        volumen = 10, 
        start_temp = 20, 
        umg_temp = 20, 
        rpm = 100, 
        wandmaterial = 'stahl', 
        wandstaerke = 5):
            
        # Speichere die Starttemperatur für reset()
        self.start_temp = start_temp

        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        self.update_stoffwerte(start_temp)

        # Geometrie des Bioreaktors (Zylinder, Annahme: H = 3 * r)
        self.volumen = volumen / 1000  # Volumen in m³ (Liter → m³)
        self.wandstaerke = wandstaerke / 1000  # Wandstärke in m (mm → m)
        self.r_i = (self.volumen / (2 * np.pi)) ** (1/3)  # Innenradius in m 
        self.r_a = self.r_i + self.wandstaerke  # Außenradius in m
        self.h_i = 3 * self.r_i  # Innenhöhe in m 
        self.h_a = self.h_i + 2 * self.wandstaerke  # Außenhöhe in m
        self.flaeche_i = 2 * np.pi * (self.r_i ** 2) + 2 * np.pi * self.r_i * self.h_i # Innenfläche in m²
        self.flaeche_a = 2 * np.pi * (self.r_a ** 2) + 2 * np.pi * self.r_a * self.h_a  # Außenfläche in m²
        self.d_ruehrer = (2 * self.r_i) / 3  # Rührerdurchmesser in m

        # Betriebsbedingungen
        self.umg_temp = umg_temp  # Umgebungstemperatur in °C
        self.ist_temp = start_temp  # Aktuelle Temperatur in °C
        self.n = rpm / 60  # Rührerdrehzahl in 1/s (1/min → 1/s)

        # Wärmeübertragung                                                              
        self.h_intern = self.berechnung_h_intern()  # Interner Wärmeübergangskoeffizient
        self.h_extern = 35  # Externer Wärmeübergangskoeffizient (natürliche Konvektion an Luft in W/(m²·K))
        self.lambda_wand = self.get_lambda(wandmaterial)  # Wärmeleitfähigkeit des Wandmaterials in W/(m·K)
        
        # Physikalische Grenzen
        self.max_leistung = 5000  # maximale Heizleistung in W
        self.min_leistung = -2000  # maximale Kühlleistung in W
        
    def update_stoffwerte(self, temp):
        """Aktualisiert die thermophysikalischen Eigenschaften des Mediums (Wasser).
        """
        try:
            temp_k = temp + 273.15
            self.spez_c = PropsSI("Cpmass", "T", temp_k, "P", 101325, "Water")  # Spezifische Wärmekapazität in J/(kg·K)
            self.dichte = PropsSI("Dmass", "T", temp_k, "P", 101325, "Water")  # Dichte in kg/m³
            self.k = PropsSI("CONDUCTIVITY", "T", temp_k, "P", 101325, "Water")  # Wärmeleitfähigkeit in W/(m·K)
            self.mu = PropsSI("VISCOSITY", "T", temp_k, "P", 101325, "Water")  # Dynamische Viskosität in Pa·s
        except Exception:
            self.spez_c = 4186  
            self.dichte = 1000
            self.k = 0.606
            self.mu = 0.001002

    def berechnung_h_intern(self):
        """Berechnung des Wärmeübergangskoeffizienten h anhand von Impeller-Drehzahl (rpm) und Rührerdurchmesser.
        """
        try:
            # Dimensionslose Kennzahlen
            Re = self.n * (self.d_ruehrer ** 2) * self.dichte / self.mu  # Reynolds-Zahl
            Pr = Prandtl(self.spez_c, self.mu, self.k)  # Prandtl-Zahl
            
            # Nusselt-Zahl basierend auf Reynolds- und Prandtl-Zahl
            if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:  
                Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)  # Übergangsbereich (4.5e3 < Re < 1e4, 0.6 < Pr < 160)                                                                              
            elif Re >= 1e4 and Pr >= 0.6:  
                Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)  # Turbulente Strömung (Re >= 10000)                                                                                                
            else:  
                Nu = 3.66  # Laminare Strömung (Re < 2300)
                
            # Berechnung des Wärmeübergangskoeffizienten h
            h = Nu * self.k / self.d_ruehrer  # Wärmeübergangskoeffizient in W/(m²·K)

            return max(h, 150)  # Mindest-Wärmeübergangskoeffizient
        except Exception:
            return 500  # Fallback-Wert bei Fehlern

    def get_lambda(self, material):
        """Wärmeleitfähigkeit für Wandmaterialien.
        """
        material = str(material).lower()

        # Datenbank für Wärmeleitfähigkeiten in W/(m·K)
        material_db = {
            'stahl': 46.0,  
            'edelstahl': 21.0,  
            'glas': 1.4,
            'kunststoff': 0.3,
            'aluminium': 230.0
        }
        return material_db.get(material, 46.0)

    def t_verlust(self):
        """Berechnet den Wärmeverlust des Reaktors basierend auf Temperaturdifferenz und Wärmeübergangskoeffizienten.
        """
        delta_t = self.ist_temp - self.umg_temp  # Temperaturdifferenz in °C
        
        # Vermeidung von Division durch Null
        if abs(delta_t) < 0.01:  
            return 0.0
        
        # Berechnung der Wärmeverluste
        R_int = 1.0 / (self.h_intern * self.flaeche_i)  # Interner Widerstand in K/W
        R_cond = self.wandstaerke / (self.lambda_wand * self.flaeche_i)  # Widerstand der Wand in K/W
        R_ext = 1.0 / (self.h_extern * self.flaeche_a)  # Externer Widerstand in K/W
        
        # Gesamtwiderstand
        R_gesamt = R_int + R_cond + R_ext  
        q_verlust = delta_t / R_gesamt  # Wärmeverlust in W
        
        return q_verlust
    
    def update_temperatur(self, 
            leistung, 
            zeitintervall = 1):
        """Aktualisiert die Reaktortemperatur basierend auf der Heizleistung, Kühlleistung, Wärmeverlust und Zeitintervall.
        """
        # Leistungsbegrenzung
        leistung = np.clip(leistung, self.min_leistung, self.max_leistung)  
        
        # Eigenschaften aktualisieren
        self.update_stoffwerte(self.ist_temp)
        self.h_intern = self.berechnung_h_intern()
        
        # Energiebilanz
        masse = self.dichte * self.volumen  # Masse des Mediums in kg
        energie_netto = (leistung - self.t_verlust()) * zeitintervall  # Nettoenergie in J
        temp_delta = energie_netto / (masse * self.spez_c)  # Temperaturänderung in K
        
        self.ist_temp += temp_delta  
        
        # Physikalische Grenzen
        self.ist_temp = max(self.ist_temp, 4)  # Praktische untere Grenze
        self.ist_temp = min(self.ist_temp, 100)  # Praktische obere Grenze
        
        return self.ist_temp
    
    def reset(self, 
        start_temp = None):
        """"Setzt die Ist-Temperatur auf die Starttemperatur zurück und aktualisiert die Stoffwerte.
        """
        if start_temp is None:
            self.ist_temp = self.start_temp
        else:
            self.ist_temp = start_temp
        
        # Aktualisiere Stoffwerte und interne Wärmeübergangskoeffizienten
        self.update_stoffwerte(self.ist_temp)

## PID-Regler-Klasse
class PID:
    """PID-Regler zur Regelung der Temperatur eines Bioreaktors.
    """
    def __init__(self, 
        kp = 0.0, 
        ki = 0.0, 
        kd = 0.0, 
        dt = 1.0, 
        output_min = -1000, 
        output_max = 1000):
    
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
        
    def run(self, sollwert, istwert):
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
        """Regler zurücksetzen
        """
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
    volumen = st.slider("Reaktorvolumen (L)", 1, 100, 10)
    start_temp = st.slider("Starttemperatur (°C)", 5, 40, 20)
    umg_temp = st.slider("Umgebungstemperatur (°C)", 5, 40, 20)
    wandmaterial = st.selectbox("Wandmaterial", ["Stahl", "Glas", "Edelstahl", "Aluminium"])
    wandstaerke = st.slider("Wandstärke (mm)", 1, 50, 5)
    rpm = st.slider("Rührerdrehzahl (1/min)", 50, 180, 100)
    
    st.header("Sollwert & Simulation")
    soll_temp = st.slider("Solltemperatur (°C)", 25, 80, 37)
    simdauer = st.slider("Simulationsdauer (min)", 10, 180, 60)
    dt = st.slider("Zeitschritt (s)", 1, 30, 5)
    
    st.header("PID-Parameter")
    kp = st.slider("Kp (P-Anteil)", 0.0, 200.0, 50.0, 1.0, help = "Höher = schneller, aber instabiler")
    ki = st.slider("Ki (I-Anteil)", 0.0, 10.0, 1.0, 0.1, help = "Eliminiert bleibende Regelabweichung")
    kd = st.slider("Kd (D-Anteil)", 0.0, 50.0, 5.0, 0.5, help = "Dämpft Schwingungen")
    
# Hauptbereich mit Tabs
tab1, tab2 = st.tabs(["📊 Simulation", "📋 Analyse"])

with tab1:
    col1, col2 = st.columns([2, 1])
    
    with col2:
        st.subheader("Reaktor-Eigenschaften")
        
        # Erstelle Reaktor für Eigenschaften-Anzeige
        reaktor_info = Bioreaktor(volumen, start_temp, umg_temp, rpm, wandmaterial.lower(), wandstaerke)
        
        st.metric("Reaktorvolumen", f"{volumen} L")
        st.metric("Reaktorradius", f"{reaktor_info.r_i*100:.1f} cm")
        st.metric("Reaktorhöhe", f"{reaktor_info.h_i*100:.1f} cm")
        st.metric("Wärmeübertr.-Fläche", f"{reaktor_info.flaeche_i:.2f} m²")
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
                
                # PID-geregeltes System
                leistung = pid.run(current_soll, reaktor_pid.ist_temp)
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

