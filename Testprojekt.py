import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

class Bioreactor:
    """Einfaches Bioreaktor-Modell mit Wärmeübertragung"""
    
    def __init__(self, volume=100, heat_capacity=4180, density=1000, 
                 heat_transfer_coeff=500, ambient_temp=20):
        self.volume = volume  # Liter
        self.heat_capacity = heat_capacity  # J/(kg*K) - Wasser
        self.density = density  # kg/m³
        self.heat_transfer_coeff = heat_transfer_coeff  # W/(m²*K)
        self.surface_area = 2 * np.pi * (volume/1000)**(2/3)  # m² (geschätzt)
        self.ambient_temp = ambient_temp  # °C
        self.temperature = ambient_temp  # Starttemperatur
        self.biological_heat = 0  # Wärme durch biologische Prozesse
        
    def update_temperature(self, heating_power, dt=1.0):
        """Temperaturänderung basierend auf Energiebilanz"""
        # Masse des Reaktorinhalts
        mass = self.volume * self.density / 1000  # kg
        
        # Wärmeverlust an Umgebung
        heat_loss = self.heat_transfer_coeff * self.surface_area * \
                   (self.temperature - self.ambient_temp)  # W
        
        # Biologische Wärmeproduktion (vereinfacht)
        bio_heat = self.biological_heat * mass  # W
        
        # Netto-Wärmezufuhr
        net_heat = heating_power + bio_heat - heat_loss  # W
        
        # Temperaturänderung
        dT = (net_heat * dt) / (mass * self.heat_capacity)  # K/s * s = K
        self.temperature += dT
        
        return self.temperature

class PIDController:
    """PID-Regler für Temperaturkontrolle"""
    
    def __init__(self, kp=1000, ki=10, kd=50):
        self.kp = kp  # Proportionalverstärkung
        self.ki = ki  # Integralverstärkung  
        self.kd = kd  # Differentialverstärkung
        self.setpoint = 37.0  # Zieltemperatur °C
        self.integral = 0
        self.last_error = 0
        self.dt = 1.0  # Zeitschritt
        
    def update(self, current_temp):
        """PID-Regelung berechnen"""
        error = self.setpoint - current_temp
        
        # Proportional-Anteil
        p_term = self.kp * error
        
        # Integral-Anteil
        self.integral += error * self.dt
        i_term = self.ki * self.integral
        
        # Differential-Anteil
        d_term = self.kd * (error - self.last_error) / self.dt
        self.last_error = error
        
        # PID-Output (Heizleistung in W)
        output = p_term + i_term + d_term
        
        # Begrenzung der Heizleistung (0-5000W)
        output = max(0, min(5000, output))
        
        return output, p_term, i_term, d_term

def simulate_system(reactor, controller, duration, with_disturbance=False):
    """Systemsimulation durchführen"""
    times = []
    temperatures = []
    setpoints = []
    heating_powers = []
    p_terms, i_terms, d_terms = [], [], []
    
    for t in range(duration):
        # Störung nach der Hälfte der Zeit
        if with_disturbance and t == duration // 2:
            reactor.biological_heat = 20  # W/kg zusätzliche Wärme
            
        # PID-Regelung
        heating_power, p, i, d = controller.update(reactor.temperature)
        
        # Reaktor-Update
        reactor.update_temperature(heating_power)
        
        # Daten sammeln
        times.append(t)
        temperatures.append(reactor.temperature)
        setpoints.append(controller.setpoint)
        heating_powers.append(heating_power)
        p_terms.append(p)
        i_terms.append(i)
        d_terms.append(d)
    
    return pd.DataFrame({
        'Zeit': times,
        'Temperatur': temperatures,
        'Sollwert': setpoints,
        'Heizleistung': heating_powers,
        'P-Anteil': p_terms,
        'I-Anteil': i_terms,
        'D-Anteil': d_terms
    })

def main():
    st.title("🧪 Bioreaktor Temperaturregelung Simulation")
    st.markdown("Simulation eines PID-geregelten Bioreaktors mit Rückkopplungsschleifen")
    
    # Sidebar für Parameter
    st.sidebar.header("📊 Simulationsparameter")
    
    # Reaktorparameter
    st.sidebar.subheader("Reaktor")
    volume = st.sidebar.slider("Volumen (L)", 50, 500, 100)
    ambient_temp = st.sidebar.slider("Umgebungstemperatur (°C)", 15, 30, 20)
    
    # PID-Parameter
    st.sidebar.subheader("PID-Regler")
    kp = st.sidebar.slider("Kp (Proportional)", 100, 2000, 1000)
    ki = st.sidebar.slider("Ki (Integral)", 1, 50, 10)
    kd = st.sidebar.slider("Kd (Differential)", 10, 100, 50)
    setpoint = st.sidebar.slider("Solltemperatur (°C)", 30, 45, 37)
    
    # Simulationseinstellungen
    st.sidebar.subheader("Simulation")
    duration = st.sidebar.slider("Dauer (Minuten)", 30, 180, 60)
    with_disturbance = st.sidebar.checkbox("Störung einschalten", False)
    
    # Simulation starten
    if st.button("🚀 Simulation starten"):
        with st.spinner("Simulation läuft..."):
            # Objekte erstellen
            reactor = Bioreactor(volume=volume, ambient_temp=ambient_temp)
            controller = PIDController(kp=kp, ki=ki, kd=kd)
            controller.setpoint = setpoint
            
            # Simulation durchführen
            data = simulate_system(reactor, controller, duration, with_disturbance)
            
            # Ergebnisse anzeigen
            st.success("✅ Simulation abgeschlossen!")
            
            # Hauptdiagramm - Temperaturverlauf
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
            
            # Temperatur
            ax1.plot(data['Zeit'], data['Temperatur'], 'b-', linewidth=2, label='Ist-Temperatur')
            ax1.plot(data['Zeit'], data['Sollwert'], 'r--', linewidth=2, label='Sollwert')
            ax1.set_ylabel('Temperatur (°C)')
            ax1.set_title('Temperaturregelung')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Heizleistung
            ax2.plot(data['Zeit'], data['Heizleistung'], 'g-', linewidth=2)
            ax2.set_ylabel('Heizleistung (W)')
            ax2.set_title('Stellgröße')
            ax2.grid(True, alpha=0.3)
            
            # PID-Anteile
            ax3.plot(data['Zeit'], data['P-Anteil'], 'r-', alpha=0.7, label='P-Anteil')
            ax3.plot(data['Zeit'], data['I-Anteil'], 'g-', alpha=0.7, label='I-Anteil') 
            ax3.plot(data['Zeit'], data['D-Anteil'], 'b-', alpha=0.7, label='D-Anteil')
            ax3.set_ylabel('PID-Anteile (W)')
            ax3.set_xlabel('Zeit (min)')
            ax3.set_title('PID-Komponenten')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            plt.tight_layout()
            st.pyplot(fig)
            
            # Statistiken
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("🎯 Endtemperatur", f"{data['Temperatur'].iloc[-1]:.1f}°C")
                
            with col2:
                steady_state_error = abs(data['Temperatur'].iloc[-10:].mean() - setpoint)
                st.metric("📊 Bleibende Regelabweichung", f"{steady_state_error:.2f}°C")
                
            with col3:
                max_overshoot = max(data['Temperatur']) - setpoint
                st.metric("📈 Max. Überschwingen", f"{max_overshoot:.2f}°C")
            
            # Datendownload
            st.subheader("📥 Simulationsdaten")
            csv = data.to_csv(index=False)
            st.download_button(
                label="CSV herunterladen",
                data=csv,
                file_name="bioreactor_simulation.csv",
                mime="text/csv"
            )
            
            # Datenvorschau
            with st.expander("📋 Datenvorschau"):
                st.dataframe(data.head(10))
    
    # Theoretischer Hintergrund
    with st.expander("📚 Theoretischer Hintergrund"):
        st.markdown("""
        ### Bioreaktor-Modell
        Das vereinfachte Modell basiert auf der Energiebilanz:
        
        **dT/dt = (Q_heiz + Q_bio - Q_verlust) / (m × cp)**
        
        - Q_heiz: Heizleistung (Stellgröße)
        - Q_bio: Biologische Wärmeproduktion  
        - Q_verlust: Wärmeverlust an Umgebung
        - m: Masse des Reaktorinhalts
        - cp: Wärmekapazität
        
        ### PID-Regler
        **u(t) = Kp×e(t) + Ki×∫e(t)dt + Kd×de(t)/dt**
        
        - Kp: Proportionalverstärkung (schnelle Reaktion)
        - Ki: Integralverstärkung (eliminiert bleibende Abweichung)
        - Kd: Differentialverstärkung (dämpft Schwingungen)
        """)

if __name__ == "__main__":
    main()