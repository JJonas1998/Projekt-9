import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Streamlit Konfiguration
st.set_page_config(page_title="Bioreaktor Temperaturregelung", page_icon="🧪", layout="wide")

class Bioreaktor:
    def __init__(self, volumen=100, start_temp=20):
        self.volumen = volumen / 1000  # L zu m³
        self.spez_c = 4186  # J/(kg*K) für Wasser
        self.dichte = 1000  # kg/m³ für Wasser
        self.flaeche = 2 * (volumen/1000)**(2/3)  # Vereinfachte Oberflächenberechnung
        self.waerme_ueber_h = 50  # W/(m²*K) - vereinfacht
        self.umg_temp = 20
        self.current_temp = start_temp
        
    def update_temperature(self, heizleistung, stoerung=0, dt=1):
        masse = self.volumen * self.dichte
        waermekapazitaet = masse * self.spez_c
        waermeverlust = self.waerme_ueber_h * self.flaeche * (self.current_temp - self.umg_temp)
        waermeaenderung = heizleistung - waermeverlust + stoerung
        dT = (waermeaenderung * dt) / waermekapazitaet
        self.current_temp += dT
        return self.current_temp

class PIDController:
    def __init__(self, kp=1.0, ki=0.1, kd=0.05):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.integral = 0.0
        self.last_error = 0.0
        self.max_output = 5000
        
    def calculate(self, setpoint, current_value, dt=1):
        error = setpoint - current_value
        self.integral += error * dt
        # Anti-Windup
        if self.integral > 1000:
            self.integral = 1000
        elif self.integral < -1000:
            self.integral = -1000
            
        proportional = self.kp * error
        integral = self.ki * self.integral
        derivative = self.kd * (error - self.last_error) / dt if dt > 0 else 0
        
        output = proportional + integral + derivative
        if output > self.max_output:
            output = self.max_output
        elif output < 0:
            output = 0
            
        self.last_error = error
        return output, proportional, integral, derivative

def run_simulation(volumen, start_temp, ambient_temp, target_temp, sim_time_min, dt, kp, ki, kd, 
                  disturbance_active, disturbance_time, disturbance_power, disturbance_duration):
    
    sim_time = sim_time_min * 60
    time_array = np.arange(0, sim_time + dt, dt)
    
    # Initialisierung
    bioreaktor = Bioreaktor(volumen, start_temp)
    bioreaktor.umg_temp = ambient_temp
    pid = PIDController(kp, ki, kd)
    
    # Ergebnis-Arrays
    temperature = []
    power = []
    error = []
    p_term = []
    i_term = []
    d_term = []
    
    for t in time_array:
        # Störung prüfen
        disturbance = 0
        if disturbance_active:
            disturbance_start = disturbance_time * 60
            disturbance_end = disturbance_start + disturbance_duration * 60
            if disturbance_start <= t <= disturbance_end:
                disturbance = disturbance_power
        
        # PID-Regelung
        current_error = target_temp - bioreaktor.current_temp
        current_power, p, i, d = pid.calculate(target_temp, bioreaktor.current_temp, dt)
        new_temp = bioreaktor.update_temperature(current_power, disturbance, dt)
        
        # Ergebnisse speichern
        temperature.append(new_temp)
        power.append(current_power)
        error.append(current_error)
        p_term.append(p)
        i_term.append(i)
        d_term.append(d)
    
    return {
        'time': time_array,
        'temperature': np.array(temperature),
        'power': np.array(power),
        'error': np.array(error),
        'p_term': np.array(p_term),
        'i_term': np.array(i_term),
        'd_term': np.array(d_term)
    }

def plot_results(results, target_temp):
    time_min = results['time'] / 60
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Temperaturverlauf', 'Heizleistung', 'PID-Komponenten', 'Regelfehler')
    )
    
    # Temperatur
    fig.add_trace(
        go.Scatter(x=time_min, y=results['temperature'], 
                  name='Temperatur', line=dict(color='blue', width=2)),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(x=time_min, y=[target_temp]*len(time_min), 
                  name='Zieltemperatur', line=dict(color='red', dash='dash')),
        row=1, col=1
    )
    
    # Heizleistung
    fig.add_trace(
        go.Scatter(x=time_min, y=results['power'], 
                  name='Heizleistung', line=dict(color='orange')),
        row=1, col=2
    )
    
    # PID-Komponenten
    fig.add_trace(
        go.Scatter(x=time_min, y=results['p_term'], 
                  name='P-Anteil', line=dict(color='green')),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(x=time_min, y=results['i_term'], 
                  name='I-Anteil', line=dict(color='purple')),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(x=time_min, y=results['d_term'], 
                  name='D-Anteil', line=dict(color='brown')),
        row=2, col=1
    )
    
    # Regelfehler
    fig.add_trace(
        go.Scatter(x=time_min, y=results['error'], 
                  name='Regelfehler', line=dict(color='red')),
        row=2, col=2
    )
    
    fig.update_layout(
        height=600, 
        title_text="Bioreaktor Temperaturregelung - Simulationsergebnisse",
        showlegend=True
    )
    
    # Achsenbeschriftungen
    fig.update_xaxes(title_text="Zeit [min]", row=1, col=1)
    fig.update_xaxes(title_text="Zeit [min]", row=1, col=2)
    fig.update_xaxes(title_text="Zeit [min]", row=2, col=1)
    fig.update_xaxes(title_text="Zeit [min]", row=2, col=2)
    
    fig.update_yaxes(title_text="Temperatur [°C]", row=1, col=1)
    fig.update_yaxes(title_text="Leistung [W]", row=1, col=2)
    fig.update_yaxes(title_text="PID-Terme [W]", row=2, col=1)
    fig.update_yaxes(title_text="Fehler [°C]", row=2, col=2)
    
    st.plotly_chart(fig, use_container_width=True)

def analyze_results(results, target_temp):
    temp = results['temperature']
    error = results['error']
    
    # Metriken berechnen
    steady_start = int(0.8 * len(temp))
    steady_error = np.mean(np.abs(error[steady_start:]))
    steady_std = np.std(temp[steady_start:])
    max_temp = np.max(temp)
    overshoot = max(0, (max_temp - target_temp) / target_temp * 100) if target_temp > 0 else 0
    avg_power = np.mean(results['power'])
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Steady-State Fehler", f"{steady_error:.2f} °C")
    with col2:
        st.metric("Standardabweichung", f"{steady_std:.3f} °C")
    with col3:
        st.metric("Überschwingen", f"{overshoot:.1f} %")
    with col4:
        st.metric("Mittlere Leistung", f"{avg_power:.0f} W")

def main():
    st.title("🧪 Bioreaktor Temperaturregelung")
    st.markdown("**PID-geregelte Temperaturkontrolle**")
    
    # Sidebar Parameter
    st.sidebar.header("⚙️ Parameter")
    
    # Bioreaktor
    st.sidebar.subheader("Bioreaktor")
    volumen = st.sidebar.slider("Volumen [L]", 10, 1000, 100)
    start_temp = st.sidebar.slider("Starttemperatur [°C]", 10.0, 50.0, 20.0)
    ambient_temp = st.sidebar.slider("Umgebungstemperatur [°C]", 5.0, 40.0, 20.0)
    
    # PID
    st.sidebar.subheader("PID-Regler")
    kp = st.sidebar.slider("Kp (Proportional)", 0.0, 10.0, 2.0, 0.1)
    ki = st.sidebar.slider("Ki (Integral)", 0.0, 2.0, 0.1, 0.01)
    kd = st.sidebar.slider("Kd (Differential)", 0.0, 2.0, 0.5, 0.1)
    
    # Simulation
    st.sidebar.subheader("Simulation")
    target_temp = st.sidebar.slider("Zieltemperatur [°C]", 20.0, 80.0, 37.0)
    simulation_time = st.sidebar.slider("Simulationszeit [min]", 1, 60, 10)
    dt = st.sidebar.selectbox("Zeitschritt [s]", [0.1, 0.5, 1.0, 2.0], index=2)
    
    # Störungen
    st.sidebar.subheader("Störungen")
    disturbance_active = st.sidebar.checkbox("Störung aktivieren")
    disturbance_time = 0
    disturbance_power = 0
    disturbance_duration = 0
    
    if disturbance_active:
        disturbance_time = st.sidebar.slider("Störung bei [min]", 1, simulation_time-1, simulation_time//2)
        disturbance_power = st.sidebar.slider("Störungsleistung [W]", -2000, 2000, -500)
        disturbance_duration = st.sidebar.slider("Störungsdauer [min]", 0.5, 5.0, 1.0)
    
    # Simulation starten
    if st.button("🚀 Simulation starten", type="primary"):
        # Progress indicator
        progress_bar = st.progress(0)
        status_text = st.empty()
        status_text.text("Simulation wird gestartet...")
        
        try:
            # Simulation ausführen
            progress_bar.progress(50)
            status_text.text("Simulation läuft...")
            
            results = run_simulation(
                volumen, start_temp, ambient_temp, target_temp, 
                simulation_time, dt, kp, ki, kd,
                disturbance_active, disturbance_time, disturbance_power, disturbance_duration
            )
            
            progress_bar.progress(100)
            status_text.text("✅ Simulation abgeschlossen!")
            
            # Ergebnisse anzeigen
            plot_results(results, target_temp)
            
            st.subheader("📈 Analyse der Ergebnisse")
            analyze_results(results, target_temp)
            
        except Exception as e:
            st.error(f"Fehler bei der Simulation: {str(e)}")
        finally:
            progress_bar.empty()
            status_text.empty()

if __name__ == "__main__":
    main()