## Streamlit-Anwendung
st.title("Projekt 9: Temperaturregelung eines Bioreaktors")
st.subheader("Simulation mit und ohne PID-Regelung")

# Benutzer-Eingaben
col1, col2 = st.columns(2)

with col1:
    ziel_temp = st.number_input("Zieltemperatur [°C]", 20.0, 100.0, 37.0)
    start_temp = st.number_input("Starttemperatur [°C]", 20.0, 60.0, 25.0)
    umgeb_temp = st.number_input("Umgebungstemperatur [°C]", 0.0, 40.0, 22.0)
    laufzeit = st.slider("Simulationszeit [s]", 10, 1000, 300)
    timestep = st.slider("Zeitschritt [s]", 1, 60, 5)

with col2:
    st.markdown("**PID-Regler Parameter**")
    kp = st.slider("Proportionalanteil Kp", 0.0, 50.0, 10.0)
    ki = st.slider("Integralanteil Ki", 0.0, 10.0, 1.0)
    kd = st.slider("Differentialanteil Kd", 0.0, 10.0, 0.5)
    max_leistung = st.slider("Max. Heizleistung [W]", 0, 3000, 1000)

# Reaktor- und Reglerinitialisierung
reaktor_geregelt = Bioreaktor(start_temp=start_temp, umg_temp=umgeb_temp)
reaktor_ungeregelt = Bioreaktor(start_temp=start_temp, umg_temp=umgeb_temp)
regler = PIDRegler(dt=timestep, kp=kp, ki=ki, kd=kd, output_max=max_leistung)

# Simulation
zeiten = np.arange(0, laufzeit, timestep)
temp_geregelt = []
temp_ungeregelt = []

for t in zeiten:
    # PID-Regelung
    leistung = regler.run(ziel_temp, reaktor_geregelt.ist_temp)
    temp_geregelt.append(reaktor_geregelt.update_temperature(leistung, dt=timestep))
    
    # Unkontrolliertes System (ohne Heizleistung)
    temp_ungeregelt.append(reaktor_ungeregelt.update_temperature(0, dt=timestep))

# Plotten
fig, ax = plt.subplots()
ax.plot(zeiten, temp_geregelt, label="Mit PID-Regelung", linewidth=2)
ax.plot(zeiten, temp_ungeregelt, label="Ohne Regelung", linestyle='--')
ax.axhline(y=ziel_temp, color='gray', linestyle=':', label="Zieltemperatur")
ax.set_xlabel("Zeit [s]")
ax.set_ylabel("Temperatur [°C]")
ax.set_title("Temperaturverlauf im Bioreaktor")
ax.legend()
ax.grid(True)

st.pyplot(fig)

# Debug-Werte anzeigen
with st.expander("Simulationsdetails"):
    st.write("Letzte geregelte Temperatur: ", round(temp_geregelt[-1], 2), "°C")
    st.write("Letzte unkontrollierte Temperatur: ", round(temp_ungeregelt[-1], 2), "°C")
    st.write("Letzte Heizleistung: ", round(leistung, 2), "W")
