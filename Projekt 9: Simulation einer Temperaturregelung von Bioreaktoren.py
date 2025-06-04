#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

# Streamlit Konfiguration  

st.set_page_config(
    page_title = "Simulation einer Temperaturregelung von Bioreaktoren",
    page_icon = "🧪",
    layout = "wide"
)

## Bioreaktor-Klasse
# Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen = 100, start_temp = 20, rpm = 100):
        self.volumen = volumen / 1000                                                                           # Umrechnung Volumen des Reaktorinhalts von l in m³
        self.spez_c = PropsSI("Cpmass", "T", (start_temp + 273.15), "P", 101325, "Water")                       # Spezifische Wärmekapazität des Mediums(Wasser) in J/(kg*K) 
        self.dichte = PropsSI("Dmass" ,"T", (start_temp + 273.15), "P", 101325, "Water")                        # Dichte des Mediums (Wasser) in kg/m³
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)                                                     # Berechnung Radius des Bioreaktors (Zylinder) in m - (Annahme: Höhe h = 2 * r)
        self.hoehe  = 2 * self.radius                                                                           # Höhe des Bioreaktors in m 
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe                    # Innenfläche des Bioreaktor in m² 
        self.umg_temp = start_temp                                                                              # Umgebungstemperatur in °C
        self.start_temp = start_temp                                                                            # Reaktor-Innentemperatur in °C (Starttemperatur)
        self.rpm = rpm                                                                                          # Drehzahl des Rührers in 1/min
        self.d_ruehrer = (2 * self.radius) / 3                                                                  # Durchmesser des Rührers in m (Annahme: Rührerdurchmesser 1/3 des Reaktordurchmessers)
        self.waerme_h = self.berechnung_h()                                                                     # Berechnung des Wärmeübergangskoeffizienten in W/(m²*K)     
        
    def berechnung_h(self):
        """
        Berechnung des Wärmeübergangskoeffizienten h anhand von Impeller-Drehzahl (rpm) und Rührerdurchmesser.
        """
        ist_temp = self.start_temp + 273.15                                                                     # Umrechnung der Starttemperatur in Kelvin                           
        k   = PropsSI("CONDUCTIVITY", "T", ist_temp, "P", 101325, "Water")                                      # Wärmeleitfähigkeit des Mediums (Wasser) in W/(m*K)    
        mu  = PropsSI("VISCOSITY" ,"T", ist_temp, "P", 101325, "Water")                                         # Dynamische Viskosität des Mediums (Wasser) in Pa*s     
        Re = ((self.rpm / 60) * (self.d_ruehrer ** 2) * self.dichte) / mu                                       # Reynolds-Zahl des Mediums (Wasser) im Bioreaktor    
        Pr = Prandtl(self.spez_c, mu, k)                                                                        # Prandtl-Zahl des Mediums (Wasser) im Bioreaktor    

        if 4.5e3 < Re < 1e4 and 0.6 < Pr < 160:
            Nu = 0.354 * (Re ** 0.714) * (Pr ** 0.260)                                                          # Nusselt-Zahl für turbulente Strömung (Impeller)                                       
        elif Re >= 1e4 and Pr >= 0.6:
            Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)                                                              # Nusselt-Zahl für turbulente Strömung (Dittus-Boelter)                                                  
        else:
            Nu = 3.66                                                                                           # Nusselt-Zahl für laminare Strömung       
                                                                                                 
        h = Nu * k / self.d_ruehrer                                                                             # Wärmeübergangskoeffizient in W/(m²*K)   
        return h       

    def update_temperature(self, heizleistung, stoerung = 0, dt = 1):

        masse = self.volumen * self.dichte                                                                      # Masse des Reaktorinhalts in kg
        waermekapazitaet = masse * self.spez_c                                                                  # Wärmekapazität des Reaktorinhalts in J/K    

        waermeverlust = self.waerme_h * self.flaeche * (self.start_temp - self.umg_temp)                        # Wärmeverlust an die Umgebung in W 
        waermeaenderung = heizleistung - waermeverlust + stoerung                                               # Änderung der Wärme im Reaktor in W    
        dT = (waermeaenderung * dt) / waermekapazitaet                                                          # Temperaturänderung in K (dT = dQ / C)
        self.current_temp += dT                                                                                 # Aktualisierung der aktuellen Temperatur in °C
        return self.current_temp                                                
        

    

         

