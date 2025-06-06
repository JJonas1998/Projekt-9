#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

## Bioreaktor-Klasse
#  Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen = 10, start_temp = 20, umg_temp = 20, rpm = 100, wandmaterial = 'steel', wandstaerke = 0.005):

        # Thermophysikalische Eigenschaften des Mediums (Wasser)
        temp_k = start_temp + 273.15                                                                # Temperatur in Kelvin
        self.volumen = volumen / 1000                                                               # Umrechnung Volumen des Reaktorinhalts von l in m³
        self.spez_c = PropsSI("Cpmass", "T", (temp_k), "P", 101325, "Water")                        # Spezifische Wärmekapazität des Mediums(Wasser) in J/(kg*K) 
        self.dichte = PropsSI("Dmass" ,"T", (temp_k), "P", 101325, "Water")                         # Dichte des Mediums (Wasser) in kg/m³

        # Geometrie des Bioreaktors (Zylinder, Annahme: H = 2 * r)
        self.radius = (self.volumen / (2 * np.pi)) ** (1/3)                                         # Radius des Bioreaktors in m 
        self.hoehe  = 2 * self.radius                                                               # Höhe des Bioreaktors in m 
        self.flaeche = 2 * np.pi * (self.radius ** 2) + 2 * np.pi * self.radius * self.hoehe        # Innenfläche des Bioreaktor in m² 
        self.wandstaerke = wandstaerke                                                              # Wandstärke des Bioreaktors in m
        self.d_ruehrer = (2 * self.radius) / 3                                                      # Durchmesser des Rührers in m (Annahme: Rührerdurchmesser 1/3 des Reaktordurchmessers) 

        # Betriebsbedingungen
        self.umg_temp = umg_temp                                                                    # Umgebungstemperatur in °C
        self.ist_temp = start_temp                                                                  # Reaktor-Innentemperatur in °C (Starttemperatur)
        self.rpm = rpm                                                                              # Drehzahl des Rührers in 1/min

        # Wärmeübertragung                                                              
        self.waerme_h = self.berechnung_h()                                                         # Berechnung des Wärmeübergangskoeffizienten in W/(m²*K)     
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
        material_db = {
            'Stahl': 21.0,                                                                          # Edelstahl V2A
            'Glas': 1.4,                                                                            # Quarzglas                                                                              
            'Aluminium': 230.0                                                                      # Aluminium
        }
        return material_db.get(material.lower(), 16.0)                                              # default: Stahl     

    def berechnung_waermeverlust(self):
        """
        Berechnung des gesamten Wärmeverlusts: 
        - durch die Wand (Fourier)
        - durch Konvektion mit h (z. B. Luft/Wasserbad)
        """
        # Temperaturdifferenz zwischen Innen- und Außentemperatur
        delta_t = self.ist_temp - self.umg_temp                                                        

        # 1. Leitung durch Wandmaterial
        q_leitung = self.lambda_wand * (delta_t * self.flaeche / self.wandstaerke)                  # Fourier'sches Gesetz: Q = λ * ( ΔT * A / d) 

        # 2. Konvektiver Verlust an Umgebung
        q_konvektion = self.waerme_h * self.flaeche * delta_t                                       # Konvektionsverlust: Q = h * A * ΔT 

        # Gesamtverlust
        q_verlust = q_leitung + q_konvektion                                                        # Gesamtwärmeverlust in W (Wärmeleitung + Konvektion)    

        return q_verlust

    def update_temperature(self, leistung, zeitintervall = 1):
        """
        Aktualisiert die Innentemperatur des Reaktors auf Basis der Energiebilanz:
        ΔT = (Q_zufuhr - Q_verlust) * dt / (m * c)
        """
        # Berechnung der Masse des Reaktorinhalts und Wärmeverlust
        masse = self.dichte * self.volumen                                                          # Masse des Reaktorinhalts in kg 
        waermeverlust = self.berechnung_waermeverlust()                                             # Gesamtwärmeverlust in W 
        energie_netto = (leistung - waermeverlust) * zeitintervall                                  # Netto-Energiezufuhr in J (Leistung - Wärmeverlust) * Zeitintervall in s
        temp_delta = energie_netto / (masse * self.spez_c)                                          # Temperaturänderung in K 
        self.ist_temp += temp_delta                                                                 # Aktualisierung der Innentemperatur des Reaktors in °C

        return self.ist_temp

## PID-Regler-Klasse
#  Implementieren Sie einen PID-Regler, um eine Zieltemperatur zu halten

class PID:
    """
    PID-Regler zur Regelung der Temperatur eines Bioreaktors.
    """
    def __init__(self, dt = 0.0, kp = 0.0, ki = 0.0, kd = 0.0, output_min = 0, output_max = 1000):
        self.dt = dt                                                                                # Zeitschritt in Sekunden
        self.kp = kp                                                                                # Proportionalanteil  
        self.ki = ki                                                                                # Integralanteil    
        self.kd = kd                                                                                # Differentialanteil
        self.output_min = output_min                                                                # Minimalwert der Stellgröße (Heizleistung)
        self.output_max = output_max                                                                # Maximalwert der Stellgröße (Heizleistung)
        self.fehler_vor = 0.0                                                                       # Vorheriger Fehler für den D-Anteil
        self.dt_vor = -1e-6                                                                         # Vorheriger Zeitschritt (initialisiert mit einem sehr kleinen Wert)       
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

        # Offset-Berechnung nur wenn Bioreaktor-Referenz vorhanden
        if self.bioreaktor is not None:                                         
            offset = self.bioreaktor.berechnung_waermeverlust()          
        else:
            offset = 0.0

        # Fehlerberechnung
        fehler = sollwert - istwert                                                                  # Berechnung des Fehlers (Sollwert - Istwert)
        
        # P-Anteil
        p_anteil = self.kp * fehler                                                                  # Proportionalanteil
        
        # I-Anteil
        self.integral += fehler * self.dt                                                            # Integralwert aktualisieren
        i_anteil = self.ki * self.integral                                                           # Integralanteil  
        
        # D-Anteil
        d_anteil = self.kd * (fehler - self.fehler_vor) / (self.dt - self.dt_vor)                    # Differentialanteil (Ableitung des Fehlers)
        
        # Ausgabe berechnen
        stellgroesse = offset + p_anteil + i_anteil + d_anteil                                       # PID-Ausgabe (Stellgröße)
        
        # Begrenzung
        if stellgroesse > self.output_max:
            output = self.output_max
            self.integral -= fehler * self.dt                                                        # Anti-Windup
        elif stellgroesse < self.output_min:
            output = self.output_min
            self.integral -= fehler * self.dt                                                        # Anti-Windup

        self.dt_vor = self.dt                                                                        # Aktuellen Zeitschritt speichern    
        self.fehler_vor = fehler                                                                     # Aktuellen Fehler speichern
        return output
    
    def reset(self):
        """Regler zurücksetzen"""
        self.fehler_vor = 0.0
        self.integral = 0.0




