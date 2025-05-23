#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from CoolProp.CoolProp import PropsSI
from fluids import Reynolds, Prandtl

## Bioreaktor-Klasse
# Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:
    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """
    def __init__(self, volumen = 100, spez_c = 4180, dichte = 1000, start_temp = 20, rpm = 100, d_ruehrer = 0,1):
        self.volumen = volumen                                                                                  # Volumen des Reaktorinhalts in Liter
        self.spez_c = PropsSI("Cpmass", "T", (start_temp + 273.15), "P", 101325, "Water")                       # spezifische Wärmekapazität des Mediums(Wasser) in J/(kg*K) - Wasser
        self.dichte = dichte                                                                                    # Dichte des Mediums in kg/m³
        self.flaeche = 2 * np.pi * (volumen / 1000) ** (2 / 3)                                                  # Fläche des Bioreaktor in m² (Zylinder-Äquivalent => geschätzt)
        self.umg_temp = start_temp                                                                              # Umgebungstemperatur in °C
        self.start_temp =  start_temp                                                                           # Reaktor-Innentemperatur in °C (Starttemperatur)
        self.rpm = rpm                                                                                          # Drehzahl des Rührers
        self.d_ruehrer = d_ruehrer                                                                              # Durchmesser des Rührers in Meter
    """
    Berechnung Wärmeübergangskoeffizient über Dittus–Boelter aus rpm und Impellerdurchmesser.
    """
    def berechnung_h_(self):
        temperatur = self.start_temp + 273.15
        k   = PropsSI("CONDUCTIVITY", "T", temperatur, "P", 101325, "Water")
        mu  = PropsSI("VISCOSITY" ,"T", temperatur, "P", 101325, "Water")
        rho = PropsSI("Dmass" ,"T", temperatur, "P", 101325, "Water")
        v = (self.rpm / 60) * np.pi * self.d_ruehrer
        Re = Reynolds(rho, v, self.d_ruehrer, mu)
        Pr = Prandtl(self.spez_z, mu, k)
        Nu = (0.037 * Re ** 0.8 * Pr) / (1 + 2.443 * Re ** (- 0.1) * (Pr ** (2 / 3) - 1))
        h = Nu * k / self.d_ruehrer
        return h       











    def set_stoerung(self, wert):
        """Setzt die Störgröße (z.B. plötzliche Wärmezufuhr, wirkt einmalig im nächsten Schritt)."""
        self.stoerung = wert
    
    def update(self, leistung, dt=1.0):
        """
        Führt einen Simulationsschritt durch.
        leistung: Heiz- oder Kühlleistung (vom Regler, wird automatisch begrenzt)
        dt: Zeitschritt (Standard 1 Sekunde)
        """
        # Begrenzung der Heiz-/Kühlleistung
        leistung = min(max(leistung, -self.max_cool), self.max_heat)

        # Temperaturdifferenz zur Umgebung (Abkühlung/Erwärmung)
        verlust = -self.waermeverlust * (self.temperatur - self.umgebung) * dt

        # Temperaturänderung durch Leistung + Störgröße
        self.temperatur += verlust + leistung * dt + self.stoerung

        # Störgröße wirkt nur einmal
        self.stoerung = 0.0
        return self.temperatur
    
    def reset(self, start_temp=None):
        """Setzt die Temperatur auf Startwert zurück."""
        if start_temp is not None:
            self.temperatur = start_temp

## PID-Regler-Klasse
# Implementieren Sie einen PID-Regler, um eine Zieltemperatur zu halten.



