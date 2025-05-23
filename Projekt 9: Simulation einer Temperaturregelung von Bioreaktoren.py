
#### Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren

### Ersteller/-in: Jonas Jahrstorfer, Johanna Niklas

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint


## Bioreaktor-Klasse
# Schreiben Sie eine Python-Funktion zur Simulation von Temperaturänderungen in einem Bioreaktor

class Bioreaktor:

    """
    Simuliert einen Bioreaktor mit Temperaturregelung. 
    Die Temperatur wird durch Heizen/Kühlen, Umgebungseinfluss und optionale Störungen verändert.
    """

    def __init__(self, start_temp = 25.0, umgebung = 22.0, waermeverlust = 0.01, max_heat = 5.0, max_cool = 5.0):
        self.temperatur = start_temp        # aktuelle Temperatur (°C)
        self.umgebung = umgebung            # Umgebungstemperatur (°C)
        self.waermeverlust = waermeverlust  # Wärmeverlustkoeffizient (1/s)
        self.max_heat = max_heat            # maximale Heizleistung (°C/s)
        self.max_cool = max_cool            # maximale Kühlleistung (°C/s)
        self.stoerung = 0.0                 # Störgröße (Wärmezufuhr/Abfuhr, wirkt einmalig)
    
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



