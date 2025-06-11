# Projekt 9: Simulation einer Temperautrregelung von Bioreaktoren.
**Ersteller/ -in**: Jonas Jahrstorfer, Johanna Niklas

**Abgabedatum**: 13.06.2025

## 1. Ziel des Programms
In diesem Projekt wird ein Bioreaktor-Modell mit einer simulierten Temperaturregelung umgesetzt. Ziel ist es, das thermische Verhalten eines Bioreaktors unter verschiedenen Betriebsbedingungen zu analysieren und durch eine PID-Regelung stabil auf eine gewünschte Solltemperatur zu bringen. Die Benutzeroberfläche wird mit Streamlit realisiert, wodurch eine interaktive Anpassung aller relevanten Parameter in Echtzeit möglich ist.

## 2. Funktionsweise der Anwendung
Die Anwendung kombiniert thermophysikalische Berechnungen mit einer Regelstrecke. Der Nutzer kann Parameter wie Reaktorgröße, Materialeigenschaften, Regelparameter und Umgebungstemperaturen anpassen. Der Temperaturverlauf des Reaktorinhalts wird daraufhin berechnet und grafisch dargestellt – mit und ohne aktive Regelung.

### PID - Regelung
Ein PID - Regler besteht aus drei Komponenten:

- **Proportional (P)**: Reagiert auf den aktuellen Unterschied zur Solltemperatur.

- **Integral (I)**: Reagiert auf die Summe aller vergangenen Abweichungen.

- **Differential (D)**: Reagiert auf die Änderungsgeschwindigkeit der Abweichung.

Diese Regelung erlaubt eine präzise Steuerung der Heizleistung, sodass Schwankungen schnell ausgeglichen werden und eine stabile Zieltemperatur erreicht wird.

## 3. Installation und Systemvoraussetzungen
Um den Code auszuführen sind zwei Vorraussetzungen notwendig:

- Python 3.9 oder höher
- Internetverbidung (für evtl. Nachinstallationen von Bibliotheken)

### Notwendige Bibliotheken:
Das Programm verwendet folgende externe Bibliotheken:
- ``` streamlit``` - UI-Komponenten für Web - App.
- ```matplotlib``` - Diagramme und Visualisieriung.
- ```CoolProp``` - Zugriff auf thermophysikalische Stoffdaten.
- ``` fluids``` - Prandtl-Zahl-Berechnung für Wärmeübertragung.

Die Installation erfolgt via: 
```bash 
pip install streamlit matplotlib CoolProp fluids
```

## 4. Start der Anwendung

1. Öffne das Terminal (z.B. in VS Code) oder die Eingabeaufforderung
2. Navigiere in den Ordner, in dem sich die Datei ... befindet
3. Starte das Programm mit dem Befehl:
```bash
streamlit run Projekt....py
```
4. Die Benutzeroberfläche öffnet sich im Standardbrowser.

## 5. Bedienungsanleitung
### Parametersteuerung
Alle Einstellungen werden über die linke Seitenleiste vorgenommen:

#### **Physikalische Parameter**
| Parameter                  | Funktion                                             |
| -------------------------- | ---------------------------------------------------- |
| **Reaktorvolumen (L)**            | Flüssigkeitsmenge im Reaktor.                        |
| **Starttemperatur (C°)**          | Starttemperatur der Flüssigkeit. |
| **Umgebungstemperatur**    | Temperatur der Luft um den Reaktor                   |
| **Wandmaterial** | Unterschiedliche Wäremleitfähigkeit (z.B. Edelstahl, Glas). Beeinflusst Wärmeleitfähigkeit und -speicherung.
   **Wandstärke (mm)**	|Bestimmt den Wärmewiderstand der Reaktorwand.  |
|**Rührerdrehzahl (1/min)** | Beeinflusst die Durchmischung und den Wärmeübergang im Reaktor. Höhere Drehzahlen fördern eine gleichmäßige Temperaturverteilung.

#### **Sollwert und Simulation**
| Parameter                  | Funktion                                             |
| -------------------------- | ---------------------------------------------------- |
| **Solltemperatur (°C)**    | Zieltemperatur, auf die geregelt werden soll.        |
| **Simulationsdauer**       | Dauer der Simulation in Minuten.                     |
|**Zeitschritte (s)**| Gibt an, in wie viele Abschnitte die Simulationszeit unterteilt wird. Höhere Werte führen zu genaueren, aber rechenintensiveren Ergebnissen.

#### **PID - Parameter**
 
| Parameter                  | Funktion                                             |
| -------------------------- | ---------------------------------------------------- |                                                                             
| **Kp (Proportional)** | Reagiert direkt auf die aktuelle Abweichung vom Sollwert. Je größer die Abweichung, desto stärker die Korrektur. |
| **Ki (Integral)**     | Summiert alle vergangenen Abweichungen. Verhindert bleibende Fehler (Offset), wirkt aber langsam.                |
| **Kd (Differential)** | Reagiert auf die Änderungsrate der Abweichung. Bremst schnelle Änderungen ab und reduziert Überschwingen.        |

#### **Erweiterte Einstellungen**
| Parameter                  | Funktion                                             |
| -------------------------- | ---------------------------------------------------- |
| **Störung aktivieren**     | Simulation einer plötzlichen Temperaturänderung      |
|**Detaillierte Analyse aktivieren**| Ermöglicht eine qualitative und quantitative Bewertung des Regelverhaltens.
|**Störungszeitpunkt (min)** | Beschreibt den Zeitpunkt, wann die Störung auftritt.
|**Störungsgröße (°C)**| Ist ein unerwarteter Einfluss von außen, z.B. plötzliche Temperaturschwankungen|

## 6. Tabs der Benutzeroberfläche

**Simulation**:

- Zeigt Temperaturverläufe mit und ohne Regelung. 
- Reagiert sofort auf Änderungen der Parameter.

**Analyse**:
 - Liefert Kennzahlen zur Regelgüte, etwa Einschwingzeit, Regelabweichung oder Energiebedarf.
 - Ermöglicht eine Bewertung der gewählten PID Einstellungen.

**Info**:
- Erläutert die theroretischen Grundlagen (z.B. PID-Regelung, Wärmeübertragung).
- Hinweise zur Interpretation der Simulation.

## 7. Interpretation der Ergebnisse

## 8. Fazit
Die Anwendung zeigt auf verständliche Weise, wie eine Temperaturregelung funktioniert, welche physikalischen Größen relevant sind und wie sich unterschiedliche Regler-Einstellungen auf das Verhalten des Systems auswirken. Sie eignet sich als Grundlage für weiterführende Projekte in Biotechnologie oder Verfahrenstechnik.


               

