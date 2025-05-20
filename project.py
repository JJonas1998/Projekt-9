"""
Projekt: Modellierung von Fermentationsdaten mit Tkinter-Interface
Autoren: Milena Rühmann (1389270) und Christian Stehno (1386026)
Datum: 14.05.2025

Test

Dieses Programm simuliert, analysiert und visualisiert Fermentationsdaten mit einer Tkinter-basierten GUI.
Es betrachtet die Auswirkungen von Temperatur, pH-Wert und Substratkonzentration
auf die Produktbildung während eines Fermentationsprozesses.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns
from scipy import stats
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from datetime import datetime

class FermentationSimulator:
    """
    Klasse zur Simulation von Fermentationsprozessen mit verschiedenen Parametern.
    """
    
    def __init__(self, duration=48, time_step=0.01):
        """
        Initialisiert den Fermentationssimulator.
        
        Parameter:
        - duration: Gesamtdauer der Fermentation in Stunden
        - time_step: Zeitintervall für die Simulation in Stunden
        """
        self.duration = duration
        self.time_step = time_step
        self.time_points = np.arange(0, duration + time_step, time_step)
        self.data = None
        
        # Standardparameter für die Simulation
        self.organism_params = {
            'CHO': {
                'max_growth_rate': 0.25,
                'substrate_affinity': 8.0,
                'yield_coefficient': 0.4,
                'maintenance_coefficient': 0.03,
                'death_rate': 0.005,
                'optimal_temp': 37.0,
                'optimal_ph': 7.2,
                'temp_sensitivity': 0.08,
                'ph_sensitivity': 1.5
            },
            'E. coli': {
                'max_growth_rate': 0.6,
                'substrate_affinity': 4.0,
                'yield_coefficient': 0.5,
                'maintenance_coefficient': 0.04,
                'death_rate': 0.01,
                'optimal_temp': 37.0,
                'optimal_ph': 7.0,
                'temp_sensitivity': 0.1,
                'ph_sensitivity': 2.0
            },
            'S. cerevisiae': {
                'max_growth_rate': 0.35,
                'substrate_affinity': 6.0,
                'yield_coefficient': 0.45,
                'maintenance_coefficient': 0.05,
                'death_rate': 0.008,
                'optimal_temp': 30.0,
                'optimal_ph': 5.5,
                'temp_sensitivity': 0.12,
                'ph_sensitivity': 2.5
            }
        }

        self.default_organism = 'E. coli'
        # default_params = Kinetik aus E. coli + Standard-Anfangswerte
        self.default_params = self.organism_params[self.default_organism].copy()
        self.default_params.update({
            'initial_substrate': 100.0,
            'initial_biomass':   0.5,
            'initial_product':   0.0
        })
    
    def simulate(self, temperature_profile, ph_profile, params=None):
        """
        Führt eine Fermentationssimulation mit gegebenen Temperatur- und pH-Profilen durch.
        
        Parameter:
        - temperature_profile: Funktion, die für jeden Zeitpunkt die Temperatur zurückgibt
        - ph_profile: Funktion, die für jeden Zeitpunkt den pH-Wert zurückgibt
        - params: Optional, Dictionary mit angepassten Modellparametern
        
        Gibt DataFrame mit simulierten Fermentationsdaten zurück.
        """
        if params is None:
            params = self.default_params
        else:
            # Kombiniere gegebene Parameter mit Standardparametern
            complete_params = self.default_params.copy()
            complete_params.update(params)
            params = complete_params
        
        # Arrays für die Simulation initialisieren
        n_steps = len(self.time_points)
        substrate = np.zeros(n_steps)
        biomass = np.zeros(n_steps)
        product = np.zeros(n_steps)
        temperature = np.zeros(n_steps)
        ph = np.zeros(n_steps)
        
        # Anfangsbedingungen setzen
        substrate[0] = params['initial_substrate']
        biomass[0] = params['initial_biomass']
        product[0] = params['initial_product']
        
        # Temperatur und pH für jeden Zeitpunkt berechnen
        for i, t in enumerate(self.time_points):
            temperature[i] = temperature_profile(t)
            ph[i] = ph_profile(t)
        
        # Monod-Kinetik mit Temperatur- und pH-Einfluss
        for i in range(1, n_steps):
            dt = self.time_step
            t = self.time_points[i]
            
            # Temperatur- und pH-Effekte auf Wachstumsrate (Gauß-Verteilung)
            temp_effect = np.exp(-params['temp_sensitivity'] * (temperature[i-1] - params['optimal_temp'])**2)
            ph_effect = np.exp(-params['ph_sensitivity'] * (ph[i-1] - params['optimal_ph'])**2)
            
            # Monod-Kinetik für das Wachstum
            growth_rate = params['max_growth_rate'] * substrate[i-1] / (params['substrate_affinity'] + substrate[i-1])
            growth_rate *= temp_effect * ph_effect  # Modifiziert durch Umgebungsbedingungen
            
            # Biomassewachstum mit Absterben
            biomass_growth = growth_rate * biomass[i-1] * dt
            biomass_death = params['death_rate'] * biomass[i-1] * dt
            biomass[i] = biomass[i-1] + biomass_growth - biomass_death
            
            # Substratverbrauch für Wachstum und Erhaltung
            substrate_consumption_growth = biomass_growth / params['yield_coefficient']
            substrate_consumption_maintenance = params['maintenance_coefficient'] * biomass[i-1] * dt
            substrate[i] = max(0, substrate[i-1] - substrate_consumption_growth - substrate_consumption_maintenance)
            
            # Produktbildung (gekoppelt an Biomassewachstum und Substratverbrauch)
            product_formation = params['yield_coefficient'] * substrate_consumption_growth
            # Zusätzlicher Produktabbau bei hohen Temperaturen
            product_degradation = 0.005 * (max(0, temperature[i-1] - 35)) * product[i-1] * dt
            product[i] = product[i-1] + product_formation - product_degradation
        
        # Ergebnisse in DataFrame speichern
        self.data = pd.DataFrame({
            'Zeit (h)': self.time_points,
            'Substrat (g/L)': substrate,
            'Biomasse (g/L)': biomass,
            'Produkt (g/L)': product,
            'Temperatur (°C)': temperature,
            'pH': ph
        })
        
        return self.data
    
    def save_data(self, filename=None):
        """
        Speichert die simulierten Daten in eine CSV-Datei.
        
        Parameter:
        - filename: Optional, Name der Datei. Wenn nicht angegeben, wird ein Zeitstempel verwendet.
        
        Gibt den Dateipfad zurück.
        """
        if self.data is None:
            raise ValueError("Keine Daten vorhanden. Führen Sie zuerst eine Simulation durch.")
        
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"fermentation_data_{timestamp}.csv"
        
        self.data.to_csv(filename, index=False)
        print(f"Daten wurden in {filename} gespeichert.")
        return filename


class FermentationAnalyzer:
    """
    Klasse zur Analyse von Fermentationsdaten.
    """
    
    def __init__(self, data=None):
        """
        Initialisiert den Fermentationsanalysator.
        
        Parameter:
        - data: Optional, DataFrame mit Fermentationsdaten
        """
        self.data = data
    
    def load_data(self, filepath):
        """
        Lädt Fermentationsdaten aus einer CSV-Datei.
        
        Parameter:
        - filepath: Pfad zur CSV-Datei
        """
        self.data = pd.read_csv(filepath)
        return self.data
    
    def calculate_correlations(self):
        """
        Berechnet die Korrelationen zwischen allen Variablen in den Daten.
        
        Gibt DataFrame mit Korrelationskoeffizienten zurück.
        """
        if self.data is None:
            raise ValueError("Keine Daten vorhanden. Laden Sie zuerst Daten.")
            
        # Ausschließen der Zeitspalte für die Korrelation
        data_for_corr = self.data.drop(columns=['Zeit (h)']) if 'Zeit (h)' in self.data.columns else self.data
        
        # Pearson-Korrelationskoeffizienten berechnen
        correlations = data_for_corr.corr()
        return correlations


class TkFermentationDashboard:
    """
    Klasse zur Erstellung eines interaktiven Tkinter-Dashboards für Fermentationsdaten.
    """
    
    def __init__(self, master):
        """
        Initialisiert das Tkinter-Dashboard.
        
        Parameter:
        - master: Das Tkinter-Root-Fenster
        """
        self.master = master
        self.master.title("Fermentations-Simulation")
        self.master.geometry("1200x800")
        self.master.minsize(1000, 700)
        
        # Simulator und Analyzer initialisieren
        self.simulator = FermentationSimulator(duration=48, time_step=0.5)
        self.analyzer = FermentationAnalyzer()
        
        # Aktuelle Abbildungen
        self.current_figures = {
            'time_series': None,
            'heatmap': None,
            'product_analysis': None
        }
        
        # GUI erstellen
        self.create_widgets()
        
        # Standardwerte initialisieren und erste Simulation ausführen
        self.init_default_values()
        self.run_simulation()
    
    def create_widgets(self):
        """
        Erstellt die GUI-Widgets für das Dashboard.
        """
        # Hauptframe erstellen
        main_frame = ttk.Frame(self.master)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Links: Parameter-Frame
        param_frame = ttk.LabelFrame(main_frame, text="Fermentationsparameter")
        param_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # Rechts: Notebook mit Tabs für Grafiken
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Tabs für die verschiedenen Grafiken
        self.time_series_tab = ttk.Frame(self.notebook)
        self.heatmap_tab = ttk.Frame(self.notebook)
        self.product_analysis_tab = ttk.Frame(self.notebook)
        
        self.notebook.add(self.time_series_tab, text="Zeitverläufe")
        self.notebook.add(self.heatmap_tab, text="Korrelationen")
        self.notebook.add(self.product_analysis_tab, text="Produktbildung")
        self.data_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.data_tab, text="Daten")
        
        ttk.Label(param_frame, text="Organismus:").grid(row=0, column=0, sticky="w", pady=5)
        self.organism_var = tk.StringVar()
        organism_combo = ttk.Combobox(param_frame, textvariable=self.organism_var,
                                    values=list(self.simulator.organism_params.keys()),
                                    state="readonly")
        organism_combo.grid(row=0, column=1, pady=5)
        organism_combo.set('E. coli')  # Default-Auswahl

        # ----- Parameter-Widgets erstellen -----
        # Substrat
        ttk.Label(param_frame, text="Anfangssubstratkonzentration (g/L):").grid(row=0, column=0, pady=5, sticky="w")
        self.substrate_var = tk.DoubleVar()
        ttk.Scale(param_frame, from_=50, to=200, orient=tk.HORIZONTAL, 
                 variable=self.substrate_var, length=200).grid(row=1, column=0, pady=5)
        ttk.Label(param_frame, textvariable=self.substrate_var).grid(row=1, column=1, pady=5)
        
        # Biomasse
        ttk.Label(param_frame, text="Anfangsbiomasse (g/L):").grid(row=2, column=0, pady=5, sticky="w")
        self.biomass_var = tk.DoubleVar()
        ttk.Scale(param_frame, from_=0.1, to=2.0, orient=tk.HORIZONTAL, 
                 variable=self.biomass_var, length=200).grid(row=3, column=0, pady=5)
        ttk.Label(param_frame, textvariable=self.biomass_var).grid(row=3, column=1, pady=5)
        
        # ----- Temperaturparameter -----
        temp_frame = ttk.LabelFrame(param_frame, text="Temperaturprofil")
        temp_frame.grid(row=4, column=0, columnspan=2, pady=10, sticky="we")
        
        # Temperaturprofil-Typ
        self.temp_profile_var = tk.StringVar()
        ttk.Radiobutton(temp_frame, text="Konstant", variable=self.temp_profile_var, 
                       value="constant").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(temp_frame, text="Rampe", variable=self.temp_profile_var, 
                       value="ramp").grid(row=1, column=0, sticky="w")
        ttk.Radiobutton(temp_frame, text="Oszillierend", variable=self.temp_profile_var, 
                       value="oscillating").grid(row=2, column=0, sticky="w")
        
        # Anfangstemperatur
        ttk.Label(temp_frame, text="Anfangstemperatur (°C):").grid(row=3, column=0, pady=5, sticky="w")
        self.initial_temp_var = tk.DoubleVar()
        ttk.Scale(temp_frame, from_=20, to=40, orient=tk.HORIZONTAL, 
                 variable=self.initial_temp_var, length=180).grid(row=4, column=0, pady=5)
        ttk.Label(temp_frame, textvariable=self.initial_temp_var).grid(row=4, column=1, pady=5)
        
        # Zieltemperatur (für Rampenprofil)
        ttk.Label(temp_frame, text="Zieltemperatur (°C):").grid(row=5, column=0, pady=5, sticky="w")
        self.target_temp_var = tk.DoubleVar()
        ttk.Scale(temp_frame, from_=20, to=40, orient=tk.HORIZONTAL, 
                 variable=self.target_temp_var, length=180).grid(row=6, column=0, pady=5)
        ttk.Label(temp_frame, textvariable=self.target_temp_var).grid(row=6, column=1, pady=5)
        
        # ----- pH-Parameter -----
        ph_frame = ttk.LabelFrame(param_frame, text="pH-Profil")
        ph_frame.grid(row=7, column=0, columnspan=2, pady=10, sticky="we")
        
        # pH-Profil-Typ
        self.ph_profile_var = tk.StringVar()
        ttk.Radiobutton(ph_frame, text="Konstant", variable=self.ph_profile_var, 
                       value="constant").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(ph_frame, text="Absinkend", variable=self.ph_profile_var, 
                       value="decreasing").grid(row=1, column=0, sticky="w")
        ttk.Radiobutton(ph_frame, text="Oszillierend", variable=self.ph_profile_var, 
                       value="oscillating").grid(row=2, column=0, sticky="w")
        
        # Anfangs-pH
        ttk.Label(ph_frame, text="Anfangs-pH:").grid(row=3, column=0, pady=5, sticky="w")
        self.initial_ph_var = tk.DoubleVar()
        ttk.Scale(ph_frame, from_=4, to=9, orient=tk.HORIZONTAL, 
                 variable=self.initial_ph_var, length=180).grid(row=4, column=0, pady=5)
        ttk.Label(ph_frame, textvariable=self.initial_ph_var).grid(row=4, column=1, pady=5)
        
        # ----- Simulationssteuerung -----
        control_frame = ttk.Frame(param_frame)
        control_frame.grid(row=8, column=0, columnspan=2, pady=20, sticky="we")
        
        ttk.Button(control_frame, text="Simulation starten", command=self.run_simulation).pack(fill=tk.X, pady=5)
        ttk.Button(control_frame, text="Daten speichern", command=self.save_current_data).pack(fill=tk.X, pady=5)
        ttk.Button(control_frame, text="Daten laden", command=self.load_data).pack(fill=tk.X, pady=5)
    
    def init_default_values(self):
        """
        Initialisiert die Standardwerte für die GUI-Elemente.
        """
        self.substrate_var.set(100)
        self.biomass_var.set(0.5)
        self.temp_profile_var.set("constant")
        self.initial_temp_var.set(30)
        self.target_temp_var.set(35)
        self.ph_profile_var.set("constant")
        self.initial_ph_var.set(7.0)
    
    def define_temp_profile(self):
        """
        Definiert das Temperaturprofil basierend auf den GUI-Einstellungen.
        
        Gibt eine Funktion zurück, die für jeden Zeitpunkt t die Temperatur berechnet.
        """
        profile_type = self.temp_profile_var.get()
        initial_temp = self.initial_temp_var.get()
        target_temp = self.target_temp_var.get()
        
        if profile_type == "constant":
            return lambda t: initial_temp
        elif profile_type == "ramp":
            return lambda t: initial_temp + (target_temp - initial_temp) * min(1.0, t / 12)
        elif profile_type == "oscillating":
            return lambda t: initial_temp + 3 * np.sin(2 * np.pi * t / 8)
        else:
            return lambda t: initial_temp
    
    def define_ph_profile(self):
        """
        Definiert das pH-Profil basierend auf den GUI-Einstellungen.
        
        Gibt eine Funktion zurück, die für jeden Zeitpunkt t den pH-Wert berechnet.
        """
        profile_type = self.ph_profile_var.get()
        initial_ph = self.initial_ph_var.get()
        
        if profile_type == "constant":
            return lambda t: initial_ph
        elif profile_type == "decreasing":
            return lambda t: max(initial_ph - 0.1 * t, 4.0)
        elif profile_type == "oscillating":
            return lambda t: initial_ph + 0.01 * np.sin(2 * np.pi * t / 12)
        else:
            return lambda t: initial_ph
    
    def run_simulation(self):
        """
        Führt die Fermentationssimulation mit den aktuellen Parametern durch
        und aktualisiert alle Grafiken.
        """
        chosen = self.organism_var.get()
        base_params = self.simulator.organism_params.get(chosen, {})
        # Überschreibe nur die vom User per Slider gewählten Startwerte
        params = base_params.copy()
        params.update({
            'initial_substrate': self.substrate_var.get(),
            'initial_biomass':  self.biomass_var.get()
        })
        
        # Profile definieren
        temp_profile = self.define_temp_profile()
        ph_profile = self.define_ph_profile()
        
        # Simulation ausführen
        try:
            data = self.simulator.simulate(temp_profile, ph_profile, params)
            self.analyzer.data = data
            
            # Grafiken aktualisieren
            self.update_time_series_plot()
            self.update_heatmap_plot()
            self.update_product_analysis_plot()
            self.update_data_table()
            
            messagebox.showinfo("Simulation", "Simulation erfolgreich durchgeführt!")
        except Exception as e:
            messagebox.showerror("Fehler", f"Fehler bei der Simulation: {str(e)}")
    
    def update_time_series_plot(self):
        """
        Aktualisiert die Zeitverlauf-Grafik im entsprechenden Tab.
        """
        # Alte Grafik entfernen, falls vorhanden
        for widget in self.time_series_tab.winfo_children():
            widget.destroy()
        
        # Neue Grafik erstellen
        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        
        # Plot 1: Substrat, Biomasse und Produkt
        axes[0].plot(self.simulator.data['Zeit (h)'], self.simulator.data['Substrat (g/L)'], 'b-', label='Substrat')
        axes[0].plot(self.simulator.data['Zeit (h)'], self.simulator.data['Biomasse (g/L)'], 'g-', label='Biomasse')
        axes[0].plot(self.simulator.data['Zeit (h)'], self.simulator.data['Produkt (g/L)'], 'r-', label='Produkt')
        axes[0].set_ylabel('Konzentration (g/L)')
        axes[0].legend()
        axes[0].grid(True)
        axes[0].set_title('Zeitverlauf von Substrat, Biomasse und Produkt')
        
        # Plot 2: Temperatur
        axes[1].plot(self.simulator.data['Zeit (h)'], self.simulator.data['Temperatur (°C)'], 'r-')
        axes[1].set_ylabel('Temperatur (°C)')
        axes[1].grid(True)
        axes[1].set_title('Temperaturprofil')
        
        # Plot 3: pH-Wert
        axes[2].plot(self.simulator.data['Zeit (h)'], self.simulator.data['pH'], 'g-')
        axes[2].set_xlabel('Zeit (h)')
        axes[2].set_ylabel('pH-Wert')
        axes[2].grid(True)
        axes[2].set_title('pH-Profil')
        
        plt.tight_layout()
        
        # Grafik im Tab anzeigen
        canvas = FigureCanvasTkAgg(fig, master=self.time_series_tab)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        canvas.draw()
        
        # Referenz speichern
        self.current_figures['time_series'] = fig
    
    def update_heatmap_plot(self):
        """
        Aktualisiert die Korrelationsheatmap im entsprechenden Tab.
        """
        # Alte Grafik entfernen, falls vorhanden
        for widget in self.heatmap_tab.winfo_children():
            widget.destroy()
        
        # Korrelationen berechnen und Heatmap erstellen
        correlations = self.analyzer.calculate_correlations()
        
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(correlations, annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0, fmt='.2f', ax=ax)
        ax.set_title('Korrelationsmatrix der Fermentationsvariablen')
        plt.tight_layout()
        
        # Grafik im Tab anzeigen
        canvas = FigureCanvasTkAgg(fig, master=self.heatmap_tab)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        canvas.draw()
        
        # Referenz speichern
        self.current_figures['heatmap'] = fig
    
    def update_product_analysis_plot(self):
        """
        Aktualisiert die Produktbildungsanalyse im entsprechenden Tab.
        """
        # Alte Grafik entfernen, falls vorhanden
        for widget in self.product_analysis_tab.winfo_children():
            widget.destroy()
        
        # Figure mit mehreren Subplots erstellen
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.flatten()
        
        # 1. Produktbildung vs. Substratkonzentration
        sns.regplot(x='Substrat (g/L)', y='Produkt (g/L)', data=self.simulator.data, ax=axes[0])
        axes[0].set_title('Produkt vs. Substrat')
        
        # 2. Produktbildung vs. Temperatur
        sns.regplot(x='Temperatur (°C)', y='Produkt (g/L)', data=self.simulator.data, ax=axes[1])
        axes[1].set_title('Produkt vs. Temperatur')
        
        # 3. Produktbildung vs. pH
        sns.regplot(x='pH', y='Produkt (g/L)', data=self.simulator.data, ax=axes[2])
        axes[2].set_title('Produkt vs. pH')
        
        # 4. Produktbildung vs. Biomasse
        sns.regplot(x='Biomasse (g/L)', y='Produkt (g/L)', data=self.simulator.data, ax=axes[3])
        axes[3].set_title('Produkt vs. Biomasse')
        
        plt.tight_layout()
        
        # Grafik im Tab anzeigen
        canvas = FigureCanvasTkAgg(fig, master=self.product_analysis_tab)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        canvas.draw()
        
        # Referenz speichern
        self.current_figures['product_analysis'] = fig

    def update_data_table(self):
            # altes Widget entfernen
            for w in self.data_tab.winfo_children():
                w.destroy()

            df = self.simulator.data
            cols = list(df.columns)

            # Treeview und Scrollbars anlegen
            tree = ttk.Treeview(self.data_tab, columns=cols, show="headings")
            vsb  = ttk.Scrollbar(self.data_tab, orient="vertical",   command=tree.yview)
            hsb  = ttk.Scrollbar(self.data_tab, orient="horizontal", command=tree.xview)
            tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

            # Layout: Treeview + Scrollbars
            tree.grid(row=0, column=0, sticky="nsew")
            vsb.grid (row=0, column=1, sticky="ns")
            hsb.grid (row=1, column=0, sticky="ew")
            self.data_tab.rowconfigure(0, weight=1)
            self.data_tab.columnconfigure(0, weight=1)

            # Spaltenüberschriften setzen
            for c in cols:
                tree.heading(c, text=c)
                tree.column(c, width=100, anchor="center")

            # Zeilen einfügen
            for _, row in df.iterrows():
                tree.insert("", "end", values=list(row))
    
    def save_current_data(self):
        """
        Speichert die aktuellen Simulationsdaten in eine CSV-Datei.
        """
        if self.simulator.data is None:
            messagebox.showwarning("Warnung", "Keine Daten zum Speichern vorhanden.")
            return
            
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV Dateien", "*.csv"), ("Alle Dateien", "*.*")],
            title="Daten speichern"
        )
        
        if file_path:
            try:
                self.simulator.save_data(file_path)
                messagebox.showinfo("Erfolg", f"Daten wurden erfolgreich gespeichert in:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Fehler", f"Fehler beim Speichern der Daten: {str(e)}")
    
    def load_data(self):
        """
        Lädt Fermentationsdaten aus einer CSV-Datei.
        """
        file_path = filedialog.askopenfilename(
            filetypes=[("CSV Dateien", "*.csv"), ("Alle Dateien", "*.*")],
            title="Daten laden"
        )
        
        if file_path:
            try:
                self.analyzer.load_data(file_path)
                self.simulator.data = self.analyzer.data
                
                # Grafiken aktualisieren
                self.update_time_series_plot()
                self.update_heatmap_plot()
                self.update_product_analysis_plot()
                self.update_data_table()
                
                messagebox.showinfo("Erfolg", f"Daten wurden erfolgreich geladen aus:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Fehler", f"Fehler beim Laden der Daten: {str(e)}")


def main():
    """
    Hauptfunktion zum Starten der Tkinter-Anwendung.
    """
    root = tk.Tk()
    root.state('zoomed')
    app = TkFermentationDashboard(root)
    root.mainloop()


if __name__ == "__main__":
    main()