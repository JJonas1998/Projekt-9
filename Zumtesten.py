import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from CoolProp.CoolProp import PropsSI
from fluids import Prandtl

t_kelvin = 293,15

k = PropsSI("CONDUCTIVITY", "T", t_kelvin, "P", 101325, "Water")  # Wärmeleitfähigkeit in W/(m·K)
mu = PropsSI("VISCOSITY", "T", t_kelvin, "P", 101325, "Water")    # Dynamische Viskosität in Pa·s

print(k)
print(mu)