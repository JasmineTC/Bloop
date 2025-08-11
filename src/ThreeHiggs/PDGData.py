"""
Make data from PDG 2025 importable.
Note: There is a PDG python API but it is missing needed data at time of
writing. So we'll use it that where possible but some things need to be added
by hand.

Important to note about the PDG API is the data set only gets updated when you
update the PDG module. 
"""
import math as m
import pdg
from scipy import constants
api = pdg.connect()
## Should this be a dict to make importing less tedious?

mHiggs = api.get_particle_by_name('H').mass
mTop = api.get_particle_by_name('t').mass
mW = api.get_particle_by_name('W+').mass
mZ = api.get_particle_by_name('Z0').mass

## This the higgs vev via the fermi constant
higgsVEV = 1/m.sqrt((m.sqrt(2)*constants.physical_constants["Fermi coupling constant"][0]))

from unittest import TestCase
class PDGUnitTests(TestCase):
    def test_HiggsMass(self):
        self.assertEqual(125.199530409718, mHiggs)
    
    def test_TopMass(self):
        self.assertEqual(172.559088345398, mTop)
        
    def test_ZMass(self):
        self.assertEqual(91.1879780919373, mZ)
    
    def test_WMass(self):
        self.assertEqual(80.377, mW)

    def test_Fermi(self):
        self.assertEqual(1.1663787e-05, constants.physical_constants["Fermi coupling constant"][0])
