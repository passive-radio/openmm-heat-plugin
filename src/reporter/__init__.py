from __future__ import absolute_import
__author__ = "Yudai Okubo"
__version__ = "0.0.0.dev1"

import re

from openmm import State, System, NonbondedForce
from openmm.app import Simulation, Topology, Atom, Residue
from openmm.unit import *

from utils import get_atoms_nonbonded

class InterAtomicQuantityReporter():
    def __init__(self, file, report_interval) -> None:
        self.__report_interval = report_interval
        
        self.__out = open(file, "+bw")
        self.__write_header()
    
    def __write_header(self):
        self.__out.write("step,donor,accepter,quantity\n".encode('UTF-8'))
        
    def describeNextReport(self, simulation: Simulation):
        step = self.__report_interval - simulation.currentStep%self.__report_interval
        
        return (step, self, True, True, True, True)
    
    def report(self, simulation: Simulation, state: State):
        """report method called by openmm.simulation instance.
        """
        
class InterAtomicForceReporter(InterAtomicQuantityReporter):
    def __init__(self, simulation: Simulation, state: State, file, report_interval, cutoff: float = None) -> None:
        super().__init__(file, report_interval)
        
        self.__cutoff = cutoff
        self.__setup()
        
    def __write_header(self):
        return super().__write_header()

    def describeNextReport(self, simulation: Simulation):
        return super().describeNextReport(simulation)
    
    def report(self, simulation: Simulation, state: State):
        return super().report(simulation, state)
    
    def __setup(self, simulation: Simulation):
        """OpenMM implemented Force C++ object
        CMAPTorsionForce
        DrudeForce
        GBSAOBCForce
        GayBerneForce
        HarmonicAngleForce
        HarmonicBondForce
        NonbondedForce
        PeriodicTorsionForce
        RBTorsionForce
        """
        
        state = simulation.context.getState(True, True, True, True, True, True, True, True, True)
        system: System = simulation.system
        topology: Topology = simulation.topology
        self.__num_atoms = topology.getNumAtoms()
        
        pattern = re.compile(r'HOH|NA|CL')
        residues = []
        for residue in topology.residues():
            residue: Residue
            if bool(pattern.match(residue.name)):
                continue
            else:
                residues.append(residue)
        
        num_atoms_protein: int = 0
        atom_ids = []
        for residue in residues:
            for atom in residue.atoms():
                atom: Atom
                atom_ids.append(atom.index)
                num_atoms_protein += 1
        
        atoms = get_atoms_nonbonded(system, atom_ids)