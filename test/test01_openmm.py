import os
import re

import openmm as mm
from openmm.app import PDBFile, ForceField, PME, HBonds, Residue, Topology, Atom, Simulation
from openmm.unit import nano, meters, pico, seconds, kelvin, atmospheres
from openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator, System

nanometers = nano * meters
picoseconds = pico * seconds

def check_all_atoms(pdbfile_name: str = '1vii-processed.pdb'):
    base_dir = os.getenv('BASE_DIR')
    pdb = PDBFile(f'{base_dir}/test/input/{pdbfile_name}')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*nanometers
    ewaldErrorTolerance = 0.0005
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 0.000001

    # Integration Options
    dt = 0.002*picoseconds
    temperature = 300*kelvin
    friction = 1.0/picoseconds
    pressure = 1.0*atmospheres
    barostatInterval = 25

    # Simulation Options
    steps = 500
    equilibrationSteps = 500
    platform = Platform.getPlatformByName('CUDA')
    platformProperties = {'Precision': 'single'}

    # Prepare the Simulation
    print('Building system...')
    topology: Topology = pdb.topology
    positions = pdb.positions
    system: System = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
        constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    simulation = Simulation(topology, system, integrator, platform, platformProperties)
    simulation.context.setPositions(positions)
        
    num_atoms = topology.getNumAtoms()
    num_bonds = topology.getNumBonds()
    
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
    
    num_paris_protein = num_atoms_protein*(num_atoms_protein-1)//2
    print('num atoms (includes waterbox, solvent):', num_atoms)
    print('num bonds (includes waterbox solvent):', num_bonds)
    print('num atoms:', num_atoms_protein)
    print('num pairs:', num_paris_protein)
    atoms_init = get_atoms_nonbonded(simulation.system, atom_ids)
    
    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    atoms_energy_minimized = get_atoms_nonbonded(simulation.system, atom_ids)

    for id in range(len(atoms_init)):
        if atoms_init[id] != atoms_energy_minimized[id]:
            print(atoms_init[id], "->", atoms_energy_minimized[id])
    
    print('Equilibrating...')
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(equilibrationSteps)
    atoms_preprocessed = get_atoms_nonbonded(simulation.system, atom_ids)
        
    for id in range(len(atoms_init)):
        if atoms_init[id] != atoms_preprocessed[id]:
            print(atoms_init[id], "->", atoms_preprocessed[id])
    
def get_atoms_nonbonded(system: System, atom_ids: list[int]):
    non_bonded_forces: mm.NonbondedForce
    for force in system.getForces():
        if type(force) == mm.NonbondedForce:
            non_bonded_forces = force
            break
        assert('No atoms recorded in NonbondedForce instance.')
    
    atoms = []
    for id in atom_ids:
        atoms.append(non_bonded_forces.getParticleParameters(id))
    
    return atoms
    
if __name__ == "__main__":
    check_all_atoms('1vii-processed.pdb')