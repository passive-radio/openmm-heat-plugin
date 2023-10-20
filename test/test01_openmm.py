import os
import re

import openmm as mm
from openmm.app import PDBFile, ForceField, PME, HBonds, Residue, Topology, Atom, Simulation
from openmm.app.topology import Bond
from openmm.unit import nano, meters, pico, seconds, kelvin, atmospheres
from openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator, System

nanometers = nano * meters
picoseconds = pico * seconds

def check_all_atoms(pdbfile_name: str = '108l-processed.pdb'):
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
    
    forces = system.getForces()
    print(forces)
        
    pattern = re.compile(r'HOH|NA|CL')
    residues = []
    for residue in topology.residues():
        residue: Residue
        if bool(pattern.match(residue.name)):
            continue
        else:
            residues.append(residue)
            
    atoms = []
    for atom in topology.atoms():
        atoms.append(atom)
            
    bonds = []
    for bond in topology.bonds():
        bonds.append(bond)
    
    bonds_ids = [i for i in range(topology.getNumBonds())]
    
    bonds_exclude_hoh = []
    for bond in bonds:
        bond: Bond
        atom1: Atom = bond.atom1
        atom2: Atom = bond.atom2
                
        if bool(pattern.match(atom1.residue.name)) or bool(pattern.match(atom2.residue.name)):
            continue
        else:
            bonds_exclude_hoh.append(bond)
    
    harmonic_bonds = get_harmonic_bond_forces(system)
    harmonic_angle_groups = get_harmonic_angle_forces(system)
    print(harmonic_bonds[0])
    print(harmonic_angle_groups[0])
    
    print('num bonds:', len(bonds_ids))
    print('num bonds exlude hoh:', len(bonds_exclude_hoh))
    print('num harmonic bonds:', len(harmonic_bonds))
    print('num harmonic angle groups:', len(harmonic_angle_groups))
    
    with open(f'{base_dir}/test/out/bonds_exclude_hoh.txt', mode='+bw') as f:
        bonds = []
        for i, bond in enumerate(bonds_exclude_hoh):
            atom1: Atom = bond.atom1
            atom2: Atom = bond.atom2
            bonds.append(f"{i}, {atom1.id}_{atom1.name}, {atom1.residue.id}_{atom1.residue.name}, {atom2.id}_{atom2.name}, {atom2.residue.id}_{atom2.residue.name}\n".encode('utf-8'))
        
        print(bonds[0])
        f.writelines(bonds)
    
    with open(f'{base_dir}/test/out/bonds_harmonic.txt', mode='+bw') as f:
        bonds = []
        for i, bond in enumerate(harmonic_bonds):
            atom1 = atoms[bond[0]]
            atom2 = atoms[bond[1]]
            bonds.append(f"{i}, {atom1.id}_{atom1.name}, {atom1.residue.id}_{atom1.residue.name}, {atom2.id}_{atom2.name}, {atom2.residue.id}_{atom2.residue.name}\n".encode('utf-8'))
        f.writelines(bonds)
        
    with open(f'{base_dir}/test/out/groups_harmonic_angle.txt', mode='+bw') as f:
        bonds = []
        for i, group in enumerate(harmonic_angle_groups):
            atom1 = atoms[group[0]]
            atom2 = atoms[group[1]]
            atom3 = atoms[group[2]]
            bonds.append(f"{i}, {atom1.id}_{atom1.name}, {atom1.residue.id}_{atom1.residue.name}, {atom2.id}_{atom2.name}, {atom2.residue.id}_{atom2.residue.name}, {atom3.id}_{atom3.name}, {atom3.residue.id}_{atom3.residue.name}\n".encode('utf-8'))
        f.writelines(bonds)
    
    # num_atoms_protein: int = 0
    # atom_ids = []
    # for residue in residues:
    #     for atom in residue.atoms():
    #         atom: Atom
    #         atom_ids.append(atom.index)
    #         num_atoms_protein += 1
    
    # num_paris_protein = num_atoms_protein*(num_atoms_protein-1)//2
    # print('num atoms (includes waterbox, solvent):', num_atoms)
    # print('num bonds (includes waterbox solvent):', num_bonds)
    # print('num atoms:', num_atoms_protein)
    # print('num pairs:', num_paris_protein)
    # atoms_init = get_atoms_nonbonded(simulation.system, atom_ids)
    
    # print('Performing energy minimization...')
    # simulation.minimizeEnergy()
    # atoms_energy_minimized = get_atoms_nonbonded(simulation.system, atom_ids)

    # for id in range(len(atoms_init)):
    #     if atoms_init[id] != atoms_energy_minimized[id]:
    #         print(atoms_init[id], "->", atoms_energy_minimized[id])
    
    # print('Equilibrating...')
    # simulation.context.setVelocitiesToTemperature(temperature)
    # simulation.step(equilibrationSteps)
    # atoms_preprocessed = get_atoms_nonbonded(simulation.system, atom_ids)
        
    # for id in range(len(atoms_init)):
    #     if atoms_init[id] != atoms_preprocessed[id]:
    #         print(atoms_init[id], "->", atoms_preprocessed[id])
    
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

def get_bonded_pairs(system: System, bond_ids: list[int]):
    bonded_forces = mm.HarmonicBondForce
    for force in system.getForces():
        if type(force) == mm.HarmonicBondForce:
            bonded_forces = force
            break
        assert('No atoms recorded in HarmonicBondForce instance.')
    
    bonds = []    
    for id in bond_ids:
        bonds.append(bonded_forces.getBondParameters(id))
    
    return bonds

def get_harmonic_bond_forces(system: System):
    bonded_forces: mm.HarmonicBondForce
    for force in system.getForces():
        if type(force) == mm.HarmonicBondForce:
            bonded_forces = force
    
    bonds = []
    num_bonds = bonded_forces.getNumBonds()
    for id in range(num_bonds):
        bonds.append(bonded_forces.getBondParameters(id))
    return bonds

def get_harmonic_angle_forces(system: System):
    bonded_forces: mm.HarmonicAngleForce
    for force in system.getForces():
        if type(force) == mm.HarmonicAngleForce:
            bonded_forces = force
    
    bonds = []
    num_bonds = bonded_forces.getNumAngles()
    for id in range(num_bonds):
        bonds.append(bonded_forces.getAngleParameters(id))
    return bonds

if __name__ == "__main__":
    check_all_atoms('1vii-processed.pdb')