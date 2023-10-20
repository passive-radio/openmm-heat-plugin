from openmm import System, NonbondedForce

def get_atoms_nonbonded(system: System, atom_ids: list[int]):
    non_bonded_forces: NonbondedForce
    for force in system.getForces():
        if type(force) == NonbondedForce:
            non_bonded_forces = force
            break
        assert('No atoms recorded in NonbondedForce instance.')
    
    atoms = []
    for id in atom_ids:
        atoms.append(non_bonded_forces.getParticleParameters(id))
    
    return atoms