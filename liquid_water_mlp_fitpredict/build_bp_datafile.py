import pathlib
import re

import numpy as np
from ase.io import read as ase_read

def main():
    test_dir = pathlib.Path(__file__).resolve().parent.parent.joinpath('data')

    liquid_frames = ase_read(test_dir.joinpath("ice_in_water_data/dataset_1000_eVAng.xyz"), index=':')

    with open(test_dir.joinpath("ice_in_water_data/dataset_1000_eVAng.xyz"), 'r') as f:
        file_content = f.read()
        f.close()

    # The forces for the liquid
    liquid_forces = np.asarray(
        re.findall(
            "[A-z][\\s]*(-?[0-9]*\.[0-9]*)[\\s]*(-?[0-9]*\.[0-9]*)[\\s]*(-?[0-9]*\.[0-9]*)[\\s]*(-?[0-9]*\.[0-9]*)[\\s]*(-?[0-9]*\.[0-9]*)[\\s]*(-?[0-9]*\.[0-9]*)", 
            file_content
        ), 
        dtype=np.float32
    )[:, 3:]

    # Only extract regex pattern of Energies in Metadata
    energy_pattern = re.compile("(TotEnergy=)(\-[0-9]*\.[0-9]*)")
    liquid_energies = np.asarray([float(energy) for buff, energy in re.findall(energy_pattern, file_content)], dtype=np.float32)

    atom_counter = 0
    with open(test_dir.joinpath('n2p2_fitting/run_pot/liquid_input.data'), 'w') as f:
        for ii_frame, frame in enumerate(liquid_frames):
            f.write('begin\n')

            lattice = frame.cell
            for ii in range(3):
                f.write('lattice\t\t')
                for jj in range(3):
                    f.write(f"{lattice[ii, jj]:.6f}\t")
                f.write('\n')
            
            for atom in frame:
                f.write('atom\t')
                f.write('\t'.join([
                    f'{atom.position[0]:.6f}', f'{atom.position[1]:.6f}', f'{atom.position[2]:.6f}',
                    atom.symbol, f'{0.:.6f}', f'{0.:.6f}',
                    f'{liquid_forces[atom_counter, 0]:.6e} {liquid_forces[atom_counter, 1]:.6e} {liquid_forces[atom_counter, 2]:.6e}\n'
                ]))
                atom_counter+=1

            f.write(f'energy\t\t{liquid_energies[ii_frame]:.6f}\n')
            f.write(f'charge\t\t{0.:.6f}\n')
            f.write('end\n')

if __name__ == "__main__":
    main()