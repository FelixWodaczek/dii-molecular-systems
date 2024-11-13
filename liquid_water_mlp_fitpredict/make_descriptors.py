import pathlib

import numpy as np
from ase.io import read as ase_read
from dscribe.descriptors import ACSF, SOAP

def get_descriptors(frames, featuriser):
    n_atoms = np.sum(np.asarray([len(frame) for frame in frames], dtype=np.int16))
    average_desc = np.zeros((len(frames), featuriser.get_number_of_features()), dtype=np.float32)
    singleo_desc = np.zeros((n_atoms, featuriser.get_number_of_features()), dtype=np.float32) # singleo_desc = np.zeros_like(average_desc)
    counter = 0
    for ii_frame, frame in enumerate(frames):
        # Get index of any oxygen in system
        symbols = frame.get_chemical_symbols()
        # first_o_ind = next(ii_symb for ii_symb, symbol in enumerate(symbols) if symbol=='O')

        frame_desc = featuriser.create(frame, n_jobs=8)
        average_desc[ii_frame, :] = np.mean(frame_desc, axis=0)
        singleo_desc[counter:counter+len(frame), :] = frame_desc # frame_desc[first_o_ind, :]
        counter+=len(frame)
        if not (ii_frame%10):
            print(ii_frame)

    return average_desc, singleo_desc

def main():
    data_dir = pathlib.Path(__file__).resolve().parent.joinpath('data')

    # If this is set to true, eta will be converted so it will be in the given range when converted to bohr radii
    make_in_bohr = True
    if make_in_bohr:
        fact = 1.8897259886**2.

    ice_and_water_dir = data_dir.joinpath('ice_and_water_data')
    target_xyz = ice_and_water_dir.joinpath('H2O-n6-l6-c6.0-g0.3-pca-d10.xyz')
    frames = ase_read(target_xyz, index=':')

    # ACSF targets
    etas = np.logspace(-3, 0.5, 6)*fact
    zetas = np.array([1, 2, 3, 4], dtype=np.float32)
    lambdas = np.array([-1, 1], dtype=np.float32)

    g2_params = np.array([[eta*fact, 0.] for eta in np.logspace(-3, 0.5, 15)], dtype=np.float32)
    g4_params = np.meshgrid(etas, zetas, lambdas)
    g4_params = np.array(g4_params).reshape((3, -1)).T
    np.savetxt(ice_and_water_dir.joinpath("g2_params.txt"), g2_params, header='eta log -3, 0.5 15')
    np.savetxt(ice_and_water_dir.joinpath("g4_params.txt"), g4_params, header='eta log -3 -1 6, zeta [1, 2, 3, 4], lambda [-1, 1]')

    acsf = ACSF(
        species=[1, 8],
        r_cut=6.0,
        g2_params=g2_params,
        g4_params=g4_params,
        periodic=True
    )

    soap = SOAP(
        n_max=6, l_max=6,
        r_cut=6.0, species=[1, 8], 
        sigma=0.3,
        periodic=True,
    )

    print("Generating ACSFs")
    average_acsf, singleo_acsf = get_descriptors(frames, acsf)
    print("Generating SOAPs")
    average_soap, singleo_soap = get_descriptors(frames, soap)

    # Save the descriptors
    storage_dir = data_dir.joinpath('water_phase_store')
    np.save(storage_dir.joinpath('average_acsf_rcut6_gridsearch_bohr_lambda.npy'), average_acsf)
    np.save(storage_dir.joinpath('singleatom_acsf_rcut6_gridsearch_bohr_lambda.npy'), singleo_acsf)
    np.save(storage_dir.joinpath('average_soap_rcut6_nmax6_lmax6_sigma03.npy'), average_soap)
    np.save(storage_dir.joinpath('singleatom_soap_rcut6_nmax6_lmax6_sigma03.npy'), singleo_soap)

if __name__ == '__main__':
    main()