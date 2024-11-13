import pathlib
import shutil

import numpy as np

class FileSystemCreator:
    def __init__(self, target_path:pathlib.Path, data_dir=None, weights_path=None):
        self.target_path: pathlib.Path = pathlib.Path(target_path).resolve()
        if data_dir is None:
            data_dir = pathlib.Path(__file__).resolve().parent.joinpath('data')
        self.data_dir = data_dir

        self.mask_dir = self.data_dir.joinpath('mlp_training/mask_files/')

        if weights_path is None:
            weights_path = self.data_dir.joinpath('water_phase_store/lasso_gammas_hartbohr_lambda.npy')
        self.weights_path = weights_path

    def copy_datafile(self):
        shutil.copy(self.mask_dir.joinpath('bc_input.data'), self.target_path.joinpath('input.data'))

    def copy_runfile(self, nfeatures:int, mode:str='nearest'):
        available_modes = {'nearest', 'exact'}
        if not (mode in available_modes):
            raise ValueError(f'Mode {mode} not available, please select from {available_modes}')
        
        if mode == 'nearest':
            available_masks = np.array(sorted([
                int(content.name) for content in self.mask_dir.glob('*') if content.name.isnumeric()
            ]))
            print(f'Found available masks: {available_masks}')
            
            # finding the closest mask could probably be done better with argpartition or something
            nearest_digit = available_masks[np.argmin(np.abs(available_masks-nfeatures))]

            print(f'For copying {nfeatures} using mask with {nearest_digit}')

            shutil.copy(
                self.mask_dir.joinpath(f'{nearest_digit}/run_desc.sh'), 
                self.target_path.joinpath(f'run_desc.sh')
            )
            shutil.copy(
                self.mask_dir.joinpath(f'{nearest_digit}/run_train.sh'),
                self.target_path.joinpath(f'run_train.sh')
            )

        elif mode == 'exact':
            if not self.mask_dir.joinpath(f'{nfeatures}').exists():
                raise FileExistsError(f'Cannot copy bash runscripts for {nfeatures} in exact mode. Use "nearest" to find next available runscripts.')
            
            shutil.copy(self.mask_dir.joinpath(f'{nfeatures}/run_desc.sh'))
            shutil.copy(self.mask_dir.joinpath(f'{nfeatures}/run_train.sh'))
        
        return

    def build_input_file(self, nfeatures, select_gammas=True):

        to_bohr = True
        mult = 1.
        if to_bohr:
            # Unit of eta is currently (1/A^2), this makes it (1/r_bohr^2)
            mult = (1./1.8897259886)**2.

        atom_types = ['H', 'O'] # Need to be in order of atomic numbers

        with open(self.mask_dir.joinpath('input_mask.txt'), 'r') as f:
            input_header = f.read()
            f.close()

        g2_params = np.loadtxt(self.data_dir.joinpath('ice_and_water_data/g2_params.txt'))
        g4_params = np.loadtxt(self.data_dir.joinpath('ice_and_water_data/g4_params.txt'))

        if (self.weights_path.suffix == '.npy') or (self.weights_path.suffix == '.np'):
            lasso_gammas = np.load(self.weights_path)[nfeatures]
        elif (self.weights_path.suffix == '.txt'):
            lasso_gammas = np.loadtxt(self.weights_path)[nfeatures]
        else:
            raise ValueError(f'Cannot parse ending {self.weights_path.suffix} for path {self.weights_path}')
            
        if not select_gammas:
            lasso_gammas = [1.0]*lasso_gammas.shape[0]
        else: # For random generation
            lasso_gammas[:] = 0.
            lasso_gammas[np.random.choice(len(lasso_gammas), (nfeatures), replace=False).astype(dtype=np.int16)] = 1.

        with open(self.target_path.joinpath('input.nn'), 'w') as f:
            f.write(input_header+'\n')

            counter = 0
            for to_at in atom_types:
                counter += 1 # add this for not included G1
                f.write(f'# radial to {to_at}\n')
                for eta, rshift in g2_params:
                    if lasso_gammas[counter]:
                        for from_at in atom_types:
                            f.write(' '.join([
                                'symfunction_short',
                                from_at, '2', to_at, # from, type of function, to
                                f'{mult*eta:.3E}', f'{rshift:.3E}', '12.00', # eta, rshift, cutoff
                                '\n'
                            ]))
                    counter+=1

            for ii_betw, betw_at in enumerate(atom_types):
                for ii_to, to_at in enumerate(atom_types):
                    if ii_to >= ii_betw:
                        f.write(f'# angular to {betw_at} and {to_at}\n')
                        for eta, zeta, lambda_ in g4_params:
                            if lasso_gammas[counter]:
                                for from_at in atom_types:
                                    f.write(' '.join([
                                        'symfunction_short',
                                        from_at, '3', betw_at, to_at, # from, type of function, middle atom, to
                                        f'{mult*eta:.3E}', f'{lambda_:.3E}', f'{zeta:.3E}', '12.00', # eta, lambda, zeta, cutoff
                                        '\n'
                                    ]))
                            counter+=1
                
        print(f"Built input file from {counter} symmetry functions using {nfeatures}.")

def main():
    nfeats = [10, 18, 25, 38, 50, 176]
    data_path = pathlib.Path(__file__).resolve().parent.joinpath('data')
    for nfeat in nfeats:
        target_path = data_path.joinpath(f'mlp_training/pot_acsf_{nfeat}')
        if not target_path.exists():
            target_path.mkdir()
        fs_creator = FileSystemCreator(
            target_path=target_path, data_dir=data_path
        )
        fs_creator.copy_datafile()
        fs_creator.copy_runfile(nfeat)
        fs_creator.build_input_file(nfeatures=nfeat)

if __name__ == '__main__':
    main()