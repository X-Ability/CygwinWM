### Installing Miniconda

1. Download the Miniconda3 Windows 64-bit installer from [here](https://docs.conda.io/en/latest/miniconda.html).

2. Install Miniconda3 with the following options:
   - Install for: Just Me
   - Do not create Start menu shortcuts or any additional options.

### Installing Pymatgen

1. Launch the Command Prompt.

2. Execute the following command to activate Miniconda:

   ```shell
   miniconda3\condabin\activate.bat
   ```

3. Install Pymatgen using Conda Forge:

   ```shell
   conda install --channel conda-forge pymatgen
   ```

4. Clean up unused packages:

   ```shell
   conda clean --all
   ```

5. You can optionally delete the `miniconda3\pkgs` folder to free up disk space.


