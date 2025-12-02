- The environment for python (```.py```) scripts will be called via anaconda from the main directory and must be called 'mri_conda_env'. These files can be run with the following command from the scripts folder: ```python -m [script_name]```. 

- The R requirements for the single R script are listed in the main directory and should be installed before running that script.

- The ```.sh``` scripts work with a slurm job manager. 
1) Manually adjust the path held in the ``$WORK`` variable based on your machine's specific paths,
1) Make each ```parallel_*.sh``` file executable by typing ```chmod +x [script_name]```,
2) Execute that file with this command: ```./[script_name]```. This will call the appropriate scripts. 

- Scripts are written specifically for the Glasser-360 atlas.

- The scripts expect a researcher-generated ```Subject_List.txt``` file. This will be a simple text file with subject IDs of those you want to run. If your cluster is similar to ours, there will be a limit on the size of the job batch you submit that will be much lower than the size of the full dataset. Thus, you will most likely need to create several and update them as you go.