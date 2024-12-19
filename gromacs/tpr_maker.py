import os
dirs = [i  for i in os.listdir() if i.startswith('Conf-')]
dirs = [i  for i in dirs if int(i.split('-')[1])>4]
for dir in dirs:
    file_gro = os.popen(f'ls {dir}/*gro*').read()
    os.system(f'gmx_mpi grompp -c {file_gro} -f {dir}/MD_U.mdp -p {dir}/topol.top -o {dir}/{dir}_PMF.tpr -r {file_gro}')
