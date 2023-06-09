#!/bin/bash
#SBATCH --nodes=1                #Número de Nós
#SBATCH --ntasks=1               #Numero total de tarefas MPI
#SBATCH -p sequana_gpu_shared                  #Fila (partition) a ser utilizada
#SBATCH -J WSM05   		    	#Nome job
#SBATCH --time=00:01:00          #Obrigatório
###SBATCH --exclusive

#Exibe os nós alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

cd  $SLURM_SUBMIT_DIR

## 1) Carrega os módulos sequana
module load sequana/current
#module load python/3.9.12_sequana
module load gcc/9.3_sequana

## 2) Carrega o Spack
workdir=/scratch/mixprecmet/roberto.souto4
version=v0.18.1
spackdir=${workdir}/spack/${version}
. ${spackdir}/share/spack/setup-env.sh
  
export SPACK_USER_CONFIG_PATH=${workdir}/.spack/${version}

#spack load netcdf-fortran@4.5.4%nvhpc@22.3
#spack load netcdf-c@4.8.1%nvhpc@22.3

## 3) Carrega o NVHPC v22.3 da NVIDIA
NVHPC_DIR=$(spack location -i nvhpc@22.3)
module load ${NVHPC_DIR}/modulefiles/nvhpc/22.3


echo -e "\n## Job iniciado em $(date +'%d-%m-%Y as %T') #####################\n"

#export PGI_ACC_TIME=1

#Configura o executavel
#EXEC=./brams-6.0_debug
EXEC=./wsm.x

#exibe informações sobre o executável
#/usr/bin/ldd $EXEC

#srun -n $SLURM_NTASKS $EXEC
#srun $EXEC
#nvprof 
time $EXEC
#mpirun -np $SLURM_NTASKS $EXEC

echo -e "\n## Job finalizado em $(date +'%d-%m-%Y as %T') ###################"
