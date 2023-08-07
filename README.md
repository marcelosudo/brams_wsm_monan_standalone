# brams_wsm_monan_standalone

O código disponível foi enviado pela equipe do Prof. Saulo e possui algumas alterações em virtude de um projeto de pesquisa relacionado a Computação Aproximada realizado com o Orientador Prof. Álvaro Fazenda (Unifesp) e colaboração do Prof. Roberto Souto (LNCC).

O arquivo wsm-test.f90 realiza a preparação dos dados e faz as chamadas para os métodos wsm do BRAMS. Nossas execuções estão configuradas para chamar o wsm5 (module_mp_wms5.f90).

A forma de execução é através do arquivo run.csh (que chama o mk.sh e que por sua vez chama o Makefile para compilação dos arquivos).
A compilação é feita pelo pgf90 com os parâmetros "-acc -Minfo=accel".

As execuções são realizadas no servidor do LNCC, em uma fila de GPU_shared, e os módulos carregados são:
- module load sequana/current
- module load gcc/9.3_sequana
- module load ${NVHPC_DIR}/modulefiles/nvhpc/22.3
Os detalhes estão no arquivo slurm_spack.srm
