# Execute primer trimming for one marker 
# V. Marques with code input from D. Kirschner 
# Last updated: 04/12/2023
# Usage: bash trimming.sh config.sh

# Source functions 
source scripts/bash_functions.sh

# Config source
# source config/config_amph.sh
source $1

# Create folder sctructure 
mkdir -p outputs/01_trim/$marker
mkdir -p outputs/02_dada2/$marker
mkdir -p outputs/03_taxo/$marker
output_trim="outputs/01_trim/"${marker}

# For primer trimming 
primer_fwd=$primer_F
primer_fwd_rc=$(reverse_complement "$primer_fwd")

primer_rev=$primer_R
primer_rev_rc=$(reverse_complement "$primer_rev")

# Remove the I files 
#find $raw_path -type f -name '*_I[12]_001.fastq.gz' | xargs rm -v

# ---------------------------------------------------------------------- # 
# DEMULTIPLEXING IF NECESSARY
# ---------------------------------------------------------------------- # 

# TBA

# ---------------------------------------------------------------------- # 
# PRIMER TRIMING
# ---------------------------------------------------------------------- # 

# Run primmer trimming 
#### now primer trimming, only if we dont do it here and not in DADA2
##dont forget to remove from both R1 and R2, thats why we have -a and -A
# Weirdly enough, reverse is not present very often so only cut the forward (and the reverse, see later or cut optionnaly)

samples=$(ls "$raw_path"| grep -v "unknown" | grep -v "Undetermined" |cut -d_ -f1 | sort | uniq)
echo $samples

# For loop is not working properly
# So here is to while loop
while read -r s; do
  echo "Procesando muestra: $s"
  
  # Buscar archivos R1 y R2 para la muestra actual
  r1_in=$(ls "${raw_path}" | grep "^${s}_R1\.fq$")
  r2_in=$(ls "${raw_path}" | grep "^${s}_R2\.fq$")
  
  # Verificar si ambos archivos existen
  if [[ -n $r1_in && -n $r2_in ]]; then
    echo "Archivos encontrados: $r1_in y $r2_in"
    
    # Ejecutar cutadapt para hacer el trimming
    cutadapt \
      -j $CORES \
      -e 0.1 \
      -m 20 \
      -g "${primer_fwd}" \
      -G "${primer_rev}" \
      --untrimmed-output $output_trim/untrimmed_${s}_R1.fq \
      --untrimmed-paired-output $output_trim/untrimmed_${s}_R2.fq \
      -o $output_trim/${s}_R1.fq -p $output_trim/${s}_R2_.fq \
      $raw_path/$r1_in $raw_path/$r2_in >> "logs/log_trim_$s.txt"
  else
    echo "Advertencia: Archivos faltantes para la muestra $s"
  fi
done <<< "$samples"

# For now, remove the untrimmed (to investigate later and save elsewhere in complete version)
rm outputs/01_trim/$marker/untrimmed*

