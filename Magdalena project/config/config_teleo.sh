# Config parameters
marker="COI" # No majuscules here
raw_path="../../../Metagenomica/Barcoding/5Bosques/COI/HN00186388/"
primer_F="TATAGCATTCCCACGAATAAATAA"
primer_R="TAAACTTCAGGGTGACCAAAAAATCA"
method_demultiplex="" # among: primer_trim / ligation_trim
convert_for_GBIF="TRUE" # values= TRUE or FALSE

# Assignment parameters
method_assignment="dada2" # among: "dada2 / decipher / ecotag" # right now, only dada2 and ecotag are ready

# For ecotag below
path_ecotag_tree="/home/shared/edna/reference_database/2023_06/teleo_custom_embl/customtaxonomy/" # Directory not file # Caution: here it needs to be the path to the tree (NCBI tree ; -t ecotag option not -d)
path_ecotag_fasta="/home/shared/edna/reference_database/2023_06/teleo_custom_embl/db_teleo_custom_and_embl.fasta"

# For dada2 below
path_dada_all="utils/teleo_crabs_dada2_all.fasta"
path_dada_species="utils/teleo_crabs_dada2_species.fasta"

# For decipher below
path_decipher_trained="utils/teleo_trained.rds"

# General
CORES=12
