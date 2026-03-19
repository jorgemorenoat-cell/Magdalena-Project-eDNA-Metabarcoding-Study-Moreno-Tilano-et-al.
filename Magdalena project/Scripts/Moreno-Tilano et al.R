###############################Process of filtering, merging, and creating ASV and abundance tables#################

# Put in the arguments to handle the paths
library(devtools)
library(dada2); packageVersion("dada2")
library(Biostrings)
library(DECIPHER)
library(dplyr)
library(patchwork)
library(ggplot2)


# Debug
marker <- "Vert01" #Switch to teleo to process Teleo data

# ---------------------------------- # 
# FUNCTIONS 

# Convert df to fasta file
df_to_fasta <- function(file, output_file_path){
  fa <- character(2 * nrow(file))
  fa[c(TRUE, FALSE)] = sprintf("> %s", file[,1])
  fa[c(FALSE, TRUE)] = as.character(file[,2])
  writeLines(fa, output_file_path)
}

extract_pattern_samplename <- function(input_string) {
  sub(".*/([^/]+)\\.[Rr][12].*", "\\1", input_string)
}

# END of FUNCTIONS
# ------------------------------ #
# Create directory
# Create directory
dir.create(paste0("outputs/Barrio Abajo (Maylin project)/", marker, "/QC/"), recursive = T)

# File parsing (modificar)
pathFR <- paste0("~/outputs/01_trim/teleo/") # CHANGE ME to the directory containing your demultiplexed  fastqs

list.files(pathFR)

filtpathFR <- file.path(pathFR, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
fastqFs <- sort(list.files(pathFR, pattern="R1.fastq"))
fastqRs <- sort(list.files(pathFR, pattern="R2.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Explore quality 
fnFs <- sort(list.files(pathFR, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(pathFR, pattern="R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

# Plot quality for QC checks
for (i in 1:length(fnFs)){
  name_sample <- extract_pattern_samplename(fnFs[i])
  p_i <- plotQualityProfile(fnFs[i]) + plotQualityProfile(fnRs[i])
  ggsave(paste0("outputs/02_dada2/", marker, "/QC/", name_sample, "_dada2QC.png"), p_i)
}

# File parsing
filtpathFR <- paste0("outputs/01_trim/",marker,"/filtered")

# Filter
filterAndTrim(fwd=file.path(pathFR, fastqFs), filt=file.path(filtpathFR, fastqFs),
              rev=file.path(pathFR, fastqRs), filt.rev=file.path(filtpathFR, fastqRs),
              maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=FALSE, verbose=TRUE, multithread=TRUE)

# FIltered reads
filtFs <- list.files(filtpathFR, pattern="R1.fastq", full.names = TRUE)
filtRs <- list.files(filtpathFR, pattern="R2.fastq", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "\\."), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "\\."), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = 20)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0("outputs/02_dada2/", marker, "/seqtab_pcr_merg.rds"))
# Merge multiple runs (if necessary)
st1 <- readRDS(paste0("outputs/02_dada2/", marker, "/seqtab_pcr_merg.rds"))
# Remove chimeras
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
saveRDS(seqtab,  paste0("outputs/02_dada2/", marker, "/seqtab_clean_pcr_merg.rds")) # CHANGE ME to where you want sequence table saved
saveRDS(t(seqtab),  paste0("outputs/02_dada2/", marker, "/seqtab_clean_t_pcr_merg.rds")) # CHANGE ME to where you want sequence table saved

write.csv(t(seqtab), paste0("outputs/02_dada2/", marker, "/seqtab_clean_t_pcr_merg.csv"), row.names = TRUE)

########################## Assign taxonomy and group ASVs by taxonomy #############################################
library(dada2); packageVersion("dada2")
library(Biostrings)
library(dplyr)

# Debug
marker <- "vert01" #Switch to teleo to process Teleo data
#Reference database path
path_dada_all <- "utils/Dada_mitofish_midori_mitofish_vert_taxonomy" # Switch to "Dada_mitofish_teleo_taxonomy" to process Teleo data
path_dada_species <- "utils/Dada_mitofish_midori_mitofish_vert_species" # Switch to "Dada_mitofish_teleo_species" to process Teleo data

cat(paste0("marker is ", marker, " \n the path dada all is ", path_dada_all, " \n and the path dada species is ", path_dada_species))

# Open seqtab - abundance table output of dada2
seqtab <- readRDS(paste0("outputs/02_dada2/", marker, "/seqtab_clean_pcr_merg.rds"))

#Taxonomic assignment (Kmers approach) with bootstrap values
taxa <- assignTaxonomy(seqtab, path_dada_all, multithread = TRUE, tryRC = TRUE, outputBootstraps = TRUE)
class(taxa)
# Add descriptive names to the elements of the resulting object
taxonomy <- taxa$tax
class(taxonomy)
bootstrap <- taxa$boot  # bootstrap Values

# Add second species assignment (Exact match approach)
taxa_plus <- addSpecies(taxonomy, path_dada_species, verbose = TRUE, allowMultiple = TRUE)

# Convert to a data frame
taxa_plus <- as.data.frame(taxa_plus)
taxa_plus$sequence <- rownames(taxa_plus)

# Include bootstrap values in additional columns
bootstrap_df <- as.data.frame(bootstrap)
colnames(bootstrap_df) <- paste0("bootstrap_", colnames(bootstrap_df))

# Join both columns
taxa_final <- cbind(taxa_plus, bootstrap_df)
#only for assignments with mitofish only 
colnames(taxa_final)[8]<-"siguientes_hit_especie"
#teleo 
write.csv(taxa_final, paste0("outputs/02_dada2/", marker, "/", "taxonomy_dada2RDP_2022_pcr_merg_database_mitofish_midori_10_04_25.csv"), row.names=F)
#vert01 
write.csv(taxa_final, paste0("outputs/02_dada2/", marker,"/", "taxonomy_dada2RDP_2022_vert01_database_mido_mito_31_04_24.csv"), row.names=F)

# Now link it to abundance table
# Open data 
samples <- read.csv(paste0("outputs/02_dada2/", marker,"/seqtab_clean_t__pcr_merg.csv"))

colnames(samples)[1] <- "sequence"
# Now link 
dada2_RDP_taxo_complete <- dplyr::left_join(taxa_final, samples)
# join ASV by taxonomy 
#When performing ASV analysis, many of these ASVs correspond to the same taxon; therefore, these ASVs are grouped together based on their identical taxonomy.
#the sequences for these ASVs are stored in the “sequence” column; the bootstrap percentages are averaged, and the number of reads per site is summed.
library(dplyr)
colnames(dada2_RDP_taxo_complete)
datos_consolidados <- dada2_RDP_taxo_complete %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(
    sequence = paste(sequence, collapse = "-"),  # join sequences with "-"
    siguientes_hit_especie =paste(siguientes_hit_especie, collapse = "-"),
    across(starts_with("SP"), sum, na.rm = TRUE), # Sum values by sites
    across(starts_with("bootstrap"), ~mean(.x, na.rm = TRUE)) # Calculates the average of bootstrap values
  ) %>%
  ungroup()
#filter by columns
columnas_taxonomia <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
columnas_bootstrap <- paste0("bootstrap_", columnas_taxonomia)
columnas_muestras <- setdiff(names(datos_consolidados), 
                             c(columnas_taxonomia, columnas_bootstrap, "siguientes_hit_especie", "sequence"))

# Automatically reorder columns
datos_consolidados_f <- datos_consolidados[, c(
  rbind(columnas_taxonomia, columnas_bootstrap),  # Alternate between taxonomy and bootstrap
  "siguientes_hit_especie",                       
  "sequence",                                     # sequences columns
  columnas_muestras                               # samples columns
)]

#Teleo
write.csv(datos_consolidados_f, paste0("~/outputs/teleo/", marker, "_dada2RDP_2022_pcr_merg_database_mitofish_31_04_24.csv"), row.names=FALSE)
#vert01
write.csv(datos_consolidados_f, paste0("~/outputs/Vert01/", marker,"_dada2RDP_2022_vert01_database_mido_mito_09_04_31.csv"), row.names=F)

#####################Combine the lists from Teleo and Vert and filter the data by number of reads#######################################
#load both data sets
Teleo<- read.csv("~/outputs/teleo/teleo_dada2RDP_2022_pcr_merg_database_mitofish_31_04_24.csv")
vert01<- read.csv("~/outputs/Vert01/Vert01_dada2RDP_2022_vert01_database_mido_mito_09_04_25.csv")

# Rename the capsules by location and merge the duplicates
New_names<-readxl::read_xlsx("~/utils/Abreviaturasf.xlsx", sheet = "Abrevi")

# Create a vector of current and new names
code_map <- setNames(New_names$Site, New_names$Capsule_code)

# Rename columns for Teleo
names(Teleo) <- ifelse(names(Teleo) %in% names(code_map),
                       code_map[names(Teleo)],
                       names(Teleo))
# Rename columns for Vert01
names(vert01) <- ifelse(names(vert01) %in% names(code_map),
                        code_map[names(vert01)],
                        names(vert01))

# Code for combining responses by Teleo sites
for (i in 1:14) {
  col1 <- paste0("S", i, "-1")
  col2 <- paste0("S", i, "-2")
  new_col <- paste0("S", i)
  
  # Set missing columns to NA if they do not exist
  if (!col1 %in% names(Teleo)) Teleo[[col1]] <- NA
  if (!col2 %in% names(Teleo)) Teleo[[col2]] <- NA
  
  # Create the new combined column by adding the existing columns together
  Teleo[[new_col]] <- rowSums(
    cbind(Teleo[[col1]], Teleo[[col2]]), 
    na.rm = TRUE
  )
}

# Code for grouping replicas by site Vert01
for (i in 1:14) {
  col1 <- paste0("S", i, "-1")
  col2 <- paste0("S", i, "-2")
  new_col <- paste0("S", i)
  
  # Inicializar columnas faltantes con NA si no existen
  if (!col1 %in% names(vert01)) vert01[[col1]] <- NA
  if (!col2 %in% names(vert01)) vert01[[col2]] <- NA
  
  # Create the new combined column by adding the existing ones together
  vert01[[new_col]] <- rowSums(
    cbind(vert01[[col1]], vert01[[col2]]), 
    na.rm = TRUE
  )
}
# Remove columns from replicas and keep only those from linked sites in Teleo
cols_to_remove <- grep("-1$|-2$|-3$|-4$", names(Teleo), value = TRUE)
Teleo <- Teleo[, !(names(Teleo) %in% cols_to_remove)]

# Remove columns from replicas and keep only those from linked sites in Vert01
cols_to_remove <- grep("-1$|-2$|-3$|-4$", names(vert01), value = TRUE)
vert01 <- vert01[, !(names(vert01) %in% cols_to_remove)]



# Define the key columns
key_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
bootstrap_cols <- paste0("bootstrap_", key_cols)
seq_col <- "sequence"
hit_col <- "siguientes_hit_especie"

# Detect site columns (all columns that are not primary, bootstrap, or special columns)
site_cols <- setdiff(
  colnames(vert01),
  c(key_cols, bootstrap_cols, seq_col, hit_col)
)

# Perform a full join between the two data frames
Mag_data <- full_join(Teleo, vert01, by = key_cols, suffix = c(".1", ".2"))

# Average bootstrap values
for (col in bootstrap_cols) {
  col1 <- paste0(col, ".1")
  col2 <- paste0(col, ".2")
  Mag_data[[col]] <- rowMeans(cbind(Mag_data[[col1]], Mag_data[[col2]]), na.rm = TRUE)
}

# Merge the "sequence" and "following_hit_especie" columns
Mag_data[[seq_col]] <- dplyr::case_when(
  is.na(Mag_data[[paste0(seq_col, ".1")]]) ~ Mag_data[[paste0(seq_col, ".2")]],
  is.na(Mag_data[[paste0(seq_col, ".2")]]) ~ Mag_data[[paste0(seq_col, ".1")]],
  TRUE ~ paste(Mag_data[[paste0(seq_col, ".1")]], Mag_data[[paste0(seq_col, ".2")]], sep = " - ")
)

Mag_data[[hit_col]] <- dplyr::case_when(
  is.na(Mag_data[[paste0(hit_col, ".1")]]) ~ Mag_data[[paste0(hit_col, ".2")]],
  is.na(Mag_data[[paste0(hit_col, ".2")]]) ~ Mag_data[[paste0(hit_col, ".1")]],
  TRUE ~ paste(Mag_data[[paste0(hit_col, ".1")]], Mag_data[[paste0(hit_col, ".2")]], sep = " - ")
)

# sum reads per site
for (site in site_cols) {
  col1 <- paste0(site, ".1")
  col2 <- paste0(site, ".2")
  Mag_data[[site]] <- rowSums(cbind(
    ifelse(is.na(Mag_data[[col1]]), 0, Mag_data[[col1]]),
    ifelse(is.na(Mag_data[[col2]]), 0, Mag_data[[col2]])
  ), na.rm = TRUE)
}

# Reorder the data frame to the desired structure
Mag_data_final <- Mag_data[, c(
  as.vector(rbind(key_cols, bootstrap_cols)), # Alternate between taxonomy and bootstrap
  "siguientes_hit_especie",                               
  "sequence",                                               # Sequence column
  site_cols                                         # Sample columns
)]

#filter the results

#Return 0 records with fewer than 10 reads
Mag_data_final_filt <- Mag_data_final %>%
  mutate(across(all_of(site_cols), ~ ifelse(.x >= 1 & .x <= 10, 0, .x)))

#Calculate the total number of ASVs detected at each site
fila_suma <- Mag_data_final_filt %>%
  summarise(across(all_of(site_cols), ~ sum(.x > 0, na.rm = TRUE))) %>%
  mutate(Kingdom = "TOTAL", .before = 1)  # Add an ID for the new row

# Identify capsules that have 5 or fewer ASVs detected in total
columnas_a_eliminar <- names(fila_suma)[fila_suma[1, ] <= 5]

#calculate the total sum of readings by ASV
suma_lecturas <- Mag_data_final_filt %>%
  mutate(total_lecturas = rowSums(across(all_of(site_cols)), na.rm = TRUE)) %>%
  pull(total_lecturas)  # Extraer solo el vector de valores

# Remove ASVs with fewer than 15 reads after the previous filtering
Mag_data_final_filt <- Mag_data_final_filt %>%
  filter(suma_lecturas > 15)
# Remove the “Kingdom” and “Kingdom_boostr” columns since they contain no data 
Mag_data_final_filt<-Mag_data_final_filt[,-(1:2)]

Mag_data_final_filt_ord <- Mag_data_final_filt %>%
  arrange(Class, Order, Family, Genus, Species)

#Save curated list by sampling sites
write.csv(Mag_data_final_filt_ord, paste0("~/outputs/Magdalena_data_2022_teleo_vert_dada2.csv"), row.names=FALSE)

######################Download the vertebrate records reported in the study area of GIBF##################
library(rgbif)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(purrr)
#load Colombia map
colombia <- ne_countries(scale = "medium", country = "Colombia", returnclass = "sf")
#load the major rivers layer and filter for the Magdalena River
rios <- ne_download(scale = "medium", type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
#select only the Magdalena River
rio_magdalena <- rios[rios$name == "Magdalena", ]
plot(rio_magdalena)
#Create a horizontal line at the desired latitude for cropping; this latitude was chosen based on the southernmost point, which is Gamarra.
poligono_corte <- st_as_sfc(st_bbox(c(ymin = 8.185939, ymax = 11.1, 
                                      xmin = -76.6, xmax = -73.5), 
                                    crs = st_crs(rio_magdalena)))
# Trim the river using the trim polygon
Bajo_magdalena <- st_intersection(rio_magdalena, poligono_corte)
# Display to confirm the cut
ggplot() +
  geom_sf(data = colombia, fill = "lightgray", color = "black") +  # Mapa de Colombia
  geom_sf(data = rio_magdalena, color = "blue", linewidth = 1) +
  geom_sf(data = Bajo_magdalena, color = "red", linewidth = 1.5) +
  geom_sf(data = poligono_corte, color = "darkred", fill = NA, linewidth = 1, linetype = "dashed") +
  ggtitle("Area de estudio articulo") +
  theme_minimal()


# Projecting onto a metric coordinate system
bajo_magdalena_m <- st_transform(Bajo_magdalena, 3116)
ggplot()+
  geom_sf(data= Bajo_magdalena)

# Create a 30 km (30,000-meter) buffer zone
buffer_magdalena <- st_buffer(bajo_magdalena_m, dist = 30000)
# Re-project the buffer back to WGS 84 (EPSG:4326)
buffer_magdalena_wgs84 <- st_transform(buffer_magdalena, crs = 4326)
# Trim the river using the trim polygon
Bajo_magdalena2 <- st_intersection(buffer_magdalena_wgs84, poligono_corte)
Bajo_magdalena2
#View the buffer on a map
ggplot() +
  geom_sf(data = colombia, fill = "gray80", color = "black") +       # Mapa base de Colombia
  geom_sf(data = Bajo_magdalena2, fill = "blue", alpha = 0.4) +     # Buffer del río Magdalena
  geom_sf(data = Bajo_magdalena, color = "darkblue", size = 1) +     # Línea del río Magdalena
  coord_sf(xlim = c(-76.6, -72.5), ylim = c(7, 11.5), expand = FALSE) +
  theme_minimal() +
  ggtitle("Polígono buffer de 30 km a cada lado del río Magdalena")
# Re-project the buffer back to WGS 84 (EPSG:4326)
Bajo_magdalena2 <- st_transform(Bajo_magdalena2, crs = 4326)
library(sf)

# Save the polygon as a shapefile
st_write(
  Bajo_magdalena2,
  "../Bajo_magdalena_buffer.shp",   
  driver = "ESRI Shapefile",
  delete_layer = TRUE           
)
st_write(
  Bajo_magdalena,
  "../Bajo_magdalena_linea.shp",   
  driver = "ESRI Shapefile",
  delete_layer = TRUE
)
#Set a password to download GIBF records 
Chordata_key <- name_backbone(name = "Chordata")$usageKey
#Generate a polygon in WKT format to download data
polygon_wkt <- st_as_text(st_union(st_geometry(buffer_magdalena_wgs84)))

# Download records with improved quality criteria for preserved and machine observation only
registros_p_m <- occ_search(
  taxonKey = Chordata_key,
  geometry = polygon_wkt,             # Area of interest
  hasCoordinate = TRUE,               # Georeferenced records only
  eventDate = "2000,2025",            # time range
  basisOfRecord = c("PRESERVED_SPECIMEN", "MACHINE_OBSERVATION"),
  occurrenceStatus = "PRESENT",       # Avoid missing entries
  limit = 15000                       # Adjust the limit as needed
)


# Create a new data frame with selected columns
Vertebrados_p <- registros_p_m$PRESERVED_SPECIMEN$data[, c("phylum", "order", 
                                                           "family", "genus", "species",
                                                           "acceptedScientificName","taxonomicStatus",
                                                           "iucnRedListCategory","decimalLatitude","decimalLongitude", "issues")]
Vertebrados_M <- registros_p_m$MACHINE_OBSERVATION$data[, c("phylum","class", "order", 
                                                            "family", "genus", "species",
                                                            "acceptedScientificName","taxonomicStatus",
                                                            "iucnRedListCategory","decimalLatitude","decimalLongitude", "issues")]
#Since the records for preserved individuals don't have a “class” field, we need to add a “class” column to avoid losing all the information from the rest of the dataset
Vertebrados_p$class<-""
#sort columns
Vertebrados_p<-Vertebrados_p[,c("phylum","class", "order", 
                                "family", "genus", "species",
                                "acceptedScientificName","taxonomicStatus",
                                "iucnRedListCategory","decimalLatitude","decimalLongitude","issues")]
#Join the two data frames
Vertebrados_p_m<- rbind(Vertebrados_p,Vertebrados_M)

#Verify the reported issues for the obtained records
table(Vertebrados_p_m$issues)
# Filter records by removing those containing ‘osu’, ‘muluriiv’, ‘inmano’, or 'osiic'
registros_filtrados <- Vertebrados_p_m[!grepl("osu|muluriiv|osiic", Vertebrados_p_m$issues), ]
# Get a list of unique species
Lis_Vert_unicoss_p_m <- registros_filtrados[!duplicated(registros_filtrados$acceptedScientificName), ]

##Perform the same procedure but download observation records
#Since there are so many records and we can only download 17,000 at a time, we've divided the study area into two parts
# Calculate the average latitude
lat_media <- (8.185939 + 11.5) / 2  # Resultado: 9.84297

# Polygon for the southern half
poligono_sur <- st_as_sfc(st_bbox(c(ymin = 8.185939, ymax = lat_media, 
                                    xmin = -76.6, xmax = -73.5), 
                                  crs = st_crs(rio_magdalena)))

# Polygon for the Northern half
poligono_norte <- st_as_sfc(st_bbox(c(ymin = lat_media, ymax = 11.5, 
                                      xmin = -76.6, xmax = -73.5), 
                                    crs = st_crs(rio_magdalena)))

# Dividing the Magdalena River into two parts
Bajo_magdalena1 <- st_intersection(rio_magdalena, poligono_sur)
Bajo_magdalena2 <- st_intersection(rio_magdalena, poligono_norte)

# Projecting onto a metric coordinate system
bajo_magdalena_m <- st_transform(Bajo_magdalena1, 3116)
bajo_magdalena_m1 <- st_transform(Bajo_magdalena2, 3116)
# Create a 30 km (30,000-meter) buffer zone
buffer_magdalena <- st_buffer(bajo_magdalena_m, dist = 30000)
buffer_magdalena1 <- st_buffer(bajo_magdalena_m1, dist = 30000)
# Re-project the buffer back to WGS 84 (EPSG:4326)
buffer_magdalena_wgs84 <- st_transform(buffer_magdalena, crs = 4326)
buffer_magdalena_wgs841 <- st_transform(buffer_magdalena1, crs = 4326)

#Set a password to download GIBF records 
Chordata_key <- name_backbone(name = "Chordata")$usageKey
#Generate a polygon in WKT format to download data
polygon_wkt <- st_as_text(st_union(st_geometry(buffer_magdalena_wgs84)))
polygon_wkt1 <- st_as_text(st_union(st_geometry(buffer_magdalena_wgs841)))

# Download records with improved quality criteria for first coordinates
registros_h_p_m1 <- occ_search(
  taxonKey = Chordata_key,
  geometry = polygon_wkt1,            
  hasCoordinate = TRUE,               
  eventDate = "2000,2025",            
  basisOfRecord = c("HUMAN_OBSERVATION"),
  occurrenceStatus = "PRESENT",      
  limit = 17000                      
)

# Download records with improved quality criteria for secondary coordinates
registros_h_p_m <- occ_search(
  taxonKey = Chordata_key,
  geometry = polygon_wkt,             
  hasCoordinate = TRUE,              
  eventDate = "2000,2025",            
  basisOfRecord = c("HUMAN_OBSERVATION"),
  occurrenceStatus = "PRESENT",       
  limit = 17000                       
)
#select the columns of interest
Listado_verteb_hm<- registros_h_p_m$data[, c("phylum","class", "order", 
                                             "family", "genus", "species",
                                             "acceptedScientificName","taxonomicStatus",
                                             "iucnRedListCategory","decimalLatitude","decimalLongitude", "issues")]
Listado_verteb_hm1<- registros_h_p_m1$data[, c("phylum","class", "order", 
                                               "family", "genus", "species",
                                               "acceptedScientificName","taxonomicStatus",
                                               "iucnRedListCategory","decimalLatitude","decimalLongitude", "issues")]
#combine the lists 
listado_vert_all<-rbind(Listado_verteb_hm,Listado_verteb_hm1)

# Get a list of unique species
Lis_Vert_unicoss_hm <- listado_vert_all[!duplicated(listado_vert_all$acceptedScientificName), ]

list_vert_all<- rbind(Lis_Vert_unicoss_hm, Lis_Vert_unicoss_p_m)
# Get a list of unique species again
Lis_Vert_unicoss_hm_f<- list_vert_all[!duplicated(list_vert_all$acceptedScientificName), ]


# Save the table to a CSV file
write.csv(Lis_Vert_unicoss_hm_f, "~/outputs/Listado de vertebrados de todas las fuentes  del area de estudio2.csv", row.names = FALSE)

#################################Figure S3#########################################
library(patchwork) 
library(cowplot)
library(ggsci)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(readODS)
library(scales)
library(readr)
library(stringdist)
library(fuzzyjoin)
library(tibble)
Asv_sub_teleo<-read.csv("~/outputs/Taxa_counts.csv")
Asv_sub_vert<-read.csv("~/outputs/Taxa_counts.csv")

# Convert “Subset” to a number (extracting the value after “sub_”)
convertir_subset_num <- function(df, col = "Subset") {
  df$Subset_num <- as.numeric(sub("sub_", "", df[[col]]))
  return(df)
}

Asv_sub_teleo<-convertir_subset_num(Asv_sub_teleo)
Asv_sub_vert<-convertir_subset_num(Asv_sub_vert)

#Function for adjusting saturation curves
ajustar_saturacion <- function(df) {
  # Logarithmic scaling
  fit_log <- try(
    nls(Num_Taxa ~ a * log(Subset_num) + b,
        data = df, start = list(a = max(df$Num_Taxa)/2, b = min(df$Num_Taxa))),
    silent = TRUE
  )
  
  # Exponential asymptotic adjustment (self-starting saturation)
  fit_asym <- try(
    nls(Num_Taxa ~ y0 + a * (1 - exp(-exp(b) * Subset_num)),
        data = df, start = list(y0 = min(df$Num_Taxa), a = max(df$Num_Taxa) - min(df$Num_Taxa), b = 0)),
    silent = TRUE
  )
  
  # Compare AIC
  modelos <- list()
  if (!inherits(fit_log, "try-error")) {
    modelos$log <- list(fit = fit_log, AIC = AIC(fit_log))
  }
  if (!inherits(fit_asym, "try-error")) {
    modelos$asym <- list(fit = fit_asym, AIC = AIC(fit_asym))
  }
  
  # Calculate ΔAIC
  AICs <- sapply(modelos, function(x) x$AIC)
  deltaAIC <- AICs - min(AICs)
  
  # Best model
  mejor <- names(which.min(AICs))
  
  # Predictions up to 1.5
  new_data <- data.frame(Subset_num = seq(0, 1.5, length.out = 200))
  new_data$pred <- predict(modelos[[mejor]]$fit, newdata = new_data)
  
  return(list(
    modelo_ganador = mejor,
    AICs = AICs,
    deltaAIC = deltaAIC,
    fit_log = fit_log,
    fit_asym = fit_asym,
    predicciones = new_data
  ))
}



#  Function to retrieve the formula as text from the template 
extraer_formula <- function(res) {
  if (res$modelo_ganador == "log") {
    # Extract parameters
    coef_log <- coef(res$fit_log)
    a <- round(coef_log["a"], 3)
    b <- round(coef_log["b"], 3)
    # Create formula  y = a*ln(x) + b
    formula_texto <- paste0("y = ", a, " * ln(x) + ", b)
    modelo_tipo <- "Logarithmic"
  } else if (res$modelo_ganador == "asym") {
    # Extract parameters
    coef_asym <- coef(res$fit_asym)
    y0 <- round(coef_asym["y0"], 3)
    a <- round(coef_asym["a"], 3)
    b <- round(coef_asym["b"], 3)
    # Fórmula tipo y = y0 + a*(1-exp(-exp(b)*x))
    formula_texto <- paste0("y = ", y0, " + ", a, " * (1 - exp(-exp(", b, ") * x))")
    modelo_tipo <- "Asymptotic / Exponential"
  } else {
    formula_texto <- NA
    modelo_tipo <- NA
  }
  
  return(list(
    formula = formula_texto,
    tipo = modelo_tipo
  ))
}
res_teleo <- ajustar_saturacion(Asv_sub_teleo)

# Extract formula and model type
formula_teleo <- extraer_formula(res_teleo)
modelo_ingles <- formula_teleo$tipo
formula_texto <- formula_teleo$formula

# Calculate delta AIC
delta_texto <- round(max(res_teleo$deltaAIC), 2)

# Define the gray area
x_pred_inicio <- 1
x_pred_fin <- max(res_teleo$predicciones$Subset_num)

# Chart
Teleo <- ggplot(Asv_sub_teleo, aes(x = Subset_num, y = Num_Taxa)) +
  # Prediction Gray Area
  geom_rect(aes(xmin = x_pred_inicio, xmax = x_pred_fin,
                ymin = -Inf, ymax = Inf),
            fill = "gray90", alpha = 0.5, inherit.aes = FALSE) +
  # Data points
  geom_point(color = "blue4", size = 2) +
  # Prediction line
  geom_line(data = res_teleo$predicciones,
            aes(x = Subset_num, y = pred),
            linetype = "dotted", color = "gray30", size = 1) +
  # Headings and labels
  labs(
    x = "Proportion of reads used",
    y = "Number of taxa recovered",
    title = "Saturation Teleo"
  ) +
  # Note with model, AIC delta, and formula
  annotate("text", x = 1, y = 0,
           label = paste("Best model:", modelo_ingles,
                         "\n", formula_texto,
                         "\nΔAIC:", delta_texto),
           hjust = 0, size = 3.5, color = "black") +
  # X-axis with all intervals
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.2), expand = c(0.01,0)) +
  # Adjusting the Y-axis
  scale_y_continuous(expand = c(0.2,0.2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, max(Asv_sub_teleo$Num_Taxa)))
###vert 01
res_vert01 <- ajustar_saturacion(Asv_sub_vert)
# Extract formula and model type
formula_vert01 <- extraer_formula(res_vert01)
modelo_ingles <- formula_vert01$tipo
formula_texto <- formula_vert01$formula

# Calculate AIC delta
delta_texto <- round(max(res_vert01$deltaAIC), 2)

# Define the gray area
x_pred_inicio <- 1
x_pred_fin <- max(res_vert01$predicciones$Subset_num)

# graph
vert01 <- ggplot(Asv_sub_vert, aes(x = Subset_num, y = Num_Taxa)) +
  # Prediction Gray Area
  geom_rect(aes(xmin = x_pred_inicio, xmax = x_pred_fin,
                ymin = -Inf, ymax = Inf),
            fill = "gray90", alpha = 0.5, inherit.aes = FALSE) +
  # Data points
  geom_point(color = "blue4", size = 2) +
  # Forecast line
  geom_line(data = res_vert01$predicciones,
            aes(x = Subset_num, y = pred),
            linetype = "dotted", color = "gray30", size = 1) +
  # titles y labels
  labs(
    x = "Proportion of reads used",
    y = "Number of taxa recovered",
    title = "Saturation Vert01"
  ) +
  # Score with model, AIC delta, and formula
  annotate("text", x = 0.7, y = 0,
           label = paste("Best model:", modelo_ingles,
                         "\n", formula_texto,
                         "\nΔAIC:", delta_texto),
           hjust = 0, size = 3.5, color = "black") +
  # X-axis with all intervals
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.2), expand = c(0.01,0)) +
  # Y-axis adjustment
  scale_y_continuous(expand = c(0.2,0.2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, max(Asv_sub_vert$Num_Taxa)))


combined_plot <- plot_grid( vert01,Teleo, 
                            ncol = 2,           # <- poner 2 columnas para que queden lado a lado
                            labels = c("a)", "b)"), 
                            label_fontfamily = "serif", 
                            label_fontface = "plain", 
                            label_size = 20,
                            align = "h"         # <- opcional, para alinear horizontalmente
)
combined_plot
###save  graph
png(filename = "~/Figure S2.png",
    
    width = 4000,height = 2000 ,units = "px", res = 300)
combined_plot
dev.off()
png(filename = "~/Figure S2.TIFF",
    
    width = 4000,height = 2000 ,units = "px", res = 300)
combined_plot
dev.off()

###################################Figure 3 y S4 #################################
Cur_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = " Curated list")
#Generate table for graphing
Table_cur <- Cur_list %>%
  mutate(
    Group = case_when(
      Class %in% c("Actinopteri", "Chondrichthyes") ~ "Fishes",
      Class == "Aves" ~ "Birds",
      Class == "Mammalia" ~ "Mammals",
      Class == "Amphibia" ~ "Amphib",
      is.na(Class) | Class == "NA" ~ "Reptiles",
      TRUE ~ "Other"
    ),
    Resolution = case_when(
      str_detect(`Final taxa`, "^[A-Za-z]{3,}_[A-Za-z]{3,}$") ~ "Specie",
      str_detect(`Final taxa`, "(sp\\.?$|_sp\\.?)") ~ "Genus",
      str_detect(`Final taxa`, "dae$") ~ "Family",
      TRUE ~ NA_character_
    )
  )
#Generate a separate data frame with only those columns and the count.
tabla_resumen <- Table_cur %>%
  filter(!is.na(Resolution)) %>%
  count(Group, Resolution, name = "Taxa") %>%
  complete(Group, Resolution, fill = list(Taxa = 0)) %>%
  arrange(Group, factor(Resolution, levels = c("Specie", "Genus", "Family")))

#Graph of curing and pre-curing of vertebrate groups
#Load verification or basic curation data
Ver_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Simple verification")
#Generate table for graphing
Table_verf <- Ver_list %>%
  mutate(
    Group = case_when(
      Class %in% c("Actinopteri", "Chondrichthyes") ~ "Fishes",
      Class == "Aves" ~ "Birds",
      Class == "Mammalia" ~ "Mammals",
      Class == "Amphibia" ~ "Amphib",
      is.na(Class) | Class == "NA" ~ "Reptiles",
      TRUE ~ "Other"
    ),
    Resolution = case_when(
      str_detect(`Final taxa`, "^[A-Za-z]{3,}_[A-Za-z]{3,}$") ~ "Specie",
      str_detect(`Final taxa`, "(sp\\.?$|_sp\\.?)") ~ "Genus",
      str_detect(`Final taxa`, "dae$") ~ "Family",
      TRUE ~ NA_character_
    )
  )

#Generate a separate data frame with only those columns and the count.
tabla_resumen_verf <- Table_verf %>%
  filter(!is.na(Resolution)) %>%
  count(Group, Resolution, name = "Taxa") %>%
  complete(Group, Resolution, fill = list(Taxa = 0)) %>%
  arrange(Group, factor(Resolution, levels = c("Specie", "Genus", "Family")))
#Load non-curing data
org_list<-readxl::read_xlsx("outputs//Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Magdalena_data_2022_teleo_vert_")

#Generate table for graphing
table_org <- org_list %>%
  mutate(
    Group = case_when(
      Class %in% c("Actinopteri", "Chondrichthyes") ~ "Fishes",
      Class == "Aves" ~ "Birds",
      Class == "Mammalia" ~ "Mammals",
      Class == "Amphibia" ~ "Amphib",
      is.na(Class) | Class == "NA" ~ "Reptiles",
      TRUE ~ "Other"
    ),
    Resolution = case_when(
      str_detect(`Final taxa`, "^[A-Za-z]{3,}_[A-Za-z]{3,}$") ~ "Specie",
      str_detect(`Final taxa`, "(sp\\.?$|_sp\\.?)") ~ "Genus",
      str_detect(`Final taxa`, "dae$") ~ "Family",
      str_detect(`Final taxa`, "NA") ~ "Higher category",
      TRUE ~ NA_character_
    )
  )
table_org<-table_org[-(256:257),]
#Generate a separate data frame with only those columns and the count.
tabla_resumen_org <- table_org %>%
  filter(!is.na(Resolution)) %>%
  count(Group, Resolution, name = "Taxa") %>%
  complete(Group, Resolution, fill = list(Taxa = 0)) %>%
  arrange(Group, factor(Resolution, levels = c("Specie", "Genus", "Family")))

#Create a table with the 3 data frames
tabla_resumen_verf$Catg<-"Verification"
tabla_resumen$Catg<-"Curated"
tabla_resumen_org$Catg <-"Non-curated"
Full_table<-rbind(tabla_resumen_verf,tabla_resumen,tabla_resumen_org)
#set order
Full_table$Catg <- factor(
  Full_table$Catg,
  levels = c("Non-curated", "Verification", "Curated")
)
Full_table$Group <- factor(
  Full_table$Group,
  levels = c("Fishes", "Mammals", "Birds", "Amphib", "Reptiles")
)

Full_table %>%
  filter(Resolution == "Genus") %>%
  summarise(suma_specie = sum(Taxa))
Full_table$Resolution <- factor(
  Full_table$Resolution,
  levels = c("Higher category", "Family", "Genus", "Specie"))

# Calculate totals by group and category
Grupos_totals <- Full_table %>%
  group_by(Group, Catg) %>%
  summarise(total_taxa = sum(Taxa), .groups = "drop")


Plots3a <- ggplot(Full_table, aes(x = Catg, y = Taxa, fill = Resolution)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(Taxa >= 5, Taxa, "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "gray4"
  ) +
  geom_text(
    data = Grupos_totals,
    aes(x = Catg, y = total_taxa, label = total_taxa),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 4,
    color = "black"
  ) +
  facet_grid(. ~ Group, space = "free_x", scales = "free_x", switch = "x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size= 11), 
    legend.position = "right",) +
  scale_fill_manual(values = c(
    "Specie" = "#6CA6CD",
    "Genus" = "#FF6A6A",
    "Family" = "#66CDAA",
    "Higher category" = "#FFB90F"
  )) +
  labs(
    x = "Vertebrate group",
    y = "Number of taxa",
    fill = "Taxonomic resolution"
  )

###general graph
library(dplyr)

# summarize by removing the Group variable
Full_table_summary <- Full_table %>%
  group_by(Catg, Resolution) %>%
  summarise(Taxa = sum(Taxa), .groups = "drop") %>%
  #sum +2 only to Pre-curing en Higher category (taxa without a category in the table )
  mutate(Taxa = ifelse(Catg == "Non-curated" & Resolution == "Higher category",
                       Taxa + 2, Taxa))

#summarize the totals by category (for the numbers above the bars)
Grupos_totals_summary <- Full_table_summary %>%
  group_by(Catg) %>%
  summarise(total_taxa = sum(Taxa), .groups = "drop")

Full_table_summary$Resolution <- factor(
  Full_table_summary$Resolution,
  levels = c("Higher category", "Family", "Genus", "Specie"))
#general graph
Plotf3a <- ggplot(Full_table_summary, aes(x = Catg, y = Taxa, fill = Resolution)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(Taxa >= 5, Taxa, "")),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "gray4"
  ) +
  geom_text(
    data = Grupos_totals_summary,
    aes(x = Catg, y = total_taxa, label = total_taxa),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 5,
    color = "black"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size =16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"), # <-- título en negrita
        legend.position = "top") +
  scale_fill_manual(values = c(
    "Specie" = "#6CA6CD",
    "Genus" = "#FF6A6A",
    "Family" = "#66CDAA",
    "Higher category" = "#FFB90F"
  )) +
  labs(
    x = NULL,
    y = "Number of taxa",
    fill = "Taxonomic resolution"
  )+
  guides(fill = guide_legend(
    title.position = "top",  # título arriba
    title.hjust = 0.5, 
    # <-- asegura que quede centrado
    nrow = 1                # 1 fila para las keys (puedes cambiar a 2 si quieres dividirlas)
  ))

#load faunal list

faunal_list<-readxl::read_xlsx("utils/Faunal list.xlsx", sheet = "Faunal list")
#Generate a list by vertebrate group
Fish1_faunal_list <- faunal_list[grepl("Actinopteri", faunal_list$class, ignore.case = TRUE), ]
Fish2_faunal_list <- faunal_list[grepl("Chondrichthyes", faunal_list$class, ignore.case = TRUE), ]
Fish_faunal_list <- rbind(Fish1_faunal_list, Fish2_faunal_list)
Anfi_faunal_list <- faunal_list[grepl("Amphibia", faunal_list$class, ignore.case = TRUE), ]
Aves_faunal_list <- faunal_list[grepl("Aves", faunal_list$class, ignore.case = TRUE), ]
Mam_faunal_list <- faunal_list[grepl("Mammalia", faunal_list$class, ignore.case = TRUE), ]
Rept_faunal_list <- faunal_list[grepl("Reptilia", faunal_list$class, ignore.case = TRUE), ]

#Load the Teleo reference database
Database_teleo <- readLines("utils/Dada_mitofish_teleo_taxonomy")

# Initialize lists
headers <- c()
sequences <- c()
seq_actual <- ""

for (line in Database_teleo) {
  if (startsWith(line, ">")) {
    if (seq_actual != "") {
      sequences <- c(sequences, seq_actual)
      seq_actual <- ""
    }
    headers <- c(headers, sub("^>", "", line))  # quitar el símbolo ">"
  } else {
    seq_actual <- paste0(seq_actual, line)
  }
}
# Add the latest sequence
sequences <- c(sequences, seq_actual)

# Separate the headers with ';'
header_parts <- strsplit(headers, ";")

# make  data frame
Database_teleo_f <- data.frame(
  kingdom   = sapply(header_parts, `[`, 1),
  phylum    = sapply(header_parts, `[`, 2),
  class     = sapply(header_parts, `[`, 3),
  order     = sapply(header_parts, `[`, 4),
  family    = sapply(header_parts, `[`, 5),
  genus     = sapply(header_parts, `[`, 6),
  species   = sapply(header_parts, `[`, 7),
  sequence  = sequences,
  stringsAsFactors = FALSE
)
#Load the vert01 reference database
Database_vert <- readLines("utils/Dada_mitofish_midori_mitofish_vert_taxonomy")

# Initialize lists
headers <- c()
sequences <- c()
seq_actual <- ""

for (line in Database_vert) {
  if (startsWith(line, ">")) {
    if (seq_actual != "") {
      sequences <- c(sequences, seq_actual)
      seq_actual <- ""
    }
    headers <- c(headers, sub("^>", "", line))  # quitar el símbolo ">"
  } else {
    seq_actual <- paste0(seq_actual, line)
  }
}
# Add the latest sequence
sequences <- c(sequences, seq_actual)

# Separate the headers with ';'
header_parts <- strsplit(headers, ";")

# make data frame
Database_vert_f <- data.frame(
  kingdom   = sapply(header_parts, `[`, 1),
  phylum    = sapply(header_parts, `[`, 2),
  class     = sapply(header_parts, `[`, 3),
  order     = sapply(header_parts, `[`, 4),
  family    = sapply(header_parts, `[`, 5),
  genus     = sapply(header_parts, `[`, 6),
  species   = sapply(header_parts, `[`, 7),
  sequence  = sequences,
  stringsAsFactors = FALSE
)
#join both databases
Databasefull<-rbind(Database_vert_f,Database_teleo_f)
# Function that calculates metrics for species, genus, and family
calcular_metricas_taxonomicas <- function(faunal_df, db_df) {
  #  SPECIES
  species_valid <- unique(faunal_df$species[!is.na(faunal_df$species)])
  species_in_db <- species_valid[species_valid %in% db_df$species]
  species_not_in_db <- species_valid[!species_valid %in% db_df$species]
  
  #  GENUS
  genus_valid <- unique(faunal_df$genus[!is.na(faunal_df$genus)])
  genus_in_db <- genus_valid[genus_valid %in% db_df$genus]
  genus_not_in_db <- genus_valid[!genus_valid %in% db_df$genus]
  
  # FAMILY
  family_valid <- unique(faunal_df$family[!is.na(faunal_df$family)])
  family_in_db <- family_valid[family_valid %in% db_df$family]
  family_not_in_db <- family_valid[!family_valid %in% db_df$family]
  
  # SUMMARY OF METRICS
  resumen <- data.frame(
    row.names = c("No faunal list", "Present in DB", "Not in DB"),
    Species = c(length(species_valid), length(species_in_db), length(species_not_in_db)),
    Genus   = c(length(genus_valid), length(genus_in_db), length(genus_not_in_db)),
    Family  = c(length(family_valid), length(family_in_db), length(family_not_in_db))
  )
  
  # TAXA LIST
  listados <- list(
    Species = list(
      valid = species_valid,
      in_db = species_in_db,
      not_in_db = species_not_in_db
    ),
    Genus = list(
      valid = genus_valid,
      in_db = genus_in_db,
      not_in_db = genus_not_in_db
    ),
    Family = list(
      valid = family_valid,
      in_db = family_in_db,
      not_in_db = family_not_in_db
    )
  )
  
  return(list(
    resumen = resumen,
    listados = listados
  ))
}


#Apply the function to the lists of fauna by vertebrate groups
resultado_fish <- calcular_metricas_taxonomicas(Fish_faunal_list, Databasefull)
resultado_Anfi <- calcular_metricas_taxonomicas(Anfi_faunal_list, Databasefull)
resultado_Mam <- calcular_metricas_taxonomicas(Mam_faunal_list, Databasefull)
resultado_Aves <- calcular_metricas_taxonomicas(Aves_faunal_list, Databasefull)
resultado_Rept <- calcular_metricas_taxonomicas(Rept_faunal_list, Databasefull)

#Extract the results table for each case
resumen_Fish <- resultado_fish$resumen
resumen_Mam <- resultado_Mam$resumen
resumen_Aves <- resultado_Aves$resumen
resumen_Rept <- resultado_Rept$resumen
resumen_Anfi <- resultado_Anfi$resumen

# Add a list-type column to each summary
resumen_Fish$type_list<- "Fishes"
resumen_Mam$type_list<- "Mammals"
resumen_Aves$type_list<- "Birds"
resumen_Rept$type_list<- "Reptiles"
resumen_Anfi$type_list<- "Amphib"
# Join tables
resumen_completo <- bind_rows(resumen_Fish, resumen_Aves,resumen_Mam,resumen_Rept, resumen_Anfi)

# Add the Presence column from the rownames
resumen_completo$Presence <- rownames(resumen_completo)
#Remove periods and numbers from the rows in the “presence” column
resumen_completo$Presence <- gsub("\\.\\.\\..*$", "", resumen_completo$Presence)

# Convert to long format
resumen_largo <- resumen_completo %>%
  pivot_longer(cols = c(Species, Genus, Family),
               names_to = "Taxonomic_level",
               values_to = "Count")

# Ensure the order of the factors 
resumen_largo$Taxonomic_level <- factor(resumen_largo$Taxonomic_level, levels = c("Species", "Genus", "Family"))
resumen_largo$type_list <- factor(resumen_largo$type_list, levels = c("Fishes", "Mammals", "Birds", "Amphib", "Reptiles"))
resumen_largo$Presence <- factor(resumen_largo$Presence, levels = c("Present in DB", "Not in DB"))
#Remove NA entries from the attendance column (unnecessary categories)
resumen_largo <- resumen_largo %>%
  filter(!is.na(Presence))

# Calculate totals by group and category
resumen_largo_total <- resumen_largo %>%
  group_by(type_list, Taxonomic_level) %>%
  summarise(total_taxa = sum(Count), .groups = "drop") %>%
  mutate(y_pos = 1.02)   # un poco arriba del 100%
PlotS3b<-ggplot(resumen_largo, aes(x = Taxonomic_level, y = Count, fill = Presence)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = ifelse(Count>=1, Count, "" )),
            position = position_fill(vjust = 0.5), size = 4, color = "gray4") +
  geom_text(data = resumen_largo_total,
            aes(x = Taxonomic_level, y = y_pos, label = total_taxa),
            inherit.aes = FALSE,
            vjust = 0,
            size = 4,
            color = "black") +
  facet_grid(. ~ type_list, space = "free_x", scales = "free_x", switch = "x") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  theme(legend.position = "right", 
        axis.text.x = element_text(size = 11))+
  labs(x = "Vertebrate group", y = "Percentage of taxa", fill = "Status in database") +
  scale_fill_manual(values = c("Present in DB" = "#458B74", "Not in DB" = "#FF6A6A"))


resumen_total <- resumen_largo %>%
  group_by(Taxonomic_level, Presence) %>%
  summarise(Count = sum(Count), .groups = "drop")

# Ensure the correct order of the factors 
resumen_total$Taxonomic_level <- factor(resumen_total$Taxonomic_level, levels = c("Species", "Genus", "Family"))
resumen_total$Presence <- factor(resumen_total$Presence, levels = c("Present in DB", "Not in DB"))
resumen_total <- resumen_total %>%
  filter(!is.na(Presence))
# Summarize the totals by category (for the numbers above the bars)
# Calculate totals by taxonomic level
Resumen_total_global <- resumen_total %>%
  group_by(Taxonomic_level) %>%
  summarise(total_taxa = sum(Count), .groups = "drop") %>%
  mutate(y_pos = 1.02)   # un poco arriba del 100%

# plot
Plotf3b <- ggplot(resumen_total, aes(x = Taxonomic_level, y = Count, fill = Presence)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = ifelse(Count >= 1, Count, "" )),
            position = position_fill(vjust = 0.5), size = 5, color = "gray4") +
  geom_text(data = Resumen_total_global,
            aes(x = Taxonomic_level, y = y_pos, label = total_taxa),
            inherit.aes = FALSE,
            vjust = 0,
            size = 5,
            color = "black") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +  # agrega espacio arriba
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size =16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"), # <-- título en negrita
        legend.position = "top") +
  labs(x = NULL , y = "Percentage of taxa", fill = "Status in database (DB)") +
  scale_fill_manual(values = c("Present in DB" = "#458B74", "Not in DB" = "#FF6A6A"))+
  guides(fill = guide_legend(
    title.position = "top",  # título arriba
    title.hjust = 0.5, 
    # <-- asegura que quede centrado
    nrow = 1                # 1 fila para las keys (puedes cambiar a 2 si quieres dividirlas)
  ))
Plotf3b
# Generate and save Figure 3
combined_plot <- plot_grid(
  Plotf3a,
  Plotf3b,
  ncol = 2, labels = c("a)", "b)"),
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_size = 20
)

png(filename = "/Figure 3.png",
    width = 3500,height = 1800 ,units = "px", res = 300)
combined_plot
dev.off()
png(filename = "./Figure 3.TIFF",
    width = 3500,height = 1800 ,units = "px", res = 300)
combined_plot
dev.off()

# Generate and save Figure S4
combined_plot2 <- plot_grid(Plots3a, 
                            plot_grid(PlotS3b, ncol = 1, labels = c("b)" ), 
                                      label_fontfamily = "serif", label_fontface = "plain", label_size = 20), 
                            ncol = 1, labels = c("a)", ""),label_fontfamily = "serif", label_fontface = "plain", label_size = 20)

png(filename = "Figure S4.png",
    width = 3500,height = 3000 ,units = "px", res = 300)
combined_plot2
dev.off()

png(filename = "Figure S4.TIFF",
    width = 3000,height = 3000 ,units = "px", res = 300)
combined_plot2
dev.off()

################################# Figure 4 ##########################################
library(dplyr)
library(stringr)
# Load curated list
Cur_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Curated list")
#filter the data frame to include only the fish
Cur_list_filtered <- Cur_list %>%
  filter(Class %in% c("Actinopteri", "Chondrichthyes"))

# Detect column names of the Sxx type
site_cols <- grep("^S\\d{1,2}$", names(Cur_list_filtered), value = TRUE)
# create a list to store results
richness_by_site <- list()

for (site in site_cols) {
  # Filter only the taxa present (≥ 1)
  filtered <- Cur_list_filtered %>%
    filter(.data[[site]] >= 1)
  
  # Filter only for species
  species_filtered <- filtered %>%
    filter(str_detect(`Final taxa`, "^[A-Za-z]{3,}_[A-Za-z]{3,}$"))
  
  # Filter only for genus 
  genus_filtered <- filtered %>%
    filter(!is.na(Genus), Genus != "NA")
  
  # Filter only Families
  family_filtered <- filtered %>%
    filter(!is.na(Family), Family != "NA")
  
  # number of species 
  n_species <- species_filtered %>%
    distinct(`Final taxa`) %>%
    nrow()
  
  # Number of genus
  n_genus <- genus_filtered %>%
    distinct(Genus) %>%
    nrow()
  
  # Number of  families
  n_family <- family_filtered %>%
    distinct(Family) %>%
    nrow()
  
  # Total number of taxa detected (rows with counts >= 1)
  n_taxones <- nrow(filtered)
  
  # Total number of read per group 
  n_reads_total  <- sum(filtered[[site]], na.rm = TRUE)
  n_reads_species <- sum(species_filtered[[site]], na.rm = TRUE)
  n_reads_genus   <- sum(genus_filtered[[site]], na.rm = TRUE)
  n_reads_family  <- sum(family_filtered[[site]], na.rm = TRUE)
  
  # save results
  richness_by_site[[site]] <- tibble(
    Site = site,
    Species = n_species,
    Species_reads = n_reads_species,
    Genus = n_genus,
    Genus_reads = n_reads_genus,
    Family = n_family,
    Family_reads = n_reads_family,
    n_taxones = n_taxones,
    n_reads_total = n_reads_total
  )
}


# combine all into a single data frame 
richness_table <- bind_rows(richness_by_site)
# Load coordinate data by location
library(geosphere)
coord<-readxl::read_xlsx("utils/Abreviaturasf.xlsx", sheet = "Coord")
coord <- coord %>%
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat))
S0_coords <- coord %>%
  filter(Site == "S0") %>%
  dplyr::select(lon, lat) %>%
  unlist() 

# Calculate the distance from S0 to each point
coord <- coord %>%
  rowwise() %>%
  mutate(
    Distance_km = distHaversine(c(lon, lat), S0_coords) / 1000
  ) %>%
  ungroup()
#Seleccionar columas de interes
Distge<-coord[,-(2:3)]
#Join this dataframe with the one created by taxon and location
final_table_cur <- left_join(richness_table, Distge, by = "Site")
#sort by distance
final_table_cur <- final_table_cur %>% arrange(desc(Distance_km))
# Make sure Site is in the correct order
final_table_cur$Site <- factor(final_table_cur$Site, levels = final_table_cur$Site)
#Generating loops with models
library(MASS)
library(broom)
library(dplyr)
library(purrr)
library(tibble)

# Responde variables
responses <- c("Species", "Genus", "Family", "n_taxones")

# Predictor variables
predictors <- list(
  "n_reads",
  "Distance_km",
  "n_reads + Distance_km"
)

#create list
resultados_cur <- list()
modelos_guardados_cur <- list()

#Iteration with automatic outlier exclusion and summary by AIC
for (resp in responses) {
  
  modelos <- list()
  
  for (pred in predictors) {
    
    # If the predictor is n_reads, use the corresponding column based on the taxonomic level
    pred_col <- pred
    
    # Replace `n_reads` with the correct column
    reads_col <- dplyr::case_when(
      resp == "Species" ~ "Species_reads",
      resp == "Genus"   ~ "Genus_reads",
      resp == "Family"  ~ "Family_reads",
      TRUE              ~ "n_reads_total"
    )
    
    if (pred == "n_reads") {
      pred_col <- reads_col
    }
    
    if (pred == "n_reads + Distance_km") {
      pred_col <- paste(reads_col, "+ Distance_km")
    }
    
    # dinamic formula
    form <- as.formula(paste(resp, "~", pred_col))
    
    # Initial model fitting
    fit <- try(glm.nb(form, data = final_table_cur), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      
      # Initial models
      resid_pearson <- residuals(fit, type = "pearson")
      cooks <- cooks.distance(fit)
      hat <- hatvalues(fit)
      
      # Identify outliers
      outliers <- which(abs(resid_pearson) > 2)
      
      # Fit without outliers, if any
      if (length(outliers) > 0) {
        data_clean <- final_table_cur[-outliers, ]
        fit_clean <- try(glm.nb(form, data = data_clean), silent = TRUE)
      } else {
        fit_clean <- fit
        data_clean <- final_table_cur
      }
      
      # Extract coefficients and p-values
      coefs <- broom::tidy(fit_clean) %>%
        dplyr::select(any_of(c("term", "estimate", "std.error", "p.value")))
      
      # Create a list of results from this predictor
      resultado_pred <- list(
        modelo_inicial = fit,
        modelo_limpio = fit_clean,
        outliers_detectados = final_table_cur$Site[outliers],
        AIC_inicial = AIC(fit),
        AIC_limpio = AIC(fit_clean),
        coeficientes = coefs,
        diagnosticos = tibble(
          Site = data_clean$Site,
          resid_pearson = residuals(fit_clean, type = "pearson"),
          cooks = cooks.distance(fit_clean),
          leverage = hatvalues(fit_clean)
        )
      )
      
      # Save to local models list
      modelos[[pred]] <- resultado_pred
      
      # save global list
      if (!resp %in% names(modelos_guardados_cur)) {
        modelos_guardados_cur[[resp]] <- list()
      }
      modelos_guardados_cur[[resp]][[pred]] <- resultado_pred
    }
  }
  
  # Save summary sorted by AIC (clean model)
  resultados_cur[[resp]] <- purrr::map_df(modelos, function(x) {
    tibble(
      Modelo = deparse(formula(x$modelo_limpio)),
      Predictor = paste(as.character(x$modelo_limpio$formula)[3], collapse = " "),
      AIC_inicial = x$AIC_inicial,
      AIC_limpio = x$AIC_limpio,
      outliers = ifelse(length(x$outliers_detectados) > 0,
                        paste(x$outliers_detectados, collapse = ", "),
                        "Ninguno"),
      p_values = paste(round(x$coeficientes$p.value, 4), collapse = "; ")
    )
  }, .id = "ID") %>% arrange(AIC_limpio)
}


# load list 
org_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Magdalena_data_2022_teleo_vert_")
#filter the data frame to include only the fish
org_list_filtered <- org_list %>%
  filter(Class %in% c("Actinopteri", "Chondrichthyes"))
# Detect column names of the Sxx type
site_cols <- grep("^S\\d{1,2}$", names(org_list_filtered), value = TRUE)

# make list to save the results 
richness_by_site <- list()

for (site in site_cols) {
  # Filter only the taxa present (≥ 1)
  filtered <- org_list_filtered %>%
    filter(.data[[site]] >= 1)
  
  # filter only species
  species_filtered <- filtered %>%
    filter(str_detect(`Final taxa`, "^[A-Za-z]{3,}_[A-Za-z]{3,}$"))
  
  # filter only genus
  genus_filtered <- filtered %>%
    filter(!is.na(Genus), Genus != "NA")
  
  #filter only families
  family_filtered <- filtered %>%
    filter(!is.na(Family), Family != "NA")
  
  # Number of species
  n_species <- species_filtered %>%
    distinct(`Final taxa`) %>%
    nrow()
  
  # Number of genus
  n_genus <- genus_filtered %>%
    distinct(Genus) %>%
    nrow()
  
  # Number of families
  n_family <- family_filtered %>%
    distinct(Family) %>%
    nrow()
  
  # Total number of taxa detected  (rows with reads >= 1)
  n_taxones <- nrow(filtered)
  
  # Total number of read by groups
  n_reads_total  <- sum(filtered[[site]], na.rm = TRUE)
  n_reads_species <- sum(species_filtered[[site]], na.rm = TRUE)
  n_reads_genus   <- sum(genus_filtered[[site]], na.rm = TRUE)
  n_reads_family  <- sum(family_filtered[[site]], na.rm = TRUE)
  
  # save results
  richness_by_site[[site]] <- tibble(
    Site = site,
    Species = n_species,
    Species_reads = n_reads_species,
    Genus = n_genus,
    Genus_reads = n_reads_genus,
    Family = n_family,
    Family_reads = n_reads_family,
    n_taxones = n_taxones,
    n_reads_total = n_reads_total
  )
}


# combine all into a single dataframe
richness_table_org <- bind_rows(richness_by_site)
#Join this dataframe with the one created by taxon and location
final_table_org <- left_join(richness_table_org, Distge, by = "Site")
#sort by distance
final_table_org <- final_table_org %>% arrange(desc(Distance_km))
# Make sure Site is in the correct order
final_table_org$Site <- factor(final_table_org$Site, levels = final_table_org$Site)

#generate models with loops
library(MASS)
library(broom)
library(dplyr)
library(purrr)
library(tibble)

# Response variables
responses <- c("Species", "Genus", "Family", "n_taxones")

# predict variables  
predictors <- list(
  "n_reads",
  "Distance_km",
  "n_reads + Distance_km"
)

# list to save results
resultados_org <- list()
modelos_guardados_org <- list()

#Iterative optimization with automatic outlier exclusion and AIC-based summary 
for (resp in responses) {
  
  modelos <- list()
  
  for (pred in predictors) {
    
    # If the predictor is n_reads, use the corresponding column based on the taxonomic level
    pred_col <- pred
    
    # Replace `n_reads` with the correct column
    reads_col <- dplyr::case_when(
      resp == "Species" ~ "Species_reads",
      resp == "Genus"   ~ "Genus_reads",
      resp == "Family"  ~ "Family_reads",
      TRUE              ~ "n_reads_total"
    )
    
    if (pred == "n_reads") {
      pred_col <- reads_col
    }
    if (pred == "n_reads + Distance_km") {
      pred_col <- paste(reads_col, "+ Distance_km")
    }
    
    # Dynamic formula
    form <- as.formula(paste(resp, "~", pred_col))
    
    # Ajuste inicial del modelo
    fit <- try(glm.nb(form, data = final_table_org), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      
      # Initial models
      resid_pearson <- residuals(fit, type = "pearson")
      cooks <- cooks.distance(fit)
      hat <- hatvalues(fit)
      
      # Identificar outliers
      outliers <- which(abs(resid_pearson) > 2)
      
      # Fit without outliers, if any
      if (length(outliers) > 0) {
        data_clean <- final_table_org[-outliers, ]
        fit_clean <- try(glm.nb(form, data = data_clean), silent = TRUE)
      } else {
        fit_clean <- fit
        data_clean <- final_table_org
      }
      
      # Extract coefficients and p-values
      coefs <- broom::tidy(fit_clean) %>%
        dplyr::select(any_of(c("term", "estimate", "std.error", "p.value")))
      
      # Create a list of results for this predictor
      resultado_pred <- list(
        modelo_inicial = fit,
        modelo_limpio = fit_clean,
        outliers_detectados = final_table_org$Site[outliers],
        AIC_inicial = AIC(fit),
        AIC_limpio = AIC(fit_clean),
        coeficientes = coefs,
        diagnosticos = tibble(
          Site = data_clean$Site,
          resid_pearson = residuals(fit_clean, type = "pearson"),
          cooks = cooks.distance(fit_clean),
          leverage = hatvalues(fit_clean)
        )
      )
      
      # Save to local template list
      modelos[[pred]] <- resultado_pred
      
      # save global list
      if (!resp %in% names(modelos_guardados_org)) {
        modelos_guardados_org[[resp]] <- list()
      }
      modelos_guardados_org[[resp]][[pred]] <- resultado_pred
    }
  }
  
  # Save summary sorted by AIC (clean model)
  resultados_org[[resp]] <- purrr::map_df(modelos, function(x) {
    tibble(
      Modelo = deparse(formula(x$modelo_limpio)),
      Predictor = paste(as.character(x$modelo_limpio$formula)[3], collapse = " "),
      AIC_inicial = x$AIC_inicial,
      AIC_limpio = x$AIC_limpio,
      outliers = ifelse(length(x$outliers_detectados) > 0,
                        paste(x$outliers_detectados, collapse = ", "),
                        "Ninguno"),
      p_values = paste(round(x$coeficientes$p.value, 4), collapse = "; ")
    )
  }, .id = "ID") %>% arrange(AIC_limpio)
}

modelo_read_distance<-modelos_guardados_org$n_taxones$Distance_km$modelo_limpio


#make plots 
#relationship chart for curated data
modelo_distance_cur_tax<-modelos_guardados_cur$n_taxones$`n_reads + Distance_km`$modelo_limpio
modelo_distance_cur_esp<-modelos_guardados_cur$Species$`n_reads + Distance_km`$modelo_limpio
modelo_distance_cur_gen<-modelos_guardados_cur$Genus$`n_reads + Distance_km`$modelo_limpio
modelo_distance_cur_fam<-modelos_guardados_cur$Family$`n_reads + Distance_km`$modelo_limpio
#  Create a distance sequence for prediction (setting the value of N_reads to the average)
new_data <- data.frame(
  Distance_km = seq(
    min(final_table_cur$Distance_km),
    max(final_table_cur$Distance_km),
    length.out = 100
  ),
  Species_reads = mean(final_table_cur$Species_reads, na.rm = TRUE),
  Genus_reads = mean(final_table_cur$Genus_reads, na.rm = TRUE),
  Family_reads = mean(final_table_cur$Family_reads, na.rm = TRUE),
  n_reads_total = mean(final_table_cur$n_reads_total, na.rm = TRUE)
)
#generate predictions
pred_species <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_cur_esp, new_data, type = "response"),
  group = "Species"
)
pred_Genus <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_cur_gen, new_data, type = "response"),
  group = "Genus"
)
pred_families <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_cur_fam, new_data, type = "response"),
  group = "Families"
)
pred_taxa <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_cur_tax, new_data, type = "response"),
  group = "Taxa (OTUs)"
)
pred_all <- bind_rows(pred_species, pred_Genus, pred_families, pred_taxa)


library(broom)

#Function to extract the p-value and pseudo-R² for the specified model only
extraer_estadisticas <- function(modelo, termino) {
  
  # Pseudo-R² de McFadden
  r2 <- 1 - as.numeric(logLik(modelo) / logLik(update(modelo, . ~ 1)))
  
  # p-value
  pval <- broom::tidy(modelo) %>%
    dplyr::filter(term == termino) %>%
    dplyr::pull(p.value)
  
  tibble(R2 = r2, p_value = pval)
}

# Extract for each model
stats_taxa <- extraer_estadisticas(modelo_distance_cur_tax, "Distance_km") %>%
  mutate(Taxonomic_level = "Taxa (OTUs)")

stats_esp <- extraer_estadisticas(modelo_distance_cur_esp, "Distance_km") %>%
  mutate(Taxonomic_level = "Species")

stats_gen <- extraer_estadisticas(modelo_distance_cur_gen, "Distance_km") %>%
  mutate(Taxonomic_level = "Genus")

stats_fam <- extraer_estadisticas(modelo_distance_cur_fam, "Distance_km") %>%
  mutate(Taxonomic_level = "Families")
# Combine into a table
stats_all <- bind_rows(stats_taxa, stats_esp, stats_gen, stats_fam) %>%
  mutate(
    R2 = round(R2, 3),
    p_value = signif(p_value, 3)
  )

#Create a table with the stats
tabla_stats <- stats_all %>%
  dplyr::mutate(
    p_value = sprintf("%.4f", p_value),
    R2 = sprintf("%.3f", R2)
  )

# Generate the table as a GROB

tabla_grob <- tableGrob(
  tabla_stats,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.9),           
      bg_params = list(fill = NA, col = "black") 
    ),
    colhead = list(
      fg_params = list(cex = 0.9, fontface = "bold"),
      bg_params = list(fill = NA, col = "black") 
    )
  )
)
# Create a graph using the included table 
Plot_curado <- ggplot() +
  geom_line(data = pred_all,
            aes(x = Distance_km, y = pred, color = group),
            size = 1.2) +
  theme_light() +
  labs(x = "Distance to the mouth (km)",
       y = "Richness",
       color = "Taxonomic level:") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5, hjust=1, size = 16),# Rotar etiquetas del eje x
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(15, 70),
                     breaks = seq(15, 70, 10)) +
  annotation_custom(
    grob = tabla_grob,
    xmin = 205,  
    xmax = 330,
    ymin = 65,
    ymax = 64
  )
Plot_curado

#Relationship chart for ORIGINAL data
modelo_distance_org_tax<-modelos_guardados_org$n_taxones$`n_reads + Distance_km`$modelo_limpio
modelo_distance_org_esp<-modelos_guardados_org$Species$`n_reads + Distance_km`$modelo_limpio
modelo_distance_org_gen<-modelos_guardados_org$Genus$`n_reads + Distance_km`$modelo_limpio
modelo_distance_org_fam<-modelos_guardados_org$Family$`n_reads + Distance_km`$modelo_limpio
# Create a distance sequence for prediction
new_data <- data.frame(
  Distance_km = seq(
    min(final_table_org$Distance_km),
    max(final_table_org$Distance_km),
    length.out = 100
  ),
  Species_reads = mean(final_table_org$Species_reads, na.rm = TRUE),
  Genus_reads = mean(final_table_org$Genus_reads, na.rm = TRUE),
  Family_reads = mean(final_table_org$Family_reads, na.rm = TRUE),
  n_reads_total = mean(final_table_org$n_reads_total, na.rm = TRUE)
)
#generate predictions
pred_species <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_org_esp, new_data, type = "response"),
  group = "Species"
)
pred_Genus <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_org_gen, new_data, type = "response"),
  group = "Genus"
)
pred_families <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_org_fam, new_data, type = "response"),
  group = "Families"
)
pred_taxa <- data.frame(
  Distance_km = new_data$Distance_km,
  pred = predict(modelo_distance_org_tax, new_data, type = "response"),
  group = "Taxa (OTUs)"
)
pred_all <- bind_rows(pred_species, pred_Genus, pred_families, pred_taxa)


library(broom)

#Function to extract the p-value and pseudo-R² for the specified model only
extraer_estadisticas <- function(modelo, termino) {
  
  # Pseudo-R² de McFadden
  r2 <- 1 - as.numeric(logLik(modelo) / logLik(update(modelo, . ~ 1)))
  
  # p-value 
  pval <- broom::tidy(modelo) %>%
    dplyr::filter(term == termino) %>%
    dplyr::pull(p.value)
  
  tibble(R2 = r2, p_value = pval)
}

# Extract for each model
stats_taxa <- extraer_estadisticas(modelo_distance_org_tax, "Distance_km") %>%
  mutate(Taxonomic_level = "Taxa (OTUs)")

stats_esp <- extraer_estadisticas(modelo_distance_org_esp, "Distance_km") %>%
  mutate(Taxonomic_level = "Species")

stats_gen <- extraer_estadisticas(modelo_distance_org_gen, "Distance_km") %>%
  mutate(Taxonomic_level = "Genus")

stats_fam <- extraer_estadisticas(modelo_distance_org_fam, "Distance_km") %>%
  mutate(Taxonomic_level = "Families")
# combine into a table
stats_all <- bind_rows(stats_taxa, stats_esp, stats_gen, stats_fam) %>%
  mutate(
    R2 = round(R2, 3),
    p_value = signif(p_value, 3)
  )


#  Create a table with the stats
tabla_stats <- stats_all %>%
  dplyr::mutate(
    p_value = sprintf("%.4f", p_value),
    R2 = sprintf("%.3f", R2)
  )

# Generate the table as a GROB

tabla_grob <- tableGrob(
  tabla_stats,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.9),          
      bg_params = list(fill = NA, col = "black") 
    ),
    colhead = list(
      fg_params = list(cex = 0.9, fontface = "bold"),
      bg_params = list(fill = NA, col = "black") 
    )
  )
)
#  Create a graph using the included table 
plot_org <- ggplot() +
  geom_line(data = pred_all,
            aes(x = Distance_km, y = pred, color = group),
            size = 1.2) +
  theme_light() +
  labs(x = "Distance to the mouth (km)",
       y = "Richness",
       color = "Taxonomic level") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5, hjust=1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(15, 70),
                     breaks = seq(15, 70, 10)) +
  annotation_custom(
    grob = tabla_grob,
    xmin = 205,  
    xmax = 330,
    ymin = 65,
    ymax = 64
  )

plot_org

# Generate and save Figure 4
combined_plot <- plot_grid(
  plot_org,
  Plot_curado,
  ncol = 2, labels = c("a)", "b)"),
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_size = 20
)

png(filename = "Figure 4.png",
    width = 3700,height = 2000 ,units = "px", res = 300)
combined_plot
dev.off()
png(filename = "Figure 4.TIFF",
    width = 3500,height = 1800 ,units = "px", res = 300)
combined_plot
dev.off()

############Spearman's correlation tests#######################

cor.test(final_table_cur$Distance_km, final_table_cur$n_reads_total, method = "spearman")

cor.test(final_table_org$Distance_km, final_table_org$n_reads_total, method = "spearman")

########################  Figure 5 ##################################
library(fields)
library(dplyr)
library(betapart)
library(ape)


Cur_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Curated list")
Cur_list_filtered <- Cur_list %>%
  filter(Class %in% c("Actinopteri", "Chondrichthyes"))
#Generate a vector containing the column names from the sites
sites <- paste0("S", 1:14)
#Generate a data frame filtered by categories
Filt_list <- Cur_list_filtered %>%
  dplyr::select(`Final taxa`, all_of(sites))
#Remove the column of names
mat_abund <- Filt_list[, -1]
mat_abund <- as.data.frame(mat_abund)
#Remove only S13 as an outlier
mat_abund<-mat_abund[,-(13)]
rownames(mat_abund) <- make.unique(Filt_list$`Final taxa`)
#Create a new data frame
pa_df <- mat_abund
#convert to binary data
pa_df[pa_df > 0] <- 1
#Calculate the dissimilarity matrix using Sorensen 
bet<-beta.pair(t(pa_df), index.family = "sorensen")
bet_cur<-beta.pair(t(pa_df), index.family = "sorensen")

sim_dist<-bet$beta.sim    
sne_dist<-bet$beta.sne    
sor_dist<-bet$beta.sor    
sor_dist_cur<-bet$beta.sor    

library(ade4)
is.euclid(bet$beta.sim)  
is.euclid(bet$beta.sne)
is.euclid(bet$beta.sor)
pcdas<-beta.multi(t(pa_df), index.family = "sorensen")
#calculate the deviation of the metrics
set.seed(123) # for reproducibility

n_iter <- 1000  
n_sites <- ncol(pa_df)

results <- replicate(n_iter, {
  # sample sites with replacement
  sampled_sites <- sample(1:n_sites, n_sites, replace = TRUE)
  sampled_matrix <- pa_df[, sampled_sites]
  
  beta_vals <- beta.multi(t(sampled_matrix), index.family = "sorensen")
  
  c(beta.SIM = beta_vals$beta.SIM,
    beta.SNE = beta_vals$beta.SNE,
    beta.SOR = beta_vals$beta.SOR)
})

# Convert to a data frame
results_df <- as.data.frame(t(results))

# Calculate means and standard deviations
summary_stats <- data.frame(
  Component = c("beta.SIM", "beta.SNE", "beta.SOR"),
  Mean = colMeans(results_df),
  SD = apply(results_df, 2, sd)
)

# Extract coordinates
pcoa_result_sor <- pcoa(sor_dist, correction = "cailliez")
pcoa_scores_sor <- as.data.frame(pcoa_result_sor$vectors[, 1:2])
pcoa_scores_sor_cur<-pcoa_scores_sor
pcoa_scores_sor$Site <- rownames(pcoa_scores_sor)
# Explained variance
Eigen <- pcoa_result_sor$values$Eigenvalues
relative_eig <- ((Eigen / sum(Eigen))*100)
pretty_Eigen<-round(relative_eig[1:2], 2)
library(glue)
labss<- c(glue("PCoA Axis 1 ({pretty_Eigen[1]}%)"),
          glue("PCoA Axis 2 ({pretty_Eigen[2]}%)"))
#Add geographic distance data to the mouth of the river
library(geosphere)
library(dplyr)

coord<-readxl::read_xlsx("../Abreviaturasf.xlsx", sheet = "Coord")
coord <- coord %>%
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat))
S0_coords <- coord %>%
  filter(Site == "S0") %>%
  dplyr::select(lon, lat) %>%
  unlist() 

# Calculate the distance from S0 to each point
coord <- coord %>%
  rowwise() %>%
  mutate(
    Distance_km = distHaversine(c(lon, lat), S0_coords) / 1000
  ) %>%
  ungroup()
#Select the columns of interest
Distge<-coord[,-(2:3)]

# Join this dataframe with the one created by taxon and location
distgeo_sor <- left_join(pcoa_scores_sor, Distge, by = "Site")

# PERMANOVA (adonis2)
permanova <- adonis2(sor_dist ~ Distance_km, data = distgeo_sor)
print(permanova)
F_val <- round(permanova$F[1], 2)
p_val <- permanova$`Pr(>F)`[1]
r_val <- round(permanova$R2 [1],2)
p_val_text <- ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 3)))
annot_text <- paste0("PERMANOVA: F = ", F_val, ", R2 = ",r_val,  ", p ", p_val_text)
# Organize a data frame of geographic distances 
orden <- c("Site", "Axis.1", "Axis.2", "Distance_km")
distgeo_sor <- distgeo_sor[, orden]
# Sort the data frame by distance
distgeo_sor <- distgeo_sor %>% arrange(Distance_km)
# Add a 10% margin to the axis range
margin_x <- 0.1 * diff(range(distgeo_sor$Axis.1))
margin_y <- 0.1 * diff(range(distgeo_sor$Axis.2))
# Create a grid on the ordering space (PCoA)
x_range <- seq(min(distgeo_sor$Axis.1) - margin_x, max(distgeo_sor$Axis.1) + margin_x, length.out = 200)
y_range <- seq(min(distgeo_sor$Axis.2) - margin_y, max(distgeo_sor$Axis.2) + margin_y, length.out = 200)

grid_vals <- expand.grid(Axis.1 = x_range, Axis.2 = y_range)

library(ggplot2)
# Calculate cell size
dx <- diff(unique(grid_vals$Axis.1))[1]
dy <- diff(unique(grid_vals$Axis.2))[1]

# Expand the boundaries of the graph by half a step on each side
xlim_exp <- range(grid_vals$Axis.1) + c(-dx/2, dx/2)
ylim_exp <- range(grid_vals$Axis.2) + c(-dy/2, dy/2)

# Create a column that indicates whether the point is in the group
puntos_grupo <- c("S5", "S6", "S7", "S8", "S9", "S10","S11","S12")
distgeo_sor$in_grupo <- ifelse(distgeo_sor$Site %in% puntos_grupo, "Grupo1", NA)

# Filter and compute the convex hull
hull_manual <- distgeo_sor %>%
  filter(in_grupo == "Grupo1") %>%
  slice(chull(Axis.1, Axis.2))

# Extracting the numeric part from the last column
distgeo_sor$Order <- as.numeric(gsub("\\D", "", distgeo_sor$Site))

# Sorting the data frame based on the "Order" column
distgeo_sor <- distgeo_sor[order(distgeo_sor$Order), ]
distgeo_sor$Order <- factor(distgeo_sor$Order, levels = sort(unique(distgeo_sor$Order)))
my_colors<- c("#FF0000", "#FF4700", "#FF8F00", "#FFD700", "#AAE42A", "#55F154", "#00FF7F", "#15CD9F", "#2B9BC0", "#4169E1",
              "#7A71E5", "#B479E9", "#EE82EE")

#plot
pcoa_sor_curing <- ggplot() +
  
  geom_polygon(data = hull_manual, aes(x = Axis.1, y = Axis.2),
               fill = NA, color = "gray30", linewidth = 0.7) +
  
  geom_point(data = distgeo_sor, aes(x = Axis.1, y = Axis.2, color = Order), size = 4) +
  
  geom_text(data = distgeo_sor, aes(x = Axis.1, y = Axis.2, label = Site), vjust = -1) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = xlim_exp, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_exp, expand = c(0, 0)) +
  
  scale_color_manual(name = "Site Order", values = my_colors) +
  
  coord_fixed() +
  theme_light() +
  
  labs(x = labss[1], y = labss[2], title = "") +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )#+
annotate("text", 
         x = max(xlim_exp) - 0.05 * diff(xlim_exp),
         y = min(ylim_exp) + 0.05 * diff(ylim_exp),
         label = annot_text, 
         hjust = 1, vjust = 0, size = 5)
pcoa_sor_curing
png(filename = "Figure 5b.png",
    width = 2500,height = 1800 ,units = "px", res = 300)
pcoa_sor_curing
dev.off()

#non curated data
library(fields)
library(dplyr)
library(betapart)
Org_list<-readxl::read_xlsx("outputs/Magdalena_data_2022_teleo_vert_dada2.xlsx", sheet = "Magdalena_data_2022_teleo_vert_")
Org_list_filtered <- Org_list %>%
  filter(Class %in% c("Actinopteri", "Chondrichthyes"))
#Generate a vector containing the names of the columns in the sites
sites <- paste0("S", 1:14)
#Generate a data frame filtered by categories
Filt_list <- Org_list_filtered %>%
  dplyr::select(`Final taxa`, all_of(sites))
#Remove the column of names
mat_abund <- Filt_list[, -1]
mat_abund <- as.data.frame(mat_abund)
mat_abund<-mat_abund[,-(13)]
rownames(mat_abund) <- make.unique(Filt_list$`Final taxa`)
#Create a new data frame
pa_df <- mat_abund
#convert to binary data
pa_df[pa_df > 0] <- 1
#Calculate the dissimilarity matrix using Sorensen 
bet<-beta.pair(t(pa_df), index.family = "sorensen")
bet_cur<-beta.pair(t(pa_df), index.family = "sorensen")

sim_dist<-bet$beta.sim    
sne_dist<-bet$beta.sne    
sor_dist<-bet$beta.sor    
sor_dist_cur<-bet$beta.sor    

library(ade4)
is.euclid(bet$beta.sim)  
is.euclid(bet$beta.sne)
is.euclid(bet$beta.sor)
pcdas<-beta.multi(t(pa_df), index.family = "sorensen")
#calculate the deviation of the metrics
set.seed(123) # for reproducibility

n_iter <- 1000  
n_sites <- ncol(pa_df)

results <- replicate(n_iter, {
  # sample sites with replacement
  sampled_sites <- sample(1:n_sites, n_sites, replace = TRUE)
  sampled_matrix <- pa_df[, sampled_sites]
  
  beta_vals <- beta.multi(t(sampled_matrix), index.family = "sorensen")
  
  c(beta.SIM = beta_vals$beta.SIM,
    beta.SNE = beta_vals$beta.SNE,
    beta.SOR = beta_vals$beta.SOR)
})

# Convert to a data frame
results_df <- as.data.frame(t(results))

# Calculate means and standard deviations
summary_stats <- data.frame(
  Component = c("beta.SIM", "beta.SNE", "beta.SOR"),
  Mean = colMeans(results_df),
  SD = apply(results_df, 2, sd)
)

# Extract coordinates
pcoa_result_sor <- pcoa(sor_dist, correction = "cailliez")
pcoa_scores_sor <- as.data.frame(pcoa_result_sor$vectors[, 1:2])
pcoa_scores_sor_cur<-pcoa_scores_sor
pcoa_scores_sor$Site <- rownames(pcoa_scores_sor)
# Explained variance
Eigen <- pcoa_result_sor$values$Eigenvalues
relative_eig <- ((Eigen / sum(Eigen))*100)
pretty_Eigen<-round(relative_eig[1:2], 2)
library(glue)
labss<- c(glue("PCoA Axis 1 ({pretty_Eigen[1]}%)"),
          glue("PCoA Axis 2 ({pretty_Eigen[2]}%)"))
#Add geographic distance data to the mouth of the river
library(geosphere)
library(dplyr)

coord<-readxl::read_xlsx("../Abreviaturasf.xlsx", sheet = "Coord")
coord <- coord %>%
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat))
S0_coords <- coord %>%
  filter(Site == "S0") %>%
  dplyr::select(lon, lat) %>%
  unlist() 

# Calculate the distance from S0 to each point
coord <- coord %>%
  rowwise() %>%
  mutate(
    Distance_km = distHaversine(c(lon, lat), S0_coords) / 1000
  ) %>%
  ungroup()
#Select the columns of interest
Distge<-coord[,-(2:3)]

# Join this dataframe with the one created by taxon and location
distgeo_sor <- left_join(pcoa_scores_sor, Distge, by = "Site")

# PERMANOVA (adonis2)
permanova <- adonis2(sor_dist ~ Distance_km, data = distgeo_sor)
print(permanova)
F_val <- round(permanova$F[1], 2)
p_val <- permanova$`Pr(>F)`[1]
r_val <- round(permanova$R2 [1],2)
p_val_text <- ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 3)))
annot_text <- paste0("PERMANOVA: F = ", F_val, ", R2 = ",r_val,  ", p ", p_val_text)
# Organize a data frame of geographic distances 
orden <- c("Site", "Axis.1", "Axis.2", "Distance_km")
distgeo_sor <- distgeo_sor[, orden]
# Sort the data frame by distance
distgeo_sor <- distgeo_sor %>% arrange(Distance_km)
# Add a 10% margin to the axis range
margin_x <- 0.1 * diff(range(distgeo_sor$Axis.1))
margin_y <- 0.1 * diff(range(distgeo_sor$Axis.2))
# Create a grid on the ordering space (PCoA)
x_range <- seq(min(distgeo_sor$Axis.1) - margin_x, max(distgeo_sor$Axis.1) + margin_x, length.out = 200)
y_range <- seq(min(distgeo_sor$Axis.2) - margin_y, max(distgeo_sor$Axis.2) + margin_y, length.out = 200)

grid_vals <- expand.grid(Axis.1 = x_range, Axis.2 = y_range)

library(ggplot2)
# Calculate cell size
dx <- diff(unique(grid_vals$Axis.1))[1]
dy <- diff(unique(grid_vals$Axis.2))[1]

# Expand the boundaries of the graph by half a step on each side
xlim_exp <- range(grid_vals$Axis.1) + c(-dx/2, dx/2)
ylim_exp <- range(grid_vals$Axis.2) + c(-dy/2, dy/2)

# Create a column that indicates whether the point is in the group
puntos_grupo <- c("S5", "S6", "S7", "S8", "S9", "S10","S11","S12")
distgeo_sor$in_grupo <- ifelse(distgeo_sor$Site %in% puntos_grupo, "Grupo1", NA)

# Filter and compute the convex hull
hull_manual <- distgeo_sor %>%
  filter(in_grupo == "Grupo1") %>%
  slice(chull(Axis.1, Axis.2))

# Extracting the numeric part from the last column
distgeo_sor$Order <- as.numeric(gsub("\\D", "", distgeo_sor$Site))

# Sorting the data frame based on the "Order" column
distgeo_sor <- distgeo_sor[order(distgeo_sor$Order), ]
distgeo_sor$Order <- factor(distgeo_sor$Order, levels = sort(unique(distgeo_sor$Order)))
my_colors<- c("#FF0000", "#FF4700", "#FF8F00", "#FFD700", "#AAE42A", "#55F154", "#00FF7F", "#15CD9F", "#2B9BC0", "#4169E1",
              "#7A71E5", "#B479E9", "#EE82EE")
#graficar
pcoa_sor_orig <- ggplot() +
  
  geom_polygon(data = hull_manual, aes(x = Axis.1, y = Axis.2),
               fill = NA, color = "gray30", linewidth = 0.7) +
  
  geom_point(data = distgeo_sor, aes(x = Axis.1, y = Axis.2, color = Order), size = 4) +
  
  geom_text(data = distgeo_sor, aes(x = Axis.1, y = Axis.2, label = Site), vjust = -1) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = xlim_exp, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_exp, expand = c(0, 0)) +
  
  
  # 🔑 Aquí usamos tu paleta personalizada
  scale_color_manual(name = "Site Order", values = my_colors) +
  
  coord_fixed() +
  theme_light() +
  
  labs(x = labss[1], y = labss[2], title = "") +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )#+

annotate("text", 
         x = max(xlim_exp) - 0.05 * diff(xlim_exp),
         y = min(ylim_exp) + 0.05 * diff(ylim_exp),
         label = annot_text, 
         hjust = 1, vjust = 0, size = 5)
pcoa_sor_orig
png(filename = "../../Proyecto Magdalena/Plots/Figuras finales/Figuras finales v2/Figure 5a.png",
    width = 2500,height = 1800 ,units = "px", res = 300)
pcoa_sor_orig
dev.off()
#############################Perform a Mantel- Procrustes tests  #############################
library(vegan)
mantel_result <- mantel(sor_dist_cur, sor_dist_org, method = "spearman", permutations = 999)

# Procrustes test
proc <- protest(pcoa_scores_sor_org, pcoa_scores_sor_cur, permutations = 999)

#save plots
# Generate and save Figure 3
combined_plot <- plot_grid(
  pcoa_sor_orig,
  pcoa_sor_curing,
  ncol = 1,               
  labels = c("a)", "b)"),
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_size = 20,
  rel_heights = c(1, 1)   # Controla altura relativa (puedes poner 2,1 si quieres que el primero sea más grande)
)

############################# Figure S5 and S6 ##################################
library(reshape2)
bet_org 
turn_mat <- as.matrix(bet_org$beta.sim)
nest_mat <- as.matrix(bet_org$beta.sne)
sor_mat <- as.matrix(bet_org$beta.sor)
sor_melted <- reshape2::melt(sor_mat)
turn_melted <- reshape2::melt(turn_mat)
nest_melted <- reshape2::melt(nest_mat)
#Set the correct order of sites (without S13)
orden_sitios <- paste0("S", c(1:12, 14))
# We create a single data frame
presi <- data.frame(
  Sitio       = rownames(sor_mat)[sor_melted$Var1],
  Comparación = rownames(sor_mat)[sor_melted$Var2],
  Sorensen     = sor_melted$value,
  Recambio    = turn_melted$value,
  Anidamiento = nest_melted$value,
  Disimilitud = rowSums(sor_mat)
)
# Convert to long format and filter out self-comparisons
derretido <- presi %>%
  filter(Sitio != Comparación) %>%
  melt() %>%
  filter(variable %in% c("Recambio", "Anidamiento")) %>%  
  mutate(
    Sitio = factor(Sitio, levels = orden_sitios),
    Comparación = factor(Comparación, levels = orden_sitios)
  )

# Plot
plots <- lapply(split(derretido, derretido$Sitio), function(df) {
  
  # Determine whether this site should display a Y-axis label
  sitio <- unique(df$Sitio)
  mostrar_y <- sitio %in% c("S1", "S5", "S9","S14")
  
  ggplot(df, aes(x = Comparación, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("#66CDAA", "#ff6a6a"),
                      name = "Components",
                      labels = c("Turnover", "Nesting")) +
    scale_y_continuous(
      name = if (mostrar_y) "Dissimilarity" else NULL,  
      limits = c(0, 0.6)
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste("", sitio),
      x = "", y = NULL
    )
})


library(patchwork)
componen_org<-wrap_plots(plots, ncol = 4) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")  

# Curated data  
bet_cur
turn_mat <- as.matrix(bet_cur$beta.sim)
nest_mat <- as.matrix(bet_cur$beta.sne)
sor_mat <- as.matrix(bet_cur$beta.sor)
sor_melted <- reshape2::melt(sor_mat)
turn_melted <- reshape2::melt(turn_mat)
nest_melted <- reshape2::melt(nest_mat)

#Set the correct order of sites (without S13)
orden_sitios <- paste0("S", c(1:12, 14))
# Create a single data frame
presi <- data.frame(
  Sitio       = rownames(sor_mat)[sor_melted$Var1],
  Comparación = rownames(sor_mat)[sor_melted$Var2],
  Sorensen     = sor_melted$value,
  Recambio    = turn_melted$value,
  Anidamiento = nest_melted$value,
  Disimilitud = rowSums(sor_mat)
)
#Convert to long format and filter out self-comparisons
derretido <- presi %>%
  filter(Sitio != Comparación) %>%
  melt() %>%
  filter(variable %in% c("Recambio", "Anidamiento")) %>%   # 🔑 Solo quedarán estas dos columnas
  mutate(
    Sitio = factor(Sitio, levels = orden_sitios),
    Comparación = factor(Comparación, levels = orden_sitios)
  )


# plot

plots <- lapply(split(derretido, derretido$Sitio), function(df) {
  
  # Determine whether this site should display a Y-axis label
  sitio <- unique(df$Sitio)
  mostrar_y <- sitio %in% c("S1", "S5", "S9")
  
  ggplot(df, aes(x = Comparación, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("#66CDAA", "#ff6a6a"),
                      name = "Components",
                      labels = c("Turnover", "Nesting")) +
    scale_y_continuous(
      name = if (mostrar_y) "Dissimilarity" else NULL,  # 👈 Condicional
      limits = c(0, 0.6)
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste("", sitio),
      x = "", y = NULL
    )
})


library(patchwork)
compone_cur<-wrap_plots(plots, ncol = 4) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")  # 🔑 Leyenda única en la parte superior

png(filename = "Figure S5.png",
    width = 3500,height = 2200 ,units = "px", res = 300)
componen_org
dev.off()
png(filename = "Figure S6.png",
    width = 3500,height = 2200 ,units = "px", res = 300)
compone_cur
dev.off()
