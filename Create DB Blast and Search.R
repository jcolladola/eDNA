## Vamos a crear una base de datos con la lista de especies que queremos.

install.packages("rentrez")
install.packages("dplyr")
library(rentrez)
library(dplyr)

# Importamos la lista de especies que queremos tener en la base de datos.

fish_species <- c("Achirus mazatlanus", "Anableps dowei", "Oxyzygonectes dovii", "Anguilla rostrata", "Cathorops fuerthii", "Cathorops steindachneri", 
                  "Sciades assimilis", "Sciades seemanni", "Atherinella argentea", "Atherinella hubbsi", "Atherinella jiloaensis", "Atherinella sardina", 
                  "Strongylura exilis", "Strongylura marina", "Oligoplites palometa", "Carcharhinus leucas", "Centropomus ensiferus", "Centropomus nigrescens", 
                  "Centropomus parallelus", "Centropomus pectinatus", "Astyanax aeneus", "Astyanax nasutus", "Astyanax nicaraguensis", "Bramocharax bransfordii", 
                  "Brycon guatemalensis", "Carlana eigenmanni", "Roeboides bouchellei", "Amatitlania nigrofasciata", "Amatitlania siquia", "Amphilophus alfari", 
                  "Amphilophus amarillo", "Amphilophus astorquii", "Amphilophus chancho", "Amphilophus citrinellus", "Amphilophus flaveolus", "Amphilophus labiatus", 
                  "Amphilophus longimanus", "Amphilophus rostratus", "Amphilophus sagittae", "Amphilophus xiloaensis", "Amphilophus zaliosus", "Archocentrus centrarchus", 
                  "Archocentrus multispinosus", "Cichlasoma urophthalmus", "Cryptoheros septemfasciatus", "Cryptoheros spilurus", "Hypsophrys nematopus", 
                  "Hypsophrys nicaraguensis", "Oreochromis aureus", "Oreochromis mossambicus", "Oreochromis niloticus niloticus", "Oreochromis urolepis hornorum", 
                  "Parachromis dovii", "Parachromis friedrichsthalii", "Parachromis loisellei", "Parachromis managuensis", "Tomocichla tuba", "Vieja maculicauda", 
                  "Dorosoma chavesi", "Lile stolifera", "Acromycter atlanticus", "Cyprinus carpio carpio", "Dormitator latifrons", "Erotelis smaragdus", 
                  "Gobiomorus maculatus", "Hemieleotris latifasciata", "Anchoa curta", "Anchoa parva", "Eucinostomus argenteus", "Eucinostomus currani", 
                  "Eucinostomus gracilis", "Eugerres plumieri", "Gerres cinereus", "Awaous banana", "Ctenogobius claytonii", "Gobioides peruanus", 
                  "Sicydium salvini", "Pomadasys bayanus", "Pomadasys crocro", "Rhamdia nicaraguensis", "Rhamdia quelen", "Kuhlia mugil", "Atractosteus tropicus", 
                  "Megalops atlanticus", "Agonostomus monticola", "Joturus pichardi", "Mugil cephalus", "Mugil curema", "Citharichthys gilberti", "Alfaro cultratus", 
                  "Alfaro huberi", "Belonesox belizanus", "Brachyrhaphis holdridgei", "Gambusia nicaraguensis", "Phallichthys amates", "Phallichthys tico", 
                  "Poecilia mexicana", "Priapichthys panamensis", "Xenophallus umbratilis", "Rivulus isthmensis", "Oncorhynchus tshawytscha", "Pseudophallus starksii")

# Obtenemos las secuencias de nucleótidos desde NCBI.

sequences <- lapply(fish_species, function(species) {
  search <- entrez_search(db="nuccore", term=species, retmax=10)
  if (length(search$ids) > 0) {
    seqs <- entrez_fetch(db="nuccore", id=search$ids, rettype="fasta")
    return(seqs)
  } else {
    return(NULL)
  }
})
names(sequences) <- fish_species

# Guardamos las secuencias en un archivo fasta.

write_fasta <- function(seqs, filename) {
  cat(unlist(seqs), file=filename, sep="\n")
}

write_fasta(sequences, "fish_sequences.fasta")

## Posteriormente habrá que crear la base de datos en formato blast.

# Primero subimos el archivo obtenido en fasta al servidor.

# makeblastdb -in fish_sequences.fasta -dbtype nucl -out fish_db

## Realizamos la búsqueda en BLAST utilizando la base de datos.

# blastn -query your_query.fasta -db fish_db -out results.out


