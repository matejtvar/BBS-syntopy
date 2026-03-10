###############################################################################*
##### Sympatry and syntopy across scales in birds                         #####*
##### V. Remes, puvodne poster na GRC 2018, viz skript "GRC2018poster.R"  #####*
##### dalsi prace pro CSPE konferenci 2019 v Olomouci                     #####*
##### dalsi prace pro rukopis, zacatek srpen 2021 - zde selekce kodu      #####*
##### ze starsich souboru (hl. "cooccurrence_Au.R")                       #####*
###############################################################################*

setwd("/Users/Lada/Documents/clanky_a_projekty/GACR_sympatry/PROJEKT_cooccurrenceAu")
setwd("/Users/vladimirremes/Documents/clanky_a_projekty/GACR_sympatry/PROJEKT_cooccurrenceAu")
require(tidyverse)
require(readxl)
require(spaa)
require(cooccur)
require(diverge)
require(visreg)
require(vegan)
require(effects)

##### SPECIES NAMES: translation table -----
names_Marki <- read_csv("results/names_Marki.csv")
names_table <- read_excel("/Users/vladimirremes/Documents/clanky_a_projekty/GACR_Australie/jmena druhu pevcu/Australia_Passerines_names_FINAL.xlsx")
names_table <- names_table %>%
  select(-c(8,11)) %>%
  separate(col = Pizzey_9e, into = c("Pizzey_genus", "Pizzey_species")) %>%
  unite(col = Pizzey_name, Pizzey_genus, Pizzey_species, sep = "_", remove = TRUE) %>%
  rename(HANZAB_name = HANZAB, BLI_gis2011_name = BLI_gis2011) %>%
  separate(col = BLI_ver6, into = c("genus", "species")) %>%
  unite(col = BLI6_name, genus, species, sep = "_", remove = TRUE)
  # pak jeste (pokud netreba ostatni jmena): select(HANZAB_Family:BLI6_name)

setdiff(names_table$BLI_gis2011_name, names_table$BLI6_name)  # 0 rozdilu
setdiff(names_table$BLI6_name, names_table$BLI_gis2011_name)  # 0 rozdilu
# totez jako (jen ukazka - muselo by se take obema smery):
anti_join(x = names_table, y = names_table, by = c("BLI_gis2011_name" = "BLI6_name"))  # 0 rozdilu

names_table_Mel <- names_table %>%
  filter(HANZAB_Family == "Meliphagidae" | HANZAB_Family == "Acanthizidae" | HANZAB_Family == "Maluridae" | HANZAB_Family == "Pardalotidae" | HANZAB_Family == "Dasyornithidae")


##### PHYLOGENY -----
require(ape)
require(geiger)
require(phangorn)
# +++ loading and plotting -----
# full tree:
tree_Marki_full <- read.nexus(file = "data/fylogeneze/marki_2017.tre")
is.ultrametric(tree_Marki_full)  # kontrola - FALSE
which(tree_Marki_full$edge.length<0)  # zadny edge uz neni zaporny
# exclude outgroups = node no. 305:
# tree_Mel_all <- extract.clade(phy = tree_Marki_full, node = 305)
# data.frame(Marki_names = tree$tip.label) %>%
#   write_excel_csv(file = "results/names_Marki.csv")

# tree for Australian Meliphagoidea, n = 142 spp
tree_Mel <- names_table %>%
  filter(!is.na(Marki_name)) %>%
  .$Marki_name %>%
  keep.tip(phy = tree_Marki_full, tip = .)
# nebo:
require(magrittr)  # different pipe that explodes variables that can be used downstream
tree_Mel <- names_table %>%
  filter(!is.na(Marki_name)) %$%
  keep.tip(phy = tree_Marki_full, tip = Marki_name)

tree_Mel_lad <- ladderize(tree_Mel)  # pro obarveni subclades (funguje divne)
quartz(width = 7, height = 20)
plot(ladderize(tree_Mel), cex = 0.5)  # Australian Meliphagoidea, n = 142 spp
nodelabels(cex = 0.5)
axisPhylo()

# tree for Australian Meliphagoidea that are in ABC data, n = 111 spp
ABCspp <- data.frame(ABC_name = names(pa_mat_Mel_red))
ABCspp <- left_join(x = ABCspp, y = names_table_Mel)
tree_MelABC <- keep.tip(phy = tree_Mel, tip = ABCspp$Marki_name)
tree_MelABC_lad <- ladderize(tree_MelABC)
quartz(width = 7, height = 20)
plot(ladderize(tree_MelABC), cex = 0.5)  # Australian Meliphagoidea, n = 111 spp
nodelabels(cex = 0.5)
axisPhylo()


# +++ phylo-distance matrix -----
phylodist_mat <- cophenetic(tree_Mel)  # matice
# overil jsem, ze cophenetic matice odpovida branching times:
branch_times <- branching.times(tree_Mel_lad)  # vector; from node to tip
hist(branch_times)

# selekce phylodist matice dle stari nodu: kdyz vydelim distance v matici dvema tak dostanu stari nodu
age <- 5
dd <- phylodist_mat
dd[(dd/2) > age] <- NA
diag(dd) <- NA
sum(!is.na(dd))/2  # 89 distances btw spp in 5my clusters (overil jsem spocitanim na obarvene fylogenezi)
ss <- colSums(dd, na.rm = TRUE)
nn <- names(ss[ss>0])  # jmena druhu, kter vystupuji alespon v jednom clusteru 5my
phylodist_age_mat <- dd[rownames(dd) %in% nn, colnames(dd) %in% nn]
# sum(!is.na(phylodist_age_mat))/2  # 89 pro 5my = souhlasi
rm(dd, ss, nn)


# +++ "VCV" phylogenetic matrix -----
# it was created in the script called "vcv_matrix.R"
res_ancestors  # matrix of ancestor node IDs for spp pairs
res_heights  # matrix of heights from root of these ancestor nodes, i.e. a "vcv" matrix for the spp-pairs dataset for use in analyses
# however, analyses with this "vcv" matrix failed

# another way: choose randomly one sp from each spp pair and make vcv matrix for this reduced dataset.
# it can be done both by constructing the matrix apriori or just parsing a phylogeny to the analysis
# 1) prune the phylogeny tree_MelABC only to spp in the final dataset dataClus, but only to include always just one sp from each spp pair. Many combinations of spp from spp pairs are possible, so I take arbitrarily just two of them: sp1_Marki and sp2_Marki in data:
tree_MelABC_sp1 <- keep.tip(phy = tree_MelABC, tip = data$sp1_Marki)
tree_MelABC_sp2 <- keep.tip(phy = tree_MelABC, tip = data$sp2_Marki)


# +++ phylo-distance df -----
# +++++ all spp -----
phylodist_df <- dist2list(as.dist(phylodist_mat))
colnames(phylodist_df) <- c("sp1_Marki","sp2_Marki","phylodist")
phylodist_df <- phylodist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = FALSE) %>%
  mutate(node_age = phylodist / 2)
# +++++ sister spp -----
require(diverge)
x <- extract_sisters(tree_MelABC_lad)
y <- extract_sisters(tree_MelABC_lad)
names(y) <- c("tree_node", "sp2", "sp1")
sister_spp_df <- bind_rows(x, y, .id = "orig_df") %>%
  rename(sp1_Marki = sp1, sp2_Marki = sp2) %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = FALSE)
rm(x, y)

sister_spp_df <- sister_spp_df %>%
  left_join(x = sister_spp_df, y = phylodist_df, by = "sp_id_Marki")  %>%
  select(-sp1_Marki.y, -sp2_Marki.y) %>%
  rename(sp1_Marki = sp1_Marki.x, sp2_Marki = sp2_Marki.x) %>%
  mutate(node_age = phylodist / 2)  # uz je v phylodist_df

# calculating allopatric spp pairs:
sister_spp_df <- left_join(x = sister_spp_df, y = overlap_df)
sister_spp_df %>%
  filter(orig_df == 1, overlap == 0) %>%
  dim()


# TAKTO SE MUSE ODSTRANIT UPPER.TRI A MA TO OBSAH JAKO OBJEKT "DIST"
xx <- phylodist_age_mat
xx[upper.tri(xx)] <- NA


# +++ obarveni clades dle vekoveho kriteria -----
# musi se pouzit uz ladderizovany strom, jinak nefunguje
# +++++ subclades by age: 142 Meliphagoidea spp -----
age <- 10
nodes_for_clades <- as.integer( names(branch_times[branch_times < age]) )
spp_for_clades_List <- Descendants(x = tree_Mel_lad, node = nodes_for_clades, type = "tips")
spp_for_clades_List %>%
  map_dbl(length) %>%
  hist(., xlab = ("No. of spp in subclade"), main="")
# extract edges to color:
edges_to_color <- spp_for_clades_List %>%
  map(function(x) which.edge(phy = tree_Mel_lad, group = x)) %>%
  flatten_int() %>%
  unique()
cls <- rep("black", Nedge(tree_Mel_lad))
cls[edges_to_color] <- "magenta"
wdt <- rep(1, Nedge(tree_Mel_lad))
wdt[edges_to_color] <- 3

quartz(width = 6, height = 12)  # figure for Supplement
plot(tree_Mel_lad, cex = 0.5, edge.color = cls, edge.width = wdt)
nodelabels(node = nodes_for_clades, pch = 19, cex = 0.7)
axisPhylo()

quartz(width = 4, height = 7)
plot(tree_Mel_lad, cex = 0.5, edge.color = cls, edge.width = wdt, show.tip.label = FALSE)  # no tip labels
axisPhylo(nodelabels(node = 160, pch = 19))

# +++++ subclades by age: 111 ABC Meliphagoidea spp -----
branch_timesABC <- branching.times(tree_MelABC_lad)  # vector; from node to tip
nodes_for_clades <- as.integer( names(branch_timesABC[branch_timesABC < age]) )
spp_for_clades_List <- Descendants(x = tree_MelABC_lad, node = nodes_for_clades, type = "tips")
spp_for_clades_List %>%
  map_dbl(length) %>%
  hist(., xlab = ("No. of spp in subclade"), main="")
# extract edges to color:
edges_to_color <- spp_for_clades_List %>%
  lapply(function(x) which.edge(phy = tree_MelABC_lad, group = x)) %>%  # map() nejak prestalo fungovat
  flatten_int() %>%
  unique()
cls <- rep("black", Nedge(tree_MelABC_lad))
cls[edges_to_color] <- "magenta"
wdt <- rep(1, Nedge(tree_MelABC_lad))
wdt[edges_to_color] <- 3
# extracting nodes to highlight:
# using phytools and getMRCA failed, so I did it manually from picture:
nodes_for_clades <- c(216,178,207,186,204,201,198,126,119,128,133,134,143,146,150,149,176,174,154,160)

quartz(width = 6, height = 12)  # figure for Supplement
plot(tree_MelABC_lad, cex = 0.5, edge.color = cls, edge.width = wdt, label.offset = 0.2)
nodelabels(node = nodes_for_clades, pch = 19, cex = 0.7)
axisPhylo()



# +++++ only sister spp pairs -----
require(diverge)
sister_spp <- extract_sisters(tree_Mel_lad)
sis_spp_List <- vector(mode = "list", length = nrow(sister_spp))
for(i in seq_along(sis_spp_List)) {
  sis_spp_List[[i]] <- as.character( sister_spp[i, 2:3] )
}
edges_to_color <- sis_spp_List %>%
  map(function(x) which.edge(tree_Mel_lad, x)) %>%
  flatten_int() %>%
  unique()
cls <- rep("black", Nedge(tree_Mel_lad))
cls[edges_to_color] <- "magenta"
wdt <- rep(1, Nedge(tree_Mel_lad))
wdt[edges_to_color] <- 3
quartz(width = 6, height = 12)
plot(tree_Mel_lad, cex = 0.5, edge.color = cls, edge.width = wdt)
axisPhylo()
quartz(width = 4, height = 7)
plot(tree_Mel_lad, cex = 0.5, edge.color = cls, edge.width = wdt, show.tip.label = FALSE)  # no tip labels
axisPhylo()


##### SYMPATRY: range overlap and symmetry -----
# using file "sympatry.R" - here also the sympatry function for calculation
# vypocet overlapu viz soubor "sympatry_prace_zde.R" v podadresari "R"
# vysledky v objektu "symp_Mel"
# +++ overlap matrix -----
overlap_mat <- symp_Mel$matrix_overlap  # upper diag. matrix w/o diagonal
# vyfiltrovat jen spp ve fylog matici def. by age (ale nejdrive Marki_names):
nam <- data.frame(rows = rownames(overlap_mat), cols = colnames(overlap_mat))
identical(nam$rows, nam$cols)  # TRUE = OK
nam <- left_join( x = nam, y = names_table, by = c("rows" = "BLI6_name") )
sum(is.na(nam))  # 0 = OK, vsechna jmena prirazena
rm(nam)

# +++ overlap matrix based on subclade age -----
overlap_age_mat <- overlap_mat
rownames(overlap_age_mat) <- nam$Marki_name
colnames(overlap_age_mat) <- nam$Marki_name
# upper diag. matrix w/o diagonal !!!:
overlap_age_mat <- overlap_age_mat[rownames(overlap_age_mat) %in% rownames(phylodist_age_mat), colnames(overlap_age_mat) %in% colnames(phylodist_age_mat)]

# +++ overlap df -----
overlap_df <- dist2list(as.dist(t(overlap_mat)))
colnames(overlap_df) <- c("sp1_BLI6","sp2_BLI6","overlap")
overlap_df <- overlap_df %>%
  unite(col = sp_id_BLI6, sp1_BLI6, sp2_BLI6, sep = "_", remove = FALSE)
# add Marki names and Marki ID:
overlap_df <- addSpID(target.df = overlap_df, names.df = names_table, guide.name = "BLI6", added.col = "Marki")

# +++ symmetry matrix -----
symmetry_mat <- symp_Mel$matrix_symmetry  # upper diag. matrix w/o diagonal
# making sure names are BLI6:
nam <- data.frame(rows = rownames(symmetry_mat), cols = colnames(symmetry_mat))
identical(nam$rows, nam$cols)  # TRUE = OK
nam <- left_join( x = nam, y = names_table, by = c("rows" = "BLI6_name") )
sum(is.na(nam))  # 0 = OK, vsechna jmena prirazena
rm(nam)

# +++ symmetry df -----
symmetry_df <- dist2list(as.dist(t(symmetry_mat)))
colnames(symmetry_df) <- c("sp1_BLI6","sp2_BLI6","symmetry")
symmetry_df <- symmetry_df %>%
  unite(col = sp_id_BLI6, sp1_BLI6, sp2_BLI6, sep = "_", remove = FALSE)
# add Marki names and Marki ID:
symmetry_df <- addSpID(target.df = symmetry_df, names.df = names_table, guide.name = "BLI6", added.col = "Marki")


##### SYNTOPY: co-occurrence from ABC data -----
# +++ read ABC data -----
pa <- read_csv("/Users/vladimirremes/Documents/clanky_a_projekty/GACR_Australie/PROJEKT_ABCdata/results/pres_abs_locs_final_LH.csv")[, -1]  # v prvnim sloupci v csv je nove poradove cislo, ktere vetsinou nepotrebuji
table(colSums(pa[, 4:ncol(pa)], na.rm = T))  # 60x nula ze 308 druhu = druh neni ani na jednom miste
# cols_to_drop = c( rep(TRUE, 3), colSums(pa[, 4:ncol(pa)]) > 0 )
# pa <- pa[, cols_to_drop]  # vysledna PA matrix bez nulovych druhu
# rm(cols_to_drop)
# sort(colSums(pa))
ABC_name <- data.frame( ABC_names = colnames(pa)[-c(1:3)] )
setdiff(ABC_name$ABC_names, names_table$BLI_gis2011_name)  # 1 rozdil: H. albispec x cinereifrons (a exoti, kteri v ABC nejsou)
pa_mat_Mel <- pa %>%
  select(contains(names_table_Mel$ABC_name))
table(colSums(pa_mat_Mel[, 1:ncol(pa_mat_Mel)], na.rm = T))  # 31x nula ze 142 spp = ani na jednom miste
pa_mat_Mel_red <- pa_mat_Mel %>%
  select(where(function(x) sum(x) > 0))  # exclude spp with 0 incidence
rownames(pa_mat_Mel_red) <- pa$Site_code
pa_mat_Mel_red_t <- t(pa_mat_Mel_red)

# +++ naive cooccur -----
# Ale tohle je uplne SPATNE, protoze bere pro druhove pary i lokality mimo ranges obou druhu! Musi se vzdy vzit jen lokality lezici ve sjednoceni (union) danych dvou arealu.
require(cooccur)
# bud s thresh=TRUE (cooc_naive_data_df) nebo FALSE (cooc_naive_data_ALL_df)
cooc_naive_Mel <- cooccur(mat = pa_mat_Mel_red_t, type = "spp_site", thresh = TRUE, spp_names = TRUE)
cooc_naive_Mel_ALL <- cooccur(mat = pa_mat_Mel_red_t, type = "spp_site", thresh = FALSE, spp_names = TRUE)
cooc_res_df <- cooc_naive_Mel$results
cooc_res_df <- cooc_naive_Mel_ALL$results
names(cooc_res_df)[10:11] <- c("sp1_ABC", "sp2_ABC")
cooc_res_df <- cooc_res_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

# +++++ summaries and graphics on cooc object -----
cooc <- cooc_naive_Mel
summary(cooc)
cooc_pairs <- pair.attributes(cooc)  # sumarizace asociaci pro kazdy druh, n = 111 spp
# plotting
require(ggplot2)
plot(cooc)  # association matrix (jen spp se signif. asociaci; od nejvic negat po nejvic pozit)
pair.profile(cooc)  # graficka sumarizace spp asociaci (tedy dat z pair.attributes funkce), all spp.
obs.v.exp(cooc) + theme_bw()

# +++++ cooc effect sizes -----
# effect size = diff btw expected and observed frequency of co-occurrence (can be standardized by dividing by N - then btw -1 and 1, where positive indicate positive associations and negative vice versa)
effect.sizes(cooc)  # ale vyloucili jsme (thresh=TRUE) ty s expected cooc<1
# pro vsechny efekty analyzu s thresh=FALSE, a pro zrychleni bez signifikanci (only_effects=TRUE):
# standardizace eff_standard=TRUE
cooc_effs_mat <- cooccur(mat = pa_mat_Mel_red_t, type = "spp_site", thresh = TRUE, spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)  # efekty v matici
cooc_effs_df <- cooccur(mat = pa_mat_Mel_red_t, type = "spp_site", thresh = FALSE, spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = FALSE)  # efekty jako pairwise table

names(cooc_effs_df)[1:2] <- c("sp1_ABC", "sp2_ABC")
cooc_effs_df <- cooc_effs_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

# +++++ creating FULL cooc data -----
cooc_naive_data_df <- full_join( x = cooc_res_df, y = select(cooc_effs_df, c(1,4)), by = "sp_id_ABC")
cooc_naive_data_ALL_df <- full_join( x = cooc_res_df, y = select(cooc_effs_df, c(1,4)), by = "sp_id_ABC")
rm(cooc_effs_df, cooc_res_df)
# add Marki names and Marki ID:
cooc_naive_data_df <- cooc_naive_data_df %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp1_ABC" = "ABC_name")) %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp2_ABC" = "ABC_name")) %>%
  unite(col = sp_id_Marki, Marki_name.x, Marki_name.y, sep = "_", remove = FALSE) %>%
  rename(sp1_Marki = Marki_name.x, sp2_Marki = Marki_name.y)


# +++ range-constrained cooccur in range INTERSECTION -----
# Toto je mask cooccur jen v intersection of ranges = mask Intersect, ale zde se zve proste "mask".

# +++++ read pools -----
pools <- read_csv("/Users/vladimirremes/Documents/clanky_a_projekty/GACR_Australie/PROJEKT_ABCdata/results/pres_abs_pools_10buff_final_LH.csv")[, -1]  # v prvnim sloupci v csv je nove poradove cislo, ktere vetsinou nepotrebuji
table(colSums(pools[, 4:ncol(pools)], na.rm = T))  # 20x nula ze 308 druhu = druh neni v poolu ani na jednom miste
# cols_to_drop = c( rep(TRUE, 3), colSums(pools[, 4:ncol(pools)]) > 0 )
# pools <- pools[, cols_to_drop]  # vysledna PA pools matrix bez nulovych druhu
# rm(cols_to_drop)
# sort(colSums(pools))
ABC_name <- data.frame( ABC_names = colnames(pools)[-c(1:3)] )
setdiff(ABC_name$ABC_names, names_table$BLI6_name)  # 1 rozdil: H. albispec x cinereifrons (a exoti, kteri v ABC nejsou)
pools_mat_Mel <- pools %>%
  select(contains(names_table_Mel$ABC_name))
min(pools_mat_Mel - pa_mat_Mel)  # min = 0, tj. druh je vzdy v poolu tam kde byl detekovan = OK
table(colSums(pools_mat_Mel[, 1:ncol(pools_mat_Mel)], na.rm = T))  # 10x nula ze 142 spp = druh neni v poolu ani na jednom miste
pools_mat_Mel_red <- pools_mat_Mel %>%
  select(contains(colnames(pa_mat_Mel_red)))
rownames(pools_mat_Mel_red) <- pa$Site_code
pools_mat_Mel_red_t <- t(pools_mat_Mel_red)

# +++++ fit mask cooccur -----
source("/Users/vladimirremes/Documents/stat_mat_a_program/R/R_functions/FUNCTIONS_spatial/cooccur/my_cooccur_functions.R")

cooc_mask_Mel <- my_cooccur(mat = pa_mat_Mel_red_t, type="spp_site", thresh = TRUE, spp_names = TRUE, site_mask = pools_mat_Mel_red_t, min_sites=20, prob="hyper", true_rand_classifier=0.1)
summary(cooc_mask_Mel)

cooc_res_df <- cooc_mask_Mel$results
names(cooc_res_df)[11:12] <- c("sp1_ABC", "sp2_ABC")
cooc_res_df <- cooc_res_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

cooc_effs_df <- my_effect.sizes(mod = cooc_mask_Mel, standardized = TRUE, matrix = FALSE)
names(cooc_effs_df)[1:2] <- c("sp1_ABC", "sp2_ABC")
cooc_effs_df <- cooc_effs_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

# +++++ creating FULL cooc data -----
cooc_mask_data_df <- full_join( x = cooc_res_df, y = select(cooc_effs_df, c(1,4)), by = "sp_id_ABC")
rm(cooc_effs_df, cooc_res_df)
# add Marki names and Marki ID:
cooc_mask_data_df <- cooc_mask_data_df %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp1_ABC" = "ABC_name")) %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp2_ABC" = "ABC_name")) %>%
  unite(col = sp_id_Marki, Marki_name.x, Marki_name.y, sep = "_", remove = FALSE) %>%
  rename(sp1_Marki = Marki_name.x, sp2_Marki = Marki_name.y)


# +++ range-constrained cooccur in range UNION -----
# Naive cooccur fituje vsech 470 mist i pro pary druhu, kdy na nekterych mistech nani ani jeden z druhoveho paru, coz je blbost. Spravny naive cooccur (tomu budu rikat UNION cooccur) musi vzdy pro druhovy par vybrat jen ta mista, kde se muze dle PA matice vyskytovat ALESPON jeden druh z analyzovaneho paru.

# +++++ fit union cooccur -----
# vypocet delam ve zvlastnim souboru "union_cooccur.R" a je ulozen v objektu "union_cooccur.RData". Duvodem je, ze nejde vytvorit N_matrix pro UNION cooccur, takze to musim delat rucne v loopu a to by zde zabiralo moc mista.
load("union_cooccur.RData")  # toto uz je v objektu "cooccurrence_MS.RData", takze se zde nemusi loadovat

# +++++ creating FULL cooc data -----
names(cooc_res_ALL_df)[11:12] <- c("sp1_ABC", "sp2_ABC")
cooc_res_ALL_df <- cooc_res_ALL_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

names(cooc_effs_ALL_df)[1:2] <- c("sp1_ABC", "sp2_ABC")
cooc_effs_ALL_df <- cooc_effs_ALL_df %>%
  unite(col = sp_id_ABC, sp1_ABC, sp2_ABC, sep = "_", remove = FALSE)

cooc_union_data_ALL_df <- full_join( x = cooc_res_ALL_df, y = select(cooc_effs_ALL_df, c(1,4)), by = "sp_id_ABC")
rm(cooc_effs_ALL_df, cooc_res_ALL_df)
# add Marki names and Marki ID:
cooc_union_data_ALL_df <- cooc_union_data_ALL_df %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp1_ABC" = "ABC_name")) %>%
  left_join( ., select(names_table, Marki_name, ABC_name), by = c("sp2_ABC" = "ABC_name")) %>%
  unite(col = sp_id_Marki, Marki_name.x, Marki_name.y, sep = "_", remove = FALSE) %>%
  rename(sp1_Marki = Marki_name.x, sp2_Marki = Marki_name.y)



##### ECOLOGY -----
ecol <- read_csv("/Users/vladimirremes/Documents/clanky_a_projekty/funkcni znaky a diverzita/DATA-Australie/data HANZAB - EK/data pro ABC projekt-specializace/data_specialization_final_LH.csv")[, -1]  # 308 rows
ecol <- ecol %>%
  unite(col = HANZAB_name, Genus_Hanzab, Species_Hanzab, sep = "_", remove = TRUE) %>%
  select(-1) %>%  # odstranen sloupec "Species" - co to bylo za jmeno?
  arrange(HANZAB_name)  # seradit dle anglicke abecedy

# overeni HANZAB names - kontrola jmen (opravil jsem nesrovnalosti, 5.8.2021):
setdiff(ecol$HANZAB_name, names_table$HANZAB_name)  # nic
setdiff(names_table$HANZAB_name, ecol$HANZAB_name)  # exoti

# create Meliphagoidea file with Marki names:
ecol_Mel <- names_table %>%
  filter(!is.na(Marki_name)) %>%
  select(HANZAB_Family, HANZAB_name, Marki_name) %>%
  left_join(x = ., y = ecol, by = "HANZAB_name") %>%
  select(-4) %>%
  arrange(Marki_name)
  
ecol_Mel %>% count(HANZAB_Family)

# +++ habitat -----
require(vegan)
class(ecol_Mel) <- "data.frame"
rownames(ecol_Mel) <- ecol_Mel$Marki_name
habdist <- vegdist(log1p(ecol_Mel[, 25:34]), method="bray")
habdist_df <- dist2list(habdist)
colnames(habdist_df) <- c("sp1_Marki","sp2_Marki","habdist")
habdist_df <- habdist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

# +++ foraging -----
# 4-7: stratum
# 8-15: substrate
# 16-24: behavior
# 25-34: habitat
# 35-42: food
require(vegan)
ecodist.str <- vegdist(ecol_Mel[, c(4:7)], method="bray")
ecodist.sub <- vegdist(ecol_Mel[, c(8:15)], method="bray")
ecodist.beh <- vegdist(ecol_Mel[, c(16:24)], method="bray")
ecodist.fud <- vegdist(ecol_Mel[, c(35:42)], method="bray")
ecodist.hab <- vegdist(ecol_Mel[, c(25:34)], method="bray")

foragdist <- vegdist(log1p(ecol_Mel[, c(4:24,35:42)]), method="bray")
foragdist_df <- dist2list(foragdist)
colnames(foragdist_df) <- c("sp1_Marki","sp2_Marki","foragdist")
foragdist_df <- foragdist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

stratumdist <- vegdist(log1p(ecol_Mel[, c(4:7)]), method="bray")
stratumdist_df <- dist2list(stratumdist)
colnames(stratumdist_df) <- c("sp1_Marki","sp2_Marki","stratumdist")
stratumdist_df <- stratumdist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

substratedist <- vegdist(log1p(ecol_Mel[, c(8:15)]), method="bray")
substratedist_df <- dist2list(substratedist)
colnames(substratedist_df) <- c("sp1_Marki","sp2_Marki","substratedist")
substratedist_df <- substratedist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

methoddist <- vegdist(log1p(ecol_Mel[, c(16:24)]), method="bray")
methoddist_df <- dist2list(methoddist)
colnames(methoddist_df) <- c("sp1_Marki","sp2_Marki","methoddist")
methoddist_df <- methoddist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

fooddist <- vegdist(log1p(ecol_Mel[, c(35:42)]), method="bray")
fooddist_df <- dist2list(fooddist)
colnames(fooddist_df) <- c("sp1_Marki","sp2_Marki","fooddist")
fooddist_df <- fooddist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)

pairs(cbind(habdist_df$habdist, foragdist_df$foragdist, stratumdist_df$stratumdist, substratedist_df$substratedist, methoddist_df$methoddist, fooddist_df$fooddist))


##### MORPHOLOGY -----
morph <- read_csv("/Users/vladimirremes/Documents/clanky_a_projekty/funkcni znaky a diverzita/DATA-Australie/morfologie_Au.csv")[, -1]  # 314 rows (jsou i "nove" druhy)
# kontrola jmen:
setdiff(morph$species_Pizzey, names_table$Pizzey_name)  # jen "nove" druhy
setdiff(names_table$Pizzey_name, morph$species_Pizzey)  # exoti
morph <- morph %>%
  select(-9) %>%  # odstranit n-ko
  rename(tarsus = tarsus_mm, mass = mass_g) %>%
  arrange(species_Pizzey)  # seradit dle anglicke abecedy

# create Meliphagoidea file with Marki names:
morph_Mel <- names_table %>%
  filter(!is.na(Marki_name)) %>%
  select(HANZAB_Family, Pizzey_name, Marki_name) %>%
  left_join( x = ., y = morph, by = c("Pizzey_name" = "species_Pizzey") ) %>%
  arrange(Marki_name)

morph_Mel %>% count(HANZAB_Family)

class(morph_Mel) <- "data.frame"
rownames(morph_Mel) <- morph_Mel$Marki_name
morphdist <- dist(x = log10(morph_Mel[, 4:10]), method = "euclidean")
morphdist_df <- dist2list(morphdist)

colnames(morphdist_df) <- c("sp1_Marki","sp2_Marki","morphdist")
morphdist_df <- morphdist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)


##### CLIMATE -----
clim_Mel <- read_csv("data/climate/Meliphagoidea_Bioclim_210813.csv")[, -1]  # file od Lenky, srpen 2021
# Bioclim_1: Annual Mean Temperature
# Bioclim_4: Temperature Seasonality (standard deviation ×100)
# Bioclim_5: Max Temperature of Warmest Month
# Bioclim_6: Min Temperature of Coldest Month
# Bioclim_7: Temperature Annual Range (BIO5-BIO6)
# Bioclim_12: Annual Precipitation
# Bioclim_13: Precipitation of Wettest Month
# Bioclim_14: Precipitation of Driest Month
# Bioclim_15: Precipitation Seasonality (Coefficient of Variation)
# min vsech promennych > 0 => lze pouzit log-transformaci (hl. u precipitation vars)
sum(is.na(clim_Mel))  # 0, nic nechybi = OK
clim_Mel <- clim_Mel %>%
  separate(col = Species, into = c("genus", "species")) %>%
  unite(col = BLI6_name, genus, species, sep = "_", remove = TRUE)
clim_Mel <- names_table_Mel %>%
  select(BLI6_name, Marki_name) %>%
  left_join(x = clim_Mel, y = ., by = "BLI6_name") %>%
  relocate(Marki_name, .after = BLI6_name)
require(GGally)
ggpairs(clim_Mel[, c(3:11)])
hist(clim_Mel[, 11])  # Bioclim 12-15 log10-transform (all have min > 0)

class(clim_Mel) <- "data.frame"
rownames(clim_Mel) <- clim_Mel$Marki_name

climdist <- dist(x = cbind( clim_Mel[, c(4:7)], log10(clim_Mel[, 8:11]) ), method = "euclidean")
climdist_df <- dist2list(climdist)
colnames(climdist_df) <- c("sp1_Marki","sp2_Marki","climdist")
climdist_df <- climdist_df %>%
  unite(col = sp_id_Marki, sp1_Marki, sp2_Marki, sep = "_", remove = TRUE)



##### RELATIONSHIPS between variables -----
require(visreg)
# prepare overall data with all variables:
cutoff_age <- 10  # max node_age is 31.32601
# + create data: naive cooccur only -----
d_age_naive <- left_join(cooc_naive_data_df, overlap_df, by = "sp_id_Marki") %>%
  left_join(., symmetry_df, by = "sp_id_Marki") %>%
  left_join(., phylodist_df, by = "sp_id_Marki") %>%
  filter(node_age <= cutoff_age)

# + create data: union cooccur only -----
d_age_union <- left_join(cooc_union_data_df, overlap_df, by = "sp_id_Marki") %>%
  left_join(., symmetry_df, by = "sp_id_Marki") %>%
  left_join(., phylodist_df, by = "sp_id_Marki") %>%
  filter(node_age <= cutoff_age)

# + create data: mask cooccur only -----
d_age_mask <- left_join(cooc_mask_data_df, overlap_df, by = "sp_id_Marki") %>%
  left_join(., symmetry_df, by = "sp_id_Marki") %>%
  left_join(., phylodist_df, by = "sp_id_Marki") %>%
  filter(node_age <= cutoff_age)

# + create data: add SES effects (see also below under "Z-scores") -----
d_age_mask <- d_age_mask %>%
  mutate(Var = ((sp1_inc*sp2_inc)/no_sites)*((no_sites-sp1_inc)/no_sites)*((no_sites-sp2_inc)/(no_sites-1)), SES = (obs_cooccur - exp_cooccur) / sqrt(Var) )

cooc_union_data_ALL_df <- cooc_union_data_ALL_df %>%
  mutate(Var = ((sp1_inc*sp2_inc)/no_sites)*((no_sites-sp1_inc)/no_sites)*((no_sites-sp2_inc)/(no_sites-1)), SES = (obs_cooccur - exp_cooccur) / sqrt(Var) )

# + create data: concatenate mask and naive/union cooccur, AGE criterion -----
dd <- left_join(x = d_age_mask, y = cooc_naive_data_ALL_df, by = "sp_id_Marki")
# NA hodnoty vznikaji tim, ze threshold odstranuje jine pary v naive vs. mask vypoctu; potvrzeno tak, ze kdys se spocita naive cooccur s thresh=FALSE, tak NA hodnoty nevzniknou
dd <- filter(dd, !is.na(effects.y))  # jen pokud NA hodnoty - viz vyse

dd_union <- left_join(x = d_age_mask, y = cooc_union_data_ALL_df, by = "sp_id_Marki")
dd_union <- filter(dd_union, !is.na(effects.y))  # jen pokud NA hodnoty - viz vyse

data_10my <- dd %>%
  select(sp_id_ABC.x, sp_id_Marki, sp_id_BLI6.x, prob_cooccur.x, effects.x, prob_cooccur.y, effects.y, overlap, symmetry, phylodist, node_age) %>%
  rename(sp_id_ABC = sp_id_ABC.x, sp_id_BLI6 = sp_id_BLI6.x, prob_cooc_mask = prob_cooccur.x, eff_mask = effects.x, prob_cooc_naive = prob_cooccur.y, eff_naive = effects.y) %>%
  left_join(., morphdist_df, by = "sp_id_Marki") %>%
  left_join(., habdist_df, by = "sp_id_Marki")  %>%
  left_join(., foragdist_df, by = "sp_id_Marki") %>%
  left_join(., climdist_df, by = "sp_id_Marki")

data_10my_union <- dd_union %>%
  select(sp_id_ABC.x, sp_id_Marki, sp1_Marki, sp2_Marki, sp_id_BLI6.x, no_sites.x, prob_cooccur.x, effects.x, p_lt.x, p_gt.x, SES.x, no_sites.y, prob_cooccur.y, effects.y, p_lt.y, p_gt.y, SES.y, overlap, symmetry, phylodist, node_age) %>%
  rename(sp_id_ABC = sp_id_ABC.x, sp_id_BLI6 = sp_id_BLI6.x, no_sites_mask = no_sites.x, prob_cooc_mask = prob_cooccur.x, eff_mask = effects.x, p_lt_mask = p_lt.x, p_gt_mask = p_gt.x, SES_mask = SES.x, no_sites_union = no_sites.y, prob_cooc_union = prob_cooccur.y, eff_union = effects.y, p_lt_union = p_lt.y, p_gt_union = p_gt.y, SES_union = SES.y) %>%
  left_join(., morphdist_df, by = "sp_id_Marki") %>%
  left_join(., habdist_df, by = "sp_id_Marki")  %>%
  left_join(., foragdist_df, by = "sp_id_Marki") %>%
  left_join(., climdist_df, by = "sp_id_Marki") %>%
  left_join(., fooddist_df, by = "sp_id_Marki") %>%
  left_join(., methoddist_df, by = "sp_id_Marki") %>%
  left_join(., substratedist_df, by = "sp_id_Marki") %>%
  left_join(., stratumdist_df, by = "sp_id_Marki")

data <- data_10my
data <- data_10my_union
identical(data_10my, data)  # TRUE = OK
identical(data_10my_union, data)  # TRUE = OK

# + create data: narrow format (from wide) -----
# for naive (this code probably no longer works; it is OK for union below):
data2_prob <- data %>%
  select(-c(sp_id_ABC, sp_id_BLI6, eff_mask, eff_naive)) %>%
  pivot_longer(cols = c(2,3), names_to = "analysis_type", values_to = "prob_cooccur")

data2_effs <- data %>%
  select(-c(sp_id_ABC, sp_id_BLI6, prob_cooc_mask, prob_cooc_naive)) %>%
  pivot_longer(cols = c(2,3), names_to = "analysis_type", values_to = "effects")

# for union
data2_prob <- data %>%
  select(-c(1,3:5,8,11,14,17)) %>%
  pivot_longer(cols = c(3,7), names_to = "analysis_type", values_to = "prob_cooccur")

data2_effs <- data %>%
  select(-c(1,3:5,7,11,13,17)) %>%
  pivot_longer(cols = c(3,7), names_to = "analysis_type", values_to = "effects")

data2_SES <- data %>%
  select(-c(1,3:5,7,8,13,14)) %>%
  pivot_longer(cols = c(5,9), names_to = "analysis_type", values_to = "SES")


# +++ comparing two different sympatry indexes -----
# new one is Overlap.both that is the proportion of overlap out of the union of areas of the two ranges
data_Overlap.both <- read.csv(file = "results/data_syntopy_Overlap.both.csv")
require(GGally)
ggpairs(data_Overlap.both[, c(20,10,9,21)], diag = list(continuous = wrap("barDiag", bins = 15, color = "black")), lower = list(continuous = "smooth" ), title = "Range sympatry and symmetry")


# +++ comparing naive vs mask cooccur -----
quartz(width = 6, height = 6)  # 1 panel
quartz(width = 10, height = 10)  # 4 panels
par(mfrow=c(2,2))
par(mar=c(7,7,4,2))
require(visreg)

# +++++ EFFECTS: -----
#visreg( lm(effects.y ~ poly(effects.x, 2), dd), xvar="effects.x", points=list(cex=0.8), xlab="Co-occurrence in overlap", ylab="Co-occurrence", cex.lab=2, main="Effect size (standardized)" )
abline(coef = c(0,1), lwd=1)

# panel A for figure:
visreg( lm(eff_union ~ poly(eff_mask, 1), data), xvar="eff_mask", points=list(cex=0.8), xlab="Co-occurrence in overlap", ylab="Co-occurrence", cex.lab=2, main="Effect size (standardized)" )
abline(coef = c(0,1), lwd=1)

# +++++ PROBS: ------
#visreg( lm(prob_cooccur.y ~ poly(prob_cooccur.x, 2), dd), xvar="prob_cooccur.x", points=list(cex=0.8), xlab="Co-occurrence in overlap", ylab="Co-occurrence", cex.lab=2, main="Probability of co-occurrence" )
abline(coef = c(0,1), lwd=1)

visreg( lm(prob_cooc_union ~ poly(prob_cooc_mask, 1), data), xvar="prob_cooc_mask", points=list(cex=0.8), xlab="Co-occurrence in overlap", ylab="Co-occurrence", cex.lab=2, main="Probability of co-occurrence" )
abline(coef = c(0,1), lwd=1)

# +++++ SIGNIFICANCES: -----
# negative
# panel C for figure
visreg( lm(p_lt_union ~ poly(p_lt_mask, 1), data), xvar="p_lt_mask", points=list(cex=0.8), xlab="Significance in overlap", ylab="Significance", cex.lab=2, main="Signif. of negative co-occurrence" )
abline(h = 0.05, lwd=1)
abline(v = 0.05, lwd=1)

# positive
# panel D for figure
visreg( lm(p_gt_union ~ poly(p_gt_mask, 1), data), xvar="p_gt_mask", points=list(cex=0.8), xlab="Significance in overlap", ylab="Significance", cex.lab=2, main="Signif. of positive co-occurrence" )
abline(h = 0.05, lwd=1)
abline(v = 0.05, lwd=1)

# diff in signif on overlap:
data_temp <- data_temp %>% mutate(p_lt_diff = p_lt_mask - p_lt_union)
data_temp <- data_temp %>% mutate(p_gt_diff = p_gt_mask - p_gt_union)

visreg( lm(p_lt_diff ~ poly(overlap, 1), data_temp), xvar="overlap", points=list(cex=0.8), xlab="Overlap", ylab="Diff. in significance (mask - union)", cex.lab=2, main="Signif. of negative co-occurrence" )
abline(h = 0)
summary(lm(p_lt_diff ~ poly(overlap, 1), data_temp))

visreg( lm(p_gt_diff ~ poly(overlap, 1), data_temp), xvar="overlap", points=list(cex=0.8), xlab="Overlap", ylab="Diff. in significance (mask - union)", cex.lab=2, main="Signif. of positive co-occurrence" )
abline(h = 0)
summary(lm(p_gt_diff ~ poly(overlap, 1), data_temp))


# +++++ NO OF SITES: -----
# panel B for figure
visreg( lm(no_sites_union ~ poly(no_sites_mask, 1), data), xvar="no_sites_mask", points=list(cex=0.8), xlab="No. of sites in overlap", ylab="No. of sites", cex.lab=2, main="No. of study sites", xlim=c(0,500), ylim=c(0,500) )
abline(coef = c(0,1), lwd=1)




###  DODELAT MAJOR AXIS REGRESSION PRO TATO SROVNANI !!!





# biases vs. difference in no. of sites:
data_temp <- data
data_temp$no_sites_diff <- data$no_sites_mask - data$no_sites_union
visreg( lm( I(eff_union-eff_mask) ~ poly(no_sites_diff, 1), data_temp), xvar="no_sites_diff", points=list(cex=0.8), xlab="Difference in no. of sites", ylab="Bias", cex.lab=2, main="Effect size (standardized)" )  # neni signif.
abline(h = 0)

visreg( lm( I(prob_cooc_union-prob_cooc_mask) ~ poly(no_sites_diff, 1), data_temp), xvar="no_sites_diff", points=list(cex=0.8), xlab="Difference in no. of sites", ylab="Bias", cex.lab=2, main="Prob. of co-occurrence (standardized)" )  # vysoce signif.
abline(h = 0)

# no. of sites and overlap:
visreg( lm( no_sites_union ~ poly(overlap, 1), data_temp), xvar="overlap", points=list(cex=0.8), cex.lab=2 )
visreg( lm( no_sites_mask ~ poly(overlap, 1), data_temp), xvar="overlap", points=list(cex=0.8), cex.lab=2 )

visreg( lm( no_sites_diff ~ poly(overlap, 1), data_temp), xvar="overlap", points=list(cex=0.8), cex.lab=2, ylab="Decrease in no. sites in overlap", xlab="Range overlap" )

cor.test(data_temp$no_sites_diff, data_temp$overlap)  # r = 0.57, p < 0.001



# +++ overlap ON phylodist -----
# +++++ subclades -----
# from matrices:
# phylodist <- orderMatrixAlphab(mat = phylodist_age_mat)
# overlap <- orderMatrixAlphab(mat = overlap_age_mat)
# phylodist <- phylodist[upper.tri(phylodist)]/2
# overlap <- overlap[upper.tri(overlap)]
# df <- data.frame(age = phylodist, overlap = overlap)
# ggplot(df, aes(x = age, y = overlap)) + 
#   geom_point() +
#   geom_smooth()
ggplot(d_age_naive, aes(x = node_age, y = overlap)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Node age", y = "Range overlap")

visreg( lm(overlap~poly(node_age,2), data = d_age_naive), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Range overlap", cex.lab=2, main="Meliphagoidea" )


# +++++ sister spp pairs -----
# only from dfs
# first, prepare sp_id based on BLI6_names:
sister_spp_df <- addSpID(target.df = sister_spp_df, names.df = names_table, guide.name = "Marki", added.col = "BLI6")

# then plot:
sister_spp_df %>%
  left_join(., overlap_df, by = "sp_id_BLI6") %>%
  filter(orig_df == "1")  %>%
  ggplot(aes(x = node_age, y = overlap)) +
  geom_point() +
  geom_smooth(method = lm)



# +++ cooccur ON phylodist -----
# +++++ naive cooccur -----
ggplot(d_age_naive, aes(x = node_age, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Node age", y = "Local co-occurrence")

quartz(width = 6, height = 6)
par(mar=c(7,7,4,2))
# ALL data:
visreg( lm(effects~poly(node_age,2), data = d_age_naive), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Local co-occurrence", cex.lab=2, main="Naive analysis" )

visreg( lm(prob_cooccur~poly(node_age,2), data = d_age_naive), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Prob. of co-occurrence", cex.lab=2, main="Naive analysis" )

# ONLY with mask data:
visreg( lm(eff_naive~poly(node_age,2), data = data), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Co-occurrence", cex.lab=2, main="Effect size" )

visreg( lm(prob_cooc_naive~poly(node_age,2), data = data), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Co-occurrence", cex.lab=2, main="Probability of co-occurrence" )


# +++++ mask cooccur -----
ggplot(d_age_mask, aes(x = node_age, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Node age", y = "Local co-occurrence")

visreg( lm(effects~poly(node_age,2), data = d_age_mask), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Local co-occurrence", cex.lab=2, main="Mask analysis" )

visreg( lm(prob_cooccur~poly(node_age,2), data = d_age_mask), xvar="node_age", points=list(cex=1), xlab="Node age", ylab="Prob. of co-occurrence", cex.lab=2, main="Mask analysis" )


# +++++ both in one plot -----
# Nnaive > Nmask
d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = node_age, y = effects, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Node age", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = node_age, y = prob_cooccur, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Node age", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nunion > Nmask
d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = node_age, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Node age", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = node_age, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Node age", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nnaive = Nmask, Nunion = Nmask podle toho, z jakych objektu vyrobim data2_effs vyse
ggplot(data = data2_effs, mapping = aes(x = node_age, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Node age", y = "Local co-occurrence", title = "Effect size")

ggplot(data = data2_prob, mapping = aes(x = node_age, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Node age", y = "Local co-occurrence", title = "Probability of co-occurrence")



# +++ cooccur ON overlap -----
# +++++ naive cooccur -----
ggplot(d_age_naive, aes(x = overlap, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Range overlap", y = "Local co-occurrence")

quartz(width = 6, height = 6)
par(mar=c(7,7,4,2))
visreg( lm(effects~poly(overlap,2), data = d_age_naive), xvar="overlap", points=list(cex=1), xlab="Range overlap", ylab="Co-occurrence", cex.lab=2, main="Effect size (standardized)" )
abline(h = 0)

visreg( lm(prob_cooccur~poly(overlap,2), data = d_age_naive), xvar="overlap", points=list(cex=1), xlab="Range overlap", ylab="Co-occurrence", cex.lab=2, main="Probability of co-occurrence" )
abline(h = 0)

# +++++ mask cooccur -----
ggplot(d_age_mask, aes(x = overlap, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Range overlap", y = "Local co-occurrence", title = "Meliphagoidea")

visreg( lm(effects~poly(overlap,2), data = d_age_mask), xvar="overlap", points=list(cex=1), xlab="Range overlap", ylab="Co-occurrence in overlap", cex.lab=2, main="Effect size (standardized)", xlim=c(0,100) )
abline(h = 0)

visreg( lm(prob_cooccur~poly(overlap,2), data = d_age_mask), xvar="overlap", points=list(cex=1), xlab="Range overlap", ylab="Co-occurrence in overlap", cex.lab=2, main="Probability of co-occurrence" )
abline(h = 0)

# +++++ both in one plot -----
# Nnaive > Nmask
d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = overlap, y = effects, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Range overlap", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = overlap, y = prob_cooccur, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Range overlap", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nunion > Nmask
d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = overlap, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range overlap", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = overlap, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range overlap", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nnaive = Nmask, Nunion = Nmask podle toho, z jakych objektu vyrobin data2_effs vyse
ggplot(data = data2_effs, mapping = aes(x = overlap, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range overlap", y = "Local co-occurrence", title = "Effect size")

ggplot(data = data2_prob, mapping = aes(x = overlap, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range overlap", y = "Local co-occurrence", title = "Probability of co-occurrence")



# +++ cooccur ON symmetry -----
# +++++ naive cooccur -----
ggplot(d_age_naive, aes(x = symmetry, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Range symmetry", y = "Local co-occurrence")

quartz(width = 6, height = 6)
par(mar=c(7,7,4,2))
visreg( lm(effects~poly(symmetry,2), data = d_age_naive), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="Co-occurrence", cex.lab=2, main="Effect size (standardized)" )
abline(h = 0)

visreg( lm(prob_cooccur~poly(symmetry,2), data = d_age_naive), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="Co-occurrence", cex.lab=2, main="Probability of co-occurrence" )
abline(h = 0)

# +++++ mask cooccur -----
ggplot(d_age_mask, aes(x = symmetry, y = effects)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,2)) +
  labs(x = "Range symmetry", y = "Local co-occurrence", title = "Meliphagoidea")

visreg( lm(effects~poly(symmetry,2), data = d_age_mask), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="Co-occurrence in overlap", cex.lab=2, main="Effect size (standardized)", xlim=c(0,50) )
abline(h = 0)

visreg( lm(prob_cooccur~poly(symmetry,2), data = d_age_mask), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="Co-occurrence in overlap", cex.lab=2, main="Probability of co-occurrence" )
abline(h = 0)

# +++++ both in one plot -----
# Nnaive > Nmask
d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = symmetry, y = effects, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Range symmetry", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = symmetry, y = prob_cooccur, color = analysis_type)) +
    geom_point() + 
    geom_smooth(method = lm, formula=y~poly(x,1)) +
    labs(x = "Range symmetry", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nunion > Nmask
d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = symmetry, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range symmetry", y = "Local co-occurrence", title = "Effect size")

d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" ) %>%
  ggplot(mapping = aes(x = symmetry, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range symmetry", y = "Local co-occurrence", title = "Probability of co-occurrence")

# Nnaive = Nmask, Nunion = Nmask podle toho, z jakych objektu vyrobin data2_effs vyse
ggplot(data = data2_effs, mapping = aes(x = symmetry, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range symmetry", y = "Local co-occurrence", title = "Effect size")

ggplot(data = data2_prob, mapping = aes(x = symmetry, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Range symmetry", y = "Local co-occurrence", title = "Probability of co-occurrence")



# +++ cooccur ON morphdist -----
# +++++ naive cooccur -----
visreg( lm(eff_naive ~ poly(morphdist,2), data = data), xvar="morphdist", points=list(cex=1), xlab="Morphol. distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(eff_naive ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence", cex.lab=2, xtrans = log10 )
visreg( lm(prob_cooc_naive ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence", cex.lab=2, xtrans = log10 )

# +++++ union cooccur -----
visreg( lm(eff_union ~ poly(morphdist,2), data = data), xvar="morphdist", points=list(cex=1), xlab="Morphol. distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(eff_union ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence", cex.lab=2, xtrans = log10 )
visreg( lm(prob_cooc_union ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence", cex.lab=2, xtrans = log10 )

# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(morphdist,2), data = data), xvar="morphdist", points=list(cex=1), xlab="Morphol. distance", ylab="Co-occurrence in overlap", cex.lab=2 )
visreg( lm(eff_mask ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence in overlap", cex.lab=2, xtrans = log10 )
visreg( lm(prob_cooc_mask ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence in overlap", cex.lab=2, xtrans = log10 )

# +++++ bias in cooccur -----
# Nnaive = Nmask
visreg( lm( I(eff_mask-eff_naive) ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence bias", cex.lab=2, xtrans = log10, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_naive) ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence bias", cex.lab=2, xtrans = log10, main = "Probability of co-occurrence" )

# Nunion = Nmask
visreg( lm( I(eff_mask-eff_union) ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence bias", cex.lab=2, xtrans = log10, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_union) ~ poly(log10(morphdist),2), data = data), xvar="morphdist", points=list(cex=1), xlab="Log Morphol. distance", ylab="Co-occurrence bias", cex.lab=2, xtrans = log10, main = "Probability of co-occurrence" )



# +++ cooccur ON habdist -----
# +++++ naive cooccur -----
visreg( lm(eff_naive ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Effect size" )
visreg( lm(prob_cooc_naive ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Prob. of co-occurrence" )

# +++++ union cooccur -----
visreg( lm(eff_union ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Effect size" )
visreg( lm(prob_cooc_union ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Prob. of co-occurrence" )

# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Effect size" )
visreg( lm(prob_cooc_mask ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence", cex.lab=2, main = "Prob. of co-occurrence" )

# +++++ both in one plot -----
# Nunion = Nmask
ggplot(data = data2_effs, mapping = aes(x = habdist, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Habitat distance", y = "Local co-occurrence", title = "Effect size")

ggplot(data = data2_prob, mapping = aes(x = habdist, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Habitat distance", y = "Prob. of co-occurrence", title = "Probability of co-occurrence")

# +++++ bias in cooccur -----
# Nnaive = Nmask
visreg( lm( I(eff_mask-eff_naive) ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_naive) ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )

# Nunion = Nmask
visreg( lm( I(eff_mask-eff_union) ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_union) ~ poly(habdist,2), data = data), xvar="habdist", points=list(cex=1), xlab="Habitat distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )



# +++ cooccur ON foragdist -----
# +++++ naive cooccur -----
visreg( lm(eff_naive ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_naive ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ union cooccur -----
visreg( lm(eff_union ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_union ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_mask ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ bias in cooccur -----
# Nnaive = Nmask
visreg( lm( I(eff_mask-eff_naive) ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_naive) ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )

# Nunion = Nmask
visreg( lm( I(eff_mask-eff_union) ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )
visreg( lm( I(prob_cooc_mask-prob_cooc_union) ~ poly(foragdist,2), data = data), xvar="foragdist", points=list(cex=1), xlab="Foraging distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )



# +++ cooccur ON foragdist - separately -----
par(mfrow=c(2,2))
hist(data$fooddist)  # sqrt transf.
hist(data$methoddist)
hist(data$substratedist)
hist(data$stratumdist)  # sqrt transf.
# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(fooddist,2), data = data), xvar="fooddist", points=list(cex=1), xlab="Food distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(eff_mask ~ fooddist, data = data), xvar="fooddist", points=list(cex=1), xlab="Food distance", ylab="Co-occurrence", cex.lab=2 )
# sqrt trans:
visreg( lm(eff_mask ~ sqrt(fooddist), data = data), xvar="fooddist", points=list(cex=1), xlab="Food distance", ylab="Co-occurrence", cex.lab=2, xtrans = sqrt )
# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(methoddist,2), data = data), xvar="methoddist", points=list(cex=1), xlab="Foraging method distance", ylab="Co-occurrence", cex.lab=2 )
# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(substratedist,2), data = data), xvar="substratedist", points=list(cex=1), xlab="Foraging substrate distance", ylab="Co-occurrence", cex.lab=2 )
# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(stratumdist,2), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(eff_mask ~ stratumdist, data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum distance", ylab="Co-occurrence", cex.lab=2 )
# sqrt trans:
visreg( lm(eff_mask ~ sqrt(stratumdist), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum distance", ylab="Co-occurrence", cex.lab=2, xtrans = sqrt )



# +++ cooccur ON climdist -----
# +++++ naive cooccur -----
visreg( lm(eff_naive ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_naive ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ union cooccur -----
visreg( lm(eff_union ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_union ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ mask cooccur -----
visreg( lm(eff_mask ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )
visreg( lm(prob_cooc_mask ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climate distance", ylab="Co-occurrence", cex.lab=2 )

# +++++ both in one plot -----
# Nnaive = Nmask
ggplot(data = data2_effs, mapping = aes(x = climdist, y = effects, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Climatic distance", y = "Local co-occurrence", title = "Effect size") +
  scale_x_log10()

ggplot(data = data2_prob, mapping = aes(x = climdist, y = prob_cooccur, color = analysis_type)) +
  geom_point() + 
  geom_smooth(method = lm, formula=y~poly(x,1)) +
  labs(x = "Climatic distance", y = "Local co-occurrence", title = "Probability of co-occurrence") +
  scale_x_log10()

# +++++ bias in cooccur -----
# Nnaive = Nmask
visreg( lm( I(eff_mask-eff_naive) ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )  # quadratic
visreg( lm( I(eff_mask-eff_naive) ~ climdist, data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )  # linear

visreg( lm( I(prob_cooc_mask-prob_cooc_naive) ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )

# Nunion = Nmask
visreg( lm( I(eff_mask-eff_union) ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )  # quadratic
visreg( lm( I(eff_mask-eff_union) ~ climdist, data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size" )  # linear

visreg( lm( I(eff_mask-eff_union) ~ log10(climdist), data = data), xvar="climdist", xtrans = log10, points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Effect size: mask - union" )  # linear

summary( lm( I(eff_mask-eff_union) ~ log10(climdist), data = data) )
# Climate is not such a strong driver of co-occurrence as might appear from naive(union) analysis.
# But not signif. when as interaction btw. climate and analysis_type, only if BIAS analyzed

visreg( lm( I(prob_cooc_mask-prob_cooc_union) ~ poly(climdist,2), data = data), xvar="climdist", points=list(cex=1), xlab="Climatic distance", ylab="Co-occurrence bias", cex.lab=2, main = "Probability of co-occurrence" )



##### simple ANALYSES -----
# +++ correlations -----
cor.test(d_age_mask$node_age, d_age_mask$overlap)  # r = 0.12, p = 0.214
cor.test(d_age_naive$node_age, d_age_naive$overlap)  # r = 0.21, p = 0.008
cor.test(d_age_union$node_age, d_age_union$overlap)  # r = 0.10, p = 0.094

cor.test(d_age_mask$node_age, d_age_mask$symmetry)  # r = 0.04, p = 0.645
cor.test(d_age_naive$node_age, d_age_naive$symmetry)  # r = 0.23, p = 0.003
cor.test(d_age_union$node_age, d_age_union$symmetry)  # r = 0.08, p = 0.177

cor.test(d_age_mask$symmetry, d_age_mask$overlap)  # r = -0.22, p = 0.018
cor.test(d_age_naive$symmetry, d_age_naive$overlap)  # r = 0.13, p = 0.090
cor.test(d_age_union$symmetry, d_age_union$overlap)  # r = 0.09, p = 0.161


# +++ preparing data -----
# all for both naive and mask => Nnaive > Nmask
data2 <- d_age_mask %>%
  select(-c(10,14)) %>%
  bind_rows( "normal" = d_age_naive, "range-constrained" = ., .id = "analysis_type" )

data2 <- d_age_mask %>%
  bind_rows( "union" = d_age_union, "intersection" = ., .id = "analysis_type" )


# +++ analysing data -----
# +++++ simple comparison of effects (Fig. 4 parts) -----
# effects (Veech)
ggplot(data2_effs, aes(x = analysis_type, y = effects)) + 
  geom_boxplot()
require(ggpubr)
p2 <- data %>%
  rename(`Range overlap` = eff_mask, `Range union` = eff_union) %>%
  ggpaired( data, cond1 = "Range overlap", cond2 = "Range union", line.color = "gray", line.size = 0.4, xlab = "Type of analysis", ylab = "Syntopy (effect size)", font.x = 20, font.y = 20 ) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20)) +
  annotate("text", x = 0.535, y = 0.14, label = "(b)", size = 5)
  
data2_effs %>%
  group_by(analysis_type) %>% 
  summarise(meanEff = mean(effects), sample = n())

summary(lm(effects ~ analysis_type, data2_effs))  # effects on average lower in union than in overlap = makes sense, see simulations

# SES
ggplot(data2_SES, aes(x = analysis_type, y = SES)) + 
  geom_boxplot()
require(ggpubr)
p2ses <- data %>%
  rename(`Range overlap` = SES_mask, `Range union` = SES_union) %>%
  ggpaired( data, cond1 = "Range overlap", cond2 = "Range union", line.color = "gray", line.size = 0.4, xlab = "Type of analysis", ylab = "Syntopy (SES)", font.x = 20, font.y = 20 ) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20)) +
  annotate("text", x = 0.535, y = 8.5, label = "(b)", size = 5)

data2_SES %>%
  group_by(analysis_type) %>% 
  summarise(meanEff = mean(SES), sample = n())

summary(lm(SES ~ analysis_type, data2_SES))  # effects on average lower in union than in overlap = makes sense, see simulations


# +++++++ comparison of syntopy effects (bias on interecept): -----
require(nlme)  # analyzy v LME:
summary(lme( I(eff_union-eff_mask) ~ 1, random = ~1|node_ID , dataClus))  # toto v RESULTS

require(phyr)  # analyzy v PHYR:
# effects (Veech):
m.inla1 <- pglmm( I(eff_union-eff_mask) ~ 1 + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( I(eff_union-eff_mask) ~ 1 + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)
sum(ranef(m.inla2)[1:3, 1])/sum(ranef(m.inla2)[, 1])

# SES:
m.inla1 <- pglmm( I(SES_union-SES_mask) ~ 1 + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( I(SES_union-SES_mask) ~ 1 + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)
sum(ranef(m.inla2)[1:3, 1])/sum(ranef(m.inla2)[, 1])

# +++++++ bias in syntopy effect vs sympatry: -----
quartz(width = 3.3, height = 4)
visreg(lm( I(eff_union-eff_mask) ~ overlap, data), xlab="Sympatry", ylab="Bias in syntopy")
summary(lm( I(eff_union-eff_mask) ~ scale(overlap), data))  # union effect biased more towards negative values with low overlap = makes sense, see simulations

# analyza v LME:
summary(lme( I(eff_union-eff_mask) ~ scale(overlap), random = ~1|node_ID , dataClus))

# analyza v PHYR:
# effects (Veech):
m.inla1 <- pglmm( I(eff_union-eff_mask) ~ scale(overlap) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( I(eff_union-eff_mask) ~ scale(overlap) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)
sum(ranef(m.inla2)[1:3, 1])/sum(ranef(m.inla2)[, 1])

# SES:
m.inla1 <- pglmm( I(SES_union-SES_mask) ~ scale(overlap) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( I(SES_union-SES_mask) ~ scale(overlap) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)
sum(ranef(m.inla2)[1:3, 1])/sum(ranef(m.inla2)[, 1])


p3 <- ggplot(data = data, aes(x = overlap, y = I(eff_union-eff_mask))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Sympatry (%)", y = "Syntopy bias") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20)) +
  annotate("text", x = 11, y = 0.07, label = "(c)", size = 5)

# the same figure with SES instead of effects:
p3ses <- ggplot(data = data, aes(x = overlap, y = I(SES_union-SES_mask))) +
  geom_point() +
#  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  labs(x = "Sympatry (%)", y = "Syntopy bias") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20)) +
  annotate("text", x = 11, y = 5, label = "(c)", size = 5)


# bias in effect vs symmetryt
visreg(lm( I(eff_union-eff_mask) ~ symmetry, data), xlab="Symmetry", ylab="Bias in syntopy")
summary(lm( I(eff_union-eff_mask) ~ symmetry, data)) 


# +++++ SYMPATRY on: habitat and resource use distances -----
# +++++++ random effects with phyr: -----
# jointly on foraging (using all four traits) AND habitat:
require(phyr)
m.inla1 <- pglmm( overlap ~ scale(foragdist) + scale(habdist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( overlap ~ scale(foragdist) + scale(habdist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# separately individual resource use traits:
m.inla1 <- pglmm( overlap ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( overlap ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# +++++++ figure k teto analyze: -----
require(visreg)
quartz(width = 7.5, height = 4.6)
par(mfrow=c(1,2), mar=c(5,5,2,0.2))
visreg( lm(overlap ~ sqrt(fooddist), data = data), xvar="fooddist", points=list(cex=1), xlab="Diet", ylab="", cex.lab=1.7, xtrans = sqrt )
title(ylab = "Sympatry (range overlap)", line=3.5, cex.lab = 1.7)
visreg( lm(overlap ~ sqrt(stratumdist), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum", ylab="", cex.lab=1.7, xtrans = sqrt )
abline(h = 0, lty = 2)


# +++++ SYNTOPY on: habitat and resource use distances (for MS) -----
# +++++++ correlations of distances: -----
require(GGally)
ggpairs(data[, c(18:25)], title="Distances in ecological traits (116 species pairs)", axisLabels="show", lower=list(continuous = wrap("cor", size=4)), diag=list(continuous="densityDiag"), upper=list(continuous=wrap("smooth_loess", colour="darkgrey", pch=21, cex=0.4)))

quartz(width = 7, height = 7)  # Fig. for Supplement
ggpairs(data[, c(19,20,22:25)], columnLabels = c("Habitat","Resource use","Diet","Method","Substrate","Stratum"), title="Distances in ecological traits (116 species pairs)", axisLabels="show", lower=list(continuous = wrap("cor", size=4)), diag=list(continuous="densityDiag"), upper=list(continuous=wrap("smooth", colour="darkgrey", pch=21, cex=0.4)))

# +++++++ mixed effects models: -----
require(nlme)
# jointly on foraging (using all four traits) AND habitat:
summary(lme( eff_mask ~ scale(foragdist) + scale(habdist), random = ~1|node_ID , dataClus))

# separately individual resource use traits:
summary(lme( eff_mask ~ scale(fooddist) + scale(methoddist) + scale(stratumdist) + scale(substratedist), random = ~1|node_ID , dataClus))
summary(lme( eff_mask ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist), random = ~1|node_ID , dataClus))  # sqrt trans where appropriate, toto v RESULTS

# +++++++ random effects with phyr: -----
# jointly on foraging (using all four traits) AND habitat:
require(phyr)
# effects (Veech):
m.inla1 <- pglmm( eff_mask ~ scale(foragdist) + scale(habdist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v Results
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( eff_mask ~ scale(foragdist) + scale(habdist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# SES:
m.inla1 <- pglmm( SES_mask ~ scale(foragdist) + scale(habdist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v Results
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( SES_mask ~ scale(foragdist) + scale(habdist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# separately individual resource use traits:
# effects (Veech):
m.inla1 <- pglmm( eff_mask ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v Results
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( eff_mask ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# SES:
m.inla1 <- pglmm( SES_mask ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v Results
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( SES_mask ~ scale(sqrt(fooddist)) + scale(methoddist) + scale(sqrt(stratumdist)) + scale(substratedist) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# +++++++ figure k teto analyze (Fig. 6): -----
require(visreg)
quartz(width = 7.9, height = 4.5)
par(mfrow=c(1,2), mar=c(5,5,2,0.5))
visreg( lm(eff_mask ~ sqrt(fooddist), data = data), xvar="fooddist", points=list(cex=1), xlab="Diet", ylab="", cex.lab=1.7, xtrans = sqrt )
abline(h = 0, lty = 2)
text("(a)", x = 0.83, y = 0.135)
title(ylab = "Syntopy (effect size)", line=3.5, cex.lab = 1.7)
visreg( lm(eff_mask ~ sqrt(stratumdist), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum", ylab="", cex.lab=1.7, xtrans = sqrt )
abline(h = 0, lty = 2)
text("(b)", x = 0.78, y = 0.135)



# +++++ SYNTOPY on: node_age, overlap and symmetry -----
m <- lm( effects ~ node_age * analysis_type, data2)
m <- lm( prob_cooccur ~ node_age * analysis_type, data2)
m <- lm( effects ~ overlap * analysis_type, data2)
m <- lm( prob_cooccur ~ overlap * analysis_type, data2)

# dohromady overlap and node_age (jejich korelace nizka):
# +++++++ EFFECTS: -----
# Nnaive > Nmask, Nunion > Nmask podle toho, z jakych objektu udelam data2 vyse:
m <- lm( effects ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(symmetry):analysis_type, data2)
# m <- lm( effects ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(overlap):scale(symmetry), data2)  # nonsignificant

# Nnaive = Nmask, Nunion = Nmask podle toho, z jakych objektu udelam data2_effs vyse:
m <- lm( effects ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(symmetry):analysis_type, data2_effs)  # signif. interakce vzdy zmizi

# Finalni model na Nunion = Nmask bez interakce:
m <- lm( effects ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type, data2_effs)

# ale je signif i interakce overlap x symmetry
m <- lm( effects ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(overlap):scale(symmetry), data2_effs)  # p = 0.026

# bud jen union nebo mask
m <- lm( eff_union ~ scale(node_age) + scale(overlap) + scale(symmetry) + scale(overlap):scale(symmetry), data)
m <- lm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + scale(overlap):scale(symmetry), data)  # p = 0.0515 
summary(m)


require(effects)
plot(predictorEffects(m))
plot( predictorEffect("symmetry", m), lines=list(multiline = TRUE), confint=list(style="auto") )
plot( predictorEffect("overlap", m), lines=list(multiline = TRUE), confint=list(style="auto") )
plot( predictorEffect("symmetry", m, residuals = TRUE), lines=list(multiline = FALSE), confint=list(style="auto"), lattice=list(layout=c(4, 1)) )
plot( predictorEffect("overlap", m, residuals = TRUE), lines=list(multiline = FALSE), confint=list(style="auto"), lattice=list(layout=c(4, 1)) )
eff <- predictorEffect("symmetry", m)

plot( predictorEffect("analysis_type", m), lines=list(multiline = TRUE), confint=list(style="auto") )
plot( predictorEffect("analysis_type", m), lines=list(multiline = FALSE), lattice=list(layout=c(5, 1)) )


# plotting effects in 3D:
require(plot3D)

# for EFF_MASK only (delano na lm, protoze s lme to nefungovalo....):
# predict values on regular xy grid
grid.lines = 50
x.pred <- seq(min(data$overlap), max(data$overlap), length.out = grid.lines)
y.pred <- seq(min(data$symmetry), max(data$symmetry), length.out = grid.lines)
xy <- expand.grid( overlap = x.pred, symmetry = y.pred)
xy$node_age <- mean(data$node_age)
z.pred <- matrix(predict(m, newdata = xy), nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(m)
# scatter plot with regression plane
quartz(width = 6, height = 6)  # TAKES AGES to draw on quartz device!
scatter3D(data$overlap, data$symmetry, data$eff_mask, pch = 18, cex = 2, 
          theta = 20, phi = 20, ticktype = "detailed",
          xlab = "Range sympatry", ylab = "Range symmetry", zlab = "Syntopy (effect size)",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "")

# for BOTH effects only:
# predict values on regular xy grid
grid.lines = 50
x.pred <- seq(min(data2_effs$overlap), max(data2_effs$overlap), length.out = grid.lines)
y.pred <- seq(min(data2_effs$symmetry), max(data2_effs$symmetry), length.out = grid.lines)
xy <- expand.grid( overlap = x.pred, symmetry = y.pred)
xy$node_age <- mean(data2_effs$node_age)
xy$analysis_type <- "eff_mask"
z.pred <- matrix(predict(m, newdata = xy), nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(m)
# scatter plot with regression plane
scatter3D(data2_effs$overlap, data2_effs$symmetry, data2_effs$effects, pch = 18, cex = 2, 
          theta = 20, phi = 20, ticktype = "detailed",
          xlab = "Range overlap", ylab = "Range symmetry", zlab = "Co-occurrence",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "Effect size")

require(rgl)
rgl::plot3d()
surface3d()

require(bamlss)
xy$eff <- predict(m, newdata = xy)
x <- xy[, c(1,2,5)]
bamlss::plot3d(x = x, residuals = r)
bamlss::plot3d(x = x, contour=T)
bamlss::plot3d(x = x, image=T)

# +++++++ EFFECTS using MIXED effects (for MS): -----
# I added cluster ID manually and thus must load the data here:
dataClus <- read_excel("results/data_10my_nodelabels.xlsx", sheet = 1)
dataClus <- dataClus %>%
  left_join(x = ., y = data[, c(2:4)], by = "sp_id_Marki") %>%
  relocate(sp1_Marki, sp2_Marki, .after = sp_id_Marki)

# ++++++++++ lme: subclade as random effect -----
table(dataClus$node_ID)  # species pairs from only 10 clusters (out of 26 original) remain!
require(nlme)
# for control using lm:
m <- lm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry), dataClus)
summary(m)
# using lme but with random effect of node_ID:
m <- lme( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry), random = ~1|node_ID , dataClus)
summary(m)

# ++++++++++ phyr: vcv matrix + subclade as random effects -----
require(phyr)
require(INLA)
# failed for spp-pairs "vcv" matrix: Failed to standardize the cov matrix
#m.inla <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp_id_Marki__), data=dataClus, family="gaussian", cov_ranef = list(sp_id_Marki = res_heights), bayes = TRUE )

# species-based vcv matrix (FUNGUJE!, ale neumim jen fylogeneticky efekt bez efektu druhu samotneho):
# BAYESIAN:
m.inla1 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = TRUE )
summary(m.inla1)

m.inla2 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = TRUE )
summary(m.inla2)

# ML (even faster than Bayesian):
# effects (Veech):
m.inla1 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# SES:
m.inla1 <- pglmm( SES_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)  # toto je v RESULTS
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( SES_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# ML without symmetry:
m.inla1 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + (1|sp1_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp1_Marki = tree_MelABC_sp1), bayes = FALSE )
summary(m.inla1)
sum(ranef(m.inla1)[1:3, 1])/sum(ranef(m.inla1)[, 1])  # % var of all random effs

m.inla2 <- pglmm( eff_mask ~ scale(node_age) + scale(overlap) + (1|sp2_Marki__) + (1|node_ID), data=dataClus, family="gaussian", cov_ranef = list(sp2_Marki = tree_MelABC_sp2), bayes = FALSE )
summary(m.inla2)

# ++++++++++ figure k teto analyze: -----
require(visreg)
quartz(width = 3.9, height = 4.6)  # for 1 panel
par(mar=c(5,5,2,0.2))
visreg( lm(eff_mask ~ overlap, data = data), xvar="overlap", points=list(cex=1), xlab="Sympatry", ylab="Syntopy (effect size)", cex.lab=1.7 )
abline(h = 0, lty = 2)

quartz(width = 7.9, height = 4.5)  # for 2 panels
par(mfrow=c(1,2), mar=c(5,5,2,0.5))
visreg( lm(eff_mask ~ overlap, data = data), xvar="overlap", points=list(cex=1), xlab="Sympatry", ylab="", cex.lab=1.7 )
abline(h = 0, lty = 2)
text("(a)", x = 97.5, y = 0.135)
title(ylab = "Syntopy (effect size)", line=3.5, cex.lab = 1.7)
visreg( lm(eff_mask ~ symmetry, data = data), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="", cex.lab=1.7 )
abline(h = 0, lty = 2)
text("(b)", x = 48, y = 0.135)



require(coxme)
# failed: Cholesky factorization failed
#m1 <- lmekin( eff_mask ~ scale(node_age) + scale(overlap) + scale(symmetry) + (1|sp_id_Marki), varlist = list( res_heights_max ), data = data )



# +++++++ PROBS: -----
# Nnaive > Nmask, Nunion > Nmask podle toho, z jakych objektu udelam data2 vyse:
m <- lm( prob_cooccur ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(overlap):analysis_type, data2)

m <- lm( log10(prob_cooccur+0.001) ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(overlap):analysis_type, data2)  # pri log10-zavisle interakce zmizi

# m <- lm( prob_cooccur ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(symmetry):analysis_type, data2)  # nonsignificant

# Nnaive = Nmask, Nunion = Nmask podle toho, z jakych objektu udelam data2_prob vyse:
m <- lm( log10(prob_cooccur) ~ scale(node_age) + scale(overlap) + scale(symmetry) + analysis_type + scale(overlap):analysis_type, data2_prob)  # signifikance interakce vzdy zmizi

summary(m)

plot( predictorEffect("overlap", m, residuals = TRUE), lines=list(multiline = FALSE), confint=list(style="auto"), lattice=list(layout=c(2, 1)) )



# +++++++ 4-panel fig for SYNTOPY (Fig. 5) -----
require(visreg)
quartz(width = 6.6, height = 7)
par(mfrow=c(2,2), mar=c(5,5,1,0.5))

# with effects (Veech)
visreg( lm(eff_mask ~ overlap, data = data), xvar="overlap", points=list(cex=1), xlab="Sympatry", ylab="", cex.lab=1.7, gg = F )
abline(h = 0, lty = 2)
text("(a)", x = 97.5, y = 0.135, cex = 1.2)
title(ylab = "Syntopy (effect size)", line=3.5, cex.lab = 1.7)
visreg( lm(eff_mask ~ symmetry, data = data), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="", cex.lab=1.7, gg = F )
abline(h = 0, lty = 2)
text("(b)", x = 48, y = 0.135, cex = 1.2)
visreg( lm(eff_mask ~ sqrt(stratumdist), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum", ylab="", cex.lab=1.7, xtrans = sqrt, gg = F )
abline(h = 0, lty = 2)
text("(c)", x = 0.78, y = 0.135, cex = 1.2)
title(ylab = "Syntopy (effect size)", line=3.5, cex.lab = 1.7)
visreg( lm(eff_mask ~ sqrt(fooddist), data = data), xvar="fooddist", points=list(cex=1), xlab="Diet", ylab="", cex.lab=1.7, xtrans = sqrt, gg = F )
abline(h = 0, lty = 2)
text("(d)", x = 0.83, y = 0.135, cex = 1.2)

# with SES
visreg( lm(SES_mask ~ overlap, data = data), xvar="overlap", points=list(cex=1), xlab="Sympatry", ylab="", cex.lab=1.7, gg = F )
abline(h = 0, lty = 2)
text("(a)", x = min(data$overlap)+3.5, y = max(data$SES_mask)-0.4, cex = 1.2)
title(ylab = "Syntopy (SES)", line=3.5, cex.lab = 1.7)
visreg( lm(SES_mask ~ symmetry, data = data), xvar="symmetry", points=list(cex=1), xlab="Range symmetry", ylab="", cex.lab=1.7, gg = F )
abline(h = 0, lty = 2)
text("(b)", x = min(data$symmetry)+2, y = max(data$SES_mask)-0.4, cex = 1.2)
visreg( lm(SES_mask ~ sqrt(stratumdist), data = data), xvar="stratumdist", points=list(cex=1), xlab="Foraging stratum", ylab="", cex.lab=1.7, xtrans = sqrt, gg = F )
abline(h = 0, lty = 2)
text("(c)", x = min(data$stratumdist)+0.025, y = max(data$SES_mask)-0.4, cex = 1.2)
title(ylab = "Syntopy (SES)", line=3.5, cex.lab = 1.7)
visreg( lm(SES_mask ~ sqrt(fooddist), data = data), xvar="fooddist", points=list(cex=1), xlab="Diet", ylab="", cex.lab=1.7, xtrans = sqrt, gg = F )
abline(h = 0, lty = 2)
text("(d)", x = min(data$fooddist)+0.04, y = max(data$SES_mask)-0.4, cex = 1.2)



## MISCELLANEOUS -----
# +++ Z-scores -----
# see paper GEB_21_Carmona_cooccurrence for dark species pools.pdf, who cite Ulrich & Gotelli 2013 for Z-scores
# they define SES = (Nobs - Nexp) / SDhypergeometric, where SDhypergeometric is sqrt(variance_hypergeometric)
# variance_hypergeometric = (ni nj / N) ((N - ni) / N) ((N - nj) / (N - 1))
d_age_mask <- d_age_mask %>%
  mutate(Var = ((sp1_inc*sp2_inc)/no_sites)*((no_sites-sp1_inc)/no_sites)*((no_sites-sp2_inc)/(no_sites-1)), SES = (obs_cooccur - exp_cooccur) / sqrt(Var) )

# do finalniho objektu data pocitam SES nahore pri vytvareni objektu data (kolem radku 590)

data$SESmask <- d_age_mask$SES

p1 <- d_age_mask %>%
  ggplot(aes(x = SES)) +
  geom_histogram() +
  geom_vline(xintercept = 0, color = "darkred")

p2 <- d_age_mask %>%
  ggplot(aes(x = effects, y = SES)) +
  geom_point() +
  geom_smooth(method = "lm")

quartz(width = 8, height = 4)
require(patchwork)
p1 | p2

# +++ conceptual figure (Fig. 3) -----
require(car)
quartz(width = 6.8, height = 4.7)
par(mfrow=c(2,3), mar=c(3,3,2,2), oma=c(2,2,3,0))
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n", main = "", cex.main = 1.5)
lines(x = c(0,1), y = c(0.2,0.8), lwd = 3, col="darkgrey")  # puvodne vsude steelblue3 (tan3, violetred3, slategray3)
lines(x = c(0.5,1), y = c(0.5,0.8), lwd = 8, col="darkgrey")
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n", main = "b1) Habitat divergence", cex.main = 1.5, font.main = 1)
lines(x = c(0,1), y = c(0.8,0.2), lwd = 3, col="darkgrey")
lines(x = c(0.5,1), y = c(0.5,0.2), lwd = 8, col="darkgrey")
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n", main = "b2) Resource divergence", cex.main = 1.5, font.main = 1)
lines(x = c(0,1), y = c(0.2,0.8), lwd = 3, col="darkgrey")
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n")
lines(x = c(0,1), y = c(0.8,0.2), lwd = 3, col="darkgrey")  # puvodne vsude tan3
lines(x = c(0,0.5), y = c(0.8,0.5), lwd = 8, col="darkgrey")
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n")
lines(x = c(0,1), y = c(0.8,0.2), lwd = 3, col="darkgrey")
lines(x = c(0.5,1), y = c(0.5,0.2), lwd = 8, col="darkgrey")
plot(c(0,1), c(0,1), type = "n", xlab="", ylab="", xaxt="n", yaxt="n")
lines(x = c(0,1), y = c(0.2,0.8), lwd = 3, col="darkgrey")

mtext("Ecological divergence of species pairs",side=1,line=-0.5,outer=TRUE,cex=1.5)
mtext("Sympatry",side=1,line=-15.7,outer=TRUE,cex=1.5)
mtext("Syntopy",side=2,line=-0.5,outer=TRUE,cex=1.5,las=0)
mtext("a) Niche conservatism", side=1, outer=TRUE, font = 2, at = 0.17, line=-32)
mtext("b) Ecological isolation", side=1, outer=TRUE, font = 2, at = 0.67, line=-32)


# +++ plot local SR Meliphagoidea -----
df <- data.frame( SR_Mel = colSums(pa_mat_Mel_red_t) )
hist(df$SR_Mel)

quartz(width = 6, height = 6)
ggplot(pa, aes(x = LONG, y = LAT)) + 
  geom_point(aes(colour = df$SR_Mel), alpha = 0.6, size = 1.5) + 
  scale_colour_gradient(high = "#132B43", low = "#56B1F7") + 
  theme(legend.position = "top", axis.text = element_text(size = 12), axis.title = element_text(size = 20)) + 
  labs(colour = "Local species richness", x = "Longitude", y = "Latitude")

require(ggmap)  # totez v ggplot
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="white")  # create a layer of borders
mp <- ggplot() + mapWorld
mp <- mp + geom_point(aes(x=pa$LONG, y=pa$LAT, colour=df$SR_Mel), alpha=0.6, size=1.5)
mp <- mp + xlim(110,155) + ylim(-45,-9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("")
mp <- mp + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "top", panel.border = element_blank())
mp <- mp + labs(colour="Species richness of Meliphagoidea") + scale_colour_gradient(high = "#132B43", low = "#56B1F7")
quartz(width = 6, height = 6)
mp


# +++ FUNCTION: order matrix alphabetically -----
orderMatrixAlphab <- function(mat) {
  abec <- mat[order(rownames(mat)), ]
  abec <- t(abec)
  abec <- abec[order(rownames(abec)), ]
  abec <- t(abec)
  return(abec)
}


# +++ FUNCTION: add species ID column for merging / joining -----
addSpID <- function(target.df, names.df, guide.name = c("HANZAB","Pizzey","ABC","Marki","BLI6"), added.col = c("HANZAB","Pizzey","ABC","Marki","BLI6")) {
  # Adds ID from genus and species that I want to add based on names already present
  # For example, df might have HANZAB_names and I want to add Marki_names and unique ID from Marki names
  # target.df: df with data and some species names present - these are used to match new names
  # names.df: df with the translation table for all names that can be used, including names to add
  # guide.name: genus and species name in target.df that are used to match new names
  # added.col: names to add to target.df from names.df
  
  require(tidyverse)
  guide.name.full <- str_c(guide.name, "name", sep = "_")
  guide.name.sp1 <- str_c("sp1", guide.name, sep = "_")
  guide.name.sp2 <- str_c("sp2", guide.name, sep = "_")
  added.name.full <- str_c(added.col, "name", sep = "_")
  added.name.sp1 <- str_c("sp1", added.col, sep = "_")
  added.name.sp2 <- str_c("sp2", added.col, sep = "_")
  sp.id.name <- str_c("sp_id", added.col, sep = "_")
  name.to.uniteX <- str_c(added.name.full, ".x")
  name.to.uniteY <- str_c(added.name.full, ".y")
  
  new.df <- target.df %>%
    left_join(., select(names.df, all_of(c(guide.name.full, added.name.full))), by = setNames(guide.name.full, guide.name.sp1)) %>%
    left_join(., select(names.df, all_of(c(guide.name.full, added.name.full))), by = setNames(guide.name.full, guide.name.sp2)) %>%
    unite(col = {{ sp.id.name }}, all_of(name.to.uniteX), all_of(name.to.uniteY), sep = "_", remove = FALSE) %>%
    rename( setNames(name.to.uniteX, added.name.sp1), setNames(name.to.uniteY, added.name.sp2) )
  return(new.df)
}



