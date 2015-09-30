#!/bin/bash

#get protein sequences for strata BLAST
#generally prefer NCBI or Ensembl, but will use VectorBase and other species-specific databases where necessary
#after getting sequences, clean up deflines to be easy to parse

#Drosophila
wget -O droMel.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dmel_r6.02/fasta/dmel-all-translation-r6.02.fasta.gz
wget -O droYak.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dyak_r1.3/fasta/dyak-all-translation-r1.3.fasta.gz
wget -O droMoj.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dmoj_r1.3/fasta/dmoj-all-translation-r1.3.fasta.gz
wget -O droVir.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dvir_r1.2/fasta/dvir-all-translation-r1.2.fasta.gz
wget -O droPse.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dpse_r3.2/fasta/dpse-all-translation-r3.2.fasta.gz
wget -O droAna.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dana_r1.3/fasta/dana-all-translation-r1.3.fasta.gz
wget -O droWil.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dwil_r1.3/fasta/dwil-all-translation-r1.3.fasta.gz
wget -O droGri.fa.gz ftp://ftp.flybase.net/releases/FB2014_05/dgri_r1.3/fasta/dgri-all-translation-r1.3.fasta.gz

#VectorBase
wget -O gloMor.fa.gz https://www.vectorbase.org/download/glossina-morsitans-yalepeptidesgmory14fagz
wget -O anoGam.fa.gz https://www.vectorbase.org/download/anopheles-gambiae-pestpeptidesagamp42fagz
wget -O aedAeg.fa.gz https://www.vectorbase.org/download/aedes-aegypti-liverpoolpeptidesaaegl32fagz
wget -O culQui.fa.gz https://www.vectorbase.org/download/culex-quinquefasciatus-johannesburgpeptidescpipj22fagz
wget -O anoDar.fa.gz https://www.vectorbase.org/download/anopheles-darlingi-coaripeptidesadarc31fagz
wget -O anoSte.fa.gz https://www.vectorbase.org/download/anopheles-stephensi-indianpeptidesastei22fagz
wget -O pedHum.fa.gz https://www.vectorbase.org/download/pediculus-humanus-usdapeptidesphumu21fagz
wget -O rhoPro.fa.gz https://www.vectorbase.org/download/rhodnius-prolixus-cdcpeptidesrproc12fagz
wget -O ixoSca.fa.gz https://www.vectorbase.org/download/ixodes-scapularis-wikelpeptidesiscaw13fagz

#NCBI
wget -O cerCap.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Ceratitis_capitata/protein/protein.fa.gz
wget -O nasVit.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Nasonia_vitripennis/protein/protein.fa.gz
wget -O apiMel.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Apis_mellifera/protein/protein.fa.gz
wget -O apiFlo.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Apis_florea/protein/protein.fa.gz
wget -O apiDor.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Apis_dorsata/protein/protein.fa.gz
wget -O bomImp.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Bombus_impatiens/protein/protein.fa.gz
wget -O bomTer.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Bombus_terrestris/protein/protein.fa.gz
wget -O megRot.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Megachile_rotundata/protein/protein.fa.gz
wget -O bomMor.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Bombyx_mori/protein/protein.fa.gz
wget -O triCas.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Tribolium_castaneum/protein/protein.fa.gz
wget -O metOcc.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Metaseiulus_occidentalis/protein/protein.fa.gz
wget -O hydMag.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Hydra_magnipapillata/protein/protein.fa.gz
wget -O acyPis.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/Acyrthosiphon_pisum/protein/protein.fa.gz

#ENSEMBL
wget -O megSca.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/megaselia_scalaris/pep/Megaselia_scalaris.Msca1.22.pep.all.fa.gz
wget -O helMel.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/heliconius_melpomene/pep/Heliconius_melpomene.Hmel1.22.pep.all.fa.gz
wget -O danPle.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/danaus_plexippus/pep/Danaus_plexippus.DanPle_1.0.22.pep.all.fa.gz
wget -O nemVec.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/nematostella_vectensis/pep/Nematostella_vectensis.GCA_000209225.1.22.pep.all.fa.gz
wget -O triAdh.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/trichoplax_adhaerens/pep/Trichoplax_adhaerens.ASM15027v1.22.pep.all.fa.gz
wget -O ampQue.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/amphimedon_queenslandica/pep/Amphimedon_queenslandica.Aqu1.22.pep.all.fa.gz
wget -O sacCer.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-22/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.22.pep.all.fa.gz
wget -O schPom.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-22/fasta/schizosaccharomyces_pombe/pep/Schizosaccharomyces_pombe.ASM294v2.22.pep.all.fa.gz
wget -O aspNid.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-22/fasta/aspergillus_nidulans/pep/Aspergillus_nidulans.ASM1142v1.22.pep.all.fa.gz
wget -O neuCra.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-22/fasta/neurospora_crassa/pep/Neurospora_crassa.ASM18292v1.22.pep.all.fa.gz
wget -O denPon.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/dendroctonus_ponderosae/pep/Dendroctonus_ponderosae.GCA_000355655.1.22.pep.all.fa.gz
wget -O tetUrt.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/tetranychus_urticae/pep/Tetranychus_urticae.GCA_000239435.1.22.pep.all.fa.gz
wget -O dapPul.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/daphnia_pulex/pep/Daphnia_pulex.GCA_000187875.1.22.pep.all.fa.gz
wget -O strMar.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/strigamia_maritima/pep/Strigamia_maritima.Smar1.22.pep.all.fa.gz
wget -O caeEle.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.22.pep.all.fa.gz
wget -O caeBre.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/caenorhabditis_brenneri/pep/Caenorhabditis_brenneri.C_brenneri-6.0.1b.22.pep.all.fa.gz
wget -O caeBri.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/caenorhabditis_briggsae/pep/Caenorhabditis_briggsae.CB4.22.pep.all.fa.gz
wget -O caeJap.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/caenorhabditis_japonica/pep/Caenorhabditis_japonica.C_japonica-7.0.1.22.pep.all.fa.gz
wget -O caeRem.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/caenorhabditis_remanei/pep/Caenorhabditis_remanei.C_remanei-15.0.1.22.pep.all.fa.gz
wget -O bruMal.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/brugia_malayi/pep/Brugia_malayi.B_malayi-3.0.22.pep.all.fa.gz
wget -O loaLoa.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/loa_loa/pep/Loa_loa.Loa_loa_V3.22.pep.all.fa.gz
wget -O oncVol.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/onchocerca_volvulus/pep/Onchocerca_volvulus.Cameroon_v3.22.pep.all.fa.gz
wget -O priPac.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/pristionchus_pacificus/pep/Pristionchus_pacificus.P_pacificus-5.0.22.pep.all.fa.gz
wget -O triSpi.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/trichinella_spiralis/pep/Trichinella_spiralis.Tspiralis1.22.pep.all.fa.gz
wget -O capTel.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/capitella_teleta/pep/Capitella_teleta.GCA_000328365.1.22.pep.all.fa.gz
wget -O helRob.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/helobdella_robusta/pep/Helobdella_robusta.GCA_000326865.1.22.pep.all.fa.gz
wget -O schMan.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/schistosoma_mansoni/pep/Schistosoma_mansoni.ASM23792v2.22.pep.all.fa.gz
wget -O craGig.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/crassostrea_gigas/pep/Crassostrea_gigas.GCA_000297895.1.22.pep.all.fa.gz
wget -O lotGig.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/lottia_gigantea/pep/Lottia_gigantea.GCA_000327385.1.22.pep.all.fa.gz
wget -O strPup.fa.gz ftp://ftp.ensemblgenomes.org/pub/release-22/metazoa/fasta/strongylocentrotus_purpuratus/pep/Strongylocentrotus_purpuratus.GCA_000002235.2.22.pep.all.fa.gz
wget -O cioInt.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/ciona_intestinalis/pep/Ciona_intestinalis.KH.75.pep.all.fa.gz
wget -O bosTau.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/bos_taurus/pep/Bos_taurus.UMD3.1.75.pep.all.fa.gz
wget -O canFam.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/canis_familiaris/pep/Canis_familiaris.CanFam3.1.75.pep.all.fa.gz
wget -O danRer.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/danio_rerio/pep/Danio_rerio.Zv9.75.pep.all.fa.gz
wget -O galGal.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/gallus_gallus/pep/Gallus_gallus.Galgal4.75.pep.all.fa.gz
wget -O homSap.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz
wget -O monDom.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/monodelphis_domestica/pep/Monodelphis_domestica.BROADO5.75.pep.all.fa.gz
wget -O musMus.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/pep/Mus_musculus.GRCm38.75.pep.all.fa.gz
wget -O oryLat.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/oryzias_latipes/pep/Oryzias_latipes.MEDAKA1.75.pep.all.fa.gz
wget -O xenTro.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.JGI_4.2.75.pep.all.fa.gz

#HGD
wget -O acrEch.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.acromyrmex/files/data/aech_OGSv3.8_pep.fa.gz
wget -O camFlo.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.camponotus/files/data/cflo_OGSv3.3_pep.fa.gz
wget -O harSal.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.harpegnathos/files/data/hsal_OGSv3.3_pep.fa.gz
wget -O linHum.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.linepithema/files/data/lhum_OGSv1.2_pep.fa.gz
wget -O pogBar.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.pogo/files/data/pbar_OGSv1.2_pep.fa.gz
wget -O solInv.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.solenopsis/files/data/sinv_OGSv2.2.3_pep.fa.gz
wget -O attCep.fa.gz http://hymenopteragenome.org/drupal/sites/hymenopteragenome.org.atta/files/data/acep_OGSv1.2_pep.fa.gz
