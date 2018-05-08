set home="/home/boliva/COLABORATION/JBONET"

python $home/archdb/scripts/PDBeChem2SQL.py -d archdb/data/pdbechem -s archdb/sql/pdbechem/pdbechem.sql.gz -v
python $home/archdb/scripts/GO2SQL.py -d archdb/data/go -s archdb/sql/go/go.sql.gz -v
python $home/archdb/scripts/TaxID2SQL.py -d archdb/data/taxid -s archdb/sql/taxid/taxid.sql.gz -v
python $home/archdb/scripts/Uniprot2SQL.py -d archdb/data/uniprot -s archdb/sql/uniprot/uniprot._.sql.gz -v
python $home/archdb/scripts/Enzyme2SQL.py -d archdb/data/enzyme -s archdb/sql/enzyme/enzyme.sql.gz - v
python $home/archdb/scripts/DrugBank2SQL.py -d archdb/data/drugbank -s archdb/sql/drugbank/drugbank.sql.gz -v
python $home/archdb/scripts/PDB2SQL.py -d archdb/data/pdbsource -q archdb/data/pdbseq -s archdb/sql/pdb -v
python $home/archdb/scripts/SCOP2SQL.py -d archdb/data/scop -s archdb/sql/scop/scop.sql.gz -v
python $home/archdb/scripts/PDBTM2SQL.py -d archdb/data/pdbtm -s archdb/sql/pdbtm/pdbtm.sql.gz -v
python $home/archdb/scripts/CDhit2SQL.py -d archdb/data/archs/PDBseq.2_5.fa.0_4.clstr -s archdb/sql/cdhit/cdhit.sql.gz -v
python $home/archdb/scripts/Arch2SQL.py -d archdb/data/archs/ -s archdb/sql/archs/ -v

#GOTO 1

#LABEL 2
/Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14 -e "select CONCAT(c.pdb,'_',c.chain,'_',lc.position), ld.type, ld.ss1L from loop_description ld, loop2chain lc, chain c WHERE lc.chain=c.nid and ld.nid=lc.loop_id and lc.assignation='D'" > loops_by_order_ss1L.tsv

# ArachType Input= loops_by_order_ss1L.tsv
# Archtype output folders = archdb/data/classification/DS & archdb/data/classification/MCL

python $home/archdb/scripts/DS2SQL.py -d archdb/data/classification/DS/ -q archdb/data/archs/loops_by_order_ss1L.tsv -s archdb/sql/classification/DS -v
python $home/archdb/scripts/MCL2SQL.py -d archdb/data/classification/MCL/ -q archdb/data/archs/loops_by_order_ss1L.tsv -s archdb/sql/classification/MCL -v

#GOTO 3

#LABEL 1
# CREATE SCHEMA IF NOT EXISTS `ArchDB14` DEFAULT CHARACTER SET latin1 COLLATE latin1_bin ;
cat ArchDB14.create.00.dbsources.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
cat ArchDB14.create.01.pdbechem.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c pdbechem/pdbechem.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.02.go.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c go/go.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.03.taxid.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c taxid/taxid.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.04.uniprot.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
for x in uniprot/*sql.gz:
gunzip -c x | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.05.enzyme.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c enzyme/enzyme.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.06.drugBank.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c drugbank/drugbank.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.07.PDB.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
for x in pdb/*/*sql.gz:
gunzip -c x | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.08.SCOP.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c scop/scop.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.09.PDBTM.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
gunzip -c pdbtm/pdbtm.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

gunzip -c cdhit/cdhit.sql.gz | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

cat ArchDB14.create.10.archs.sql  | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
for x in archs/*/*sql.gz:
gunzip -c x | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

#GOTO 2

#LABEL 3
cat ArchDB14.create.11.loopclassification.sql | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14
for x in classificaton/*/*sql.gz:
gunzip -c x | /Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14

# Postprocessing
/Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14 -e " SELECT CONCAT(m.name,'.', cs.name) AS subclass_name, ld.loop_id AS loop_name, lc.seq_ali AS sequence FROM loop_description ld, loop2cluster lc, cluster_subclass cs, cluster_class cc, method m WHERE cs.nid=lc.cluster_nid AND cc.nid=cs.class_nid AND m.nid=cc.method and ld.nid = lc.loop_nid" > archdb_alignments.tsv
--> stamp
/Applications/MAMP/Library/bin/mysql --host=localhost -uroot -proot ArchDB14 -e "SELECT DISTINCT(CONCAT(c.pdb,'_',c.chain)) AS pdbID FROM chain c, loop2chain lc WHERE lc.chain=c.nid AND lc.assignation='D'" > loop_sources_list.txt
