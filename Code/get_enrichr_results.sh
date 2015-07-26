
for f in */*.nodes;
do
  echo "doing $f..."
  python ../../Code/queryEnrichr.py -g $f -d KEGG_2015 > KEGG_2015_table.txt;
  python ../../Code/queryEnrichr.py -g $f -d GO_Biological_Process > GO_biological_process.txt;
  python ../../Code/queryEnrichr.py -g $f -d MSigDB_Computational > MSigDB_computational.txt ;
  python ../../Code/queryEnrichr.py -g $f -d Drug_Perturbations_From_GEO > Drug_Perturbations_From_GEO.txt;
  python ../../Code/queryEnrichr.py -g $f -d CMAP_up > CMAP_up.txt;
  python ../../Code/queryEnrichr.py -g $f -d CMAP_down > CMAP_down.txt;
  echo "...done"
done;
