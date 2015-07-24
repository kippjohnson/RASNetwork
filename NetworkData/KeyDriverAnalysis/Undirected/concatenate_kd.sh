rm -f keydriver_undirected.xls # delete if exists
touch keydriver_undirected.xls # create keydriver.xls

# Cat every file into keydriver.xls
for f in */*RASopathy_gene_keydriver.xls; do
    echo "$f" >>keydriver_undirected.xls ;
    cat $f >> keydriver_undirected.xls ;
done
