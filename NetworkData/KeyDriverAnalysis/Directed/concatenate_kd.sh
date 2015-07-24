rm -f keydriver_directed.xls # delete if exists
touch keydriver_directed.xls # create keydriver.xls

# Cat every file into keydriver.xls
for f in */*RASopathy_gene_keydriver.xls; do
    echo "$f" >>keydriver_directed.xls ;
    cat $f >> keydriver_directed.xls ;
done
