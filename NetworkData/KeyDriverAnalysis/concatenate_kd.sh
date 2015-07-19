rm -f keydriver.xls # delete if exists
touch keydriver.xls # create keydriver.xls

# Cat every file into keydriver.xls
for f in */*RASopathy_gene_keydriver.xls; do
    echo "$f" >>keydriver.xls ;
    cat $f >> keydriver.xls ;
done
