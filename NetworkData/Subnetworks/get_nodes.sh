#rm -f keydriver_directed.xls # delete if exists
#touch keydriver_directed.xls # create keydriver.xls

# Cat every file into keydriver.xls
for f in ./*.sif; do
    python get_unique_nodes.py $f >> "$f".nodes ;
done
