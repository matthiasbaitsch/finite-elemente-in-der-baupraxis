rm -f *.zip
for i in *-app.nb; do
    bn=$(basename $i .nb)
    zip $bn.zip $i
done