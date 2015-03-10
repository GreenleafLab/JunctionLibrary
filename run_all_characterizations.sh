# /usr/bin/env bash

python ./allKTjunctionsNB.06_dataframe.py
python ./threebythree.00_standard_dataframe.py
python ./threeBythree.01_helixcontext_dataframe.py
python ./threeBythree.02_alonghelix_dataframe.py
python ./threebythree.03_centralregion_dataframe.py
python ./threebythree.04_doubledouble_dataframe.py
python ./threeBythree.05_receptorLoop_dataframe.py

# remove headers
wd=/Users/Sarah/Dropbox/HJH_project/HJH_project/libraries/v1_characterization
for file in $wd/*characterization.txt;
do
    tail -n+2 $file > $(echo $file | sed 's/characterization.txt/.noheader.characterization/g')
done
head -n1 $file > header

# cat and add header
files=$(find $wd -name "*.noheader.characterization")
cat header $files > $wd/allJunctions.characterization