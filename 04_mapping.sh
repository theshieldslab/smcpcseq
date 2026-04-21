
#!/usr/bin/bash

scripts=${PWD}
cd /blue/robert.shields/arwalker/fastqc

mkdir images 

echo "<!DOCTYPE html>" > index.html
echo "<html>" >> index.html
echo "<head>" >> index.html
echo "<meta charset="utf-8">" >> index.html
echo "<title>UA159 Tn-Seq Analysis</title>" >> index.html
echo "</head>" >> index.html
echo "<body>" >> index.html


for i in $(cat ../samples.txt | cut -f 1 )
    do 
        unzip -p ../fastqc/${i}_fastqc.zip ${i}_fastqc/Images/per_base_quality.png > ./images/${i}.png

        echo -e "<p>"${i}"</p>" >> index.html
        echo -e "<figure class="half" style="display:flex">\"" >> index.html
        echo -e "<img style="width:600px" src= \"./images/${i}.png\" >" >> index.html
        echo -e "</figure>" >> index.html


done



echo "</body>" >> index.html

echo "</html>" >> index.html


cat index.html | sed 's\<figure class=half style=display:flex>"\<figure class=half style=display:flex>\g' > temp.html

mv temp.html samples_fastqc.html

rm index.html


cd ${scripts}