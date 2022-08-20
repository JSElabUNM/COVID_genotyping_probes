for FILE in Multicore/scratch/*_012_R*
do 
  echo ${FILE}
  /usr/local/bin/Rscript replacement_B3.R ~/Dropbox/JeremyEdwards/COVID/Manuscript3/Data/$FILE
done






