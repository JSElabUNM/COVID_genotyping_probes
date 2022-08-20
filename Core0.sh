for i in 30 32 38 40 45 48 59 64 
do
  echo "Looping ... number $i"

  for j in 0.5s 1s 4s
  do
    mkdir ${j}
    mkdir ${j}/${i}
    echo "Heatmaps/${j}/${i}_A*heatmap.csv | perl ReadHeatmap-Core0.pl > ${j}/$i/WholeGenome-int.json" 
    cat Heatmaps/${j}/${i}_A*heatmap.csv | perl ReadHeatmap-Core0.pl > ${j}/$i/WholeGenome-int.json
  done

done
