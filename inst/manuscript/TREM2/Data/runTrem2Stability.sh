nodes=(14 15 16 1 2 3 5 6 8 9 10 11)
for i in `seq 1 10`
do 
cat > runTest.sh << EOF
Rscript /data/dcosorioh/knkStability/trem2Stability.R $i
EOF
qsub -l hostname=compute-1-${nodes[$i]} -e /data/dcosorioh/trem2e${i}.log -o /data/dcosorioh/trem2o${i}.log runTest.sh
done
