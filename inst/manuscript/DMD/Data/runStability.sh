nodes=(0 1 2 3 5 6 8 9 10 11 12 14)
for i in `seq 1 10`
do 
cat > runTest.sh << EOF
Rscript /data/dcosorioh/knkStability/testStability.R $i
EOF
qsub -l hostname=compute-1-${nodes[$i]} -e /data/dcosorioh/e${i}.log -o /data/dcosorioh/o${i}.log runTest.sh
done
