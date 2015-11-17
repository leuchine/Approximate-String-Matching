for i in 1 2 5 10 20
do
    ./psearch -k $i ../data/data-dblp.txt ../data/dblp-query-new.txt >topk_rlts_appgram_sample_$i.rlt
done
echo "topk query finish!"
