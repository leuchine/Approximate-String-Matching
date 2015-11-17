for i in 1 2 5 10 20
do
    ./psearch -k $i ../data/author.data ../data/author-query >topk_rlts_appgram_aurthor_$i.rlt
done
echo "topk query finish!"
