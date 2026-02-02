#### Generate sample list for demultiplexing

donor=read.table("TOB_982_ind_list.txt",header=T)
pool=read.table("atac_pool_donor_match_list.txt",header=F,sep="\t")
id=read.table("TOB_phase1_976_full_ind_list.txt",header=T)
colnames(pool)=c("pool","experiment","id")
pool=pool[order(pool$pool),]

for(i in 1:nrow(pool)){

pool_name = pool$pool[i]
split_ids <- unlist(strsplit(pool$id[i], ", "))

if(all(split_ids %in% donor$TOB_ID2)){

onek1k_id = donor[match(split_ids,donor$TOB_ID2),"id"]
write.table(sort(onek1k_id),paste0(pool_name,"_sample_list.txt"),row.names=F,col.names=F,quote=F)
print(pool_name)

}else{

N_missing = 8-length(split_ids %in% donor$TOB_ID2)
print(paste0(pool_name, " is missing ", N_missing," donors. Double check the ids."))

}


}


####
