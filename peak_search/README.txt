mapping.pl
	各サンプルごとにtrimmomatic, bwaでmapping, samtoolsでsort&index
	trimmomatic > sample/trimed1.fastq or sample/trimed2.fastq
	bwa&satools sort > sample/sample.bam

peak_by_sliding_window.pl
	どれか一つのサンプルでも25bp windowの平均depthが100以上となっているpeakを抽出
	> sample/sample_allpeak.tsv.gz

peak_connect_count.pl
	上で作ったsample_allpeak.tsv.gzをまとめてlist化。複数sampleで共有しているpeakをまとめる際100bpのgapを許している。
	頭の$peak_depthの値を変えると平均depth $peak_depth以上のlistを吐き出すことが可能
	> 100_peak_list.tsv ($peak_depht = 100の場合)

extract_each_peak_depth.pl
	各peakでのsampleごとのdepthをsamtools bedcovを使ってtidyに書き出し。(最初はsamtools dpethでやっていたのでそのサブルーチンが最後に残っている)
	> all_peak_table_tidy.tsv

extract_each_peak_depth_of_plasmid.pl
	各peakにおいてreadのpaired-endやchimeric alignment先がGFPかmCherryなどのplasmidの配列のものを抽出し、それぞれreadのみでの平均depthを書き出し。
	> paired_plasmid_table_tidy.tsv (paired-endがplasmid sequenceのもの)
	> chimeric_align_table_tidy.tsv (chimeric alignmentがplasmid sequenceのもの)

search_peak_homologous_region.pl
	各sampleごとにsignificant(p_value <10**-5)なpeakをまとめてfastaを書き出し
	> fasta_dir_signif/sample_peak.fa
	これらの中でkmer=25で完全相同な領域がある場合にlist up
	> homologous_peak_lsit_signif/sample.tsv

