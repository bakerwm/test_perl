#! /bin/bash

## GRID-seq pipeline 
## version: 2016-11-3

genome="/home/wangming/work/gridseq/genome/"

# genome code for MACS2
gg="dm"

# genome version
gx="dm3"
gsize="$genome/${gx}.chrom.sizes"
gene="$genome/${gx}.ensGene.bed"

for k in S2; do 
## maping: RNA/DNA to genome
## filtering unimapped RNA mapq > 2
	for f in data/gridseq.${k}_rep*.fq.gz; do
		bam=${f/data/mapped}
		bam="${bam/.fq.gz/}.${gx}.bam";

		bowtie2 -p 32 --local -x $genome/idx_bowtie2/${gx} -U $f | \
		samtools view -Sbq 2 - > $bam
	done

## GRIDseq: unique readmates
	for i in 1 2; do 
		DNA="mapped/gridseq.${k}_rep${i}.gDNA.${gx}.bam"
		RNA="mapped/gridseq.${k}_rep${i}.cDNA.${gx}.bam"
	
		bedu="mapped/gridseq.${k}_rep${i}.rna_dna.bedu.gz"
		tmp=${bedu/bedu/tmp}
		
		bedtools bamtobed -i $RNA | \
		awk 'BEGIN{OFS="\t"} {$6 = $6=="+"? "-":"+"; print}' > $tmp
		bedtools bamtobed -i $DNA >> $tmp
	
		awk 'BEGIN { OFS = "\t" } {if(K[$4]) print K[$4], $0; else K[$4]=$0;}' $tmp | \
		awk 'length($1)<6 && length($7)<6{
			x = $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $12;
			D[x] = D[x]? D[x] "," $4 : $4 
		}
		END { OFS = "\t"; for(x in D) print x, D[x]} ' | gzip -c > $bedu
	
		rm $tmp
	done


## RNA Peakcalling
	bedu="mapped/gridseq.${k}_rep*.rna_dna.bedu.gz"
	bed1="mapped/gridseq.${k}.RNA+.bed"
	bed0="mapped/gridseq.${k}.RNA-.bed"
	pkrna="peaks/gridseq.${k}.RNA.peaks.bed"
	pkdna="peaks/gridseq.${k}.DNA.bin1k.bed"

	zcat $bedu | awk 'BEGIN{OFS="\t"} $4=="+"{print $1,$2,$3,$9,0,$4}' | sort -k1,1 -k2,2n > $bed1
	zcat $bedu | awk 'BEGIN{OFS="\t"} $4=="-"{print $1,$2,$3,$9,0,$4}' | sort -k1,1 -k2,2n > $bed0
	
	macs2 callpeak --verbose 0 --broad --nolambda --nomodel --extsize 100 -g $gg -n ${bed1/.bed/} -t $bed1
	awk 'BEGIN{OFS="\t"} $9>3{print $1,$2,$3,$1 ":" $2 "-" $3 ",+", $9, "+"}' ${bed1/.bed/_peaks.broadPeak} > ${bed1/RNA+/RNA.pktmp}

	macs2 callpeak --verbose 0 --broad --nolambda --nomodel --extsize 100 -g $gg -n ${bed0/.bed/} -t $bed0
	awk 'BEGIN{OFS="\t"} $9>3{print $1,$2,$3,$1 ":" $2 "-" $3 ",-", $9, "-"}' ${bed0/.bed/_peaks.broadPeak} >> ${bed1/RNA+/RNA.pktmp}

	sort -k1,1 -k2,2n ${bed1/RNA+/RNA.pktmp} > $pkrna
	rm ${bed1/RNA+/RNA.pktmp} $bed1 $bed0

## DNA bin=1k
	awk -v bin=1000 'BEGIN { OFS="\t" }
	!/chrM/{
		for(i=0;i<$2/bin-1;i++) {
			print $1, i*bin, (i+1)*bin, $1 ":" (i*bin) "-" ((i+1)*bin), 0, "+"
		}
	}' $gsize > $pkdna


## Gross coverage of gene tx
	ginfo="peaks/${gx}.ensGene_wo_chrM.bed"
	gcvg="peaks/gridseq.${k}.genecvg.bed"

if [[ ! -e $ginfo ]]; then
	awk -v g=$gsize ' BEGIN{OFS="\t"; while((getline < g)>0) G[$1] = $2}
	$1!="chrM" && G[$1]{print}' $gene | sort -k1,1 -k2,2n > $ginfo
fi

	zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3,$9,0,$4}' | \
	sort -k1,1 -k2,2n | bedtools coverage -sorted -g $gsize -counts -s -a $ginfo -b - | \
	awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$7,$6}' > $gcvg


## Active gene and Hit genes
	cf=$( Rscript gridRsub.getGeneHitCutoff.R $gcvg )

	pkgene="peaks/gridseq.${k}.RNA.peakgene.bed"
	pkcvg="peaks/gridseq.${k}.RNA.peakgene.cvg"


## Filtering the genes with peaks
	awk '$5>5 && !/chrM/' $pkrna | sort -k1,1 -k2,2n | \
	bedtools closest -s -d -a stdin -b $ginfo | \
	awk 'BEGIN{OFS="\t"} 
		$13!=0 && $3-$2>1000{x++; print $1,$2,$3,"Tx." x,$5,$6}
		$13==0{print $7,$8,$9,$10,$11,$12}' | \
	sort -k1,1 -k2,2n | uniq > $pkgene

## Calculating the gross coverage of gene
	zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3,$9,0,$4}' | 
	sort -k1,1 -k2,2n | bedtools coverage -sorted -g $gsize -counts -s -a $pkgene -b - | \
	awk -v cf=$cf 'BEGIN{OFS="\t"} $7>cf{print $1,$2,$3,$4,$7,$6}' > $pkcvg
	mv $pkcvg $pkgene
	
## Calculating the point coverage of gene
	zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
	sort -k1,1 -k2,2n | bedtools closest -s -d -a - -b $pkgene | \
	awk 'BEGIN{OFS="\t"} 
	$13==0{	split($4,X,"@"); 
		print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" $10 "@" $11 "@" $12 ,0,X[4]}' | \
	sort -k1,1 -k2,2n | bedtools closest -d -a - -b $pkgene | \
	awk '$13==0{rid=$4; did=$7 "@" $8 "@" $9 "@" $10 "@" $11 "@" $12; 
		X[rid "~" did]++; T[rid]++
		S[rid] += rid == did ? 1:0
	}
	END{ OFS="\t"; 
		for(i in X) {split(i, K, "~"); M[K[1]] = M[K[1]] > X[i] ? M[K[1]] : X[i] }
		for(k in T) {
			x = k; gsub("@","\t",x); 
			kt = T[k]
			ks = S[k] > 0 ? S[k] : 0
			km = M[k] > 0 ? M[k] : 0
			print x, kt, ks, km
		} 
	}' | sort -k1,1 -k2,2n > $pkcvg 


## GeneMids and GeneHits (gene)
	gmask="peaks/GeneMask.${k}.bed"
	gmids="peaks/GeneMids.${k}.bed"
	ghits="peaks/GeneHits.${k}.bed"

	bedtools intersect -sorted -s -wo -a $pkgene -b $pkgene | \
	awk '$4!=$10&&$3-$2==$13&&$4!~/protein_coding/' | \
	cut -f1-6 | sort -uk1,1 -k2,2n > $gmask

	awk '$4~/protein_coding/' $pkgene > $gmids

	awk -v fmsk=$gmask '
	BEGIN {OFS="\t"; trpk=100; mrpk=10; while((getline < fmsk)>0) G[$4]++}
	!G[$4] && ($7/($3-$2)*1000>=trpk || $9/($3-$2)*1000>=mrpk){print $1,$2,$3,$4,$5,$6}
	' $pkcvg > $ghits



## RNA-DNA network (gene)
	netg="matrix/gridseq.${k}.ghits.net.gz"
	netb="matrix/gridseq.${k}.gmids.net.gz"

	zcat $bedu | awk 'BEGIN{OFS="\t"} {print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
	sort -k1,1 -k2,2n | bedtools intersect -sorted -s -v -f 0.5 -wa -a - -b $gmask | \
	bedtools intersect -sorted -s -wa -wb -a - -b $ghits | \
	awk 'BEGIN{OFS="\t"} {split($4,X,"@"); split($10,Y,"|"); 
		print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" Y[1] "@" $11 "@" $12 ,$11,X[4]}' | \
	sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $pkdna | \
	awk '{rid=$4; did=$7 "@" $8; D[rid "~" did]++} 
	END{ OFS="\t"; 
		for(k in D) {
			split(k, X, "~"); 
			gsub("@", "\t", X[1]); gsub("@", "\t", X[2]);
			print X[1], X[2], D[k]
		}
	}' | sort -k1,1 -k2,2n -k7,7 -k8,8n | gzip -c > $netg 
	

	zcat $bedu | awk 'BEGIN{OFS="\t"} {print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
	sort -k1,1 -k2,2n | bedtools intersect -sorted -s -v -f 0.5 -wa -a - -b $gmask | \
	bedtools intersect -sorted -s -wa -wb -a - -b $gmids | \
	awk 'BEGIN{OFS="\t"} {split($4,X,"@"); split($10,Y,"|"); 
		print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" Y[1] "@" $11 "@" $12 ,$11,X[4]}' | \
	sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $pkdna | \
	awk '{rid=$4; did=$7 "@" $8; D[rid "~" did]++} 
	END{ OFS="\t"; 
		for(k in D) {
			split(k, X, "~"); 
			gsub("@", "\t", X[1]); gsub("@", "\t", X[2]);
			print X[1], X[2], D[k]
		}
	}' | sort -k1,1 -k2,2n -k7,7 -k8,8n | gzip -c > $netb


## Chromosome Backgrounds by mid genes
## Smooth bin (dm 10k, hs 100k)
	w=100

	bgdna="matrix/gridseq.${k}.DNA.bin1k.bg.bed"

	zcat $netb | awk -v pkdna=$pkdna -v g=$gsize -v w=$w '
	BEGIN{ OFS="\t"; 
		while((getline < g)>0) { C[$1] = $2 }; close(g)
	}
	$1 != $7{
		did = $7 ":" $8; D[did] += $9; N[$7] += $9
	}
	END{
		while((getline < pkdna)>0) {
			x = 0
			for(i=-w/2; i<w/2; i++) {
				j = ($2+i*1000)
				did = $1 ":" j
				v = D[did]/N[$1] * C[$1]/1e3
				x += v/w
			}
			print $1, $2, $3, $4, x, $6
		}
	}' > $bgdna


## RNA-DNA network normalize and peakfilter
	pknetg="matrix/gridseq.${k}.ghits.pkbin.net.gz"

	zcat $netg | awk -v g=$gsize -v fb=$bgdna -v win=10 -v nf=2 -v mwin=3 '
	BEGIN{ OFS="\t"; 
		while((getline < g)>0) {nc += $2; G[$1]=$2}; close(g)
		while((getline < fb)>0) B[$1 ":" $2] = $5; close(fb)
	}
	gid != $4 || $7 != chr || $8-pe>win*2000 {
		delete D
		if (NR>1) {
			for(p=ps; p<=pe; p+=1000) {
				nwl = 0; nwm = 0; nwr = 0
				for(i=-win; i<-win/2; i++) nwl += R[p+i*1000] > nf
				for(i=-win/2; i<win/2; i++) nwm += R[p+i*1000] > nf
				for(i=win/2; i<win; i++) nwr += R[p+i*1000] > nf

				R[p] = nwl > mwin || nwm > mwin || nwr > mwin ? R[p]:0
			}
			for(p=ps-win*500; p<pe+win*500; p+=1000) {
				vf = 0
				for(i=-win/2; i<win/2; i++) 
					vf += R[p+i*1000]/win
				D[p] = vf
			}
			for(p=ps-win*1000; p<pe+win*1000; p+=1000) {
				vf = 0
				for(i=-win/2; i<win/2; i++) 
					vf += D[p+i*1000]/win
				if(vf > 0 && p>=0 && p<G[chr]) 
					print ginfo, chr, p, vf
			}
		}
		ps = $8; delete R
	}
	{ 
		gid = $4; chr = $7; pe = $8+1000
		ginfo = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
		x = ($9/$5 * nc/1e3 + 1) / (B[$7 ":" $8] + 1)
		if(x > nf) R[$8] = x
	}
	END {
		delete D
		for(p=ps; p<=pe; p+=1000) {
			nwl = 0; nwm = 0; nwr = 0
			for(i=-win; i<-win/2; i++) nwl += R[p+i*1000] > nf
			for(i=-win/2; i<win/2; i++) nwm += R[p+i*1000] > nf
			for(i=win/2; i<win; i++) nwr += R[p+i*1000] > nf

			R[p] = nwl > mwin || nwm > mwin || nwr > mwin ? R[p]:0
		}
		for(p=ps-win*500; p<pe+win*500; p+=1000) {
			vf = 0
			for(i=-win/2; i<win/2; i++) 
				vf += R[p+i*1000]/win
			D[p] = vf
		}
		for(p=ps-win*1000; p<pe+win*1000; p+=1000) {
			vf = 0
			for(i=-win/2; i<win/2; i++) 
				vf += D[p+i*1000]/win
			if(vf > 0 && p>=0 && p<G[chr]) 
				print ginfo, chr, p, vf
		}
	}
	' | gzip -c > $pknetg

## peak regions
	pkinfo="matrix/gridseq.${k}.ghits.dnapeak.info"
	gflts="matrix/GeneHits.${k}.flt.bed"

	zcat $pknetg | awk -v gap=10 '
	gid != $4 || $7 != chr || $8-pe > gap*1000 {
		vp = chr ":" ps ":" pe
		if(NR > 1) D[gid] = D[gid] ? D[gid] "," vp : vp
		ps = $8	
	}
	{ 	gid = $4; chr = $7; pe = $8+1000; 
		G[gid] = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
	}
	END {	OFS = "\t"
		vp = chr ":" ps ":" pe
		D[gid] = D[gid] ? D[gid] "," vp : vp

		for(k in D) print G[k], D[k]
	}
	' | sort -k1,1 -k2,2n -k4,4 | \
	tee $pkinfo | cut -f1-6 > $gflts


### require >100G memory
### gene hits network --> matrix
#	mtx=${pknetg/net.gz/mtx.gz}
#
#	zcat $pknetg | sort -k7,7 | awk -v g=$gsize -v pkinfo=$pkinfo '
#	BEGIN{ OFS="\t"; 
#		while((getline < g)>0) C[$1] = $2; close(g)
#		while((getline < pkinfo)>0) { gids = gids ? gids "\t" $4 : $4 }; close(hit)
#		print "Chrom", gids
#		ylen = split(gids, G, "\t")
#	}
#	NR>1 && gchr != $7 {
#		xlen = int(C[gchr]/1000)
#
#		for(i=0; i<xlen; i++) {
#			vstr = gchr ":" (i*1000)
#			for(j=1; j<=ylen; j++) {
#				k = G[j] "|" i
#				x = ( V[k] ? int(V[k]*100)/100 : 0 )
#				vstr = vstr "\t" x
#			}
#			print vstr 
#		}
#
#		delete V
#	}
#	{ V[$4 "|" ($8/1000)] = $9; gchr = $7 }
#	END {
#		xlen = int(C[gchr]/1000)
#
#		for(i=0; i<xlen; i++) {
#			vstr = gchr ":" (i*1000)
#			for(j=1; j<=ylen; j++) {
#				k = G[j] "|" i
#				x = ( V[k] ? int(V[k]*100)/100 : 0 )
#				vstr = vstr "\t" x
#			}
#			print vstr 
#		}
#	}
#	' | gzip -c > $mtx
done 
