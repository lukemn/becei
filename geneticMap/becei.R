# becei genome/genetic map

require(ggplot2)
require(qtl)


gmap <- function(){
  
  getMarkerPos <- function(cross, mpos, nameSplitChar='.', lgs=1:6){
    
    gmap <- pull.map(cross, chr = lgs)
    mapn = unlist(lapply(gmap, names))
    mapdf <- data.frame(genetic=as.numeric(unlist(gmap)), 
                        ix = unlist(sapply(gmap, function(i) 1:length(i))),
                        marker = gsub('^[1-9].', '', mapn),
                        lg = tstrsplit(mapn, '.', fixed = T)[[1]],
                        stringsAsFactors = F
    )
    mapdf <- merge(mapdf, mpos[,c('marker', 'start', 'end')], sort=F)
    mapdf$scaf <- unlist(lapply(strsplit(mapdf$marker, nameSplitChar, fixed = T), '[[', 1))
    mapdf
  }
  
  # hifiasm (haploid mode)
  # 538 markers, 17 primary contigs, 3 dupes
  # all LGs concordant, BUT huge gap between contigs 3 (ptg000022l) and 13 (ptg000018l) on LG3
  # potential missassembly: HiC shows contacts between start of 18 and both start and, more weakly, ~0.5Mb end of 22
  # 22 end is 1Mb contig_36 (23) in m10 (and shows 2 small duplications around the 0.5Mb mark where HiC contacts change in hifiasm 22)
  # rest of 22 is scaffold_77 (4) in m10
  # also 8cM gap at very end of LG2 within Ctg016 (hifiasm 12), though ordered correctly. hifiasm 12 ~colinear with m10 contig_72 (some small deviations, but latter half is clean)
  # HiC 3 splits hifiasm 12 after 2.05Mb as repetitive, with the remainder (partially duplicated and inverted) in HiC 4 (which is otherwise == hifiasm 22)
  # 94.1Mb span
  load('~/Documents/github/becei/geneticMap/hifiasm/rqtl/gmap_hifi_q99_rimap_3.rda', verbose = T)
  pull.map(cross, chr=4)
  
  # TAKE AS BACKBONE
  # do manual stitching, then map hifi reads to each LG, take consensus with deepvariant > bcftools
  bb = reo = fread('~/Documents/github/becei/geneticMap/hifiasm/gmap_hifi_q99_orderedContigs.csv')
  dscaf = tstrsplit(names(dupes), '.', fixed=T)[[1]]
  # add dupes... no new sequences
  bbdupe = lapply(bb$scaf, function(x) {print(x); xx = unlist(dupes[dscaf==x]); if(!is.null(xx)) xscaf = unique(tstrsplit(xx, '.', fixed=T)[[1]]); xscaf[xscaf!=x]})
  bbdupe = cbind(bb$scaf[sapply(bbdupe, length)>0], sapply(bbdupe[sapply(bbdupe, length)>0], function(x) paste(unlist(x), collapse='_')))
  
  reo = reo[,c('scaf', 'fa', 'lg', 'strand')]
  reo$asm = 'hifiasm'
  reos = split(reo, reo$lg)
  # use m10 scaffold_8 to link hifiasm 54/6 
  reos[[1]] = rbind(reos[[1]], data.frame(scaf = 'Ctg014', fa = 'scaffold_8', lg=1, strand=NA, asm='m10'))
  
  # use m10 contig_68 to link hifiasm 37/9 (9 extends a bit further)
  reos[[2]] = rbind(reos[[2]], data.frame(scaf = 'Ctg002', fa = 'contig_68', lg=2, strand=NA, asm='m10'))
  
  # hic scaffolds 3 & 4 link hifiasm 22,12,18 (just). 
  # 12 sus. 
  # use hic scaffolds to thread m10 backbone (contig_2, contig_72, scaffold_77)
  # then take consensus
  reos[[3]] = rbind(data.frame(scaf = NA, fa = 'HiC_scaffold_3', lg=3, strand=NA, asm='hifiasm_hic'),
                    data.frame(scaf = NA, fa = 'HiC_scaffold_4', lg=3, strand=NA, asm='hifiasm_hic'),
                    data.frame(scaf = 'Ctg015', fa = 'contig_2', lg=3, strand='-', asm='m10'),
                    data.frame(scaf = 'Ctg018', fa = 'contig_72', lg=3, strand='-', asm='m10'),
                    data.frame(scaf = 'Ctg004', fa = 'scaffold_77', lg=3, strand='+', asm='m10'),
                    data.frame(scaf = 'Ctg023', fa = 'contig_36', lg=3, strand='+', asm='m10')
  )
                          
  
  # HiC_scaffold_11 links all. hifiasm 68 suspect (both HiC and m10 fragment it)
  # take consensus for both hifiasm hic 11 and m10 HiC 1, and compare mapping scores
  reos[[4]] = rbind(data.frame(scaf='NA', fa = 'HiC_scaffold_11', lg=4, strand=NA, asm='hifiasm_hic'),
                    data.frame(scaf='NA', fa = 'HiC_scaffold_1', lg=4, strand=NA, asm='m10_hic')
  )
  
  # hifiasm (and both HiCs) fragged. use overlapping? m10 Ctg003 (contig_166) + Ctg032 (contig190)
  reos[[5]] = rbind(data.frame(scaf='Ctg003', fa = 'contig_166', lg=5, strand=NA, asm='m10'),
                    data.frame(scaf='Ctg032', fa = 'contig_190', lg=5, strand=NA, asm='m10')
  )
  
  # HiC 9 links 43,3, HiC 10 links 3,62
  reos[[6]] = rbind(data.frame(scaf='NA', fa = 'HiC_scaffold_9', lg=6, strand=NA, asm='hifiasm_hic'),
                    data.frame(scaf='NA', fa = 'HiC_scaffold_10', lg=6, strand=NA, asm='hifiasm_hic')
  )
  do.call(rbind, reos)
  
  
  # flye m10
  load('~/Documents/github/becei/geneticMap/m10/gmap_hifi_q99_rimap_3.rda', verbose = T)
  load('~/Documents/github/becei/geneticMap/m10/gmap_hifi_q99_bestOrders_LL.rda', verbose = T)
  # 488 markers, 27 primary, 7 dupes
  # 93.29 Mb
  # drop 2 terminal singleton markers
  cross <- drop.markers(cross, c("Ctg036.1", "Ctg035.1"))
  
  # Ctg026.101 alone in LG7, two other contig markers are on LG5
  # but addition causes gaps
  tcross <- addmarker(cross, genotypes = cross$geno[[7]]$data, markername = "Ctg026.101", chr=5, pos=3.25)
  tcross <- tcross[1:6]
  tmap <- est.map(tcross, map.function = 'morgan')
  tcross <- replace.map(tcross, tmap)
  plotMap(cross)
  # 1 37:27:2 
  #     almost entirely Ctg002 (contig_68), clean split over hifiasm 9/37, meta contig_275/435
  # 2 17:14:7:16:20
  #     17 (contig_71) is merged with 14 (scaffold_8) in mapped_m9 contig_279 & hifiasm 6
  #     7 (contig_73) is merged with 16 (contig_70) in mapped_m9 scaffold_334 & hifiasm 30
  # 3 == Ctg003 (contig_166). clean split of hifiasm 24/38, mapped_m9 contig_572/1083
  # 4 15:18:4:23
  #     15 (contig_2 == hifiasm 18) is merged with 18 (contig_72 == hifiasm 12) in mapped_m9 contig_460
  #     4 (scaffold_77) is joined with 23 (contig_36) in hifiasm 22 (but gappy, see above)
  # 5 19:24:26:28:11:5:10:12
  #     19 (contig_3) == hifiasm 35 (redundant across all ass)
  #     24 (contig_35) contained by hifiasm 41 (with contig_37,236 : dupes in the right place)
  #     26 (contig_214) contained by hifiasm 68 (with several others), which is linked to contig_38 (11), part of hifiasm 28 (contains contig_210 [5], 84 [10], 9 [12])
  # 6 13:6:9 (13 ooo)
  #     13 (contig_213) == hifiasm 43 (redundant across all ass)
  #     6 (contig_193) contained by hifiasm 3, which links to 9 (contig_211, which contains hifiasm 62)
  
  # canu (all Q20 reads)
  # 482 markers, 82 primary, 19 dupes
  # 100.3 Mb span
  load('~/Documents/github/becei/geneticMap/canu/gmap_hifi_q99_rimap_3.rda', verbose = T)
  # adds ~no new joins. tig 195 links hifiasm 9 and 37

  # hic assembly (hifiasm / m10) adds a little, but cannot be trusted: 
  # frequent small scale rearrangements that are consistent across hifi assemblies, and greater discordance with genetic data
  # use links only
  # hifiasm: 
  #     HiC 1 links hifiasm 87,38 /2 linked by hifiasm 38
  #     HiC 5 links 9,37(frag), Hifiasm 37 links HiC 5,6
  #     HiC 7 links 6,54,30, Hifiasm 30 links HiC 7,8
  #     HiC 8 links 30,17 
  #     HiC 9 links 43,3, Hifiasm 3 links HiC 9,10
  #     HiC 10 links 3,62
  #     HiC 11 links 35,107,41,68(frag),43,28
  # m10:
  #     HiC 1 links contig_3,97,37,35,236,118*,120,128,214*,38(with terminal inversion introduced),210,84,9
  #     HiC 2 ~= contig_166 + HiC 32 et al
  #     HiC 3 links contig_68 and many others
  # others are more fragmented than input
  
  # m10 hic v hifiasm hic
  #     m10 1 ~= hifiasm 11
  #     m10 2 joins hifiasm 1,2
  #     m10 3 joins hifiasm 3,4*,5*,6,10*,9,8,7 
  #     but, 
  #       hifiasm 4 split over m10 3 and 46
  #       hifiasm 5 split over m10 3 and 12, and tiny 10,11
  #       hifiasm 8 split over m10 3 and 67
  #       hifiasm 7 split over 3 and 67 and also partially inverted wrt m10 3
  #       hifiasm 10 split over 3 and 60
  
  # m10 hic v m10 (both on map ctg)
  #     hic 1 merges many (major contig_3,37,35,236..214,38,210,84,9 = Ctg 19,22,24,25..26,11,5,10,12)
  unique(mapdf$scaf[mapdf$lg==5])
  # == m10 LG5
  #     hic 2 splits contig_166 (Ctg003) into several pieces. HiC is choppy, but genetic data looks fine
  # == m10 LG3
  unique(mapdf$scaf[mapdf$lg==3])
  rle(sign(diff(as.numeric(tstrsplit(names(pull.map(cross, chr=3)[[1]]), '.', fixed=T)[[2]]))))
  plotMap(cross, chr=3)
  x = subset(mm, lg==3)
  unique(x$scaf)
  ggplot(x, aes(cpos, genetic)) + geom_point()
  #   hic 3 merges many (contig_72,2,scaffold_77*,36,68*,211,193*,213,70,192,71,scaffold_8,73*)
  #     18,15,4,23,2,9,6,13,6,20,17,14,7
  #     but * = partial, and contig_68 inverted
  # == m10 LG2
  unique(tstrsplit(names(pull.map(cross, chr=2)[[1]]), '.', fixed=T)[[1]])
  names(pull.map(cross, chr=1)[[1]])
  names(pull.map(cross, chr=2)[[1]])
  plotMap(cross, chr=2)
  x = subset(mm, lg==2)
  unique(x$scaf)
  ggplot(x, aes(cpos, genetic)) + geom_point()
  
  
  # flye meta
  # 541 markers, 26 primary contigs, 6 dupes
  # concordant but for LG2 (both termini)
  # span 93.39 Mb
  load('~/Documents/github/becei/geneticMap/meta/gmap_hifi_q95_rimap_3.rda', verbose = T)
  
  # flye mapped reads to m10_120x genetic map contigs assembled at m9
  # 545 markers, 23 primary contigs, 7 dupes
  # concordant but for singleton contig marker on LG5 (Ctg033.1 within Ctg011)
  # 93.38 Mb span
  load('~/Documents/github/becei/geneticMap/mapped_m9/gmap_hifi_q95_rimap_3.rda', verbose = T)
  
  # flye mapped reads to m10_120x genetic map contigs (default settings)
  # 548 markers, 23 primary contigs, 5 dupes
  # concordant but for singleton contig marker on LG5 (Ctg034.1 within Ctg010)
  # 93.15 Mb span
  load('~/Documents/github/becei/geneticMap/mapped/gmap_hifi_q95_rimap_3.rda', verbose = T)
  
  ####
  # stitch1 pseudochromosomes (hifiasm backbone, m10 and hifiasm hic)
  # 6 LGs, 93.21Mb 
  load('~/Documents/github/becei/geneticMap/stitch1/gmap_hifi_q95_rimap_3.rda', verbose = T)
  # clean but for mangled terminal marker block in LG2 (pseudochrom LG4, a hifiasm hic scaffold)
  # good from 4737164-
  # 159  Ctg001.12851  37.984 73  2  4547985  4575899 Ctg001 73
  # 160  Ctg001.13451  38.516 74  2  4737164  4747942 Ctg001 74
  # 161   Ctg001.8701  39.048 75  2  3254000  3274847 Ctg001 75
  # 162   Ctg001.5651  39.575 76  2  2334130  2341158 Ctg001 76
  # 163   Ctg001.4051  41.217 77  2  1452699  1463685 Ctg001 77
  # 164      Ctg001.1  41.736 78  2     1194    67242 Ctg001 78
  # 165    Ctg001.401  41.736 79  2   545185   588052 Ctg001 79
  # 166  Ctg001.14201  49.451 80  2  4930221  4942465 Ctg001 80
  
  mapi = pull.map(cross, chr=2)[[1]]
  rle(sign(diff(as.numeric(tstrsplit(names(mapi), '.', fixed=T)[[2]]))))
  o = oo = 1:nmar(cross)[2]
  ### q95 no different. assuming gterr reasonable
  oo[8:9] = oo[9:8]
  oo[14:15] = oo[15:14]
  oo[42:43] = oo[43:42]
  oo[49:50] = oo[50:49]
  # oo[70:71] = oo[71:70]
  ### q99
  # 1 switch at same genetic pos
  oo[21:22] = oo[22:21]
  # 1 switch, discordant 0.19cM (LOD ~0)
  oo[60:61] = oo[61:60]
  # 1 switch, discordant 0.54cM (LOD -3)
  # oo[70:71] = oo[71:70]
  # invert 3 markers, LOD -6.5
  # oo[72:74] = oo[74:72]
  mapi[oo]
  compareorder(cross, 2, oo)
  cross <- switch.order(cross, 2, oo)
  # leaves first 69 markers concordant (up to 5.456173 Mb), terminal 10 markers not
  # hifiasm 28 runs from 5.41Mb > 
  # first 5 Mb is hifiasm 35,107,41,68 (fragmented in HiC, and 5 contigs in m10)
  # m10 contig 38 has small inversion that looks to be correct (4.96-5.19Mb in HiC 11)
  tail(subset(mapdf, lg==2), 12)
  
  ####
  # stitch2 pseudochromosomems (as above, hifiasm LG1 and 4 redone)
  # 6 LGs, 93.86Mb
  load('~/Documents/github/becei/geneticMap/stitch2/gmap_hifi_q99_rimap_3.rda', verbose = T)
  
  (cliks = sapply(1:6, function(i) attr(cross$geno[[i]]$map, 'loglik')))
  (bliks = sapply(1:6, function(i) attr(rcross$geno[[i]]$map, 'loglik')))
  for(i in 1:6) if(bliks[i]>cliks[i]) {print(i); cross$geno[[i]]$map = rcross$geno[[i]]$map; cross$geno[[i]]$data = rcross$geno[[i]]$data}
  
  summary(cross)
  # gmap <- est.map(cross, error.prob=gterr, map.function = 'morgan')
  # cross <- replace.map(cross, gmap)
  cross <- est.rf(cross)
  dim(table(mapdf$scaf, mapdf$lg))
  dim(table(c(mm$scaf, mm$marker_scaf), c(mm$lg, mm$lg)))
  sum(fai$len[fai$scaf %in% c(mm$scaf, mm$marker_scaf)])
  
  summaryMap(cross)
  plotMap(cross)
  
  # elegans homologs
  # reorder cross
  fai = fai[order(fai$fa),]
  fai$chrom = c(5, 2, 1, 6, 3, 4)
  fai = fai[order(fai$chrom),]
  fai$chromr = romanChroms()
  fai = fai[order(fai$scaf),]
  names(cross$geno) = fai$chromr
  cross$geno = cross$geno[order(fai$chrom)]
  
  # RF
  cross <- est.rf(cross)
  plotRF(cross)
  rf = melt(pull.rf(cross))
  names(rf)[1] = 'marker'
  rf = merge(rf, mpos[,c('contig', 'marker', 'start')])
  names(rf)[1] = 'marker.x'
  names(rf)[2] = 'marker'
  rf = merge(rf, mpos[,c('contig', 'marker', 'start')], by = 'marker')
  rf = rf[order(rf$contig.x),]
  rf$chrom.x = factor(rf$contig.x, labels = fai$chrom)
  rf$chrom.x = factor(rf$chrom.x, levels=1:6)
  rf = rf[order(rf$contig.y),]
  rf$chrom.y = factor(rf$contig.y, labels = fai$chrom)
  rf$chrom.y = factor(rf$chrom.y, levels=1:6)
  rf = subset(rf, !is.na(rf$value))
  p <- ggplot(rf, aes(start.x/1e6, start.y/1e6, col = value)) + geom_point(stroke=0, size=1, alpha=0.9, shape=15) + facet_grid(chrom.y~chrom.x, scales='free') + scale_color_viridis_c('Rec. Fraction', direction = -1) + mythbw + labs(x='Mb', y='Mb') +
    scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n=3)) + scale_x_reverse(expand=c(0,0), breaks = scales::pretty_breaks(n=3)) +
    theme(legend.position = 'top', strip.background = element_rect(color='NA'))
  ggsave('~/Documents/github/becei/geneticMap/geno_RF.pdf', h=6, w=6)
  
  # seg dist
  gt <- geno.table(cross); gt <- gt[order(gt$P),]
  gt$marker = row.names(gt)
  gt = merge(gt, mpos)
  ggplot(gt, aes(start/1e6, -log10(P.value))) + geom_point() + facet_grid(.~chr, scales='free')
  gt$aa = gt$AA/(gt$AA+gt$BB+gt$missing)
  gt$bb = gt$BB/(gt$AA+gt$BB+gt$missing)
  gt$miss = gt$missing/(gt$AA+gt$BB+gt$missing)
  names(gt)[2] = 'chrom'
  ggplot(gt, aes(start/1e6, aa)) + geom_point() + facet_grid(.~chrom, scales='free')
  gtprop = melt(as.data.table(gt[,c('aa', 'bb', 'miss', 'start', 'chrom')]), 4:5)
  gtprop$variable = factor(gtprop$variable, labels = c('QG2082', 'QG2083', 'NA'))
  p <- ggplot(gtprop, aes(start/1e6, y=value, fill = variable)) + geom_area(alpha=0.8, position=position_fill()) + facet_grid(.~chrom, scales ='free') + theme_classic() + labs(x='Physical distance (Mb)', y = 'Genotype frequency') + theme(legend.position = 'top', strip.background = element_rect(color=NA), axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n=3)) + coord_cartesian(ylim=c(0,1)) + geom_hline(aes(yintercept=0.5), col='white') + scale_x_continuous(expand=c(0.1,0), breaks = scales::pretty_breaks(n=3)) +
    scale_fill_brewer('', palette ='Dark2')
  ggsave('~/Documents/github/becei/geneticMap/geno_props.pdf', h=3, w=6)
  
  # marey
  cross <- flip.order(cross, chr=c('I', 'III', 'V'))
  mapdf = merge(getMarkerPos(cross, mpos, lgs = romanChroms()), fai)
  p <- ggplot(mapdf, aes(start/1e6, genetic)) + geom_point(size=1, alpha=.75, stroke=0) + facet_grid(.~chrom, scales = 'free', space='free') + theme_classic() + 
    labs(x='Physical distance (Mb)', y='Genetic distance (cM)') + scale_x_continuous(n.breaks = 3) +
    theme(legend.position = 'top', axis.line = element_blank(), strip.background = element_blank(), panel.border = element_rect(fill = NA), panel.spacing = unit(1, 'mm'), legend.margin = margin(0, 0, -5, 0))
  ggsave('~/Documents/github/becei/geneticMap/marey.pdf', h=4, w=8)
  
  p <- ggplot(mapdf, aes(start/1e6, genetic)) + geom_point(size=1, alpha=.75, stroke=0) + facet_grid(.~chrom, scales = 'free') + theme_classic() + 
    labs(x='Physical distance (Mb)', y='Genetic distance (cM)') + scale_x_continuous(n.breaks = 3) +
    theme(legend.position = 'top', axis.line = element_blank(), strip.background = element_blank(), panel.border = element_rect(fill = NA), panel.spacing = unit(1, 'mm'), legend.margin = margin(0, 0, -5, 0))
  ggsave('~/Documents/github/becei/geneticMap/marey_fixed.pdf', h=4, w=8)
  
  # cf selfers
  load('~/Documents/ctrop/selfer_maps.RData')
  mapn = do.call(rbind, lapply(split(mapdf, mapdf$chrom), function(x) cbind(x, gnorm = x$genetic/max(x$genetic), pnorm = x$start/max(x$start))))
  mapn = rbind(cbind(mapn[,c('pnorm', 'gnorm')], chrom=mapn$chromr, sp='becei'), gcf[,c('chrom', 'pnorm', 'gnorm', 'sp')])
  p <- ggplot(mapn, aes(pnorm, gnorm, col=sp)) + geom_point(size=1, alpha=.5, stroke=0) + geom_line() + facet_grid(.~chrom, scales = 'free', space='free') + theme_classic() + 
    labs(x='Normalized physical distance', y='Normalized genetic distance') + scale_x_continuous(expand=c(0.01,0.01), breaks = c(0,1)) + 
    scale_y_continuous(expand=c(0.01,0.01), breaks = c(0,1)) + scale_color_manual('', values = c('black', selfer_pal$col), labels = c('outcrosser', expression(italic(C.~briggsae)), expression(italic(C.~elegans)), expression(italic(C.~tropicalis)))) + 
    theme(legend.position = 'top', axis.text = element_blank(), axis.line = element_blank(), strip.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), panel.spacing = unit(1, 'mm'), legend.margin = margin(0, 0, -5, 0)) + 
    guides(color = guide_legend(override.aes = list(alpha=1, size=3, shape=19)))
  ggsave('~/Dropbox/Apps/Overleaf/lille/marey_fixed.png', h=4, w=8, dpi=300)
  
  
  save(cross, mapdf, dupes, file = '~/Documents/github/becei/geneticMap/becei_geneticMap.RData')
  
}

segmentDomains <- function(){
  
  load('~/Documents/github/becei/geneticMap/becei_geneticMap.RData', verbose=T)
  
  require(strucchange)
  require(segmented)
  
  o = mapdf[order(mapdf$chrom, mapdf$start),]
  # standardize to markers per chrom equally spaced in physical distance
  nMarkersPerChrom=200
  sampleMarkers <- function(x){
    ox = data.frame(chrom = x$chrom[1], start = as.integer(seq(1, max(x$start), length.out = nMarkersPerChrom)))
    ox$genetic = approx(x$start, x$genetic, xout = ox$start, rule = 2)$y
    ox
  }
  o = do.call(rbind, lapply(split(o, o$chrom), function(x) sampleMarkers(x)))
  
  doms =  do.call(rbind, lapply(1:6, function(i) {
    subo = subset(o, chrom==i)
    subo = subset(subo, genetic>0 & genetic<max(subo$genetic))
    bp = breakpoints(genetic~start, data=subo, breaks=2)
    oo = data.frame(apply(confint(bp, breaks=2, level=0.999)$conf+1, 2, function(x) subo$start[x]))
    names(oo) = c('bp.l', 'bp', 'bp.r')
    oo = cbind(chrom = i, j = c('l', 'r'), oo, method='struc')
    
    fit <- lm(genetic~start, subo)
    segfit <- segmented(fit, seg.Z = ~start, npsi=2)
    segfit = as.data.frame(confint.segmented(segfit, level = 0.999))
    names(segfit) = c('bp', 'bp.l', 'bp.r')
    oo = rbind(oo, cbind(chrom = i, j = c('l', 'r'), segfit, method='seg'))
  }))
  row.names(doms) = NULL
  
  # bps are ~identical except for X. 
  # take seg method
  o$chromr = factor(o$chrom, labels = romanChroms())
  doms$chromr = factor(doms$chrom, labels = romanChroms())
  ggplot(o, aes(start/1e6, genetic)) + geom_point(stroke=0, alpha=0.5) + theme_classic() + facet_grid(.~chromr, scales='free') +
    geom_vline(data = doms, aes(xintercept=bp/1e6, col = method)) + 
    geom_vline(data = doms, aes(xintercept=bp.l/1e6, col = method), linetype=2) + 
    geom_vline(data = doms, aes(xintercept=bp.r/1e6, col = method), linetype=2) + 
    theme(legend.position = 'top', strip.background = element_blank()) + 
    labs(x = 'Physical position (Mb)', y = "Genetic position (cM)") +
    scale_x_continuous(n.breaks = 3)
  ggsave('~/Documents/github/becei/geneticMap/marey_breakpoints.pdf', h=4, w=8)
  fwrite(doms, file = '~/Documents/github/becei/geneticMap/domains.csv')
  
  
  
  
  
  
}