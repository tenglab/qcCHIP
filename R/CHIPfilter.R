#' Extract potential CHIP from an annotated txt file with different filters.
#'
#' A function to extract CHIP with different filters
#'
#' @details
#' This function will return a txt file of CHIP candidates based on the filter settings.
#'
#' @param input a data frame of text annotated file with proper column names.
#' @param max_percent maximum percentage of sample size allowed for one variant.
#' @param max_percent_locus maximum percentage of sample size allowed for one locus.
#' @param DP_min minimum read depth.
#' @param SOR_min minimum SOR value.
#' @param Qual_min minimum Qual(TLOD) value.
#' @param Alt_AD_min minimum ALT AD value.
#' @param VAF_min minimum VAF value.
#' @param VAF_max maximun VAF value.
#' @param SAF_min minimum SAF value.
#' @param SAR_min minimum SAR value.
#' @param tumor_sample if there are paired tumor samples.
#' @param tumor_VAF_min minimum VAF value in paired tumor sample. (Available when tumor_sample is used)
#' @param tumor_muti factor used for VAF comparison between blood and tumor samples. (Available when tumor_sample is used)
#' @param germline_VAF_min minimum VAF value in blood sample when there is paired tumor sample. (Available when tumor_sample is used)
#' @param gnomad if gnomad reference are used.
#' @param gnomad_min minimum non cancer AF popmax from gnomad reference file. (Available when gnomad is used)
#' @param dpGAP if dpGAP reference are used.
#' @param dpGAP_min minimum ALT dpGAP PopFreq value in blood sample.
#' @param variant_effect_VAF_min minimum VAF value in blood sample for variant effect filter. (Available when dpGAP or gnomad reference are used)
#' @param not_nonsynonymous if including not nonsynonymous option.
#' @param blacklist_f a bed file contain blacklist regions to exclude.
#' @param info if print processing massage.
#'
#' @return
#' A text file with CHIP candidates.
#'
#' @import data.table
#' @import GenomicRanges
#'
#' @export
#' @examples
#' # input file
#' in_df <- fread("data/demo_input.txt")
#' # run default setting
#' out_1 <- CHIPfilter(in_df)
#'
#' # change different metrics
#' out_2 <- CHIPfilter(in_df,max_percent=0.02,DP_min = 40,VAF_min=0.002)
#'
#' # with paired tumor sample
#' out_3 <- CHIPfilter(in_df,tumor_sample = T,tumor_VAF_min = 0.02)
#'
#' # with gnomad or dpGAP reference file
#' out_4 <- CHIPfilter(in_df,gnomad = F,dpGAP = T)
#'
#' # with blacklist region
#' bl_f <- fread("blacklist.bed")
#' out_5 <- CHIPfilter(in_df,blacklist_f = bl_f)

CHIPfilter <- function(input,
                      # population parameter
                      max_percent=0.1,
                      max_percent_locus=0.3,
                      # normal sample filter parameter
                      DP_min=20,
                      SOR_min=3,
                      Qual_min=6.3,
                      Alt_AD_min=1,
                      VAF_min=0.02,
                      VAF_max=0.4,
                      SAF_min=0,
                      SAR_min=0,

                      # tumor
                      tumor_sample=F,
                      tumor_VAF_min=0.35,
                      tumor_muti=1,
                      germline_VAF_min=0.35,

                      # variant effect
                      gnomad=T,
                      gnomad_min=0.005,
                      dpGAP=F,
                      dpGAP_min=0.005,
                      variant_effect_VAF_min=0.25,
                      not_nonsynonymous=F,

                      # blacklist region
                      blacklist_f=F,
                      info=T) {

  # change column structure
  input$TLOD <- as.numeric(input$TLOD)
  input$SOR <- as.numeric(input$SOR)

  # add mutation name
  input$v_name <- paste(input$Chr,input$Start,input$Ref,input$Alt,sep=":")
  input$loci <- paste(input$Chr,input$Start,sep=":")

  ########
  # population
  ########
  if (info==T) {
    print("Perform population metrics")
  }

  # get sample size
  sampleSize <- length(unique(input$SampleID))


  # count number of samples have variant
  # variant freq
  tmp1 <- unique(input[,c("v_name","SampleID")])
  tmp2 <- unique(input[,c("loci","SampleID")])

  v_freq <- data.frame(table(tmp1$v_name))

  # loci freq
  v_loci_freq <- data.frame(table(tmp2$loci))

  #---------------
  # check filters exclude 1 or 2
  # 1. found in > 10% of normal samples
  v_name_f1 <- as.character(v_freq$Var1[which(v_freq$Freq > max_percent*sampleSize)])
  # 2. any variants in same locus in > 30% normal samples
  v_name_f3_tmp <- as.character(v_loci_freq$Var1[which(v_loci_freq$Freq > max_percent_locus*sampleSize)])
  v_name_f3 <- unique(input$v_name[which(input$loci %in% v_name_f3_tmp)])

  pop_exclude_v <- unique(c(v_name_f1,v_name_f3))

  pop_keep <- input[!(input$v_name %in% pop_exclude_v),]

  ########
  # normal sample filter (include if 1 AND 2 AND 3 AND 4 AND 5 AND 6 AND 7 AND8)
  ########
  if (info==T) {
    print("Perform technique metrics")
  }

  # add mutation-sample pair
  pop_keep$mut_sample <- paste0(pop_keep$v_name,"|",pop_keep$SampleID)

  normal_keep <- pop_keep[which(pop_keep$DP>DP_min &
                                  (pop_keep$SOR<SOR_min | is.na(pop_keep$SOR)) &
                                  (pop_keep$TLOD>Qual_min| is.na(pop_keep$TLOD)) &
                                  pop_keep$AD_alt >Alt_AD_min &
                                  pop_keep$AF>=VAF_min &
                                  pop_keep$AF<=VAF_max &
                                  pop_keep$SAF>SAF_min&
                                  pop_keep$SAR>SAR_min),]

  ########
  # tumor pair
  ########
  if (tumor_sample==T) {
    if (info==T) {
      print("Perform individual metrics with paired tumor sample")
    }

    normal_keep$tumor_AF <- as.numeric(normal_keep$tumor_AF)
    q1 <- which(normal_keep$AF > germline_VAF_min & normal_keep$tumor_AF > tumor_VAF_min)
    q2 <- which(normal_keep$AF > tumor_muti*normal_keep$tumor_AF)
    q3 <- which(normal_keep$tumor_AF==0)
    q4 <- which(normal_keep$AF > germline_VAF_min &
                  normal_keep$tumor_AF > tumor_VAF_min &
                  normal_keep$cosmic70!=".")

    # included if (not q1) AND (2 OR 3 OR 4)
    tmp_i <- unique(c(q2,q3,q4))
    index_keep <- tmp_i[! tmp_i %in% q1]

    sample_pair_keep <- normal_keep[index_keep,]
  } else {
    if (info==T) {
      print("No paired tumor sample, skip")
    }
    sample_pair_keep <- normal_keep
  }

  ########
  # variantEffect
  ########
  if (info==T) {
    print("Perform functional metrics")

  }
  # part 1: include
  # 1. include located in exonic region AND Not synonymous
  if (not_nonsynonymous==T) {
    if (info==T) {
      print("Perform not nonsunonymous metrics")

    }
    varianteffect_1 <- sample_pair_keep[which(sample_pair_keep$Func.refGene=="exonic" &
                                                sample_pair_keep$ExonicFunc.refGene!="synonymous SNV" &
                                                sample_pair_keep$ExonicFunc.refGene!="nonsynonymous SNV"),]
  } else {
    if (info==T) {
      print("No not nonsunonymous metrics find, skip")

    }
    varianteffect_1 <- sample_pair_keep[which(sample_pair_keep$Func.refGene=="exonic" &
                                                sample_pair_keep$ExonicFunc.refGene!="synonymous SNV"),]
  }


  # part 2: exclude (gnomAD non-cancer > 0.005 | dpGAP > 0.005) and AF>0.25
  if (gnomad==T & dpGAP==T) {
    if (info==T) {
      print("Perform gnomad and dpGAP metrics")
    }
    veffect_exclude <- varianteffect_1$mut_sample[which((as.numeric(varianteffect_1$non_cancer_AF_popmax) > gnomad_min |
                                                           as.numeric(varianteffect_1$Alt_dpGAP_PopFreq)> dpGAP_min) &
                                                          varianteffect_1$AF> variant_effect_VAF_min)]

  } else if (gnomad==T & dpGAP==F) {
    if (info==T) {
      print("Perform gnomad metrics only")
    }

    veffect_exclude <- varianteffect_1$mut_sample[which(as.numeric(varianteffect_1$non_cancer_AF_popmax) > gnomad_min &
                                                          varianteffect_1$AF> variant_effect_VAF_min)]

  } else if (gnomad==F & dpGAP==T) {
    if (info==T) {
      print("Perform dpGAP metrics only")
    }
    veffect_exclude <- varianteffect_1$mut_sample[which(as.numeric(varianteffect_1$Alt_dpGAP_PopFreq) > dpGAP_min &
                                                          varianteffect_1$AF> variant_effect_VAF_min)]
  } else {
    if (info==T) {
      print("No gnomad and dpGAP metrics find, skip")
    }
    veffect_exclude <- NA
  }

  varianteffect_keep <- varianteffect_1[!(varianteffect_1$mut_sample %in% veffect_exclude),]

  ########
  # remove blacklist regions
  ########
  if (is.null(nrow(blacklist_f))) {
    if (info==T) {
      print("No blacklist region bed file find, skip")

    }
    output_f <- varianteffect_keep
  } else {
    if (info==T) {
      print("Perform blacklist region excluding")

    }
    dust_gr <- GRanges(seqnames=blacklist_f[,1],IRanges(start=as.integer(blacklist_f[,2]),
                                                  end=as.integer(blacklist_f[,3])))

    in_gr <- GRanges(seqnames = varianteffect_keep$Chr,IRanges(start=varianteffect_keep$Start,
                                                               end=varianteffect_keep$End))


    # find overlap
    ol_tmp <- findOverlaps(in_gr,dust_gr)

    repet_exclude <- varianteffect_keep$mut_sample[unique(queryHits(ol_tmp))]

    output_f <- unique(varianteffect_keep[!(varianteffect_keep$mut_sample %in% repet_exclude),])

  }

  return(output_f)
}

