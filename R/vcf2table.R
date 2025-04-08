vcf2table <- function(vcf_in,
                      tumor=F,
                      tumor_in=NA,

                      # annotation ref gene
                      refGene=F,
                      refGene_func_name=NA,
                      refGene_gene_name=NA,
                      refGene_Exonicfunc_name=NA,
                      refGene_AAchange_name=NA,
                      refGene_GeneDetail_name=NA,

                      # annotation cosmic
                      cosmic=F,
                      cosmic_name=NA,

                      # annotation genomad
                      gnomad=F,
                      gnomad_name=NA,

                      # annotation dpGAP
                      dpGAP=F,
                      dpGAP_name=NA) {

  vcf <- read.vcfR(vcf_in,verbose = F)
  out_tmp <- data.frame(SampleID=colnames(vcf@gt)[2],
                        Chr=vcf@fix[,"CHROM"],
                        Start=as.integer(vcf@fix[,"POS"]),
                        End=as.integer(vcf@fix[,"POS"])+nchar(vcf@fix[,"REF"])-1,
                        Ref=vcf@fix[,"REF"],
                        Alt=vcf@fix[,"ALT"],
                        TLOD=extract.info(vcf,element="TLOD",as.numeric = T),
                        SOR=extract.info(vcf,element="SOR",as.numeric = T),
                        AD_alt=as.numeric(sapply(strsplit(extract.gt(vcf,element="AD"),","),function(x) x[2])),
                        AF=extract.gt(vcf,element="AF",as.numeric = T)[,1],
                        DP=extract.gt(vcf,element="DP",as.numeric = T)[,1],
                        SAF=as.numeric(sapply(strsplit(extract.gt(vcf,element="SB"),","),function(x) x[3])),
                        SAR=as.numeric(sapply(strsplit(extract.gt(vcf,element="SB"),","),function(x) x[4])))

  if (tumor==T) {
    tumor_vcf <- read.vcfR(tumor_in,verbose = F)
    tumor_df <- data.frame(SampleID_t=colnames(tumor_vcf@gt)[2],
                           Chr=tumor_vcf@fix[,"CHROM"],
                           Start=as.integer(tumor_vcf@fix[,"POS"]),
                           Ref=tumor_vcf@fix[,"REF"],
                           Alt=tumor_vcf@fix[,"ALT"],
                           tumor_AF=extract.gt(tumor_vcf,element="AF",as.numeric = T)[,1])

    out_tmp <- merge(out_tmp,tumor_df,by=c("Chr","Start","Ref","Alt"),all.x=T)
  }

  # check refGene databse
  if (refGene==T) {
    if (!is.na(refGene_func_name)) {
      out_tmp$Func.refGene <- extract.info(vcf,element=refGene_func_name)
    }

    if (!is.na(refGene_gene_name)) {
      out_tmp$Gene.refGene <- extract.info(vcf,element=refGene_gene_name)
    }

    if (!is.na(refGene_GeneDetail_name)) {
      out_tmp$GeneDetail.refGene <- extract.info(vcf,element=refGene_GeneDetail_name)
    }

    if (!is.na(refGene_Exonicfunc_name)) {
      out_tmp$ExonicFunc.refGene <- extract.info(vcf,element=refGene_Exonicfunc_name)
    }

    if (!is.na(refGene_AAchange_name)) {
      out_tmp$AAChange.refGene <- extract.info(vcf,element=refGene_AAchange_name)
    }
  }


  # check cosmic databse
  if (gnomad==T) {
    out_tmp$cosmic70 <- extract.info(vcf,element=cosmic_name)
  }
  # check genomad
  if (gnomad==T) {
    out_tmp$non_cancer_AF_popmax <- extract.info(vcf,element=gnomad_name,as.numeric = T)
  }
  # check dpGAP
  if (gnomad==T) {
    out_tmp$Alt_dpGAP_PopFreq <- extract.info(vcf,element=dpGAP_name,as.numeric = T)
  }


  return(out_tmp)

}
