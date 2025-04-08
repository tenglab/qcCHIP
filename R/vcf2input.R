vcf2input <- function(vcf_path,
                      sample_list,
                      tumor=F,
                      tumor_path,
                      tumor_list,

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

  sample_list_full <- paste0(vcf_path,"/",sample_list)

  # no tumor
  vcf_list <- lapply(sample_list_full, function(file) {
    vcf2table(file,
              # annotation ref gene
              refGene=refGene,
              refGene_func_name=refGene_func_name,
              refGene_gene_name=refGene_gene_name,
              refGene_Exonicfunc_name=refGene_Exonicfunc_name,
              refGene_AAchange_name=refGene_AAchange_name,
              refGene_GeneDetail_name=refGene_GeneDetail_name,

              # annotation cosmic
              cosmic=cosmic,
              cosmic_name=cosmic_name,

              # annotation genomad
              gnomad=gnomad,
              gnomad_name=gnomad_name,

              # annotation dpGAP
              dpGAP=dpGAP,
              dpGAP_name=dpGAP_name)
  })

  # with tumor sample
  if (tumor==T) {
    tumor_full <- paste0(tumor_path,"/",tumor_list)
    paired_vectors <- Map(list, sample_list_full, tumor_full)
    vcf_list <- lapply(paired_vectors,function(pair)
      vcf2table(pair[[1]],
                tumor=tumor,
                tumor_in=pair[[2]],
                # annotation ref gene
                refGene=refGene,
                refGene_func_name=refGene_func_name,
                refGene_gene_name=refGene_gene_name,
                refGene_Exonicfunc_name=refGene_Exonicfunc_name,
                refGene_AAchange_name=refGene_AAchange_name,
                refGene_GeneDetail_name=refGene_GeneDetail_name,

                # annotation cosmic
                cosmic=cosmic,
                cosmic_name=cosmic_name,

                # annotation genomad
                gnomad=gnomad,
                gnomad_name=gnomad_name,

                # annotation dpGAP
                dpGAP=dpGAP,
                dpGAP_name=dpGAP_name))

  }

  out_final <- do.call(rbind,vcf_list)
  rownames(out_final) <- NULL

  # replace all NA to "."
  out_final[is.na(out_final)] <- "."

  return(out_final)
}
