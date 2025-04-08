#' Convert individual VCF file to qcCHIP input table.
#'
#' A function to convert individual VCF file to table for qcCHIP
#'
#' @details
#' This function will convert individual VCF file to qcCHIP input file.
#'
#' @param vcf_path Path of VCF files of individual samples.
#' @param sample_list VCF file names of individual samples.
#' @param tumor If there are paired tumor samples.
#' @param tumor_path Path of VCF files of paired tumor samples.
#' @param tumor_list VCF file names of paired tumor samples.
#' @param refGene If there is annotation INFO of refGene.
#' @param refGene_func_name Name that contains the value function of gene (e.g. exonic, UTR3, ncRNA_intronic, etc.).
#' @param refGene_gene_name Name that contains the value name of gene.
#' @param refGene_Exonicfunc_name Name that contains the value exonic function of gene (e.g. stopgain, synonymous_SNV, nonsynonymous_SNV, etc.).
#' @param refGene_AAchange_name Name that contains the value AAchange of gene.
#' @param refGene_GeneDetail_name Name that contains the value Genedetail of gene.
#' @param cosmic If there is annotation from cosmic database.
#' @param cosmic_name Name that contains the cosmic value.
#' @param gnomad If there is annotation from gnomad database.
#' @param gnomad_name Name that contains the gnomad value.
#' @param dpGAP If there is annotation from dpGAP database.
#' @param dpGAP_name Name that contains the dpGAP value.
#'
#' @return
#' a table for qcCHIP
#'
#' @import vcfR
#'
#' @export
#' @examples
#'
#' # No paired tumor sample and no annotation of other database
#' out_1 <- vcf2input(vcf_path = vcf_path,
#'                    sample_list=sample_list)
#'
#' # With paired tumor sample and no annotation of other database
#' out_2 <- vcf2input(vcf_path = vcf_path,
#'                    sample_list=sample_list,
#'                    tumor=T,
#'                    tumor_path=tumor_path,
#'                    tumor_list=tumor_list)
#'
#' # Annotated VCF of sample and paired tumor sample. Annotated with refGene, cosmic, and gnomad.
#' out_3 <- vcf2input(vcf_path = vcf_path_annot,
#'                    sample_list=sample_list_annot,
#'                    tumor=T,
#'                    tumor_path=tumor_path_annot,
#'                    tumor_list=tumor_list_annot,
#'                    refGene=T,
#'                    refGene_func_name="Func.refGeneWithVer",
#'                    refGene_gene_name="Gene.refGeneWithVer",
#'                    refGene_Exonicfunc_name="ExonicFunc.refGeneWithVer",
#'                    refGene_AAchange_name=NA,
#'                    refGene_GeneDetail_name=NA,
#'                    cosmic=T,
#'                    cosmic_name="cosmic70",
#'                    gnomad=T,
#'                    gnomad_name="non_cancer_AF_popmax")

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
