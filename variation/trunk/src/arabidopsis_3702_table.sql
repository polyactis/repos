set quoted_identifier  OFF 
GO

CREATE TABLE [AlleleFreqBySsPop] (
	[subsnp_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[allele_id] [int] NOT NULL ,
	[source] [char] (2) NOT NULL ,
	[cnt] [real] NOT NULL ,
	[freq] [real] NOT NULL ,
	[last_updated_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [Batch] (
	[batch_id] [int] NOT NULL ,
	[handle] [varchar] (20) NOT NULL ,
	[loc_batch_id] [varchar] (64) NOT NULL ,
	[loc_batch_id_upp] [varchar] (64) NOT NULL ,
	[batch_type] [char] (3) NOT NULL ,
	[status] [int] NOT NULL ,
	[simul_sts_status] [tinyint] NOT NULL ,
	[moltype] [varchar] (8) NOT NULL ,
	[method_id] [int] NOT NULL ,
	[samplesize] [int] NULL ,
	[synonym_type] [varchar] (255) NULL ,
	[submitted_time] [smalldatetime] NOT NULL ,
	[linkout_url] [varchar] (255) NULL ,
	[pop_id] [int] NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[success_rate_int] [int] NULL ,
	[build_id] [int] NULL ,
	[tax_id] [int] NOT NULL 
)
GO


CREATE TABLE [BatchCita] (
	[batch_id] [int] NOT NULL ,
	[position] [int] NOT NULL ,
	[pub_id] [int] NOT NULL ,
	[citation] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchCommLine] (
	[batch_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchCultivar] (
	[batch_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchMeExLine] (
	[batch_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchPrivLine] (
	[batch_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchStrain] (
	[batch_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [BatchValCode] (
	[batch_id] [int] NOT NULL ,
	[validation_status] [tinyint] NOT NULL 
)
GO


CREATE TABLE [Contact] (
	[batch_id] [int] NOT NULL ,
	[handle] [varchar] (20) NOT NULL ,
	[name] [varchar] (255) NOT NULL ,
	[fax] [varchar] (255) NULL ,
	[phone] [varchar] (255) NULL ,
	[email] [varchar] (255) NULL ,
	[lab] [varchar] (255) NULL ,
	[institution] [varchar] (255) NULL ,
	[address] [varchar] (255) NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [FreqSummaryBySsPop] (
	[subsnp_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[source] [varchar] (1) NOT NULL ,
	[chr_cnt] [int] NOT NULL ,
	[ind_cnt] [int] NOT NULL ,
	[non_founder_ind_cnt] [int] NOT NULL ,
	[chisq] [real] NULL ,
	[df] [tinyint] NULL ,
	[hwp] [real] NULL ,
	[het] [real] NULL ,
	[het_se] [real] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [GeneIdToName] (
	[locus_id] [int] NOT NULL ,
	[locus_symbol] [varchar] (64) NOT NULL ,
	[locus_name] [varchar] (255) NULL ,
	[tax_id] [int] NOT NULL ,
	[last_update_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [GtyFreqBySsPop] (
	[subsnp_id] [int] NULL ,
	[pop_id] [int] NULL ,
	[unigty_id] [int] NULL ,
	[source] [varchar] (1) NULL ,
	[cnt] [real] NULL ,
	[freq] [real] NULL ,
	[last_updated_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [HapSetSnpList] (
	[hapset_id] [int] NOT NULL ,
	[snp_order] [int] NOT NULL ,
	[subsnp_id] [int] NOT NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [IndGrpCode] (
	[code] [tinyint] NOT NULL ,
	[name] [varchar] (32) NOT NULL ,
	[descrip] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [IndivAlias] (
	[src_id] [int] NOT NULL ,
	[src_ind_id] [varchar] (64) NOT NULL ,
	[alias] [varchar] (64) NOT NULL ,
	[create_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [IndivBySource] (
	[ind_id] [int] NOT NULL ,
	[src_id] [int] NOT NULL ,
	[src_ind_id] [varchar] (64) NOT NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[src_ind_grp] [varchar] (32) NULL 
)
GO


CREATE TABLE [IndivSourceCode] (
	[code] [int] NOT NULL ,
	[name] [varchar] (22) NOT NULL ,
	[descrip] [varchar] (255) NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[src_type] [varchar] (10) NULL ,
	[display_order] [tinyint] NULL 
)
GO


CREATE TABLE [Individual] (
	[ind_id] [int] NOT NULL ,
	[descrip] [varchar] (255) NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[tax_id] [int] NULL ,
	[ind_grp] [tinyint] NULL 
)
GO


CREATE TABLE [Log_delBatch] (
	[issue_time] [smalldatetime] NOT NULL ,
	[issue_user] [varchar] (20) NOT NULL ,
	[cmd] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [OmimVarLocusIdSNP] (
	[omim_id] [int] NOT NULL ,
	[locus_id] [int] NOT NULL ,
	[omimvar_id] [char] (4) NOT NULL ,
	[locus_symbol] [char] (10) NOT NULL ,
	[var1] [char] (2) NOT NULL ,
	[aa_position] [int] NOT NULL ,
	[var2] [char] (2) NOT NULL ,
	[var_class] [int] NOT NULL ,
	[snp_id] [int] NOT NULL 
)
GO


CREATE TABLE [Pedigree] (
	[ped_id] [numeric](7, 0) NOT NULL ,
	[curator] [varchar] (12) NOT NULL ,
	[curator_ped_id] [varchar] (12) NOT NULL ,
	[create_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [PedigreeIndividual] (
	[ped_id] [decimal](7, 0) NOT NULL ,
	[ind_id] [int] NOT NULL ,
	[ma_ind_id] [int] NULL ,
	[pa_ind_id] [int] NULL ,
	[sex] [char] (1) NULL ,
	[create_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [PipelineIn] (
	[snp_type] [char] (2) NOT NULL ,
	[snp_id] [int] NOT NULL ,
	[input_time] [datetime] NOT NULL ,
	[inproc_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [PopLine] (
	[pop_id] [int] NOT NULL ,
	[line_num] [int] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [PopMandLine] (
	[pop_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [Population] (
	[pop_id] [int] NOT NULL ,
	[handle] [varchar] (20) NOT NULL ,
	[loc_pop_id] [varchar] (64) NOT NULL ,
	[loc_pop_id_upp] [varchar] (64) NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[src_id] [int] NULL 
)
GO


CREATE TABLE [RsMergeArch] (
	[rsHigh] [int] NOT NULL ,
	[rsLow] [int] NOT NULL ,
	[build_id] [smallint] NOT NULL ,
	[orien] [tinyint] NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[last_updated_time] [smalldatetime] NOT NULL ,
	[rsCurrent] [int] NULL ,
	[orien2Current] [tinyint] NULL 
)
GO


CREATE TABLE [SNP] (
	[snp_id] [int] NOT NULL ,
	[avg_heterozygosity] [real] NULL ,
	[het_se] [real] NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[CpG_code] [tinyint] NULL ,
	[tax_id] [int] NOT NULL ,
	[validation_status] [tinyint] NULL ,
	[exemplar_subsnp_id] [int] NOT NULL ,
	[univar_id] [int] NOT NULL ,
	[cnt_subsnp] [tinyint] NOT NULL ,
	[map_property] [tinyint] NULL 
)
GO


CREATE TABLE [SNP2Omim] (
	[omim_id] [int] NOT NULL ,
	[var_id] [char] (4) NOT NULL ,
	[locus_id] [int] NOT NULL ,
	[snp_id] [int] NOT NULL ,
	[aa1] [char] (1) NOT NULL ,
	[aa_pos] [int] NOT NULL ,
	[aa2] [char] (1) NOT NULL 
)
GO


CREATE TABLE [SNP3D] (
	[snp_id] [int] NOT NULL ,
	[protein_acc] [char] (10) NOT NULL ,
	[master_gi] [int] NOT NULL ,
	[neighbor_gi] [int] NOT NULL ,
	[aa_position] [int] NOT NULL ,
	[var_res] [char] (1) NOT NULL ,
	[contig_res] [char] (1) NOT NULL ,
	[neighbor_res] [char] (1) NOT NULL ,
	[neighbor_pos] [int] NOT NULL ,
	[var_color] [int] NOT NULL ,
	[var_label] [int] NOT NULL 
)
GO


CREATE TABLE [SNPAlleleFreq] (
	[snp_id] [int] NOT NULL ,
	[allele_id] [int] NOT NULL ,
	[chr_cnt] [real] NOT NULL ,
	[freq] [real] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [SNPAncestralAllele] (
	[snp_id] [int] NOT NULL ,
	[ancestral_allele_id] [int] NOT NULL 
)
GO


CREATE TABLE [SNPGtyFreq] (
	[snp_id] [int] NOT NULL ,
	[unigty_id] [int] NOT NULL ,
	[ind_cnt] [float] NULL ,
	[freq] [float] NULL ,
	[last_updated_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [SNPHWProb] (
	[snp_id] [int] NOT NULL ,
	[df] [tinyint] NULL ,
	[chisq] [real] NULL ,
	[hwp] [real] NULL ,
	[ind_cnt] [smallint] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [SNPHistory] (
	[snp_id] [int] NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL ,
	[history_create_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [SNPIndGtyFlag] (
	[snp_id] [int] NOT NULL ,
	[ind_id] [smallint] NOT NULL ,
	[gty_flag] [tinyint] NOT NULL ,
	[last_updated_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [SNPSubSNPLink] (
	[subsnp_id] [int] NOT NULL ,
	[snp_id] [int] NOT NULL ,
	[substrand_reversed_flag] [bit] NOT NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[build_id] [int] NULL ,
	[comment] [varchar] (255) NULL 
)
GO


CREATE TABLE [SNPSubSNPLinkHistory] (
	[subsnp_id] [int] NOT NULL ,
	[snp_id] [int] NOT NULL ,
	[build_id] [int] NOT NULL ,
	[history_create_time] [smalldatetime] NOT NULL ,
	[link_create_time] [datetime] NULL ,
	[link_last_updated_time] [datetime] NULL ,
	[orien] [tinyint] NULL ,
	[build_id_when_history_made] [int] NULL 
)
GO


CREATE TABLE [SNPVal] (
	[batch_id] [int] NOT NULL ,
	[snp_id] [int] NOT NULL 
)
GO


CREATE TABLE [SnpInSts] (
	[sts_uid] [int] NOT NULL ,
	[tax_id] [int] NULL ,
	[snp_id] [int] NOT NULL ,
	[create_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [SnpPidAcc] (
	[pid] [int] NOT NULL ,
	[accession] [varchar] (15) NULL ,
	[accession_ver] [int] NULL ,
	[tax_id] [int] NULL ,
	[create_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [SubInd] (
	[batch_id] [smallint] NOT NULL ,
	[subsnp_id] [int] NOT NULL ,
	[submitted_ind_id] [smallint] NOT NULL ,
	[submitted_strand_code] [int] NULL ,
	[allele_flag] [tinyint] NULL ,
	[gty_id] [int] NULL ,
	[submitted_rs] [int] NULL 
)
GO


CREATE TABLE [SubPop] (
	[batch_id] [int] NOT NULL ,
	[subsnp_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[type] [char] (3) NOT NULL ,
	[samplesize] [int] NOT NULL ,
	[submitted_strand_code] [tinyint] NULL ,
	[submitted_rs] [int] NULL ,
	[allele_flag] [tinyint] NULL ,
	[ambiguity_status] [tinyint] NULL ,
	[sub_heterozygosity] [real] NULL ,
	[est_heterozygosity] [real] NULL ,
	[est_het_se_sq] [real] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL ,
	[observed] [varchar] (255) NULL ,
	[sub_het_se_sq] [real] NULL ,
	[subpop_id] [int] NOT NULL 
)
GO


CREATE TABLE [SubPopAllele] (
	[batch_id] [int] NOT NULL ,
	[subsnp_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[allele] [char] (1) NOT NULL ,
	[other] [varchar] (255) NULL ,
	[freq] [real] NULL ,
	[cnt_int] [int] NULL ,
	[freq_min] [real] NULL ,
	[freq_max] [real] NULL ,
	[data_src] [varchar] (6) NULL ,
	[type] [char] (3) NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[allele_flag] [tinyint] NULL ,
	[cnt] [real] NULL ,
	[allele_id] [int] NULL ,
	[subpop_id] [int] NOT NULL 
)
GO


CREATE TABLE [SubPopGty] (
	[subpop_id] [int] NOT NULL ,
	[gty_id] [int] NOT NULL ,
	[gty_str] [varchar] (64) NOT NULL ,
	[cnt] [real] NULL ,
	[freq] [real] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [SubPopRsMerge] (
	[subpop_id] [int] NOT NULL ,
	[submitted_rs] [int] NOT NULL ,
	[current_rs] [int] NOT NULL ,
	[orien] [bit] NOT NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [SubPop_new] (
	[batch_id] [int] NOT NULL ,
	[subsnp_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[type] [char] (3) NOT NULL ,
	[samplesize] [int] NOT NULL ,
	[submitted_strand_code] [tinyint] NULL ,
	[submitted_rs] [int] NULL ,
	[allele_flag] [tinyint] NULL ,
	[ambiguity_status] [tinyint] NULL ,
	[sub_heterozygosity] [real] NULL ,
	[est_heterozygosity] [real] NULL ,
	[est_het_se_sq] [real] NULL ,
	[last_updated_time] [smalldatetime] NOT NULL ,
	[observed] [varchar] (255) NULL ,
	[sub_het_se_sq] [real] NULL ,
	[subpop_id] [int] IDENTITY (1, 1) NOT NULL 
)
GO


CREATE TABLE [SubSNP] (
	[subsnp_id] [int] NOT NULL ,
	[known_snp_handle] [varchar] (20) NULL ,
	[known_snp_loc_id] [varchar] (64) NULL ,
	[known_snp_loc_id_upp] [varchar] (64) NULL ,
	[batch_id] [int] NOT NULL ,
	[loc_snp_id] [varchar] (64) NULL ,
	[loc_snp_id_upp] [varchar] (64) NULL ,
	[synonym_names] [varchar] (255) NULL ,
	[loc_sts_id] [varchar] (64) NULL ,
	[loc_sts_id_upp] [varchar] (64) NULL ,
	[segregate] [char] (1) NOT NULL ,
	[indiv_homozygosity_detected] [char] (1) NULL ,
	[PCR_confirmed_ind] [char] (1) NULL ,
	[gene_name] [varchar] (64) NULL ,
	[sequence_len] [int] NULL ,
	[samplesize] [int] NULL ,
	[EXPRESSED_SEQUENCE_ind] [char] (1) NULL ,
	[SOMATIC_ind] [char] (1) NULL ,
	[sub_locus_id] [int] NULL ,
	[create_time] [smalldatetime] NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[ancestral_allele] [varchar] (255) NULL ,
	[CpG_code] [tinyint] NULL ,
	[variation_id] [int] NULL ,
	[top_or_bot_strand] [char] (1) NULL ,
	[validation_status] [tinyint] NULL ,
	[snp_id] [int] NULL ,
	[tax_id] [int] NOT NULL ,
	[chr_id] [tinyint] NULL 
)
GO


CREATE TABLE [SubSNPAcc] (
	[subsnp_id] [int] NOT NULL ,
	[acc_type_ind] [char] (1) NOT NULL ,
	[acc_part] [varchar] (16) NOT NULL ,
	[acc_ver] [int] NULL 
)
GO


CREATE TABLE [SubSNPCommLine] (
	[subsnp_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [SubSNPDeletedBySubmitter] (
	[subsnp_id] [int] NOT NULL ,
	[handle] [varchar] (20) NOT NULL ,
	[loc_snp_id_upp] [varchar] (64) NOT NULL ,
	[comment_id] [smallint] NOT NULL ,
	[reload_time] [datetime] NULL ,
	[create_time] [datetime] NOT NULL ,
	[delete_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [SubSNPMdFailLn] (
	[subsnp_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [SubSNPNoVariSeq] (
	[subsnp_id] [int] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [SubSNPSeq3] (
	[subsnp_id] [int] NOT NULL ,
	[type] [tinyint] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [SubSNPSeq5] (
	[subsnp_id] [int] NOT NULL ,
	[type] [tinyint] NOT NULL ,
	[line_num] [tinyint] NOT NULL ,
	[line] [varchar] (255) NOT NULL 
)
GO


CREATE TABLE [SubSNPSeqPos] (
	[subsnp_id] [int] NOT NULL ,
	[contig_acc] [varchar] (20) NOT NULL ,
	[contig_pos] [int] NOT NULL ,
	[chr] [varchar] (2) NULL ,
	[upstream_len] [int] NOT NULL ,
	[downstream_len] [int] NOT NULL ,
	[last_update_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [SubSNP_top_or_bot] (
	[subsnp_id] [int] NOT NULL ,
	[top_or_bot] [char] (1) NULL ,
	[step] [tinyint] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [SubmittedIndividual] (
	[submitted_ind_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[loc_ind_id_upp] [varchar] (64) NOT NULL ,
	[ind_id] [int] NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[last_updated_time] [smalldatetime] NULL ,
	[tax_id] [int] NOT NULL ,
	[loc_ind_alias] [varchar] (64) NULL ,
	[loc_ind_id] [varchar] (64) NULL ,
	[loc_ind_grp] [varchar] (32) NULL 
)
GO


CREATE TABLE [Synonym] (
	[subsnp_id] [int] NOT NULL ,
	[type] [varchar] (64) NOT NULL ,
	[name] [varchar] (64) NULL 
)
GO


CREATE TABLE [UnigeneSnp] (
	[gi] [int] NOT NULL ,
	[snp_id] [int] NOT NULL ,
	[unigene_cid] [int] NOT NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [dnSummaryByType] (
	[type] [varchar] (20) NOT NULL ,
	[id] [int] NOT NULL ,
	[build_id] [int] NOT NULL ,
	[ss_cnt] [int] NOT NULL ,
	[rs_cnt] [int] NOT NULL ,
	[create_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [dn_IND_batchCount] (
	[batch_id] [int] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[ss_cnt] [int] NOT NULL ,
	[rs_cnt] [int] NOT NULL ,
	[ind_cnt] [int] NOT NULL ,
	[create_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [dn_IND_batch_pop] (
	[batch_id] [smallint] NOT NULL ,
	[pop_id] [int] NOT NULL ,
	[update_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [dn_batchCount] (
	[batch_id] [int] NOT NULL ,
	[ss_cnt] [int] NOT NULL ,
	[rs_cnt] [int] NOT NULL ,
	[rs_validated_cnt] [int] NOT NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[pop_cnt] [int] NULL ,
	[ind_cnt] [int] NULL 
)
GO


CREATE TABLE [dn_chr_locusList] (
	[locus_id] [int] NOT NULL ,
	[contig_chr] [char] (2) NOT NULL 
)
GO


CREATE TABLE [dn_crossHandle] (
	[inHandle] [varchar] (20) NULL ,
	[inHandleSource] [char] (1) NULL ,
	[sbid] [int] NULL ,
	[pbid] [int] NULL ,
	[ibid] [int] NULL ,
	[common_ss_cnt] [int] NULL ,
	[common_rs_cnt] [int] NULL ,
	[create_time] [smalldatetime] NULL ,
	[sorder] [int] NULL ,
	[porder] [int] NULL ,
	[iorder] [int] NULL ,
	[last_updated_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [dn_gty_rsCnt_byChr] (
	[ind_id] [int] NULL ,
	[contig_chr] [char] (2) NOT NULL ,
	[rsCnt] [int] NULL 
)
GO


CREATE TABLE [dn_gty_rsCnt_byLocus] (
	[ind_id] [int] NULL ,
	[locus_id] [int] NOT NULL ,
	[top_level_class] [char] (5) NOT NULL ,
	[rsCnt] [int] NULL 
)
GO


CREATE TABLE [dn_handleCount] (
	[handle] [varchar] (20) NOT NULL ,
	[batch_type] [char] (3) NOT NULL ,
	[ss_cnt] [int] NOT NULL ,
	[rs_cnt] [int] NULL ,
	[rs_validated_cnt] [int] NULL ,
	[create_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [dn_snpFxnCnt] (
	[build_id] [int] NOT NULL ,
	[fxn_class] [tinyint] NULL ,
	[snp_cnt] [int] NOT NULL ,
	[gene_cnt] [int] NOT NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[last_updated_time] [smalldatetime] NOT NULL ,
	[tax_id] [int] NOT NULL 
)
GO


CREATE TABLE [dn_summary] (
	[tax_id] [int] NOT NULL ,
	[build_id] [int] NOT NULL ,
	[type] [varchar] (20) NOT NULL ,
	[cnt] [int] NOT NULL ,
	[create_time] [smalldatetime] NOT NULL ,
	[last_updated_time] [smalldatetime] NOT NULL 
)
GO


CREATE TABLE [dn_table_rowcount] (
	[tabname] [varchar] (64) NOT NULL ,
	[row_cnt] [int] NOT NULL ,
	[build_id] [int] NOT NULL ,
	[update_time] [datetime] NOT NULL 
)
GO


CREATE TABLE [dn_validationSummary] (
	[tax_id] [int] NOT NULL ,
	[snp_class] [tinyint] NOT NULL ,
	[validation_status] [tinyint] NOT NULL ,
	[rs_cnt] [int] NULL ,
	[build_id] [int] NULL ,
	[create_time] [smalldatetime] NULL 
)
GO


CREATE TABLE [strainConv] (
	[strain_curated] [varchar] (255) NULL ,
	[strain_entry] [varchar] (255) NULL ,
	[cell_line] [varchar] (255) NULL ,
	[tax_id] [int] NOT NULL ,
	[crdate] [smalldatetime] NULL ,
	[modate] [smalldatetime] NULL 
)
GO


CREATE TABLE [uniStsLoc] (
	[sts_uid] [int] NULL ,
	[sts_name] [varchar] (50) NULL ,
	[tax_id] [int] NULL ,
	[sts_weight] [int] NULL ,
	[contig_acc] [varchar] (15) NULL ,
	[contig_ver] [tinyint] NULL ,
	[contig_from] [int] NULL ,
	[contig_to] [int] NULL ,
	[create_time] [smalldatetime] NULL 
)
GO


