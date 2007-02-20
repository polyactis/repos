set quoted_identifier  OFF 
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


SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Allele" as select * from dbSNP_main..Allele

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "AlleleFlagCode" as select * from dbSNP_main..AlleleFlagCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "AlleleMotif" as select * from dbSNP_main..AlleleMotif

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "AllocIds" as select * from dbSNP_main..AllocIds

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Author" as select * from dbSNP_main..Author

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "BaseURLList" as select * from dbSNP_main..BaseURLList

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "BuildIdToAllocId" as select * from dbSNP_main..BuildIdToAllocId

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "ChiSqPValue" as select * from dbSNP_main..ChiSqPValue

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "ChiSqPValueLookUp" as select * from dbSNP_main..ChiSqPValueLookUp

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "CpGCode" as select * from dbSNP_main..CpGCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "FtpValidation" as select * from dbSNP_main..FtpValidation

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "GenderCode" as select * from dbSNP_main..GenderCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "GtyAllele" as select * from dbSNP_main..GtyAllele

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "LoadHistory" as select * from dbSNP_main..LoadHistory

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "MapLinkCode" as select * from dbSNP_main..MapLinkCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Method" as select * from dbSNP_main..Method

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "MethodClass" as select * from dbSNP_main..MethodClass

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "MethodLine" as select * from dbSNP_main..MethodLine

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Moltype" as select * from dbSNP_main..Moltype

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Motif" as select * from dbSNP_main..Motif

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "ObsGenotype" as select * from dbSNP_main..ObsGenotype

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "ObsVariation" as select * from dbSNP_main..ObsVariation

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "OrganismTax" as select * from dbSNP_main..OrganismTax

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "PopClass" as select * from dbSNP_main..PopClass

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "PopClassCode" as select * from dbSNP_main..PopClassCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Publication" as select * from dbSNP_main..Publication

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "SnpClassCode" as select * from dbSNP_main..SnpClassCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "SnpFunctionCode" as select * from dbSNP_main..SnpFunctionCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "SnpValidationCode" as select * from dbSNP_main..SnpValidationCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "StrandCode" as select * from dbSNP_main..StrandCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "SubSNPDelComm" as select * from dbSNP_main..SubSNPDelComm

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "SubSNPSeqTypeCode" as select * from dbSNP_main..SubSNPSeqTypeCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "Submitter" as select * from dbSNP_main..Submitter

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "TemplateType" as select * from dbSNP_main..TemplateType

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "UniGty" as select * from dbSNP_main..UniGty

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "UniVariAllele" as select * from dbSNP_main..UniVariAllele

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "UniVariation" as select * from dbSNP_main..UniVariation

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "UniVariationSrcCode" as select * from dbSNP_main..UniVariationSrcCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "VarFlagCode" as select * from dbSNP_main..VarFlagCode

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "VariAllele" as select * from dbSNP_main..VariAllele

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "dn_Allele_rev" as select * from dbSNP_main..dn_Allele_rev

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "dn_Motif_rev" as select * from dbSNP_main..dn_Motif_rev

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "dn_UniGty_allele" as select * from dbSNP_main..dn_UniGty_allele

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "dn_UniGty_rev" as select * from dbSNP_main..dn_UniGty_rev

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "dn_UniVariation_rev" as select * from dbSNP_main..dn_UniVariation_rev

GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER ON 
GO
SET ANSI_NULLS ON 
GO

CREATE VIEW "idRange" as select * from dbSNP_main..idRange


GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS OFF 
GO

create view vwIndivBySource
as
select distinct c.code, ibs.src_id, ibs.ind_id, c.name, ibs.src_ind_id, ibs.src_ind_grp,
  case when left(ibs.src_ind_id, 4) = 'CEPH' then  substring(  ibs.src_ind_id, 5, charindex('.', ibs.src_ind_id) -5 )
else  ibs.src_ind_id
end ceph_ped_id, a.alias
from IndivBySource ibs left join IndivAlias a on ibs.src_id = a.src_id and ibs.src_ind_id = a.src_ind_id
join IndivSourceCode c on ibs.src_id = c.code


GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS OFF 
GO

create view vwIndivBySource_coriell_min_id
as
select ind_id, src_id, min(src_ind_id) src_ind_id from IndivBySource where src_id = 1
group by ind_id, src_id 
union
select ind_id,src_id,  src_ind_id from IndivBySource where src_id != 1


GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS OFF 
GO


create view vwSNP_max_len (snp_id, max_len) /* with SCHEMABINDING */
     as     
select l.snp_id, max(sequence_len) max_len
		from dbo.SubSNP s, dbo.SNPSubSNPLink l 
		where  s.subsnp_id = l.subsnp_id 
		group by l.snp_id  		



GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS OFF 
GO




create view SNPSeqPick(snp_id, subsnp_id) /* with SCHEMABINDING */
     as     
	
	select l.snp_id, min(s.subsnp_id) subsnp_id
		from	dbo.vwSNP_max_len t, dbo.SNPSubSNPLink l, dbo.SubSNP s
		where t.snp_id = l.snp_id
		and l.subsnp_id = s.subsnp_id
		and s.sequence_len = t.max_len 	
		group by l.snp_id



GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS OFF 
GO

create view vwPedigreeIndividual
as
/* Note: Special treatment for ceph(src_id=2):
 *         when display ceph pedigree, only show ma and pa within the same curator_ped_id.
 *        Since PedigreeIndividual table is on ind_id.and one ind_id can have several src_ind_id in ceph.
 * 
 */
select distinct pi.ped_id, pi.ind_id, p.curator, p.curator_ped_id, ibs.src_id, ibs.src_ind_id, ibs.src_ind_grp,
 maibs.src_ind_id ma_src_ind_id, 
 paibs.src_ind_id pa_src_ind_id, maibs.ind_id ma_ind_id ,paibs.ind_id pa_ind_id
from PedigreeIndividual pi join IndivBySource ibs on pi.ind_id = ibs.ind_id and ibs.src_id != 2
join Pedigree p on p.ped_id = pi.ped_id
left join vwIndivBySource_coriell_min_id maibs on pi.ma_ind_id = maibs.ind_id and ibs.src_id = maibs.src_id
left join vwIndivBySource_coriell_min_id paibs on pi.pa_ind_id = paibs.ind_id and ibs.src_id = paibs.src_id
union
select distinct pi.ped_id, pi.ind_id, p.curator, p.curator_ped_id, ibs.src_id, ibs.src_ind_id,  ibs.src_ind_grp,
 maibs.src_ind_id ma_src_ind_id, 
 paibs.src_ind_id pa_src_ind_id, maibs.ind_id ma_ind_id ,paibs.ind_id pa_ind_id
from PedigreeIndividual pi join vwIndivBySource ibs on pi.ind_id = ibs.ind_id and ibs.src_id =2
join Pedigree p on p.ped_id = pi.ped_id and ibs.ceph_ped_id = p.curator_ped_id
left join vwIndivBySource maibs on pi.ma_ind_id = maibs.ind_id and ibs.src_id = maibs.src_id and ibs.ceph_ped_id = maibs.ceph_ped_id
left join vwIndivBySource paibs on pi.pa_ind_id = paibs.ind_id and ibs.src_id = paibs.src_id and ibs.ceph_ped_id = paibs.ceph_ped_id





GO
SET QUOTED_IDENTIFIER OFF 
GO
SET ANSI_NULLS ON 
GO

