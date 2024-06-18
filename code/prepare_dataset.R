.libPaths( c('.Rpackages',.libPaths() ) )

library(ewastools)
library(stringi)
library(data.table)
library(magrittr)
library(purrr)

# Download phenotype data (copy URL in browser to see what the requested file looks like)
pheno = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116339&targ=gsm&form=text&view=brief"
pheno = readLines(pheno)

# Split into individual samples
pheno = split(pheno,cumsum(pheno %like% "^\\^SAMPLE = GSM"))

# Extract GSM accessions
names(pheno) = map(pheno,1) %>% stri_match_first(regex="GSM\\d+")

# Parse pheno data
imap(pheno,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> pheno

pheno = rbindlist(pheno)

# Keep only information on sample characteristics and supplementary files
pheno = pheno[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = pheno[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = pheno$value[i] %>% stri_split(fixed=": ")
pheno$variable[i] = map_chr(ch,1)
pheno$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
pheno[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
pheno[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

# Reshape data.table from long to wide format
pheno = dcast(pheno, gsm ~ variable)

# Select and parse the relevant variables
pheno = pheno[,.(gsm,age,gender,red,grn)]
pheno = rbind(pheno
	,list(
	 "GSM8002189"
	,"66"
	,"Male"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8002nnn/GSM8002189/suppl/GSM8002189_205125590022_R06C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8002nnn/GSM8002189/suppl/GSM8002189_205125590022_R06C01_Grn.idat.gz")
	,list(
	 "GSM5322148"
	,"54"
	,"Male"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5322nnn/GSM5322148/suppl/GSM5322148_204237140013_R03C01_Red.idat.gz"
	,"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5322nnn/GSM5322148/suppl/GSM5322148_204237140013_R03C01_Grn.idat.gz")
	 )


# Select the samples
setkey(pheno,"gsm")
pheno = pheno[c(
	## 15 smokers
  "GSM3229078","GSM3228657","GSM3228869", "GSM3229168",
  "GSM3229167","GSM3228683","GSM3229086","GSM3228891",
  "GSM3228837","GSM3228698","GSM3228928","GSM3228838",
  "GSM3229054","GSM3229161","GSM3228886",
	
	##  15 non-smokers
  "GSM3228804","GSM3228829","GSM3229169","GSM3228981",
  "GSM3228768","GSM3228739", "GSM3228889", "GSM3229117",
  "GSM3229074","GSM3228947","GSM3229087","GSM3228648",
  "GSM3228943","GSM3229147","GSM3229108"
  

# ,"GSM2260543" # same person as GSM2260485 (from 450K processing)
# ,"GSM2260653" # this is the potentially contaminated sample (from 450K processing)
,"GSM5322148" # unrelated sample from another GSE, CD4T instead of whole blood
,"GSM8002189" # unrelated sample of brain tissue
,"GSM3229163" # sample for which we'll change sex
)]

dir.create("data",showWarnings=FALSE)

# Adding chip and position information
pheno$Sentrix_ID = substring(pheno$red, 79, 90)
pheno$row = substring(pheno$red, 92, 94)
pheno$col = substring(pheno$red, 95, 97)

# Download .idat files
map2(pheno$red,"data/" %s+% pheno$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(pheno$grn,"data/" %s+% pheno$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
pheno$red = NULL; pheno$grn = NULL

# Import the methylation data
meth = read_idats("data/" %s+% pheno$gsm)
pheno[,c("X","Y"):=check_sex(meth)]
pheno[,sex:=ifelse(X>1.,"f","m")]

# add predicted smoking activity
smoke <- read.csv("full_phenotypes.csv")
smoke <- smoke %>% dplyr::select(geo_accession, Smoking) %>% 
  dplyr::rename(gsm = geo_accession, smoker = Smoking)
pheno <- dplyr::left_join(pheno, smoke, by = "gsm")

pheno = pheno[,.(gsm,age,sex,smoker,Sentrix_ID,row,col)]
pheno[gsm=="GSM3229163",sex:="m"]

write.csv(pheno,file="data/pheno.csv",row.names=FALSE)
write.csv(pheno[1:30],file="data/pheno_clean.csv",row.names=FALSE)
