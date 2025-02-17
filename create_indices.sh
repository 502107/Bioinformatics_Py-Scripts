#!/bin/bash
#SBATCH --mem 2G
#SBATCH -J create_indices
#SBATCH -p jic-long,nbi-long,RG-Diane-Saunders
#SBATCH -o create_indices_%j.out
#SBATCH -e create_indices_%j.err

# example: sbatch /jic/scratch/groups/Diane-Saunders/your_folder/your_reference/create_indices.sh \
# --refPath /jic/scratch/groups/Diane-Saunders/your_folder/your_reference \
# --refGenomeName your_reference.fasta \
# --refDescriptionName your_reference.gff3


while [[ $# -gt 1 ]]
do
    key=$1
    case $key in
        -r|--refPath)
        refPath="$2"
        shift
        ;;
        -g|--refGenomeName)
        refGenomeName="$2"
        shift
        ;;
        -d|--refDescriptionName)
        refDescriptionName="$2"
        shift
        ;;
        *) echo "Unkown option: $1 " >&2
        exit 1
        ;;
    esac
    shift
done

echo "parameters:"
echo "refPath               $refPath"
echo "refGenomeName         $refGenomeName"
echo "refDescriptionName    $refDescriptionName"

if [ -z $refPath ]
then
    echo "ERROR: refPath not given. refPath (path to reference directory) must be given. This is needed to find the genome and description files, as well as for outputting the created indices."
    exit 1
fi
if [ -z $refGenomeName ]
then
    echo "ERROR: refGenomeName not given. refGenomeName (name of genome sequence file, e.g. pst104e.fasta) must be given. This should exist in refPath"
    exit 1
fi
if [ -z $refDescriptionName ]
then
    echo "ERROR: refDescriptionName not given. refDescriptionName (name of genome description file, e.g. pst104e.gff3) must be given. This should exist in refPath"
    exit 1
fi

refGenomePath="$refPath"/"$refGenomeName"
refDescriptionPath="$refPath"/"$refDescriptionName"

starOutPath="$refPath"/for_star
bwaOutPath="$refPath"/for_bwa
mkdir "$starOutPath"
mkdir "$bwaOutPath"

logsPath="$refPath"/create_indices_logs
mkdir "$logsPath"
cd "$logsPath"

parentFeature="ID"

printf "\nderived parameters: \n"
echo "refGenomePath         $refGenomePath"
echo "refDescriptionPath    $refDescriptionPath"
echo "starOutPath           $starOutPath"
echo "parentFeature         $parentFeature"
echo "bwaOutPath            $bwaOutPath"
echo "logsPath              $logsPath"

# Copy the description and genome files to the star and bwa directories
cp "$refDescriptionPath" "$refGenomePath" "$starOutPath"
cp "$refDescriptionPath" "$refGenomePath" "$bwaOutPath"


# Create .fai file from genome file

# samtools - 1.10
source package 758be80b-33cc-495a-9adc-11882ab145b1
faidxJobId=$(sbatch --parsable -J faidx --mem 4G -p jic-long,nbi-long,RG-Diane-Saunders -o faidx.out -e faidx.err --wrap "samtools faidx $refGenomePath")
# Once the .fai file is made, copy it to bwa and star directories
sbatch --dependency=afterok:"$faidxJobId" -J copy_fai --mem 2G -p jic-long,nbi-long,RG-Diane-Saunders --wrap "cp $refPath/*.fai $starOutPath && cp $refPath/*.fai $bwaOutPath"



# Create files for STAR

source STAR-2.5.a
starJobId=$(sbatch -J star_index --mem 32G -c 8 -N 1 -p jic-long,nbi-long,RG-Diane-Saunders -o star_index.out -e star_index.err --wrap "STAR --runMode genomeGenerate --genomeDir $starOutPath --genomeFastaFiles $refGenomePath --sjdbGTFtagExonParentTranscript $parentFeature --sjdbGTFfile $refDescriptionPath --runThreadN 8 --sjdbOverhang 99")
# picard - 2.21.2
source package ce4daee0-abd9-4cc6-8b13-5fe8ede3b149
# Create the dictionary file in star directory using the genome file that was copied to the star directory
dictionaryJobId=$(sbatch -J dictionary --mem 4G -c 8 -N 1 -p jic-long,nbi-long,RG-Diane-Saunders -o dictionary.out -e dictionary.err --wrap "picard CreateSequenceDictionary R=$starOutPath/$refGenomeName")



# Create files for BWA

source bwa-0.7.5
bwaJobId=$(sbatch --parsable -J bwa_index --mem 8G -c 1 -N 1 -p jic-long,nbi-long,RG-Diane-Saunders -o bwa_index.out -e bwa_index.err --wrap "bwa index $refGenomePath")
# Once the files for BWA are made, move them to the bwa directory
sbatch --dependency=afterok:"$bwaJobId" -J move_bwa --mem 2G -p jic-long,nbi-long,RG-Diane-Saunders --wrap "mv $refPath/*.ann $refPath/*.amb $refPath/*.bwt $refPath/*.pac $refPath/*.sa $bwaOutPath"

