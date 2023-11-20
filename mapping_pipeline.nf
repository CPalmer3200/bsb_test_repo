nextflow.enable.dsl=2 // Force nextflow to run DSL2

// Process 'mapping_to_bam' parameters
params.genome_index="/home/biouser/bsb_test/index_dir/chr21.fasta" // Directory containing the indexed chromosome files
params.fastq="/home/biouser/bsb_test/data/fastqs/Sample*.fastq.gz" // Path to the samples in .fastq format


process mapping_to_bam {

input:
 val genome_index // Expects indexed genome file location as a value channel
 tuple val(name), path(fastq) // Expects a tuple containing sample name and path to the file

output:
 tuple val(name), path ("${name}.bam") // Emits a tuple containing sample name and the work directory path to the .bam file

script:  // Map the input .fastq to the indexed genome with bwa mem, pipe to samtools and convert to a .bam format - Requires bwa and bedtools to be installed in the active environment
"""
bwa mem ${genome_index} ${fastq} | samtools view -h -b -o ${name}.bam -
"""

}


process bam_to_bed {

input:
 tuple val(name), path(bam_file) // Expects a tuple containing sample name and path to the .bam file

output:
 tuple val(name), path("${name}.bed") // Emits a tuple containing sample name and the work directory path to the .bed file

script: // Convert the input .bam file into a .bed file
"""
bedtools bamtobed -i ${bam_file} > ${name}.bed
"""

}


process trim_reads {

input:
 tuple val(name), path(input_path) // Expects a tuple containing sample name and path to the .bed file

script: // Runs trim_reads.py to trim the reads to include just the break sites for each sample
"""
python3 /home/biouser/bsb_test/trim_reads.py ${input_path} ${name} /home/biouser/bsb_test/trimmed_beds/
"""

}


workflow {

index_ch=Channel.value(params.genome_index) // Define the indexed genome location as a value channel - needs to be consumed multiple times

fastq_ch =Channel
                .fromPath( params.fastq ) // Define as a path channel
                .map { tuple( it.getBaseName(2), it ) } // Define as a tuple to access sample name and path

mapping_to_bam(index_ch, fastq_ch) // Call Process 'mapping_to_bam' with index_ch and fastq_ch as parameters

bam_to_bed(mapping_to_bam.out) // Call Process 'bam_to_bed' with the 'mapping_to_bam' output as the parameter

trim_reads(bam_to_bed.out) // Call Process 'trim_reads' with the 'bam_to_bed' output as the parameter

}