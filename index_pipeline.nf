nextflow.enable.dsl=2 // Force nextflow to run DSL2

// Process 'Index' parameters
params.ref="./data/chr21/chr21.fasta" // Path to the unprocessed chromosome .fasta file
params.index_dir="./index_dir" // Publish directory for the emitted files


process index {

publishDir("${params.index_dir}", mode: 'move') // Emit files to the publish directory - 'move' as 'copy' occasionally fails to copy file(s)

input:
 path genome // Expects the path to the unprocessed chromosome .fasta file

output:
 path "*" // Emit all generated files

script: // Index the input .fasta file - BWA must be in installed in the active environment
"""
bwa index ${genome}
"""

}


workflow {

ref_ch=Channel.fromPath(params.ref) // Define the reference parameter as a path channel

index(ref_ch) // Call the index process with ref_ch as the parameter

index.out.view() // Echo the files emitted by process 'index'

}