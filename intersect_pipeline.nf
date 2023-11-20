nextflow.enable.dsl=2 // Force nextflow to run DSL2

// find_intersections parameters
params.trimmed_beds_dir="/home/biouser/bsb_test/trimmed_beds/trimmed_Sample*.bed" // Path to the trimmed sample.bed files
params.asiSI_sites_dir="/home/biouser/bsb_test/data/chr21_AsiSI_sites.t2t.bed" // Path to the asiSI_sites_file.bed file
params.intersected_beds_dir="./intersected_beds" // Publish directory for the emitted files


process find_intersections {

publishDir("${params.intersected_beds_dir}", mode: 'move') // Emit files to the publish directory - 'move' as 'copy' occasionally fails to copy file(s)

input:
 tuple val(name), path(trimmed_bed_file) // Expects a tuple containing sample name and path to the trimmed .bed file
 path asiSI_sites_file // Expects the path to the asiSI_sites_file.bed file

output:
 path "intersected_${name}.bed" // Emits the path to the intersected_trimmed_Sample*.bed file

script: // Perform bedtools intersection and save the result to a .bed file
"""
bedtools intersect -a ${asiSI_sites_file} -b ${trimmed_bed_file} -c > intersected_${name}.bed
"""

}


workflow {

trimmed_beds_ch = Channel
                        .fromPath( params.trimmed_beds_dir ) // Define as a path channel
                        .map { tuple(it.getBaseName(), it ) } // Define as a tuple to access sample name and path

asiSI_sites_ch = Channel
                        .value( params.asiSI_sites_dir ) // Define the asiSI_sites_file.bed location as a value channel - needs to be consumed multiple times

find_intersections(trimmed_beds_ch, asiSI_sites_ch) // Call Process 'find_intersections' with the trimmed_beds_ch and asiSI_sites_ch channels as the parameters

}