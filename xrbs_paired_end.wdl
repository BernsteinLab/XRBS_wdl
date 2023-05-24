version 1.0

workflow XRBS {
	String pipeline_ver = 'dev'
	String image_id = sub(pipeline_ver, "dev", "latest")

	input {
        String sample_id
		File fastq1
		File fastq2
		File reference
        File reference_index
        File reference_positions
	}

	call fastqc { input:
        fastq1 = fastq1,
        fastq2 = fastq2,
        sample_id = sample_id
    }

    call trimming {input: 
                    fastq1_in = fastq1, 
                    fastq2_in = fastq2,
                    sample_id = sample_id
                }

    call fqconv {input: 
                    fastq1 = trimming.fastq1, 
                    fastq2 = trimming.fastq2,
                    sample_id = sample_id
                    }

    call align {input: 
                    fastq = fqconv.fastq,
                    reference_index = reference_index,
                    sample_id = sample_id
                }

	call bconv {input: 
		bam = align.bam_out,
        sample_id = sample_id
	}

	call filter {input: 
		sorted_bam = bconv.sorted_bam_out,
        sorted_bam_index = bconv.sorted_bam_index,
        sample_id = sample_id
	}

	call methylation {input: 
		filtered_bam = filter.bam,
        #filtered_bam = bconv.sorted_bam_out,
		reference_pos = reference_positions,
        sample_id = sample_id
	}

    call qc_stats {input:
        fastq1 = fastq1,
        trimmed_file = trimming.fastq1,
        stats_file = filter.stats,
        sample_id = sample_id,
        methylation_call = methylation.calls
    }

	output {
		File meth_calls = methylation.calls
        File meth_metrics = methylation.metrics
        File flagstat = bconv.flagstat
        File bconv_metrics = bconv.metrics
        File fastqc_out = fastqc.fastqc_out
        File cutadapt_report = trimming.cutadapt_report

        String reads_number = qc_stats.reads_number
        String trimmed_number = qc_stats.trimmed_number
        String mapped_unique = qc_stats.mapped_unique
        String msp1Pos_filtered = qc_stats.msp1Pos_filtered
        String size_filtered = qc_stats.size_filtered
        String dup_filtered = qc_stats.dup_filtered
        String cpg_sites = qc_stats.cpg_sites
	}
}

task fastqc {
    input {
        File fastq1
        File fastq2
        String sample_id
    }
    command {
        mkdir ${sample_id}_fastqc_out
        fastqc -o ${sample_id}_fastqc_out ${fastq1} ${fastq2}
        zip -r ${sample_id}_fastqc_out.zip ${sample_id}_fastqc_out
    }
    runtime {
        docker: "salvacasani/fastqc"
    }
    output {
        File fastqc_out = "${sample_id}_fastqc_out.zip"
    }
}

task trimming {
    input {
        File fastq1_in
        File fastq2_in
        String sample_id
    }
    command {
        cutadapt --discard -a GCTCTTCCGATCT -o ${sample_id}_R2.cutadapt.fastq.gz -p ${sample_id}_R1.cutadapt.fastq.gz ${fastq2_in} ${fastq1_in} > ${sample_id}_cutadapt_report.txt
        /TrimGalore-0.6.10/trim_galore --path_to_cutadapt cutadapt --paired --illumina --nextseq 20 ${sample_id}_R1.cutadapt.fastq.gz ${sample_id}_R2.cutadapt.fastq.gz
    }
    runtime {
        docker: "salvacasani/trimming:latest"
    }
    output {
        File cutadapt_report = "${sample_id}_cutadapt_report.txt"
        File fastq1 = "${sample_id}_R1.cutadapt_val_1.fq.gz"
        File fastq2 = "${sample_id}_R2.cutadapt_val_2.fq.gz"
    }
}


task fqconv {
    input {
        File fastq1
        File fastq2
        String sample_id
    }
    command {
        /methylCtools/methylCtools fqconv -1 ${fastq1} -2 ${fastq2} ${sample_id}.conv.fq

    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
		memory: "20G"
		disks: "local-disk " + "500" + " SSD"
    }
    output {
        File fastq = "${sample_id}.conv.fq"
    }
}


task align {
    input {
        File reference_index
        File fastq
        String sample_id
    }
    command {
        mkdir genome_index
        tar zxvf ${reference_index} -C genome_index
        # Get genome index fasta name. 
        BWT=$(find genome_index -name '*.bwt')  
        GENOME_INDEX_FA="$(dirname $BWT)"/"$(basename $BWT .bwt)"
        echo "Using bwa index: $GENOME_INDEX_FA"

        bwa mem -t 4 -p -M -T 0 $GENOME_INDEX_FA ${fastq} > align.sam
        samtools view -Sb align.sam > ${sample_id}.reads.conv.bam

    }
    runtime {
        docker: "salvacasani/bwa:latest"
        bootDiskSizeGb: 40
		memory: "20GB"
		disks: "local-disk " + "500" + " SSD"
    }
    output {
        File bam_out = "${sample_id}.reads.conv.bam"
    }
}


task bconv {
    input {
        File bam
        String sample_id
    }
    command {
        /methylCtools/methylCtools bconv ${bam} -m ${sample_id}.human.conv.sort.metrics.txt - | samtools sort -T ${sample_id}.human.sort -@ 4 - > ${sample_id}.sorted.bam
        samtools index ${sample_id}.sorted.bam
        samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.sorted.bam.flagstat
    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
        memory: "20GB"
        disks: "local-disk " + "500" + " SSD"
    }
    output {
        File sorted_bam_out = "${sample_id}.sorted.bam"
        File sorted_bam_index = "${sample_id}.sorted.bam.bai"
        File metrics = "${sample_id}.human.conv.sort.metrics.txt"
        File flagstat = "${sample_id}.sorted.bam.flagstat"
    }
}


task filter {
    input {
        File sorted_bam
        File sorted_bam_index
        String sample_id
    }
    command {
        Rscript /home/filter.vh20200112.R ${sorted_bam}
    }
    runtime {
        docker: "salvacasani/r_filter:latest"
        bootDiskSizeGb: 40
        memory: "200GB"
        disks: "local-disk " + "1000" + " SSD"
    }
    output {
        File bam = "${sample_id}.sorted.bam.filter.bam"
        File stats = "${sample_id}.sorted.bam.filter.stats"
    }
}

task methylation {
    input {
        File filtered_bam
        File reference_pos
        String sample_id
    }
    command {
        unzip ${reference_pos}
        /methylCtools/methylCtools bcall --trimPE --metrics ${sample_id}.human.sort.filter.call.metrics GRCh38.pos.gz ${filtered_bam} - | bgzip > ${sample_id}.call.gz

    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
        memory: "20GB"
        disks: "local-disk " + "500" + " SSD"
    }
    output {
        File metrics = "${sample_id}.human.sort.filter.call.metrics"
        File calls = "${sample_id}.call.gz"
    }
}


task qc_stats {
    input {
        File fastq1
        File trimmed_file
        File stats_file
        String sample_id
        File methylation_call
    }
    command {

        reads_number=$(echo $(zcat ${fastq1}|wc -l)/4|bc);

        trimmed=$(echo $(zcat ${trimmed_file}|wc -l)/4|bc);

        mapped_unique=$(echo $(cat ${stats_file} | sed 1d | cut -f 2));

        msp1Pos_filtered=$(echo $(cat ${stats_file} | sed 1d | cut -f 5));

        size_filtered=$(echo $(cat ${stats_file} | sed 1d | cut -f 6));

        dup_filtered=$(echo $(cat ${stats_file} | sed 1d | cut -f 7));

        cpg_sites=$(echo $(zcat ${methylation_call} | wc -l));

    }
    runtime {
        docker: "salvacasani/r_filter:latest"
        bootDiskSizeGb: 40
        memory: "20GB"
        disks: "local-disk " + "500" + " SSD"
    }
    output {
        String reads_number = "${reads_number}"
        String trimmed_number = "${trimmed_number}"
        String mapped_unique = "${mapped_unique}"
        String msp1Pos_filtered = "${msp1Pos_filtered}"
        String size_filtered = "${size_filtered}"
        String dup_filtered = "${dup_filtered}"
        String cpg_sites = "${cpg_sites}"
    }
}


