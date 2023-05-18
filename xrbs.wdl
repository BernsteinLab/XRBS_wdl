workflow XRBS {
	String pipeline_ver = 'dev'
	String image_id = sub(pipeline_ver, "dev", "latest")

	input {
        String sample_id
		File fastq1
		File? fastq2
		File reference
        File reference_index
		String adapter
        String library
        Boolean paired
        File reference_positions
	}

	call fastqc { input:
		fastq1 = fastq1,
        fastq2 = fastq2,
        sample_id = sample_id
	}

    if (paired){
        call trimming {input: 
                        fastq1 = fastq1, 
                        fastq2 = fastq2, 
                        adapter = adapter,
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
    }

    if (!paired){
        call trimming2 {input: 
                        fastq1 = fastq1,
                        adapter = adapter,
                        sample_id = sample_id
                    }

        call fqconv2 {input: 
                        fastq1 = trimming2.fastq1,
                        sample_id = sample_id
                    }

        call align {input: 
                        fastq = fqconv2.fastq,
                        reference_index = reference_index,
                        sample_id = sample_id
                    }
    }

	call bconv {input: 
		bam = align.bam_out,
        sample_id = sample_id
	}

	call filter {input: 
		bam = bconv.sorted_bam_out,
        sample_id = sample_id
	}

	call methylation {input: 
		bam = filter.filter_bam_out,
		pos = reference_positions
	}

	output {
		File meth_calls = methylation.calls
        File meth_metrics = methylation.metrics
        File flagstat = bconv.flagstat
        File bconv_metrics = bconv.metrics
        File fastqc_out = fastqc.fastqc_out
        if (paired){
            File cutadapt_report = trimming.cutadapt_report
        }
        if (!paired){
            File cutadapt_report = trimming2.cutadapt_report
        }
	}


}

task fastqc {
    input {
        File fastq1
        File? fastq2
        String sample_id
    }
    command {
        fastqc -o ${sample_id}_fastqc_out fastq1 fastq2
        zip ${sample_id}_fastqc_out
    }
    runtime {
        docker: "pegi3s/fastqc"
        bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastqc_out = ${sample_id}_fastqc_out.zip
    }
}

task trimming {
    input {
        File fastq1
        File fastq2
        String adapter
        String platform
        String sample_id
    }
    command {
        ~/.local/bin/cutadapt --discard -a ${adapter} -o ${sample_id}_R2.cutadapt.fastq.gz -p ${sample_id}_R1.cutadapt.fastq.gz ${fastq2} ${fastq1} > ${sample_id}_cutadapt_report.txt
        software/TrimGalore-0.6.5/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt --paired --illumina --nextseq 20 ${sample_id}_R1.cutadapt.fastq.gz ${sample_id}_R2.cutadapt.fastq.gz
    }
    runtime {
        docker: "salvacasani/trimming:latest"
        bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
    }
    output {
        File cutadapt_report = cutadapt_report.txt
        File fastq1 = ${sample_id}_R1.cutadapt_val_1.fq.gz
        File fastq2 = ${sample_id}_R2.cutadapt_val_2.fq.gz
    }
}

task trimming2 {
    input {
        File fastq1
        String adapter
        String platform
        String sample_id
    }
    command {
        cutadapt --discard -a ${adapter} -o ${sample_id}.cutadapt.fastq.gz ${fastq1} > ${sample_id}_cutadapt_report.txt
        TrimGalore-0.6.5/trim_galore --path_to_cutadapt cutadapt --illumina --nextseq 20 ${sample_id}.cutadapt.fastq.gz
    }
    runtime {
        docker: "salvacasani/trimming:latest"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }
    output {
        File cutadapt_report = ${sample_id}_cutadapt_report.txt
        File fastq1 = ${sample_id}.cutadapt_val_1.fq.gz
    }
}


task fqconv {
    input {
        File fastq1
        File fastq2
        String sample_id
    }
    command {
        methylCtools fqconv -1 ${fastq1} -2 ${fastq2} ${sample_id}.conv.fq

    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastq = ${sample_id}.conv.fq
    }
}


task fqconv2 {
    input {
        File fastq1
        String sample_id
    }
    command {
        methylCtools fqconv -1 ${fastq1} ${sample_id}.conv.fq

    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastq = ${sample_id}.conv.fq
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

        bwa mem -t 4 -p -M -T 0 $GENOME_INDEX_FA ${fastq} | samtools view -Sb - > ${sample_id}.reads.conv.bam

    }
    runtime {
        docker: "salvacasani/bwa:latest"
        bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastq = ${sample_id}.conv.fq
    }
}


task bconv {
    input {
        File bam
        String sample_id
    }
    command {
        methylCtools bconv bam -m ${$sample_id}.human.conv.sort.metrics.txt - | samtools sort -T ${sample_id}.human.sort -@ 4 - > ${sample_id}.sorted.bam
        samtools index ${sample_id}.sorted.bam
        samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.sorted.bam.flagstat
    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastq = ${sample_id}.sorted.bam
        File metrics = ${sample_id}.human.conv.sort.metrics.txt
        File flagstat = ${sample_id}.sorted.bam.flagstat
    }
}



task filter {
    input {
        File sorted.bam
        String sample_id
    }
    command {
        Rscript filter.vh20200112.R  sorted.bam

    }
    runtime {
        docker: "salvacasani/r_filter:latest"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }
    output {
        File fastq = ${sample_id}.filter.sorted.bam
    }
}


task methylation {
    input {
        File filtered.bam
        File reference_pos
        String sample_id
    }
    command {
        methylCtools bcall --trimPE --metrics ${sample_id}.human.sort.filter.call.metrics reference_pos filtered.bam - | bgzip > ${sample_id}.call.gz

    }
    runtime {
        docker: "salvacasani/methylctools:latest"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }
    output {
        File metrics = ${sample_id}.human.sort.filter.call.metrics
        File calls = ${sample_id}.call.gz
    }
}