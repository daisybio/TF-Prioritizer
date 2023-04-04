#! /usr/bin/make -f

ROOT = $(realpath ehmm_demo)
DATA_DIR = $(ROOT)/bam

ATAC-seq_DIR = $(DATA_DIR)/ATAC-seq
H3K27AC_DIR = $(DATA_DIR)/H3K27AC
H3K4ME1_DIR = $(DATA_DIR)/H3K4ME1
H3K4ME3_DIR = $(DATA_DIR)/H3K4ME3
TREATMENT_DIR = $(DATA_DIR)/ATAC-seq $(DATA_DIR)/H3K27AC $(DATA_DIR)/H3K4ME1 $(DATA_DIR)/H3K4ME3

ATAC-seq_ENC = $(ATAC-seq_DIR)/ENCFF550NVA
H3K27AC_ENC = $(H3K27AC_DIR)/ENCFF524ZFV
H3K4ME1_ENC = $(H3K4ME1_DIR)/ENCFF788JMC
H3K4ME3_ENC = $(H3K4ME3_DIR)/ENCFF211WGC
HM_ENC = $(H3K27AC_ENC) $(H3K4ME1_ENC) $(H3K4ME3_ENC)
HM_ENC_BAMS = $(addsuffix .bam,$(HM_ENC))
ATAC-seq_ENC_BAM = $(addsuffix .bam,$(ATAC-seq_ENC))

TARGETS = $(addsuffix .bam,$(TREATMENT_DIR))

# ------------------------------------------------------------------------------

all: $(TARGETS)

# ------------------------------------------------------------------------------

$(TARGETS): %.bam:
	ln -s "$<" "$@"
	samtools index $@

$(ATAC-seq_DIR).bam: $(ATAC-seq_ENC).sort.rmdup.bam
$(H3K27AC_DIR).bam: $(H3K27AC_ENC).sort.rmdup.bam
$(H3K4ME1_DIR).bam: $(H3K4ME1_ENC).sort.rmdup.bam
$(H3K4ME3_DIR).bam: $(H3K4ME3_ENC).sort.rmdup.bam

# ------------------------------------------------------------------------------

$(HM_ENC_BAMS): %.bam:
	mkdir -p $(@D)
	wget https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(basename $(notdir $@)).bam -O $@

$(ATAC-seq_ENC_BAM): %.bam:
	mkdir -p $(@D)
	wget https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(basename $(notdir $@)).bam -O $@
	samtools view -hq 30 $@ | grep -v 'XA:Z:' | grep -v 'chrM' | samtools view -Sb - > $@.tmp.bam
	mv $@.tmp.bam $@
	samtools sort -n $@ > $@.tmp.bam
	mv $@.tmp.bam $@
	samtools fixmate $@ $@.tmp.bam
	mv $@.tmp.bam $@  

%.sort.bam: %.bam
	samtools sort $< > $@

%.sort.rmdup.bam: %.sort.bam
	samtools rmdup -s $< $@

# keep intermediate files and delete files whenever an error occurs and set 'all' target to be phony
# ------------------------------------------------------------------------------
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all 
