#include "BamParser.hpp"
#include <stdexcept>

BamParser::BamParser(const std::string& filepath, int num_threads) {
    in_file.reset(sam_open(filepath.c_str(), "r"));
    if (!in_file) {
        throw std::runtime_error("Not possible to open BAM file: " + filepath);
    }

    if (num_threads > 1)
        hts_set_threads(in_file.get(), num_threads);

    header.reset(sam_hdr_read(in_file.get()));
    if (!header) {
        throw std::runtime_error("Not possible to read BAM header.");
    }

    record.reset(bam_init1());
}

bool BamParser::get_next_contact(std::string& contig_a, std::string& contig_b) {
    while (sam_read1(in_file.get(), header.get(), record.get()) >= 0) {
        const bam1_core_t* core = &record->core;

        bool is_paired = (core->flag & BAM_FPAIRED) != 0;
        bool is_unmapped = (core->flag & BAM_FUNMAP) != 0;
        bool is_mate_unmapped = (core->flag & BAM_FMUNMAP) != 0;
        
        bool is_read1 = (core->flag & BAM_FREAD1) != 0;
        bool is_secondary = (core->flag & BAM_FSECONDARY) != 0;
        bool is_supplementary = (core->flag & BAM_FSUPPLEMENTARY) != 0;

        if (is_paired && !is_unmapped && !is_mate_unmapped && is_read1 && !is_secondary && !is_supplementary) {
            if (core->tid >= 0 && core->mtid >= 0) {
                contig_a = header->target_name[core->tid];
                contig_b = header->target_name[core->mtid];
                
                return true;
            }
        }
    }
    return false;
}