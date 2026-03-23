#pragma once

#include <string>
#include <memory>
#include <htslib/sam.h>

struct SamFileDeleter { 
    void operator()(samFile* p) const { if (p) sam_close(p); } 
};
struct BamHdrDeleter { 
    void operator()(bam_hdr_t* p) const { if (p) bam_hdr_destroy(p); } 
};
struct Bam1Deleter { 
    void operator()(bam1_t* p) const { if (p) bam_destroy1(p); } 
};

class BamParser {
private:
    std::unique_ptr<samFile, SamFileDeleter> in_file;
    std::unique_ptr<bam_hdr_t, BamHdrDeleter> header;
    std::unique_ptr<bam1_t, Bam1Deleter> record;

public:
    explicit BamParser(const std::string& filepath);

    bool get_next_contact(std::string& contig_a, std::string& contig_b);
};