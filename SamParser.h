/* ReadType here means if the read is unalignable, alignable or aligned too much. It is NOT siheaderngle read or paired-end read */
#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"

#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"
#include "SingleHit.h"
#include "PairedEndHit.h"

#include "Transcripts.h"
#include "TranscriptsGenome.h"

class SamParser {
public:
	SamParser(const char* inpF, const char* aux, TranscriptsGenome& transcripts, const char* imdName, bool useGenome = false);
	~SamParser();

	/**
	 * return value
	 * -1 : no more alignment
	 * 0 : new read , type 0
	 * 1 : new read , type 1 with alignment
	 * 2 : new read , type 2
	 * 5 : new alignment but same read
	 */
	int parseNext(SingleRead& read, SingleHit& hit);
	int parseNext(SingleReadQ& read, SingleHit& hit);
	int parseNext(PairedEndRead& read, PairedEndHit& hit);
	int parseNext(PairedEndReadQ& read, PairedEndHit& hit);

	static void setReadTypeTag(const char* tag) {
		strcpy(rtTag, tag);
	}

private:
	samFile *sam_in;
	bam_hdr_t *header;
	bam1_t *b, *b2;

	TranscriptsGenome& transcripts;

	int n_warns; // Number of warnings
	
	//tag used by aligner
	static char rtTag[STRLEN];

	//0 ~ N0, 1 ~ N1, 2 ~ N2
	int getReadType(const bam1_t* b) {
	  if (bam_is_mapped(b)) return 1;
	  if (!strcmp(rtTag, "")) return 0;
	  uint8_t *p = bam_aux_get(b, rtTag);
	  return (p == NULL || bam_aux2i(p) <= 0) ? 0 : 2;
	}

	// for paired-end reads
	int getReadType(const bam1_t* b, const bam1_t* b2) {
	  if (bam_is_mapped(b) && bam_is_mapped(b2)) return 1;
	  if (!strcmp(rtTag, "")) return 0;
	  
	  uint8_t *p = bam_aux_get(b, rtTag);
	  if (p != NULL && bam_aux2i(p) > 0) return 2;
	  
	  p = bam_aux_get(b2, rtTag);
	  if (p != NULL && bam_aux2i(p) > 0) return 2;
	  
	  return 0;
	}

	void initSamFile(const char* inpF, const char* aux) {
		sam_in = sam_open(inpF, "r");
		general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

		if (aux != NULL) hts_set_fai_filename(sam_in, aux);
		header = sam_hdr_read(sam_in);
		general_assert(header != 0, "Fail to parse sam header!");
	}

	void closeSamFile() {
		sam_close(sam_in);
		sam_in = NULL;
		header = NULL;
	}

	int getHitPosition(bam1_t* b, bool isReversed) {
		if (transcripts.getIsUsingGenomeFile()) {
			Transcript t = transcripts.getTranscriptAt(
				transcripts.getInternalSid(b, header->target_name[b->core.tid])
			);
			if (isReversed) {
				return t.getLength() - (b->core.pos - t.getPosition()) - b->core.l_qseq;
			} else {
				return b->core.pos - t.getPosition();
			}
		} else if (isReversed) {
			return header->target_len[b->core.tid] - b->core.pos - b->core.l_qseq;
		} else {
			return b->core.pos;
		}
	}
};

char SamParser::rtTag[STRLEN] = ""; // default : no tag, thus no Type 2 reads

// aux, if not 0, points to the file name of fn_list
SamParser::SamParser(const char* inpF, const char* aux, TranscriptsGenome& transcripts, const char* imdName, bool useGenome)
	: transcripts(transcripts), n_warns(0)
{
	initSamFile(inpF, aux);

	if (!useGenome) {
		printf("Building mappings from non-genome file\n");
		transcripts.buildMappings(header->n_targets, header->target_name, imdName);
	} else {
		printf("Building mappings from genome file\n");
		
		transcripts.buildMappingsGenome(sam_in, header, imdName);
		// Close and re-open the file, as the genome file has been read through by the above function
		closeSamFile();
		initSamFile(inpF, aux);
	}

	b = bam_init1();
	b2 = bam_init1();
}

SamParser::~SamParser() {
	if (n_warns > 0) fprintf(stderr, "Warning: Detected %d lines containing read pairs whose two mates have different names.\n", n_warns);
	
  	bam_hdr_destroy(header);
	sam_close(sam_in);
	bam_destroy1(b);
	bam_destroy1(b2);
}

// If sam_read1 returns 0 , what does it mean?
//Assume b.core.tid is 0-based
int SamParser::parseNext(SingleRead& read, SingleHit& hit) {
	int val; // return value

	if (sam_read1(sam_in, header, b) < 0) return -1;

	std::string name = bam_get_canonical_name(b);
	
	general_assert(!bam_is_paired(b), "Read " + name + ": Find a paired end read in the file!");

	int readType = getReadType(b);
	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		read = SingleRead(name, bam_get_read_seq(b));
	}
	else {
		general_assert(read.getReadLength() == b->core.l_qseq, "Read " + name + " has alignments with inconsistent read lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!\n");
	  if (bam_is_rev(b)) {
	    hit = SingleHit(
			-transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, true)
		);
	  }
	  else {
	    hit = SingleHit(
			transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, false)
		);
	  }
	}

	return val;
}

int SamParser::parseNext(SingleReadQ& read, SingleHit& hit) {
	int val;

	if (sam_read1(sam_in, header, b) < 0) return -1;

	std::string name = bam_get_canonical_name(b);
	
	general_assert(!bam_is_paired(b), "Read " + name + ": Find a paired end read in the file!");

	int readType = getReadType(b);
	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		read = SingleReadQ(name, bam_get_read_seq(b), bam_get_qscore(b));
	}
	else {
		general_assert(read.getReadLength() == b->core.l_qseq, "Read " + name + " has alignments with inconsistent read lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!\n");
	  if (bam_is_rev(b)) {
	    hit = SingleHit(
			-transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, true)
		);
	  }
	  else {
	    hit = SingleHit(
			transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, false)
		);
	  }
	}

	return val;
}

//Assume whether aligned or not , two mates of paired-end reads are always get together
int SamParser::parseNext(PairedEndRead& read, PairedEndHit& hit) {
	int val;

	if ((sam_read1(sam_in, header, b) < 0) || (sam_read1(sam_in, header, b2) < 0)) return -1;

	if (!bam_is_read1(b)) { bam1_t * tmp = b; b = b2; b2 = tmp; }
	std::string name = bam_get_canonical_name(b);
	
	general_assert(bam_is_paired(b) && bam_is_paired(b2), "Read " + name + ": One of the mate is not paired-end! (RSEM assumes the two mates of a paired-end read should be adjacent)");
	general_assert((bam_is_read1(b) && bam_is_read2(b2)), "Read " + name + ": The adjacent two lines do not represent the two mates of a paired-end read! (RSEM assumes the two mates of a paired-end read should be adjacent)");
	general_assert((bam_is_mapped(b) && bam_is_mapped(b2)) || (!bam_is_mapped(b) && !bam_is_mapped(b2)), "Read " + name + ": RSEM currently does not support partial alignments!");
	
	std::string name2 = bam_get_canonical_name(b2);	
	if (name != name2) 
	  if (++n_warns <= MAX_WARNS) 
	    fprintf(stderr, "Warning: Detected a read pair whose two mates have different names--%s and %s!\n", name.c_str(), name2.c_str());

	int readType = getReadType(b, b2);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		SingleRead mate1(name, bam_get_read_seq(b));
		SingleRead mate2(name2, bam_get_read_seq(b2));
		read = PairedEndRead(mate1, mate2);
	}
	else {
		general_assert(read.getMate1().getReadLength() == b->core.l_qseq && read.getMate2().getReadLength() == b2->core.l_qseq, "Paired-end read " + name + " has alignments with inconsistent mate lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b) && bam_check_cigar(b2), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!");
	  general_assert(b->core.tid == b2->core.tid, "Read " + name + ": The two mates do not align to a same transcript! RSEM does not support discordant alignments.");
	  if (bam_is_rev(b)) {
	    hit = PairedEndHit(
			-transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, true),
			b->core.pos + b->core.l_qseq - b2->core.pos
		);
	  }
	  else {
	    hit = PairedEndHit(
			transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, false),
			b2->core.pos + b2->core.l_qseq - b->core.pos
		);
	  }
	}

	return val;
}

int SamParser::parseNext(PairedEndReadQ& read, PairedEndHit& hit) {
	int val;
	
	if ((sam_read1(sam_in, header, b) < 0) || (sam_read1(sam_in, header, b2) < 0)) return -1;

	if (!bam_is_read1(b)) { bam1_t *tmp = b; b = b2; b2 = tmp; } // swap if the first read is not read 1
	std::string name = bam_get_canonical_name(b);
	
	general_assert(bam_is_paired(b) && bam_is_paired(b2), "Read " + name + ": One of the mate is not paired-end! (RSEM assumes the two mates of a paired-end read should be adjacent)");
	general_assert(bam_is_read1(b) && bam_is_read2(b2), "Read " + name + ": The adjacent two lines do not represent the two mates of a paired-end read! (RSEM assumes the two mates of a paired-end read should be adjacent)");
	general_assert((bam_is_mapped(b) && bam_is_mapped(b2)) || (!bam_is_mapped(b) && !bam_is_mapped(b2)), "Read " + name + ": RSEM currently does not support partial alignments!");
	
	std::string name2 = bam_get_canonical_name(b2);	
	if (name != name2)
	  if (++n_warns <= MAX_WARNS)
	    fprintf(stderr, "Warning: Detected a read pair whose two mates have different names--%s and %s!\n", name.c_str(), name2.c_str());

	int readType = getReadType(b, b2);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		SingleReadQ mate1(name, bam_get_read_seq(b), bam_get_qscore(b));
		SingleReadQ mate2(name2, bam_get_read_seq(b2), bam_get_qscore(b2));
		read = PairedEndReadQ(mate1, mate2);
	}
	else {
		general_assert(read.getMate1().getReadLength() == b->core.l_qseq && read.getMate2().getReadLength() == b2->core.l_qseq, "Paired-end read " + name + " has alignments with inconsistent mate lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b) && bam_check_cigar(b2), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!");
	  general_assert(b->core.tid == b2->core.tid, "Read " + name + ": The two mates do not align to a same transcript! RSEM does not support discordant alignments.");
	  if (bam_is_rev(b)) {
	    hit = PairedEndHit(
			-transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, true),
			b->core.pos + b->core.l_qseq - b2->core.pos
		);
	  }
	  else {
	    hit = PairedEndHit(
			transcripts.getInternalSid(b, header->target_name[b->core.tid]),
			getHitPosition(b, false),
			b2->core.pos + b2->core.l_qseq - b->core.pos
		);
	  }
	}
	
	return val;
}

#endif /* SAMPARSER_H_ */
