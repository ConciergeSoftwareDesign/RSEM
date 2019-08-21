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
	 * 99: skip read
	 */
	template<class ReadType, class HitType>
	int parseNext(ReadType&, HitType&);

	// The following overload is used in the case where we are parsing reads 
	// aligned to the genome, and we want to output both the hit relative to a transcript
	// as well as a hit relative to the genome
	template<class ReadType, class HitType>
	int parseNext(ReadType& read, HitType& hit, HitType& genomeHit);

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

	/**
	 *  Get the hit position for a certain read, specifying if the read has the incorrect 
	 * 	orientation (isReversed). Returns -1 for reads where there is no valid position 
	 */
	int getHitPosition(bam1_t* b, int internalSid, bool isReversed, bool calculateTranscriptPos) {
		if (calculateTranscriptPos) {
			if (internalSid == -1) {
				return -1;
			}

			//TODO: This is a bit of strange place to do this
			transcripts.setHasAppeared(internalSid, true);

			Transcript t = transcripts.getTranscriptAt(internalSid);
			if (isReversed) {
				return t.getLength() - (b->core.pos - t.getPosition()) - b->core.l_qseq - 1;
			} else {
				return b->core.pos - t.getPosition() + 1;
			}
		} else if (isReversed) {
			return header->target_len[b->core.tid] - b->core.pos - b->core.l_qseq;
		} else {
			return b->core.pos;
		}
	}

	/**
	 * Get the length of a hit given the two reads as a pair, the hit position, 
	 * and whether or not it's reversed
	 */
	int getHitLength(bam1_t* b1, bam1_t* b2, int hitPosition, bool isReversed, bool calculateTranscriptPos) {
		int length;
		if (isReversed) {
			length = b->core.pos + b->core.l_qseq - b2->core.pos;
		} else {
			length = b2->core.pos + b2->core.l_qseq - b->core.pos;
		}

		// In the case where we are searching for overlap between a read's position
		// and a transcript's position (when we are using the genome file), the "hitPosition"
		// of a read can be negative if it starts before the transcript. We subtract the 
		// position of the read from the length, and later set the negative position to
		// zero to reflect the read's position in the transcript.
		if (calculateTranscriptPos && hitPosition < 0) {
			length += hitPosition;
		}

		return length;
	}

	template<typename ReadType>
	int parseRead(ReadType&, PairedEndHit&, int&, bool);

	template<typename ReadType>
	int parseRead(ReadType&, SingleHit&, int&, bool);

	void createHit(SingleHit&, int, bool);
	void createHit(PairedEndHit&, int, bool);
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

template<class ReadType, class HitType>
int SamParser::parseNext(ReadType& read, HitType& hit) {
	int sid;
	return parseRead(read, hit, sid, transcripts.getIsUsingGenomeFile());
}

template<class ReadType, class HitType>
int SamParser::parseNext(ReadType& read, HitType& hit, HitType& genomeHit) {
	int sid;
	int result = parseRead(read, hit, sid, true);
	createHit(genomeHit, sid, false);

	return result;
}

// Parse read in "single" read type
template<typename ReadType>
int SamParser::parseRead(ReadType& read, SingleHit& hit, int& sid, bool calculateTranscriptPos) {
	int val;

	if (sam_read1(sam_in, header, b) < 0) return -1;

	std::string name = bam_get_canonical_name(b);
	
	general_assert(!bam_is_paired(b), "Read " + name + ": Find a paired end read in the file!");

	int readType = getReadType(b);
	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		createSingleRead(read, name, b);
	}
	else {
		general_assert(read.getReadLength() == b->core.l_qseq, "Read " + name + " has alignments with inconsistent read lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!\n");
		sid = transcripts.getInternalSid(b, header->target_name[b->core.tid]);
		if (sid == -1) return 99;

		createHit(hit, sid, calculateTranscriptPos);
	}

	return val;
}

// Parse read in "paired" read type
template<typename ReadType>
int SamParser::parseRead(ReadType& read, PairedEndHit& hit, int& sid, bool calculateTranscriptPos) {
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
		createPairedRead(read, name, name2, b, b2);
	}
	else {
		general_assert(read.getMate1().getReadLength() == b->core.l_qseq && read.getMate2().getReadLength() == b2->core.l_qseq, "Paired-end read " + name + " has alignments with inconsistent mate lengths!");
		val = 5;
	}

	if (readType == 1) {
	  general_assert(bam_check_cigar(b) && bam_check_cigar(b2), "Read " + name + ": RSEM currently does not support gapped alignments, sorry!");
	  general_assert(b->core.tid == b2->core.tid, "Read " + name + ": The two mates do not align to a same transcript! RSEM does not support discordant alignments.");

		sid = transcripts.getInternalSid(b, header->target_name[b->core.tid]);
		if (sid == -1) return 99;

		createHit(hit, sid, calculateTranscriptPos);
	}
	
	return val;
}

void SamParser::createHit(SingleHit &hit, int sid, bool calculateTranscriptPos) {
	if (bam_is_rev(b)) {
		hit = SingleHit(
			-sid,
			getHitPosition(b, sid, true, calculateTranscriptPos)
		);
	}
	else {
		hit = SingleHit(
			sid,
			getHitPosition(b, sid, false, calculateTranscriptPos)
		);
	}
}

void SamParser::createHit(PairedEndHit &hit, int sid, bool calculateTranscriptPos) {
	bool isReversed = bam_is_rev(b);
	int position = getHitPosition(b, sid, isReversed, calculateTranscriptPos);
	int length = getHitLength(b, b2, position, isReversed, calculateTranscriptPos);

	hit = PairedEndHit(
		isReversed ? -sid : sid,
		max(position, 0), // make sure the position is "bounded" to the transcript
		length
	);
}

void createSingleRead(SingleRead& read, std::string name, bam1_t* b) {
	read = SingleRead(name, bam_get_read_seq(b));
}

void createSingleRead(SingleReadQ& read, std::string name, bam1_t* b) {
	read = SingleReadQ(name, bam_get_read_seq(b), bam_get_qscore(b));
}

void createPairedRead(PairedEndReadQ& read, std::string name1, std::string name2, bam1_t* b1, bam1_t* b2) {
	SingleReadQ mate1(name1, bam_get_read_seq(b1), bam_get_qscore(b1));
	SingleReadQ mate2(name2, bam_get_read_seq(b2), bam_get_qscore(b2));
	read = PairedEndReadQ(mate1, mate2);
}

void createPairedRead(PairedEndRead& read, std::string name1, std::string name2, bam1_t* b1, bam1_t* b2) {
	SingleRead mate1(name1, bam_get_read_seq(b1));
	SingleRead mate2(name2, bam_get_read_seq(b2));
	read = PairedEndRead(mate1, mate2);
}

#endif /* SAMPARSER_H_ */
