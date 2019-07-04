#include "htslib/sam.h"
#include "sam_utils.h"

class TranscriptsGenome : public Transcripts {
public:
	bool getIsUsingGenomeFile() {
		return isUsingGenomeFile;
	}

	int getInternalSid(bam1_t* b, std::string target_name) {
		if (getIsUsingGenomeFile()) {
			return e2iString[buildExternalKey(b, target_name)];
		} else {
			return Transcripts::getInternalSid(b->core.tid);
		}
	}

	void buildMappings(int, char**, const char*);
    void buildMappingsGenome(samFile*, bam_hdr_t* , const char* = NULL);
	std::map<std::string, ITNode*> createIntervalTrees(std::vector<Transcript>, int);

private:
	bool isUsingGenomeFile = false;
	std::map<std::string, int> e2iString;

	std::string buildExternalKey(bam1_t* b, std::string target_name) {
		return target_name + itos(b->core.pos);
	}
};

void TranscriptsGenome::buildMappings(int n_targets, char** target_name, const char* imdName) {
	isUsingGenomeFile = false;
	Transcripts::buildMappings(n_targets, target_name, imdName);
}

void TranscriptsGenome::buildMappingsGenome(samFile* sam_in, bam_hdr_t* header, const char* imdName) {
	std::vector<bool> appeared;
	i2e.assign(M + 1, 0);
	appeared.assign(M + 1, false);
	e2iString.clear();
	isUsingGenomeFile = true;

	// TODO: Build interval trees. Each interval tree should be for a specific
	// chromosome. chr -> intervalTree
	std::map<std::string, ITNode*> intervalTrees = createIntervalTrees(transcripts, M);

	bam1_t *b = bam_init1();
	while(sam_read1(sam_in, header, b) > -1) {
		// We sequentially read each genome mapping from the sam file, creating an 
		// index from that mapping (i) to our internal ids (in the dict)
		// TODO: We cannot just use target_name[i] like the other function, as target_name[i]
		// is actually the name of the chromosome. Instead, we need to use an interval tree (?)
		// to look up the transcript name for the particular genome mapping `b`
		if (!bam_is_proper(b)) {
			continue;
		}

		std::string chromosome = header->target_name[b->core.tid];
		std::string externalId = buildExternalKey(b, chromosome);

		ITInterval interval = {};
		interval.low = b->core.pos;
		interval.high = interval.low + b->core.l_qseq;

		ITNode* node = overlapSearch(intervalTrees[chromosome], interval);
		if (node == NULL) {
			continue; // TODO: Is this an error condition?
		}

		Transcript* transcript = &transcripts[node->transcriptId];
		std::string transcriptName = isAlleleSpecific() ? transcript->getSeqName() : transcript->getTranscriptID();

		e2iString[externalId] = node->transcriptId;
		appeared[node->transcriptId] = true;
	}

	if (imdName != NULL) {
	  char omitF[STRLEN];
	  sprintf(omitF, "%s.omit", imdName);
	  FILE *fo = fopen(omitF, "w");
	  for (int i = 1; i <= M; i++) 
	    if (!appeared[i]) fprintf(fo, "%d\n", i);
	  fclose(fo);
	}
}

std::map<std::string, ITNode*> TranscriptsGenome::createIntervalTrees(std::vector<Transcript> transcripts, int tLength) {
	std::map<std::string, ITNode*> trees;
	for (int i = 1; i <= tLength; i++) {
		Transcript* transcript = &transcripts[i];
		std::string chromosome = transcript->getSeqName();
		std::vector<Interval> structure = transcript->getStructure();

		for(std::vector<Interval>::iterator it = structure.begin(); it != structure.end(); ++it) {
			ITInterval interval = {};
			interval.low = (*it).start;
			interval.high = (*it).end;

			trees[chromosome] = insert(
				trees.find(chromosome) == trees.end() ? NULL : trees[chromosome],
				interval,
				i
			);
		}
	}

	return trees;
}
