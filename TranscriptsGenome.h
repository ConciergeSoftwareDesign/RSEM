#include "htslib/sam.h"
#include "sam_utils.h"

class TranscriptsGenome : public Transcripts {
public:
	bool getIsUsingGenomeFile() {
		return isUsingGenomeFile;
	}

	int getInternalSid(bam1_t* b, std::string target_name) {
		if (getIsUsingGenomeFile()) {
			int internalId = searchInternalId(b, target_name);
			if (internalId < 0) {
				return INFINITY;
			} else {
				return internalId;
			}
		} else {
			return Transcripts::getInternalSid(b->core.tid);
		}
	}

	int searchInternalId(bam1_t* b, std::string chromosome) {
		ITInterval interval = {};
		interval.low = b->core.pos;
		interval.high = interval.low + b->core.l_qseq;

		ITNode* node = overlapSearch(intervalTrees[chromosome], interval);
		if (node == NULL) {
			return -1;
		}

		return node->transcriptId;
	}

	void buildMappings(int, char**, const char*);
  void buildMappingsGenome(samFile*, bam_hdr_t* , const char* = NULL);
	void setHasAppeared(int, bool);
	void writeOmitFile(const char*);

	std::map<std::string, ITNode*> createIntervalTrees(std::vector<Transcript>, int);

private:
	bool isUsingGenomeFile = false;
	std::vector<bool> appeared;
	std::map<std::string, ITNode*> intervalTrees;
	std::map<std::string, int> externalKeyToInternalId;

	std::string buildExternalKey(bam1_t* b, std::string target_name) {
		return target_name + itos(b->core.pos);
	}
};

void TranscriptsGenome::buildMappings(int n_targets, char** target_name, const char* imdName) {
	isUsingGenomeFile = false;
	Transcripts::buildMappings(n_targets, target_name, imdName);
}

void TranscriptsGenome::buildMappingsGenome(samFile* sam_in, bam_hdr_t* header, const char* imdName) {
	i2e.assign(M + 1, 0);
	appeared.assign(M + 1, false);
	isUsingGenomeFile = true;

	// Build a map of chromosome -> interval tree
	intervalTrees = createIntervalTrees(transcripts, M);
}

void TranscriptsGenome::setHasAppeared(int index, bool hasAppeared) {
	appeared[index] = hasAppeared;
}

void TranscriptsGenome::writeOmitFile(const char* filename) {
	if (filename != NULL && isUsingGenomeFile) {
	  char omitF[STRLEN];
	  sprintf(omitF, "%s.omit", filename);
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
