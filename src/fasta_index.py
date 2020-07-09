from src.po import FastaSequence


class FastaIndex:

    # TODO factory methods: from_bam / from_fai / from_vcf

    def __init__(self, bam_header: dict):
        """
        :param bam_header: bam header parsed by pysam
        """
        indices = []
        start = 0
        counter = 0

        for record in bam_header['SQ']:
            # ref_name = record['SN'] if keep_chr else strip_chr(record['SN'])
            ref_name = record['SN']
            indices.append(FastaSequence(index=counter, name=ref_name, start=start, length=record['LN']))
            start += record['LN']
            counter += 1

        if len(indices) == 0:
            raise ValueError("Empty BAM sequence header")

        self._indices = indices
        self._dict = dict((reference.name, reference) for reference in self._indices)

    # def __getitem__(self, index: int):
    #     return self._indices[index]

    def __iter__(self):
        for item in self._indices:
            yield item

    def __len__(self):
        return len(self._indices)

    # def __next__(self) -> Optional[FaiRecord]:
    #     try:
    #         return next(self._indices)
    #     except StopIteration:
    #         return None

    def resolve_start_index(self, start_ref_name, start_ref_pos):
        if start_ref_name is None and start_ref_pos is not None:
            raise ValueError("Start reference must be supplied along with start position.")

        if start_ref_name is None:
            start_ref_name = self.first_ref().name

        if start_ref_pos is None:
            start_ref_pos = 0

        return self.pos2index(start_ref_name, start_ref_pos)

    def resolve_end_index(self, end_ref_name, end_ref_pos):
        if end_ref_name is None and end_ref_pos is not None:
            raise ValueError("End reference must be supplied along with end position.")

        if end_ref_name is None:
            end_ref_name = self.last_ref().name

        if end_ref_pos is None:
            end_ref_pos = self.last_ref().length - 1

        return self.pos2index(end_ref_name, end_ref_pos)

    def first_ref(self) -> FastaSequence:
        return self._indices[0]

    def last_ref(self) -> FastaSequence:
        return self._indices[len(self._indices) - 1]

    def first_index(self):
        return self.first_ref().start

    def last_index(self):
        fai_ref = self.last_ref()
        return fai_ref.start + fai_ref.length - 1

    def resolve_start_pos(self, start_index) -> (str, int):
        if start_index is None:
            return self.index2pos(self.first_index())
        else:
            return self.index2pos(start_index)

    def resolve_end_pos(self, end_index) -> (str, int):
        if end_index is None:
            return self.index2pos(self.last_index())
        else:
            return self.index2pos(end_index)

    def ref_id(self, ref_name):
        for i in range(len(self._indices)):
            if self._indices[i].name == ref_name:
                return i

        raise ValueError("Reference name not found in BAM sequence header.")

    def ref_name(self, ref_id):
        return self._indices[ref_id].name

    def index2pos(self, index) -> (str, int):
        """
        Convert absolute position (index) on genome to reference position.
        :param index: 0-based position on genome
        :return: reference position as tuple (<reference name>, <0-based position>)
        """
        ref_pos = index
        for ref in self._indices:
            new_ref_pos = ref_pos - ref.length
            if new_ref_pos < 0:
                return ref.name, ref_pos
            else:
                ref_pos = new_ref_pos
        raise ValueError("reference position for index %d not found" % index)

    def pos2index(self, ref_name, ref_pos) -> int:
        """
        Convert reference position to absolute position (index) on genome.
        :param ref_name: reference name
        :param ref_pos: 0-based reference position
        :return: 0-based position on genome
        """
        return self._dict[ref_name].start + ref_pos
