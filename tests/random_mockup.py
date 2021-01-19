from src.very_random import VeryRandom


class VeryRandomMockup(VeryRandom):
    def __init__(self):
        super().__init__()
        self._multirand_counter = 0
        self._random_counter = -1

    def random(self) -> float:
        self._random_counter += 1
        return self._random_counter % 2

    def rand_int(self, **kwargs):
        return 0

    def multirand_index(self, p_dist: list) -> int:
        """
        return non zero index in order from _dist between multiple calls
        :param p_dist:
        :return:
        """
        not_zero_indices = []
        for i in range(len(p_dist)):
            if p_dist[i] > 0:
                not_zero_indices.append(i)

        curr_index = self._multirand_counter
        self._multirand_counter += 1
        if self._multirand_counter == len(not_zero_indices):
            self._multirand_counter = 0

        return not_zero_indices[curr_index]
