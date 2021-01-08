import random

import numpy as np


class VeryRandom:
    @staticmethod
    def create(seed: int = None):
        """
        :param seed: optional random number generator seed
        Providing seed is not cryptographically secure and is intended mainly for testing.
        """
        if seed is None:
            # use cryptographycally secure random generator
            rnd = random.SystemRandom()
            np_rnd = np.random.seed(rnd.randint(0, 2 ** 32))
        else:
            # use seedable random generator
            rnd = random.Random()
            rnd.seed(seed)
            np_rnd = np.random.RandomState(seed)

        return VeryRandom(rnd, np_rnd)

    @staticmethod
    def seed_rng(seed: int) -> random.Random:
        return random.Random(seed)

    # TODO remove optional parameters
    def __init__(self, rnd: random.Random = None, np_rnd: np.random.RandomState = None):
        """
        Use create() method to create new instance. Direct call to constructor is intended for testing.
        :param rnd:
        :param np_rnd:
        """
        if rnd is None:
            self._rnd = random.Random()
        else:
            self._rnd = rnd

        if np_rnd is None:
            self._np_rnd = np.random.RandomState()
        else:
            self._np_rnd = np_rnd

    def random(self) -> float:
        return self._rnd.random()

    def rand_bytes(self, n) -> bytes:
        """
        :param n: number of random bytes to return
        :return:
        """
        return bytes([self._rnd.getrandbits(8) for _ in range(n)])

    def rand_int(self, low=0, high=2 ** 32) -> int:
        return self._rnd.randint(low, high)

    def rand_ints(self, low, high, n) -> np.ndarray:
        """
        Uses Numpy random generator for maximum performance.
        
        timeit.timeit("sorted(random.sample(range(int(3e9)), int(1e5)))", "import random", number=10)
        timeit.timeit("np.sort(np.random.randint(0, int(3e9), int(1e5)))", "import numpy as np", number=10)
        
        :param low:
        :param high:
        :param n: number of random integers to return
        :return:
        """
        return self._np_rnd.randint(low, high, n)

    def multirand_index(self, p_dist: list) -> int:
        """
        Draw index from multinomial probability distribution.
        :param p_dist: Probability distribution. When all probabilities are zero, each outcome has equal probability.
        Each value represents probability of one outcome.
        :return: Always returns index of outcome in p_dist array.
        """
        if len(p_dist) == 1:
            # only one choice
            return 0

        if sum(p_dist) == 0:
            # each outcome has equal probability
            return self._rnd.randint(0, len(p_dist) - 1)

        p_level = 0
        rnd_value = self._rnd.random()

        # make relative
        p_value = rnd_value * sum(p_dist)
        for i in range(len(p_dist)):
            p_level += p_dist[i]
            if p_level > p_value:
                # random outcome has been reached
                return i

        raise ValueError(
            "sum of probability distribution %d must be greater then the probability value %d" % (sum(p_dist), p_value))

    def shuffle(self, x: list):
        self._rnd.shuffle(x, self._rnd.random)
